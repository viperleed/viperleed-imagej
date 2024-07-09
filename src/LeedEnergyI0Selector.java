import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.io.*;
import ij.plugin.*;
import ij.text.*;
import ij.measure.ResultsTable;
import ij.util.Tools;
import ij.plugin.frame.Recorder;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.File;

/** This class reads metadata such as energies, I0 and I00 (if present) from
 *  the 'Set Energies, I0, t' user dialog and/or the slice labels of an image stack */

/** This code is part of the ViPErLEED package for LEED I(V) analysis.
 *  Licensed under GNU General Public License v3.0 or later (GPL-3.0-or-later),
 *  https://www.gnu.org/licenses/gpl-3.0.html
 *  The authors may decide later to put part of the auxiliary code in this work into the public domain,
 *  to allow incorporation into ImageJ if desired (ImageJ is in the public domain).
 *  When using and/or modifying this program for scientific work, please cite
 *  the paper describing it:
 *  M. Schmid, F. Kraushofer, A. M. Imre, T. Kißlinger, L. Hammer, U. Diebold, and M. Riva,
 *  ViPErLEED package II: Spot tracking, extraction and processing of I(V) curves,
 *  Phys. Rev. Research, 2024. 
 *  @author Michael Schmid, IAP/TU Wien, 2020-2024
 */


public class LeedEnergyI0Selector implements DialogListener {
    static final int N_DIALOG_DATA = 4; //number of data columns: ENERGY, TIME, CURRENTI0, CURRENTI00; we don't handle the others
    static final int ENERGY=LEED_Spot_Tracker.ENERGY, TIME=LEED_Spot_Tracker.TIME,
            CURRENTI0=LEED_Spot_Tracker.CURRENTI0, CURRENTI00=LEED_Spot_Tracker.CURRENTI00; //indices of the output arrays
    static final int[] keyParams = new int[] {          //these entries of LeedParams give the search strings for the columns
            LeedParams.ENERGYKEYS, LeedParams.TIMEKEYS, LeedParams.CURRENTI0KEYS, LeedParams.CURRENTI00KEYS,
            LeedParams.TEMPERATUREKEYS, LeedParams.AUXKEYS, -1, -1, -1, -1};
    static final String[] LONG_NAMES = new String[] {"Energies", "Frame Times", "Beam currents I0", "Beam-off currents I00"};
    static final int NONE=0, STACK=1, TABLE=2, MANUAL=3, PREV=4, DARKI0=5, DARKI0_LIN=6;  //sources of the data
    static final String[] SOURCE_NAMES = new String[] {
            "None", "Input Stack", "Table", "Manual input", "Previous data", "Dark frame I0", "Dark frame I0, linear fit"};
    static final String[] SOURCE_SHORT = new String[] {"[?]", "", "[t]", "[!]", "[?]", "[D]", "[d]"};
    static final String TABLE_PREFIX = "Table: ";

    static final double MIN_SMOOTHING = 1.5;            //minimum I0 smoothing (points of comparable moving average filter)

    static int[] dataSources = new int[N_DIALOG_DATA];  //data sources actually used
    static String[] lastTableColumn = new String[N_DIALOG_DATA]; //last name in a table

    //global for use during dialog callback
    int nSlices;
    Hashtable<String, ResultsTable> tableList;
    double[] values1 = new double[N_DIALOG_DATA];       //'manual' values for each data field
    double[] values2 = new double[N_DIALOG_DATA];
    double smoothingI0Points = Double.NaN;              //equivalent smoothing points of moving average
    double[][] dataFromStack, previousData, newData;
    double[] darkI0 = null;
    int[] dialogDataSources = numberToDataSources((int)LeedParams.get(LeedParams.METADATASOURCES));
    int xAxisVariable;
    String[] macroSources = new String[N_DIALOG_DATA];  //Strings with data sources for macro call
    LeedLabeledField[] valueFields = new LeedLabeledField[2*N_DIALOG_DATA];
    LeedLabeledField smoothingI0Field;
    Label messageLabel;
    Checkbox showPlotCbx;
    Choice xAxisChoice;


    /** The dialog asking for energies etc. First tries to read the data from
     *  the slice labels of the stack; the user may change these data and/or
     *  provide new input.
     *  The result, as far as provided, is written to energiesEtc; for
     *  data not provided, the respective array in energiesEtc remains unchanged.
     *  The data in energiesEtc remain unchanged on cancel or error. */
    public void runSelectionDialog(LEED_Spot_Tracker spotTracker, ImageStack stack, ImageStack darkStack, double[][] energiesEtc) {
        try {
            double[][] dataFromStack = decode(/*spotTracker=*/null, stack); //also sets 'dataSources' from LeedParams
            if (stack == null)  {
                spotTracker.enableAndHighlightComponents(true);
                return;  //leaves it unchanged
            }
            if (darkStack != null)
                darkI0 = decode(/*spotTracker=*/null, darkStack)[CURRENTI0];
            nSlices = stack.size();
            this.dataFromStack = dataFromStack;
            this.previousData = energiesEtc;
            if (LeedParams.getBoolean(LeedParams.LEEMMODE) && //in LEEM mode, energy was saved as E_LEEM; handle as ENERGY again
                    previousData[LEED_Spot_Tracker.E_LEEM] != null && previousData[LEED_Spot_Tracker.ENERGY] == null) {
                previousData[LEED_Spot_Tracker.ENERGY] = previousData[LEED_Spot_Tracker.E_LEEM];
                previousData[LEED_Spot_Tracker.E_LEEM] = null;
            }
            newData = new double[LEED_Spot_Tracker.N_DATA_E_ETC][];
            tableList = getTables();                //find open ResultsTables
            String[] tableNames = tableList.keySet().toArray(new String[0]);
            Arrays.sort(tableNames);
            GenericDialog gd = new GenericDialog(LEED_Spot_Tracker.PLUGIN_NAME+" - Energies, I0, Time");
            for (int d=0; d<N_DIALOG_DATA; d++) {
                if (d <= CURRENTI00) { //in case of more than 4 items: we have v1, v2 only for the first 4 in LeedParams
                    values1[d] = LeedParams.get(LeedParams.EMANUALFIRST +2*d);
                    values2[d] = LeedParams.get(LeedParams.EMANUALLAST +2*d);
                }
                if (!isCurrent(d) && dialogDataSources[d] != MANUAL)
                    makeValues1Values2(d, NONE);
                ArrayList<String> options = new ArrayList<String>();
                boolean offerPrevious = previousData != null && previousData[d] != null;    //we want a 'use previous' option
                if (dataFromStack[d] != null) {
                    options.add(SOURCE_NAMES[STACK]);
                    if (dataSources[d] == STACK) offerPrevious = false;
                }
                if (d==CURRENTI00 && darkI0 != null) {
                    options.add(SOURCE_NAMES[DARKI0]);
                    if (dataSources[d] == DARKI0) offerPrevious = false;
                }
                if (d==CURRENTI00 && darkI0 != null && darkI0.length==stack.size()) {
                    options.add(SOURCE_NAMES[DARKI0_LIN]);
                    if (dataSources[d] == DARKI0_LIN) offerPrevious = false;
                }
                if (offerPrevious)
                    options.add(SOURCE_NAMES[PREV]);
                for (int i=0; i<tableNames.length; i++)
                    options.add(TABLE_PREFIX+tableNames[i]);
                options.add(SOURCE_NAMES[MANUAL]);
                options.add(SOURCE_NAMES[NONE]);
                String[] optionArray = options.toArray(new String[0]);
                gd.addChoice(LONG_NAMES[d]+" from:", optionArray, defaultOption(d, optionArray, dialogDataSources[d]));
                gd.addToSameRow();
                gd.addNumericField(isCurrent(d) ? "I@E=0" : "start", values1[d], isCurrent(d) ? 4 : 1);
                valueFields[2*d] = LeedLabeledField.numeric(gd);
                gd.addToSameRow();
                gd.addNumericField(isCurrent(d) ? "slope" : "last", values2[d], isCurrent(d) ? 6 : 1);
                valueFields[2*d+1] = LeedLabeledField.numeric(gd);
            }
            double smoothI0points = LeedParams.get(LeedParams.SMOOTHI0POINTS);
            if (smoothI0points < MIN_SMOOTHING) smoothI0points = 0;
            gd.addNumericField("I0 smoothing (points)", smoothI0points, 0);
            smoothingI0Field = LeedLabeledField.numeric(gd);
            gd.setInsets(0,20,0); //top, left, bottom
            gd.addCheckbox("Plot I0", false);
            showPlotCbx = (Checkbox)(gd.getCheckboxes().lastElement());
            gd.setInsets(10,20,0); //top, left, bottom
            gd.addCheckbox("LEEM mode", LeedParams.getBoolean(LeedParams.LEEMMODE));
            Checkbox leemModeCbx = (Checkbox)(gd.getCheckboxes().lastElement());
            String[] xAxisChoices = getXAxisChoices(dataFromStack);
            int defaultXAxis = (int)LeedParams.get(LeedParams.XAXISVARIABLE);
            gd.addChoice("X axis variable", xAxisChoices, LEED_Spot_Tracker.ENERGY_ETC_NAMES[defaultXAxis]);
            xAxisChoice = (Choice)gd.getChoices().lastElement();

            gd.addMessage("");
            messageLabel = (Label)gd.getMessage();
            for (int d=0; d<N_DIALOG_DATA; d++)
                adjustDialog(gd, (Choice)(gd.getChoices().get(d)));  //prepare dialog according to all choice states
            gd.addDialogListener(this);
            gd.addHelp(LeedSpotTrackerHelp.getEnergyI0Help());
            //spotTracker.setCurrentFrontDialog(gd);

            gd.showDialog();
            spotTracker.enableAndHighlightComponents(true);

            if (gd.wasCanceled()) return;

            for (int d=0; d<N_DIALOG_DATA; d++) {
                if (dialogDataSources[d] == MANUAL &&
                        !Double.isNaN(values1[d]) && !Double.isNaN(values2[d])) {
                    if (LeedParams.get(LeedParams.EMANUALFIRST + 2*d) != values1[d] ||
                            LeedParams.get(LeedParams.EMANUALLAST  + 2*d) != values2[d])
                        spotTracker.setPostTrackingChanges(LONG_NAMES[d]);
                    LeedParams.set(LeedParams.EMANUALFIRST + 2*d, values1[d]);
                    LeedParams.set(LeedParams.EMANUALLAST  + 2*d, values2[d]);
                }
            }

            if (LeedParams.get(LeedParams.SMOOTHI0POINTS) != smoothingI0Points)
                spotTracker.setPostTrackingChanges("I0 smoothing");
            if (smoothingI0Field.isEnabled() && !Double.isNaN(smoothingI0Points))
                LeedParams.set(LeedParams.SMOOTHI0POINTS, smoothingI0Points);

            if (LeedParams.getBoolean(LeedParams.LEEMMODE) != leemModeCbx.getState())
                spotTracker.changeStatus(0, LEED_Spot_Tracker.TRACK_OK);
            LeedParams.set(LeedParams.LEEMMODE, leemModeCbx.getState());

            if (xAxisVariable < 0) {    //no suitable column, selected 'index': create aux column
                newData[LEED_Spot_Tracker.AUX] = LeedUtils.createSequence(stack.size());
                xAxisVariable = LEED_Spot_Tracker.AUX;
            }
            int dataSourcesNumber = dataSourcesToNumber(dialogDataSources);
            if (LeedParams.get(LeedParams.XAXISVARIABLE) != xAxisVariable ||
                    LeedParams.get(LeedParams.METADATASOURCES) != dataSourcesNumber)
                spotTracker.changeStatus(0, LEED_Spot_Tracker.TRACK_OK);            
            LeedParams.set(LeedParams.XAXISVARIABLE, xAxisVariable);
            LeedParams.set(LeedParams.METADATASOURCES, dataSourcesNumber);

            tempDataToFinal(newData, energiesEtc, stack.size(), dialogDataSources);
            if (LeedParams.getBoolean(LeedParams.LEEMMODE)) {
                setLeemMode(energiesEtc);
                if (xAxisVariable == LEED_Spot_Tracker.ENERGY && energiesEtc[LEED_Spot_Tracker.E_LEEM] != null)
                    xAxisVariable = LEED_Spot_Tracker.E_LEEM;
            }
            spotTracker.setXAxisVariable(xAxisVariable);
            spotTracker.updateEnergiesEtc(/*checkProcessedI0=*/true, /*showPlot=*/energiesEtc[CURRENTI0] != null && showPlotCbx.getState());

            if (Recorder.record) {
                for (int d=0; d<N_DIALOG_DATA; d++) //record SetArray from table column
                    if (macroSources[d] != null)
                        LEED_Spot_Tracker.recordMacro(LEED_Spot_Tracker.MACRO_COMMAND_SETARRAY,
                                '"'+LEED_Spot_Tracker.ENERGY_ETC_NAMES[d]+"\","+macroSources[d]);
                LEED_Spot_Tracker.recordMacro(LEED_Spot_Tracker.MACRO_COMMAND_UPDATE, "'metadata'");
            }
        } catch (Exception e) { IJ.handleException(e); }
        return;
    }

   /** This callback method is called when the user changes choices or text fields the dialog.
    *  It also populates the newData array, which will be used at the end (except for <cancel>) */
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        try {
            Object src = e == null ? null : e.getSource();  //e=null on final call, when dialog is closed
            if (src instanceof Choice && src != xAxisChoice)
                adjustDialog(gd, (Choice)src);              //enable/disable & set numeric fields, ask in case of table
            for (int d=0; d<N_DIALOG_DATA; d++) {
                String option = gd.getNextChoice();
                int iOption = LeedUtils.arrayIndexOf(SOURCE_NAMES, option);
                double[] data = null;
                double v1 = gd.getNextNumber();             //always read it (to stay in sync), even if we don't use it
                double v2 = gd.getNextNumber();
                if (iOption < 0) {                          //iOption<0 if data from a table, these have been read already into newData[d]
                    data = newData[d];
                } else if (iOption == STACK) {
                    data = dataFromStack[d];
                    macroSources[d] = null;
                } else if (iOption == PREV) {
                    data = previousData[d];
                    macroSources[d] = null;
                } else if (iOption == MANUAL) {
                    if (Double.isNaN(v1) || Double.isNaN(v2)) {
                        messageLabel.setText("BAD NUMBER: "+LONG_NAMES[d]);
                        return false;
                    }
                    if (d == ENERGY) {
                        if (!energyOK(v1) || !energyOK(v2)) {
                            messageLabel.setText("BAD ENERGY");
                            return false;
                        }
                    }
                    values1[d] = v1;
                    values2[d] = v2;
                    Object dataOrString = makeLinearRamp(d, nSlices, v1, v2, newData[ENERGY]);
                    if (dataOrString instanceof String)  //error
                        messageLabel.setText((String)dataOrString);
                    else
                        data = (double[])dataOrString;
                    macroSources[d] = null;
                } else if (iOption == DARKI0 || iOption == DARKI0_LIN) {
                    Object dataOrString = getI00FromDark(nSlices, darkI0, iOption);
                    if (dataOrString instanceof String)  //error
                        messageLabel.setText((String)dataOrString);
                    else
                        data = (double[])dataOrString;
                }
                newData[d] = data;
                if (data == null)
                    macroSources[d] = "newArray(0)"; //none
                //IJ.log(d+" src="+iOption+" macroSrc="+macroSources[d]);
            }
            for (int d=N_DIALOG_DATA; d<Math.min(newData.length, dataFromStack.length); d++)
                newData[d] = dataFromStack[d];              //the remaining data always come from the stack

            smoothingI0Points = gd.getNextNumber();
            if (smoothingI0Field.isEnabled() && !(smoothingI0Points >= 0)) {
                messageLabel.setText("ERROR: I0 smoothing points must be >= 0");
                return false;
            }
            String[] xAxisChoices = getXAxisChoices(newData);   //data may have changed, update options for x axis
            LeedUtils.setChoiceItems(xAxisChoice, xAxisChoices);
            String xAxisStr = gd.getNextChoice();
            xAxisVariable = LeedUtils.arrayIndexOf(LEED_Spot_Tracker.ENERGY_ETC_NAMES, xAxisStr); //maybe -1 for 'index'

            messageLabel.setText(""); //no error
        } catch (Exception ex) { IJ.handleException(ex); return false; }
        return true;
    }    

    /** Converts the data sources to a single number for the LeedParams:
     *  Least-siginificant digit for energy, then t, I0, I00; 0-4 for none/stack/table/man/prev */
    private static int dataSourcesToNumber(int[] dataSources) {
        int out = 0;
        for (int i=dataSources.length-1; i>=0; i--)
            out = 10*out + dataSources[i];
        return out;
    }

    /** Converts a single number from the LeedParams to an array of data sources. Number format:
     *  Least-siginificant digit for energy, then t, I0, I00; 0-4 for none/stack/table/man/prev */
    private static int[] numberToDataSources(int num) {
        int[] out = new int[N_DIALOG_DATA];
        for (int i=0; i<N_DIALOG_DATA; i++) {
            int div10 = num/10;
            out[i] = num - 10*div10; //modulo
            num = div10;
        }
        return out;
    }

    /** Sets the data source type for the given column, for documentation purposes.
     *  Does not modify the data.
     *  Used when set via a macro to an array (source=PREV) or null (source=NONE). */
    public static void setDataSource(int dataColumn, int source) {
        dataSources[dataColumn] = source;
    }

    /** Creates a linear ramp for the given dataColumn.
     *  For I0 and I00, v1 and v2 are intercept and slope of the function w.r.t. energy
     *  (This requires that the energy column is set unless the slope v2 is 0).
     *  Otherwise v1 & v2 are start and final values.
     *  Returns the data array if successful, otherwise an error String starting with 'ERROR:' */
    private static Object makeLinearRamp(int dataColumn, int nSlices, double v1, double v2, double[] energies) {
        if (Double.isNaN(v1) || Double.isNaN(v2))
            return "ERROR: Invalid NaN value for "+LEED_Spot_Tracker.ENERGY_ETC_SHORTNAMES[dataColumn]+" linear ramp";
        double[] data = new double[nSlices];
        if (isCurrent(dataColumn)) {
            if (energies == null && v2 != 0) {
                return "ERROR: no energy for "+LEED_Spot_Tracker.ENERGY_ETC_SHORTNAMES[dataColumn]+"(E) = ";
            } else {
                for (int i=0; i<nSlices; i++) {
                    data[i] = v1;
                    if (v2 != 0) data[i] += energies[i]*v2; //if v2==0, energies may be null
                }
            }
        } else {
            for (int i=0; i<nSlices; i++)
                data[i] = v1 + i*(v2-v1)/(nSlices-1);
        }
        return data;
    }

    /** Returns the I00 values from 'darkI0', the I0 of the dark frame, where 'method' can be
     *  'DARKI0' (value-by-value if the size of darkI0 fits 'size' the main stack 'size', average otherwise)
     *  or 'DARKI0_LIN' (linear fit, requires that the size of darkI0 fits 'size' the main stack 'size'.
     *  Instead returns an error String starting with 'ERROR:' if data are not appropriate */
    private static Object getI00FromDark(int size, double[] darkI0, int method) {
        if (darkI0 == null) return "ERROR: Dark frame has no I0, cannot use it as I00";
        double[] i00Values = null;
        if (darkI0.length==size) {
            if (method == DARKI0) {
                i00Values = darkI0;
            } else if (method == DARKI0_LIN) {
                double[] x = LeedUtils.createSequence(size);
                LeedLinearRegression lin = new LeedLinearRegression();
                lin.addPoints(x, darkI0);
                i00Values = lin.getFitValues(x);
            }
        } else if (method == DARKI0) {
            i00Values = new double[size];
            Arrays.fill(i00Values, LeedUtils.getArraySum(darkI0)/darkI0.length);
        }
        if (i00Values == null)
            return "ERROR: cannot create I00 from "+SOURCE_NAMES[method];
        return i00Values; //success
    }

    /** Returns whether this data column is a current (which may be specified by offset&slope w.r.t. energy) */
    private static boolean isCurrent(int dataColumn) {
        return dataColumn==CURRENTI0 ||dataColumn==CURRENTI00;
    }

    /** Returns a String array with the choices for the x axis if the energy is not a suitable x axis.
     *  Does not check whether the energy is a suitable axis. */
    private static String[] getXAxisChoices(double[][] data) {
        ArrayList<String> choiceList = new ArrayList<String>(10);
        for (int d=0; d<data.length; d++) {
            if (d != LEED_Spot_Tracker.CURRENTI00 && d != LEED_Spot_Tracker.PROCESSEDI0 &&
                    data[d] != null && usableValues(data[d]))
                choiceList.add(LEED_Spot_Tracker.ENERGY_ETC_NAMES[d]);

            //IJ.log(LEED_Spot_Tracker.ENERGY_ETC_NAMES[d]+": data="+data[d]+" usable="+usableValues(data[d]));
        }
        if (choiceList.size() == 0)
            choiceList.add("index");
        return choiceList.toArray(new String[0]);
    }

    /** Whether the values are not all the same */
    private static boolean usableValues(double[] values) {
        if (values == null) return false;
        double[] minMax = Tools.getMinMax(values);
        return minMax[1] > minMax[0];
    }

    /** Updates the dialog according to the Choice for a dataColumn
     *  This is called initially and after user changes of the Choice item.
     *  Also asks for the column and reads the data if a table is provided */
    private void adjustDialog(GenericDialog gd, Choice choice) {
        int dataColumn = gd.getChoices().indexOf(choice);
        String option = choice.getSelectedItem();
        int iOption = LeedUtils.arrayIndexOf(SOURCE_NAMES, option);
        //IJ.log("choice d="+dataColumn+" selSrc="+iOption);
        if (iOption >= 0) {
            dialogDataSources[dataColumn] = iOption;
            setNumericFields(gd, dataColumn);
        } else {                    //data from a table requested
            if (!option.startsWith(TABLE_PREFIX))
                throw new RuntimeException("Unrecognized option: "+option);
            String tableName = option.substring(TABLE_PREFIX.length());
            ResultsTable rt = tableList.get(tableName);
            double[] data = askForTableColumn(gd, dataColumn, tableName, rt);
            if (data != null) {
                newData[dataColumn] = data;
                dialogDataSources[dataColumn] = TABLE;
                makeValues1Values2(dataColumn, TABLE);
                setNumericFields(gd, dataColumn);
            } else {                //non-successful read of column
                iOption = dialogDataSources[dataColumn];
                choice.select(SOURCE_NAMES[iOption]);
            }
        }
        setNumericFieldsEnabled(gd, dataColumn, iOption==MANUAL);
        if (!isCurrent(dataColumn) && (iOption==STACK || iOption==PREV)) {
            makeValues1Values2(dataColumn, iOption);
            setNumericFields(gd, dataColumn);
        }
        if (dataColumn == CURRENTI0) {
            boolean weHaveI0 = iOption != NONE;
            smoothingI0Field.setEnabled(weHaveI0);
            showPlotCbx.setEnabled(weHaveI0);
        }
    }

    /** Sets the numericFields for a data column from the 'firstValues' and 'increments' */
    private void setNumericFields(GenericDialog gd, int dataColumn) {
        Vector numericFields = gd.getNumericFields();
        double v1 = values1[dataColumn];
        double v2 = values2[dataColumn];
        if (Double.isNaN(v1)) {
            ((TextField)(numericFields.get(2*dataColumn))).setText("");
            ((TextField)(numericFields.get(2*dataColumn+1))).setText("");
        } else {
            ((TextField)(numericFields.get(2*dataColumn))).setText(IJ.d2s(v1,1));
            ((TextField)(numericFields.get(2*dataColumn+1))).setText(IJ.d2s(v2, isCurrent(dataColumn) ? 4 : 1));
        }
    }

    /** Enables or disables the numeric fields for a data colunm */
    private void setNumericFieldsEnabled(GenericDialog gd, int dataColumn, boolean b) {
        //IJ.log(dataColumn+" Enable="+b);
        Vector numericFields = gd.getNumericFields();
        ((TextField)(numericFields.get(2*dataColumn))).setEnabled(b);
        ((TextField)(numericFields.get(2*dataColumn+1))).setEnabled(b);
        if (valueFields[2*dataColumn] != null) valueFields[2*dataColumn].setEnabled(b);
        if (valueFields[2*dataColumn+1] != null) valueFields[2*dataColumn+1].setEnabled(b);
    }

    /** Returns a list of all open ResultsTables with suitable size.
     *  The window title is the key String. */
    private Hashtable<String, ResultsTable> getTables() {
        Hashtable<String, ResultsTable> list = new Hashtable<String, ResultsTable>();
        Frame[] nonImageWindows = WindowManager.getNonImageWindows();
        for (Frame win : nonImageWindows) {
            if (!(win instanceof TextWindow)) continue;
            TextPanel tp = ((TextWindow)win).getTextPanel();
            ResultsTable rt = tp.getResultsTable();
            if (rt == null || rt.getColumnHeadings().trim().length()==0) continue;
            if (rt.size() == nSlices)
                list.put(win.getTitle(), rt);
        }
        return list;
    }

    /** Populates the 'firstValues' and 'increments' for a non-current data column,
     *  with source = NONE for automatic (STACK, if possible). */
    private void makeValues1Values2(int dataColumn, int source) {
        if (isCurrent(dataColumn)) return;
        double[] sourceArray = null;
        if (source == NONE)
            sourceArray = dataFromStack[dataColumn] != null ? dataFromStack[dataColumn] :
                (previousData != null ? previousData[dataColumn] : null);
        else if (source == STACK)
            sourceArray = dataFromStack[dataColumn];
        else if (source == PREV && previousData != null)
            sourceArray = previousData[dataColumn];
        else if (source == TABLE)
            sourceArray = newData[dataColumn];
        if (sourceArray != null) {
            double[] firstAndInc = LeedUtils.getFirstAndIncrement(sourceArray);
            values1[dataColumn] = firstAndInc[0];  //NaNs if nonlinear
            values2[dataColumn] = firstAndInc[0] + firstAndInc[1]*(nSlices-1);
            //IJ.log("Preset d="+dataColumn+" src="+source+" v1="+(float)values1[dataColumn]+" v2="+(float)values2[dataColumn]);
            return;
        }
    }

    /** Returns the default option, with a given array of options for data sources.
     *  If the previous choice is not available, selects the first available of 'stack' or 'none' */
    private String defaultOption(int dataColumn, String[] optionArray, int previousSource) {
        String option = null;
        int index = LeedUtils.arrayIndexOf(optionArray, SOURCE_NAMES[previousSource]);
        if (index >= 0) return optionArray[index];   //keep previous choice if available
        index = LeedUtils.arrayIndexOf(optionArray, SOURCE_NAMES[PREV]);
        if (index >= 0 && previousSource==TABLE) return optionArray[index];   //after having used a table, default is 'previous'
        return SOURCE_NAMES[NONE];
    }

    /** Displays a dialog asking for the column of a data table and returns the array */
    private double[] askForTableColumn(GenericDialog gd, int dataColumn, String tableName, ResultsTable rt) {
        String[] headings = rt.getColumnHeadings().split("\t");
        int nColOk = 0;                         //determine which columns are ok (numeric)
        boolean[] colOk = new boolean[headings.length];
        int lastRow = rt.size() - 1;
        for (int i=0; i<headings.length; i++) {
            if (!rt.columnExists(headings[i])) continue;   //ignore empty columns and row label column
            double v = rt.getValue(headings[i], 0);
            colOk[i] = !Double.isNaN(rt.getValue(headings[i], 0)) && !Double.isNaN(rt.getValue(headings[i], lastRow));
            if (colOk[i]) nColOk++;
        }
        if (nColOk == 0) {
            IJ.error(LEED_Spot_Tracker.PLUGIN_NAME, "ERROR: Table does not contain any fully numeric column");
            return null;
        }
        
        String[] columns = new String[nColOk];
        int defaultCol = -1;
        for (int i=0, iOut=0; i<headings.length; i++) {
            if (!colOk[i]) continue;
            columns[iOut] = headings[i];
            if (columns[iOut].equalsIgnoreCase(lastTableColumn[dataColumn]))
                defaultCol = iOut;
            iOut++;
        }
        if (defaultCol < 0)                     //no column with the previous name? try keywords
            defaultCol = getDefaultColumnIndex(dataColumn, columns);

        GenericDialog gd1 = new GenericDialog(LEED_Spot_Tracker.PLUGIN_NAME+" "+SOURCE_NAMES[dataColumn]);
        gd1.addMessage("Select data column in "+tableName+"\nfor "+LEED_Spot_Tracker.ENERGY_ETC_NAMES[dataColumn]);
        gd1.addChoice("Column", columns, columns[defaultCol]);
        gd1.showDialog();
        if (gd1.wasCanceled()) return null;
        String col = gd1.getNextChoice();
        double[] data = rt.getColumn(col);
        if (data != null && data.length != nSlices) {
            IJ.error(LEED_Spot_Tracker.PLUGIN_NAME, "ERROR: Table has wrong length (changed?)");
            return null;
        }
        for (int i=0; i<data.length; i++) {
            if (Double.isNaN(data[i])) {
                IJ.error(LEED_Spot_Tracker.PLUGIN_NAME, "ERROR: Non-numeric value in table row "+i);
                return null;
            }
            if (dataColumn == ENERGY && !energyOK(data[i])) {
                IJ.error(LEED_Spot_Tracker.PLUGIN_NAME, "ERROR: bad energy in table row "+i+": "+IJ.d2s(data[i]));
                return null;
            }
        }
        lastTableColumn[dataColumn] = col;
        macroSources[dataColumn] = "Table.getColumn(\""+col+"\", \""+tableName+"\")";
        return data;
    }

    /** Finds the column with the name best matching one of the keywords for the given data type */
    int getDefaultColumnIndex(int dataType, String[] columns) {
        String[] searchKeys = (String[])(LeedParams.getStringArray(keyParams[dataType]).clone());
        for (int i=0; i<searchKeys.length; i++)
            if (searchKeys[i].endsWith("=") || searchKeys[i].endsWith("|"))
                searchKeys[i] = searchKeys[i].substring(0, searchKeys[i].length()-1);  //remove trailing '=' and ':'
        for (int iOut=0; iOut<columns.length; iOut++)
            for (String key : searchKeys)
                if (columns[iOut].equalsIgnoreCase(key))    //try case-insensitive match first
                    return iOut;
        for (int iOut=0; iOut<columns.length; iOut++)
            for (String key : searchKeys)
                if (columns[iOut].contains(key))            //next try case-sensitive substring
                    return iOut;
        return 0;   //default if nothing found
    }

    /** Updates the arrays according to the LeedParams: This can set linear ramps for I0 or I00
     *  or an input to 'none' or the value from the input stack. It cannot set metadata from a table.
     *  Returns null if successful or an error string starting with 'ERROR:'
     *  Call spottracker.updateEnergiesEtc thereafter. */
    public static String updateFromParams(LEED_Spot_Tracker spotTracker, ImageStack stack, ImageStack darkStack) {
        int[] newDataSources = numberToDataSources((int)LeedParams.get(LeedParams.METADATASOURCES));
        boolean leemMode = LeedParams.getBoolean(LeedParams.LEEMMODE);
        int xAxisVariable = (int)LeedParams.get(LeedParams.XAXISVARIABLE);
        spotTracker.setXAxisVariable(xAxisVariable); //can be needed for spotTracker.getEnergyArray() below
        double[][] dataFromStack = decode(/*spotTracker=*/null, stack);
        double[] darkI0 = null;
        if (darkStack != null)
            darkI0 = decode(/*spotTracker=*/null, darkStack)[CURRENTI0];

        for (int d=0; d<N_DIALOG_DATA; d++) {
            double[] data = null;  //remains null is source==NONE
            int source = newDataSources[d];
            if (source == STACK) {
                if (dataFromStack[d] == null) continue; //leave it untouched if no data
                data = dataFromStack[d];
            } else if (source == TABLE) {
                continue; //can't set from table this way; leave everything untouched
            } else if (source == MANUAL) {
                Object dataOrString = makeLinearRamp(d, stack.size(),
                    LeedParams.get(LeedParams.EMANUALFIRST + 2*d), LeedParams.get(LeedParams.EMANUALLAST + 2*d), spotTracker.getEnergyArray());
                if (dataOrString instanceof String)
                    return (String)dataOrString; //error
                else
                    data = (double[])dataOrString;
            } else if (source == PREV) continue; //leave all untouched
            if (leemMode && d == ENERGY && xAxisVariable == ENERGY) {
                spotTracker.energiesEtc[ENERGY] = null;
                spotTracker.energiesEtc[LEED_Spot_Tracker.E_LEEM] = data;
            } else if (source == DARKI0 || source == DARKI0_LIN) {
                if (d != CURRENTI00)
                    return "ERROR: Invalid metaDataSources, only I0 can be read from dark frame I0";
                Object dataOrString = getI00FromDark(stack.size(), darkI0, source);
                if (dataOrString instanceof String)
                    return (String)dataOrString; //error
                else
                    data = (double[])dataOrString;
            } else
                spotTracker.energiesEtc[d] = data;
            dataSources[d] = data == null ? NONE : source;
        }
        return null;
    }


    /** Decodes the energy, I0 and I00 (measured beam current and offset of beam current measurements)
     *  as well as the time from the slice labels of a stack and darkStack (the latter when supplied
     *  and if I00 should be taken from the dark stack acccording to the LeedParams).
     *  Returns four arrays with these quantities.
     *  For quantities that could not be read the respective array is null.
     *  When spotTracker is non-null, can enable LEEM mode and sets the xAxisVariable of SpotTracker */
    public static double[][] decode(LEED_Spot_Tracker spotTracker, ImageStack stack, ImageStack darkStack) {
        double[][] output = decode(spotTracker, stack);
        int[] theSources = numberToDataSources((int)LeedParams.get(LeedParams.METADATASOURCES));
        int sourceI00 = theSources[CURRENTI00];
        if (darkStack != null && (sourceI00 == DARKI0 || sourceI00 == DARKI0_LIN)) {
            double[] darkI0 = decode(/*spotTracker=*/null, darkStack)[CURRENTI0];
            Object dataOrString = getI00FromDark(stack.size(), darkI0, sourceI00);
            if (dataOrString instanceof String) {
                dataSources[CURRENTI00] = NONE;
                LeedUtils.logError("Error reading I00 from dark: "+(String)dataOrString);
            } else {
                output[CURRENTI00] = (double[])dataOrString;
                dataSources[CURRENTI00] = sourceI00;
            }
        }
        return output;
    }


    /** Decodes the energy, I0 and I00 (measured beam current and offset of beam current measurements)
     *  as well as the time from the slice labels of a stack.
     *  Returns four arrays with these quantities for the stack;
     *  For quantities that could not be read the respective array is null.
     *  When spotTracker is non-null, can enable LEEM mode and sets the xAxisVariable of SpotTracker */
    public static double[][] decode(LEED_Spot_Tracker spotTracker, ImageStack stack) {
        int size = stack == null ? 0 : stack.size();
        if (size == 0)                              //no stack or stack closed? empty output
            return new double[LEED_Spot_Tracker.N_DATA_E_ETC][];
        double[][] output = new double[LEED_Spot_Tracker.N_DATA_E_ETC][size];
        boolean[] hasError = new boolean[LEED_Spot_Tracker.N_DATA_E_ETC];
        String[][] searchKeys = new String[LEED_Spot_Tracker.N_DATA_E_ETC][];
        for (int d=0; d<LEED_Spot_Tracker.N_DATA_E_ETC; d++)
            if (keyParams[d] >= 0)
                searchKeys[d] = LeedParams.getStringArray(keyParams[d]);
        for (int d=0; d<LEED_Spot_Tracker.N_DATA_E_ETC; d++)
            Arrays.fill(output[d], Double.NaN);

        for (int i=0; i<stack.size(); i++) {        //decode all values
            for (int d=0; d<LEED_Spot_Tracker.N_DATA_E_ETC; d++) {
                if (searchKeys[d] == null) continue; //we don't read PROCESSEDI0 from the stack
                String slicelabel = stack.getSliceLabel(i+1);
                output[d][i] = slicelabel == null ? Double.NaN : LeedUtils.getNumberFromString(slicelabel, searchKeys[d]);
                if (Double.isNaN(output[d][i])) {
                    if (IJ.debugMode && !hasError[d]) IJ.log("Error reading "+LEED_Spot_Tracker.ENERGY_ETC_NAMES[d]+" for stack slice "+(i+1)+":\n"+slicelabel);
                    hasError[d] = true;
                }
            }
        }

        if (!hasError[ENERGY]) {                    //check energy values & steps
            double minStep=Double.MAX_VALUE, maxStep=-Double.MAX_VALUE;
            for (int i=0; i<stack.size(); i++) {
                double energy = output[ENERGY][i];
                if (!energyOK(energy))
                    hasError[ENERGY] = true;
                if (i > 0) {
                    double step = energy - output[ENERGY][i-1];;
                    if (step < minStep) minStep = step;
                    if (step > maxStep) maxStep = step;
                }
            }
            if (minStep * maxStep < 0)              //energy not monotonic?
                hasError[ENERGY] = true;
        }
        for (int d=0; d<LEED_Spot_Tracker.N_DATA_E_ETC; d++) {
            if (hasError[d] || (output[d] != null && Double.isNaN(output[d][0])))
                output[d] = null;                   //invalid data arrays become null
        }

        int xAxisVariable = (int)LeedParams.get(LeedParams.XAXISVARIABLE);
        if (xAxisVariable > output.length || output[xAxisVariable] == null)
            xAxisVariable = LEED_Spot_Tracker.ENERGY;

        if (spotTracker != null) {
            if (LeedParams.getBoolean(LeedParams.LEEMMODE)) {
                setLeemMode(output);
                if (xAxisVariable == LEED_Spot_Tracker.ENERGY && output[LEED_Spot_Tracker.E_LEEM] != null)
                    xAxisVariable = LEED_Spot_Tracker.E_LEEM;
            }
            spotTracker.setXAxisVariable(xAxisVariable);
        }

        if (IJ.debugMode)
            for(int d=0; d<output.length; d++)
                if(output[d]!=null)
                    IJ.log("Data from stack: "+LEED_Spot_Tracker.ENERGY_ETC_NAMES[d]+": "+output[d][0]+"..."+output[d][output[d].length-1]);
        return output;
    }

    /** Copies the arrays from tempData to 'energiesEtc'.
     *  The corresponding data source codes can be provided as 'tempSource';
     *  If tempSource is null, the source is assumed to be 'STACK', and only non-null
     *  arrays are copied.
     *  For simplicity, the data source codes are kept as static class variables. */
    public static void tempDataToFinal(double[][] tempData, double[][] energiesEtc, int stackSize, int[] tempSource) {
        for (int d=0; d<LEED_Spot_Tracker.N_DATA_E_ETC; d++) {
            if (tempData[d] != null || tempSource != null) {
                energiesEtc[d] = tempData[d];
                if (d < dataSources.length)
                    dataSources[d] = tempSource == null ? STACK :
                        (energiesEtc[d] == null ? NONE : tempSource[d]);
                    //if (d < dataSources.length) IJ.log(LEED_Spot_Tracker.ENERGY_ETC_SHORTNAMES[d]+" copiedToFinal, src="+dataSources[d]);
            } else if (d < dataSources.length && energiesEtc[d] == null) {     //new input stack, old data kept (nothing new)
                dataSources[d] = PREV;
            }
            if (energiesEtc[d] != null && energiesEtc[d].length != stackSize) {
                energiesEtc[d] = null;  //delete items that do not fit the stack size
                if (d < dataSources.length) {
                    dataSources[d] = NONE;
                    if (IJ.debugMode) IJ.log(LEED_Spot_Tracker.ENERGY_ETC_SHORTNAMES[d]+": bad length");
                }
            }
        }
    }

    /** Returns an array with the processed beam current I0 (i.e. after subtraction
     *  of I00 and smoothing, if selected), or null if none.
     *  Even if no modification of I0 is required, always returns a new array */
    public static double[] getProcessedI0(double[][] energiesEtc) {
        double[] output = energiesEtc[CURRENTI0];
        if (output == null) return null;
        if (energiesEtc[CURRENTI00] != null) {
            output = Arrays.copyOf(output, output.length);
            for (int i=0; i<output.length; i++)
                output[i] -= energiesEtc[CURRENTI00][i];
        }
        double i0Smoothing = LeedParams.get(LeedParams.SMOOTHI0POINTS);
        if (i0Smoothing >= MIN_SMOOTHING) {
            int kernelHalfWidth = LeedSmoother.invNoiseGainSquareToM(i0Smoothing);
            output = (new LeedSmoother(kernelHalfWidth)).smooth(output);
        }
        if (output == energiesEtc[CURRENTI0])
            output = Arrays.copyOf(output, output.length);
        return output;
    }

    /** When I0 smoothing is on, augments the processed I0 by the high-frequency components
     *  from the slice-dependent background of the LEED images.
     *  The correction is applied by modifying the processedI0 array */
    public static void applyI0BackgroundVariations(double[] processedI0, double[] background) {
        if (processedI0 == null || background == null) return;
        double i0Smoothing = LeedParams.get(LeedParams.SMOOTHI0POINTS);
        if (i0Smoothing < MIN_SMOOTHING) return;

        int kernelHalfWidth = LeedSmoother.invNoiseGainSquareToM(i0Smoothing);
        double[] modBg = (new LeedSmoother(kernelHalfWidth)).smooth(background);
        for (int i=0; i<modBg.length; i++) {
            modBg[i] = background[i]/modBg[i]; //correction factor, a bit like high-pass filtered
            if (!(background[i] > 0) || !(modBg[i] > 0.1 && modBg[i] < 10)) {
                IJ.error(LEED_Spot_Tracker.PLUGIN_NAME,
                        "I0 correction from background intensity not performed.\nSlice "+(i+1)+
                        " has background "+LeedUtils.d2s(background[i], 4)+
                        ", correction factor would be "+LeedUtils.d2s(modBg[i], 3));
                return;
            }
        }
        for (int i=0; i<modBg.length; i++)
            processedI0[i] *= modBg[i];
    }

    /** Returns a plot of I0 and processed I0 (if different), or null if none.
     *  This is called with a title containing 'error' if a negative (processed) I0 is detected */
    public static Plot getProcessedI0Plot(String title, double[][] energiesEtc, int xAxisVariable) {
        double[] processedI0 = getProcessedI0(energiesEtc);
        if (processedI0 == null) return null;
        if (energiesEtc[xAxisVariable] == null) return null;
        Plot plot = new Plot(title, LEED_Spot_Tracker.ENERGY_ETC_SHORTNAMES[xAxisVariable], LEED_Spot_Tracker.ENERGY_ETC_SHORTNAMES[CURRENTI0]);
        plot.addPoints(energiesEtc[xAxisVariable], energiesEtc[CURRENTI0], Plot.LINE);
        String legend = LEED_Spot_Tracker.ENERGY_ETC_SHORTNAMES[CURRENTI0];
        if (energiesEtc[CURRENTI00] != null) {
            plot.setColor(Color.GRAY);
            plot.addPoints(energiesEtc[xAxisVariable], energiesEtc[CURRENTI00], Plot.LINE);
            legend += '\t'+LEED_Spot_Tracker.ENERGY_ETC_SHORTNAMES[CURRENTI00];
        }
        if (!processedI0.equals(energiesEtc[CURRENTI0])) {
            plot.setColor(Color.RED);
            plot.setLineWidth(2);
            plot.addPoints(energiesEtc[xAxisVariable], processedI0, Plot.LINE);
            legend += "\tfinal";
        }
        if (title.toLowerCase().contains("error") && processedI0 != null) {
            double[] minmax = Tools.getMinMax(processedI0);
            plot.setLimits(Double.NaN, Double.NaN,      //xmin, xmax auto
                    1.1*minmax[0]-0.1*minmax[1],        //y range: encompass processedI0 and y=0, + 10% padding
                    Math.max(1.1*minmax[1]-0.1*minmax[0],0.1*Math.abs(minmax[0])));
        } else
            plot.setLimits(Double.NaN, Double.NaN, 0, Double.NaN);
        plot.setColor(Color.BLACK);
        plot.setLineWidth(1);
        if (legend.indexOf('\t') > 0)
            plot.addLegend(legend);
        return plot;
    }

    /** Returns the x-axis type, ENERGY or TIME (the latter if energy is constant). */
    /*public static int getXAxisType(double[][] energiesEtc) {
        if (energiesEtc[TIME] == null) return ENERGY;
        if (energiesEtc[ENERGY] == null) return TIME;
        double[] energyMinMax = Tools.getMinMax(energiesEtc[ENERGY]);
        return energyMinMax[0] == energyMinMax[1] ? TIME : ENERGY;
    }*/

    /** Returns a short summary for the field next to the 'Energy Range' button */
    public static String getStatusText(double[][] energiesEtc, int xAxisVariable) {
        double[] minMax = null;
        if (energiesEtc[xAxisVariable] != null)
            minMax = Tools.getMinMax(energiesEtc[xAxisVariable]);
        String str = LEED_Spot_Tracker.ENERGY_ETC_SHORTNAMES[xAxisVariable];
        if (xAxisVariable != ENERGY)
            str += "!";
        if (minMax != null && (xAxisVariable != ENERGY || minMax[0] > 0)) {
            str += LeedUtils.d2s(minMax[0], 4)+"-"+LeedUtils.d2s(minMax[1], 4);
            if (xAxisVariable < dataSources.length)
                str += SOURCE_SHORT[dataSources[xAxisVariable]];
        } else {
            str += "???";
        }
        double smoothingI0 = LeedParams.get(LeedParams.SMOOTHI0POINTS);
        for (int d=0; d<N_DIALOG_DATA; d++)
            if (d != xAxisVariable && energiesEtc[d] != null && !isAllZero(energiesEtc[d])) {
                if (d == CURRENTI00 && energiesEtc[CURRENTI0] == null)
                    continue;
                str += ", " + LEED_Spot_Tracker.ENERGY_ETC_SHORTNAMES[d] + SOURCE_SHORT[dataSources[d]];  //adds, e.g. ', I0[s]
                if (d == CURRENTI0 && smoothingI0 > 0)
                    str += "s"+(int)smoothingI0;
            }
            if (energiesEtc[CURRENTI0] == null)
                str += "; no I0";
        return str;
    }

    /** Sets LEEM mode by setting the E_LEEM array to ENERGIES and the ENERGIES array to null.
     *  Does nothing if the ENERGIES array is null, so it can be called twice if desired */
    static void setLeemMode(double[][] energiesEtc) {
        if (energiesEtc[LEED_Spot_Tracker.ENERGY] != null) {
            energiesEtc[LEED_Spot_Tracker.E_LEEM] = energiesEtc[LEED_Spot_Tracker.ENERGY];
            energiesEtc[LEED_Spot_Tracker.ENERGY] = null;
        }
    }

    /** Checks whether an array has only zero entries */
    static boolean isAllZero(double[] a) {
        for (double v : a)
            if (v != 0) return false;
        return true;
    }

    /** Checks whether a given energy value is valid; returns false if invalid or NaN */
    public static boolean energyOK(double energy) {
        return energy > 0.001 && energy < 1e4;
    }
}
