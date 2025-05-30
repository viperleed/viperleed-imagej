import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.io.*;
import ij.measure.*;
import ij.text.*;
import ij.plugin.*;
import ij.util.Tools;
import ij.util.FloatArray;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.File;

/** This class contains the dialog starting the interactive I(V) Curve Editor
 *  (The Curve editor itself in class LeedCurvePlotEditor)
 * 
 *  Implementation note: In contrast to the Spot Tracker, which is a singleton, there
 *  can be multiple instances of the I(V) Curve Editor (the LeedCurvePlotEditor class).
 *  Therefore, the LeedCurvePlotEditor must not read the parameter values
 *  from the LeedParams or ImageJ Prefs but one has to pass the parameters
 *  to the constructor of the LeedCurvePlotEditor.
 */

/** This ImageJ plugin is part of the ViPErLEED package for LEED I(V) analysis.
 *  Licensed under GNU General Public License v3.0 or later (GPL-3.0-or-later),
 *  https://www.gnu.org/licenses/gpl-3.0.html
 *  The authors may decide later to put part of the auxiliary code in this work into the public domain,
 *  to allow incorporation into ImageJ if desired (ImageJ is in the public domain).
 *  The help text is licensed under the Creative Commons Attribution 4.0 (CC BY 4.0) license.
 *  When using and/or modifying this program for scientific work, please cite
 *  the paper describing it:
 *  M. Schmid, F. Kraushofer, A. M. Imre, T. Kißlinger, L. Hammer, U. Diebold, and M. Riva,
 *  ViPErLEED package II: Spot tracking, extraction and processing of I(V) curves,
 *  Phys. Rev. Research, 2024. 
 *  @author Michael Schmid, IAP/TU Wien, 2020-2024
 */

public class LEED_Curve_Editor implements PlugIn, DialogListener, ActionListener {

    static final String PLUGIN_NAME = LeedCurvePlotEditor.PLUGIN_NAME;
    static final String PATTERN_FILE_BUTTON = "Spot Pattern File...";
    static final String EDIT_FILE_BUTTON = "Previous Edit File...";
    static final String  NONE = "*NONE*"; //edit file label is none
    //prefs keys for parameters
    static final String PREFS_KEY_PREFIX = "leed.editor_";
    static final String PREFS_KEY_SMOOTHEV = "smoothev";
    static final String PREFS_KEY_MINSPANEV = "minspanev";
    static final String PREFS_KEY_CURVESTARTEV = "curvestartev";
    static final String PREFS_KEY_CURVEENDEV = "curveendev";
    static final String PREFS_KEY_NOISELIMIT = "noiselimit";
    //default parameter values that we don't share with the Spot Tracker (eV except noise limit)
    static final double DEFAULT_SMOOTHEV = 5;
    static final double DEFAULT_MINSPANEV = 30;
    static final double DEFAULT_CURVESTARTEV = 50;
    static final double DEFAULT_CURVEENDEV = 1000;
    static final double DEFAULT_NOISELIMIT = 0.05;

    private String pathPrefix;
    private String patternFilePath;
    private String editFilePath;
    LeedSpotPattern spotPattern;
    LeedIVData ivData;

    GenericDialog gd;
    Label patternLabel, editLabel;
    Checkbox autoSelectCbx;
    private static final String HELP_STRING =
            "<html>"+
            "<h1>Parameters and Settings for the ViPErLEED I(V) Curve Editor</h1>"+
            "<ul>"+
            "<li><b>V0i</b>: Imaginary part of the inner potential  (strictly speaking, its absolute value), needed for calculating the <i>R</i> factor. "+
            "When averaging symmetry-equivalent curves that do not all include the full energy range, "+
            "V0i also determines the width of the transition zone with smooth fade in/fade out. "+
            "If V0i is not known, use, e.g., a value of 5 eV.</li>"+
            "<li><b>Default energy range to smooth</b>: How much the curves should be smoothed; this is the "+
            "energy range of a moving average with same noise suppression, "+
            "typically 0.6&nbsp;V0i (very good data or data including low energies &lt; 50 eV) to 1.3&nbsp;V0i (noisy data). "+
            "If a <a href='#editFile'>previous edit file</a> is used, "+
            "this setting has no effect and the value from the previous session is used.</li>"+
            "<li><b>Min. data range per curve</b>: Curves with fewer contiguous data points than this are ignored. "+
            "Note that you cannot make ignored curves (with too few points) visible any more in the current editing session. "+
            "You have to close and open the Curve Editor to get them back.</li>"+
            "<li><b>Lowest energy to use</b> and <b>Highest energy to use</b>: "+
            "The default energy range  if a beam is selected for the first time. "+
            "These parameters also limit the range for automatic (noise-dependent) selection of energy ranges. "+
            "The limits do not restrict manual selection in the I(V) curve editor. "+
            "If a <a href='#editFile'>previous edit file</a> is used, "+
            "the settings for the default energy range from the previous session are used, "+
            "also for curves not selected in that session.</li>"+
            "<li><b>Select curves/ranges automatically</b>: Automatically selects which groups of symmetry-equivalent <i>I</i>(<i>V</i>) curves "+
            "and which energy ranges will be selected. This selection is mainly based on an estimation of the noise and "+
            "its impact on the <i>R</i> factor.</li>"+
            "<li><b>Noise limit</b>: Automatic selection is done such that the estimated impact on the <i>R</i> factor is less than "+
            "the noise limit given. Typical values are around 0.05; lower values are more selective, higher values "+
            "lead to a larger total energy range (larger 'database'), at the cost of increased noise.</li>"+
            "<li><b>Pattern file</b>: Usually, the information on symmetry-equivalent beams is read from the headings "+
            "of the file with the <i>I</i>(<i>V</i>) curves (spot groups in square brackets). If this is not the case or you want to use a different "+
            "symmetry, a spot pattern file can be specified. "+
            "A spot pattern file is a comma-separated list (.csv file) of beam names, <i>h</i>, <i>k</i>, <i>gx</i>, <i>gy</i> "+
            "(reciprocal lattice vector in Cartesian coordinates), and beam group number (equal for symmetry-equivalent beams). "+
            "Such a file can be created with the pattern simulator of the ViPErLEED GUI.</li>"+
            "<li><a name='editFile'><b>Previous edit file</b></a>: If these <i>I</i>(<i>V</i>) data have been edited previously, one may select the file from the "+
            "previous session to continue with the same selection of curves, energy ranges and smoothing. "+
            "If an edit file is selected and you want a fresh start (not continue the session stored in a that edit file) "+
            "press the button to open the file dialog and then &lt;Cancel&gt;.<br>&nbsp;&nbsp;&nbsp;&nbsp;"+
            "If 'Select curves/ranges automatically' is checked, automatic selection will supersede "+
            "the selection in the edit file; only the smoothing (and deselection of individual curves from the average) "+
            "will be read from the edit file.</li>"+
            "</ul>"+
            "<h2><a name='ivEditorLicense'>License</a></h2>"+
            "<p>The code is licensed under <a href='http://www.gnu.org/licenses/gpl-3.0.html'>GNU General Public License v3.0</a> "+
            "or later (GPL-3.0-or-later).</p>"+
            "<p>&nbsp;&nbsp;&nbsp;&nbsp;The ViPErLEED ImageJ plugin collection is free software: you can redistribute it and/or modify it "+
            "under the terms of the GNU General Public License as published by the Free Software Foundation, "+
            "either version 3 of the License, or (at your option) any later version.<br>&nbsp;&nbsp;&nbsp;&nbsp;"+
            "The ViPErLEED ImageJ plugin collection is distributed in the hope that it will be useful, "+
            "but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. "+
            "See the GNU General Public License for more details.<br>&nbsp;&nbsp;&nbsp;&nbsp;"+
            "You should have received a copy of the GNU General Public License along with these plugins. "+
            "If not, see <a href='https://www.gnu.org/licenses/'>https://www.gnu.org/licenses/</a>.</p>"+
            "<p>The authors may decide later to put part of the auxiliary code in this work into the public domain, "+
            "to allow incorporation into ImageJ if desired (ImageJ is in the public domain).</p>"+
            "<p>The documentation, including the help texts, is licensed under the "+
            "<a href='http://creativecommons.org/licenses/by/4.0/'>Creative Commons Attribution 4.0</a> "+
            "(CC BY 4.0) license.</p>"+
            "<p>When using this program (in its original or modified form) for scientific work, "+
            "please cite the paper describing the program [<a href='#paper'>1</a>].</p>"+
            "<p>A copy of these license texts and the source code is included in the jar archive holding this plugin "+
            "(use an unzip utility to view its contents).</p>"+
            "<p><a name='paper'>[1]</a> M. Schmid, F. Kraushofer, A. M. Imre, T. Kißlinger, L. Hammer, U. Diebold, and M. Riva, "+
            "<i>ViPErLEED package II: Spot tracking, extraction and processing of I(V) curves</i>, "+
            "Phys. Rev. Research, 2024. <a href='https://arxiv.org/abs/2406.18413/'>arXiv:2406.18413</a></p>"+
            "</html>";
    // The following is for Compile&Run in ImageJ 1.x only, to ensure everything is updated
    // if it is not directly called here.
    // (we do not check for classes called by the LEED_Spot_Tracker, assuming these are up to date)
    LeedPlotter dummy1; LeedUtils dummy2; LeedIntegerArray dummy3; LeedPlotter dummy4;
    LeedRFactor dummy5; LeedSmoother dummy6; LEED_Spot_Tracker dummy7; LeedIVData dummy8;
    LeedInterpolator dummy9; LeedCurveAverager dummy10;

    /** This constructor may be used in the future */
    /* unused public LEED_Curve_Editor(String pathPrefix, LeedSpotPattern spotPattern, LeedIVData ivData) {
        this.pathPrefix = pathPrefix;
        this.spotPattern = spotPattern;
        this.ivData = ivData;
    } /**/

    /** Empty constructor needed for use as ImageJ plugin [before calling run(String arg)] */
    public LEED_Curve_Editor() {}

    public void run(String arg) {
        if (IJ.versionLessThan("1.53i")) return;

        // P R E P A R E
        String directory = null;
        String filename = null;
        if (ivData == null) {        //Ask for files with intensities
            OpenDialog od = new OpenDialog("Select ..._IVcurves.csv File with I(V) Curves");
            directory = od.getDirectory();
            filename = od.getFileName();
            if (filename == null) return;  //cancelled
            if (!filename.endsWith(".csv") && !IJ.showMessageWithCancel("Leed Curve Editor WARNING",
                    filename+" is not named .csv\nOpen anyhow?")) return;
            IJ.showStatus("Reading data...");
            IJ.showProgress(0.1);
            ivData = LeedIVData.fromFile(directory+filename, null,
                LeedIVData.E_ASCENDING|LeedIVData.ZERO_IS_NAN); //IJ.error on failure
            IJ.showProgress(1.0);
            if (ivData == null) return;
            boolean ok = ivData.trimData();
            if (!ok) {
                IJ.error(PLUGIN_NAME, "Error: no valid data in file\n"+filename);
                return;
            }
            ivData.trimEnergies();

            pathPrefix = directory + File.separator+LeedUtils.removeExtension(filename);
            if (pathPrefix.endsWith("_IVcurves"))
                pathPrefix = pathPrefix.substring(0, pathPrefix.length()-9);
            else if (pathPrefix.endsWith("_Int"))
                pathPrefix = pathPrefix.substring(0, pathPrefix.length()-4);

            spotPattern = new LeedSpotPattern(ivData.spotNames, /*groupsRequired=*/true);

            if (spotPattern.size() <= 0) spotPattern = null;
            if (spotPattern == null) {      //could not get spotPattern from groups in []; try getting it from the log file
                String logFilePath = pathPrefix+"_log.txt";
                if (LeedUtils.fileOk(logFilePath)) {            //open log file to extract path to spot pattern file
                    String log = IJ.openAsString(logFilePath);
                    if (log != null && !log.startsWith("Error:")) {
                        String[] logLines = Tools.split(log,"\n");
                        final String key="leed.tracker_s0=";
                        for (String logLine : logLines)
                            if (logLine.startsWith(key)) {
                                patternFilePath = logLine.substring(key.length());
                                break;
                            }
                    }
                    if (!LeedUtils.fileOk(patternFilePath))
                        patternFilePath = null;
                }
            }
        }

        //default parameter values
        double v0i = LeedParams.get(LeedParams.V0I);
        double smoothEv = Prefs.get(PREFS_KEY_PREFIX+PREFS_KEY_SMOOTHEV, DEFAULT_SMOOTHEV);
        double minSpanEv = Prefs.get(PREFS_KEY_PREFIX+PREFS_KEY_MINSPANEV, DEFAULT_MINSPANEV);
        double curveStartEv = Prefs.get(PREFS_KEY_PREFIX+PREFS_KEY_CURVESTARTEV, DEFAULT_CURVESTARTEV);
        double curveEndEv = Prefs.get(PREFS_KEY_PREFIX+PREFS_KEY_CURVEENDEV, DEFAULT_CURVEENDEV);
        double noiseLimit = Prefs.get(PREFS_KEY_PREFIX+PREFS_KEY_NOISELIMIT, DEFAULT_NOISELIMIT);

        final double eStep = LeedUtils.getEnergyStep(ivData.energies);
        if (v0i < eStep) v0i = eStep;


        if (patternFilePath != null)
            readPatternFile(false);

        IJ.showStatus("");
        File editFile = getEditFile(pathPrefix);
        editFilePath = editFile==null ? null : editFile.getPath();

        // I N I T I A L   D I A L O G
        gd = new GenericDialog("LEED I(V) Curve Editor");
        File dir = (new File(pathPrefix)).getParentFile();
        if (!dir.canWrite())
            gd.addMessage("WARNING: No write access to directory of input file:\n"+dir+"\nSaving a new edit file will fail!");
        gd.setInsets(5, 0, 5); //top left bottom
        gd.addMessage("I(V) Curves from file:\n"+filename);
        gd.addNumericField("V0i", v0i, 1, 6, "eV (min."+IJ.d2s(eStep+0.0499999,1)+")");
        gd.addNumericField("Default energy range to smooth", smoothEv, 1, 6, "eV * **");
        gd.addNumericField("Min. data span per curve", minSpanEv, 0, 6, "eV (min. 2*V0i)");
        gd.addNumericField("Lowest energy to use", curveStartEv, 0, 6, "eV *");
        gd.addNumericField("Highest energy to use", curveEndEv, 0, 6, "eV *");
        gd.addCheckbox("Select curves/ranges automatically", editFilePath==null && !filename.toLowerCase().contains("theobeams"));
        autoSelectCbx = (Checkbox)(gd.getCheckboxes().get(0));
        gd.addNumericField("Noise limit for automatic selection", noiseLimit, 3, 6, " (typ. 0.02-0.08)");
        gd.addMessage("* Used only if no previous edit file is selected\n"+
                "    (for energy range, see 'Help' for details)\n"+
                "** Energy range of moving average with same noise suppression,\n"+
                "    typically 0.6*V0i (very good data or energies < 50 eV) to 1.3*V0i (noisy data).");
        patternLabel = addButtonAndLabel(gd, PATTERN_FILE_BUTTON, "");
        setPatternLabel();
        editLabel = addButtonAndLabel(gd, EDIT_FILE_BUTTON, editFile==null ? NONE : editFile.getName());
        gd.addDialogListener(this);
        gd.addHelp(HELP_STRING);
        // S H O W   D I A L O G
        gd.showDialog();
        if (gd.wasCanceled()) return;

        v0i = Math.abs(gd.getNextNumber());
        if (v0i >= eStep)
            LeedParams.set(LeedParams.V0I, v0i);
        else
            v0i = Math.max(LeedParams.get(LeedParams.V0I), eStep);
        smoothEv = gd.getNextNumber();
        if (smoothEv>=0)
            Prefs.set(PREFS_KEY_PREFIX+PREFS_KEY_SMOOTHEV, smoothEv);
        else
            smoothEv = 0;
        minSpanEv = gd.getNextNumber();
        if (minSpanEv >= 2*v0i)
            Prefs.set(PREFS_KEY_PREFIX+PREFS_KEY_MINSPANEV, minSpanEv);
        else
            minSpanEv = Math.max(Prefs.get(PREFS_KEY_PREFIX+PREFS_KEY_MINSPANEV, DEFAULT_MINSPANEV), 2*v0i);
        curveStartEv = gd.getNextNumber();
        if (curveStartEv >= 0)
            Prefs.set(PREFS_KEY_PREFIX+PREFS_KEY_CURVESTARTEV, curveStartEv);
        else
            curveStartEv = 0;
        curveEndEv = gd.getNextNumber();
        if (curveEndEv >= curveStartEv + minSpanEv)
            Prefs.set(PREFS_KEY_PREFIX+PREFS_KEY_CURVEENDEV, curveEndEv);
        else
            curveEndEv = 1000000;
        boolean autoSelect = gd.getNextBoolean();
        noiseLimit = gd.getNextNumber();
        if (noiseLimit >= 0.001 && noiseLimit <= 0.3)
            Prefs.set(PREFS_KEY_PREFIX+PREFS_KEY_NOISELIMIT, noiseLimit);

        if (spotPattern != null)
            ivData.setSpotPattern(spotPattern, /*showWarning=*/true);         //creates ivData.spotIndices

        while (ivData.spotIndices == null) {            //spotPattern does not fit or we have no spotPattern so far
            askForFile(/*spotPatternWanted=*/true);
            if (spotPattern == null)
                return;
            ivData.setSpotPattern(spotPattern, /*showWarning=*/true);
        }
        if (ivData.spotIndices == null) return;
        fillShortGaps(ivData, LeedCurvePlotEditor.INTERPOLATE_OVER_V0I*v0i/eStep);
        removeShortData(ivData, spotPattern, minSpanEv/eStep, LeedCurvePlotEditor.MIN_RANGE_OVER_V0I*v0i/eStep);
        if (ivData.spotIndices.length == 0) {
            IJ.error(PLUGIN_NAME, "No data with contiguous energy range >="+
                    IJ.d2s(Math.max(minSpanEv,LeedCurvePlotEditor.MIN_RANGE_OVER_V0I*v0i),1));
            return;
        }

        // S T A R T   P L O T   E D I T O R
        LeedCurvePlotEditor curveEditor = new LeedCurvePlotEditor(pathPrefix, spotPattern, ivData, editFilePath,
                v0i, smoothEv, minSpanEv, curveStartEv, curveEndEv, noiseLimit, autoSelect); //no single-instance-listener, so LeedCurveEditor must not read prefs
        new Thread(curveEditor, "LEED_Curve_Editor").start();
    }

    /** Adds a panel with a button and a label, returns a reference to the label */
    Label addButtonAndLabel(GenericDialog gd, String buttonLabel, String labelText) {
        Panel panel = new Panel();
        panel.setLayout(new FlowLayout(FlowLayout.LEFT, 0, 0));
        Button button = new Button(buttonLabel);
        button.addActionListener(this);
        panel.add(button);
        Label theLabel = new Label(labelText);
        panel.add(theLabel);
        gd.addPanel(panel);
        return theLabel;
    }

    /** This callback method is called when the user changes choices or text fields the dialog. */
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        try {
            if (e != null && e.getSource().equals(autoSelectCbx))
                setEditLabel();
            double v0i = Math.abs(gd.getNextNumber());
            double smoothEv = gd.getNextNumber();
            double minSpanEv = gd.getNextNumber();
            double curveStartEv = gd.getNextNumber();
            double curveEndEv = gd.getNextNumber();
            double noiseLimit = gd.getNextNumber();
            if (!(v0i > 1)) return false;
            if (smoothEv < 0) return false;
            if (!(minSpanEv >= 2*v0i)) return false;
            if (!(curveStartEv > 0)) return false;
            if (!(curveEndEv > curveStartEv)) return false;
            if (!(noiseLimit >= 0.001 && noiseLimit <= 0.3)) return false;
            return true;
        } catch (Exception ex) { IJ.handleException(ex); return false; }
    }

    /** Button pressed: asks for Pattern file or Edit file */
    public void actionPerformed(ActionEvent e) {
        boolean spotPatternWanted = e.getActionCommand().equals(PATTERN_FILE_BUTTON);
        askForFile(spotPatternWanted);
    }

    /** Asks for the spotPattern file */
    void askForFile(boolean spotPatternWanted) {
        String path = spotPatternWanted ? patternFilePath : editFilePath;
        if (spotPatternWanted && (path == null || path.length() == 0))
            path = LeedParams.getString(LeedParams.PATTERNFILE);
        String dialogHeading = spotPatternWanted ? "Select Spot Pattern .csv File" : "Select Previous Edit File (cancel=none)";
        OpenDialog od;
        if (path != null && path.length() > 0) {
            File file = new File(path);
            od = new OpenDialog(dialogHeading, file.getParent(), file.getName());
        } else
            od = new OpenDialog(dialogHeading);
        if (spotPatternWanted) {
            patternFilePath = od.getPath();
            if (patternFilePath == null)
                spotPattern = null;
            else
                readPatternFile(true);
            setPatternLabel();
        } else {
            boolean wasEditFile = editFilePath != null;
            editFilePath = od.getPath();
            boolean isEditFile = editFilePath != null;
            if (isEditFile != wasEditFile)
                autoSelectCbx.setState(!isEditFile);
            setEditLabel();
        }
    }

    /** Sets the label next to the edit file button */
    void setEditLabel() {
        String editLabelText = editFilePath == null ?
                NONE : (new File(editFilePath)).getName();
        if (autoSelectCbx.getState() && editFilePath != null)
            editLabelText += " (selection superseded by 'auto')";
        editLabel.setText(editLabelText);
        gd.pack();
    }

    /** Gets the spotPattern from the pattern file */
    void readPatternFile(boolean showErrors) {
        spotPattern = new LeedSpotPattern(patternFilePath);
        if (spotPattern.size() < 1 || spotPattern.getErrorMessage() != null) {
            if (showErrors)
                IJ.error("Pattern invalid or empty\n"+spotPattern.getErrorMessage());
            spotPattern = null;
            patternFilePath = null;
        } else if (showErrors && LeedUtils.fileOk(patternFilePath))      //interactively specified pattern file successfully read
            LeedParams.setString(LeedParams.PATTERNFILE, patternFilePath);
    }

    /** Sets the label next to the button for the spotPattern file */
    void setPatternLabel() {
        if (patternLabel == null) return;
        String label = spotPattern == null ? "<<< SELECT" :
                (patternFilePath == null ? "<read from I(V) curves>" : (new File(patternFilePath)).getName());
        patternLabel.setText(label);
        gd.pack();
    }

    /** Fills short gaps in the data by interpolation */
    static void fillShortGaps(LeedIVData ivData, double interpolateOverPoints) {
        int maxInterpolate = (int)Math.ceil(interpolateOverPoints+1e-6);
        LeedLinearRegression regression = new LeedLinearRegression();
        for (int ic=0; ic<ivData.spotIndices.length; ic++) {
            if (ivData.spotIndices[ic] < 0) continue;
            double[] data = ivData.data[ic];
            int gapLength = 1000000;  //at the start, we can't interpolate, so assume a large gap
            for (int iE=0; iE<ivData.energies.length; iE++) {
                if (!Double.isNaN(data[iE])) {   //valid data
                    if (gapLength <= maxInterpolate) {  //...after a short gap: interpolate
                        for (int jE=Math.max(iE-gapLength-maxInterpolate/2-1, 0); jE<Math.min(iE+maxInterpolate/2+1, ivData.energies.length); jE++)
                            if (!Double.isNaN(data[iE]))
                                regression.addPoint(jE, data[jE]);
                        double slope = regression.getSlope();
                        double offset = regression .getOffset();
                        for (int jE=iE-gapLength; jE<iE; jE++)
                            data[jE] = offset + jE * slope;
                        regression.clear();
                    }
                } else
                    gapLength++;
            }
        }
    }

    /** Deletes columns where we don't have enough contiguous data spanning at least 'minRangeInColumn'
     *  within the column and 'minRangeInGroup' (in points), taking all columns of the group together.
     *  (any data within the group of symmetry-equivalent spots counts, no check for overlap)
     *  Replaces the intData and spotIndices arrays if columns have to be deleted.
     *  Minimum requirement is 2 points. */
    static void removeShortData(LeedIVData ivData, LeedSpotPattern spotPattern, double minRangeInGroup, double minRangeInColumn) {
        int[] spotIndices = ivData.spotIndices;
        int minRangeInCol = (int)Math.ceil(minRangeInColumn + 1e-6);
        if (minRangeInCol < 2) minRangeInCol = 2;
        int minRangeInGrp = (int)Math.ceil(minRangeInGroup + 1e-6);
        for (int ic=0; ic<spotIndices.length; ic++) {  // Check minERangeInColumn
            double[] data = ivData.data[ic];
            if (spotIndices[ic] < 0) continue;
            int nOk = 0, nBestOk = 0;
            for (int iE=0; iE<ivData.energies.length; iE++) {
                if (!Double.isNaN(data[iE]))
                    nOk++;
                else if (nOk > nBestOk) {
                    nBestOk = nOk;
                    nOk = 0;
                }
            }
            if (nOk > nBestOk)
                nBestOk = nOk;
            if (nBestOk < minRangeInCol)
                spotIndices[ic] = -1;
        }
        if (minRangeInGrp > minRangeInCol) {        // Check minERangeInGroup
        boolean[] columnChecked = new boolean[spotIndices.length];
            for (int ic=0; ic<spotIndices.length; ic++) {
                if (spotIndices[ic] < 0) continue;
                if (columnChecked[ic]) continue;
                int group = spotPattern.getGroup(spotIndices[ic]);
                int[] groupIndices = spotPattern.getAllSpotsForGroup(group);
                LeedIntegerArray groupColumns = new LeedIntegerArray(groupIndices.length);
                for (int spotInGroup : groupIndices) {
                    int spotColumn = LeedUtils.arrayIndexOf(spotIndices, spotInGroup);
                    if (spotColumn >= 0)
                        groupColumns.add(spotColumn);
                }
                int nOk = 0, nBestOk = 0;
                for (int iE=0; iE<ivData.energies.length; iE++) {
                    boolean ok = false;
                    for (int icg=0; icg<groupColumns.size(); icg++) {
                        double[] dataColumn = ivData.data[groupColumns.get(icg)];
                        if (!Double.isNaN(dataColumn[iE])) {
                            ok = true; break;
                        }
                    }
                    if (ok)
                        nOk++;
                    else if (nOk > nBestOk) {
                        nBestOk = nOk;
                        nOk = 0;
                    }
                }
                if (nOk > nBestOk)
                    nBestOk = nOk;
                for (int icg=0; icg<groupColumns.size(); icg++) {
                    int column = groupColumns.get(icg);
                    columnChecked[column] = true;
                    if (nBestOk < minRangeInGrp)
                        spotIndices[column] = -1;
                }
            }
        }
        ivData.deleteIfSpotIndex(-1);
    }

    /** Returns the best guess for a file with previous editing results, null if we have one */
    File getEditFile(String pathPrefix) {
        File prefixFile = new File(pathPrefix);
        String prefix = prefixFile.getName()+"_Edit";
        File directoryF = prefixFile.getParentFile();
        if (directoryF == null) return null;
        File[] files = directoryF.listFiles();
        String nameFound = "";
        File fileFound = null;
        for (File file : files) {
            String name = file.getName();
            if (name.startsWith(prefix) && name.compareTo(nameFound) > 0) {
                fileFound = file;
                nameFound = name;
            }
        }
        return fileFound;
    }
}
