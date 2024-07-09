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
import java.lang.reflect.Array;

/** Compares two data files and creates a list of mutual R factors.
 *
 *  Interpolates to the finer of the two energy grids or 1 eV, whichever is finest.
 *  Symmetry-equivalent beams are averaged before calculating the R factors; if
 *  averaging is performed, the energy range is limited to the range where all
 *  symmetry-equivalent beams are present.
 *  
 *  Forbidden spots are ignored.
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
 *  @author Michael Schmid, IAP/TU Wien, 2021-2024
 */


public class LEED_R_Factor_Between_Datasets implements PlugIn {
    static final String PLUGIN_NAME = "R Factor Between Data Sets";
    static final double MAX_E_STEP = 0.5;       //energy step for R factor comparison no larger than this
    static final int R_FACTOR_TYPE = LeedRFactor.R_PENDRY;
    //static final int R_FACTOR_TYPE = LeedRFactor.R_PMOD;
    /** The following should be remembered until ImageJ is closed */
    static String[] filePaths = new String[2];
    static String spotPatternPath = LeedParams.getString(LeedParams.PATTERNFILE);
    static double v0i = LeedParams.get(LeedParams.V0I);
    static boolean averageEquivalent = true;
    static boolean allowShift = true;
    static boolean showPlot;
    private static final String HELP_STRING =
            "<html>"+
            "<h1>"+PLUGIN_NAME+"</h1>"+
            "<p>This plugin compares two <i>I</i>(<i>V</i>) data files (.csv format).</p>"+
            "<p>It provides a list with Pendry's R factor for each beam as well as the total <i>R</i> factor, "+
            "which is the weighted average of all beams, with weights according to the energy span of each beam.<br />"+
            "In addition to the <i>R</i> factor, the energy span (where both files contain valid data) "+
            "and ratio between the average intensity per beam of the two files are reported. "+
            "The ratio is normalized such that the overall ratio is 1.0.<br />"+
            "In case of a superstructure (i.e., if at least one non-integer beam is present), "+
            "separate summaries are given for the integer and superstructure beams. "+
            "(If the intensity ratios for superstructure beams strongly differ from integer beams, "+
            "this may be an indication for partial coverage by the superstructure, or, "+
            "when comparing experiment and calculation, for a wrong structure model.)</p>"+
            "<p>Note that the <i>R</i> factors are calculated in a rather simple way, ignoring the end points "+
            "(where the derivatives are poorly defined). Results may slightly differ from those of other <i>R</i> factor programs.</p>"+
            "<p>The option <b>Average symmetry-equivalent beams</b> is mainly meant for averaging calculated beams. "+
            "The experimental <i>I</i>(<i>V</i>) curves should be averaged with the <i>I</i>(<i>V</i>) Curve Editor. "+
            "In contrast to the <i>I</i>(<i>V</i>) Curve Editor, the averaging algorithm used here does not use smooth fade-in/fade-out "+
            "if the curves start/end at different energies; "+
            "the average is rather limited to the common energy range of the input data. "+
            "If a curve is defined over an energy range that is too short (if would reduce the energy range to less than half) it is ignored.<br />"+
            "If averaging is off, comparison is done on a beam-per-beam basis. "+
            "If both of the input files contain only one beam per group of symmetry-equivalent beams, "+
            "comparison will be done for beam groups. In this case, beams will be also compared if the group is the same, "+
            "but the beam indices are different.</p>"+
            "<p><b>Allow shifting</b> searches for the best offset of the energy axes between the two input files "+
            "(the offset is the same for all beams). "+
            "When comparing calculated and experimental <i>I</i>(<i>V</i>) curves, this essentially corresponds to a variation of "+
            "the real part of the inner potential (V_0r) or a variation of the filament work function.</p>"+
            "<p><b>Show I(V) curve plot</b> creates a plot stack, with one plot for each beam, "+
            "showing the respective <i>I</i>(<i>V</i>) curves of the two input files. The curves of the first plot are red, the second black.</p>"+
            "<p>If the groups of symmetry-equivalent beams cannot be read from both input files "+
            "(or if the group numbers in these two files are inconsistent) a <b>spot pattern file</b> is required. "+
            "This is a .csv file containing a list of the beams with additional information on them, "+
            "such as the group number, which defines groups of symmetry-equivalent beams. "+
            "A spot pattern file can be also supplied to override the symmetry assumed in the input data files "+
            "(i.e., the group numbers in square brackets).<br />"+
            "A spot pattern file can be created with the ViPErLEED gui (select 'export').</p>"+
            "&nbsp;<br /><hr />"+
            "<p>The code is licensed under <a href='http://www.gnu.org/licenses/gpl-3.0.html'>GNU General Public License v3.0</a> "+
            "or later (GPL-3.0-or-later). "+
            "The authors may decide later to put part of the auxiliary code in this work into the public domain, "+
            "to allow incorporation into ImageJ if desired (ImageJ is in the public domain).</p>"+
            "<p>This documentation is licensed under the <a href='http://creativecommons.org/licenses/by/4.0/'>Creative Commons Attribution 4.0</a> "+
            "(CC BY 4.0) license.</p>"+
            "<p>When using this program (in its original or modified form) for scientific work, "+
            "please cite the paper describing the program [<a href='#paper'>1</a>].</p>"+
            "<p>You should find a copy of these license texts and the source code in the zip/jar archive holding this plugin "+
            "(use an unzip utility to view its contents).</p>"+
            "<p><a name='paper'>[1]</a> M. Schmid, F. Kraushofer, A. M. Imre, T. Kißlinger, L. Hammer, U. Diebold, and M. Riva, "+
            "<i>ViPErLEED package II: Spot tracking, extraction and processing of <i>I</i>(<i>V</i>) curves</i>, "+
            "Phys. Rev. Research, 2024.</p>"+

            "</html>";

    public void run(String arg) {
        if (IJ.versionLessThan("1.52u")) return;

        // I N I T I A L   D I A L O G
        if (spotPatternPath==null)
            spotPatternPath = LeedParams.getString(LeedParams.PATTERNFILE);
        GenericDialog gd = new GenericDialog("R Factor Between Sets of IV curves");
        gd.addFileField("File_1", filePaths[0], 50);
        gd.addFileField("File_2", filePaths[1], 50);
        gd.addFileField("Spot pattern file", spotPatternPath, 50);
        gd.addMessage(".csv files required; you can drag and drop files into the fields");
        gd.addNumericField("V0i", v0i, 1, 6, "eV");
        gd.addCheckbox("Average symmetry-equivalent beams", averageEquivalent);
        gd.addCheckbox("Allow shifting (V0r variation)", allowShift);
        gd.addCheckbox("Show I(V) curve plot", showPlot);

        DialogListener dialogListener = new DialogListener() {  // for checking whether the files exist and the weight is ok
            public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
                int nFiles = 0;
                for (int i=0; i<3; i++) {
                    String path = gd.getNextString();
                    if (path.trim().length()==0)
                        return i>=2;                            // spot pattern file may be empty
                    if (!LeedUtils.fileOk(path))
                        return false;
                }
                double v0i = Math.abs(gd.getNextNumber());
                if (!(v0i > 0)) return false;
                return true;
            }
        };
        gd.addDialogListener(dialogListener);
        Button okButton = gd.getButtons()[0];
        okButton.setEnabled(dialogListener.dialogItemChanged(gd, null));
        gd.addHelp(HELP_STRING);
        gd.showDialog();
        if (gd.wasCanceled()) return;

        filePaths[0] = gd.getNextString();
        filePaths[1] = gd.getNextString();
        spotPatternPath = gd.getNextString();
        v0i = gd.getNextNumber();
        averageEquivalent = gd.getNextBoolean();
        allowShift = gd.getNextBoolean();
        showPlot = gd.getNextBoolean();

        if (!(v0i > 0)) {
            IJ.error("Error: invalid V0i, positive number required");
            return;
        }

        // P R E P A R E   D A T A
        LeedIVData[] ivData = new LeedIVData[filePaths.length];
        for (int i=0; i<filePaths.length; i++) {
            IJ.showProgress(0.1+0.3*i);
            ivData[i] = LeedIVData.fromFile(filePaths[i], null, LeedIVData.E_ASCENDING | LeedIVData.ZERO_IS_NAN);
            if (ivData[i] == null) {
                IJ.showProgress(1.0);
                return;
            }
        }
        IJ.showProgress(1.0);
        String spotPatternError = null;
        LeedSpotPattern spotPattern = null;
        if (spotPatternPath.trim().length() == 0) {
            LeedSpotPattern[] spotPatterns = new LeedSpotPattern[ivData.length];
            for (int i=0; i<ivData.length; i++) {
                spotPatterns[i] = new LeedSpotPattern(ivData[i].spotNames, /*groupsRequired=*/averageEquivalent);
                spotPatternError = spotPatterns[i].getErrorMessage();
                if (spotPatternError != null)
                    break;
            }
            if (spotPatternError == null) {         //merge the spot patterns read from the files
                spotPattern = new LeedSpotPattern(spotPatterns[0], spotPatterns[1], /*groupsRequired=*/averageEquivalent);
                spotPatternError = spotPattern.getErrorMessage();
            }
            if (spotPatternError != null)
                spotPattern = null;
        }
        do {                 //make sure we have a valid SpotPattern
            if (spotPatternPath != null && spotPatternPath.trim().length() > 0) {
                spotPattern = new LeedSpotPattern(spotPatternPath);
                if (spotPattern.size() < 2) {
                    IJ.error(PLUGIN_NAME, "Error reading Spot Pattern from file:\n"+spotPattern.getErrorMessage());
                    spotPattern = null;
                    spotPatternPath = null;
                }
            } else if (spotPatternError != null) {
                IJ.error(PLUGIN_NAME, "Cannot read information on symmetry-equivalent groups from I(V) curve files:\n"+
                        spotPatternError+
                        "\nYou will be asked to specify the spot pattern file");
                spotPatternError = null;
            }
            if (spotPattern == null) {
                spotPattern = LeedSpotPattern.openWithDialog(LeedParams.getString(LeedParams.PATTERNFILE));
                if (spotPattern == null) return;
            }
            if (spotPattern != null) {
                for (LeedIVData ivd : ivData) {
                    ivd.setSpotPattern(spotPattern, /*showWarning=*/false);
                    if (ivd.spotIndices == null) {
                        spotPattern = null;
                        spotPatternPath = null;
                        break;
                    }
                }
            }
        } while (spotPattern == null);
        if (spotPattern.getPath() != null)
            spotPatternPath = spotPattern.getPath();    //remember for the next invocation


        // C A L C U L A T E
        IJ.showProgress(0.7);
        if (averageEquivalent && spotPattern.hasEquivalentBeams()) {
            for (int i=0; i<ivData.length; i++)
                averageEquivalentBeams(ivData[i], spotPattern);
        }
        Object ivDataOrError = getCommonIVData(ivData, spotPattern);
        if (ivDataOrError instanceof String) {
            IJ.error((String)ivDataOrError);
            IJ.showProgress(1.0);
            return;
        }
        LeedIVData[] commonIVData = (LeedIVData[])ivDataOrError;
        Object rFactorDataOrError = getRFactorData(commonIVData[0], commonIVData[1], v0i, R_FACTOR_TYPE, allowShift, /*summaryOnly=*/false);
        if (rFactorDataOrError instanceof String) {
            IJ.error((String)rFactorDataOrError);
            IJ.showProgress(1.0);
            return;
        }
        double[][] rFactorData = (double[][])rFactorDataOrError;
        IJ.showProgress(0.8);

        double sumRatio = 0;                                            //normalize intensity ratio 2/1
        double sumRange = 0;
        for (int ic=0; ic<commonIVData[0].spotIndices.length; ic++)
            if (rFactorData[ic][LeedRFactor.RATIO] > 0) {
                double overlap = rFactorData[ic][LeedRFactor.N_OVERLAP];
                sumRatio += rFactorData[ic][LeedRFactor.RATIO] * overlap;
                sumRange += overlap;
            }
        for (int ic=0; ic<rFactorData.length; ic++)
            rFactorData[ic][LeedRFactor.RATIO] *= sumRange/sumRatio;

        double[][] sumIntSup = new double[2][LeedRFactor.OUTPUT_SIZE];  //[0] integer, [1] superstructure
        int[] nIntSup = new int[2];
        if (spotPattern.isSuperstructure()) {                           //create separate statistics for integer & superstructure
            for (int ic=0; ic<commonIVData[0].spotIndices.length; ic++) {
                int intOrSup = spotPattern.isSuperstructure(commonIVData[0].spotIndices[ic]) ?
                        1 : 0;
                for (int i=0; i<LeedRFactor.OUTPUT_SIZE; i++)
                    sumIntSup[intOrSup][i] += rFactorData[ic][i];
                nIntSup[intOrSup]++;
            }
        }

        // O U T P U T
        String title = ivData[0].title+"_vs_"+ivData[1].title;
        ResultsTable rt = new ResultsTable();
        Plot plot = showPlot ? new Plot(title, "E", "intensity", Plot.DEFAULT_FLAGS & (~Plot.Y_NUMBERS)) : null;
        if (showPlot) {
            normalizeBeams(commonIVData);
            double xLeft = Math.max(commonIVData[0].energies[0], commonIVData[1].energies[0]);
            double xRight = Math.min(LeedUtils.lastElement(commonIVData[0].energies), LeedUtils.lastElement(commonIVData[1].energies));
            plot.setLimits(1.01*xLeft - 0.01*xRight, 1.01*xRight - 0.01*xLeft, -0.01, 1.01);
        }

        for (int ic=0; ic<rFactorData.length; ic++) {                   //add data to output table&plot
            double[] rFactorEtc = rFactorData[ic];
            if (rFactorEtc[LeedRFactor.N_OVERLAP] == 0) continue;
            boolean isSpotData = ic<commonIVData[0].spotIndices.length; //not the summary ('overall' result)
            String label = isSpotData ? commonIVData[0].spotNames[ic] : "OVERALL";
            rt.incrementCounter();
            rt.addLabel(label);
            rt.addValue(LeedRFactor.R_FACTOR_NAMES[R_FACTOR_TYPE], rFactorEtc[R_FACTOR_TYPE]);
            rt.addValue("E_Range", rFactorEtc[LeedRFactor.N_OVERLAP]);
            rt.addValue("Ratio2/1", rFactorEtc[LeedRFactor.RATIO]);

            if (showPlot && isSpotData) {
                plot.setColor(Color.RED);
                plot.addPoints(commonIVData[0].energies, commonIVData[0].data[ic], Plot.LINE);
                plot.setColor(Color.BLACK);
                plot.addPoints(commonIVData[1].energies, commonIVData[1].data[ic], Plot.LINE);
                plot.setJustification(Plot.RIGHT);
                label += "\nR="+IJ.d2s(rFactorEtc[R_FACTOR_TYPE], 3)+
                        "\nratio2/1="+IJ.d2s(rFactorEtc[LeedRFactor.RATIO]);
                plot.addLabel(0.99, 0.05, label);
                plot.addToStack();
                if ((ic&0x1f)==0)
                    IJ.showProgress(0.8+0.2*ic/rFactorData.length);
            }
        }
        if (spotPattern.isSuperstructure() && nIntSup[0]>0 && nIntSup[1]>0) {
            for (int intOrSup=0; intOrSup<2; intOrSup++) {
                rt.incrementCounter();
                rt.addLabel(intOrSup == 0 ? "Integer" : "Superstr.");
                rt.addValue(LeedRFactor.R_FACTOR_NAMES[R_FACTOR_TYPE], sumIntSup[intOrSup][R_FACTOR_TYPE]/nIntSup[intOrSup]);
                rt.addValue("E_Range", sumIntSup[intOrSup][LeedRFactor.N_OVERLAP]);
                rt.addValue("Ratio2/1", sumIntSup[intOrSup][LeedRFactor.RATIO]/nIntSup[intOrSup]);
            }
        }
        if (allowShift) {
            rt.incrementCounter();
            rt.addLabel("E Shift");
            double shift = LeedUtils.lastElement(rFactorData)[LeedRFactor.OUTPUT_SIZE];
            rt.addValue(LeedRFactor.R_FACTOR_NAMES[R_FACTOR_TYPE], IJ.d2s(shift)+" eV");
            rt.addValue("E_Range", "");
            rt.addValue("Ratio2/1", "");
        }
        rt.setPrecision(5);
        int eRangeColumn = rt.getColumnIndex("E_Range");
        if (eRangeColumn > 0)
            rt.setDecimalPlaces(eRangeColumn, 1);
        rt.show(title);

        if (showPlot) {
            plot.show();
            plot.getImagePlus().setTitle(title);
        }
        IJ.showProgress(1.0);
    }

    /** Based on an ivData[2] input array, returns an ivData[2] array of ivData,
     *  limited to the beams common to both datasets and interpolated to the same
     *  energy step (fine enough for an R factor calculation).
     *  Note that spotPattern may be null; then beams to compare must have exactly the same indices
     *  Returns an error String on failure */
    public static Object getCommonIVData(LeedIVData[] ivData, LeedSpotPattern spotPattern) {
        boolean groupsMode = spotPattern.hasGroups() &&
                !hasEquivalentBeams(ivData[0], spotPattern) &&
                !hasEquivalentBeams(ivData[1], spotPattern);  //we got one beam per group, so let's compare groups
        /* duplicate (so we can modify the data) and remove unused data */
        ivData = new LeedIVData[] {ivData[0].duplicate(), ivData[1].duplicate()};
        for (LeedIVData ivd : ivData) {
            boolean ok = ivd.trimEnergies();
            if (!ok)
                return("Error: No valid data"+(ivd.title==null ? "" : " in "+ivd.title));
            ivd.trimData();
        }
        if (groupsMode) {
            useFirstGroupBeam(ivData[0], spotPattern);
            useFirstGroupBeam(ivData[1], spotPattern);
        }
        /* find the common beams and delete the others */
        deleteNonCommonBeams(ivData[0], ivData[1], spotPattern);
        if (ivData[0].spotIndices.length == 0)
            return("Error: Can't calculate R factor:\nData sets have no common beams");

        /* interpolate to a fine energy grid */
        double[] firstInc0 = LeedUtils.getFirstAndIncrement(ivData[0].energies);
        double[] firstInc1 = LeedUtils.getFirstAndIncrement(ivData[1].energies);
        double step0 = firstInc0[1];
        double step1 = firstInc1[1];
        if (!(step0>0 && step1>0))
            return("Error: Energies are not evenly spaced with positive steps");

        double step = Math.min(step0, step1);
        if (step > MAX_E_STEP) step = MAX_E_STEP; //we always calculate the R factor with these energy steps, or less
        if (Math.abs(step0 - step) > 1e-6)
            ivData[0] = LeedInterpolator.interpolate(ivData[0], step);
        if (Math.abs(step1 - step) > 1e-6)
            ivData[1] = LeedInterpolator.interpolate(ivData[1], step);
        return ivData;
    }

    /** Returns the R factor between two data sets that have been prepared with getCommonIVData.
     *  In case of an error, returns NaN, and displays an error message if showError=true */
    public static double getRFactor(LeedIVData ivData0, LeedIVData ivData1, double v0i, int rFactorType, boolean allowShift, boolean showError) {
        Object rFactorData = getRFactorData(ivData0, ivData1, v0i, rFactorType, allowShift, /*summaryOnly=*/true);
        if (rFactorData instanceof String) {
            if (showError) IJ.error((String)rFactorData);
            return Double.NaN;
        }
        return LeedUtils.lastElement((double[][])rFactorData)[rFactorType];
    }

    /** Returns double[][] data on the R factor for each common beam as given by LeedRFactor.getRFactor,
     *  and, as a last array, the overall result, which is the weighted average
     *  except for overlap in points, which is the sum (the total overlap).
     *  With summaryOnly=true, the result will contain only one subarray,
     *  that for the total energy range and R factor.
     *  Returns an error String on failure */
    static Object getRFactorData(LeedIVData ivData0, LeedIVData ivData1, double v0i, int rFactorType, boolean allowShift, boolean summaryOnly) {
        double step = LeedUtils.getEnergyStep(ivData0.energies);
        int eShift = (int)Math.round((ivData0.energies[0] - ivData1.energies[0])/step);

        Object rFactorData = allowShift ?
                findBestRFactorData(ivData0.data, ivData1.data, v0i, step, eShift, rFactorType, summaryOnly) :
                getRFactorData(ivData0.data, ivData1.data, v0i, step, eShift, rFactorType, summaryOnly);
        return rFactorData;
    }

    /** Returns the individual R factor data for each beam group, and the total R factors as the last array.
     *  Columns of data0, data1 must correspond the same beam if the column index is the same.
     *  The output has one array as provided by LeedRFactor for each spot present in both data sets and,
     *  as a last array, the overall result, which is the weighted average
     *  except for overlap in points, which is the sum (the total overlap).
     *  With summaryOnly=true, the result will contain only one subarray,
     *  that for the total energy range and R factor. */
    static double[][] getRFactorData(double[][] data0, double[][] data1, double v0i, double step,
            int shift, int rFactorType, boolean summaryOnly) {
        double[][] output = new double[summaryOnly ? 1 : data0.length+1][];
        double[] sums = new double[LeedRFactor.OUTPUT_SIZE];
        for (int ic=0; ic<data0.length; ic++) {
            double[] rFactorEtc = LeedRFactor.getRFactor(data0[ic], data1[ic], shift, v0i/step, null);
            rFactorEtc[LeedRFactor.N_OVERLAP] *= step;
            if (!summaryOnly)
                output[ic] = rFactorEtc; //we must have a consistent number of data sets even if overlap = 0 at some shift
            if (rFactorEtc[LeedRFactor.N_OVERLAP] > 0) {
                for (int i=0; i<LeedRFactor.OUTPUT_SIZE; i++)
                    sums[i] += (i == LeedRFactor.N_OVERLAP) ?
                            rFactorEtc[i] :
                            rFactorEtc[i]*rFactorEtc[LeedRFactor.N_OVERLAP];
            }
        }
        double totalOverlap = sums[LeedRFactor.N_OVERLAP];
        for (int i=0; i<LeedRFactor.OUTPUT_SIZE; i++)
            if (i != LeedRFactor.N_OVERLAP)
                sums[i] *= 1.0/totalOverlap;
        output[output.length-1] = sums;
        return output;
    }

    /** Finds the optimum R factor, using 3-point parabolic interpolation.
     *  Returns as a double[][] array a list with the individual R factor data
     *  for each beam group, and the total R factors.
     *  Columns of data0, data1 must correspond the same beam if the column index is the same.
     *  The output has one array as provided by LeedRFactor for each spot present in both data sets and,
     *  as a last array, the overall result, which is the weighted average 
     *  except for overlap in points, which is the total overlap.
     *  In addition, the last array has the shift added at the end.
     *  With summaryOnly=true, the result will contain only one subarray,
     *  that for the total energy range and R factor.
     *  Returns an error String if no minimum found. */
    static Object findBestRFactorData(double[][] data0, double[][] data1, double v0i, double step,
            int eShift, int rFactorType, boolean summaryOnly) {
        int maxShift = (int)Math.round(v0i/step*2);
        double[][][] rFactorData = new double[2*maxShift+1][][];
        int searchDirection = 1; //+1 or -1
        int shift = 0;
        double[][] theData = getRFactorData(data0, data1, v0i, step, shift+eShift, rFactorType, summaryOnly);
        rFactorData[shift+maxShift] = theData;
        double bestR = theData[theData.length-1][rFactorType];
        int bestShift = shift;
        if (IJ.debugMode) IJ.log("shift "+shift+": R="+IJ.d2s(bestR, 7));
        boolean first = true;
        while (true) {
            shift += searchDirection;
            if (Math.abs(shift) > maxShift)
                return "Error: No minimum of R vs. shift found within "+IJ.d2s((maxShift-1)*step);
            theData = getRFactorData(data0, data1, v0i, step, shift+eShift, rFactorType, summaryOnly);
            rFactorData[shift+maxShift] = theData;
            double rFactor = theData[theData.length-1][rFactorType];
            if (IJ.debugMode) IJ.log("shift "+shift+": R="+IJ.d2s(rFactor, 7));
            if (rFactor < bestR) {
                bestR = rFactor;
                bestShift = shift;
            } else if (first) {  //we are searching in the wrong direction or shift=0 is best. Do the other direction.
                searchDirection = -searchDirection;
                shift = 0;
            } else { //R goes up: we passed by the minimum, so we can fit a parabola
                double rLeft = LeedUtils.lastElement(rFactorData[bestShift+maxShift-1])[rFactorType];
                double rMid  = LeedUtils.lastElement(rFactorData[bestShift+maxShift])[rFactorType];
                double rRight= LeedUtils.lastElement(rFactorData[bestShift+maxShift+1])[rFactorType];
                double twob = rRight - rLeft;                   //parabola: a+b*dx+c*dx^2, has extremum at x=-b/2c
                double c = 0.5*(rRight + rLeft) - rMid;
                double dx = -0.25*twob/c;
                double[][] output = new double[summaryOnly ? 1 : data0.length+1][];
                for (int ic=0; ic<output.length; ic++) {        //write best R factor to the table
                    if (rFactorData[bestShift+maxShift][ic] == null) continue;
                    int outputSize = LeedRFactor.OUTPUT_SIZE;
                    if (ic == output.length-1) outputSize++;    //last array gets shift added at end
                    if (!summaryOnly || ic == output.length-1)
                        output[ic] = new double[outputSize];
                    for (int i=0; i<LeedRFactor.OUTPUT_SIZE; i++) {
                        if (i == LeedRFactor.N_OVERLAP) {
                            output[ic][i] = rFactorData[bestShift+maxShift][ic][LeedRFactor.N_OVERLAP];
                        } else {
                            double left = rFactorData[bestShift+maxShift-1][ic][i];
                            double mid  = rFactorData[bestShift+maxShift][ic][i];
                            double right= rFactorData[bestShift+maxShift+1][ic][i];
                            double bDx  = 0.5*dx*(right - left);    //parabola: a+b*dx+c*dx^2
                            double cDx2 = dx*dx*(0.5*(right + left) - mid);
                            output[ic][i] = mid + bDx + cDx2;
                        }
                        if (ic == output.length-1)
                            output[ic][output[ic].length-1] = (bestShift + dx)*step; //shift of parabola minimum
                    }
                }
                return output;
            }
            first = false;
        }
    }

    /** Returns whether the ivData contain symmetry-equivalent beams.
     *  ivData must have a spotIndices array for the given spotPattern. */
    static boolean hasEquivalentBeams(LeedIVData ivData, LeedSpotPattern spotPattern) {
        LeedIntegerArray groups = new LeedIntegerArray(ivData.data.length);
        for (int ic=0; ic<ivData.data.length; ic++) {
            if (ivData.data[ic] == null) continue;
            if (ivData.spotIndices[ic] < 0) continue;
            int group = spotPattern.getGroup(ivData.spotIndices[ic]);
            if (groups.indexOf(group) >= 0) return true;
            if (group != LeedSpotPattern.UNKNOWN_GROUP) groups.add(group);
        }
        return false;
    }

    /** Replaces the spotIndices with those of the first beam of each group
     *  of symmetry-equivalent beams, to ensure that symmetry-equivalent beams
     *  of both files have the same name and spot index.
     *  Does not care about the spotNames (this will be done later) */
    static void useFirstGroupBeam(LeedIVData ivData, LeedSpotPattern spotPattern) {
        for (int ic=0; ic<ivData.data.length; ic++) {
            if (ivData.data[ic] == null) continue;
            int group = spotPattern.getGroup(ivData.spotIndices[ic]);
            int[] groupIndices = spotPattern.getAllSpotsForGroup(group);
            ivData.spotIndices[ic] = groupIndices[0];
        }
    }

    /** Averages the symmetry-equivalent beams.
     *  One beam per group of symmetry-equivalent beams will be replaced
     *  by the average, the others will be deleted from  ivData.
     *  ivData must have a spotIndices array for the given spotPattern. */
    void averageEquivalentBeams(LeedIVData ivData, LeedSpotPattern spotPattern) {
        double energyStep = LeedUtils.getEnergyStep(ivData.energies);
        boolean deletion = false;
        for (int ic=0; ic<ivData.data.length; ic++) {
            if (ivData.data[ic] == null) {
                deletion = true;
                continue;
            }
            int group = spotPattern.getGroup(ivData.spotIndices[ic]);
            int[] gSpots = spotPattern.getAllSpotsForGroup(group);
            if (gSpots.length <=1) continue;                    //skip if no equivalent spots in spotPattern
            ArrayList<double[]> dataForAvg = new ArrayList<double[]>(12);
            for (int iSpot : gSpots) {
                int iData = LeedUtils.arrayIndexOf(ivData.spotIndices, iSpot);
                if (iData >= 0) {
                    dataForAvg.add(ivData.data[iData]);
                    if (iData != ic)
                        ivData.data[iData] = null;              //data merged in will become obsolete
                }
            }
            if (dataForAvg.size() <= 1) continue;               //skip if no equivalent spots in ivData
            ivData.data[ic] = calculateAverage(dataForAvg);
            deletion = true;
        }
        if (deletion)
            ivData.trimData();
    }

    /** Deletes the beams not common to both ivData sets, sorts the sets and
     *  and, if we have a spotPattern, sets the spotNames in both ivData sets to the canonical form.
     *  Both ivData sets must have a spotIndices array. */
    static void deleteNonCommonBeams(LeedIVData ivData0, LeedIVData ivData1, LeedSpotPattern spotPattern) {
        boolean[] inUse1 = new boolean[ivData1.data.length];
        for (int ic=0; ic<ivData0.data.length; ic++) {
            if (ivData0.data[ic] == null) {
                ivData0.spotIndices[ic] = -1;
                continue;
            }
            int spotIndex = ivData0.spotIndices[ic];
            if (spotIndex < 0) continue;
            int arrayIndex1 = LeedUtils.arrayIndexOf(ivData1.spotIndices, spotIndex);
            if (arrayIndex1 >= 0 && ivData1.data[arrayIndex1] != null)
                inUse1[arrayIndex1] = true;
            else
                ivData0.spotIndices[ic] = -1;   //no counterpart in ivData1, mark for deletion
        }
        for (int ic1=0; ic1<inUse1.length; ic1++)
            if (!inUse1[ic1])
                ivData1.spotIndices[ic1] = -1;  //mark unused ivData1 columns for deletion
        ivData0.deleteIfSpotIndex(-1);
        ivData1.deleteIfSpotIndex(-1);
        ivData0.sort();
        ivData1.sort();
        ivData0.setCanonicalNames(spotPattern);
        ivData1.setCanonicalNames(spotPattern);
        if (!Arrays.equals(ivData0.spotIndices, ivData1.spotIndices))
            throw new RuntimeException("Internal error, beam sets should have the same beams now, but that's not the case");
    }

    /** Normalizes all beams such that the maximum intensity in the overlap region
     *  of ivData[0] and ivData[1] becomes 1.0. Both ivData sets must have the
     *  same number of beams (number of elements in 'data') with 1:1 correspondence
     *  of the beams, and the same energy step. */
    static void normalizeBeams(LeedIVData[] ivData) {
        double step = LeedUtils.getEnergyStep(ivData[0].energies);
        int eShift =                            //positive if 0 starts later, first energies of 1 are unused
                (int)Math.round((ivData[0].energies[0] - ivData[1].energies[0])/step);
        for (int ic=0; ic<ivData[0].data.length; ic++) {
            int[] startEnd0 = LeedUtils.getStartEndNum(ivData[0].data[ic]);
            int[] startEnd1 = LeedUtils.getStartEndNum(ivData[1].data[ic]);
            if (startEnd0[0]<0 || startEnd1[0]<0)
                continue;                       //empty data, should not happen since we have used trimData()
            int start0 = Math.max(startEnd0[0], startEnd1[0]-eShift);
            int end0 = Math.min(startEnd0[1], startEnd1[1]-eShift);            
            for (int id=0; id<ivData.length; id++) {
                int start = id == 0 ? start0 : start0 + eShift;
                int end   = id == 0 ? end0   : end0   + eShift;
                double max = LeedUtils.getArrayMax(ivData[id].data[ic], start, end);
                if (max > 0)
                    for (int i=0; i<ivData[id].data[ic].length; i++)
                        ivData[id].data[ic][i] *= 1.0/max;
            }
        }
    }

    /** Calculates the average over equivalent beams, without fade in/fade out.
     *  Input data that would lead to a reduction of the usable range by more than 50% are ignored.
     *  Returns null if no usable input data are present.
     *  This type of averaging is used, e.g., for calculated I(V) curves where the ranges differ slightly
     *  if slightly off-normal incidence is used to simulate a converging beam. */
    public static double[] calculateAverage(ArrayList<double[]> inData) {
        inData = new ArrayList<double[]>(inData); //copy the original list since we will remove the elements from the list
        int length = -1;
        /* eliminate unusable input arrays */
        for (int i=inData.size()-1; i>=0; i--) {
            if (inData.get(i) == null || LeedUtils.countNonNaN(inData.get(i)) < 3) {
                inData.remove(i);
                continue;
            }
            if (length < 0)
                length = inData.get(i).length;
            else if (inData.get(i).length != length)
                throw new IllegalArgumentException("Arrays for abveraging have different lengths: "+
                        inData.get(i).length+" vs "+length);
        }
        if (inData.isEmpty())
            return null;
        if (inData.size() == 1)
            return inData.get(0);
        /* find the array with the best (longest) total data range */
        int bestI = -1;
        int bestNonNaN = Integer.MIN_VALUE;
        for (int i=0; i<inData.size(); i++) {
            int nNonNaN = LeedUtils.countNonNaN(inData.get(i));
            if (nNonNaN > bestNonNaN) {
                bestNonNaN = nNonNaN;
                bestI = i;
            }
        }
        double[] out = (double[])(inData.get(bestI).clone());
        inData.remove(bestI);
        int nAveraged = 1;
        while (!inData.isEmpty()) {
            double[] bestOverlapData = null;
            int bestOverlapPoints = Integer.MIN_VALUE;
            for (int i=inData.size()-1; i>=0; i--) {
                int[] overlap = LeedRFactor.getOverlap(inData.get(i), out);    //(returns number of points as the last element)
                if (overlap == null || 2*overlap[overlap.length-1] < bestNonNaN) {
                    inData.remove(i);               //overlap too short, ignore data column
                    continue;
                }
                int overlapPoints = overlap[overlap.length-1];
                if (overlapPoints > bestOverlapPoints) {
                    bestOverlapPoints = overlapPoints;
                    bestOverlapData = inData.get(i);
                }
            }
            if (bestOverlapData == null) break;
            for (int i=0; i<length; i++)
                out[i] += bestOverlapData[i];
            inData.remove(bestOverlapData);
            nAveraged++;
        }
        for (int j=0; j<length; j++)
            out[j] *= 1.0/nAveraged;
        return out;
    }
}
