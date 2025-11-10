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
 *  This means that the collection of ViPErLEED ImageJ plugins is free software:
 *  you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation,
 *  either version 3 of the License, or (at your option) any later version.
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  The collection of ViPErLEED ImageJ plugins is distributed
 *  in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  The authors may decide later to put part of the auxiliary code in this work into
 *  the public domain, to allow incorporation into ImageJ if desired.
 *  (ImageJ is in the public domain.)
 *
 *  The help text is licensed under the Creative Commons Attribution 4.0 (CC BY 4.0) license.
 *
 *  When using and/or modifying this program for scientific work, please cite
 *  the paper describing it:
 *  M. Schmid, F. Kraushofer, A. M. Imre, T. Kißlinger, L. Hammer, U. Diebold, and M. Riva,
 *  ViPErLEED package II: Spot tracking, extraction and processing of I(V) curves,
 *  Phys. Rev. Research 7, 013006 (2025).
 *  @author Michael Schmid, IAP/TU Wien, 2019-2025
 */


public class LEED_R_Factor_Between_Datasets implements PlugIn {
    static final String PLUGIN_NAME = "R Factor Between Data Sets";
    static final double MAX_E_STEP = 0.5;       //energy step for R factor comparison no larger than this
    /** Dialog parameters, should be remembered until ImageJ is closed */
    static String[] filePaths = new String[2];
    static String spotPatternPath = LeedParams.getString(LeedParams.PATTERNFILE);
    static double v0i = LeedParams.get(LeedParams.V0I);
    static int rFactorType = (int)LeedParams.get(LeedParams.RFACTORTYPE);
    static boolean averageEquivalent = true;
    static boolean allowShift = true;
    static final int NO_PLOT = 0, PLOT_I = 1, PLOT_I_AND_DY = 2; //plot types
    static final String[] PLOT_TYPES = new String[] {
            "nothing", "Intensities", "Intensities & Y Difference"};
    static int plotType = 0;
    /** Remember results for call in a macro */
    static double resultR = Double.NaN, resultOverlap = Double.NaN, resultRatio = Double.NaN;
    /** Shorthand notation for LeedRFactor conatants */
    static final int NUMERATOR = LeedRFactor.NUMERATOR;
    static final int DENOMINATOR = LeedRFactor.DENOMINATOR;
    private static final String HELP_STRING =
            "<html>"+
            "<h1>"+PLUGIN_NAME+"</h1>"+
            "<p>This plugin compares two <i>I</i>(<i>V</i>) data files (.csv format).</p>"+
            "<p>It provides a list with the R factor for each beam as well as the total <i>R</i> factor.<br />"+
            "In addition to the <i>R</i> factor, the energy span (where both files contain valid data) "+
            "and ratio between the average intensity per beam of the two files are reported. "+
            "The ratio is normalized such that the overall ratio is 1.0.<br />"+
            "In the case of a superstructure (i.e., if at least one non-integer beam is present), "+
            "separate summaries are given for the integer and superstructure beams. "+
            "(If the intensity ratios for superstructure beams strongly differ from integer beams, "+
            "this may be an indication for partial coverage by the superstructure, or, "+
            "when comparing experiment and calculation, for a wrong structure model.)</p>"+
            "<p>Note that the <i>R</i> factors are calculated in a rather simple way, ignoring the end points "+
            "(where the derivatives are poorly defined). Results may slightly differ from those of other <i>R</i> factor programs.</p>"+
            "<p><b>V0i</b> is the absolute value of the imaginary part of the inner potential, required for "+
            "calculating the R factors (typically 3.5 to 6 eV).</p>"+
            "<p><b>R factor type</b>: Currently supported are Pendry's R factor <i>R</i><sub>P</sub> and "+
            "the 'smooth' R factor <i>R</i><sub>S</sub> <a href='#Rsmooth'>[2]</a>.</p>"+
            "<p>The option <b>Average symmetry-equivalent beams</b> is mainly meant for averaging calculated beams. "+
            "The experimental <i>I</i>(<i>V</i>) curves should be averaged with the <i>I</i>(<i>V</i>) Curve Editor. "+
            "In contrast to the <i>I</i>(<i>V</i>) Curve Editor, the averaging algorithm used here does not use smooth fade-in/fade-out "+
            "if the curves start/end at different energies; "+
            "the average is rather limited to the common energy range of the input data. "+
            "If a curve is defined over an energy range that is too short (if would reduce the energy range to less than half) "+
            "it is omitted in the average.<br />"+
            "If averaging is off, comparison is done on a beam-per-beam basis. "+
            "If both of the input files contain only one beam per group of symmetry-equivalent beams, "+
            "comparison will be done for beam groups. In this case, beams will be also compared if the group is the same, "+
            "but the beam indices are different.</p>"+
            "<p><b>Allow shifting</b> searches for the best offset of the energy axes between the two input files "+
            "(the offset is the same for all beams). "+
            "When comparing calculated and experimental <i>I</i>(<i>V</i>) curves, this essentially corresponds to a variation of "+
            "the real part of the inner potential (V_0r) or a variation of the filament work function.</p>"+
            "<p><b>Plot:</b> When 'Intensities' is selected, creates a plot stack, with one plot for each beam, "+
            "showing the respective <i>I</i>(<i>V</i>) curves of the two input files. "+
            "The curves of the first plot are red, the second black. "+
            "With 'Intensities &amp; Y Difference', there is an additional blue curve, with the difference between "+
            "the two Y functions of the R factor. Since the integral of the square of this difference determines the "+
            "R factor, its highest deviations from the zero line (center of the plot) show which sections of the curves "+
            "are cause a high (='bad') R factor.</p>"+
            "<p>If the groups of symmetry-equivalent beams cannot be read from both input files "+
            "(or if the group numbers in these two files are inconsistent) a <b>spot pattern file</b> is required. "+
            "This is a .csv file containing a list of the beams with additional information on them, "+
            "such as the group number, which defines groups of symmetry-equivalent beams. "+
            "A spot pattern file can be also supplied to override the symmetry assumed in the input data files "+
            "(i.e., the group numbers in square brackets).<br />"+
            "A spot pattern file can be created with the ViPErLEED gui (select 'export').</p>"+
            "<p><b>For use with the ImageJ Macro Language</b>, use the Macro Recorder to see how to call the plugin. "+
            "You can then access the overall result not only from the table created but also via the macro command<br>"+
            "<tt>&nbsp;&nbsp;result=call('LEED_R_Factor_Between_Datasets.getResult', 'R');</tt><br>"+
            "The result is a String that can be converted to a number with <tt>parseFloat(result)</tt>. "+
            "The argument 'R' may be replaced with 'overlap' to access the total energy range or 'ratio' "+
            "for the average intensity ratio between file 2 (black) and file 1 (red in the plot, if a plot is created).</p>"+
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
            "<p><a name='paper'>[1]</a> M. Schmid, F. Kraushofer, A. M. Imre, T. Ki&szlig;linger, L. Hammer, U. Diebold, and M. Riva<br>"+
            "<i>ViPErLEED package II: Spot tracking, extraction and processing of <i>I</i>(<i>V</i>) curves</i><br>"+
            "<a href='https://doi.org/10.1103/PhysRevResearch.7.013006'>Phys. Rev. Research 7, 013006 (2025)</a>.</p>"+
            "<p><a name='Rsmooth'>[2]</a> A. M. Imre, L. Hammer, U. Diebold, M. Riva, M. Schmid<br>"+
            "<i>An improved reliability factor for quantitative low-energy electron diffraction</i><br>"+
            "<a href='https://arxiv.org/abs/2511.05448/'>arXiv:2511.05448</a>.</p>"+
            "</html>";

    public void run(String arg) {
        if (IJ.versionLessThan("1.52u")) return;

        // I N I T I A L   D I A L O G
        if (spotPatternPath==null)
            spotPatternPath = LeedParams.getString(LeedParams.PATTERNFILE);
        GenericDialog gd = new GenericDialog("R Factor Between Sets of IV curves");
        gd.addFileField("File_1 (red)", filePaths[0], 50);
        gd.addFileField("File_2 (black)", filePaths[1], 50);
        gd.addFileField("Spot pattern file", spotPatternPath, 50);
        gd.addMessage(".csv files required; you can drag and drop files into the fields");
        gd.addNumericField("V0i", v0i, 2, 6, "eV");
        gd.addChoice("R factor type", LeedRFactor.R_FACTOR_NAMES, LeedRFactor.R_FACTOR_NAMES[rFactorType]);
        gd.addCheckbox("Average symmetry-equivalent beams", averageEquivalent);
        gd.addCheckbox("Allow shifting (V0r variation)", allowShift);
        gd.addChoice("Plot:", PLOT_TYPES, PLOT_TYPES[plotType]);

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
        rFactorType = gd.getNextChoiceIndex();

        averageEquivalent = gd.getNextBoolean();
        allowShift = gd.getNextBoolean();
        plotType = gd.getNextChoiceIndex();

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
        double step = LeedUtils.getEnergyStep(commonIVData[0].energies);
        Object rFactorDataOrError = getRFactorData(commonIVData[0], commonIVData[1], rFactorType, v0i, step, allowShift, /*summaryOnly=*/false);
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
        if (spotPattern.isSuperstructure()) {                           //create separate statistics for integer & superstructure
            for (int ic=0; ic<commonIVData[0].spotIndices.length; ic++) {
                if (rFactorData[ic][LeedRFactor.N_OVERLAP] == 0) continue;
                int intOrSup = spotPattern.isSuperstructure(commonIVData[0].spotIndices[ic]) ?
                        1 : 0;
                for (int i=0; i<LeedRFactor.OUTPUT_SIZE; i++) {
                    if (i==LeedRFactor.N_OVERLAP)
                        sumIntSup[intOrSup][i] += rFactorData[ic][i];
                    else
                        sumIntSup[intOrSup][i] += rFactorData[ic][i]*rFactorData[ic][LeedRFactor.N_OVERLAP]; //weighted sum
                }
            }
        }
        double[] energies1 = commonIVData[1].energies;
        if (allowShift && plotType != NO_PLOT) {
            double shift = LeedUtils.lastElement(rFactorData)[LeedRFactor.OUTPUT_SIZE];
            energies1 = new double[commonIVData[1].energies.length];
            for (int i=0; i<energies1.length; i++)
                energies1[i] = commonIVData[1].energies[i] - shift;
        }
        int eShift = (int)Math.round((commonIVData[0].energies[0] - commonIVData[1].energies[0])/step); //>0 if 0 starts later

        // O U T P U T
        String title = ivData[0].title+"_vs_"+ivData[1].title;
        ResultsTable rt = new ResultsTable();
        Plot plot = plotType == NO_PLOT ? null : new Plot(title, "E", "intensity", Plot.DEFAULT_FLAGS & (~Plot.Y_NUMBERS));
        if (plotType != NO_PLOT) {
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
            rt.addValue(LeedRFactor.R_FACTOR_NAMES[rFactorType], rFactorEtc[NUMERATOR]/rFactorEtc[DENOMINATOR]);
            rt.addValue("E_Range", rFactorEtc[LeedRFactor.N_OVERLAP]);
            rt.addValue("Ratio2/1", rFactorEtc[LeedRFactor.RATIO]);
            //IJ.log(label+", denom/dE,"+IJ.d2s(rFactorEtc[DENOMINATOR]/rFactorEtc[LeedRFactor.N_OVERLAP]));//shows that denominator is roughly proportional to E range
            if (!isSpotData) {                                          //remember data for macro call
                resultR = rFactorEtc[NUMERATOR]/rFactorEtc[DENOMINATOR];
                resultOverlap = rFactorEtc[LeedRFactor.N_OVERLAP];
                resultRatio = rFactorEtc[LeedRFactor.RATIO];
            }

            if (plotType != NO_PLOT && isSpotData) {
                plot.setColor(Color.RED);
                plot.addPoints(commonIVData[0].energies, commonIVData[0].data[ic], Plot.LINE);
                plot.setColor(Color.BLACK);
                plot.addPoints(energies1, commonIVData[1].data[ic], Plot.LINE);
                plot.setJustification(Plot.RIGHT);
                label += "\nR="+IJ.d2s(rFactorEtc[NUMERATOR]/rFactorEtc[DENOMINATOR], 3)+
                        "\nratio2/1="+IJ.d2s(rFactorEtc[LeedRFactor.RATIO]);
                plot.addLabel(0.99, 0.05, label);
                if (plotType == PLOT_I_AND_DY) {
                    double[] data1forY = LeedInterpolator.getInterpolated(energies1, commonIVData[1].energies, commonIVData[1].data[ic]);
                    double[] y0 = LeedRFactor.getYcurve(commonIVData[0].data[ic], rFactorType, v0i/step, null);
                    double[] y1 = LeedRFactor.getYcurve(data1forY, rFactorType, v0i/step, null);
                    
                    for (int i0=0; i0<y0.length; i0++) {        //y0 becomes difference of y functions
                        int i1 = i0 + eShift;
                        if (i1 >= 0 && i1 < y1.length)
                            y0[i0] = 0.5 + 0.25*(y1[i1] - y0[i0]);
                        else
                         y0[i0] = Double.NaN;
                    }
                    plot.setColor(Color.BLUE);
                    plot.addPoints(commonIVData[0].energies, y0, Plot.LINE);
                }
                plot.addToStack();
                if ((ic&0x1f)==0)
                    IJ.showProgress(0.8+0.2*ic/rFactorData.length);
            }
        }
        if (spotPattern.isSuperstructure() && sumIntSup[0][LeedRFactor.N_OVERLAP]>0 && sumIntSup[1][LeedRFactor.N_OVERLAP]>0) {
            for (int intOrSup=0; intOrSup<2; intOrSup++) {
                rt.incrementCounter();
                rt.addLabel(intOrSup == 0 ? "Integer" : "Superstr.");
                rt.addValue(LeedRFactor.R_FACTOR_NAMES[rFactorType],
                        LeedRFactor.toRFactor(sumIntSup[intOrSup]));
                rt.addValue("E_Range", sumIntSup[intOrSup][LeedRFactor.N_OVERLAP]);
                rt.addValue("Ratio2/1", sumIntSup[intOrSup][LeedRFactor.RATIO]/sumIntSup[intOrSup][LeedRFactor.N_OVERLAP]);
            }
        }
        if (allowShift) {
            rt.incrementCounter();
            rt.addLabel("E Shift");
            double shift = LeedUtils.lastElement(rFactorData)[LeedRFactor.OUTPUT_SIZE];
            rt.addValue(LeedRFactor.R_FACTOR_NAMES[rFactorType], IJ.d2s(shift)+" eV");
            rt.addValue("E_Range", "");
            rt.addValue("Ratio2/1", "");
        }
        rt.setPrecision(5);
        int eRangeColumn = rt.getColumnIndex("E_Range");
        if (eRangeColumn > 0)
            rt.setDecimalPlaces(eRangeColumn, 1);
        rt.show(title);

        if (plotType != NO_PLOT) {
            plot.show();
            plot.getImagePlus().setTitle(title);
        }
        IJ.showProgress(1.0);
    }


    /** Returns the results of the last invocation as a String.
     *  The argument (case-insensitive) determines which value to return:
     *  'overlap' for the overlap, 'ratio' for ratio between data in 2 and 1. Otherwise returns the R factor.  */
    public static String getResult(String what) {
        what = what.toLowerCase();
        double result;
        if (what.contains("over")) result = resultOverlap;
        else if (what.contains("rat")) result = resultRatio;
        else result = resultR;
        return Double.toString(result);
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
    public static double getRFactor(LeedIVData ivData0, LeedIVData ivData1,  int rFactorType, double v0i, double step, boolean allowShift, boolean showError) {
        Object rFactorData = getRFactorData(ivData0, ivData1, rFactorType, v0i, step, allowShift, /*summaryOnly=*/true);
        if (rFactorData instanceof String) {
            if (showError) IJ.error((String)rFactorData);
            return Double.NaN;
        }
        return LeedRFactor.toRFactor(LeedUtils.lastElement((double[][])rFactorData));
    }

    /** Returns double[][] data on the R factor for each common beam as given by LeedRFactor.getRFactorData,
     *  and, as a last array, the overall result, which is the weighted average
     *  except for overlap in points, which is the sum (the total overlap).
     *  With summaryOnly=true, the result will contain only one subarray,
     *  that for the total energy range and R factor.
     *  Returns an error String on failure */
    static Object getRFactorData(LeedIVData ivData0, LeedIVData ivData1, int rFactorType, double v0i, double step, boolean allowShift, boolean summaryOnly) {
        int eShift = (int)Math.round((ivData0.energies[0] - ivData1.energies[0])/step);

        Object rFactorData = allowShift ?
                findBestRFactorData(ivData0.data, ivData1.data, rFactorType, v0i, step, eShift, summaryOnly) :
                getRFactorData(ivData0.data, ivData1.data, rFactorType, v0i, step, eShift, summaryOnly);
        return rFactorData;
    }

    /** Returns the individual R factor data for each beam group, and the total R factors as the last array.
     *  Columns of data0, data1 must correspond the same beam if the column index is the same.
     *  The output has one array as provided by LeedRFactor for each spot present in both data sets and,
     *  as a last array, the overall result, which is the weighted average
     *  except for overlap in points, which is the sum (the total overlap).
     *  With summaryOnly=true, the result will contain only one subarray,
     *  that for the total energy range and R factor. */
    static double[][] getRFactorData(double[][] data0, double[][] data1, int rFactorType, double v0i, double step,
            int shift, boolean summaryOnly) {
        double[][] output = new double[summaryOnly ? 1 : data0.length+1][];
        double[] sums = new double[LeedRFactor.OUTPUT_SIZE];
        for (int ic=0; ic<data0.length; ic++) {
            double[] rFactorEtc = LeedRFactor.getRFactorData(data0[ic], data1[ic], shift, rFactorType, v0i/step, null);
            rFactorEtc[LeedRFactor.N_OVERLAP] *= step;
            if (!summaryOnly)
                output[ic] = rFactorEtc; //we must have a consistent number of data sets even if overlap = 0 at some shift
            if (rFactorEtc[LeedRFactor.N_OVERLAP] > 0) {
                for (int i=0; i<LeedRFactor.OUTPUT_SIZE; i++)
                    sums[i] += (i > LeedRFactor.N_OVERLAP) ? //for average I and ratio, weighted sum
                            rFactorEtc[i]*rFactorEtc[LeedRFactor.N_OVERLAP] :
                            rFactorEtc[i];
            }
        }
        double totalOverlap = sums[LeedRFactor.N_OVERLAP];
        for (int i=0; i<LeedRFactor.OUTPUT_SIZE; i++)
            if (i > LeedRFactor.N_OVERLAP)
                sums[i] *= 1.0/totalOverlap;
        output[output.length-1] = sums;
        return output;
    }

    /** Finds the optimum R factor, using 3-point parabolic interpolation between energy shifts.
     *  On success, returns as a double[][] array a list with the individual R factor data
     *  for each beam group, and the total R factors.
     *  Columns of data0, data1 must correspond the same beam if the column index is the same.
     *  The output has one array as provided by LeedRFactor for each spot present in both data sets and,
     *  as a last array, the overall result, which is the weighted average 
     *  except for overlap in points, which is the total overlap.
     *  In addition, the last array has the shift added at the end.
     *  With summaryOnly=true, the result will contain only one subarray,
     *  that for the total energy range and R factor.
     *  Returns an error String if no minimum found. */
    static Object findBestRFactorData(double[][] data0, double[][] data1, int rFactorType, double v0i, double step,
            int eShift, boolean summaryOnly) {
        int maxShift = (int)Math.round(v0i/step*2);
        double[][][] rFactorData = new double[2*maxShift+1][][];
        int searchDirection = 1; //+1 or -1
        int shift = 0;
        double[][] theData = getRFactorData(data0, data1, rFactorType, v0i, step, shift+eShift, summaryOnly);
        rFactorData[shift+maxShift] = theData;
        double bestR = LeedRFactor.toRFactor(theData[theData.length-1]);
        int bestShift = shift;
        if (IJ.debugMode) IJ.log("shift "+shift+": R="+IJ.d2s(bestR, 7));
        boolean first = true;
        while (true) {
            shift += searchDirection;
            if (Math.abs(shift) > maxShift)
                return "Error: No minimum of R vs. shift found within "+IJ.d2s((maxShift-1)*step);
            theData = getRFactorData(data0, data1, rFactorType, v0i, step, shift+eShift, summaryOnly);
            rFactorData[shift+maxShift] = theData;
            double rFactor = LeedRFactor.toRFactor(theData[theData.length-1]);
            if (IJ.debugMode) IJ.log("shift "+shift+": R="+IJ.d2s(rFactor, 7));
            if (rFactor < bestR) {
                bestR = rFactor;
                bestShift = shift;
            } else if (first) {  //we are searching in the wrong direction or shift=0 is best. Do the other direction.
                searchDirection = -searchDirection;
                shift = 0;
            } else { //R goes up: we passed by the minimum, so we can fit a parabola
                double rLeft = LeedRFactor.toRFactor(LeedUtils.lastElement(rFactorData[bestShift+maxShift-1]));
                double rMid  = LeedRFactor.toRFactor(LeedUtils.lastElement(rFactorData[bestShift+maxShift]));
                double rRight= LeedRFactor.toRFactor(LeedUtils.lastElement(rFactorData[bestShift+maxShift+1]));
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
     *  by the average, the others will be deleted from ivData.
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
            //IJ.log(ivData.title+": group "+group+" n="+gSpots.length);
            if (gSpots.length <=1) continue;                    //skip if no equivalent spots in spotPattern
            ArrayList<double[]> dataForAvg = new ArrayList<double[]>(12);
            for (int iSpot : gSpots) {
                int iData = LeedUtils.arrayIndexOf(ivData.spotIndices, iSpot);
                if (iData >= 0) {
                    dataForAvg.add(ivData.data[iData]);
                    //IJ.log(ivData.title+": avg for curve="+ic+" includes curve="+iData+" name="+spotPattern.getName(ivData.spotIndices[iData])+" 1st="+ivData.data[iData][0]);

                    if (iData != ic)
                        ivData.data[iData] = null;              //data merged in will become obsolete
                }
            }
            if (dataForAvg.size() <= 1) continue;               //skip if no equivalent spots in ivData
            ivData.data[ic] = calculateAverage(dataForAvg);
            //IJ.log("averaging "+dataForAvg.size()+" curves; avg 1st="+ivData.data[ic][0]);
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
