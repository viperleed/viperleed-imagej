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
import java.io.FileReader;
import java.io.BufferedReader;

/** This ImageJ plugin provides statistics on the R factors between symmetry-equivalent beams
 *  and other statistics concerning the data quality, as a table or plot. */

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


public class LEED_Data_Quality_Statistics implements PlugIn, DialogListener {

    public static final int PLOT=0, LIST=1;

	private static final String PLUGIN_NAME   = "Leed Data Quality Statistics";
    private static final int    R_FACTOR_TYPE = LeedRFactor.R_PENDRY; //R factor type for plot
    private static final String[] TASKS = new String[] {
            "Plot","Data table"};
    private static final String HELP_STRING =
            "<html>"+
            "<h1>"+PLUGIN_NAME+"</h1>"+
            "<p>This ImageJ plugin determines R Factors between I(V) curves of symmetry-equivalent beams "+
            "and other statistics useful to assess the quality of experimental I(V) data.</p>"+
            "<dl>"+
            "<dt>I(V) Input Table/File:</dt><dd>File with the raw (unsmoothed) I(V) curves (csv Format).</dd>"+
            "<dt>Pattern File:</dt><dd>Usually, the <em>spot pattern</em> defining the symmetry-equivalent beams will be determined from the "+
            "heading of the I(V) curves (based on the group numbers in square brackets). "+
            "You can also select a spot pattern file of your choice. "+
            "See the Help of the ViPErLEED Spot Tracker for more information on spot pattern files.</dd>"+
            "<dt>What to do:</dt><dd>This plugin can ..."+
            "<ul><li> Create plot.<br />"+
            "This plot includes a scatter plot of R factors between equivalent beams vs average intensity "+
            "of the data section. There is also a line with the R factor against the total energy range "+
            "of all curve section pairs with an R factor better than the given number. "+
            "(If there are many symmetry-equivalent beams, this is far more than "+
            "the energy range available up to that R factor because there are many pairs for one group of beams.) "+
            "Furthermore, if there are beams with negative intensities, the absolute values of the most negative intensities "+
            "of these beams are plotted (with the energy range in keV of the negative intensities on the y axis).</li>"+
            "<li> Create a table of the mutual R factors and additional information (similar to the plot).</li></ul></dd>"+
            "<dt>V0i:</dt><dd>Imaginary part of the inner potential, needed to calculate Pendry's R factor "+
            "(strictly speaking, the absolute value; usually about 4&ndash;5&nbsp;eV).</dd>"+
            "<dt>Split energy range:</dt><dd>The I(V) curves are split into energy ranges (sections) of approximately this length. "+
            "Thus, energy regions with substantially different intensities will be plotted or listed separately. "+
            "Leave this field empty or set it to 0 to avoid dividing the I(V) curves into sections.</dd>"+
            "<dt>Energy range to smooth:</dt><dd>Determines the smoothing applied before calculating the R factor. "+
            "Typically 0.8&times;V0i to 1.3&times;V0i. Stronger smoothing leads to better noise suppression, and, thus, lower R factors. "+
            "Thus, when comparing different data sets, V0i and the smoothing parameter should be the same for all data sets.</dd>"+
            "</dl>"+
            "<h2>License</h2>"+
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
            "<i>ViPErLEED package II: Spot tracking, extraction and processing of I(V) curves</i>, "+
            "Phys. Rev. Research, 2024.</p>"+
            "</html>";
    LeedSpotPattern spotPattern;
    //dialog options
    private static String prevFilePath;
    private static String prevPatternPath;
    private static double v0i = Double.NaN;
    private static double eVsmooth = Double.NaN;
    private static double eVseparate = Double.NaN;
    private static int task = PLOT;

    GenericDialog gd;
    boolean inputHasGroups;   //whether the table of file has spot groups
    String patternLabelText;  //for spot pattern albel, if there is no user-supplied file
    Frame[] resultsWindows;   //the windows with results tables currently open


    public void run(String arg) {
        if (IJ.versionLessThan("1.52s")) return;

        // P R E P A R E
        LeedParams.initialize();
        if (Double.isNaN(v0i))
            v0i = LeedParams.get(LeedParams.V0I);
        if (Double.isNaN(eVsmooth))
            eVsmooth = LeedParams.get(LeedParams.SMOOTHEV);
        if (Double.isNaN(eVseparate))
            eVseparate = LeedParams.get(LeedParams.SPLITEV);

        // C R E A T E   D I A L O G
        gd = new GenericDialog(PLUGIN_NAME);
        gd.addFileField("I(V) input file", prevFilePath, 60);
        gd.addFileField("Spot pattern file*", prevPatternPath, 60);
        gd.setInsets(5, 0, 5); //top left bottom
        gd.addChoice("What to do", TASKS, TASKS[task]);
        gd.addNumericField("V0i", v0i, 1, 6, "eV");
        gd.addNumericField("Energy range to smooth", eVsmooth, 1, 6, "eV **");
        if (eVseparate > 1e6) eVseparate = 0;
        gd.addNumericField("Split energy range into", eVseparate, 0, 6, "eV sections ***");
        gd.addMessage("*   Spot pattern file may be omitted if groups of equivalent "+
                "beams are given in the I(V) input file [in square brackets].\n"+
                "**  Energy range of moving average with same noise suppression, "+
                "typically 0.6 V0i (very low noise) to 1.3 V0i (noisy data).\n"+
				"*** Enter 0 or leave empty to always keep the full range "+
                "per pair in one section.");

        gd.addDialogListener(this);
        gd.addHelp(HELP_STRING);
        Button okButton = gd.getButtons()[0];
        okButton.setEnabled(dialogItemChanged(gd, null));
        // S H O W   D I A L O G
        gd.showDialog();
        if (gd.wasCanceled()) return;

        String inputFilePath = gd.getNextString();
        String patternFilePath = gd.getNextString();
        int task = gd.getNextChoiceIndex();
        v0i = Math.abs(gd.getNextNumber());
        eVsmooth = gd.getNextNumber();
        if (Double.isNaN(eVsmooth)) eVsmooth = 0;
        eVseparate = gd.getNextNumber();
        if (!(eVseparate > 3)) eVseparate = 1e7;

        LeedIVData ivData = LeedIVData.fromFile(inputFilePath, /*spotPattern=*/null, LeedIVData.E_ASCENDING|LeedIVData.ZERO_IS_NAN);
        if (ivData == null) return;
        prevFilePath = inputFilePath;
        LeedSpotPattern spotPattern = null;
        if (patternFilePath.trim().length() > 0) {
            spotPattern = new LeedSpotPattern(patternFilePath);
            if (spotPattern.size() < 1) {
                IJ.error(PLUGIN_NAME, "Error: Invalid spot pattern file\n"+spotPattern.getErrorMessage());
                spotPattern = null;
            } else {
                ivData.setSpotPattern(spotPattern, /*showWarning=*/true);
                if (ivData.spotIndices != null)
                    prevPatternPath = patternFilePath;
            }
        } else {
            spotPattern = new LeedSpotPattern(ivData.spotNames, /*groupsRequired=*/true);
            if (spotPattern.size() > 1)
                ivData.setSpotPattern(spotPattern, /*showWarning=*/true);
        }
        while (ivData.spotIndices == null) {
            String defaultPatternPath = LeedParams.getString(LeedParams.PATTERNFILE);
            String savedDir = OpenDialog.getDefaultDirectory();
            if (defaultPatternPath != null && defaultPatternPath.length() > 0)
                OpenDialog.setDefaultDirectory((new File(defaultPatternPath)).getParent());
            OpenDialog od = new OpenDialog("Select Spot Pattern .csv File");
            String spotPatternPath = od.getPath();
            if (savedDir != null)
                OpenDialog.setDefaultDirectory(savedDir);
            if (spotPatternPath == null)
                return;
            spotPattern = new LeedSpotPattern(spotPatternPath);
            if (spotPattern.size() > 0)
                ivData.setSpotPattern(spotPattern, /*showWarning=*/true);
            else
                IJ.error(PLUGIN_NAME, "Error: Invalid spot pattern file\n"+spotPattern.getErrorMessage());
        }
        analyze(task, ivData, spotPattern, v0i, eVsmooth, eVseparate);
        IJ.showProgress(1.0);
    }

    /** This callback method is called when the user changes choices or text fields the dialog. */
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        String inputFilePath = gd.getNextString();
        if (!LeedUtils.fileOk(inputFilePath))
            return false;
        String patternFilePath = gd.getNextString();
        if (patternFilePath.trim().length() > 0 && !LeedUtils.fileOk(patternFilePath))
            return false;
        double v0i = Math.abs(gd.getNextNumber());
        double eVsmooth = gd.getNextNumber();
        double eVseparate = gd.getNextNumber();
        if (!(v0i > 0) || !(eVsmooth >= 0 && eVsmooth < 5*v0i) || eVseparate < 0) return false;
        return true;
    }

    /** Analyzes the equivalent spots in the ivData */
    static void analyze(int task, LeedIVData ivData, LeedSpotPattern spotPattern, double v0i, double eVsmooth, double eVseparate) {
        String title = ivData.title+"_Quality_Statistics";
        Object output = analyze(title, task, ivData.energies, ivData.spotIndices, ivData.data, spotPattern, v0i, eVsmooth, eVseparate);

        if (output instanceof Plot)
            ((Plot)output).show();
        else if (output instanceof ResultsTable)
            ((ResultsTable)output).show(title);
        else
            IJ.error(PLUGIN_NAME, "Error: no symmetry-equivalent curves with sufficient overlap");
    }

    /** Analyzes the equivalent spots in a arrays and, depending on the task,
     *  returns a plot or ResultsTable, or null if no overlaps found.
     *  The 'indices' array contains the spot index for each column in the intensities 'intData' */
    public static Object analyze(String title, int task, double[] energies, int[] indices, double[][] intData,
            LeedSpotPattern spotPattern, double v0i, double eVsmooth, double eVseparate) {

        if (indices.length != intData.length)
            throw new IllegalArgumentException("Inconsistent array lengths "+indices.length+"!="+intData.length);

        double eStep = LeedUtils.getEnergyStep(energies);
        double v0iOverStep = v0i/eStep;
        int rangePoints = (int)Math.round(eVseparate/eStep);
        int smoothRadius = LeedSmoother.invNoiseGainSquareToM(eVsmooth/eStep);
        double actualSmooth = LeedSmoother.mToInvNoiseGainSquare(smoothRadius);
        String paramString = "V0i = "+IJ.d2s(v0i, 1) + " smooth = "+IJ.d2s(actualSmooth*eStep, 1);
        if (rangePoints < 1e7) paramString += " split E to "+IJ.d2s(eVseparate,0);

        //smooth the data and use the smoothed data from now on
        if (smoothRadius > 0) {
            LeedSmoother smoother = new LeedSmoother(smoothRadius);
            double[][] smoothedData = new double[intData.length][];
            for (int i=0; i<indices.length; i++)
                smoothedData[i] = smoother.smooth(intData[i]);
            intData = smoothedData;
        }

        //The following expandable arrays will contain the data for each pair of beams and energy section:
        FloatArray intArray = new FloatArray(indices.length/3 + 10);                //average intensity
        FloatArray rFArray  = new FloatArray(indices.length/3 + 10);                //R_Factor
        FloatArray ovlpArray = new FloatArray(indices.length/3 + 10);               //overlap points
        FloatArray eAvgArray = new FloatArray(indices.length/3 + 10);               //average energy
        LeedIntegerArray groupArray = new LeedIntegerArray(indices.length/3 + 10);  //spot group
        ArrayList<String> spotName1 = new ArrayList<String>(indices.length/3 + 10); //first beam
        ArrayList<String> spotName2 = new ArrayList<String>(indices.length/3 + 10); //2nd beam
        ArrayList<String> rangeStrs = new ArrayList<String>(indices.length/3 + 10); //energy range in readable form
        //expandable arrays for negative intensities (per beam)
        FloatArray iNegativeArray = new FloatArray(indices.length/10 + 10);          //max negative intensity (abs)
        FloatArray eNegativeArray = new FloatArray(indices.length/10 + 10);          //energy range/beam (keV) with negative intensity
        FloatArray iNegativeArray1 = new FloatArray(indices.length/3 + 10);          //max negative intensity (abs) beam 1 per range
        FloatArray eNegativeArray1 = new FloatArray(indices.length/3 + 10);          //energy range/beam (keV) with negative intensity
        FloatArray iNegativeArray2 = new FloatArray(indices.length/3 + 10);          //max negative intensity (abs) beam 2
        FloatArray eNegativeArray2 = new FloatArray(indices.length/3 + 10);          //energy range/beam (keV) with neg. intensity. beam2
        double[] rAndOverlap = new double[LeedRFactor.OUTPUT_SIZE];
        double maxIntensity = 0;

        for (int i=0; i<indices.length; i++) {
            if (indices[i] < 0) continue;
            int group = spotPattern.getGroup(indices[i]);
            String spotI = spotPattern.getName(indices[i]);
            for (int j=i+1; j<indices.length; j++) {
                if (indices[j] < 0 || spotPattern.getGroup(indices[j]) != group) continue;
                int[] ranges = LeedRFactor.getOverlap(intData[i], intData[j]);
                if (ranges == null) continue;
                for (int r=0; r<ranges.length/2; r++) {
                    int iStart = ranges[2*r];
                    int iEnd = ranges[2*r+1];
                    if (iEnd - iStart < 10) continue;        //ignore overlaps < 10 data points
                    int nSubRanges = (int)Math.round((iEnd - iStart)/(double)rangePoints);
                    if (nSubRanges <=0) nSubRanges = 1;
                    for (int s=0; s<nSubRanges; s++) {
                        int iS = iStart + (int)Math.round(s*(iEnd - iStart)/nSubRanges);
                        int iE = iStart + (int)Math.round((s + 1)*(iEnd - iStart)/nSubRanges);
                        rAndOverlap = LeedRFactor.getRFactor(intData[i], intData[j], iS, iE, /*shift=*/0, v0iOverStep, rAndOverlap);
                        rFArray.add((float)rAndOverlap[R_FACTOR_TYPE]);
                        intArray.add((float)rAndOverlap[LeedRFactor.AVG_INTENSITY]);
                        ovlpArray.add((float)rAndOverlap[LeedRFactor.N_OVERLAP]);
                        groupArray.add(group);
                        spotName1.add(spotPattern.getName(indices[i]));
                        spotName2.add(spotPattern.getName(indices[j]));
                        eAvgArray.add((float)(energies[iS]+energies[iE-1])*0.5f);
                        rangeStrs.add(IJ.d2s(energies[iS],1)+"-"+IJ.d2s(energies[iE-1],1));
                        double mostNegative1 = 0, mostNegative2 = 0;
                        int nNegative1 = 0, nNegative2 = 0;
                        for (int p=iS; p<iE; p++) {
                            if (intData[i][p] < 0) {
                                if (intData[i][p] < mostNegative1)      //most negative in the overlap
                                    mostNegative1 = intData[i][p];
                                nNegative1++;
                            }
                            if (intData[j][p] < 0) {
                                if (intData[j][p] < mostNegative2)
                                    mostNegative2 = intData[j][p];
                                nNegative2++;
                            }
                        }
                        iNegativeArray1.add(-(float)mostNegative1);     //abs value
                        eNegativeArray1.add((float)(nNegative1));
                        iNegativeArray2.add(-(float)mostNegative2);
                        eNegativeArray2.add((float)(nNegative2));
                    }
                }
            }
            double mostNegative = 0;
            int nNegative = 0;
            for (int p=0; p<intData[i].length; p++) {       //intensity statistics...
                if (intData[i][p] > maxIntensity)
                    maxIntensity = intData[i][p];           //overall maximum intensity
                if (task == PLOT && intData[i][p] < 0) {
                    if (intData[i][p] < mostNegative)       //most negative per beam (for plot)
                        mostNegative = intData[i][p];
                    nNegative++;
                }
            }
            if (nNegative != 0) {
                iNegativeArray.add(-(float)mostNegative);   //abs value
                eNegativeArray.add((float)(nNegative));
            }
        }
		if (groupArray.size() == 0) return null;


        float[] intensities = intArray.toArray();
        float[] rFactors = rFArray.toArray();
        float[] overlaps = ovlpArray.toArray();
        float[] eAvgs = eAvgArray.toArray();
        float[] iNegatives = iNegativeArray.toArray();
        float[] eNegatives = eNegativeArray.toArray();
        float[] iNegatives1 = task == LIST ? iNegativeArray1.toArray() : null;
        float[] eNegatives1 = task == LIST ? eNegativeArray1.toArray() : null;
        float[] iNegatives2 = task == LIST ? iNegativeArray2.toArray() : null;
        float[] eNegatives2 = task == LIST ? eNegativeArray2.toArray() : null;
        int[] groups = groupArray.toArray();
        //normalize intensities to 0-1000 for average
        float normFactor = (float)(1000./maxIntensity);
        for (int i=0; i<intensities.length; i++)
            intensities[i] *= normFactor;

        if (task == PLOT) {
            return plotRFactorOfEquivalent(title, intensities, rFactors, overlaps,
                    iNegatives, eNegatives, eStep, paramString);
        } else if (task == LIST) {
            ResultsTable outRt = tabulateRFactorOfEquivalent(title, intensities, rFactors, overlaps,
                    eAvgs, groups, spotName1, spotName2, rangeStrs, iNegatives1, eNegatives2, iNegatives2, eNegatives2);
            return outRt;
        } else return null;
    }

    /** Creates a plot with statistics, e.g. the mutual R factors vs. mean intensity
     *  Splits the ranges into ranges of roughly rangePoints
     *  Returns null if there are no overlapping spots */
    static Plot plotRFactorOfEquivalent(String title, float[] intensities, float[] rFactors,
            float[] overlaps, float[] iNegatives, float[] eNegatives, double eStep, String paramString) {
        int n = intensities.length;
        String rFactorName = LeedRFactor.R_FACTOR_NAMES[R_FACTOR_TYPE];
        double[] minMax = Tools.getMinMax(intensities);
        double plotXMin = minMax[0]*0.95;
        double plotXMax = Math.max(minMax[1]*1.05, 100.);
        //normalize sum of overlaps to 100
        double sumOverlaps = 0;
        for (int i=0; i<n; i++)
            sumOverlaps += overlaps[i];
		if (sumOverlaps == 0)
			return null;

        //get cumulated range vs. R
        double[] rsDouble = Tools.toDouble(rFactors);
        int[] rRanks = Tools.rank(rsDouble);
        float[] rSumOverlaps = new float[n];
        float[] sortedRs = new float[n];
        double rSum = 0.;
        for (int i=0; i<n; i++) {
            int index = rRanks[i];
            rSum += overlaps[index];
            rSumOverlaps[i] = (float)(rSum*eStep*1e-3); //overlap in keV
            sortedRs[i] = rFactors[index];
        }

        Plot plot = new Plot(title,
                "Intensity | Overlap (keV)", rFactorName);
        plot.setColor(new Color(0x4020FF), new Color(0xc0c0FF));
        plot.addPoints(intensities, rFactors, null, Plot.CIRCLE, rFactorName+" vs. avg. intensity");
        plot.setColor(new Color(0x000080));
        plot.addPoints(rSumOverlaps, sortedRs, null, Plot.LINE, rFactorName+" vs. total pair E overlap");
        if (iNegatives.length > 0) {
            for (int i=0; i<eNegatives.length; i++)
                eNegatives[i] *= (float)(eStep*0.001);  //to keV, to make it fit the R factor axis
            plot.setColor(Color.RED, Color.RED);
            String legendLabel = "Negative intensities (" +
                    (eNegatives.length > 0 ? "y: keV range)" : "none)");
            plot.addPoints(iNegatives, eNegatives, null, Plot.TRIANGLE, "Negative intensities (y: keV range)");
            double[] minMaxINeg = Tools.getMinMax(iNegatives);
            if (minMaxINeg[0] < plotXMin)               //if there are negative intensities, x range should include at least some
                plotXMin = Math.max(minMaxINeg[0], 0.1*plotXMin);
            if (minMaxINeg[1] < plotXMin)
                plotXMin = 0.995*minMaxINeg[1];
            if (minMaxINeg[1] > plotXMax)
                plotXMax = minMaxINeg[1];
        }
        plot.setColor(Color.BLACK, null);
        plot.setLimits(plotXMin, plotXMax, 0, 1);
        plot.setLogScaleX();
        plot.setJustification(Plot.LEFT);
        String label = paramString;
        label += "\nTotal pair overlap: "+IJ.d2s(sumOverlaps*eStep, 0)+" eV";
        plot.addLabel(0.01, 0.08, label);
        plot.setLegend(null, Plot.AUTO_POSITION);
        return plot;
    }

    /** Creates a ResultsTable with the mutual R factors */
    static ResultsTable tabulateRFactorOfEquivalent(String title, float[] intensities, float[] rFactors,
            float[] overlaps, float[] eAvgs, int[] groups, ArrayList<String> spotName1, ArrayList<String> spotName2, ArrayList<String> rangeStrs,
            float[] iNegatives, float[] eNegatives, float[] iNegatives2, float[] eNegatives2) {

        ResultsTable outRt = new ResultsTable();
        for (int i=0; i<intensities.length; i++) {
            outRt.incrementCounter();
            outRt.addValue("Group", groups[i]);
            outRt.addValue("Beam1", spotName1.get(i));
            outRt.addValue("Beam2", spotName2.get(i));
            outRt.addValue("Range", rangeStrs.get(i));
            outRt.addValue("E_Avg", eAvgs[i]);
            outRt.addValue(LeedRFactor.R_FACTOR_NAMES[R_FACTOR_TYPE], rFactors[i]);
            outRt.addValue("Intensity", intensities[i]);
            outRt.addValue("nPoints", overlaps[i]);
            if (iNegatives.length > 0) {
                outRt.addValue("negative I(1)",   -iNegatives[i]);
                outRt.addValue("neg. I points(1)", eNegatives[i]);
            }
            if (iNegatives.length > 0) {
                outRt.addValue("negative I(2)",   -iNegatives2[i]);
                outRt.addValue("neg. I points(2)", eNegatives2[i]);
            }
        }
        outRt.setDecimalPlaces(4, 1); //E_Avg
        outRt.setDecimalPlaces(5, 3); //RPe
        outRt.setDecimalPlaces(6, 3); //RMod
        outRt.setDecimalPlaces(7, 3); //Intensity
        for (String colStr : new String[] {"negative I(1)", "negative I(2)"}) {
            int col = outRt.getColumnIndex(colStr);
            if (col > 0)
                outRt.setDecimalPlaces(col, 3);
        }
        return outRt;
    }


}
