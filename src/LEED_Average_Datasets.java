import ij.*;
import ij.gui.*;
import ij.io.*;
import ij.measure.*;
import ij.text.*;
import ij.plugin.*;
import ij.util.Tools;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;

/**
 *  Averages up to five datasets, i.e., sets of I(V) curves.
 *  IV curves present in only one of the datasets are copied.
 *  Uses soft fade in/fade out if the curves start or end at
 *  energies that are further apart than V0i; otherwise sets
 *  the limit accroding to the common energy range.
 */

/** This ImageJ plugin is part of the ViPErLEED package for LEED I(V) analysis.
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

public class LEED_Average_Datasets implements PlugIn {
    static final String PLUGIN_NAME = "Average LEED I(V) Curves";
    static final int N_MAX = 6;                     //we have up to six curves
    static final int MIN_OVERLAP = 2;               //minimum 2 points overlap
    static String[] filePaths = new String[N_MAX];  //remember for next invocation
    static double v0i = LeedParams.get(LeedParams.V0I);

    public void run(String arg) {
        if (IJ.versionLessThan("1.52u")) return;

        // I N I T I A L   D I A L O G
        GenericDialog gd = new GenericDialog("Average Sets of IV curves");
        for (int i=0; i<N_MAX; i++)
            gd.addFileField("File_"+(i+1), filePaths[i], 50);
        gd.addNumericField("V0i", v0i, 1, 4, " eV");
        gd.addMessage("V0i determines the energy scale for smooth fade in/fade out");

        DialogListener dialogListener = new DialogListener() {  // for checking whether the files exist and the weight is ok
            public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
                int nFiles = 0;
                for (int i=0; i<N_MAX; i++) {
                    String path = gd.getNextString();
                    if (path.trim().length()>0) {
                        if (!LeedUtils.fileOk(path))
                            return false;
                        nFiles++;
                    }
                }
                if (nFiles < 2) return false;
                double v0i = Math.abs(gd.getNextNumber());
                if (!(v0i > 0)) return false;
                return true;
            }
        };
        gd.addDialogListener(dialogListener);
        Button okButton = gd.getButtons()[0];
        okButton.setEnabled(dialogListener.dialogItemChanged(gd, null));

        gd.showDialog();
        if (gd.wasCanceled()) return;

        for (int i=0; i<N_MAX; i++)
            filePaths[i] = gd.getNextString();
        v0i = Math.abs(gd.getNextNumber());

        // R E A D   F I L E S
        String directory = null;
        ArrayList<LeedIVData> allIVData = new ArrayList<LeedIVData>(N_MAX);
        for (int i=0; i<N_MAX; i++) {
            IJ.showProgress(0.1+i*0.7/N_MAX);
            if (filePaths[i].trim().length() == 0) continue;
            LeedIVData ivd = LeedIVData.fromFile(filePaths[i], null, LeedIVData.E_ASCENDING | LeedIVData.ZERO_IS_NAN);
            if (ivd == null) {
                IJ.showProgress(1.0);
                return;
            } else
                allIVData.add(ivd);
            if (directory == null)
                directory = (new File(filePaths[i])).getParent();
        }
        if (allIVData.size() < 2) {
            IJ.error(PLUGIN_NAME, "Error: Only "+allIVData.size()+" data sets, nothing to average");
            IJ.showProgress(1.0);
            return;
        }
        IJ.showProgress(0.8);

        // A N A L Y Z E   E N E R G I E S ,   I N T E R P O L A T E   I F   N E E D E D
        double minE=Double.MAX_VALUE, maxE=-Double.MAX_VALUE;
        double minStep = Double.MAX_VALUE;
        boolean differentEnergyAxes = false;
        boolean first = true;
        for (LeedIVData ivData : allIVData) {
            double[] energies = ivData.energies;
            double[] minMax = Tools.getMinMax(energies);
            if (!first && (minE != minMax[0]) || maxE != minMax[1])
                differentEnergyAxes = true;

            if (minMax[0] < minE)
                minE = minMax[0];
            if (minMax[1] > maxE)
                maxE = minMax[1];

            double step = LeedUtils.getEnergyStep(energies);
            if (!first && step != minStep)
                differentEnergyAxes = true;
            if (step < minStep)
                minStep = step;

            first = false;
        }
        if (differentEnergyAxes)
            for (int i=0; i<allIVData.size(); i++)
                allIVData.set(i, LeedInterpolator.interpolate(allIVData.get(i), minE, maxE, minStep));

        double[] energies = allIVData.get(0).energies;      //we have a common energy axis now
        int nEnergies = energies.length;
        double energyStep = LeedUtils.getEnergyStep(energies);
        IJ.showProgress(0.82);


        // A V E R A G E   F O R   A L L   S P O T S ;   S O R T E D   B Y   G R O U P
        ArrayList<LeedSpotPattern> allSpotPatterns = new ArrayList<LeedSpotPattern>(allIVData.size());
        int firstGroup = Integer.MAX_VALUE, lastGroup = Integer.MIN_VALUE;
        for (LeedIVData ivData : allIVData) {
            LeedSpotPattern spotPattern = new LeedSpotPattern(ivData.spotNames, /*groupsRequired=*/false);
            allSpotPatterns.add(spotPattern);
            for (int i=0; i<spotPattern.size(); i++) {
                int g = spotPattern.getGroup(i);
                if (g < firstGroup) firstGroup = g;
                if (g > lastGroup)  lastGroup = g;
            }
        }
        int initialOutputSize = allIVData.get(0).spotNames.length*2;
        ArrayList<String> namesOut  = new ArrayList<String>(initialOutputSize);
        ArrayList<double[]> dataOut = new ArrayList<double[]>(initialOutputSize);
        LeedIntegerArray iStarts = new LeedIntegerArray(12);
        LeedIntegerArray iEnds    = new LeedIntegerArray(12);
        for (int group=firstGroup; group<=lastGroup; group++) {
            IJ.showProgress(0.82 + 1.17/(lastGroup-firstGroup)*(group-firstGroup));
            LinkedHashSet<String> spotNames = new LinkedHashSet<String>(12);
            for (int id=0; id<allIVData.size(); id++) {     //merge spot names of all input data into one set per group
                LeedIVData ivData = allIVData.get(id);
                LeedSpotPattern spotPattern = allSpotPatterns.get(id);
                int[] spots = spotPattern.getAllSpotsForGroup(group);
                if (spots.length == 0) continue;
                for (int spotIndex : spots)
                    spotNames.add(spotPattern.getNameWithGroup(spotIndex, /*replaceComma=*/true));
            }
            if (spotNames.isEmpty()) continue;
            int iStart=Integer.MAX_VALUE, iEnd = Integer.MIN_VALUE;
            for (String spotName : spotNames) {
                ArrayList<double[]> dataToAverage = new ArrayList<double[]>(allIVData.size());
                for (int id=0; id<allIVData.size(); id++) { //find data for the given spot name in all input files
                    LeedIVData ivData = allIVData.get(id);
                    LeedSpotPattern spotPattern = allSpotPatterns.get(id);
                    int spotIndex = spotPattern.getIndex(spotName);
                    if (spotIndex < 0) continue;
                    if (ivData.spotIndices != null) {       //if there is no 1:1 correspondence between spotNames array and spotPattern
                        spotIndex = LeedUtils.arrayIndexOf(ivData.spotIndices, spotIndex);
                        if (spotIndex < 0) continue;
                    }
                    int[] startEnd = LeedUtils.getStartEndNum(ivData.data[spotIndex]);
                    if (startEnd[0] < 0) continue;          //ignore empty data columns
                    dataToAverage.add(ivData.data[spotIndex]);
                    iStarts.add(startEnd[0]);               //find start & end of data ranges
                    iEnds.add(startEnd[1]);
                    if (startEnd[0] < iStart)
                        iStart = startEnd[0];
                    if (startEnd[1] > iEnd)
                        iEnd = startEnd[1];
                }
                if (iStarts.size() == 0) continue;          //ignore spots without any valid data
                                                           //avoid fade in/out if a small change of the limits can avoid it
                int iStart1 = snapStartOrEnd(iStarts, iStart, (int)Math.round(v0i/energyStep));
                int iEnd1   = snapStartOrEnd(iEnds,   iEnd,   (int)Math.round(v0i/energyStep));
                double[] average = LeedCurveAverager.calculateAverage(dataToAverage, iStart1, iEnd1, MIN_OVERLAP, v0i/energyStep);
                if (average != null && LeedUtils.countNonNaN(average) > 0) {
                    if (iStart1 != iStart || iEnd1 != iEnd)
                        average = LeedUtils.restrictRange(average, iStart1, iEnd1);
                    namesOut.add(spotName);                 //averaging done, add to output
                    dataOut.add(average);
                }
                iStarts.clear();
                iEnds.clear();
            }
        }
        IJ.showProgress(1.0);
        if (namesOut.size() == 0) {
            IJ.error (PLUGIN_NAME, "Error: averaging results in no valid data");
            return;
        }
        String[] headingsOut = namesOut.toArray(new String[namesOut.size()]);
        LeedSpotPattern spotPatternOut = new LeedSpotPattern(headingsOut, /*groupsRequired=*/false);  //check spot names for duplicate spots and validity
        if (spotPatternOut.size() != headingsOut.length) {
            IJ.error(PLUGIN_NAME, "Error: "+spotPatternOut.getErrorMessage());
            IJ.showProgress(1.0);
            return;

        }
        LeedIVData ivDataOut = new LeedIVData("averaged", energies, dataOut.toArray(new double[dataOut.size()][]),
                headingsOut, /*spotIndices=*/null);

        // O U T P U T
        String commonName = null;
        for (LeedIVData ivData : allIVData)
            commonName = commonName==null ? ivData.title :
                LeedUtils.getCommon(commonName, ivData.title);

        String outName = commonName+"_avg.csv";
        IJ.showStatus("Save average as?");
        SaveDialog sd = new SaveDialog("Save average as", outName, ".csv");
        IJ.showStatus("");
        String outDirectory = sd.getDirectory();
        String outFileName = sd.getFileName();
        if (outFileName != null)
            ivDataOut.save(outDirectory+outFileName);
    }

    /** Returns iBound snapped to the most different value of the 'iBounds' array,
     *  as long as that value is within 'plusMinus' from the original value */
    int snapStartOrEnd(LeedIntegerArray iBounds, int iBound, int plusMinus) {
        int delta = 0;
        for (int i=0; i<iBounds.size(); i++) {
            int diff = iBounds.get(i) - iBound;
            if (Math.abs(diff) <= plusMinus && Math.abs(diff) > Math.abs(delta))
                delta = diff;
        }
        return iBound+delta;
    }

}
