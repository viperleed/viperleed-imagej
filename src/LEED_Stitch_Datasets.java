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
 *  Stitches two datasets, i.e., two sets of I(V) curves.
 *  There must be some overlap between the IV curves; in the overlap region the
 *  curves are averaged with weights that lead to a smooth transition.
 *  Curves only present in one dataset are copied and filled with NaNs in
 *  the energy range where no data are present.
 *  Normalization is done by multiplying the second dataset (i.e., for ascending
 *  energies, the one with the higher energies) to have equal average intensity in
 *  the overlap region. The scale factor is the same for all beams (and calculated
 *  form all beams where both input files have values at the same energy).
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
 *  @author Michael Schmid, IAP/TU Wien, 2021-2024
 */


public class LEED_Stitch_Datasets implements PlugIn {
    static final String PLUGIN_NAME = "Stitch LEED I(V) Curves";
    static String[] filePaths = new String[2];

    public void run(String arg) {
        if (IJ.versionLessThan("1.52u")) return;


        // I N I T I A L   D I A L O G
        GenericDialog gd = new GenericDialog("Stitch Data Sets of IV curves");
        gd.addFileField("File_1", filePaths[0], 35);
        gd.addFileField("File_2", filePaths[1], 35);
        gd.addMessage("Note: energy steps must be equal for both;\nthere must be some overlap between the data sets.\n"+
                "You can drag&drop .csv files onto the 'File' fields.");

        DialogListener dialogListener = new DialogListener() {  //only for checking whether the files exist
            public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
                for (int i=0; i<filePaths.length; i++) {
                    String path = gd.getNextString();
                    if (!LeedUtils.fileOk(path))
                        return false;
                }
                return true;
            }
        };
        gd.addDialogListener(dialogListener);
        Button okButton = gd.getButtons()[0];
        okButton.setEnabled(dialogListener.dialogItemChanged(gd, null));

        gd.showDialog();
        if (gd.wasCanceled()) return;

        filePaths[0] = gd.getNextString();
        filePaths[1] = gd.getNextString();

        String[] fileNames = new String[2];
        String[] directories = new String[2];
        for (int i=0; i<fileNames.length; i++) {
            File file = new File(filePaths[i]);
            fileNames[i] = file.getName();
            directories[i] = file.getParent();
        }

        // R E A D   F I L E S   &   A N A L Y Z E
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

        double[] firstInc0 = LeedUtils.getFirstAndIncrement(ivData[0].energies);
        double[] firstInc1 = LeedUtils.getFirstAndIncrement(ivData[1].energies);
        double step0 = firstInc0[1];
        double step1 = firstInc1[1];
        if (!(Math.abs((step1-step0)/step1) < 1e-2)) {   //NaN if not evenly spaced
            IJ.error(PLUGIN_NAME, "Error: Energies must be evenly spaced with the same step in both files:\n"+
                    fileNames[0]+": "+(Double.isNaN(step0)?"not evenly spaced":"step="+IJ.d2s(step0, 4))+"\n"+
                    fileNames[1]+": "+(Double.isNaN(step1)?"not evenly spaced":"step="+IJ.d2s(step1, 4)));
            return;
        }

        double step = 0.5*(step0 + step1);
        double first0 = firstInc0[0];
        double first1 = firstInc1[0];
        double last0 = firstInc0[0] + step0*(ivData[0].energies.length-1);
        double last1 = firstInc1[0] + step1*(ivData[1].energies.length-1);
        boolean firstIs0 = (first1 - first0)*step > 0;  // whether dataset 0 is first (earlier start)
        boolean lastIs0 = (last0 - last1)*step > 0;     // whether dataset 0 is last (later end)
        if (firstIs0 == lastIs0) {
            IJ.error(PLUGIN_NAME, "Error: Not two successive energy ranges:\n"+
                    fileNames[0]+": "+(float)first0+"-"+(float)last0+"\n"+
                    fileNames[1]+": "+(float)first1+"-"+(float)last1);
            return;
        }
        int firstIndex = firstIs0 ? 0 : 1;
        int secondIndex = firstIs0 ? 1 : 0;

        double overlapERange = firstIs0 ? last0 - first1 : last1 - first0;
        double nOverlapD = overlapERange/step + 1;
        int nOverlap = (int)Math.round(nOverlapD);
        if (nOverlap < 1) {
            IJ.error(PLUGIN_NAME, "Error: No energy overlap:\n"+
                    fileNames[0]+": "+(float)first0+"-"+(float)last0+"\n"+
                    fileNames[1]+": "+(float)first1+"-"+(float)last1);
            return;
        }
        int firstOverlapStart = ivData[firstIndex].energies.length - nOverlap;
        int nEnergies = firstOverlapStart + ivData[secondIndex].energies.length;

        // M E R G E   S P O T   P A T T E R N S   (works also without spot groups)
        LeedSpotPattern[] spotPatterns = new LeedSpotPattern[ivData.length];
        for (int i=0; i<ivData.length; i++) {
            spotPatterns[i] = new LeedSpotPattern(ivData[i].spotNames, /*groupsRequired=*/false);
            String error = spotPatterns[i].getErrorMessage();
            if (error != null) {
                 IJ.error(PLUGIN_NAME, error);
                 return;
            }
        }
        LeedSpotPattern spotPattern = new LeedSpotPattern(spotPatterns[0], spotPatterns[1], /*groupsRequired=*/false);
        String error = spotPattern.getErrorMessage();
        if (error != null) {
            IJ.error(PLUGIN_NAME, error);
            return;
        }
        for (int id=0; id<ivData.length; id++) {
            ivData[id].setSpotPattern(spotPattern, /*showWarning=*/true);
            if (ivData[id].spotIndices == null)
                return;
        }
        int[][] col = new int[2][spotPattern.size()];
        for (int iSpot=0; iSpot<spotPattern.size(); iSpot++) {
            for (int id=0; id<ivData.length; id++)
                col[id][iSpot] = LeedUtils.arrayIndexOf(ivData[id].spotIndices, iSpot);
        }
        String[] spotNamesArray = spotPattern.getAllSpotNamesWithGroup();

        // P R E P A R E   N O R M A L I Z A T I O N
        double firstSum=0, secondSum=0;                 // for comparing intensities in the overlap region
        for (int iSpot=0; iSpot<spotPattern.size(); iSpot++) {
            if (col[0][iSpot] < 0 || col[1][iSpot] < 0) continue;
            int firstCol  = col[firstIndex][iSpot];
            int secondCol = col[secondIndex][iSpot];
            for (int i=0; i<nOverlap; i++) {
                double firstV  = ivData[firstIndex].data[firstCol][firstOverlapStart+i];
                double secondV = ivData[secondIndex].data[secondCol][i];
                if (firstV > 0 && secondV > 0) {
                    firstSum  += firstV;
                    secondSum += secondV;
                }
            }
        }
        if (firstSum == 0) {
            IJ.error(PLUGIN_NAME, "Error: No positive values in overlap");
            return;
        }
        double firstOverSecond = firstSum/secondSum;

        // C A L C U L A T E   O U T P U T
        double[] outEnergies = new double[nEnergies];
        for (int i=0; i<nEnergies; i++)
            outEnergies[i] = i < ivData[firstIndex].energies.length - nOverlap/2 ?
                    ivData[firstIndex].energies[i] :
                    ivData[secondIndex].energies[i-firstOverlapStart];

        double[][] outData = new double[spotPattern.size()][nEnergies];
        for (int iSpot=0; iSpot<spotPattern.size(); iSpot++) {
            int firstCol  = col[firstIndex][iSpot];
            int secondCol = col[secondIndex][iSpot];
            if (firstCol >= 0 && secondCol >= 0) {  //common column, smooth transition
                for (int i=0; i<firstOverlapStart; i++)
                    outData[iSpot][i] = ivData[firstIndex].data[firstCol][i];
                for (int i=firstOverlapStart, j=0; i<ivData[firstIndex].energies.length; i++, j++) {
                    double weight = (j+1)*(1.0/(nOverlap+1));   //e.g. for nOverlap = 4: weights 1/4, 2/4, 3/4
                    outData[iSpot][i] = ivData[firstIndex].data[firstCol][i]*(1.-weight)
                            + ivData[secondIndex].data[secondCol][j]*weight*firstOverSecond;
                }
                for (int i=ivData[firstIndex].energies.length, j=nOverlap; i<nEnergies; i++,j++)
                    outData[iSpot][i] = ivData[secondIndex].data[secondCol][j]*firstOverSecond;
            } else if (firstCol >= 0) {
                for (int i=0; i<ivData[firstIndex].energies.length; i++)
                    outData[iSpot][i] = ivData[firstIndex].data[firstCol][i];
                for (int i=ivData[firstIndex].energies.length; i<nEnergies; i++)
                    outData[iSpot][i] = Double.NaN;
            } else if (secondCol >= 0) {
                for (int i=0; i<firstOverlapStart; i++)
                    outData[iSpot][i] = Double.NaN;
                for (int i=firstOverlapStart, j=0; i<nEnergies; i++, j++)
                    outData[iSpot][i] = ivData[secondIndex].data[secondCol][j]*firstOverSecond;
            } else
                throw new RuntimeException("Stitch datasets: Internal error, column '"+spotPattern.getName(iSpot)+"' not found");
        }
        LeedIVData ivDataOut = new LeedIVData(/*title=*/null, outEnergies, outData, spotNamesArray, null);

        // O U T P U T   F I L E
        String firstFileName = new File(filePaths[firstIndex]).getName();
        String secondFileName = new File(filePaths[secondIndex]).getName();
        int commonBeginLength = LeedUtils.nMatchingBeginCharacters(firstFileName, secondFileName);
        if (commonBeginLength < 2)
            commonBeginLength = 0;
        int commonEndLength = LeedUtils.nMatchingEndCharacters(firstFileName, secondFileName);
        String outName = firstFileName.substring(0, firstFileName.length()-commonEndLength)+"+"+secondFileName.substring(commonBeginLength);
        String savedDir = null;
        if (directories[0].equals(directories[1])) { //if both are in the same directory, suggest this dir for output
            savedDir = OpenDialog.getDefaultDirectory();
            OpenDialog.setDefaultDirectory(directories[0]);
        }
        IJ.showStatus("Save stitched as?");
        SaveDialog sd = new SaveDialog("Save stitched data as", outName, ".csv");
        IJ.showStatus("");
        String outDirectory = sd.getDirectory();
        String outFileName = sd.getFileName();
        if (savedDir != null)
            OpenDialog.setDefaultDirectory(savedDir);
        if (outFileName == null) return;
        ivDataOut.save(outDirectory+outFileName);
    }

}
