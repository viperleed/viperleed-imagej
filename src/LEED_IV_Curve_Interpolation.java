import ij.*;
import ij.gui.*;
import ij.io.*;
import ij.measure.*;
import ij.process.*;
import ij.text.*;
import ij.plugin.*;
import ij.util.Tools;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;

/**
 *  This class provides the user interface for interpolating I(V) curves to new energy steps.
 *  The actual interplation is done in the LeedInterpolator class.
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


public class LEED_IV_Curve_Interpolation implements PlugIn {
    /** Interpolation type: Cubic splines as implemented in ImageJ */
    public static final int CUBIC_SPLINE = 0;
    /** Interpolation type: Catmull-Rom, i.e., smooth cubic curves like bilinear interpolation in images */
    public static final int CUBIC_SMOOTH = 1;
    /** Interpolation type: Cubic curves by piecewise fitting */
    public static final int CUBIC_FIT = 2;
    /** Interpolation type: nearest neighbor */
    public static final int NEAREST = 3;
    /** Interpolation type: windowed sinc with 6 points input. The sinc interpolations are for testing and do not properly handle the boundaries */
    public static final int SINC6 = 4;
    /** Interpolation type: windowed sinc with 12 points input, mimicking a cubic spline. No proper handling of bopundaries. */
    public static final int SINC12 = 5;
    /** Interpolation type: windowed sinc with 20 points input, mimicking quintic spline. No proper handling of boundaries. */
    public static final int SINC20 = 6;
    static String filePathS;
    static double newStepS=0.5;

    int interpolationType = CUBIC_SPLINE;

    public void run(String arg) {
        if (IJ.versionLessThan("1.52u")) return;

        // I N I T I A L   D I A L O G
        String filePath = filePathS;
        double newStep = newStepS;
        GenericDialog gd = new GenericDialog("LEED IV Curve Interpolation");
        gd.addFileField("File (*.csv)", filePath, 35);
        gd.addNumericField("New step size", newStep, 2);
        gd.addMessage("You can drag&drop files onto the 'File' field.\n");

        DialogListener dialogListener = new DialogListener() {  //only for checking whether the file exist
            public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
                String path = gd.getNextString();
                if (!LeedUtils.fileOk(path))
                    return false;
                double newStep = gd.getNextNumber();
                if (!(newStep > 0))
                    return false;
                return true;
            }
        };
        gd.addDialogListener(dialogListener);
        Button okButton = gd.getButtons()[0];
        okButton.setEnabled(dialogListener.dialogItemChanged(gd, null));

        gd.showDialog();
        if (gd.wasCanceled()) return;

        filePath = gd.getNextString();
        newStep = gd.getNextNumber();
        filePathS = filePath;
        newStepS = newStep;

        // P R O C E S S
        File file = new File(filePath);
        String filename  = file.getName();
        String directory = file.getParent();

        IJ.showStatus("Reading data");
        LeedIVData ivData = LeedIVData.fromFile(filePath, null, LeedIVData.ZERO_IS_NAN);
        if (ivData == null || ivData.energies == null)
            return;
      
        LeedIVData newData = LeedInterpolator.interpolate(ivData, newStep);
        
        IJ.showStatus("Saving...");
        int energyDigits = LeedUtils.getEnergyDigits(newData.energies);
        String newFilename = LeedUtils.removeExtension(filename);
        int stepIndex = newFilename.indexOf("_step");
        if (stepIndex > 0)                  //remove previous step size
            newFilename = newFilename.substring(0, stepIndex);
        newFilename += "_step_"+IJ.d2s(newStep, energyDigits)+".csv";
        String newPath = directory + File.separator + newFilename;
        newData.save(newPath);
        IJ.showStatus("");
    }
}
