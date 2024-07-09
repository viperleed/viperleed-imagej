import ij.*;
import ij.IJ;
import ij.gui.*;
import ij.process.*;
import ij.measure.Measurements;
import ij.plugin.filter.MaximumFinder;
import ij.plugin.filter.GaussianBlur;
import ij.util.Tools;
import java.awt.*;
import java.awt.event.*;
import java.util.Vector;
import java.util.Arrays;
/**
 *  Part of the LEED I/V package
 *
 * The dialog for selecting the slice where spot indices will be defined
 * and the significance limit for detection of spots.
 * In interactive operation, called from the LeedIndexSelector.
 *
 */

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
 *  @author Michael Schmid, IAP/TU Wien, 2019-2024
 */

public class LeedIndexSliceSelector  implements DialogListener, ImageListener {
    public static final int MIN_N_MAXIMA = 3;        //we need three maxima for a basis
    LEED_Spot_Tracker spotTracker;
    LeedIndexSelector indexSelector;
    ImagePlus stackImp;
    ImagePlus maskImp;
    Roi maskRoi;
    double[] energies;
    double[] screenFitterArray;     // previous relation between k space and screen coordinates (if any)
    double minSignificance;
    double[][] xyMax;               // pixel coordinates of maxima (spots) found (one array each for x & y)
    int nMaxima;                    // number of maxima in xMax, yMax
    LeedScreenFitter screenFitter;  // refined relation between k space and screen coordinates (if any)
    NonBlockingGenericDialog gd;
    boolean blockOkButton;


    /** Asks the user to set the slice and minSignificance and returns the maxima found,
     *  or, if a previous screenFitterArray is provided and the user says so, the
     *  LeedScreenFitter describing the relation between k-space and screen coordinates.
     *  When returning the maxima as double[][] array (index=0, 1 for x, y), these
     *  maxima are sorted in the sequence of decreasing pixel value.
     *  Returns null on error or if cancelled. */
    public Object getMaxima(LEED_Spot_Tracker spotTracker, LeedIndexSelector indexSelector,
            ImagePlus stackImp, ImagePlus maskImp, Roi maskRoi, double[] energies, double[] screenFitterArray) {
        try {
            this.spotTracker = spotTracker;
            this.indexSelector = indexSelector;
            this.stackImp = stackImp;
            this.maskImp = maskImp;
            this.maskRoi = maskRoi;
            this.energies = energies;
            this.screenFitterArray = screenFitterArray;

            minSignificance = LeedParams.get(LeedParams.MINSIGNIFICANCEINDEX);
            stackImp.getWindow().toFront();
            gd = new NonBlockingGenericDialog("Select Stack Slice");
            gd.setInsets(5, 0, 0); //top, left, bottom
            gd.addMessage("Select a stack slice where as many spots as possible\n"+
                    "are visible all over the screen, also near the edges\n"+
                    "(minimum "+MIN_N_MAXIMA+").\n"+
                    "You should be able to identify some of the spots.");
            gd.setInsets(20, 5, 0); //top, left, bottom
            gd.addSlider("Noise rejection *", 1.5, 5.0, minSignificance);
            //addNumericField("Noise rejection", minSignificance, 2, 6, "(1-10, usually ~2) *");
            gd.setInsets(0, 0, 0); //top, left, bottom
            gd.addMessage("* 0.5-10, typically ~2.\n"+
                    "Set high enough to suppress most noise maxima\n"+
                    "(you may want to check different stack frames).\n"+
                    "In case of problems, it may help to cancel and\n"+
                    "adjust the radius to better fit the spots.");

            if (screenFitterArray != null)
                gd.enableYesNoCancel("Labels are OK", "Name Spots...");

            gd.addHelp(LeedSpotTrackerHelp.getSetIndices1Help());
            gd.addDialogListener(this);

            findAndShowSpotMaxima();    //enables/disables the OK button
            ImagePlus.addImageListener(this);
            spotTracker.setCurrentFrontDialog(gd);
            gd.showDialog();
            spotTracker.setCurrentFrontDialog(null);
            ImagePlus.removeImageListener(this);

            if (gd.wasCanceled()) return null;

            minSignificance = gd.getNextNumber();
            LeedParams.set(LeedParams.MINSIGNIFICANCEINDEX, minSignificance);

            if (nMaxima < MIN_N_MAXIMA) {
                return null;
            } else {
                if (screenFitterArray != null && gd.wasOKed()) { //user said "Labels are OK"
                    return screenFitter;
                } else {                                         //user has to name the spots
                    LeedOverlay.removeNameOverlays(stackImp);
                    return xyMax;
                }
            }
        } catch (Exception e) {IJ.handleException(e); return null;}
    }

    /** This callback method is called when the user changes fields the dialog. */
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        try{
            if (e == null) return true;                     //ignore call when closing the dialog
            minSignificance = gd.getNextNumber();
            blockOkButton = gd.invalidNumber() || minSignificance < 0.5 || minSignificance > 10;
            if (blockOkButton || stackImp == null)
                return false;
            LeedParams.set(LeedParams.MINSIGNIFICANCEINDEX, minSignificance);
            findAndShowSpotMaxima();                        //enables/disables ok button
            Button okButton = gd.getButtons()[0];
            return okButton.isEnabled();
        } catch(Exception ex) {IJ.handleException(ex);return false;}
    }

    /** required by ImageListener interface, not used */
    public void imageOpened(ImagePlus imp) {}

    public void imageClosed(ImagePlus imp) {
        if (imp == stackImp || imp == maskImp) {
            stackImp = null;
            Button okButton = gd.getButtons()[0];
            okButton.setEnabled(false);         //no 'OK' when the stack or mask disappears.
            okButton.setVisible(false);
            if (screenFitterArray != null) {    //we have 'Labels OK' and No='Name Spots' buttons
                Button noButton = gd.getButtons()[2];
                noButton.setEnabled(false);
                noButton.setVisible(false);
            }
            gd.toFront();
        }
    }

    /** When the stack or mask changes, find the maxima again */
    public void imageUpdated(ImagePlus imp) {
        try {
            if (imp == stackImp || imp == maskImp)
                findAndShowSpotMaxima();
        } catch(Exception ex) {IJ.handleException(ex);}
    }

    /** Detects and shows the maxima in the current slice and sets the xMax, yMax and nMaxima.
     *  Enables the 'OK' button if there are enough maxima. */
    void findAndShowSpotMaxima() {
        if (stackImp == null) return;

        int stackSlice = stackImp.getCurrentSlice();
        double energy = LEED_Spot_Tracker.sliceToEnergy(energies, stackSlice);
        double radius = LeedRadiusSelector.radius(energy, false);
        int spotBackgrShape = (int)LeedParams.get(LeedParams.BACKGROUNDTYPE);
        double azBlurRadians = Math.toRadians(LeedParams.get(LeedParams.AZIMUTHBLURANGLE));

        xyMax = findSpotMaxima(stackImp, maskImp, maskRoi, spotBackgrShape, radius, azBlurRadians, minSignificance);
        nMaxima = xyMax[0].length;
        screenFitter = screenFitterArray == null ? null :
                indexSelector.useScreenFitterArray(xyMax[0], xyMax[1], stackSlice, energy, /*searchRadius:default*/Double.NaN, /*showLabels=*/true);

        showMaxOverlay(stackImp, stackSlice, xyMax, radius);

        Button okButton = gd.getButtons()[0];
        if (screenFitterArray != null) {    //we have 'Labels OK' and No='Name Spots' buttons
            Button noButton = gd.getButtons()[2];
            noButton.setEnabled(nMaxima > MIN_N_MAXIMA && !blockOkButton);
            okButton.setEnabled(nMaxima > MIN_N_MAXIMA && screenFitter != null && !blockOkButton);
        } else
            okButton.setEnabled(nMaxima > MIN_N_MAXIMA && !blockOkButton);
    }

    /** Shows circles at the spots found */
    public static void showMaxOverlay(ImagePlus stackImp, int stackSlice, double[][] xyMax, double radius) {
        LeedOverlay.removeSpotOverlays(stackImp);
        for (int i=0; i<xyMax[0].length; i++)    //show circles at the maxima
            LeedOverlay.add(stackImp, stackSlice, xyMax[0][i], xyMax[1][i], radius, null);
        stackImp.draw();
    }

    /** Detects the maxima within the mask roi and returns their coordinates */
    public static double[][] findSpotMaxima(ImagePlus stackImp, ImagePlus maskImp,
            Roi maskRoi, int spotBackgrShape, double radius, double azBlurRadians, double minSignificance) {
        if (maskImp == null || maskRoi == null)
            throw new RuntimeException("Error: Finding maxima (spots) requires a mask");
        Rectangle maskRect = maskRoi.getBounds();
        double maskCenterX = maskRect.x + 0.5*maskRect.width;
        double maskCenterY = maskRect.y + 0.5*maskRect.height;

        FloatProcessor stackIp = (FloatProcessor)stackImp.getProcessor();
        FloatProcessor ip = (FloatProcessor)stackIp.duplicate(); //used temporarily to find maxima
        //assuming sigma is radius/3, smoothing with sigma=radius/6 makes almost no change
        (new GaussianBlur()).blurGaussian(ip, /*sigma=*/0.167*radius);
        ip.setRoi(maskRoi);
        ImageStatistics is = ImageStatistics.getStatistics(ip, Measurements.MEAN+Measurements.MIN_MAX, stackImp.getCalibration());
        ip.setThreshold(is.min, is.mean, ImageProcessor.NO_LUT_UPDATE);     //needed for better resolution of histogram
        ImageStatistics is2 = ImageStatistics.getStatistics(ip, Measurements.LIMIT+Measurements.MODE, stackImp.getCalibration());
        double minProminence = 0.02*(is2.dmode-is.min)+5e-5*(is.max-is.min);//rather sensitive settings for MaximumFinder
        if (IJ.debugMode) new ImagePlus("findMaximaHere"+LeedUtils.d2s(minProminence,3), ip).show();
        Polygon maxCoords = (new MaximumFinder()).getMaxima(ip, minProminence, /*strict=*/true, /*excludeOnEdges=*/true);
        int nMax = maxCoords.npoints;
        if (IJ.debugMode) IJ.log("findSpotMaxima minProminence="+LeedUtils.d2s(minProminence,3)+" dmode="+IJ.d2s(is2.dmode)+" (in "+IJ.d2s(is2.histMin)+"-"+IJ.d2s(is2.histMax)+") min-max="+IJ.d2s(is.min)+"-"+IJ.d2s(is.max)+" nMax="+nMax);
        double[] xMax = new double[nMax];
        double[] yMax = new double[nMax];
        ByteProcessor maskIp = LeedUtils.getMaskIp(maskImp);
        double[] spotData = new double[LeedSpotAnalyzer.OUTPUT_SIZE];
        int nMaxima = 0;
        for (int i=0; i<nMax; i++) { //we start at the highest maximum
            spotData = LeedSpotAnalyzer.centerAndAnalyzeSpot(maxCoords.xpoints[i], maxCoords.ypoints[i], maskCenterX, maskCenterY,
                spotBackgrShape, radius, azBlurRadians, minSignificance, stackIp, maskIp, spotData);
            //if (i<10)IJ.log("max@"+maxCoords.xpoints[i]+","+maxCoords.ypoints[i]+" -> "+IJ.d2s(spotData[LeedSpotAnalyzer.X])+","+IJ.d2s(spotData[LeedSpotAnalyzer.Y]));
            if (Double.isNaN(spotData[LeedSpotAnalyzer.X])) continue; //no maximum found at this position
            xMax[nMaxima] = spotData[LeedSpotAnalyzer.X];
            yMax[nMaxima] = spotData[LeedSpotAnalyzer.Y];
            boolean maxOk = true;
            for (int j=0; j<nMaxima; j++) {
                double distanceSqr = sqr(xMax[j] - xMax[nMaxima]) + sqr(yMax[j] - yMax[nMaxima]);
                if (distanceSqr < sqr(2.5*radius)) { //too close to a previously found spot?
                maxOk = false;
                break;
                }
            }
            if (maxOk) nMaxima++;
        }
        IJ.showStatus(nMaxima+" spots");
        double[][] xyMax = new double[2][nMaxima];
        System.arraycopy(xMax, 0, xyMax[0], 0, nMaxima);
        System.arraycopy(yMax, 0, xyMax[1], 0, nMaxima);
        return xyMax;
    }

    static double sqr(double x) {return x*x;}
}
