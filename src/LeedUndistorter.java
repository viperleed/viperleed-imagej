import ij.*;
import ij.IJ;
import ij.gui.*;
import ij.process.*;
import ij.measure.Measurements;
import ij.plugin.filter.MaximumFinder;
import ij.plugin.filter.GaussianBlur;
import ij.util.Tools;
import java.lang.ref.*;
import java.awt.*;
import java.awt.event.*;
import java.util.Vector;
import java.util.Arrays;
/**
 *  Part of the LEED I/V package
 *
 *  This class can undistort an image based on the distortions determined
 *  when fitting the basis.
 *  It can either undistort a single image or the whole stack; in the latter case,
 *  the output is a Virtual Stack.
 *  When undistorting the whole stack, one can have either
 *   (1) a fixed scale factor w.r.t the original images or
 *   (2) a fixed scale factor between k space and the resulting images;
 *       then a given spot is at the same position in all images.
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
 *  @author Michael Schmid, IAP/TU Wien, 2024
 */


public class LeedUndistorter extends VirtualStack implements Runnable {
    // What task to do: undostort current slice, full stack with (1) E-dependent scale, or (2) fixed reciprocal scale
    final static int SLICE=0, STACK_IN_SCALE=1, STACK_K_SCALE=2;
    // Dialog designations for the tasks
    final static String[] TASK_NAMES = new String[] {
            "Currently displayed image", "Stack, scale like input", "Stack, fixed k-space scale"};
    final static int MIN_SIZE = 40, MAX_SIZE = 2000;   //limits of output size (pixels)
    int task;           //dialog parameters used in this vartual stack
    int size;
    double scale;       //k-space value (gx, gy of the SpotPattern) per output pixel
    boolean normalizeProcI0;
    String title;       //title of the input stack (w/o extension)
    LEED_Spot_Tracker spotTracker;
    ImageStack stack;
    Roi maskRoi;
    int slice;
    LeedScreenFitter screenFitter;
    LeedIVAnalyzer ivAnalyzer;
    double[][] energiesEtc;
    ByteProcessor mask; //from the maskRoi
    Reference<float[]>[] cacheRefs; //references for caching calculated pixels

    /** Presents a dialog, undistorts the image with the given slice or the stack and shows the result.
     *  Note that 'slice' is 1-based, i.e., 1 <= slice < n.
     *  Currently, only floating-point stacks are supported.
     *  If an ivAnalyzer is supplied (non-null) and has a fit to the energy-dependent overall
     *  x/y offset, scale and angle, this correction is taken into account. */
    public void undistort(String title, LEED_Spot_Tracker spotTracker, ImageStack stack, Roi maskRoi,
            int slice, double[][] energiesEtc, LeedScreenFitter screenFitter, LeedIVAnalyzer ivAnalyzer) {
        this.title = title;
        this.spotTracker = spotTracker;
        this.stack = stack;
        this.maskRoi = maskRoi;
        this.slice = slice;
        this.energiesEtc = energiesEtc;
        this.screenFitter = screenFitter;
        this.ivAnalyzer = ivAnalyzer;
        mask = new ByteProcessor(stack.getWidth(),stack.getHeight());
		mask.setColor(255);
		mask.fill(maskRoi);
        boolean leemMode = energiesEtc[LEED_Spot_Tracker.ENERGY] == null && energiesEtc[LEED_Spot_Tracker.E_LEEM] != null;
        String[] taskNames = leemMode ?     //the last task is not available in LEEM mode
                Arrays.copyOf(TASK_NAMES, TASK_NAMES.length-1) : TASK_NAMES;
        task = (int)LeedParams.get(LeedParams.UNDISTORTTASK);
        if (task >= taskNames.length) task = taskNames.length-1;
        size = (int)LeedParams.get(LeedParams.UNDISTORTSIZE);
        normalizeProcI0 = LeedParams.getBoolean(LeedParams.UNDISTORTNORMALIZE);
        GenericDialog gd = new GenericDialog("Undistort...");
        gd.addChoice("Undistort what:", taskNames, TASK_NAMES[task]);
        gd.addNumericField("Output image size:", size, 0, 6, "pixels ("+MIN_SIZE+"-"+MAX_SIZE+")");
        if (energiesEtc[LEED_Spot_Tracker.PROCESSEDI0] != null)
            gd.addCheckbox("Normalize by processed I0", normalizeProcI0);
        if (ivAnalyzer == null)
            gd.addMessage("Stack undistortion should be run after spot tracking!");
        DialogListener dialogListener = new DialogListener() {
            public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
                double size = gd.getNextNumber();
                if (!(size>=MIN_SIZE && size<=MAX_SIZE && size==(int)size))
                    return false;
                return true;
            }
        };
        Button okButton = gd.getButtons()[0];
        okButton.setEnabled(dialogListener.dialogItemChanged(gd, null));
        gd.addDialogListener(dialogListener);
        gd.addHelp(LeedSpotTrackerHelp.getDarkFlatProcessingHelp());
        gd.showDialog();

        if (gd.wasCanceled()) return;
        task = gd.getNextChoiceIndex();
        LeedParams.set(LeedParams.UNDISTORTTASK, task);                       //save as default for the next invocation
        if (energiesEtc[LEED_Spot_Tracker.PROCESSEDI0] != null) {
            normalizeProcI0 = gd.getNextBoolean();
            LeedParams.set(LeedParams.UNDISTORTNORMALIZE, normalizeProcI0);
        }
        size = (int)gd.getNextNumber();
        if (!(size>=MIN_SIZE))
            size = MIN_SIZE;
        else if (!(size<=MAX_SIZE))
            size = MAX_SIZE;
        else
            LeedParams.set(LeedParams.UNDISTORTSIZE, size);

        // prepare for processing
        int sliceForScale = task==SLICE ? slice : stack.size();
        double energy = getEnergy(energiesEtc, sliceForScale);
        double[] xyOffsets = getXYOffsets(sliceForScale, ivAnalyzer);
        double maxKradius = maxAbsK(screenFitter, 1./Math.sqrt(energy), maskRoi, xyOffsets[0],  xyOffsets[1]);
        if (task == STACK_IN_SCALE)
            maxKradius *= 1./Math.sqrt(energy);
        scale = maxKradius/(0.5*(size - 1));
        //DEBUG IJ.log("E="+IJ.d2s(energy)+" maxKradius="+(float)maxKradius+" scale="+(float)scale);
        if (task != SLICE) {                //also check scale of first slice and make sure that both fit
            energy = getEnergy(energiesEtc, 1);
            xyOffsets = getXYOffsets(1, ivAnalyzer);
            maxKradius = maxAbsK(screenFitter, 1./Math.sqrt(energy), maskRoi, xyOffsets[0],  xyOffsets[1]);
            if (task == STACK_IN_SCALE)
                maxKradius *= 1./Math.sqrt(energy);
            double scale1 = maxKradius/(0.5*(size - 1));
            //DEBUG IJ.log("E="+IJ.d2s(energy)+" scale="+(float)scale+" scale1="+(float)scale1);
            scale = Math.max(scale, scale1);
        }
        spotTracker.enableAndHighlightComponents(false);

        Thread thread = new Thread(this, "LEED Undistorter");
        thread.start();
    }

    /** Undistortion happens in a separate thread, so we can wait until the input stack is prepared (if necessary) */
    public void run() {
        if (stack instanceof LeedDarkFlatVirtualStack) {
            if (!spotTracker.waitForStackUpdate()) {
                prepareExit();
                return;
            }
        }
        ImagePlus imp = null;
        title += "-undist";
        title = WindowManager.getUniqueName(title);
        if (task == SLICE) {
            ImageProcessor ip = getProcessor(slice);
            if (ip == null) {   // input asynchronously closed
                prepareExit();
                return;
            }
            imp = new ImagePlus("undistorted", ip);
        } else {
            cacheRefs = new Reference[stack.size()];
            imp = new ImagePlus("undistorted", this);
        }
        imp.show();
        prepareExit();
    }

    /** Things to do before exit */
    void prepareExit() {
        spotTracker.enableAndHighlightComponents(true);
    }

    /** Returns the undistorted image for the given slice */
    public ImageProcessor getProcessor(int n) {
        Reference<float[]> ref = cacheRefs == null ? null : cacheRefs[n-1];
        if (ref != null) {
            float[] cachedPixels = ref.get();
            if (cachedPixels != null)
                return new FloatProcessor(size, size, cachedPixels);
        }
        ImageProcessor ip = stack != null && stack.size() >= n ? stack.getProcessor(n) : null;
        if (ip == null)                         // asynchronously closed
            return null;
        double energy = getEnergy(energiesEtc, n);
        double invSqrtEnergy = 1./Math.sqrt(energy);
        float normalizationFactor = normalizeProcI0 && energiesEtc[LEED_Spot_Tracker.PROCESSEDI0] != null ?
                (float)(1./energiesEtc[LEED_Spot_Tracker.PROCESSEDI0][n-1]) : 1.0f;

        double xOffset = 0, yOffset = 0;        // actual (0,0) spot coordinates minus screenFitter (if available)
        double scaleCorr = 1.;
        double cos = 1, sin = 0;                // rotation matrix
        double[][] regressionResults2D = ivAnalyzer != null ? ivAnalyzer.getRegressionResults2D() : null;
        if (regressionResults2D != null && !Double.isNaN(regressionResults2D[LeedIVAnalyzer.REG2D_DELTAX][n-1])) {
            xOffset = regressionResults2D[LeedIVAnalyzer.REG2D_DELTAX][n-1];
            yOffset = regressionResults2D[LeedIVAnalyzer.REG2D_DELTAY][n-1];
            scaleCorr = 1. + 0.01*regressionResults2D[LeedIVAnalyzer.REG2D_SCALE_FIT][n-1];
            double angle = Math.toRadians(regressionResults2D[LeedIVAnalyzer.REG2D_ANGLE_FIT][n-1]);
            cos = Math.cos(angle);
            sin = Math.sin(angle);
        }
        double r0 = 0.5*(size-1);
        double scaleFactor = scale * scaleCorr;
        if (task == STACK_IN_SCALE)
            scaleFactor /= invSqrtEnergy;

        double[] xy10and01 = screenFitter.getLinearMatrix(invSqrtEnergy);
        double pxlToKscale = Math.sqrt(Math.abs(determinantOf2x2Matrix(xy10and01)));
        double downsize = pxlToKscale*scaleFactor;  // > 1 when downsizing (considering the linear transformation only)
        if (downsize > 1.05) {                  // when downsizing, smooth to avoid undersampling
            ip = ip.duplicate();
            LeedImageSmoother.smooth((FloatProcessor)ip, mask, downsize - 1);
        }

        FloatProcessor outIp = new FloatProcessor(size, size);
        double[] xy = new double[2];
        for (int y=0; y<size; y++) {
            double ky = scaleFactor * (r0 - y); // inverted sign: we want ky up; the image has y down
            for (int x=0; x<size; x++) {
                double kx = scaleFactor * (x - r0);
                if (sin == 0) {
                    xy = screenFitter.screenCoordinates(kx, ky, invSqrtEnergy, xy);
                } else {
                    double kxr = cos*kx - sin*ky;
                    double kyr = sin*kx + cos*ky;
                    xy = screenFitter.screenCoordinates(kxr, kyr, invSqrtEnergy, xy);
                }
                float pixelValue = mask.getPixel((int)Math.round(xy[0]+xOffset), (int)Math.round(xy[1]+yOffset)) != 0 ?
                        (float)ip.getInterpolatedValue(xy[0]+xOffset, xy[1]+yOffset) :
                        Float.NaN;  //also when out-of-image (then mask.getPixel returns 0)
                outIp.setf(x, y, pixelValue*normalizationFactor);
            }
        }
        outIp.setSliceNumber(n);
        if (cacheRefs != null)
            cacheRefs[n-1] = new SoftReference<float[]>((float[])outIp.getPixels());
        return outIp;
    }

    /** Returns the image width of the virtual stack */
    public int getWidth() {
        return size;
    }

    /** Returns the image height of the virtual stack */
    public int getHeight() {
        return size;
    }

    /** Returns the number of images in this virtual stack (if it is one) */
    public int getSize() {
        return stack.size();
    }

    /** Returns the label of the specified slice in this virtual stack (if it is one). */
    public String getSliceLabel(int n) {
        return stack.getSliceLabel(n);
    }

    /** Adds a slice to this virtual stack, unsupported */
    public void addSlice(String sliceLabel, Object pixels) { addSlice(""); }
    public void addSlice(String sliceLabel, ImageProcessor ip) { addSlice(""); }
    public void addSlice(String sliceLabel, ImageProcessor ip, int n) { addSlice(""); }
    public void addSlice(String fileName) {
        throw new RuntimeException("Adding Slices NOT Supported. Duplicate the stack, then try again.");
    }

    /** Deletes the specified image from this virtual stack, unsupported */
    public void deleteSlice(int n) {
        throw new RuntimeException("Deleting Slices NOT Supported. Duplicate the stack, then try again.");
    }

    
    /** Returns the maximum maximum value of kx, ky on the roi boundary.
     *  This determines the image size required for the undistorted image. */
    static double maxAbsK(LeedScreenFitter screenFitter, double invSqrtEnergy,
            Roi maskRoi, double xOffset, double yOffset) {
        double[] xy10and01 = screenFitter.getLinearMatrix(invSqrtEnergy); // dx/dkx, dy/dkx, dx/dky, dy/dky
        double[] screenToKmatrix = inverted2x2Matrix(xy10and01);
        //DEBUG IJ.log("matrix=["+IJ.d2s(xy10and01[0])+","+IJ.d2s(xy10and01[1])+"; "+IJ.d2s(xy10and01[2])+","+IJ.d2s(xy10and01[3])+"]");
        //DEBUG IJ.log("Invmat=["+IJ.d2s(screenToKmatrix[0],5)+","+IJ.d2s(screenToKmatrix[1],5)+"; "+IJ.d2s(screenToKmatrix[2],5)+","+IJ.d2s(screenToKmatrix[3],5)+"]");

        //double[] xy00 = screenFitter.screenCoordinates(0, 0, invSqrtEnergy, new double[2]);
        double x00 = screenFitter.getX00();
        double y00 = screenFitter.getY00();

        double maxRadius = 0;
        double dx = 0, dy = 0;      //as a starting position for the new point, we take the previous
        double kx = 0, ky = 0;
        double[] xy = new double[2];
        FloatPolygon fp = maskRoi.getFloatConvexHull();     //find ky, ky for all points on the convex hull of the mask
        for (int i=0; i<fp.npoints; i++) {
            double dxP = fp.xpoints[i] - (x00 + xOffset);
            double dyP = fp.ypoints[i] - (y00 + yOffset);
            for (int iter=0; iter<10; iter++) {             //iteratively refine kx, ky
                kx += screenToKmatrix[0]*(dxP - dx) + screenToKmatrix[1]*(dyP - dy);
                ky += screenToKmatrix[2]*(dxP - dx) + screenToKmatrix[3]*(dyP - dy);
                //DEBUG IJ.log("guess kx,ky="+IJ.d2s(kx,5)+","+IJ.d2s(ky,5)+" from "+IJ.d2s(dxP - dx) +","+IJ.d2s(dyP - dy));
                xy = screenFitter.screenCoordinates(kx, ky, invSqrtEnergy, xy);
                if (Double.isNaN(xy[0])) {                  //distortion at the border too bad, skip the point
                    dx = 0; dy = 0;                         //and try a fresh start
                    kx = 0; ky = 0;
                    break;
                }
                dx = xy[0] - x00;
                dy = xy[1] - y00;
                //DEBUG IJ.log("xP,yP="+IJ.d2s(dxP)+","+IJ.d2s(dyP)+" dx,dy="+IJ.d2s(dx)+","+IJ.d2s(dy)+" error x,y="+IJ.d2s(dx - dxP) +","+IJ.d2s(dy - dyP));
                double deviation = sqr(dx - dxP) + sqr(dy - dyP);
                if (deviation < 0.1) {
                    double r = Math.max(Math.abs(kx), Math.abs(ky));
                    if (r > maxRadius)
                        maxRadius = r;
                    break;
                }
            }
        }
        return maxRadius;
    }

   /** Returns the offsets in x&y in pixels, i.e., the position devviation of the (0,0) spot
     *  from the ScreenFitter model, as determined via the fit over energy.
     *  The offsets are returned as a two-element array.
     *  If the values cannot be determined, i.e., if there is no valid fit, both offsets returned are 0. */
    static double[] getXYOffsets(int slice, LeedIVAnalyzer ivAnalyzer) {
        double xOffset = 0, yOffset = 0;  // actual (0,0) spot coordinates minus screenFitter (if available)
        double[][] regressionResults2D = ivAnalyzer != null ? ivAnalyzer.getRegressionResults2D() : null;
        if (regressionResults2D != null && !Double.isNaN(regressionResults2D[LeedIVAnalyzer.REG2D_DELTAX][slice-1])) {
            xOffset = regressionResults2D[LeedIVAnalyzer.REG2D_DELTAX][slice-1];
            yOffset = regressionResults2D[LeedIVAnalyzer.REG2D_DELTAY][slice-1];
        }
        return new double[] {xOffset, yOffset};
    }

    /** Returns the energy for the given stack slice in normal LEED mode, otherwise a fixed value */
    static double getEnergy(double[][] energiesEtc, int slice) {
        double energy = energiesEtc[LEED_Spot_Tracker.ENERGY] != null ?
                energiesEtc[LEED_Spot_Tracker.ENERGY][slice-1] : LeedRadiusSelector.UNKNOWN_ENERGY;
        double invSqrtEnergy = 1./Math.sqrt(energy);
        return energy;
    }

     /** Returns the maximum of the current maximum and the square of the
     *  distance between points (x0,y0) and (x1,y1) */
    static double maxOfDistanceSqr(double currentMax, double x0, double y0, double x1, double y1) {
        double d = sqr(x1 - x0) + sqr (y1 - y0);
        return Math.max(currentMax, d);
    }

    /** Returns the inverse of a 2x2 matrix {a1, a12, a21, a22} */
    static double[] inverted2x2Matrix(double[] a) {
        double invDet = 1./determinantOf2x2Matrix(a);
        return new double[] {a[3]*invDet, -a[1]*invDet, -a[2]*invDet, a[0]*invDet};
    }

    /** Returns the determinant of a 2x2 matrix {a1, a12, a21, a22} */
    static double determinantOf2x2Matrix(double[] a) {
        return a[0]*a[3] - a[1]*a[2];
    }

    static double sqr(double x) { return x*x; }
}
