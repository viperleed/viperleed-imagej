import ij.*;
import ij.process.*;
import ij.plugin.*;
import ij.plugin.filter.*;
import ij.plugin.filter.ThresholdToSelection;
import ij.IJ;
import ij.util.Tools;
import ij.gui.*;
import java.awt.*;


/**
 *  This class fits a linear function of z at each pixel position, and
 *  provides a stack of the images resulting from the fit.
 *  This is useful for noise suppression of dark-frame stacks with slow and
 *  linear energy dependence of the brightness, as it can happen for
 *  (voltage-dependent) glow due to field emission.
 *  The input can be also the difference of two stacks, inStack minus DarkStack.
 *  For application to flat field image stacks, polynomial fitting (inside a mask)
 *  as in the LeedDarkFlatVirtualStack class is also supported.
 *  In this case, the polynomial fit is applied to the two images that would be
 *  returned for the first and last slice, and the images in between are interpolated.
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
 *  @author Michael Schmid, IAP/TU Wien, 2021-2024
 */


public class Linear_Fit_LEED_Stack extends VirtualStack implements PlugInFilter {
    ImageStack inStack, darkStack;
    ImagePlus imp;              //for use as PlugInFilter only
    int polynomialOrder;        //for polynomial fitting (no fit if <=0)
    ImagePlus maskImp;          //determines the area of polynomial fitting
    float[] pixels1, pixelsN;
    int width, height, nSlices;
    boolean showProgress = true;

    /** Creator for use by Java code.
     *  @param inStack The input image stack
     *  @param darkStack Dark frame images. When non-null, will be subtracted form inStack.
     *          This stack must have either exactly one slice or exactly the same number of slices as inStack.
     *  @param polynomialOrder Order of a polynomial. When 0, normalization to an average of 1.0 inside
     *          the mask area is done. When >0, the polynomial fit is done in the mask area
     *          in the logarithmic domain (for details see LeedDarkFlatVirtualStack).
     *          No fitting is done when polynomialOrder < 0.
     *  @param maskImp For polynomial fitting, defines the area of fitting (must be a binary 8-bit image) */
    public Linear_Fit_LEED_Stack(ImageStack inStack, ImageStack darkStack, int polynomialOrder, ImagePlus maskImp) {
        this.inStack = inStack;
        this.darkStack = darkStack;
        this.polynomialOrder = polynomialOrder;
        this.maskImp = maskImp;
        setup();
    }

    /** No-argument creator required for an ImageJ PlugInFilter */
    public Linear_Fit_LEED_Stack() {}

    /** When called as a PlugInFilter, this method is called first */
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        return DOES_8G | DOES_16 | DOES_32 | STACK_REQUIRED | NO_CHANGES;
    }

    /** When called as a PlugInFilter, this method is called second */
    public void run(ImageProcessor ip) {
        this.inStack = imp.getStack();
        setup();
        String title = WindowManager.getUniqueName("FIT_"+imp.getTitle());
        getImagePlus(title).show();
    }

    /** Creates an ImagePlus showing this stack */
    public ImagePlus getImagePlus(String title) {
        ImagePlus outImp = new ImagePlus(title, this);
        if (imp != null) outImp.setCalibration(imp.getCalibration());
        return outImp;
    }

    /** Here we do the linear fit; and the polynomial fit if required */
    void setup () {
        this.nSlices = inStack.size();
        this.width = inStack.getWidth();
        this.height = inStack.getHeight();
        int darkSize = darkStack == null ? 0 : darkStack.size();
        int size = this.width*this.height;
        double[] sumY = new double[size];
        double[] sumXY = new double[size];
        double sumX=0, sumX2=0;
        long t0 = System.currentTimeMillis() + 400;
        boolean progressShown = false;
        for (int i=0; i<nSlices; i++) {
            Object oPixels = inStack.getPixels(i+1);  //stack slices start with 1
            ImageProcessor darkIp = darkStack == null ? null : darkStack.getProcessor(darkSize == 1 ? 1 : i+1);
            if (oPixels instanceof byte[]) {
                byte[] pixels = (byte[])oPixels;
                for (int p=0; p<size; p++) {
                    double v = pixels[p]&0xff;
                    if (darkIp != null) v -= darkIp.getf(p);
                    sumY[p] += v;
                    sumXY[p] += v*i;
                }
            } else  if (oPixels instanceof short[]) {
                short[] pixels = (short[])oPixels;
                for (int p=0; p<size; p++) {
                    double v = pixels[p]&0xffff;
                    if (darkIp != null) v -= darkIp.getf(p);
                    sumY[p] += v;
                    sumXY[p] += v*i;
                }
            } else  if (oPixels instanceof float[]) {
                float[] pixels = (float[])oPixels;
                for (int p=0; p<size; p++) {
                    double v = pixels[p];
                    if (darkIp != null) v -= darkIp.getf(p);
                    sumY[p] += v;
                    sumXY[p] += v*i;
                }
            }
            sumX += i;
            sumX2 += i*(double)i;
            long t = System.currentTimeMillis();
            if (t > t0 + 100) {
                IJ.showProgress(i+1,nSlices);
                t0 = t;
                progressShown = true;
            }
        }
        this.pixels1 = new float[size];
        this.pixelsN = new float[size];
        for (int p=0; p<size; p++) {
            double slope = nSlices > 1 ? (sumXY[p]-sumY[p]*(sumX/nSlices))*(1./(sumX2-sumX*(sumX/nSlices))) : 0;
            double offset = (sumY[p] - slope*sumX)*(1./nSlices);
            this.pixels1[p] = (float)offset;
            this.pixelsN[p] = (float)(offset + (nSlices-1)*slope);
        }
        if (maskImp != null) {
            if (polynomialOrder == 0) {             //normalization to average=1.0 inside mask area
                LeedFlatFitter.normalize(this.pixels1, maskImp);
                LeedFlatFitter.normalize(this.pixelsN, maskImp);
            } else if (polynomialOrder > 0) {       //polynomial fit inside mask area
                double[] xyFitScale = LeedFlatFitter.getXYFitScale(LeedUtils.getMaskIp(maskImp));
                LeedFlatFitter.correctFlatWithPoly(this.pixels1, maskImp, xyFitScale, polynomialOrder, null);
                LeedFlatFitter.correctFlatWithPoly(this.pixelsN, maskImp, xyFitScale, polynomialOrder, null);
            }
        }
        if (progressShown)
            IJ.showProgress(1.0);
    }

    public ImageProcessor getProcessor(int n) {
        return new FloatProcessor(width, height, getPixels(n));
    }

    /** Returns the pixels for slice n, from linear interpolation between the ends.
     *  The output array is always newly created and may therefore be modified. */
    public float[] getPixels(int n) {
        float[] pixels = new float[width*height];
        float weight1 = (float)((nSlices-n)/(double)(nSlices-1));
        float weightN = (float)((n-1)/(double)(nSlices-1));
        for (int p=0; p<pixels.length; p++)
            pixels[p] = pixels1[p]*weight1 + pixelsN[p]*weightN;
        return pixels;
    }

    /** Returns the image width of the virtual stack */
    public int getWidth() {
        return width;
    }

    /** Returns the image height of the virtual stack */
    public int getHeight() {
        return height;
    }

    /** Returns the number of images in this virtual stack */
    public int getSize() {
        return nSlices;
    }

    /** Returns the label of the specified slice in this virtual stack (if it is one). */
    public String getSliceLabel(int n) {
        return inStack.getSliceLabel(n);
    }

    /** Adds a slice to this virtual stack, unsupported */
    public void addSlice(String sliceLabel, Object pixels) { addSlice(""); }
    public void addSlice(String sliceLabel, ImageProcessor ip) { addSlice(""); }
    public void addSlice(String sliceLabel, ImageProcessor ip, int n) { addSlice(""); }
    public void addSlice(String fileName) {
        throw new RuntimeException("Adding Slices NOT Supported");
    }

    /** Deletes the specified image from this virtual stack, unsupported */
    public void deleteSlice(int n) {
        throw new RuntimeException("Deleting Slices NOT Supported");
    }

}
