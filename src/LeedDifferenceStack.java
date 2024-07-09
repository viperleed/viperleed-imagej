import ij.*;
import ij.process.*;
import ij.IJ;
import ij.util.Tools;
import java.lang.ref.*;


/**
 *  This class is a virtual stack providing the difference of two stacks, which can either
 *  have the same size or the second stack can have a size of 1 whereas the first one does not.
 *  It also provides normalization or polynomial fit for the result.
 *  This class does only weak caching; only the fit coefficients are kept.
 *  Stronger chaching (SoftReference) is done by the LeedDarkFlatVirtualStack.
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
 *  @author Michael Schmid, IAP/TU Wien, 2020-2024
 */


public class LeedDifferenceStack extends VirtualStack {
    ImageStack inStack, darkStack;
    int width, height, nSlices;
    float[] oneSlicePixels;         //keeps the result if we have only one slice
    ImagePlus maskImp;              //determines the area of polynomial fitting
    int polynomialOrder;            //for polynomial fitting (no fit if <=0)
    double[] xyFitScale;            //for polynomial fitting, see LeedDarkFlatVirtualStack
    double[][] polyCoefficients;    //first index is order in x&y, second index the slice (0 unused, n>0 if a stack, n=1 if flat is not a stack)
    Reference<float[]>[] cacheRefs; //keeps non-garbage-collected output

    /** Creates a new VirtualStack with the (possibly processed) difference of the two stacks.
     *  @param inStack The input image stack
     *  @param darkStack Dark frame images. When non-null, will be subtracted form inStack.
     *          This stack must have either exactly one slice or exactly the same number of slices as inStack.
     *  @param smoothing: Determines the amout of smoothing. Noise suppression is
     *          roughly comparable to unweighted averaging over 'smoothing' adjacent slices.
     *          The actual effective kernel is larger than this by about 1/3.
     *  @param polynomialOrder Order of a polynomial. When 0, normalization to an average of 1.0 inside
     *          the mask area is done. When >0, the polynomial fit is done in the mask area
     *          in the logarithmic domain (for details see LeedDarkFlatVirtualStack).
     *          No fitting is done when polynomialOrder < 0.
     *  @param maskImp For polynomial fitting, defines the area of fitting (must be a binary 8-bit image) */
    public LeedDifferenceStack(ImageStack inStack, ImageStack darkStack, int polynomialOrder, ImagePlus maskImp) {
        this.inStack = inStack;
        this.darkStack = darkStack;
        this.polynomialOrder = polynomialOrder;
        this.maskImp = maskImp;
        this.nSlices = darkStack == null ?
                inStack.size() :
                Math.max(inStack.size(), darkStack.size());
        this.width = inStack.getWidth();
        this.height = inStack.getHeight();
        if (maskImp == null)
            this.polynomialOrder = -1;
        if (this.polynomialOrder > 0) {
            xyFitScale = LeedFlatFitter.getXYFitScale(LeedUtils.getMaskIp(maskImp));
            polyCoefficients = LeedFlatFitter.getNewPolyCoeffs(this.polynomialOrder, nSlices);
        }
        cacheRefs = new Reference[nSlices];
    }

    public ImageProcessor getProcessor(int n) {
        ImageProcessor ip = new FloatProcessor(width, height, getPixels(n));
        ip.setSliceNumber(n);
        return ip;
    }

    /** Returns the pixels for slice n (1 <= n <= nSlices). Do not modify the output array. */
    @SuppressWarnings("unchecked")
    public float[] getPixels(int n) {
        if (n<1 || n>nSlices) throw new IllegalArgumentException("Slice index "+n+" out of range 1-"+nSlices);
        if (nSlices == 1 && oneSlicePixels != null)
            return oneSlicePixels;
        Reference<float[]> ref = cacheRefs[n-1];
        if (ref != null) {              //maybe we still have a previous array for this slice
            float[] cachedPixels = (float[])ref.get();
            if (cachedPixels != null) return cachedPixels;
        }

        ImageProcessor inIp = inStack == null || inStack.size() == 0 ?
                null : inStack.getProcessor(inStack.size() == 1 ? 1 : n);
        if (inIp == null) return null;
        ImageProcessor darkIp = darkStack == null || darkStack.size() == 0 ?
                null : darkStack.getProcessor(darkStack.size() == 1 ? 1 : n);
        float[] pixels = new float[width*height];
        byte[] maskPixels = polynomialOrder == 0 ? LeedUtils.getMaskPixels(maskImp) : null;
        double sum = 0;
        int count = 0;
        for (int p=0; p<pixels.length; p++) {
            float value = inIp.getf(p);
            if (darkIp != null) value -= darkIp.getf(p);
            pixels[p] = value;
            if (polynomialOrder == 0 && maskPixels[p] != 0 && value > 0) {
                sum += value;
                count++;
            }
        }
        if (polynomialOrder == 0)       //normalize
            for (int p=0; p<pixels.length; p++)
                pixels[p] *= (float)(count/sum);
        else if (polynomialOrder > 0)   //polynomial correction, calculates polyCoefficients if required
            LeedFlatFitter.correctFlatWithPoly(pixels, maskImp, xyFitScale, polynomialOrder, polyCoefficients[n-1]);

        if (nSlices == 1)
            oneSlicePixels = pixels;    //if a single frame: remember the result
        else
            cacheRefs[n-1] = new WeakReference(pixels);
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
