import ij.*;
import ij.process.*;
import ij.plugin.*;
import ij.plugin.filter.*;
import ij.plugin.filter.ThresholdToSelection;
import ij.IJ;
import ij.util.Tools;
import ij.gui.*;
import java.awt.*;
import java.util.*;
import java.lang.ref.*;


/**
 *  This class provides a stack smoothed in z direction.
 *  It starts with averaging groups of adjacent slices and then runs a triangular
 *  smoothing kernel acting on four adjacent group averages. In contrast to applying
 *  a kernel directly to the data, the pre-averaging into groups reduces both memory
 *  and computing time (especially when a polynomial fit in x&y is applied to the
 *  group averages).
 *  The disadvantage of this approach: The center of peaks in the data is not
 *  accurately preserved.
 *
 *  Noise suppression for white noise by convolution with the triangular kernel
 *  is between sqrt(5/16) and sqrt(3/8) (since the kernel also interpolates
 *  between the group averages, the value depends on the interpolation position).
 *  The average of the squared noise gain (over all interpolation positions) is 1/3.
 *  Including averaging within groups, noise gain is reduced by a further factor of sqrt(1/n),
 *  where n is the number of data points (stack slices) per group.
 *
 *  "Amount of smoothing": This parameter ('smoothing') determines the amout of smoothing.
 *  Noise suppression is roughly comparable to unweighted averaging over 'smoothing'
 *  adjacent slices. The actual effective kernel is about 4/3 times larger, however.
 *  Note that the 'smoothing' parameter is not exact, mainly due to less smoothing at the
 *  ends, but also as the noise suppression depends on the position of the output slice
 *  with respect to the group it belongs to (see noise gain of the kernel, above).
 *
 *  Handling of near-end slices: We use the same kernel, but truncated where
 *  it falls outside the region of input data. This is a compromise between
 *  preserving linear trends like a 1st-order Savitzky-Golay-like filter (which
 *  has strongly increasing noise towards the end) and mofifying the kernel
 *  to keep the noise suppresion the same at the ends (at the cost of strongly
 *  flattening out the curve of pixel value vs slice number near the end).
 *
 *  The input can be also the difference of two stacks, inStack minus DarkStack.
 *  For application to flat field image stacks, polynomial fitting (inside a mask)
 *  as in the LeedDarkFlatVirtualStack class is also supported.
 *  In this case, the polynomial fit is applied to the group averages (when iterating
 *  thorugh all stack slices, this is less computational effort than fitting each
 *  output slice).
 *
 *  Caching: This class also works if is not possible to hold all group averages in
 *  memory. The group averages are held via Java SoftReferences, i.e., they are deleted
 *  only if memory is required and the garbage collector cannot free any other memory.
 *  In case of insufficient memory, the preformance will be poor, however, due to
 *  frequent garbage collection as well as recalculating group averages (especially when
 *  iterating through the stack more than once and polynomial fitting is used).
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

public class Averaging_LEED_Stack extends VirtualStack implements PlugInFilter {
    public static final int MIN_SMOOTHING = 10;
    public static final int K_SIZE = 4; //we average over 4 groups
    final double GROUP_SIZE_PER_SMOOTH = 1./3.; //roughly the average kernel noise gain
    static double lastSmoothing = 50;  //default for use as PlugInFilter
    ImagePlus imp;              //for use as PlugInFilter only
    ImageStack inStack;         //the main input
    ImageStack darkStack;       //must be subtracted from main input, may be null
    int width, height, nSlices;
    int polynomialOrder = -1;   //for polynomial fitting
    ImagePlus maskImp;          //determines the area of polynomial fitting
    double[] xyFitScale;        //for polynomial fitting, see LeedDarkFlatVirtualStack
    double[][] polyCoefficients;//remember fitting coefficients, in case the processed group averages get gc'd
    double kernelLength;        //in units of original slices
    Reference[] groupPxlAvgs;   //holds the averages over groups
    int[] groupStartSliceIs;    //for each group, first slice number contained minus 1 (slices start at 1)
    int[] groupEndSliceIs;      //for each group, last slice number contained minus 1 + 1
    Reference<float[]>[] cacheRefs; //keeps non-garbage-collected output
    boolean showProgress = true;//set to false to fully suppress updating the ImageJ progress bar

    /** Creator for use by Java code.
     *  @param inStack The input image stack
     *  @param darkStack Dark frame images. When non-null, will be subtracted form inStack.
     *          This stack must have either exactly one slice or exactly the same number of slices as inStack.
     *  @param smoothing: Determines the amount of smoothing. Noise suppression is
     *          roughly comparable to unweighted averaging over 'smoothing' adjacent slices.
     *          The actual effective kernel is larger than this by about 1/3.
     *  @param polynomialOrder Order of a polynomial. When 0, normalization to an average of 1.0 inside
     *          the mask area is done. When >0, the polynomial fit is done in the mask area
     *          in the logarithmic domain (for details see LeedDarkFlatVirtualStack).
     *          No fitting is done when polynomialOrder < 0.
     *  @param maskImp For polynomial fitting, defines the area of fitting (must be a binary 8-bit image) */
    public Averaging_LEED_Stack(ImageStack inStack, ImageStack darkStack, double smoothing, int polynomialOrder, ImagePlus maskImp) {
        this.imp = imp;
        this.inStack = inStack;
        this.darkStack = darkStack;
        this.polynomialOrder = polynomialOrder;
        this.maskImp = maskImp;
        setup(smoothing);
    }

    /** No-argument creator required for an ImageJ PlugInFilter */
    public Averaging_LEED_Stack() {}

    /** When called as a PlugInFilter, this method is called first */
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        return DOES_8G | DOES_16 | DOES_32 | STACK_REQUIRED | NO_CHANGES;
    }

    /** When called as a PlugInFilter, this method is called second */
    public void run(ImageProcessor ip) {
        double smoothing = IJ.getNumber("Amount of smoothing (min "+MIN_SMOOTHING+"):", lastSmoothing);
        if (smoothing == IJ.CANCELED) return;
        if (!(smoothing >= MIN_SMOOTHING)) {
            IJ.error ("invalid smoothing: "+IJ.d2s(smoothing, 1)+"<"+MIN_SMOOTHING);
            return;
        }
        this.inStack = imp.getStack();
        setup(smoothing);
        String title = WindowManager.getUniqueName("AVG"+IJ.d2s(smoothing, 0)+'_'+imp.getTitle());
        getImagePlus(title).show();
    }

    /** Creates an ImagePlus showing this stack */
    public ImagePlus getImagePlus(String title) {
        ImagePlus outImp = new ImagePlus(title, this);
        if (imp != null) outImp.setCalibration(imp.getCalibration());
        return outImp;
    }

    /** Determines whether to use the ImageJ progress bar */
    public void setShowProgress(boolean b) {
        showProgress = b;
    }

    /** Preparation only; all computationally expensive calculations are done in getPixels().
     *  @param smoothing Determines the amout of smoothing. Noise suppression is
     *          roughly comparable to unweighted averaging over 'smoothing' adjacent slices.
     *          The actual effective kernel is larger than this by about 1/3. */
    @SuppressWarnings("unchecked")
    void setup (double smoothing) {
        this.nSlices = inStack.size();
        if (nSlices == 0) return;   //asynchronously closed
        this.width = inStack.getWidth();
        this.height = inStack.getHeight();
        double avgGroupSize = smoothing*GROUP_SIZE_PER_SMOOTH;
        int nGroups = (int)Math.round(nSlices/avgGroupSize);
        if (nGroups < 1) nGroups = 1;
        groupPxlAvgs = new SoftReference[nGroups];
        groupStartSliceIs = new int[nGroups];
        groupEndSliceIs = new int[nGroups];
        kernelLength = K_SIZE*avgGroupSize;
        int nextStart = 0;
        for (int ig=0; ig<nGroups; ig++) {
            groupStartSliceIs[ig] = nextStart;
            groupEndSliceIs[ig] = (int)Math.round((ig+1.0)/nGroups*nSlices+1e-10);
            double fourGroupCenterSpan = ig < K_SIZE ? Double.MAX_VALUE : groupCenter(ig) - groupCenter(ig-K_SIZE);
            if (kernelLength > fourGroupCenterSpan)
                kernelLength = fourGroupCenterSpan;
            nextStart = groupEndSliceIs[ig];
            //IJ.log(ig+": from "+groupStartSliceIs[ig]+" to "+groupEndSliceIs[ig]+" 4gSpan="+fourGroupCenterSpan+" kLen="+kernelLength);
        }
        if (maskImp != null && polynomialOrder > 0) {
            xyFitScale = LeedFlatFitter.getXYFitScale(LeedUtils.getMaskIp(maskImp)); //(null if asynchronously closed)
            polyCoefficients = LeedFlatFitter.getNewPolyCoeffs(this.polynomialOrder, nGroups);
        }
        cacheRefs = new Reference[nSlices];
    }

    public ImageProcessor getProcessor(int n) {
        ImageProcessor ip = new FloatProcessor(width, height, getPixels(n));
        ip.setSliceNumber(n);
        return ip;
    }

    /** Returns the pixels for slice n (1 <= n <= nSlices), as obtained by
     *  smoothing in z and normalization/averaging if requested.
     *  Do not modify the output array.
     *  This call may be interrupted and then returns an array with zero values. */
    public float[] getPixels(int n) {
        if (n<1 || n>nSlices) throw new IllegalArgumentException("Slice index "+n+" out of range 1-"+nSlices);
        Reference<float[]> ref = cacheRefs[n-1];
        if (ref != null) {                  //maybe we still have a previous array for this slice
            float[] cachedPixels = (float[])ref.get();
            if (cachedPixels != null) return cachedPixels;
        }
        float[] pixels = new float[width*height];
        int iSlice = n-1;
        int nGroups = groupPxlAvgs.length;
        int iGroup = (nGroups*iSlice)/nSlices-K_SIZE/2-1;   //first group that may be relevant, slightly too low (to be on the safe side)
        if (iGroup < 0) iGroup = 0;
        int firstGroupInside = -1;
        float[] kernel = new float[K_SIZE];
        double sumKernel = 0;
        for (int ik=0; ik<K_SIZE && iGroup<nGroups; ) {
            double dx = (groupCenter(iGroup) - iSlice)/kernelLength;  //distance from kernel center
            if (dx <= (-0.5 + 1e-10)) {     //we need not process this group (it is below or at the kernel start):
                iGroup++;                   //  continue without incrementing kernel index ik
                continue;
            } else if (dx >= (0.5 - 1e-10)) //this group is already at or beyond the end of kernel -- done
                break;
            if (Thread.currentThread().isInterrupted())
                return pixels;
            kernel[ik] = (float)(0.5-Math.abs(dx)) * (groupEndSliceIs[iGroup] - groupStartSliceIs[iGroup]); //groups with more slices get more weight
            sumKernel += kernel[ik];
            if (firstGroupInside < 0) firstGroupInside = iGroup;
            ik++; iGroup++;
        }
        for (int ik=0; ik<K_SIZE; ik++)
            kernel[ik] *= (float)(1./sumKernel);
        if (showProgress && !Thread.currentThread().getName().contains("refetch"))
            IJ.showProgress(0.0);
        for (int ik=0, ig=firstGroupInside; ik<K_SIZE && kernel[ik]!=0; ik++, ig++) {
            float[] sumPxls = getGroupPxlSum(ig);
            if (showProgress && !Thread.currentThread().getName().contains("refetch"))
                IJ.showProgress((ik+0.9)*(1.0/K_SIZE));
            for (int p=0; p<pixels.length; p++)
                pixels[p] += sumPxls[p]*kernel[ik];
        }
        cacheRefs[n-1] = new WeakReference(pixels);
        if (showProgress && !Thread.currentThread().getName().contains("refetch"))
            IJ.showProgress(1.0);
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

    /** central index (slice number-1) of group ig */
    double groupCenter(int ig) {
        return 0.5*(groupStartSliceIs[ig]+groupEndSliceIs[ig]-1);
    }

    /** Returns the pixel sums (summed over the slices, for each x, y) of the i-th group. */
    @SuppressWarnings("unchecked")
    float[] getGroupPxlSum(int ig) {
        Reference ref = groupPxlAvgs[ig];
        if (ref != null) {
            float[] pxlSums = (float[])ref.get();
            if (pxlSums != null )           //group pixel sums calculated already and still available
                return pxlSums;
        }
        int darkSize = darkStack == null ? 0 : darkStack.size();

        float[] pxlSums = new float[width*height];
        for (int is=groupStartSliceIs[ig]; is<groupEndSliceIs[ig]; is++) {
            ImageProcessor darkIp = darkStack == null ? null : darkStack.getProcessor(darkSize == 1 ? 1 : is+1);
            Object pixels = inStack.getPixels(is+1);
            if (pixels instanceof byte[])
                for (int p=0; p<pxlSums.length; p++) {
                    float v = ((byte[])pixels)[p] & 0xff;
                    if (darkIp != null) v -= darkIp.getf(p);
                    pxlSums[p] += v;
                }
            else if (pixels instanceof short[])
                for (int p=0; p<pxlSums.length; p++) {
                    float v = ((short[])pixels)[p] & 0xffff;
                    if (darkIp != null) v -= darkIp.getf(p);
                    pxlSums[p] += v;
                }
            if (pixels instanceof float[])
                for (int p=0; p<pxlSums.length; p++) {
                    float v = ((float[])pixels)[p];
                    if (darkIp != null) v -= darkIp.getf(p);
                    pxlSums[p] += v;
                }
        }
        for (int p=0; p<pxlSums.length; p++)            //sum to average
            pxlSums[p] *= (float)(1./(groupEndSliceIs[ig]-groupStartSliceIs[ig]));

        if (maskImp != null) {
            if (polynomialOrder == 0 || xyFitScale == null)  //normalization
                LeedFlatFitter.normalize(pxlSums, maskImp);
            else if (polynomialOrder > 0)              //polynomial fit
                LeedFlatFitter.correctFlatWithPoly(pxlSums, maskImp, xyFitScale, polynomialOrder, polyCoefficients[ig]);
        }
        groupPxlAvgs[ig] = new SoftReference(pxlSums);  //keep the result (at least as long we have enough memory)
        return pxlSums;
    }
}
