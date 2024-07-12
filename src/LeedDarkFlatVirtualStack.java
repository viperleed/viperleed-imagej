import ij.plugin.PlugIn;
import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.util.Tools;
import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import java.lang.ref.*;
import java.util.concurrent.atomic.*;

/** This class takes the main input for the given stack slice, corrects if with dark frame and flat fields (if given),
 *  and provides the corrected image as output. The flat field can be also processed by fitting a polynomial
 *  to its logarithm (with suitable weights) and subtracting the fit from the logarithm.
 *  Note that dark (when given) is subtracted from the source, dark2 (when given) is subtracted from the flat field.
 *  Usually, dark and dark2 will be identical.
 *  All these images must have the same size; dark & flats may be single images or stacks of the same size as the source.
 *  The output is a virtualStack, which is calculated on the fly if required.
 *  The stack is chached as cached in memory via the Java WeakReference (if fast to re-compute)
 *  or SoftReference (if more computational effort) mechanism.
 */
 /*<pre>
 Data Flow:

 mainStack -----------------\
                             > subtract ---------------------------------------------\
 darkStack  -> [averaging] -/                                                         \
                                                                                       > subtract ---- output
 flatStack -----------------\            / [averaging & fitting/normalization] -\     /
                             > subtract <                                        >---<
 dark2Stack -> [averaging] -/            \ [fitting/normalization] -------------/     \--------------- show processed flat

 [averaging]: - averaging of all slices by averageSlices method, or
              - Linear_Fit_LEED_Stack for linear fit over energy, or
              - Averaging_LEED_Stack for smoothing, or
              - no processing (input directly used)
 [fitting/normalization] (when flat has averaging, done on intermediate images, not on each image)
              - fitting of polynomial in x&y (order 2 or 4) in the logarithmic domain, inside the mask
              - normalization to average 1 inside the mask
              - no fitting/normalization
 </pre>*/

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

public class LeedDarkFlatVirtualStack extends VirtualStack {
    static final int[] POLY_ORDERS = new int[] {-1, 0, 2, 4};   //allowed polynomial orders (-1 for none)
    static final String[] POLY_NAMES = new String[] {           //corresponding names
            "None","Normalize average intensity", "2nd-order polynomial fit", "4th-order polynomial fit"};
    static final int NONE=0, AVG_ALL=1, LINEAR=2, SMOOTH=3;     //averaging types
    static final String[] AVG_NAMES = new String[] {            //corresponding names
            "No Averaging", "Average All", "Linear Fit", "Smoothing"};
    static final String[] AVG_SHORTNAMES = new String[] {       //corresponding short names
            "-", "avg", "lin", "s"};
    static final int MIN_SMOOTHING = Averaging_LEED_Stack.MIN_SMOOTHING;
    static final String PROCESSED_FLAT_INFOKEY = "ViPErLEED_processedFlat";
    ImagePlus mainImp, darkImp, flatImp, dark2Imp, maskImp;   // the input images
    ImageStack mainStack, darkStack, flatStack, dark2Stack;   // may be averaged or linear fit vs slice.
    int width, height, nSlices;         //dimensions of source input stack
    int polynomialOrder = -1;           //for fitting the flat, -1 = no fit, 0 = normalize only
    int darkAveraging, flatAveraging, dark2Averaging;
    int sliceWithBadFlat = -1;          //if a non-positive processed flat value appears, slice where it was seen
    boolean doStrongCaching;            //whether to keep a SoftReference cache for the output (otherwise WeakReference), for faster access
    Reference<float[]>[] cacheRefs;     //SoftReference (gets garbage collected only if space is needed) or WeakReference (gc'd normally) to slices
    static WeakHashMap<ImageStack,ImageStack> averagedStacks = new WeakHashMap<ImageStack,ImageStack>(); //caches averaged stacks
    static WeakHashMap<ImageStack,ImageStack> linearFitStacks = new WeakHashMap<ImageStack,ImageStack>(); //caches linear fit stacks
    Thread prefetchThread;
    Object cacheSynchronizer = new Object(); //for synchronizing the cache
    public volatile boolean isChanging; //set true externally, will be set to false after changes have been done.
    public static volatile boolean prefetchAllowed;  //SpotTracker sets this to true while it is active
    AtomicInteger nWorkingThreads = new AtomicInteger();   //tells prefetch to wait while another thread is getting data
    boolean showProgress = true;        //progress bar for stack slice averaging

    /** Creates a basic stack based on the main source stack, without any further processing. */
    public LeedDarkFlatVirtualStack(ImagePlus source) {
        this.mainImp = source;
        width = mainImp.getWidth();
        height = mainImp.getHeight();
        nSlices = mainImp.getImageStackSize();
        mainStack = mainImp.getStack();
    }

    /** Call this method first:
     *  It sets the main source image stack, the dark field applied to the source,
     *  the flat field and the dark frame ('dark2') applied to the flat field.
     *  Except for the source, the others may be single images (valid for all stack slices),
     *  stacks with a different size as the source (then they are averaged)
     *  or stacks with the same size as the source (for slice-by-slice application).
     *  The mask is used for determining the are for fitting or normalizing the flat field,
     *  it should contain a binary imaeg with foreground (pixel value 255) where the pixel values are used.
     *  Except for the source, these images may be also null; then the respective correction is not applied.
     *  Averaging: NONE=0 for no averaging, AVG_ALL=1 for averaging all, LINEAR=2 for linear over energy.
     *  TODO: larger numbers for smoothing.
     *  The polynomialOrder for fitting the flat field can be -1 to leave the flat field untouched,
     *  0 for a normalization of each flat field slice in slice-by-slice flats (same stack size as the source),
     *  2 or 4 for fitting a polynomial of the given order to the logarithm of the flat field.
     *  When done, resets 'isChanging' to false unless the Thread is interrupted.
     *  If stackImp is not null, calls update() thereafter.  */
    public void setInputDarkFlatDark2Mask(ImagePlus source, ImagePlus dark, ImagePlus flat, ImagePlus dark2,
            ImagePlus mask, int darkAveraging, int flatAveraging, int dark2Averaging, int polynomialOrder, ImagePlus stackImp) {
        //DEBUG IJ.log("setInputDarkFlatDark2Mask stack="+Integer.toHexString(this.hashCode())+" t="+(System.currentTimeMillis()-LEED_Spot_Tracker.t0));
        this.mainImp = source;
        this.darkImp  = dark;
        this.flatImp  = flat;
        this.dark2Imp = flat != null ? dark2 : null;
        width = mainImp.getWidth();
        height = mainImp.getHeight();
        nSlices = mainImp.getImageStackSize();
        if (dark != null && (dark.getWidth() != width || dark.getHeight() != height))
            throw new IllegalArgumentException("Dark & Source Image: Size Mismatch");
        if (flat != null && (flat.getWidth() != width || flat.getHeight() != height))
            throw new IllegalArgumentException("Flat & Source Image: Size Mismatch");
        if (dark2 != null && (dark2.getWidth() != width || dark2.getHeight() != height))
            throw new IllegalArgumentException("Dark2 & Source Image: Size Mismatch");

        mainStack = mainImp.getStack();
        setAveragingAndPoly(darkAveraging, flatAveraging, dark2Averaging, polynomialOrder, mask, stackImp); // also does resetCache();
        if (IJ.debugMode && darkImp != null)IJ.log("dark: "+darkImp.getTitle()+ ": "+darkStack);
        if (IJ.debugMode && flatImp != null) IJ.log("flat: "+flatImp.getTitle()+ ": "+flatStack);
        if (IJ.debugMode && dark2Imp != null) IJ.log("dark2: "+dark2Imp.getTitle()+ ": "+dark2Stack);
    }

    /** Sets the averaging of dark & flat: NONE=0 for no averaging, AVG_ALL=1 for averaging all,
     *  LINEAR=2 for linear over energy (dark, dark2 only), or larger numbers for smoothing.
     *  Also sets the polynomial order for fitting the flat and the mask for the fit.
     *  If the polynomial order is 0, the flat field (after background subtraction,
     *  if appropriate) will be normalized, i.e., divided by its average intensity
     *  within the mask. Orders 2 and 4 are also supported; with these the
     *  logarithm of the flat will be fitted with a polynomial of the given order;
     *  e.g. for 2nd order the flat field will be divided by a 2D Gaussian.
     *  Call with polynomial order = -1 to disable fitting and normalization.
     *  When done, resets 'isChanging' to false unless the Thread is interrupted.
     *  If stackImp is not null, calls update() thereafter.
     *  This method creates the stacks darkStack, flatStack, dark2Stack.  */
    public void setAveragingAndPoly(int darkAveraging, int flatAveraging, int dark2Averaging,
            int polynomialOrder, ImagePlus maskImp, ImagePlus stackImp) {
        //DEBUG IJ.log("setAveragingAndPoly stack="+Integer.toHexString(this.hashCode())+" t="+(System.currentTimeMillis()-LEED_Spot_Tracker.t0));
        this.darkAveraging = darkAveraging;
        this.flatAveraging = flatAveraging;
        this.dark2Averaging = dark2Averaging;
        this.polynomialOrder = polynomialOrder;
        int optionIndex = LeedUtils.arrayIndexOf(POLY_ORDERS, polynomialOrder);
        if (optionIndex < 0 && flatImp != null) {       //should never happen unless we get a bad macro parameter
            IJ.error("Macro parameter error or internal error:\nFlat-field fit polynomial order not supported: "+polynomialOrder);
            polynomialOrder = 0;
        };
        sliceWithBadFlat = -1;
        this.maskImp = maskImp;
        if (maskImp == null) {
            this.polynomialOrder = -1;
        }
        if (this.polynomialOrder >= 0 && (maskImp.getWidth() != width || maskImp.getHeight() != height))
            throw new IllegalArgumentException("Mask & Source Image: Size Mismatch: "+
                    maskImp.getWidth()+"x"+maskImp.getHeight()+" != "+width+"x"+height);
        if (this.polynomialOrder >= 0 && !(maskImp.getProcessor() instanceof ByteProcessor))
            throw new IllegalArgumentException("Wrong mask Type: "+maskImp.getProcessor());
        if (Thread.currentThread().isInterrupted())
            return;
        darkStack = null;
        if (darkImp != null) {                          //we need a dark frame, possibly averaged/smoothed
            if (darkImp.getImageStackSize() == 1) {
                darkStack = darkImp.getStack();
            } else {
                boolean mustAverageDark = darkImp.getImageStackSize() != nSlices;
                if (darkAveraging == AVG_ALL || mustAverageDark)
                    darkStack = averageSlices(darkImp);
                else if (darkAveraging == LINEAR) {     //pixel-wise linear fit (as a function of slice)
                    ImageStack stack = linearFitStacks.get(darkImp.getStack());
                    if (stack != null) {
                        darkStack = stack;              //use cached version if possible
                    } else {
                        darkStack = new Linear_Fit_LEED_Stack(darkImp.getStack(), null, -1, null);
                        if (darkStack.size() == nSlices && (darkStack instanceof VirtualStack))
                            linearFitStacks.put(darkImp.getStack(), darkStack);
                    }
                } else if (darkAveraging > LINEAR)      //smoothing
                    darkStack = new Averaging_LEED_Stack(darkImp.getStack(), null, darkAveraging, -1, null);
                else
                    darkStack = darkImp.getStack();
            }
        } else
            darkStack = null;
        //DEBUG IJ.log("setAveragingAndPoly got dark stack="+Integer.toHexString(this.hashCode())+" t="+(System.currentTimeMillis()-LEED_Spot_Tracker.t0));

        if (Thread.currentThread().isInterrupted())
            return;
        dark2Stack = null;
        if (dark2Imp != null && flatImp != null) {      //we need a dark2 (for flat field), possibly averaged/smoothed
            if (dark2Imp.getImageStackSize() == 1) {
                dark2Stack = dark2Imp.getStack();
            } else {
                boolean mustAverageDark2 = flatImp.getImageStackSize() != nSlices ||
                        flatAveraging == AVG_ALL || dark2Imp.getImageStackSize() != nSlices;
                int theDark2Averaging = mustAverageDark2 ? AVG_ALL : dark2Averaging; //this is what we really have to do with dark2
                if (dark2Imp == darkImp && theDark2Averaging == darkAveraging)  //we can use (averaged) dark
                    dark2Stack = darkStack;
                else if (theDark2Averaging == AVG_ALL)
                    dark2Stack = averageSlices(dark2Imp);
                else if (theDark2Averaging == LINEAR) {
                    ImageStack stack = linearFitStacks.get(dark2Imp.getStack());
                    if (stack != null) {
                        dark2Stack = stack;              //use cached version if possible
                    } else {
                        dark2Stack = new Linear_Fit_LEED_Stack(dark2Imp.getStack(), null, -1, null);
                        if (dark2Stack.size() == nSlices && (dark2Stack instanceof VirtualStack)) //(unless closed)
                            linearFitStacks.put(dark2Imp.getStack(), dark2Stack);
                    }
                } else if (theDark2Averaging > LINEAR)    //smoothing
                    dark2Stack = new Averaging_LEED_Stack(dark2Imp.getStack(), null, dark2Averaging, -1, null);
                else
                    dark2Stack = dark2Imp.getStack();
            }
        } else
            dark2Stack = null;

        if (Thread.currentThread().isInterrupted())
            return;
        if (flatImp != null) {                          //we need a flat field, possibly averaged/smoothed, possibly normalized/with fit
            if (flatImp.getImageStackSize() == 1) {
                flatStack = new LeedDifferenceStack(flatImp.getStack(), dark2Stack, this.polynomialOrder, maskImp);
            } else {
                boolean mustAverageFlat = flatImp.getImageStackSize() != nSlices;
                flatStack = flatImp.getStack();
                if (flatStack == null || flatStack.size() == 0 ||    //asynchronously closed
                        (maskImp != null && maskImp.getProcessor() == null))
                    return;
                if (flatAveraging == AVG_ALL || mustAverageFlat) {
                    ImageStack flatAvgStack = averageSlices(flatImp);
                    if (flatAvgStack != null && dark2Stack != null && dark2Stack.size() > 0)
                        flatStack = new LeedDifferenceStack(flatAvgStack, dark2Stack, this.polynomialOrder, maskImp);
                } else if (flatAveraging == LINEAR) {
                    flatStack = new Linear_Fit_LEED_Stack(flatStack, dark2Stack, this.polynomialOrder, maskImp);
                } else if (flatAveraging > LINEAR) {
                    flatStack = new Averaging_LEED_Stack(flatStack, dark2Stack, flatAveraging, this.polynomialOrder, maskImp);
                } else {
                    flatStack = new LeedDifferenceStack(flatStack, dark2Stack, this.polynomialOrder, maskImp);
                }
            }
        } else
            flatStack = null;
        if (Thread.currentThread().isInterrupted())
            return;

        sliceWithBadFlat = -1;
        doStrongCaching =   //if calculation is computationally expensive, try to cache the results as far as memory allows
                polynomialOrder > 0 || darkAveraging > LINEAR || flatAveraging >LINEAR || dark2Averaging >LINEAR ||
                (mainStack instanceof VirtualStack);
        //DEBUG IJ.log("setAveragingAndPoly pre-update stack="+Integer.toHexString(this.hashCode())+" t="+(System.currentTimeMillis()-LEED_Spot_Tracker.t0));

        if (stackImp != null)
            update(stackImp);   //includes resetCache
        else
            resetCache();       //also initializes polyCoefficients[0] to NaN, enforcing calculation
            
        if (!Thread.currentThread().isInterrupted()) { //unless we have to update again
            isChanging = false;
            IJ.showProgress(1.0);
        }
        //DEBUG IJ.log("DONE setAveragingAndPoly stack="+Integer.toHexString(this.hashCode())+" t="+(System.currentTimeMillis()-LEED_Spot_Tracker.t0));
    }

    /** Creates an ImagePlus showing this stack */
    public ImagePlus getImagePlus(String title) {
        if (mainImp == null) return null;
        ImagePlus imp = new ImagePlus(title, this);
        imp.setCalibration(mainImp.getCalibration());
        return imp;
    }

    /** Shows the dialog for dark & flat-field processing options;
     *  returns 0 if there are no processing options or the dialog was cancelled.
     *  Otherwise, the calling code should thereafter call
     *  setDarkFlatDark2AveragingFlatPoly with the values from the LeedParams.
     *  With a return code of 2, also the processed flat field should be created. */
    public int doProcessingDialog() {
        boolean canProcessDark = darkImp != null && darkImp.getImageStackSize() == nSlices;
        boolean canProcessFlat = flatImp != null;
        boolean canProcessDark2 = flatImp != null && flatImp.getImageStackSize() == nSlices &&
                dark2Imp != null && dark2Imp.getImageStackSize() == nSlices;
        if (!(canProcessDark || canProcessFlat || canProcessDark2)) return 0; //should be the same as !hasProcessingOptions() below
        boolean canAverageFlat = flatImp != null && flatImp.getImageStackSize() == nSlices;

        final HashMap<CheckboxGroup, LeedLabeledField> smoothFieldForRadioButton = new HashMap<CheckboxGroup, LeedLabeledField>(3);
        String title = "";
        if (canProcessDark || canProcessDark2) title += "Dark Frame ";
        if (canProcessFlat) title += "Flat Field ";
        GenericDialog gd = new GenericDialog(title+" Processing");
        if (canProcessDark) {
            gd.setInsets(5, 0, 0); //top left bottom
            gd.addRadioButtonGroup("DARK FRAME:", AVG_NAMES, 0, AVG_NAMES.length, AVG_NAMES[Math.min(darkAveraging, AVG_NAMES.length-1)]);
            gd.addNumericField("Smooth amount (min "+MIN_SMOOTHING+")", darkAveraging<MIN_SMOOTHING ? 50 : darkAveraging, 0);
            CheckboxGroup radioBG = (CheckboxGroup)gd.getRadioButtonGroups().lastElement();
            smoothFieldForRadioButton.put(radioBG, LeedLabeledField.numeric(gd));
        }
        int optionIndex = LeedUtils.arrayIndexOf(POLY_ORDERS, polynomialOrder);
        if (optionIndex < 0) { IJ.error("internal error: polynomial order="+polynomialOrder); optionIndex = 0; }; //should never happen
        if (canProcessFlat) {
            if (canAverageFlat) {
                gd.setInsets(canProcessDark ? 15 : 5, 0, 0); //top left bottom
                gd.addRadioButtonGroup("FLAT FIELD:", AVG_NAMES, 0, AVG_NAMES.length, AVG_NAMES[Math.min(flatAveraging, AVG_NAMES.length-1)]);
                gd.addNumericField("Smooth amount (min "+MIN_SMOOTHING+")", flatAveraging<MIN_SMOOTHING ? 50 : flatAveraging, 0);
                CheckboxGroup radioBG = (CheckboxGroup)gd.getRadioButtonGroups().lastElement();
                smoothFieldForRadioButton.put(radioBG, LeedLabeledField.numeric(gd));
            }
            gd.addChoice("Flat field processing type", POLY_NAMES, POLY_NAMES[optionIndex]);
        }
        if (canProcessDark2) {
            gd.setInsets(canProcessFlat ? 15 : 5, 0, 0); //top left bottom
            gd.addRadioButtonGroup("DARK FRAME 2 (for flat field):", AVG_NAMES, 0, AVG_NAMES.length, AVG_NAMES[Math.min(dark2Averaging, AVG_NAMES.length-1)]);
            gd.addNumericField("Smooth amount (min "+MIN_SMOOTHING+")", dark2Averaging<MIN_SMOOTHING ? 50 : dark2Averaging, 0);
            CheckboxGroup radioBG = (CheckboxGroup)gd.getRadioButtonGroups().lastElement();
            smoothFieldForRadioButton.put(radioBG, LeedLabeledField.numeric(gd));
        }
        if (canProcessFlat)
            gd.addCheckbox("Show processed flat field", false);

        DialogListener dialogListener = new DialogListener() {
            public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
                for (CheckboxGroup radioBG : smoothFieldForRadioButton.keySet()) {
                    Checkbox cbx = radioBG.getSelectedCheckbox();
                    if (cbx == null) return false;
                    boolean smoothSelected = cbx.getLabel().equals(AVG_NAMES[SMOOTH]);
                    LeedLabeledField numField = smoothFieldForRadioButton.get(radioBG);
                    numField.setEnabled(smoothSelected);
                    if (smoothSelected) {
                        String numStr = numField.getText();
                        if (!(Tools.parseDouble(numStr) >= MIN_SMOOTHING)) return false;
                    }
                }
                return true;
            }
        };
        Button okButton = gd.getButtons()[0];
        okButton.setEnabled(dialogListener.dialogItemChanged(gd, null));
        gd.addDialogListener(dialogListener);
        gd.addHelp(LeedSpotTrackerHelp.getDarkFlatProcessingHelp());

        gd.showDialog();
        if (gd.wasCanceled()) return 0;
        boolean showFlat = false;

        if (canProcessDark) {
            String optionStr = gd.getNextRadioButton();
            int darkAveraging = LeedUtils.arrayIndexOf(AVG_NAMES, optionStr);
            int smoothAmount = (int)gd.getNextNumber();
            if (darkAveraging == SMOOTH) darkAveraging = smoothAmount;
            LeedParams.set(LeedParams.DARKAVERAGING, darkAveraging);
        }
        if (canProcessFlat) {
            if (canAverageFlat) {
                String optionStr = gd.getNextRadioButton();
                int flatAveraging = LeedUtils.arrayIndexOf(AVG_NAMES, optionStr);
                int smoothAmount = (int)gd.getNextNumber();
                if (flatAveraging == SMOOTH) flatAveraging = smoothAmount;
                LeedParams.set(LeedParams.FLATAVERAGING, flatAveraging);
            }
            optionIndex = gd.getNextChoiceIndex();
            int polynomialOrder = POLY_ORDERS[optionIndex];
            LeedParams.set(LeedParams.FLATFITORDER, polynomialOrder);
        }
        if (canProcessDark2) {
            String optionStr = gd.getNextRadioButton();
            int dark2Averaging = LeedUtils.arrayIndexOf(AVG_NAMES, optionStr);
            int smoothAmount = (int)gd.getNextNumber();
            if (dark2Averaging == SMOOTH) dark2Averaging = smoothAmount;
            LeedParams.set(LeedParams.DARK2AVERAGING, dark2Averaging);
        }
        if (canProcessFlat)
            showFlat = gd.getNextBoolean();
        return showFlat ? 2 : 1;
    }

    /** Returns whether we have dark/flat/dark2 options to ask for in doProcessingDialog() */
    public boolean hasProcessingOptions() {
        boolean canProcessDark = darkImp != null && darkImp.getImageStackSize() == nSlices;
        boolean canProcessFlat = flatImp != null;
        boolean canProcessDark2 = flatImp != null && flatImp.getImageStackSize() == nSlices &&
                dark2Imp != null && dark2Imp.getImageStackSize() == nSlices;
        return canProcessDark || canProcessFlat || canProcessDark2;
    }

    /** Returns a reference to the source ImaegPlus */
    public ImagePlus getMainImp() {
        return mainImp;
    }

    /** Returns a FloatProcessor for the specified slice of this virtual stack
     *  where 1<=n<=nslices. Returns null on error. */
    public ImageProcessor getProcessor(int n) {
        float[] pixels = getPixels(n);
        if (pixels == null) return new FloatProcessor(width, height); //dummy if interrupted
        ImageProcessor ip = new FloatProcessor(width, height, pixels);
        ip.setSliceNumber(n);
        return ip;
    }

    /** Returns the pixels array for the specified slice of this virtual stack
     *  where 1<=n<=nslices. Returns null on error. */
    public float[] getPixels(int n) {
        Reference[] cacheRefs = this.cacheRefs;
        Reference<float[]> ref = cacheRefs == null ? null : cacheRefs[n];
        if (ref != null) {
            float[] cachedPixels = ref.get();
            if (cachedPixels != null) return cachedPixels;
        }
        ImageStack mainStack = this.mainStack; // work on local copies, because these can change asynchronously
        ImageStack darkStack = this.darkStack;
        ImageStack flatStack = this.flatStack;
        //DEBUG IJ.log("synchronized got stacks stack="+Integer.toHexString(this.hashCode())+" t="+(System.currentTimeMillis()-LEED_Spot_Tracker.t0));
        if (Thread.currentThread().isInterrupted()) return zeroArray();
        nWorkingThreads.incrementAndGet();
        ImageProcessor sourceIp = mainStack.getProcessor(n);
        if (Thread.currentThread().isInterrupted() || sourceIp == null)
            return zeroArray();
        ImageProcessor darkIp = null;
        if (darkStack != null && darkStack.size() > 0) {
            darkIp = darkStack.getProcessor(darkStack.size() == 1 ? 1 : n);
            if (Thread.currentThread().isInterrupted() || darkIp == null)
            return zeroArray();
        }
        float[] flatPixels = null;
        if (flatStack != null && flatStack.size() > 0) {
            Object fPixels = flatStack.getPixels(flatStack.size() == 1 ? 1 : n);
            flatPixels = fPixels instanceof float[] ?
                    (float[])fPixels :
                    LeedUtils.toFloat(fPixels);
            if (Thread.currentThread().isInterrupted() || flatPixels == null)
                return zeroArray();
        }

        float[] pixels = null;
        if (flatPixels != null || darkIp != null) {
            pixels = new float[width*height];
            ImageProcessor maskIp = LeedUtils.getMaskIp(maskImp);
            byte[] maskPixels = maskIp==null ? null : (byte[])maskIp.getPixels();
            for (int p=0; p<pixels.length; p++) {
                float result = sourceIp.getf(p);
                if (darkIp != null) result -= darkIp.getf(p);
                if (flatPixels != null) {
                    if (!(flatPixels[p] > 0) && (maskPixels==null || maskPixels[p] != 0))
                        sliceWithBadFlat = n;
                    result /= flatPixels[p];    //up to here, 'pixels' holds the corrected flat
                    if (!(result >= -(float)1e30 && result < (float)1e30)) result = 0f;
                }
                pixels[p] = result;
            }
            if (pixels != null && cacheRefs != null)
                cacheRefs[n] = doStrongCaching ? new SoftReference<float[]>(pixels) : new WeakReference<float[]>(pixels);
        } else {
            pixels = (float[])sourceIp.convertToFloat().getPixels();
        }

        nWorkingThreads.decrementAndGet();
        return pixels;
    }

    /** Returns an ImagePlus with the processed flat field. */
    public ImagePlus getProcessedFlat() {
        ImageStack flatStack = this.flatStack;          //work on local copies, to avoid asynchronous changes
        if (flatStack == null) return null;
        ImagePlus flatImp = this.flatImp;
        if (flatImp == null) return null;
        String title = flatImp.getTitle()+"_processed_"+polynomialOrder;
        String flatInfo = flatImp.getInfoProperty();   //transfer 'info' property from original flat field
        if (flatInfo == null) flatInfo = "";
        flatInfo += PROCESSED_FLAT_INFOKEY + ": polynomialOrder="+polynomialOrder+'\n';  //retrieve with imp.getProp(PROCESSED_FLAT_INFOKEY)

        ImagePlus outImp = new ImagePlus(title, flatStack);
        outImp.setProperty("Info", flatInfo);
        return outImp;
    }

    /** If interrupted, returns this array to avoid a NullPointerException */
    private float[] zeroArray() {
            return new float[width*height];
    }

    /** Returns the size of the flat field (stack or single image) in bytes */
    public long getFlatSizeInBytes() {
        if (flatStack == null)
            return 0;
        else
            return width * height * 4L * flatStack.size() ;
    }

    /** Returns a String describing processing of the dark frames & flat field */
    public String getFlatProcessingString() {
        if (!hasProcessingOptions())
            return "";
        String str = "";
        if (darkImp != null && darkStack != darkImp.getStack()) {
            if (darkStack != null && darkStack.size() == nSlices) {
                str += "dark:"+AVG_SHORTNAMES[Math.min(darkAveraging, AVG_SHORTNAMES.length-1)];
                if (darkAveraging>=SMOOTH) str+=IJ.d2s(darkAveraging,0);
            } else str += "dark:avg";
        }
        if (flatImp != null) {
            if (str.length() > 0) str += ' ';
            str += "flat:";
            if (flatAveraging != NONE && flatImp.getImageStackSize() == nSlices && flatImp.getStack() != flatStack) {
                str += AVG_SHORTNAMES[Math.min(flatAveraging, AVG_SHORTNAMES.length-1)];
                if (flatAveraging>=SMOOTH) str+=IJ.d2s(flatAveraging,0);
                str += ' ';
            }
            if (polynomialOrder < 0) {
                boolean isProcessedFlat = flatImp.getProp(PROCESSED_FLAT_INFOKEY) != null;
                if (flatImp.getImageStackSize() > 1)
                    str+= isProcessedFlat ? "not norm." : "NOT normalized";
            } else if (polynomialOrder == 0)
                str+= "norm.";
            else //if (polynomialOrder > 0)
                str+= "fit"+polynomialOrder;

            if (flatStack != null && dark2Stack != null && dark2Imp != null && dark2Stack != dark2Imp.getStack()) {
                if (str.length() > 0) str += ' ';
                if (dark2Stack.size() == nSlices && flatStack.size() == nSlices) {
                    str += "dark2:"+AVG_SHORTNAMES[Math.min(dark2Averaging, AVG_SHORTNAMES.length-1)];
                    if (dark2Averaging>=SMOOTH) str+=IJ.d2s(dark2Averaging,0);
                } else
                    str += "dark2:avg";
            }
        }
        return str;
    }

    /** Returns the image width of the virtual stack */
    public int getWidth() {
        return mainImp.getWidth();
    }

    /** Returns the image height of the virtual stack */
    public int getHeight() {
        return mainImp.getHeight();
    }

    /** Returns the number of images in this virtual stack (if it is one) */
    public int getSize() {
        return mainImp.getImageStackSize();
    }

    /** Returns the label of the specified slice in this virtual stack (if it is one). */
    public String getSliceLabel(int n) {
        return mainStack.getSliceLabel(n);
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

    /** Updates the display of the ImagePlus showing the stack. Also resets the cache. */
    public void update(final ImagePlus imp) {
        if (imp == null) return;
        if ((mainStack != null && mainStack.size() == 0) ||
                (darkStack != null && darkStack.size() == 0) ||
                (flatStack != null && flatStack.size() == 0) ||
                (dark2Stack != null && dark2Stack.size() == 0))
            return;         //if an input stack closed asynchronously
        resetCache();
        int slice = imp.getCurrentSlice();
        ImageProcessor impIp = imp.getProcessor();
        if (impIp == null)
            return;         //if closed asynchronously
        float[] imagePixels = (float[])impIp.getPixels();
        float[] pixels = (float[])getProcessor(slice).getPixels();
        if (pixels == null || imagePixels == null)
            return;         //if closed asynchronously
        System.arraycopy(pixels, 0, imagePixels, 0, pixels.length);
        EventQueue.invokeLater(new Runnable() { //imp.updateAndDraw() is synchronized; avoid potential deadlock with zSelector
                public void run() {
                    imp.updateAndDraw();
                }});
    }

    /** If a non-positive flat was encountered returns its slice number, otherwise -1 */
    public int getSliceWithBadFlat() {
        return sliceWithBadFlat;
    }

    /** Resets the cache of the processed images, to enforce recalculation.
     *  Does not recalculate the currently displayed image; use update() for this.
     *  This method may also start a low-priority thread for prefetch (filling the cache) in the background.
     *  Synchronization since it may be called from different threads */
    @SuppressWarnings("unchecked")
    public void resetCache() {
        synchronized(cacheSynchronizer) {
            if (cacheRefs == null)
                cacheRefs = new Reference[nSlices+1];
            else
                Arrays.fill(cacheRefs, null);
            if (prefetchThread != null) terminatePrefetch();
            if (doStrongCaching && prefetchAllowed) {
                prefetchThread = new Thread(new Runnable() {
                    public void run() {                     //Prefetch all slices of output 'Spot tracking' stack
                        try {
                            Thread.sleep(10000);             //start prefetch with 10 sec delay
                        } catch (InterruptedException e) { return; }
                        try {
                            long lastTime = System.currentTimeMillis();
                            for (int n=1; n<=nSlices; n++) {
                                if (Thread.currentThread().isInterrupted()) return;
                                if (!prefetchAllowed) return;
                                long currentMemory = IJ.currentMemory();
                                long maxMemory = IJ.maxMemory();
                                long freeMemory = maxMemory - currentMemory;
                                long wantedMemory = width * height * 4L *200;   //space for >=200 frames wanted
                                int retry = 0;
                                while (freeMemory < wantedMemory) {
                                    if (retry > 3) return;  //probably we don't have enough memory to keep everything in the cache
                                    System.gc();
                                    IJ.wait(1000);
                                    retry++;
                                }
                                long time = System.currentTimeMillis();
                                if (time - lastTime > 1000) //after a second of work, yield to other processes that want to access the data
                                    try {                   //(the low priority of this thread does not apply to I/O, only to CPU)
                                        Thread.sleep(500);
                                        lastTime = System.currentTimeMillis();
                                    } catch (InterruptedException e) { return; }
                                while (nWorkingThreads.get() > 0)   //wait until no one else wants data
                                    try {
                                        Thread.sleep(500);
                                    } catch (InterruptedException e) { return; }
                                Reference[] cr = cacheRefs; //may become null asynchronously
                                if (cr != null && cr.length>n && cr[n] == null) {
                                    Object o = getPixels(n);//fill the cache
                                    if (o == null)          //error, e.g. input closed asynchronously
                                        return;
                                }
                            }
                        }
                        catch (Exception e) {IJ.handleException(e);}
                    }
                }, "LEED_Spot_Tracker_Prefetch");
                prefetchThread.setPriority(Math.max(Thread.MIN_PRIORITY, Thread.currentThread().getPriority()-3));
                prefetchThread.start();
            }
        }
    }

    /** Terminates prefetch. Call only if this stack will be discarded.
     *  Blocks until prefetch has ended, to avoid a previous prefetch
     *  overwriting the results of a prefetch thread started later. */
    public void terminatePrefetch() {
        Thread prefetchThread = this.prefetchThread;    //safe if it changes asynchronously
        if (prefetchThread == null) return;
        try {
            prefetchThread.interrupt();
            if (prefetchThread.isAlive())
                prefetchThread.join();
        } catch (Exception e) {}
    }

    /** Normalizes the pixels such that the average inside the mask is 1 */
    public static void normalize(float[] pixels, ImagePlus maskImp) {
        byte[] maskPixels = LeedUtils.getMaskPixels(maskImp);
        float factor = getFlatNormalizationFactor(pixels, maskPixels);
        for (int p=0; p<pixels.length; p++)
            pixels[p] *= factor;
    }

    /** Returns the normalization factor for a float pixels array
     *  inside the foreground area of the mask, also given by a pixels array */
    static float getFlatNormalizationFactor(float[] flatPixel, byte[] maskPixel) {
        double sum = 0;
        int n = 0;
        for (int i=0; i<flatPixel.length; i++)
            if (maskPixel[i] != 0 && flatPixel[i] > 0) {
                sum += flatPixel[i];
                n++;
            }
        return (float)(sum/n);
    }

    /** Determines whether progress bars will be shown by the substacks that support it
     *  (currently only relevant for substacks of type Averaging_LEED_Stack, which may take
     *  some time for the calculation) */
    public void setShowProgress(boolean b) {
        for (Object stack : new Object[] {mainStack, darkStack, flatStack, dark2Stack}) {
            if (stack instanceof Averaging_LEED_Stack)
                ((Averaging_LEED_Stack)stack).setShowProgress(b);
        }
        this.showProgress = b;
    }

    /** Returns the average of the flat field stack, corrected by the average of dark2, when present.
     *  Returns null if there is no flat or an input is anynchronously closed. */
    public ImageProcessor getAverageFlat() {
        if (flatImp == null) return null;
        ImageStack flatAvgStack = averageSlices(flatImp);
        if (flatAvgStack == null) return null;
        ImageStack dark2AvgStack = null;
        if (dark2Imp != null)
            dark2AvgStack = averageSlices(dark2Imp);
        ImageProcessor flatAvgIp = flatAvgStack.size() == 1 ?
                flatAvgStack.getProcessor(1) : null;
        if (flatAvgIp == null || flatAvgIp.getPixels() == null) return null;
        ImageProcessor dark2AvgIp = dark2AvgStack != null && dark2AvgStack.size() == 1 ?
                dark2AvgStack.getProcessor(1) : null;
        if (dark2AvgIp != null) {
            ImageStack differenceStack = new LeedDifferenceStack(flatAvgStack, dark2AvgStack,
                    /*polynomialOrder=*/0, /*maskImp=*/null);
            flatAvgIp = differenceStack.getProcessor(1);
        }
        return flatAvgIp;
    }

    /** Returns an ImageStack with averaged slices of the input.
     *  Returns null if interrupted or an input is asynchronously closed. */
    ImageStack averageSlices(ImagePlus imp) {
        int width = imp.getWidth();
        int height = imp.getHeight();
        ImageStack inStack = imp.getStack();
        int nSlices = inStack == null ? 0 : inStack.size();
        if (nSlices == 0)
            return null;
        else if (nSlices == 1)
            return inStack;
        ImageStack cachedStack = averagedStacks.get(inStack);
        if (cachedStack != null)
            return cachedStack;
        FloatProcessor outIp = new FloatProcessor(width, height);
        float[] outPixels = (float[])outIp.getPixels();
        for (int n=1; n<=inStack.size(); n++) {
            if (Thread.currentThread().isInterrupted()) return null;
            if (showProgress)           //accessing slices of a virtual stack can be slow
                IJ.showProgress((n-1)/(double)nSlices);
            ImageProcessor inIp = inStack.getProcessor(n);
            for (int i=0; i<width*height; i++)
                outPixels[i] += inIp.getf(i);
        }
        for (int i=0; i<width*height; i++)
            outPixels[i] *= 1.0/nSlices;
        ImageStack outStack = new ImageStack(width, height, 1);
        outStack.setPixels(outPixels, 1);
        outStack.setSliceLabel("AVG_"+imp.getTitle(), 1);
        if (showProgress)
            IJ.showProgress(1.0);
        if (inStack.size() > 1 && (inStack instanceof VirtualStack))  //unless the input stack has been closed
            averagedStacks.put(inStack, outStack);
        return outStack;
    }
}




