import ij.*;
import ij.process.*;
import ij.measure.*;
import ij.IJ;
import ij.util.Tools;
import ij.plugin.ZProjector;
import ij.plugin.filter.RankFilters;
import ij.plugin.filter.ThresholdToSelection;
import ij.gui.*;
import java.awt.*;
import java.util.Arrays;


/**
 *  This class helps in creating the mask of valid pixels in a LEED I(V) stack.
 *  If targetImp has a LeedDarkFlatVirtualStack, it checks for a flat field therein.
 *  If there is a flat field, it uses the average of the flat field; otherwise
 *  it calculates the pixel-by-pixel stddev of the input stack as an input.
 *  Then it tries to find an ellipse outline of the LEED screen, which limits
 *  the mask to pixels inside that ellipse (unless deselected).
 *  The mask is created by an intensity threshold of the flat field average or stddev.
 *  Unles this option is deselected, it also corrects for the radially decreasing intensity
 *  of the flat field average or stddev image.
 *  It asks the user to find a suitable threshold for the inside, and then creates a mask
 *  from this threshold. The mask may be modified by shrinking or expanding (by a few pixels)
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

public class LeedMaskCreator implements Runnable, DialogListener {
    LEED_Spot_Tracker spotTracker;
    ImagePlus inStackImp;               //the input image stack, if we don't take it from the flat of targetImp
    ImagePlus targetImp;                //here we show the mask outline
    ImageProcessor inputIp;             //the image (stddev of input or average of flat) that will can be thresholded
    ImageProcessor inputIpWithFit;      //the same, corrected for intensity decreasing with radius
    double hMin, hMax;                  //histogram min & max
    static boolean useEllipse = true;   //whether to restrict to elliptical outline, remember as default during the session
    static boolean useRadialCorr = true; //whether to fit radial decrease, remember as default during the session
    Roi ellipseOutlineRoi;
    boolean usingFlatField;
    int initialIntThreshold;

    public LeedMaskCreator(LEED_Spot_Tracker spotTracker, ImagePlus inStackImp, ImagePlus targetImp) {
        this.spotTracker = spotTracker;
        this.inStackImp = inStackImp;
        this.targetImp = targetImp;
    }

    /** Here the actual work is done*/
    public void run() {
        try {
            if (inStackImp == null || targetImp == null) return;    //should never happen
            spotTracker.enableAndHighlightComponents(false);

            IJ.showStatus("preparing for mask...");
            
            ImageStack mainStack = targetImp.getStack();
            if (mainStack instanceof LeedDarkFlatVirtualStack) {    //if possible, try to get flat-field average
                inputIp = ((LeedDarkFlatVirtualStack)mainStack).getAverageFlat();
                usingFlatField = inputIp != null;
            }
            if (inputIp == null) {
                ImagePlus stddevImp = ZProjector.run(inStackImp, "sd");
                inputIp = stddevImp.getProcessor();  //otherwide use std deviation of input stack along z
            }
            if (inputIp == null) throw new RuntimeException("Mask creation failed, could not get input");
            ellipseOutlineRoi = LeedEllipseEdgeFinder.getEllipse(inputIp);

            inputIp.resetRoi();  //(just to make sure)
            inputIp.resetMinAndMax();
            //DEBUG thImp = new ImagePlus("maskThresholdInput", inputIp); thImp.setRoi(ellipseOutlineRoi); thImp.show();
            ImageStatistics stats = inputIp.getStatistics();
            int[] histogram = stats.histogram;
            this.hMin = stats.histMin;
            this.hMax = stats.histMax;
            int nPxl = 0;
            for (int i=0; i<histogram.length; i++)
                nPxl += histogram[i];
            int minNBgPxl = (int)Math.round(0.2*nPxl);               //we want at least 20% of background area for the initial threshold
            int maxN0Pxl = (int)Math.round(0.05*nPxl);               //we take the 5% darkest pixels as zero intenisty ("darkPixelValue")
            int minIniThreshold = 0;
            double darkPixelValue = 0;
            for (int i=0, n=0; i<histogram.length/2; i++) {
                minIniThreshold = i;
                if (n <= maxN0Pxl)
                    darkPixelValue = hMin + (hMax - hMin)/(histogram.length)*i;
                if (n > minNBgPxl) break;
                n += histogram[i];
            }
            AutoThresholder thresholder = new AutoThresholder();
            int iThreshold = thresholder.getThreshold(AutoThresholder.Method.Huang, histogram);
            if (iThreshold > histogram.length/2) iThreshold = histogram.length/2;
            int iOfMax=0, histoMax=0;           //get histogram peak (of foreground pixels)
            for (int i=iThreshold; i<histogram.length; i++)
                if (histogram[i] > histoMax) {
                    histoMax = histogram[i];
                    iOfMax = i;
                }
            int iOfMin=0, histoMin=histoMax;    //get histogram valley (between foreground & background)
            for (int i=iOfMax-1; i>0; i--) {
                if (histogram[i] < histoMax/10) {
                    iOfMin = i;                 //1/10 of peak is deep enough, use this as default threshold
                    break;
                } else if (histogram[i] < histoMin) {
                    histoMin = histogram[i];
                    iOfMin = i;
                } else if (histogram[i] > histoMax || histogram[i] > 4*histoMin)
                    break;                      //looks like we have reached the background
            }
            this.initialIntThreshold = Math.max(iOfMin, minIniThreshold);
            //we use negative log for the sliders
            double logHMinOverIniTh = -Math.log(250./(this.initialIntThreshold+0.1));
            double logHMaxOverIniTh = -Math.log(0.3/(this.initialIntThreshold+0.1));
            //DEBUG IJ.log("histo: "+(float)hMin+"-"+(float)hMax+" dark="+(float)darkPixelValue+" iniThr="+(float)(hMin + (hMax-hMin)/255.*this.initialIntThreshold));
            if (ellipseOutlineRoi != null)
                inputIpWithFit = getRadialDecreaseCorrected(inputIp, ellipseOutlineRoi, darkPixelValue);

            //DEBUG IJ.log("hMin,Max="+(float)hMin+","+(float)hMax+" huang:"+iThreshold+" max@"+iOfMax+" iOfMin="+iOfMin+" slider:"+(float)logHMinOverIniTh+","+logHMaxOverIniTh);
            targetImp.getWindow().toFront();
            IJ.showProgress(1.0);

            GenericDialog gd = new NonBlockingGenericDialog(LEED_Spot_Tracker.PLUGIN_NAME+" - MaskCreation");
            if (ellipseOutlineRoi != null)
                gd.addCheckbox("Limit to elliptical fit", useEllipse);
            if (inputIpWithFit != null)
                gd.addCheckbox("Correct radius-dependent intensity", useRadialCorr);
            gd.addSlider("Mask threshold adjustment", logHMinOverIniTh, logHMaxOverIniTh, 0, 0.01);
            gd.addSlider("Shrink/grow (pxls)", -10, 10, 0, 0.5);
            gd.addMessage("Adjust the sliders or values for the best outline\n"+
                    "of the useable LEED screen area.\n"+
                    "You can modify the this outline afterwards by\n"+
                    "modifying the mask image that will be created.");

            gd.addDialogListener(this);
            dialogItemChanged(gd, null);        //show threshold with default value
            spotTracker.setCurrentFrontDialog(gd);
            gd.showDialog();

            spotTracker.setCurrentFrontDialog(null);
            if (gd.wasOKed() && targetImp.getRoi() != null) {
                ByteProcessor maskBp = targetImp.createRoiMask();
                String title = inStackImp.getTitle();
                title = LeedUtils.removeExtension(title)+"_Mask";
                title = WindowManager.getUniqueName(title);
                if (!Prefs.blackBackground) maskBp.invertLut();
                new ImagePlus(title, maskBp).show();

                if (ellipseOutlineRoi != null)
                    useEllipse = gd.getNextBoolean();
                if (inputIpWithFit != null)
                    useRadialCorr = gd.getNextBoolean();
            }
        } catch (Exception e) { IJ.handleException(e); }
        spotTracker.enableAndHighlightComponents(true);
        targetImp.resetRoi();
    }

    /** This callback method is called when the user changes fields the dialog. */
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        try {
            boolean useEllipse = ellipseOutlineRoi != null ? gd.getNextBoolean() : false;
            boolean useRadialCorr = inputIpWithFit != null ? gd.getNextBoolean() : false;
            double logThreshold = gd.getNextNumber();
            double growRadius = gd.getNextNumber();
            if (Double.isNaN(logThreshold) || Double.isNaN(growRadius)) return false;
            double intThreshold = Math.exp(-logThreshold)*(initialIntThreshold + 0.1);  //int index of histogram analyzed above
            double threshold = hMin + (hMax - hMin)/255. * intThreshold;
            //DEBUG IJ.log("intTh="+(float)intThreshold+" thr="+threshold);
            ImageProcessor inputIp = useRadialCorr ? this.inputIpWithFit : this.inputIp;
            inputIp.setThreshold(threshold, Float.MAX_VALUE, ImageProcessor.NO_LUT_UPDATE);
            int width = inputIp.getWidth();
            int height = inputIp.getHeight();
            int x=0, y=0;
            for (int i=0; i<Math.min(width/2, height/2); i++) {  //starting from center, find first thresholded pixel
                y=height/2;
                x = width/2 + i;
                if (inputIp.getPixelValue(x,y) > threshold) break;
                x = width/2 - i;
                if (inputIp.getPixelValue(x,y) > threshold) break;
                x = width/2;
                y = height/2 + i;
                if (inputIp.getPixelValue(x,y) > threshold) break;
                y = height/2 - i;
                if (inputIp.getPixelValue(x,y) > threshold) break;
            }
            Wand wand = new Wand(inputIp);
            wand.autoOutline(x, y, threshold, Float.MAX_VALUE, Wand.EIGHT_CONNECTED);
            if (wand.npoints>0) {
                Roi roi = new PolygonRoi(wand.xpoints, wand.ypoints, wand.npoints, Roi.TRACED_ROI);
                if (growRadius != 0 || useEllipse) {
                    ByteProcessor ip = new ByteProcessor(inStackImp.getWidth(), inStackImp.getHeight());
                    ip.setColor(255);
                    ip.fill(roi);
                    if (growRadius != 0)
                        (new RankFilters()).rank(ip, Math.abs(growRadius), growRadius > 0 ? RankFilters.MAX : RankFilters.MIN);
                    if (useEllipse) {
                        ip.setValue(0);
                        ip.fillOutside(ellipseOutlineRoi);
                        ip.setValue(255);
                    }
                    ip.setThreshold(254, 255, ImageProcessor.NO_LUT_UPDATE);
                    ThresholdToSelection tts = new ThresholdToSelection();
                    roi = tts.convert(ip);
                }
                targetImp.setRoi(roi);
                return true;
            } else {
                targetImp.killRoi();
                return false;
            }
        } catch(Exception ex) {IJ.handleException(ex);return false;}
    }

    /** Creates an ImageProcessor with compensation of the radial decrease
     *  of diffuse intensity typical for LEED images.
     *  The center is assumed to lie in the center of the elliptical roi.
     *  Elliptical distortion (if the ellipse is not a circle) is taken into account.
     *  The fit function for the radial decrease is
     *      y = a*exp(1-b*x*x*x/(c*c+x*x))
     *  Note that this function is always positive, so we can divide
     *  the intensity (minus the background, given as darkPixelValue) by the fit.
     *  The scale of the output is roughly preserved, so that the threshold values are still useful.
     */
    static ImageProcessor getRadialDecreaseCorrected(ImageProcessor ip, Roi ellipseOutlineRoi, double darkPixelValue) {
        if (!(ellipseOutlineRoi instanceof EllipseRoi)) return null;
        int width = ip.getWidth();
        int height = ip.getHeight();
        double[] roiParams = ((EllipseRoi)ellipseOutlineRoi).getParams();
        //IJ.log("Ellipse w,h="+width+","+height+";roiParams="+Arrays.toString(roiParams));
        if (LeedUtils.countNonNaN(roiParams) != roiParams.length) return null;
        double xC = 0.5*(roiParams[0] + roiParams[2]);
        double yC = 0.5*(roiParams[1] + roiParams[3]);
        double vLength = Math.sqrt(sqr(roiParams[0] - roiParams[2]) + sqr(roiParams[1] - roiParams[3]));
        double semiMajor = 0.5*vLength;
        double cosMajor = (roiParams[0] - roiParams[2])/vLength;
        double sinMajor = (roiParams[1] - roiParams[3])/vLength;
        double aspectRatio = roiParams[4];
        int nRadii = (int)Math.round(semiMajor);
        //for fitting, get radial increase towards the center, but without any decrease
        double[] rMajor = new double[nRadii];
        double[] yValue =  new double[nRadii];
        double maxValue = 0;
        double xCorner = 0;     //where intensities level off
        for (int i=0; i<nRadii; i++) {
            rMajor[i] = semiMajor - i;
            double scale = rMajor[i]/semiMajor;
            if (scale > 0.05) { //we don't fit the inner 5%; too few points (unreliable)
                EllipseRoi roi = new EllipseRoi(
                    xC + rMajor[i]*cosMajor, yC + rMajor[i]*sinMajor,
                    xC - rMajor[i]*cosMajor, yC - rMajor[i]*sinMajor,
                    aspectRatio);
                yValue[i] = LeedEllipseEdgeFinder.getIntensityOnRoiLine(ip, roi) - darkPixelValue;
                if (yValue[i] < maxValue) {
                    yValue[i] = maxValue;
                } else {
                    maxValue = yValue[i];
                    xCorner = rMajor[i];
                }
            } else {            //in the inner 5%, assume constant intensity
                yValue[i] = maxValue;
                if (xCorner == 0) xCorner = rMajor[i];
            }
        }
        //fit radial decrease with y = a*exp(1-b*x*x*x/(c*c+x*x))
        CurveFitter cf = new CurveFitter(rMajor, yValue);
        UserFunction userFunction = new UserFunction() {
                public double userFunction(double[] params, double x) {
                    return params[0]*Math.exp(1-params[1]*x*x*x/(params[2]*params[2]+x*x));
                }
        };
        cf.setOffsetMultiplySlopeParams(/*offsetParam=*/-1, /*multiplyParam=*/0, /*slopeParam=*/-1);
        double[] initialParams = new double[] {yValue[yValue.length-1]/Math.exp(1), 5./semiMajor, 0.2*xCorner};
        cf.doCustomFit(userFunction, /*numParams=*/3, "radial decrease: a*exp(1-b*x³/(c²+x²))", initialParams, null, /*showSettings=*/false);
        //DEBUG cf.getPlot(nRadii).show();
        //examine the fit result (we don't care about 'inaccurate result' errors)
        if (cf.getStatus() != Minimizer.SUCCESS && cf.getStatus() != Minimizer.MAX_RESTARTS_EXCEEDED) {
            IJ.log("Cannot fit radial decrease: "+cf.getStatusString());
            return null;
        }
        double[] fitParams = cf.getParams();
        double midIntensity = userFunction.userFunction(fitParams, 0.5*semiMajor);
        fitParams[0] *= 1./midIntensity;    //normalize the fit to values around 1
        //apply the fit: divide intensity minus darkPixelValue by the fit
        FloatProcessor fp = new FloatProcessor(width, height);
        for (int y=0; y<height; y++) {
            double dy = y - yC;
            for (int x=0; x<width; x++) {
                double dx = x - xC;
                double dMajor = dx*cosMajor + dy*sinMajor;
                double dMinor = dx*sinMajor - dy*cosMajor;
                double r = Math.sqrt(sqr(dMajor) + sqr(dMinor/aspectRatio));    //distance, taking elliptical distortion into account
                float value = ip.getf(x, y);
                value = (float)((value - darkPixelValue)/userFunction.userFunction(fitParams, r) + darkPixelValue);
                fp.setf(x, y, (float)value);
            }
        }
        //DEBUG new ImagePlus("fitted", fp).show();
        return fp;
    }

    static double sqr(double x) {
        return x*x;
    }
}
