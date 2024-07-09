import ij.*;
import ij.process.*;
import ij.plugin.*;
import ij.plugin.filter.*;
import ij.IJ;
import ij.util.Tools;
import ij.gui.*;
import java.awt.*;


/**
 *  This class helps in creating the mask of valid pixels in a LEED I(V) stack.
 *  It calculates the pixel-by-pixel stddev of the input stack,
 *  asks the user to find a suitable threshold for it,
 *  and then creates a mask.
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
    ImagePlus inStackImp;   //the input image stack
    ImagePlus targetImp;    //here we show the mask outline
    ImagePlus threshImp;    //this image (stddev of input) will be thresholded
    double hMin, hMax;      //histogram min & max
    int initialIntThreshold;

    public LeedMaskCreator(LEED_Spot_Tracker spotTracker, ImagePlus inStackImp, ImagePlus targetImp) {
        this.spotTracker = spotTracker;
        this.inStackImp = inStackImp;
        this.targetImp = targetImp;
    }

    /** Here the actual work is done*/
    public void run() {
        try {
            if (inStackImp == null || targetImp == null) return;       //should never happen
            spotTracker.enableAndHighlightComponents(false);

            IJ.showStatus("preparing for mask...");
            threshImp = ZProjector.run(inStackImp, "sd");
            ImageProcessor threshIp = threshImp.getProcessor(); //std deviation of input stack along z
            threshIp.resetMinAndMax();
            //DEBUG new ImagePlus("stddev", threshIp).show();
            ImageStatistics stats = threshIp.getStatistics();
            int[] histogram = stats.histogram;
            int nPxl = 0;
            for (int i=0; i<histogram.length; i++)
                nPxl += histogram[i];
            int minBgPxl = (int)Math.round(0.2*nPxl);           //we want at least 20% of background area
            int minIniThreshold = 0;
            for (int i=0, n=0; i<histogram.length/2; i++) {
                minIniThreshold = i;
                if (n > minBgPxl) break;
                n += histogram[i];
            }
            this.hMin = stats.histMin;
            this.hMax = stats.histMax;
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
            double logHMinOverIniTh = Math.log(1./(this.initialIntThreshold+0.1));
            double logHMaxOverIniTh = Math.log(255./(this.initialIntThreshold+0.1));

            //DEBUG IJ.log("hMin,Max="+(float)hMin+","+(float)hMax+" huang:"+iThreshold+" max@"+iOfMax+" iOfMin="+iOfMin+" slider:"+(float)logHMinOverIniTh+","+logHMaxOverIniTh);
            targetImp.getWindow().toFront();
            IJ.showProgress(1.0);

            GenericDialog gd = new NonBlockingGenericDialog(LEED_Spot_Tracker.PLUGIN_NAME+" - MaskCreation");
            gd.addSlider("Mask threshold adjustment", logHMinOverIniTh, logHMaxOverIniTh, 0, 0.01);
            gd.addSlider("Shrink/grow (pxls)", -10, 10, 0, 0.5);
            gd.addMessage("Adjust the sliders or values for the best outline\n"+
                    "of the useable LEED screen area.\n"+
                    "You can modify the this outline afterwards\n"+
                    "(modify the mask image that will be created).");

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
            }
        } catch (Exception e) { IJ.handleException(e); }
        spotTracker.enableAndHighlightComponents(true);
        targetImp.resetRoi();
    }

    /** This callback method is called when the user changes fields the dialog. */
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        try{
            double logThreshold = gd.getNextNumber();
            double growRadius = gd.getNextNumber();
            if (Double.isNaN(logThreshold) || Double.isNaN(growRadius)) return false;
            double intThreshold = Math.exp(-logThreshold)*(initialIntThreshold + 0.1);  //with respect to previous histogram
            double threshold = hMin + (hMax - hMin)/255. * intThreshold;
            //DEBUG IJ.log("intTh="+(float)intThreshold+" thr="+threshold);
            ImageProcessor threshIp = threshImp.getProcessor();
            threshIp.setThreshold(threshold, Float.MAX_VALUE, ImageProcessor.NO_LUT_UPDATE);
            int width = threshIp.getWidth();
            int height = threshIp.getHeight();
            int x=0, y=0;
            for (int i=0; i<Math.min(width/2, height/2); i++) {  //starting from center, find first thresholded pixel
                y=height/2;
                x = width/2 + i;
                if (threshIp.getPixelValue(x,y) > threshold) break;
                x = width/2 - i;
                if (threshIp.getPixelValue(x,y) > threshold) break;
                x = width/2;
                y = height/2 + i;
                if (threshIp.getPixelValue(x,y) > threshold) break;
                y = height/2 - i;
                if (threshIp.getPixelValue(x,y) > threshold) break;
            }
            Wand wand = new Wand(threshIp);
            wand.autoOutline(x, y, threshold, Float.MAX_VALUE);
            if (wand.npoints>0) {
                Roi roi = new PolygonRoi(wand.xpoints, wand.ypoints, wand.npoints, Roi.TRACED_ROI);
                if (growRadius != 0) {
                    ByteProcessor ip = new ByteProcessor(inStackImp.getWidth(), inStackImp.getHeight());
                    ip.setColor(255);
                    ip.fill(roi);
                    (new RankFilters()).rank(ip, Math.abs(growRadius), growRadius > 0 ? RankFilters.MAX : RankFilters.MIN);
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

}
