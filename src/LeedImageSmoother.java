import ij.plugin.filter.PlugInFilter;
import ij.*;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.process.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

/** This class smooths an image with a (roughly) triangular kernel in x&y,
 *  as a preprocessing step before sampling with a spacing that would otherwise
 *  lead to undersampling. Used by the LeedUndistorter.
 *  If a mask is supplied, image pixels are considered only the corresponding
 *  mask pixel is nonzero. Output pixels are NaN if there are no input data to average.
 *  The input must be a floating-point image and a binary mask of the same size.
 *  This class is meant as utility class; the user interface (implementing
 *  PlugInFilter) is for testing only.
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

public class LeedImageSmoother implements PlugInFilter {
    private final static String[] TYPES = new String[] {
        "mean", "border-limited mean", "median", "minimum", "maximum", "eliminate maxima", "eliminate minima",
        "background from minima", "background from maxima", "background from median"
    };
    private static double radius = 2;           // The kernel radius, may be non-integer
    private ByteProcessor maskIp;
    private final int flags = DOES_ALL|CONVERT_TO_FLOAT|SNAPSHOT;
    /**
     * When used with user interface, this method is called by ImageJ for initialization.
     * @param arg Unused here. For plugins in a .jar file this argument string can
     *            be specified in the plugins.config file of the .jar archive.
     * @param imp The ImagePlus containing the image (or stack) to process.
     * @return    The method returns flags (i.e., a bit mask) specifying the
     *            capabilities (supported formats, etc.) and needs of the filter.
     *            See PlugInFilter.java and ExtendedPlugInFilter in the ImageJ
     *            sources for details.
     */
    public int setup(String arg, ImagePlus imp) {
        if (IJ.versionLessThan("1.54a"))        // generates an error message for older versions
            return DONE;
        GenericDialog gd = new GenericDialog("Triangular-Kernel Smooth...");
        gd.addNumericField("Radius", radius, 1);
        gd.showDialog();                        // display the dialog; preview runs in the background now
        if (gd.wasCanceled()) return DONE;
        radius = gd.getNextNumber();
        Roi roi = imp.getRoi();
        if (roi != null && roi.isArea()) {      //create mask from current Roi (if any)
            maskIp = new ByteProcessor(imp.getWidth(),imp.getHeight());
            maskIp.setColor(255);
            maskIp.fill(roi);
        }
        return flags;
    }

    /** This method is called by ImageJ when used with a user interface (after 'setup') */
    public void run(ImageProcessor ip) {
        smooth((FloatProcessor)ip, maskIp, radius);
    }
    
    /** Filters a float image fp. The mask determines which pixels are considered valid input.
     *  All pixels are considered valid if the mask is null.
     *  For integer values of the kernel radius, the kernel is a triangular function with
     *  the given radius, i.e. decaying to zero at the given radius+1.
     *  Thus, no smoothing is performed with radius <=0.
     *  When using as preprocessing step for downsampling, the radius should be set to 1/scale - 1
     *  where 'scale' is the factor between the size of the downsampled image and the original image. */
    public static void smooth(FloatProcessor ip, ByteProcessor maskIp, double radius) {
        if (radius <= 0) return;
        int width = ip.getWidth();
        int height = ip.getHeight();
        float[] fPixels = (float[])ip.getPixels();
        byte[] mPixels = maskIp == null ? null : (byte[])maskIp.getPixels();
        float[] fLine = new float[Math.max(width, height)];     //image orw or column
        double[] cache = new double[Math.max(width, height)];     //work array
        for (int y=0, p=0; y<height; y++, p+=width)
            filterLine(fPixels, mPixels, cache, p, 1, width, radius);
        for (int x=0, p=0; x<width; x++, p++)
            filterLine(fPixels, mPixels, cache, p, width, height, radius);
    }

    /** Filter one row (inc=1) or column (inc=width).
     *  The first pixel to handle is given by 'pointer' and the number of pixels
     *  in the row or column ny nData.
     *  Filtering with a triangular kernel is achieved by filtering twice with a rectangular kernel.
     *  Strictly speaking, the kernel is triangular (from 2x rectangular) for integer radius only. */
    private static void filterLine(float[] fPixels, byte[]mPixels, double[] cache, int pointer, int inc, int nData, double radius) {
        int r1 = (int)Math.floor(radius+1);   //filtering for arbitrary radius is a weighted sum of two filters with integer radius
        int r2 = r1 + 1;
        int dp1 = r1*inc, dp2 = r2*inc;     //offset in arrays to previous point that we have to subtract
        double weight2 = radius+1 - r1;
        double weight1 = 1 - weight2;
        { //first pass with a rectangular kernel, smooth towards the right
            double sum1 = 0, sum2 = 0;      //sums for the two rectangular filters
            int n1 = 0, n2 = 0;             //number of pixels in the sum
            int validFrom = 0;              //valid data start at i=0 (unless out-of-mask)
            for (int i=0, p = pointer; i<nData; i++, p+=inc) {
                boolean isInside = mPixels==null || mPixels[p]!=0;
                if (isInside) {
                    sum1 += fPixels[p];
                    sum2 += fPixels[p];
                    n1++;
                    n2++;
                }
                if (i-r1 >=validFrom && (mPixels==null || mPixels[p-dp1]!=0)) {
                    sum1 -= fPixels[p-dp1];
                    n1--;
                }
                if (i-r2 >=validFrom && (mPixels==null || mPixels[p-dp2]!=0)) {
                    sum2 -= fPixels[p-dp2];
                    n2--;
                }
                if (!isInside) {
                    sum1 = 0; sum2 = 0;
                    n1 = 0; n2 = 0;
                    validFrom = i + 1;      //at best, valid data start with the next point
                }
                cache[i] = weight1*sum1/n1 + weight2*sum2/n2;
            }
        }
        { //second pass with a rectangular kernel, smooth towards the left
            double sum1 = 0, sum2 = 0;      //sums for the two rectangular filters
            int n1 = 0, n2 = 0;             //number of pixels in the sum
            int validEnd = nData;           //valid data end where the array ends (or earlier if restricted by the mask)
            for (int i=nData-1, p = pointer + (nData-1)*inc; i>=0; i--, p-=inc) {
                boolean isInside = mPixels==null || mPixels[p]!=0;
                if (isInside) {
                    sum1 += cache[i];
                    sum2 += cache[i];
                    n1++;
                    n2++;
                }
                if (i+r1 < validEnd && (mPixels==null || mPixels[p+dp1]!=0)) {
                    sum1 -= cache[i+r1];
                    n1--;
                }
                if (i+r2 < validEnd && (mPixels==null || mPixels[p+dp2]!=0)) {
                    sum2 -= cache[i+r2];
                    n2--;
                }
                if (!isInside) {
                    sum1 = 0; sum2 = 0;
                    n1 = 0; n2 = 0;
                    validEnd = i;           //no valid data
                }
                fPixels[p] = (float)(weight1*sum1/n1 + weight2*sum2/n2);
            }
        }
    }
}
