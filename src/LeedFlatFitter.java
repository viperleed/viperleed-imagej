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

/** This class contains static methods for normalization and fitting the flat field.
 *  The class is also used for fitting a nonlinear background
 *  Normalization and fitting only takes these pixels into account where the corresponding
 *  ask pixel is nonzero.
 *  Polynomial fitting (order 2 or 4) is done not with the original data but with the logarithm
 *  of the data (with suitable weights) and the fit is subtracted from the logarithm.
 *  (or, equivalently, the value is divided by exp(fit_value). Note that a 2nd-order polynomial
 *  in the logarithmic domain then corresponds to a 2D Gaussian.
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


public class LeedFlatFitter {
    public static final int X0=0, Y0=1, XFACT=2, YFACT=3; //indices in xyFitScale array for x0, y0, xFact, yFact

    /** Normalizes the pixels such that the average inside the mask is 1.0 */
    public static void normalize(float[] pixels, ImagePlus maskImp) {
        byte[] maskPixels = (byte[])LeedUtils.getMaskIp(maskImp).getPixels();
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

    /** Returns the scaling coefficients for the polynomial as a 4-element array.
     *  The polynomial uses scaled coordinates x, y with a mean around 0 and typical absolute order of magnitude one,
     *  to improve stability of the fitting. The variables of the polynomial are (x-x0)*xFact and (y-y0)*yFact,
     *  where x and y are the pixel coordinates. */
    public static double[] getXYFitScale(ImageProcessor maskIp) {
        if (maskIp == null)
            return null;            //if asynchronously closed
        byte[] maskPixels = (byte[])maskIp.getPixels();
        if (maskPixels == null)
            return null;
        double sumX=0, sumX2=0, sumY=0, sumY2=0;
        int count=0;
        int width = maskIp.getWidth();
        int height = maskIp.getHeight();
        for (int y=0, p=0; y<height; y++)
            for (int x=0; x<width; x++,p++)
                if (maskPixels[p] != 0) {
                    sumX += x;  sumX2 += x*x;
                    sumY += y;  sumY2 += y*y;
                    count++;
                }
        double[] xyFitScale = new double[4];
        double x0 = sumX/count;
        double y0 = sumY/count;
        xyFitScale[X0] = x0;
        xyFitScale[Y0] = y0;
        xyFitScale[XFACT] = 1./Math.sqrt(sumX2/count - x0*x0); //scale factor is 1/stddev
        xyFitScale[YFACT] = 1./Math.sqrt(sumY2/count - y0*y0);
        return xyFitScale;
    }

    /** Calculates the polynomial-corrected flat field for slice n. The (dark-field corrected) pixels of the
     *  flat field must be supplied and will be modified.
     *  Calculates the polynomial fit coefficients if required and enters them in the coefficients array
     *  Note that weights of the polynomial fit to log(pixel value) change from slice to slice, thus the matrix
     *  must be calculated and inverted for each slice separately (for the same reason, also an approach with
     *  orthogonal basis functions common to all slices would not work)
     *  */
    public static void correctFlatWithPoly(float[] flatPixels, ImagePlus maskImp,
            double[] xyFitScale, int polynomialOrder, double[] coefficients) {
        int width = maskImp.getWidth();
        int height = maskImp.getHeight();
        double x0 = xyFitScale[X0];
        double y0 = xyFitScale[Y0];
        double xFact = xyFitScale[XFACT];
        double yFact = xyFitScale[YFACT];
        ImageProcessor maskIp = LeedUtils.getMaskIp(maskImp);
        byte[] maskPixels = maskIp == null ? null : (byte[])maskIp.getPixels();
        if (maskPixels == null) return; //(may have become null asynchronously, e.g. on 'Close All')
        if (coefficients == null) {
            coefficients = new double[nPolyCoeffs(polynomialOrder)];
            coefficients[0] = Double.NaN;
        }
        if (Double.isNaN(coefficients[0])) {  //we have to calculate the polynomial fit first
            boolean fitOk = calculateCoefficients(coefficients, width, height, flatPixels, maskPixels, xyFitScale, polynomialOrder);
            if (!fitOk) {
                LeedUtils.logError("Flat field fit (matrix inversion) failed");
                if (Double.isNaN(coefficients[0])) LeedUtils.logError("Flat (minus dark') has non-positive pixel values!");
            }
        }
        int nCoeff = coefficients.length;
        //String s="coeff:"; for (int c=0; c<nCoeff; c++) s+=IJ.d2s(coefficients[c],-3)+", "; IJ.log(s);
        for (int y=0, p=0; y<height; y++) {  //apply the polynomial
            double yP = (y-y0)*yFact;
            for (int x=0; x<width; x++,p++) {
                if (maskPixels[p] != 0) {
                    double xP = (x-x0)*xFact;
                    double fitValue = 0;
                    double polyTermY = 1;
                    for (int yOrder = 0, c=0; yOrder<=polynomialOrder; yOrder++, polyTermY *= yP) {
                        double polyTerm = polyTermY;
                        for (int xOrder = 0; xOrder<=polynomialOrder-yOrder; xOrder++, c++, polyTerm *= xP)
                            fitValue += polyTerm*coefficients[c];    //term with x^xOrder * y^yOrder
                    }
                    flatPixels[p] *= (float)Math.exp(-fitValue);            //the equivalent of subtracting fitValue from the logarithm
                } else
                    flatPixels[p] = 1f;
            } // for x
        } //for y
    }

    /** Calculates the fit coefficients and enters them in the array provided.
     *  Otherwise sets the first value to the average and leaves the others untouched. */
    static boolean calculateCoefficients(double[] coefficients, int width, int height,
            float[] fPixels, byte[] maskPixels, double[] xyFitScale, int polynomialOrder) {
        double x0 = xyFitScale[X0];
        double y0 = xyFitScale[Y0];
        double xFact = xyFitScale[XFACT];
        double yFact = xyFitScale[YFACT];
        int nCoeff = nPolyCoeffs(polynomialOrder);
        double[][] matrix = new double[nCoeff][nCoeff];
        double[] vector = new double[nCoeff];
        double[] polyTerms = new double[nCoeff];
        for (int y=0, p=0; y<height; y++) {
            double yP = (y-y0)*yFact;
            for (int x=0; x<width; x++,p++) {
                if (maskPixels != null && maskPixels[p] == 0) continue;
                double pixelValue = fPixels[p];
                if (pixelValue <= 0) continue;              //since we fit log(pixelValue)
                double xP = (x-x0)*xFact;
                double polyTermY = 1;
                for (int yOrder = 0, c=0; yOrder<=polynomialOrder; yOrder++, polyTermY *= yP) {
                    double polyTerm = polyTermY;
                    for (int xOrder = 0; xOrder<=polynomialOrder-yOrder; xOrder++, c++, polyTerm *= xP)
                        polyTerms[c] = polyTerm;            //this is x^xOrder * y^yOrder
                }
                double value = Math.log(pixelValue);        //we fit log(pixelValue) with weights pixelValue
                double weight = pixelValue;
                for (int c=0; c<nCoeff; c++) {
                    matrix[c][c] += polyTerms[c]*(polyTerms[c]*weight);
                    for (int c1=c+1; c1<nCoeff; c1++) {     //off-diagonal terms
                        matrix[c1][c] += polyTerms[c1]*(polyTerms[c]*weight);
                        matrix[c][c1] += polyTerms[c1]*(polyTerms[c]*weight);
                    }
                    vector[c] += value*(polyTerms[c]*weight);
                }
            } //for x
        } //for y
        double sumWeights = matrix[0][0];
        //for (int r=0; r<nCoeff; r++) {String s=""; for (int c=0; c<nCoeff; c++) s+=IJ.d2s(matrix[r][c], -3)+", "; IJ.log(s);}
        boolean matrixOk = LeedUtils.invertSymmetricMatrix(matrix);
        Arrays.fill(coefficients, 0.0);
        if (matrixOk) {
            for (int c=0; c<nCoeff; c++)
                for (int c1=0; c1<nCoeff; c1++)
                    coefficients[c1] += matrix[c][c1] * vector[c];
        } else          //matrix inversion failed, use normalization (i.e., only 0-th order polynomial)
            coefficients[0] = vector[0]/sumWeights;

        return matrixOk;
    }

    /** Creates an image corresponding to the polynomial fit in the log domain. For debug purposes only. */
    public static FloatProcessor createFitImage(ImageProcessor maskIp,
            double[] xyFitScale, int polynomialOrder, double[] coefficients) {
        int width = maskIp.getWidth();
        int height = maskIp.getHeight();
        double x0 = xyFitScale[X0];
        double y0 = xyFitScale[Y0];
        double xFact = xyFitScale[XFACT];
        double yFact = xyFitScale[YFACT];
        byte[] maskPixels = maskIp == null ? null : (byte[])maskIp.getPixels();
        int nCoeff = coefficients.length;
        FloatProcessor outFp = new FloatProcessor(width, height);
        float[] outPixels = (float[])outFp.getPixels();
        for (int y=0, p=0; y<height; y++) {  //use the polynomial
            double yP = (y-y0)*yFact;
            for (int x=0; x<width; x++,p++) {
                if (maskPixels == null || maskPixels[p] != 0) {
                    double xP = (x-x0)*xFact;
                    double fitValue = 0;
                    double polyTermY = 1;
                    for (int yOrder = 0, c=0; yOrder<=polynomialOrder; yOrder++, polyTermY *= yP) {
                        double polyTerm = polyTermY;
                        for (int xOrder = 0; xOrder<=polynomialOrder-yOrder; xOrder++, c++, polyTerm *= xP)
                            fitValue += polyTerm*coefficients[c];    //term with x^xOrder * y^yOrder
                    }
                    outPixels[p] = (float)Math.exp(fitValue);
                }
            } // for x
        } //for y
        return outFp;
    }

    /** Returns the number of fit coefficients required for a polynomial order */
    public static int nPolyCoeffs(int polynomialOrder) {
        return (polynomialOrder + 1)*(polynomialOrder + 2)/2;
    }

    /** Returns an array of n sets of polynomial fit coefficients, ready for use (set to NaN) */
    public static double[][] getNewPolyCoeffs(int polynomialOrder, int n) {
        double[][] polyCoefficients = new double[n][nPolyCoeffs(polynomialOrder)];
        for (int i=0; i<polyCoefficients.length; i++)
            polyCoefficients[i][0] = Double.NaN;
        return polyCoefficients;
    }
}




