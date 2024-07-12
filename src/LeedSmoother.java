import ij.IJ;
import java.util.Arrays;

/**
 *  Part of the LEED I/V package
 *
 *  Simple implementation of the Modified Sinc (MS) smoother of degree
 *  4 as described in
 *  M. Schmid, D. Rath and U. Diebold,
 *  'Why and how Savitzky-Golay filters should be replaced', in
 *  ACS Measurement Science Au 2, 185 (2022).
 *  doi:10.1021/acsmeasuresciau.1c00054
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

public class LeedSmoother {
    /** Smalles kernel halfwidth that can be used for smoothing. Everything below leads to no smoothing. */
    public static final int MIN_HALFWIDTH = 4;
    /** This class implements only smoothing degree 4 */
    private static final int DEGREE = 4;
    /** The kernel for filtering */
    private double[] kernel;
    /** The weights for linear fitting for extending the data at the boundaries */
    private double[] fitWeights;

    /** Constructor for testing only. The result should be
     *  [-2.358052962436155E-17, 0.15190721432800858, 0.0010247965627336176, -0.0566986446979822,
     *    0.005934455804735752, 0.21808020781696172, -0.0448991641611013, -1.0164302357649069,
     *    0.1724817848366889, 5.227319949237206, 9.045810323567304, 7.258346561237802,
     *    3.162495697496988, 0.7964266384865009, -0.13024084634616687]
     */
    /*public LeedSmoother() {
        int m = 7;
        kernel = makeKernel(m);
        fitWeights = makeFitWeights(m);
        double[] data = new double[] //arbitrary test data
                {0, 1, -2, 3, -4, 5, -6, 7, -8, 9, 10, 6, 3, 1, 0};
        double[] out = smooth(data);
        IJ.log("Filtered data:");
        IJ.log(Arrays.toString(out));
    }/**/


    /** Constructor, creates the kernel and the weights for near-boundary points,
     *  as used for linear extrapolation at the boundaries.
     *  Note that no smoothing will be performed if the kernelHalfwidth m is less then MIN_HALFWIDTH = 4. */
    public LeedSmoother(int kernelHalfwidth) {
        if (kernelHalfwidth < MIN_HALFWIDTH) return;
        kernel = makeKernel(kernelHalfwidth);
        fitWeights = makeFitWeights(kernelHalfwidth);
    }

    /**
     * Calculates the kernel halfwidth ("radius") m best suited for obtaining a given noise gain.
     * @param invNoiseGainSqr The square of the reciprocal value of the noise gain
     *         for white noise. (The RMS value of white noise would be multiplied by
     *         the noise gain.) invNoiseGainSqr is also the width of a moving-averages
     *         filter with equal white noise gain. If invNoiseGainSqr is n points,
     *         white noise will be reduced by a factor of 1/sqrt(n).
     * @return The kernel halfwidth m required.
     */
    public static int invNoiseGainSquareToM(double invNoiseGainSqr) {
        if (invNoiseGainSqr < 1.3) return 0; //no smoothing; 1.42 corresponds to the minimum m=4
        final double exponent = -2.5-0.8*DEGREE;
        final double a = 1.494 + 0.4965*DEGREE;
        final double b = 0.52;
        double m = -1 + a*invNoiseGainSqr +
                     b*Math.pow(invNoiseGainSqr, exponent);
        return (int)Math.round(m);
    }

    /** Calculates the square of the reciprocal white noise gain for a given kernel halfwidth ("radius") m.
     *  The result of this function can be also interpreted as the number of points in a
     *  moving-average filter with the same suppression of white noise. */
    public static double mToInvNoiseGainSquare(int m) {
        if (m < MIN_HALFWIDTH) return 0;
        final double exponent = -2.5-0.8*DEGREE;
        final double a = 1.494 + 0.4965*DEGREE;
        final double b = 0.52;
        double ings = (m + 1)*(1.0/a); //first approximation
        ings = (m + 1 - b*Math.pow(ings, exponent))*(1.0/a); //2nd approximation
        return ings;
    }

    /**
     * Calculates the kernel halfwidth m that comes closest to the desired
     * band width, i.e., the frequency where the response decreases to
     * -3 dB, i.e., 1/sqrt(2).
     * @param bandwidth The desired band width, with respect to the sampling frequency.
     *         The value of <code>bandwidth</code> must be less than 0.5
     *         (the Nyquist frequency).
     * @return The kernel halfwidth m.
     */
    public static int bandwidthToM(double bandwidth) {
        final int degree = 4;
        if (bandwidth <= 0 || bandwidth >= 0.5)
            throw new IllegalArgumentException("Invalid bandwidth value: "+bandwidth);
        double radius = (0.74548 + 0.24943*degree)/bandwidth - 1.0;
        return (int)Math.round(radius);
    }

    /** Returns the smoothed data in a new array.
     *  If the data contain NaN values, these are boundaries of a range for smoothing.
     */
    public double[] smooth(double[] data) {
        if (kernel == null)
            return (double[])(data.clone()); // kernelHalfwidth < MIN_HALFWIDTH, don't smooth
        int radius = kernel.length - 1; //kernel halfwidth m; we need m extrapolated points
        int[] rangeLimits = LeedUtils.getRangeLimits(data, 0, -1);
        int extendLeft = radius - rangeLimits[0];
        if (extendLeft < 0) extendLeft = 0;
        int extendRight = radius - data.length + rangeLimits[rangeLimits.length-1];
        if (extendRight < 0) extendRight = 0;
        double[] extendedData = new double[data.length + extendLeft + extendRight];
        System.arraycopy(data, 0, extendedData, extendLeft, data.length);
        Arrays.fill(extendedData, 0, extendLeft, Double.NaN);
        Arrays.fill(extendedData, extendLeft+data.length, extendedData.length, Double.NaN);
        for (int i=0; i<rangeLimits.length/2; i++) {
            int start = i==0 ? 0 : rangeLimits[2*i-1]+extendLeft;   //range to fill with extrapolated data: [start, end[
            int end = rangeLimits[2*i]+extendLeft;
            if (end - start > radius) start = end - radius;
            int inputEnd = rangeLimits[2*i+1]+extendLeft;           //fit input: [end, inputEnd[
            if (inputEnd > end + fitWeights.length)
                inputEnd = end + fitWeights.length;
            extrapolateLeft(extendedData, start, end, inputEnd);
        }
        for (int i=0; i<rangeLimits.length/2; i++) {
            boolean isLast = i == rangeLimits.length/2-1;
            int start = rangeLimits[2*i+1]+extendLeft;              //range to fill with extrapolated data: [start, end[
            int end = isLast ? extendedData.length : rangeLimits[2*(i+1)]+extendLeft;
            int nOverlap = isLast ? 0 : 2*radius - (end-start);     //number of points in overlap with previous extrapolateLeft
            if (nOverlap > end - start) nOverlap = end - start;
            if (end - start > radius) end = start + radius;
            int inputStart = rangeLimits[2*i]+extendLeft;           //fit input: [inputStart, start[
            if (start - inputStart > fitWeights.length)
                inputStart = start - fitWeights.length;
            extrapolateRight(extendedData, start, end, inputStart, nOverlap);
        }
        double[] smoothed = convolve(extendedData);

        double[] out = null;
        if (extendLeft == 0 && extendRight == 0) {
            out = smoothed;
        } else {
            out = new double[data.length];
            System.arraycopy(smoothed, extendLeft, out, 0, out.length);
        }
        int p = 0;
        int nBridgeGap = fitWeights.length/2;   //for short gaps, close the gap by smoothing
        for (int i=0; i<rangeLimits.length/2; i++) {
            if (i > 0 && rangeLimits[2*i] - p < nBridgeGap)
                p = rangeLimits[2*i];       //short gap; skip setting to NaN
            for (; p<rangeLimits[2*i]; p++) //set NaN till start of valid range
                out[p] = Double.NaN;
            p = rangeLimits[2*i+1];         //continue after end of this range
        }
        for (; p<out.length; p++)
            out[p] = Double.NaN;
        return out;
    }

    /** Extrapolates to the left:
     *  Fits a linear function to 'data' in the range [end, inputEnd[ and applies the fit
     *  to extendedData in the range [start, end[. */
    private void extrapolateLeft(double[] extendedData, int start, int end, int inputEnd) {
        LeedLinearRegression linreg = new LeedLinearRegression();
        int fitLength = inputEnd - end;
        if (fitLength < 1) throw new RuntimeException("ExtrapolateLeft: Cannot fit over "+fitLength+" points");
        for (int i=0; i<fitLength; i++)
            linreg.addPointW(i, extendedData[end+i], fitWeights[i]);
        double offset = linreg.getOffset();
        double slope  = linreg.getSlope();
        for (int i=-1, p=end-1; p>=start; i--,p--)
            extendedData[p] = offset + slope*i;
    }

    /** Extrapolates data to the right:
     *  Fits a linear function to 'data' in the range [inputStart, start[ and applies the fit
     *  to the range [start, end[ . */
    private void extrapolateRight(double[] extendedData, int start, int end, int inputStart, int nOverlap) {
        LeedLinearRegression linreg = new LeedLinearRegression();
        int fitLength = start - inputStart;
        if (fitLength < 1) throw new RuntimeException("ExtrapolateRight: Cannot fit over "+fitLength+" points");
        for (int i=0, p=start-1; i<fitLength; i++,p--)
            linreg.addPointW(i, extendedData[p], fitWeights[i]);
        double offset = linreg.getOffset();
        double slope  = linreg.getSlope();
        for (int i=1, p=start; p<end; i++, p++) {
            double extrapolated = offset - slope*i;
            if (p < end-nOverlap) {
                extendedData[p] = extrapolated;
            } else {
                double w = nOverlap == 1 ? 0.5 : (end-start-i)/(double)(nOverlap-1);  //weight of extrapolated vs. previous extrapolateLeft
                extendedData[p] = w*extrapolated + (1.0 - w)*extendedData[start-1+i];
            }
        }
    }

    /**
     * Smooths the data with the kernel created in the constructor,
     * except for the near-boundary points.
     * @param data The input data.
     * @return  The smoothed data. If the kernel halfwidth is m,
     *          the first and last m points remain zero.
     */
    double[] convolve(double[] data) {
        double[] out = new double[data.length];
        int radius = kernel.length - 1; //kernel halfwidth m
        for (int i=radius; i<data.length-radius; i++) {
            double sum = kernel[0]*data[i];
            for (int j=1; j<kernel.length; j++) {
                sum += kernel[j]*(data[i-j]+data[i+j]);
            }
            out[i] = sum;
        }
        return out;
    }

    /**
     * Creates a kernel and returns it.
     * @param isMS1  if true, calculates the kernel for the MS1 variant.
     *          Otherwise, standard MS kernels are used.
     * @param degree The degree n of the kernel, i.e. the polynomial degree of a Savitzky-Golay filter
     *          with similar passband, must be 2, 4, ... MAX_DEGREE.
     * @param m The half-width of the kernel (the resulting kernel
     *          has <code>2*m+1</code> elements).
     * @param coeffs Correction parameters for a flatter passband, or
     *          null for no correction (used for degree 2).
     * @return  One side of the kernel, starting with the element at the
     *          center. Since the kernel is symmetric, only one side with
     *          <code>m+1</code> elements is needed.
     */
    static double[] makeKernel(int m) {
        if (m < 4) throw new IllegalArgumentException("Internal error: Kernel halfwidth m="+m+" too low");
        double[] kernel = new double[m+1];
        double sum = 0;
        for (int i=0; i<=m; i++) {
            double x = i*(1./(m+1)); //0 at center, 1 at zero
            double sincArg = Math.PI*0.5*(DEGREE+4)*x;
            double k = i==0 ? 1 : Math.sin(sincArg)/sincArg;
            final double decay = 4;  //decay alpha
            k *= Math.exp(-x*x*decay) + Math.exp(-(x-2)*(x-2)*decay) +
                    Math.exp(-(x+2)*(x+2)*decay) - 2*Math.exp(-decay) - Math.exp(-9*decay);
            kernel[i] = k;
            sum += k;
            if (i > 0) sum += k;    //off-center kernel elements appear twice
        }
        for (int i=0; i<=m; i++)
            kernel[i] *= 1./sum;    //normalize the kernel to sum = 1
        return kernel;
    }

    /**
     * Returns the weights for the linear fit used for linear extrapolation
     * at the end. The weight function is a Hann (cos^2) function. For beta=1
     * it would decay to zero at the position of the
     * first zero of the sinc function in the kernel. Larger beta values lead
     * to stronger noise suppression near the edges, but the smoothed curve
     * does not follow the input as well as for lower beta (for high degrees,
     * also leading to more ringing near the boundaries).
     * @param m The half-width of the kernel (the full kernel
     *          has 2*m+1 elements).
     * @return  The fit weights, with array element [0] corresponding
     *          to the data value at the very end.
     */
    static double[] makeFitWeights(int m) {
        double firstZero = (m+1)/(1.5+0.5*DEGREE);
        double beta = 0.70 + 0.14*Math.exp(-0.60*(DEGREE-4));
        int fitLength = (int)Math.ceil(firstZero*beta);
        double[] weights = new double[fitLength];
        for (int p=0; p<fitLength; p++)
            weights[p] = sqr(Math.cos(0.5*Math.PI/(firstZero*beta)*p));
        return weights;
    }

    /** Returns the square of a number */
    static double sqr(double x) {return x*x;}

}
