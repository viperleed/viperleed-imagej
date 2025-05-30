import ij.IJ;
import ij.util.Tools;
import ij.gui.Plot;

/** This class calculates of Pendry's R factor.
 *  It contains provisions to extend it to calculate a modified R factor,
 *  which must be also based on a Y function.
 *  Code for calculating the energy-dependent Y function is also provided.
 *  The derivatives are calculated in a very simple way, from the difference of
 *  the two values above and below the current point.
 *  Note that the viperleed.calc/TensErLEED R factor code calculates the derivatives
 *  with finite difference coefficients over 7 points (for interior points;
 *  fewer near the boundary). In case of large energy steps, this is more
 *  accurate than the simple difference done here.
 *  https://en.wikipedia.org/wiki/Finite_difference_coefficient
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


public class LeedRFactor {
    public static final int R_PENDRY = 0, R_PMOD = 1;   //R factor types; what follows are indices in output array of other outputs
    public static final int MAX_1 = 2, MAX_2 = 3;       //highest intensity in overlap region of the two curves (for normalization when plotting)
    public static final int AVG_INTENSITY = 4;          //average intensity of the curves (actually sqrt(average intensity 1 * average intensity 2)
    public static final int N_OVERLAP = 5;              //number of points in common region of the curves
    public static final int RATIO = 6;                  //ratio of data2/data1 average
    public static final int OUTPUT_SIZE = 7;            //the length of the output array 0...6
    public static final String[] R_FACTOR_NAMES = new String[] {"R_Pendry", "R_PeMod"};


    /** Returns the R factors in the range where the two curves overlap, the average intensity
     *  of the curves and maximum intensity of each curve in the overlap region, as well as
     *  the number of points of the overlap.
     *  The input arrays 'data1', 'data2' should be the intensities over energy (both with the
     *  same energy steps; the 'shift' parameter should be used if the start energies are different).
     *  When 'shift' is nonzero, the 'data2' array is shifted by 'shift' points to the left.
     *  'v0iOverStep' is the imaginary part of the inner potential divided by the energy step.
     *  If the data contain negative points, an offset is added to all points
     *  of this data array to make the data non-negative.
     *  For simplicity, differentiation is done by taking the differenc of the data values at i+1 and i-1
     *  The output array contains the R factors R_Pendry and R_PMod, the max intensities of data1 and data2
     *  as well as the average intensity in the overlap range, and the overlap range in points.
     *  An array for the output may be provided or null; the output is always in the array returned.
     *  The output array contains [0] the R factor and [1] the number of points used for the comparison
     *  (note that the two end points of each range where both data arrays are non-NaN are discarded) */
    public static double[] getRFactor(double[] data1, double[] data2, int shift, double v0iOverStep, double[] output) {
        return getRFactor(data1, data2, 0, data1.length, shift, v0iOverStep, output);
    }

    /** Returns the R factors in the range where the two curves overlap, the average intensity
     *  of the curves and maximum intensity of each curve in the overlap region, as well as
     *  the number of points of the overlap.
     *  The input arrays 'data1', 'data2' should be the intensities over energy (both with the
     *  same energy steps; the 'shift' parameter should be used if the start energies are different).
     *  The calculation is restricted to an range between indices iStart and iEnd  of data1
     *  (including iStart, excluding iEnd).
     *  When 'shift' is nonzero, the 'data2' is shifted by 'shift' points to the left.
     *  'v0iOverStep' is the imaginary part of the inner potential divided by the energy step.
     *   If the data contain negative points, an offset is added to all points
     *  of this data array to make the data non-negative.
     *  For simplicity, differentiation is done by taking the differenc of the data values at i+1 and i-1
     *  The output array contains the R factors R_Pendry and R_PMod, the max intensities of data1 and data2
     *  as well as the average intensity in the overlap range, and the overlap range in points.
     *  An array for the output may be provided or null; the output is always in the array returned.
     *  The output array contains [0] the R factor and [1] the number of points used for the comparison
     *  (note that the two end points of each range where both data arrays are non-NaN are discarded) */
    public static double[] getRFactor(double[] data1, double[] data2, int iStart, int iEnd, int shift, double v0iOverStep, double[] output) {
        if (output == null || output.length < 2)
            output = new double[OUTPUT_SIZE];

        double offset1 = 0, offset2 = 0;        //in case of negative data, find the most negative values in overlap range
        int maxEnd = Math.min(data1.length, data2.length-shift);
        for (int i=Math.max(iStart, -shift); i<Math.min(iEnd, maxEnd); i++) {
            double v1 = data1[i], v2 = data2[i+shift];
            if (!Double.isNaN(v1) && !Double.isNaN(v2)) {
                if (v1 < offset1) offset1 = v1;
                if (v2 < offset2) offset2 = v2;
            }
        }
        offset1 = -offset1; offset2 = -offset2; //these offsets will be added

        double nPoints = 0;
        double sumNominP = 0, sumDenomP = 0;    //nominator and denominator in Pendry R factor
        double sumNominM = 0, sumDenomM = 0;    //nominator and denominator in Modified Pendry R factor
        double maxIntensity1 = 0, maxIntensity2 = 0;  //maximum intensity in overlap range
        double sumIntensity1 = 0, sumIntensity2 = 0;  //sum over intensity in overlap range
        int nPreviousOk = 0;
        double mid1 = 0, mid2 = 0;              //mid data values for the two curves
        double left1 = 0, left2 = 0;            //data values before mid
        iStart--;                               //when not starting at 0, need to read the previous point for differentiating
        if (iStart < 0) iStart = 0;
        if (iStart < -shift) iStart = -shift;   //for data2, we access index i+shift
        iEnd++;                                 //when the end is not array end, we need to read the next point as well
        if (iEnd > data1.length) iEnd = data1.length;
        if (iEnd > data2.length-shift) iEnd = data2.length-shift;
        for (int i=iStart; i<iEnd; i++) {
            double right1 = data1[i] + offset1, right2 = data2[i+shift] + offset2;
            if (Double.isNaN(right1) || Double.isNaN(right2)) {
                nPreviousOk = 0;                //no valid data, so we can't do anything
            } else {                            //valid data at least at this point
                if (nPreviousOk >= 2) {         //enough data to differentiate at the mid point
                    double y1 = calculateY(left1, mid1, right1, v0iOverStep, R_PENDRY);   //factor 0.5 because we have twice the derivative
                    double y2 = calculateY(left2, mid2, right2, v0iOverStep, R_PENDRY);   //  (i.e., difference between points two steps apart)
                    sumNominP += sqr(y2 - y1);
                    sumDenomP += sqr(y2) + sqr(y1);
                    /*
                    y1 = calculateY(left1, mid1, right1, v0iOverStep, R_PMOD);
                    y2 = calculateY(left2, mid2, right2, v0iOverStep, R_PMOD);
                    sumNominM += sqr(y2 - y1);
                    sumDenomM += sqr(y2) + sqr(y1);
                    */
                    if (mid1 > maxIntensity1) maxIntensity1 = mid1;
                    if (mid2 > maxIntensity2) maxIntensity2 = mid2;
                    sumIntensity1 += mid1;
                    sumIntensity2 += mid2;
                    nPoints++;
                }
                left1 = mid1; left2 = mid2;
                mid1 = right1; mid2 = right2;
                nPreviousOk++;
            }
        }
        double rFactorP = sumNominP/sumDenomP;
        //double rFactorM = sumNominM/sumDenomM;
        output[R_PENDRY] = rFactorP;
        //output[R_PMOD] = rFactorM;
        output[MAX_1] = maxIntensity1;
        output[MAX_2] = maxIntensity2;
        output[AVG_INTENSITY] = Math.sqrt(sumIntensity1*sumIntensity2)/nPoints;
        output[N_OVERLAP] = nPoints;
        output[RATIO] = sumIntensity2/sumIntensity1;
        return output;
    }

    /** Returns the range(s) where two curves overlap (i.e., both are not NaN
     *  for at least 3 points) in the form:
     *  start1, end1, start2, end2, ... nOverlap,
     *  where the ranges include the 'start', but not the 'end' points,
     *  and the last element 'nOverlap' is the total number of points in the overlap regions.
     *  Returns null if there is no overlap. */
    public static int[] getOverlap(double[] data1, double[] data2) {
        LeedIntegerArray ranges = new LeedIntegerArray();
        int nOverlap = 0;
        boolean overlap = false;  //whether we are currently in an overlap region
        boolean foundNonZero1 = false, foundNonZero2 = false; //TensErLEED has instead of NaN
        int iEnd = Math.min(data1.length, data2.length);
        for (int i=0; i<iEnd; i++) {
            if (!foundNonZero1 && data1[i] != 0) foundNonZero1 = true;
            if (!foundNonZero2 && data2[i] != 0) foundNonZero2 = true;
            if (!foundNonZero1 || !foundNonZero2) continue;
            if (!overlap && !Double.isNaN(data1[i]) && !Double.isNaN(data2[i])) {
                ranges.add(i);
                overlap = true;
            }
            if (overlap && (Double.isNaN(data1[i]) || Double.isNaN(data2[i]))) {
                if (i - ranges.getLast() >= 3) {
                    nOverlap += i - ranges.getLast();
                    ranges.add(i);                  //range ends here
                } else
                    ranges.removeLast(1);           //range too short, delete previous start
                overlap = false;
            }
        }
        if (overlap) {
            if (iEnd - ranges.getLast() >= 3) {
                nOverlap += iEnd - ranges.getLast();
                ranges.add(iEnd);                   //range ends here
            } else
                ranges.removeLast(1);               //range too short, delete previous start
        }
        if (nOverlap == 0) return null;
        ranges.add(nOverlap);
        return ranges.toArray();
    }

    /** Returns the "Y" normalized to be within -0.99...+0.99 for an I(V) curve */
    public static double[] getYcurve(double[] data, int type, double v0iOverStep, double[] output) {
        if (output == null || output.length < data.length)
            output = new double[data.length];
        double offset = getOffset(data, 0, data.length);        //offset to add in case of negative data

        double nPoints = 0;
        double sumNomin = 0, sumDenom = 0;      //nominator and denominator in R factor
        int nPreviousOk = 0;
        double mid = 0;                        //mid data value
        double left = 0;                       //data value before mid
        for (int i=0; i<data.length; i++) {
            output[i] = Double.NaN;
            double right = data[i] + offset;
            if (Double.isNaN(right)) {
                nPreviousOk = 0;                //no valid data, so we can't do anything
            } else {                            //valid data at least at this point
                if (nPreviousOk >= 2)           //enough data to differentiate at the mid point
                    output[i-1] = calculateY(left, mid, right, v0iOverStep, type);
                left = mid;
                mid = right;
                nPreviousOk++;
            }
        }
        double[] minMax = Tools.getMinMax(output);
        double normalizeFactor = 0.99 / Math.max(minMax[1], -minMax[0]);
        for (int i=0; i< output.length; i++)
            output[i] *= normalizeFactor;
        return output;
    }

    /** Returns the R factor between two Y curves */
    public static double getRFactor(double[] y1, double[] y2) {
        double sumNomin=0, sumDenom=0;
        for (int i=0; i<Math.min(y1.length, y2.length); i++) {
            if (Double.isNaN(y1[i]) || Double.isNaN(y2[i])) continue;
            sumNomin += sqr(y2[i] - y1[i]);
            sumDenom += sqr(y2[i]) + sqr(y1[i]);
        }
        return sumNomin/sumDenom;
    }

    /** Returns the 'Y' value from three successive points */
    static double calculateY(double left, double mid, double right, double v0iOverStep, int type) {
        if (type == R_PENDRY) {
            double twoL = (right - left)/(mid + 1e-100);       //2*logarithmic derivative
            double y = twoL/(1 + sqr(0.5*v0iOverStep * twoL)); //factor 2 does not hurt, will be normalized
            return y;
        } else { // R_PMOD, not impelemented yet
            return 0;
        }
    }

    /** Returns the offset that has to be added to make the array non-negative. */
    static double getOffset(double[] data, int iStart, int iEnd) {
        double offset = 0;
        for (int i=iStart; i<iEnd; i++)
            if (data[i] < offset)
                offset = data[i];
        return -offset;
    }

    static double sqr(double x) {
        return x*x;
    }
}
