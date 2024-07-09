import ij.*;
import ij.IJ;
import ij.gui.*;
import ij.process.*;
/**
 *  Fits an affine transformation
 *    y1 = a1 + b1*x1 + c1*x2
 *    y2 = a2 + b2*x1 + c2*x2
 *  for the points with measured coordinates (y1, y2)
 *  and calculated coordinates (x1, x2), and given weights for each point.
 *  (Note: Here, naming of x, y is not for the coordinates but for the variables
 *   of the fitting problem in each direction)
 *  The fit output contains:
 *  - The offsets a1, a2 (array indices X_OFFSET, Y_OFFSET)
 *  - The determinant of the [b1,c1; b2,c2] matrix, the rotation angle.
 *  Assuming a conformal transtormation, the determinant is the the square of the
 *  scale factor.
 *  - The rotation angle in degrees of y1 vx x1, assuming a conformal transformation.
 *  - The weight, which is an indication of the significance (non-degeneracy)
 *  in arbitrary units.
 *  - The sum of the statistical weights of the input data.
 *
 *  When applied to measured vs calculated spot positions in LEED I(V),
 *  the energy-dependent offsets a1, a2 are an indication of residual fields
 *  or poor alignment of the electron source; the energy-dependent rotation
 *  angle of the LEED pattern indicates longitudinal magnetic fields [these
 *  are irrelevant for LEED I(V)], and the energy dependence of the determinant
 *  (scale factor) indicates an energy offset, e.g. due to the difference
 *  of the vacuum levels of the filament and the region between sample and
 *  LEED screen, or charging of the sample, or field distortions due to the
 *  fringe field corrector electrodes of MCP-LEED optics.
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

public class Leed2DRegression {
    /** length of output array for the coefficients */
    public static final int OUTPUT_SIZE=6;
    /** indices in output array */
    public static final int X_OFFSET=0, Y_OFFSET=1, DETERMINANT=2, ANGLE=3, WEIGHT=4, SUM_DATA_WEIGHTS=5;
    double sumW, sumx1, sumx2, sumx1Sq, sumx2Sq, sumx1x2,
            sumy1, sumy2, sumx1y1, sumx1y2, sumx2y1, sumx2y2;

    /** Adds a point with measured screen coordinates (y1Exp, y2Exp) and
     *  calculated screen coordinates (x1Calc, x2Calc). */
    public void addPoint(double y1Exp, double y2Exp, double x1Calc, double x2Calc, double weight) {
        sumW    += weight;
        sumx1   += weight*x1Calc;
        sumx2   += weight*x2Calc;
        sumx1x2 += weight*x1Calc*x2Calc;
        sumx1Sq += weight*x1Calc*x1Calc;
        sumx2Sq += weight*x2Calc*x2Calc;
        sumy1   += weight*y1Exp;
        sumy2   += weight*y2Exp;
        sumx1y1 += weight*x1Calc*y1Exp;
        sumx1y2 += weight*x1Calc*y2Exp;
        sumx2y1 += weight*x2Calc*y1Exp;
        sumx2y2 += weight*x2Calc*y2Exp;
    }

    /** Clears all data */
    public void clear() {
        sumW   = 0;
        sumx1  = 0;
        sumx2  = 0;
        sumx1x2 = 0;
        sumx1Sq = 0;
        sumx2Sq = 0;
        sumy1  = 0;
        sumy2  = 0;
        sumx1y1 = 0;
        sumx1y2 = 0;
        sumx2y1 = 0;
        sumx2y2 = 0;
    }

    /** If x1 is not NaN, applies the fit result to x1, x2 and
     *  returns the values as a 2-element array.
     *  If x1 is NaN, calculates the fit coefficients and returns an array of
     *  {a1, a2, determinant, weight, sumDataWeights} where 'determinant'
     *  is the square of the scale factor and 'weight' an indication
     *  of the certainty of the result. With sum of all data weights,
     *  one can check for a non-degenerate 2D range of the 'calc' data:
     *  weight/sumDataWeights should be comparable order of magnitude
     *  as the 'calc' data.
     *  The output array can be supplied or be null. */
    public double[] getFit(double x1, double x2, double[] output) {
        double sumx1x2r = sumx1x2 - (1./sumW)*sumx1*sumx2;
        double sumx1Sqr = sumx1Sq - (1./sumW)*sumx1*sumx1;
        double sumx2Sqr = sumx2Sq - (1./sumW)*sumx2*sumx2;
        double sumx1y1r = sumx1y1 - (1./sumW)*sumx1*sumy1;
        double sumx1y2r = sumx1y2 - (1./sumW)*sumx1*sumy2;
        double sumx2y1r = sumx2y1 - (1./sumW)*sumx2*sumy1;
        double sumx2y2r = sumx2y2 - (1./sumW)*sumx2*sumy2;
        double denom = sumx1Sqr*sumx2Sqr - sumx1x2r*sumx1x2r;
        double b1 = (sumx2Sqr*sumx1y1r - sumx1x2r*sumx2y1r)*(1./denom);
        double c1 = (sumx1Sqr*sumx2y1r - sumx1x2r*sumx1y1r)*(1./denom);
        double b2 = (sumx2Sqr*sumx1y2r - sumx1x2r*sumx2y2r)*(1./denom);
        double c2 = (sumx1Sqr*sumx2y2r - sumx1x2r*sumx1y2r)*(1./denom);
        double a1 = (sumy1 - b1*sumx1 - c1*sumx2)*(1./sumW);
        double a2 = (sumy2 - b2*sumx1 - c2*sumx2)*(1./sumW);
        if (Double.isNaN(x1)) { //we have to return the offsets, determinant etc
            double determinant = b1*c2 - b2*c1;
            double weight = Math.sqrt(Math.abs(denom));
            if (output == null || output.length < OUTPUT_SIZE)
                output = new double[OUTPUT_SIZE];
            output[X_OFFSET] = a1;
            output[Y_OFFSET] = a2;
            output[DETERMINANT] = b1*c2 - b2*c1;
            output[ANGLE] = Math.toDegrees(output[DETERMINANT] > 0 ?
                    Math.atan2(c1-b2, b1+c2) : Math.atan2(c1+b2, b1-c2));
            output[WEIGHT] = Math.sqrt(Math.abs(denom));
            output[SUM_DATA_WEIGHTS] = sumW;
        } else {                //we have to calculate the fit result for given x1, x2
            if (output == null || output.length < 2)
                output = new double[2];
            output[0] = a1 + b1*x1 + c1*x2;
            output[1] = a2 + b2*x1 + c2*x2;
        }
        return output;
    }

    /** Constructor for testing only; comment out for actual use!
     *  Results with the first data set should be offset&slopes
     *     -4.103581, 0.086409, 0.087602 for expDa1,
     *     -1.595220, 0.051618, 0.119580 for expDa2,
     * determinant = 0.0058109 *//**
    public Leed2DRegression() {
        //double[] expDa1 = new double[] {1,2,1,3,2,3,3,4,4,3,5,3,3,2,4,5,5,5,4,3};
        //double[] expDa2 = new double[] {2,3,2,4,3,4,4,5,5,4,9,8,8,7,5,6,6,6,5,4};
        //double[] x1Calc = new double[] {40,45,38,50,48,55,53,55,58,40,55,48,45,55,60,60,60,65,50,58};
        //double[] x2Calc = new double[] {25,20,30,30,28,30,34,36,32,34,38,28,30,36,34,38,42,38,34,38};
        for (int i=0; i<expDa1.length; i++)
            addPoint(expDa1[i], expDa2[i], x1Calc[i], x2Calc[i], 0.12345);
        double[] out = getFit(Double.NaN, Double.NaN, null);
        IJ.log("offs a1="+IJ.d2s(out[X_OFFSET]) + " sumy1="+sumy1 + " sumx1="+sumx1 + " sumx2="+sumx2);
        IJ.log("offs a2="+IJ.d2s(out[Y_OFFSET]) + " sumx1y1="+sumx1y1 + " sumx2y1="+sumx2y1 + " sumx1x2="+sumx1x2);
        IJ.log("det="+IJ.d2s(out[DETERMINANT],8) + " sumx1Sq="+sumx1Sq + " sumx2Sq="+sumx2Sq);
    } /**/
}
