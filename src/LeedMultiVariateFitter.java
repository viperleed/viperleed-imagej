import ij.*;
import ij.IJ;
import ij.gui.*;
import ij.process.*;
import java.util.Arrays;
/**
 *  Fits a linear multi-parameter function to data
 *    y = b0*x0 + b1*x1 + ...
 *  where x0, x1, ... are n independent variables
 *  and b0, b1, ... are the parameters; these should be determined.
 * 
 *  Note that numeric stability is limited; especially for large numbers
 *  of parameters the x values should be scaled such that they have
 *  root-mean-square values around 1.
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

public class LeedMultiVariateFitter {
    private int n;                      //number of parameters
    private int count;                  //number of data points
    private double[][] sumXiXj;
    private double[] sumXiY;
    private double[] parameters = null; //the final parameters

    /** Constructor, sets the number of parameters */
    public LeedMultiVariateFitter(int nParams) {
        this.n = nParams;
        sumXiXj = new double[n][n];
        sumXiY  = new double[n];
    }

    /** Clears everything, for a new fit */
    public void clear() {
        count = 0;
        Arrays.fill(sumXiY, 0);
        for (int i=0; i<n; i++)
            Arrays.fill(sumXiXj[i], 0);
        parameters = null;
    }

    /** Adds a data point with a given experimental value and the x values.
     *  The array size of x must be equal to the number of parameters */
    public void addPoint(double expY, double[] x) {
        for (int i=0; i<n; i++) {
            sumXiY[i] += expY*x[i];
            sumXiXj[i][i] += x[i]*x[i];
            for (int j=0; j<i; j++)
                sumXiXj[i][j] += x[i]*x[j];
        }
        count++;
        parameters = null;
    }

    /** Returns the number of data points */
    public int getNPoints() {
        return count;
    }

    /** Returns the fit parameters. */
    public double[] getFitParameters() {
        if (parameters == null) {
            parameters = new double[n];
            for (int i=0; i<n; i++)     //first complete the symmetric matrix
                for (int j=0; j<i; j++)
                    sumXiXj[j][i] = sumXiXj[i][j];

            boolean success = LeedUtils.invertSymmetricMatrix(sumXiXj);
            if (success) {
                for (int i=0; i<n; i++)
                    for (int j=0; j<n; j++)
                        parameters[i] += sumXiXj[i][j]*sumXiY[j];
            } else
                Arrays.fill(parameters, Double.NaN);
            return parameters;
        }
        return parameters;
    }

	/** Returns the fit value for given x values.
     *  The array size of x must be equal to the number of parameters */
	public double getFitValue(double[] x) {
        if (parameters == null)
            getFitParameters();
        double result = 0;
        for (int i=0; i<n; i++)
            result += parameters[i]*x[i];
        return result;
    }

    /** This constructor is for testing only.
     *  The result should be -4.103581, 0.086409, 0.087602 *//**
    public LeedMultiVariateFitter() {
        LeedMultiVariateFitter fit = new LeedMultiVariateFitter(3);
        double[] expY = new double[] {1,2,1,3,2,3,3,4,4,3,5,3,3,2,4,5,5,5,4,3};
        double[] x1 = new double[] {40,45,38,50,48,55,53,55,58,40,55,48,45,55,60,60,60,65,50,58};
        double[] x2 = new double[] {25,20,30,30,28,30,34,36,32,34,38,28,30,36,34,38,42,38,34,38};
        for (int i=0; i<expY.length; i++)
            fit.addPoint(expY[i], new double[] {1, x1[i], x2[i]});
        double[] params = fit.getFitParameters();
        for (int i=0; i<3; i++) IJ.log((char)('a'+i)+"="+(float)params[i]);
    }/**/


}
