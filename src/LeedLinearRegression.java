import ij.*;
import ij.IJ;
import ij.gui.*;
import ij.process.*;
import ij.gui.Plot;
/**
 *  This class implements linear regression for a single function of one variable
 *  with (optional) weight factors.
 *  Note that the equations related to the errors (e.g. getFitValueWithMinSlope)
 *  are valid only in the limit of a large number of points, as we do not correct
 *  for the number of degrees of freedom (which would make sense with counters
 *  only, not with weights).
 *
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
 
public class LeedLinearRegression {
    /** number of data points, or sum of weights*/
    protected double counter = 0;
    /** sum of all x values*/
    protected double sumX = 0;
    /** sum of all y values*/
    protected double sumY = 0;
    /** sum of all x*y products*/
    protected double sumXY = 0;
    /** sum of all squares of x*/
    protected double sumX2 = 0;
    /** sum of all squares of y*/
    protected double sumY2 = 0;
    /** result of the regression: offset */
    protected double offset = Double.NaN;
    /** result of the regression: residuals */
    protected double rmsResiduals = Double.NaN;
    /** result of the regression: slope */
    protected double slope = Double.NaN;
    /** result: reciprocal of (n^2 * squared x_stddev) */
    protected double oneOverN2xRmsSq;
    /** whether the results (offset, slope) have been calculated already */
    protected boolean calculated = false;

    /** Clear the linefit */
    public void clear() {
        counter = 0;
        sumX = 0;
        sumY = 0;
        sumXY = 0;
        sumX2 = 0;
        sumY2 = 0;
        offset = Double.NaN;
        slope = Double.NaN;
        calculated = false;
    }

    /** Add a point x,y with weight unless any value is NaN*/
    public void addPointW(double x, double y, double weight) {
        if (Double.isNaN(x) || Double.isNaN(y) || Double.isNaN(weight)) return;
        counter += weight;
        sumX += weight*x;
        sumY += weight*y;
        sumXY += weight*x*y;
        sumX2 += weight*x*x;
        sumY2 += weight*y*y;
        calculated = false;
    }

    /** Add a point x,y with weight 1; neither x nor y must be NaN for having a NaN result*/
    public void addPoint(double x, double y) {
        counter++;
        sumX += x;
        sumY += y;
        sumXY += x*y;
        sumX2 += x*x;
        sumY2 += y*y;
        calculated = false;
    }

    /** Adds the x, y points from the arrays */
    public void addPoints(double[] x, double[] y) {
        for (int i=0; i<x.length; i++)
            addPoint(x[i], y[i]);
    }

    /** Do the actual regression calculation - only necessary if data are added after getting fit results or to
     * undo a previous setSlope command */
    public void reCalculate() {
        double stdX2TimesN = sumX2-sumX*sumX*(1/counter);
        double stdY2TimesN = sumY2-sumY*sumY*(1/counter);
        oneOverN2xRmsSq = 1./(stdX2TimesN*counter);    //needed for error estimate
        if (counter>0) {
            slope=(sumXY-sumX*sumY*(1/counter))/stdX2TimesN;
            if (Double.isNaN(slope)) slope=0;           //slope 0 if fit was unsuccessful
        } else
            slope = Double.NaN;
        offset = (sumY-slope*sumX)/counter;
        double r = (sumXY - sumX*sumY*(1/counter)) / Math.sqrt(stdX2TimesN*stdY2TimesN);
        if (counter > 1e-100 && Double.isNaN(r)) r = 1; //if stddev of y = 0;
        rmsResiduals = stdY2TimesN*(1/counter)*(1 - r*r);
        if (rmsResiduals < 0) rmsResiduals = 0;         //numerical errors could cause rmsResiduals<0
        // alternative equation: (sumY2 + slope*slope*sumX2 + offset*offset*counter - 2*slope*sumXY - 2*offset*sumY + 2*offset*slope*sumX) / counter;
        calculated = true;
		//DEBUG IJ.log("off="+(float)offset+" sl="+(float)slope+" r2="+(float)(r*r)+" rmsR="+(float)rmsResiduals);

    }

    /** Do the regression calculation with a forced slope. Returns the offset. */
    public double setSlope(double slope) {
        if (counter>0) {
            this.slope=slope;
            offset=(sumY-slope*sumX)/counter;
        }
        calculated = true;
        return offset;
    }

    /** Returns the offset (intersection on y axis) of the fit */
    public double getOffset() {
        if (!calculated) reCalculate();
        return offset;
    }

    /** Returns the slope of the line of the fit */
    public double getSlope() {
        if (!calculated) reCalculate();
        return slope;
    }

    /** Returns the fit function value at a given x */
    public double getFitValue(double x) {
        if (!calculated) reCalculate();
        return offset + x*slope;
    }

    /** Returns the values of the fit function value at the x avlues of the array */
    public double[] getFitValues(double[] x) {
        double[] out = new double[x.length];
        for (int i=0; i<x.length; i++)
            out[i] = getFitValue(x[i]);
        return out;
    }

    /** Returns a rough estimate of the fit function value at a given x
     *  assuming the smallest slope compatible with the data, if the
     *  y values have an uncertainty of sqrt(yErrSqr).
     *  Also works with weights if there are n values with equal weight */
    public double getFitValueWithMinSlope(double x, double yErrSqr, int n) {
        if (!calculated) reCalculate();
        if (n < 2) return offset;
		double stdX2 = (sumX2-sumX*sumX*(1/counter))*(1/counter);
		yErrSqr = Math.max(rmsResiduals, yErrSqr);
		double varSlope = Math.sqrt(yErrSqr/(stdX2*(n-1)));
        double minSlope = Math.abs(slope) > varSlope ? slope - Math.signum(slope) * varSlope : 0;
        double dx = x - sumX*(1/counter);
        return offset + sumX*(1/counter)*slope + dx*minSlope;
    }


    /** Returns a measure of the fit weight at a given point (proportional to 1/error^2)
     *  Note that the error would be very roughly 1/sqrt(weight * n),
     *  where n are the number of points minus 2 */
    public double getFitWeight(double x) {
        if (!calculated) reCalculate();
        return 1./(rmsResiduals*(1 + (x*counter - sumX)*(x*counter - sumX)*oneOverN2xRmsSq)+1e-50);
    }

    /** Returns the mean X value */
    public double getMeanX() {return sumX/counter;}

    /** Returns the mean Y value */
    public double getMeanY() {return sumY/counter;}

    /** Returns the mean square deviation from average x value */
    public double getMeanDX2() {return sumX2/counter-sumX*sumX/counter/counter;}

    /** Returns the mean square deviation from average y value */
    public double getMeanDY2() {return sumY2/counter-sumY*sumY/counter/counter;}

    /** Returns the mean square deviation from regression y value */
    public double getRmsResiduals() {
        if (!calculated) reCalculate();
        return rmsResiduals;
    }

    /** Get counter (or sum of weights) value */
    public double getCounter() {return counter;}

    /** Whether we have any data at all */
    public boolean dataPresent() {return counter>0;}

    /** String representation */
    public String toString() {
        String out = "LineFit of "+counter+" points";
        if (calculated)
            out += "; offset="+(float)offset+" slope="+(float)slope;
        else
            out += ", not evaluated yet";
        return out;
    }

    /** Constructor for testing only, plots fit & errors from 1./getFitWeight(x) */ /*
    public LeedLinearRegression() {
        double[] xdata = new double[]{3,4,5,6,7,11};
        double[] ydata = new double[]{4,4,5,6,5,10};
        //double[] ydata = new double[]{4,4,4,4,4,4};
        addPoints(xdata,ydata);
        Plot plot = new Plot("regression test","x","y");
        plot.addPoints(xdata,ydata,Plot.CIRCLE);
        double offs=getOffset();
        double slop=getSlope();
        plot.drawLine(0, offs, 15, offs+15*slop);
        double[] xe = new double[30], yePlus = new double[30], yeMinus = new double[30];
        for (int i=0; i<xe.length; i++) {
            double x = i-5; xe[i] = x;
            double err = Math.sqrt(1./getFitWeight(x));
            yePlus[i]  = offs+slop*x+err;
            yeMinus[i] = offs+slop*x-err;
        }
        plot.addPoints(xe, yePlus, Plot.LINE);
        plot.addPoints(xe, yeMinus,Plot.LINE);
        plot.show();
    }/**/
}
