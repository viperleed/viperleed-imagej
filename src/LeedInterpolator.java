import ij.*;
import ij.gui.*;
import ij.io.*;
import ij.measure.*;
import ij.process.*;
import ij.text.*;
import ij.plugin.*;
import ij.util.Tools;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;

/**
 *  This class interpolates LEED I(V) curves to a finer energy grid.
 *  It uses the ImageJ SplineFitter class.
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
 *  @author Michael Schmid, IAP/TU Wien, 2021-2024
 */


public class LeedInterpolator {

    /** Returns a LeedIVData object created by interpolation to a different (finer) energy grid
     *  with spacing 'newStep', in the same energy range as the input.
     *  Data are interpolated using cubic splines.
     *  Note that the spotNames and spotIndices arrays are kept, not duplicated; modification
     *  of these arrays affects the original and the interpolated LeedIVData objects. */
    public static LeedIVData interpolate(LeedIVData ivData, double newStep) {
        double[] minMax = Tools.getMinMax(ivData.energies);
        return interpolate(ivData, minMax[0], minMax[1], newStep);
    }

    /** Returns a LeedIVData object created by interpolation to a different (finer) energy grid.
     *  The output range will be firstEnergy, lastEnergy (inclusive, except for cases where these
     *  energies are not multiples of the fine grid spacing 'newStep')
     *  Data are interpolated using cubic splines.
     *  Note that the spotNames and spotIndices arrays are kept, not duplicated; modification
     *  of these arrays affects the original and the interpolated LeedIVData objects. */
    public static LeedIVData interpolate(LeedIVData ivData, double firstEnergy, double lastEnergy, double newStep) {
        if (!(newStep > 0)) throw new IllegalArgumentException ("Step "+newStep+" not positive");
        double[] newEnergies = getNewE(firstEnergy, lastEnergy, newStep);
        double[][] newData = new double[ivData.data.length][];
        for (int ic=0; ic<ivData.data.length; ic++) {
            if (ivData.data[ic] == null) continue;
            newData[ic] = getInterpolated(ivData.energies, newEnergies, ivData.data[ic]);
        }
        LeedIVData newLeedIVData = ivData.shallowClone();
        newLeedIVData.energies = newEnergies;
        newLeedIVData.data = newData;
        return newLeedIVData;
    }


    /** For interpolation to a new energy step, returns the new independent energy axis */
    public static double[] getNewE(double firstEnergy, double lastEnergy, double newStep) {
        double newFirst = Math.ceil(firstEnergy/newStep - 1e-8)*newStep;
        double newLast = Math.floor(lastEnergy/newStep + 1e-8)*newStep;
        int nNew = (int)Math.round((newLast-newFirst)/newStep+1);
        double[] newE = new double[nNew];

        for (int i=0; i<newE.length; i++)
            newE[i] = newFirst + newStep*i;
        return newE;
    }

    /** Returns data arrays interpolated from the old to the new x axis. */
    public static double[] getInterpolated(double[] oldX, double[] newX, double[] oldData) {
        double newStep = (newX[newX.length-1] - newX[0])/(newX.length-1);
        double[] newData = new double[newX.length];
        Arrays.fill(newData, Double.NaN);
        int[] rangeLimits = LeedUtils.getRangeLimits(oldData, 0, oldData.length);
        if (rangeLimits[0] < 0)                     //input is only NaN
            return newData;
        for (int iRange=0; iRange<rangeLimits.length/2; iRange++) {
            int rStart = rangeLimits[2*iRange];
            int rEnd = rangeLimits[2*iRange+1];
            int newStart = (int)Math.ceil((oldX[rStart]-newX[0])/newStep - 1e-10);
            int newEnd = (int)Math.floor((oldX[rEnd-1]-newX[0])/newStep + 1e-10) + 1;
            if (newStart < 0) newStart = 0;
            if (newEnd > newX.length) newEnd = newX.length;
            if (newEnd - newStart < 1) continue;    //nothing to do in output array
            if (rEnd - rStart >= 2)                 //spline needs at least two points
                splineInterpolate(oldX, newX, oldData, newData, rStart, rEnd, newStart, newEnd);
            else
                Arrays.fill(newData, newStart, newEnd, oldData[rStart]);
        }
        return newData;
    }

    /** Spline interpolation, where valid indices of old quantities are firstValid to endValid-1 and
     *  output points in newData are written between newStart to newEnd-1.  */
    static void splineInterpolate(double[] oldX, double[] newX, double[] oldData, double[] newData,
            int firstValid, int endValid, int newStart, int newEnd) {
        float[] oldXF = extractFloat(oldX, firstValid, endValid);
        float[] oldDF = extractFloat(oldData, firstValid, endValid);
        SplineFitter spline = new SplineFitter(oldXF, oldDF, endValid-firstValid);
        for (int i=newStart; i<newEnd; i++) {
            double x = newX[i];
            newData[i] = spline.evalSpline(x);
        }
    }

    /** Extracts a subarray between 'start (inclusive) and 'end' (exclusive)
     *  as a float array */
    static float[] extractFloat(double[] in, int start, int end) {
        float[] out = new float[end - start];
        for (int iOut=0, iIn=start; iOut<out.length; iIn++,iOut++)
            out[iOut] = (float)in[iIn];
        return out;
    }

}
