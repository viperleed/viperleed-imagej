import ij.IJ;
import java.util.Arrays;
import java.util.ArrayList;

/**
 *  This class averages over LEED I(V) curves, with a smooth fade-in and fade-out
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


public class LeedCurveAverager {
    /* We blend over this energy range (in V_0i units): */
    public static final double BLEND_OVER_V0I = LeedCurvePlotEditor.BLEND_OVER_V0I;
    /* If intensity ratio between two beams is larger than the following, beams averaging does not take the ratio into account */
    static final double MAX_EQUIV_I_RATIO = LeedCurvePlotEditor.MAX_EQUIV_I_RATIO; 

    /** Calculates the average over equivalent beams, with a soft transition
     *  (fade in/fade out) where data for beams end, for use in the energy
     *  range given by indices iEStart (inclusive) to iEEnd (exclusive).
     *  Data outside that range are NOT set to NaN.
     *  Requires 0 <= iEStart and iEEnd <= length of data arrays in inData.
     *  The minimum overlap (in data points) is 'minOverlap'.
     *  Returns null if no input data are present. */
    public static double[] calculateAverage(ArrayList<double[]> inData, int iEStart, int iEEnd, int minOverlap, double v0iOverStep) {
        inData = new ArrayList<double[]>(inData); //copy the original list since we will remove the elements from the list
        int length = -1;
        /* eliminate unusable input arrays */
        for (int i=inData.size()-1; i>=0; i--) {
            if (inData.get(i) == null) {
                inData.remove(i);
                continue;
            }
            if (length < 0)
                length = inData.get(i).length;
            else if (inData.get(i).length != length)
                throw new IllegalArgumentException("Arrays for averaging have different lengths: "+
                        inData.get(i).length+" vs "+length);
            int[] startEnd = LeedUtils.getStartEndNum(inData.get(i));
            if (startEnd[0] < 0 || startEnd[0]>=iEEnd || startEnd[1]<=iEStart) {
                inData.remove(i);
                continue;
            }
        }
        if (inData.isEmpty()) return null;

        /* find the pair with the largest overlap, or if none, the curve with the longest contiguous range */
        int[] bestOverlapColumnIndices = bestOverlapIndices(inData, iEStart, iEEnd, minOverlap);
        if (bestOverlapColumnIndices == null) {
            int bestColumnIndex = bestRangeIndex(inData, iEStart, iEEnd);
            if (bestColumnIndex >= 0)
                return inData.get(bestColumnIndex);
            else
                return null;
        }
        int nCol = inData.size();
        /* merge the data and set those merged in to null in inData */
        double[] sumData = (double[])inData.get(bestOverlapColumnIndices[0]).clone();
        double[] weights = new double[length];
        for (int i=0; i<length; i++)
            if (!Double.isNaN(sumData[i])) weights[i] = 1.0;
        if (IJ.debugMode) IJ.log("start average with #"+bestOverlapColumnIndices[0]+" E indices="+iEStart+"-"+iEEnd);
        inData.set(bestOverlapColumnIndices[0], null);      //mark as used    
        int nextIndexToProcess = bestOverlapColumnIndices[1];
        int blendLength = (int)Math.round(BLEND_OVER_V0I*v0iOverStep+1e-6);
        do {
            if (IJ.debugMode) IJ.log("merge with #"+nextIndexToProcess);
            mergeSmoothly(sumData, weights, inData.get(nextIndexToProcess), blendLength, iEStart, iEEnd);
            inData.set(nextIndexToProcess, null);           //mark as used
            nextIndexToProcess = bestOverlapIndex(inData, sumData, iEStart, iEEnd, minOverlap);
        } while (nextIndexToProcess >= 0);
        for (int i=0; i<length; i++)
            sumData[i] /= weights[i];
        return sumData;
    }

    /** Sums 'sumData' and 'newData' into 'sumData'. It is assumed that each point of 'sumData' contains 'weights'
     *  entries (non-integer in the transition zone). Each item of 'weights' is increased by the weight of the point
     *  added.
     *  At the ends (if different), uses a soft transition over a range blendLength */
    static void mergeSmoothly(double[] sumData, double[] weights, double[] newData, int blendLength, int iEStart, int iEEnd) {

        int first0=-1, first1=-1;                                   //first valid data for sumData (0) and newData (1)
        for (int iE=0; iE<=sumData.length; iE++) {
            boolean valid0 = iE<sumData.length && !Double.isNaN(sumData[iE]);
            boolean valid1 = iE<sumData.length && !Double.isNaN(newData[iE]);
            if (valid0)
                if (first0 < 0) first0=iE;
            if (valid1)
                if (first1 < 0) first1=iE;
            if ((!valid0 || !valid1) && first0>=0 && first1>=0) {    //end of common data (or end of all data)
                if (iE >= iEStart) {
                    int firstCommon = Math.max(first0, first1);
                    boolean blend0Start = first0 > first1 && first0 > iEStart;
                    boolean blend1Start = first1 > first0 && first1 > iEStart;
                    boolean blend0End = !valid0 && valid1 && iE < iEEnd;
                    boolean blend1End = valid0 && !valid1 && iE < iEEnd;
                    boolean blendAtStart = blend0Start || blend1Start;
                    boolean blendAtEnd = blend0End || blend1End;
                    int nOverlap = iE - firstCommon;                //number of points in overlap region
                    boolean useTrend =  nOverlap > 5*blendLength;   //for long overlap ranges (20*V0i, ~100 eV), adjust ratio at the ends separately if we blend both
                    if (IJ.debugMode) IJ.log("end common @index "+iE+"... blendAtStart="+blendAtStart+" blendAtEnd="+blendAtEnd+" useTrend="+useTrend+" overlap="+nOverlap);
                    double newOverAvgStart = blendAtStart && useTrend ?  //at start & end, get ratio over two blend lengths (8*V0i)
                            getNewOverAvgRatio(sumData, weights, newData, firstCommon, firstCommon+2*blendLength) :
                            Double.NaN;
                    double newOverAvgEnd = blendAtEnd && useTrend ?
                            getNewOverAvgRatio(sumData, weights, newData, iE-2*blendLength, iE) :
                            Double.NaN;
                    double newOverAvgAll = Double.isNaN(newOverAvgStart) || Double.isNaN(newOverAvgEnd) ?
                            getNewOverAvgRatio(sumData, weights, newData, firstCommon, iE) : //ratio in complete overlap range
                            Double.NaN;
                    double slopeOfLog = 0;                          //log(new/old) slope vs energy index
                    double offsetOfLog = Double.NaN;                //log(new/old) at index 'firstCommon'
                    String where="neither";
                    if (!Double.isNaN(newOverAvgStart) && !Double.isNaN(newOverAvgEnd)) {
                        slopeOfLog = Math.log(newOverAvgEnd/newOverAvgStart)/(nOverlap-2*blendLength);
                        offsetOfLog = Math.log(newOverAvgStart) - slopeOfLog*(blendLength-0.5);
                        where="both";
                    } else if (!Double.isNaN(newOverAvgStart)) {
                        slopeOfLog = Math.log(newOverAvgAll/newOverAvgStart)/(0.5*(nOverlap-blendLength));
                        offsetOfLog = Math.log(newOverAvgStart) - slopeOfLog*(blendLength-0.5);
                        where="start";
                    } else if (!Double.isNaN(newOverAvgEnd)) {
                        slopeOfLog = Math.log(newOverAvgEnd/newOverAvgAll)/(0.5*(nOverlap-blendLength));
                        offsetOfLog = Math.log(newOverAvgAll) - slopeOfLog*(0.5*(nOverlap-1));
                        where="end";
                    } else {
                        offsetOfLog = Math.log(newOverAvgAll);
                    }
                    if (Double.isNaN(offsetOfLog) && nOverlap > 0) {
                        slopeOfLog = 0;
                        offsetOfLog = 0;
                    }
                    if (IJ.debugMode) IJ.log("ratio ok at "+where+". New/avg offsetOfLog="+(float)offsetOfLog+" slopeOfLog="+(float)slopeOfLog);
                    if (!Double.isNaN(offsetOfLog)) {               //apply the correction and average
                        double sumW0 = 0, sumW1 = 0;                //get average weights in overlap region
                        for (int i=firstCommon; i<iE; i++) {
                            double w0 = weights[i];
                            double w1 = 1.0;                        //these weights get reduced in the blending zone
                            if (i-firstCommon<blendLength && blend0Start)
                                w0 *= (i-firstCommon+1)*(1.0/(blendLength+1));
                            if (i-firstCommon<blendLength && blend1Start)
                                w1 *= (i-firstCommon+1)*(1.0/(blendLength+1));
                            if (iE-i <= blendLength && blend0End)
                                w0 *= (iE-i)*(1.0/(blendLength+1));
                            if (iE-i <= blendLength && blend1End)
                                w1 *= (iE-i)*(1.0/(blendLength+1));
                            sumW0 += w0;
                            sumW1 += w1;
                        }
                        double rWeight0 = sumW0 / (sumW0 + sumW1);   //relative average weights of 'sumData' aka 0
                        double rWeight1 = sumW1 / (sumW0 + sumW1);   //... and 'newData' aka 1 in the overlap zone
                        double factor0 = Math.exp(offsetOfLog*rWeight1);  //correction factors at the start of the overlap zone
                        double factor1 = Math.exp(-offsetOfLog*rWeight0);
                        for (int i=0; i<firstCommon; i++) {         //correct or supply data before the common range
                            if (Double.isNaN(sumData[i])) {
                                if (!Double.isNaN(newData[i])) {
                                    sumData[i] = newData[i] * factor1;
                                    weights[i] = 1.0;
                                }
                            } else
                                sumData[i] *= factor0;
                        }
                        if (IJ.debugMode) IJ.log("factors start: "+(float)factor0+","+(float)factor1+" log f0/f1="+(float)Math.log(factor0/factor1));
                        for (int i=firstCommon; i<iE; i++) {
                            double w0 = 1.0;
                            double w1 = 1.0;                        //these get reduced in the blending zone
                            if (i-firstCommon<blendLength && blend0Start)
                                w0 *= (i-firstCommon+1)*(1.0/(blendLength+1));
                            if (i-firstCommon<blendLength && blend1Start)
                                w1 *= (i-firstCommon+1)*(1.0/(blendLength+1));
                            if (iE-i <= blendLength && blend0End)
                                w0 *= (iE-i)*(1.0/(blendLength+1));
                            if (iE-i <= blendLength && blend1End)
                                w1 *= (iE-i)*(1.0/(blendLength+1));
                            double norm = 1./(w0*weights[i] + w1);
                            rWeight0 = w0*weights[i] * norm; rWeight1 = w1 * norm;
                            factor0 *= Math.exp(slopeOfLog*rWeight1);
                            factor1 *= Math.exp(-slopeOfLog*rWeight0);
                            sumData[i] = w0*factor0*sumData[i] + w1*factor1*newData[i];
                            weights[i] = w0*weights[i] + w1;
                        }
                        if (IJ.debugMode) IJ.log("factors end: "+(float)factor0+","+(float)factor1+" log f0/f1="+(float)Math.log(factor0/factor1));
                        boolean moreSumData = false;
                        for (int i=iE; i<sumData.length; i++) {     //correct or supply data beyond the common range
                            if (!Double.isNaN(sumData[i])) {
                                sumData[i] *= factor0;
                                moreSumData = true;
                            } else if (!moreSumData) {
                                if (!Double.isNaN(newData[i])) {
                                    sumData[i] = newData[i] * factor1;
                                    weights[i] = 1.0;
                                    iE = i;                         //do not apply further processing to the new data
                                }
                            }
                        }
                    }
                }
                if (iE > iEEnd) break;
                if (!valid0) first0 = -1;
                if (!valid1) first1 = -1;
            } //at end of common data
        } // for iE
    }

    /** Returns the ratio, as determined from the average of newData and oldData/weights
     *  in the region given by the indices start-end.
     *  Returns NaN if no valid points or the ratio would be outside the [1/3, 3] range */
    static double getNewOverAvgRatio(double[] sumData, double[] weights, double[] newData, int start, int end) {
        if (start < 0) start = 0;
        if (end > sumData.length) end = sumData.length;
        boolean validDataFound = false;
        for (int i=start; i<end; i++) {
            if (Double.isNaN(sumData[i]) || Double.isNaN(newData[i])) {
                if (validDataFound) end = i;
                else start = i;
            } else
                validDataFound = true;
        }
        if (end <= start) return Double.NaN;
        double sumOld = 0, sumNew = 0;
        for (int i=start; i<end; i++) {
            sumNew += newData[i];
            sumOld += sumData[i] / weights[i];
        }
        double ratio = sumNew/sumOld;
        if (ratio > MAX_EQUIV_I_RATIO || ratio < 1./MAX_EQUIV_I_RATIO)
            ratio = Double.NaN;    //a high factor is unlikely, probably bad data
        return ratio;
    }

    /** Returns the index among the data array in 'inData' given with the longest
     *  contiguous range of data in the energy range between indices iEStart, iEEnd.
     *  Returns -1 if there are no valid data at all. */
    static int bestRangeIndex(ArrayList<double[]> inData, int iEStart, int iEEnd) {
        int bestArrayIndex = -1;
        int bestOverlap = -1;
        for (int ic=0; ic<inData.size(); ic++) {
            double[] icData = inData.get(ic);
            if (icData == null) continue;
            int nOverlap = 0;
            for (int i=iEStart; i<=iEEnd; i++) {
                boolean valid = i<iEEnd ? !Double.isNaN(icData[i]) : false;
                if (!valid) {
                    if (nOverlap > bestOverlap) {
                        bestOverlap = nOverlap;
                        bestArrayIndex = ic;
                        nOverlap = 0;
                    }
                } else {
                    nOverlap++;
                }
            }
        }
        return bestArrayIndex;
    }

    /** Returns the index of the inData column with longest contiguous overlap to otherData
     *  in the energy range between indices iEStart, iEEnd.
     *  Columns with inData = null are ignored.
     *  Returns -1 if no data or no overlap */
    static int bestOverlapIndex(ArrayList<double[]> inData, double[] otherData, int iEStart, int iEEnd, int minOverlap) {
        int bestColIndex = -1;
        int bestOverlap = -1;
        for (int ic=0; ic<inData.size(); ic++) {
            double[] icData = inData.get(ic);
            if (icData == null) continue;
            int nOverlap = 0;
            for (int i=iEStart; i<=iEEnd; i++) {
                boolean valid0 = i<iEEnd ? !Double.isNaN(icData[i]) : false;
                boolean valid1 = i<iEEnd ? !Double.isNaN(otherData[i]) : false;
                if (!valid0 || !valid1) {
                    if (nOverlap > bestOverlap) {
                        bestOverlap = nOverlap;
                        bestColIndex = ic;
                        nOverlap = 0;
                    }
                } else {
                    nOverlap++;
                }
            }
        }
        if (bestOverlap < minOverlap)
            return -1;
        else
            return bestColIndex;
    }


    /** Returns the pair of column indices among those in 'inData' with the
     *  longest contiguous mutual overlap in the energy range between indices iEStart, iEEnd.
     *  Returns or null if none or only one column, or no overlap of at least 'minOverlap' points */
    static int[] bestOverlapIndices(ArrayList<double[]> inData, int iEStart, int iEEnd, int minOverlap) {
        int[] bestPair = new int[] {-1, -1};
        int bestOverlap = -1;
        for (int ic=0; ic<inData.size()-1; ic++) {
            double[] icData = inData.get(ic);
            if (icData == null) continue;
            for (int jc = ic+1; jc<inData.size(); jc++) {
                double[] jcData = inData.get(jc);
                if (jcData == null) continue;
                int nOverlap = 0;
                for (int i=iEStart; i<=iEEnd; i++) {
                    boolean valid0 = i<iEEnd ? !Double.isNaN(icData[i]) : false;
                    boolean valid1 = i<iEEnd ? !Double.isNaN(jcData[i]) : false;
                    if (!valid0 || !valid1) {
                        if (nOverlap > bestOverlap) {
                            bestOverlap = nOverlap;
                            bestPair[0] = ic; bestPair[1] = jc;
                            nOverlap = 0;
                        }
                    } else {
                        nOverlap++;
                    }
                }
            }
        }
        if (bestOverlap < minOverlap)
            return null;
        else
            return bestPair;
    }

}
