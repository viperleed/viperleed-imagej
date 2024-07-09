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
 *  This class contains I(V) data:
 *  Energies, intensities and the name for each intensity column
 * 
 *  The class includes methods to read LeedIVData from a file and
 *  write them to a file, as well as methods related to the spot names,
 *  limiting the energy range, and deleting beams.
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
 *  @author Michael Schmid, IAP/TU Wien, 2023-2024
 */


public class LeedIVData implements Cloneable {
    /** Flag: ensure that the energy axis is all positive, equally spaced, and ascending */
    public static final int E_ASCENDING = 1;
    /** Flag: zero values at start or end of a data column are NaN */
    public static final int ZERO_IS_NAN = 2;
    /** The array with the energies */
    public double[] energies;
    /** The arrays with the intensities for each column (=beam) */
    public double[][] data;
    /** The array with the spot names corresponding to the data columns (=beams) */
    public String[] spotNames;
    /** The array with the indices in the spotPattern corresponding to the data columns (beams). 
     *  Can be null (if no spotPattern was supplied). */
    public int[] spotIndices;
    /** Typically the name of the file or ResultsTable where the data were read */
    public String title;

    /** Creates a LeedIVData object from the given arrays, where the columns of 'data'
     *  and the spotNames must have a 1:1 correspondence. 'spotIndices' may be null. */
    public LeedIVData(String title, double[] energies, double[][] data, String[] spotNames, int[] spotIndices) {
        this.title = title;
        this.energies = energies;
        this.data = data;
        this.spotNames = spotNames;
        this.spotIndices = spotIndices;
    }

    /** Creates and returns a LeedIVData object read from a .csv file.
     *  If the spotPattern is given, also the spotIndices array will be populated with the
     *  indices of the spots for the data columns (beams).
     *  Displays an error message (IJ.error) and returns null if not successful. */
    public static LeedIVData fromFile(String filePath, LeedSpotPattern spotPattern, int flags) {
        String str = IJ.openAsString(filePath);
        if (str==null) return null;                 //cancel
        if (str.startsWith("Error:")) {    // failed to open
            IJ.error("I(V) Data", str+"\n"+filePath);
            return null;
        }

        LeedIVData ivData = null;
        String title = (new File(filePath)).getName();
        if (title.toLowerCase().endsWith(".csv"))
            title = LeedUtils.removeExtension(title);

        String[] lines = Tools.split(str, "\n");
        int n = -1;                                 //number of data lines, start with minus one for the header
        for (int l=0; l<lines.length; l++) {
            lines[l] = lines[l].trim();
            if (LeedUtils.isCommentLine(lines[l]))
                lines[l] = null;
            else
                n++;
        }
        boolean waitForHeader = true;               //the first non-comment line will be the header
        int iData=0;
        int nCol = -1;
        int energyColumn = -1;
        for (int l=0; l<lines.length; l++) {
            String line = lines[l];
            if (line == null) continue;
            String[] items = line.split(",");
            if (waitForHeader) {                    //header line found
                energyColumn = LeedUtils.getEnergyColumn(items);
                if (energyColumn < 0) {
                    IJ.error("I(V) Data", "ERROR: No energy column in "+title);
                    return null;
                }
                nCol = items.length-1;
                String[] spotNames = new String[nCol];
                int ic = 0;
                for (int i=0; i<items.length; i++)
                    if (i != energyColumn)
                        spotNames[ic++] = items[i];

                double[][] data = new double[nCol][n];
                for (double[] col : data)
                    Arrays.fill(col, Double.NaN);
                double[] energies = new double[n];
                ivData = new LeedIVData(title, energies, data, spotNames, null);
                waitForHeader = false;
            } else {                            //normal data line
                if (items.length <= energyColumn) {
                    IJ.error("I(V) Data", "Line "+(l+1)+" in "+title+" too short:\n"+line);
                    return null;
                }
                int ic = 0;
                for (int i=0; i<items.length; i++)
                    if (i == energyColumn)
                        ivData.energies[iData] = Tools.parseDouble(items[i].trim());
                    else
                        ivData.data[ic++][iData] = Tools.parseDouble(items[i].trim());
                iData++;
            }
        }
        double[] firstAndInc = LeedUtils.getFirstAndIncrement(ivData.energies);
        if (LeedUtils.flagSet(flags, E_ASCENDING)) {
            if (!(ivData.energies[0] > 0) || !(firstAndInc[1] > 0)) {
                IJ.error("I(V) Data", "ERROR: Energies in "+title+" are not positive, evenly spaced and ascending");
                return null;
            }
        } else if (!isMonotonous(ivData.energies)) {
                IJ.error("I(V) Data", "ERROR: Energies in "+title+" are not monotonously ascending or descending");
                return null;            
        }
        if (LeedUtils.flagSet(flags, ZERO_IS_NAN))
            for (double[] d : ivData.data)
                convertZeroToNaN(d);

        if (spotPattern != null) {
            ivData.setSpotPattern(spotPattern, /*showWarning=*/true);
            if (ivData.spotIndices == null)
                return null;
        }
        return ivData;
    }

    /** Creates the spotIndices array, i.e., the array of indices in the
     *  spotPattern for each column (beam) of the data.
     *  When called with a null argument, the spotPattern is created from the
     *  spotNames array.
     *  If there are beams not found in the spotPattern provided as an argument,
     *  the spotIndices array is set to null and an error message displayed. 
     *  'showWarning' determines whether to show a warning if there is a disagreement
     *  between the current groups and the spot pattern file. */
    public void setSpotPattern(LeedSpotPattern spotPattern, boolean showWarning) {
        if (spotPattern == null)
            spotPattern = new LeedSpotPattern(spotNames, false);
        spotIndices = spotPattern.getIndices(spotNames);
        String error = spotPattern.getErrorMessage();
        String sourceName = title==null ? "" : " in "+title+"\n";
        if (LeedUtils.countNonPositives(spotIndices) > 0) {
            IJ.error("I(V) data", "ERROR: Spots"+sourceName+" missing in spot pattern file"+(error == null ? "" : ("\n"+error)));
            spotIndices = null;
            return;
        }
        if (showWarning && error != null)
            IJ.error("I(V) Data", "Warning: Disagreement with spot pattern file\n"+sourceName+error);
    }

    /** Saves this IVData as a comma delimited text file.
     *  Displays an error message and returns 'false' if there is an error.
     *  Values of 0 are written as 2e-45 (the smallest 32-bit float number)
     *  since some LEED programs treat 0 as NaN.
     *  Writes an empty field if the data do not exist or are NaN */
    public boolean save(String path) {
        if (spotNames.length != data.length)
            throw new RuntimeException("Number of headings different from number of data arrays: "+spotNames.length+" != "+data.length);
        int energyDigits = LeedUtils.getEnergyDigits(energies);
        PrintWriter pw = null;
        try {
            FileOutputStream fos = new FileOutputStream(path);
            BufferedOutputStream bos = new BufferedOutputStream(fos);
            pw = new PrintWriter(bos);
            pw.print("E");
            for (int ic=0; ic<spotNames.length; ic++) {
                pw.print(",");
                pw.print(spotNames[ic]);
            }
            pw.println();

            for (int i=0; i<energies.length; i++) {
                pw.print(IJ.d2s(energies[i], energyDigits));
                for (int ic=0; ic<data.length; ic++) {
                    pw.print(',');
                    if (data[ic] == null || Double.isNaN(data[ic][i])) continue;
                    double v = data[ic][i];
                    if ((float)v == 0f)             //avoid zero, TensErLEED may read 0 as "no valid data"
                        v = Math.copySign(2e-38, v==0 ? 1 : v);
                    pw.print(LeedUtils.d2s(v,7));   //7 significant digits
                }
                pw.println();
            }
            pw.close();
            return true;
        } catch (Exception e) {
            IJ.error("Error writing file "+path+"\n"+e);
            if (pw != null)
                try {pw.close();} catch (Exception e2) {}
            return false;
        }
    }

    /** Sets the names according the the 'canonical' name in the spotPattern,
     *  with the a vertical bar '|' between the h and k indices, not a comma,
     *  and including the spot group in square brackets. */
    public void setCanonicalNames(LeedSpotPattern spotPattern) {
        for (int ic=0; ic<spotIndices.length; ic++) {
            int spotIndex = spotIndices[ic];
            if (spotIndex >= 0)
                spotNames[ic] = spotPattern.getNameWithGroup(spotIndex, true);
        }
    }

    /** Deletes all columns (beams) where the spotIndex is the given value */
    public void deleteIfSpotIndex(int spotIndexToDelete) {
        int nDelete = 0;
        for (int i=0, iNew=0; i<spotIndices.length; i++)
            if (spotIndices[i] == spotIndexToDelete)
                nDelete++;
        if (nDelete == 0) return;
        int nGoodIndices = spotIndices.length - nDelete;
        double[][] newData = new double[nGoodIndices][];
        String[] newSpotNames = new String[nGoodIndices];
        int[] newSpotIndices = new int[nGoodIndices];
        for (int i=0, iNew=0; i<spotIndices.length; i++) {
            if (spotIndices[i] != spotIndexToDelete) {
                newSpotNames[iNew] = spotNames[i];
                newData[iNew] = data[i];
                newSpotIndices[iNew] = spotIndices[i];
                iNew++;
            }
        }
        data = newData;
        spotNames = newSpotNames;
        spotIndices = newSpotIndices;
    }

    /** Deletes energies at the start and end of the energy range,
     *  if there are no valid (non-NaN) data at this energy.
     *  Returns false (and does not modify the data) if all data are NaN. */
    public boolean trimEnergies() {
        int start = Integer.MAX_VALUE, end = 0;
        for (int ic=0; ic<data.length; ic++) {
            if (data[ic] == null) continue;
            for (int i=0; i<Math.min(start, data[ic].length); i++) {
                if (!Double.isNaN(data[ic][i])) {
                    start = i;                  //this column starts earlier
                    break;
                }
            }
            for (int i=data[ic].length-1; i>=end; i--) {
                if (!Double.isNaN(data[ic][i])) {
                    end = i+1;                  //this column ends later
                    break;
                }
            }
        }
        return cropEnergyRange(start, end);
    }

    /** Reduces the energy range to the limits given by the two array elements.
     *  The first must be the lower limit. If the energy is equal to any limit it is included.
     *  Returns false (and does not modify the data) if no energies would remain. */
    public boolean cropEnergyRange(double[] range) {
        int start = Integer.MAX_VALUE, end = energies.length;
        for (int i=0; i<energies.length; i++) {
            if (energies[i] >= range[0]) {
                start = i;
                break;
            }
        }
        for (int i=energies.length-1; i>=0; i--) {
            if (energies[i] > range[1])
                end = i;
            else
                break;
        }
        return cropEnergyRange(start, end);
    }

    /** Reduces the data range to indices start (inclusive) to end (exclusive) */
    private boolean cropEnergyRange(int start, int end) {
        if (start >= end) return false;
        if (start == 0 && end == energies.length) return true;
        energies = Arrays.copyOfRange(energies, start, end);
        for (int ic=0; ic<data.length; ic++)
            data[ic] = Arrays.copyOfRange(data[ic], start, end);
        return true;
    }

    /** Deletes data columns that are null or contain only NaN values.
     *  Returns false if there are no non-NaN at all.
     *  In this case, no modifications are done. */
    public boolean trimData() {
        boolean noSpotIndices = spotIndices == null;
        if (noSpotIndices) spotIndices = new int[data.length];
        boolean goodDataFound = false;
        for (int ic=0; ic<data.length; ic++)
            if (data[ic] == null || LeedUtils.countNonNaN(data[ic]) == 0)
                spotIndices[ic] = Integer.MIN_VALUE;
            else
                goodDataFound = true;
        if (!goodDataFound) return false;
        deleteIfSpotIndex(Integer.MIN_VALUE);
        if (noSpotIndices) spotIndices = null;
        return true;
    }

    /** Sorts the entries in the sequence of ascending spotIndices. */
    public void sort() {
        if (LeedUtils.isSorted(spotIndices)) return;
        int[] sortedIndices = LeedUtils.sortedIndices(spotIndices);
        LeedIVData copy = duplicate();
        for (int ic=0; ic<spotIndices.length; ic++) {
            data[ic] = copy.data[sortedIndices[ic]];
            spotNames[ic] = copy.spotNames[sortedIndices[ic]];
            spotIndices[ic] = copy.spotIndices[sortedIndices[ic]];
        }
    }

    /** Changes zero values at the beginning and end of the array to NaN */
    static void convertZeroToNaN(double[] a) {
        for (int i=0; i<a.length; i++) {
            if (a[i] == 0)
                a[i] = Double.NaN;
            else
                break;
        }
        for (int i=a.length-1; i>0; i--) {
            if (a[i] == 0)
                a[i] = Double.NaN;
            else
                break;
        }
    }

    /** Returns a shallow clone; no array contents are cloned.
     *  Thus, no not modify arrays of a clone, only replace them (at top level
     *  for the data[][] array). The duplicate() function returns a somewhat
     *  deeper clone (see there for details). */
    public Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

    /** Returns a shallow clone; no array contents are cloned.
     *  Thus, no not modify arrays of a clone, only replace them (at top level
     *  for the data[][] array). The duplicate() function returns a somewhat
     *  deeper clone (see there for details). */
    public LeedIVData shallowClone() {
        try {
            return (LeedIVData)this.clone();
        } catch (CloneNotSupportedException e) {return null;}
    }

    /** Returns a copy of this LeedIVData object.
     *  Note that the data subarrays and the energy array are not duplicated;
     *  any changes to them also affects the original! */
    public LeedIVData duplicate() {
        LeedIVData ivData = (LeedIVData)this.shallowClone();
        ivData.data = (double[][])data.clone();
        ivData.spotNames = (String[])spotNames.clone();
        ivData.spotIndices = (int[])spotIndices.clone();
        return ivData;
    }

    /** Returns a copy of this LeedIVData object.
     *  Only the 'data' arrays are duplicated (deep clone),
     *  all other arrays are those of the original and must NOT
     *  be modified in the copy */
    public LeedIVData duplicateData() {
        LeedIVData ivData = (LeedIVData)this.shallowClone();
        ivData.data = LeedUtils.duplicate(data);
        return ivData;
    }

    /** Returns whether the array is strictly monotonous (increasing or decreasing).
     *  An array with less than two elements is considered monotonous. */
    static boolean isMonotonous(double[] energies) {
        if (energies.length < 2)
            return true;
        double stepSign = Math.signum(energies[1] - energies[0]);
        if (!(Math.abs(stepSign) == 1))
            return false;
        for (int i=2; i<energies.length; i++)
            if (!(Math.signum(energies[i] - energies[i-1]) == stepSign))
                return false;
        return true;
    }

    /** Returns an ImageJ ResultsTable for this LeedIVData.
     *  If 'zeroReplacement' is nonzero, values with an absolute value less
     *  than this number are replaced by this number (with the correct sign,
     *  if nonzero). This should avoid reading values of zero as undefined
     *  by some LEED programs. */
    /* UNUSED
    public ResultsTable toResultsTable(double zeroReplacement) {
        if (spotNames.length != data.length)
            throw new RuntimeException("Number of headings different from number of data arrays: "+spotNames.length+" != "+data.length);
        ResultsTable rt = new ResultsTable(energies.length);

        for (int i=0; i<energies.length; i++)
            rt.setValue("E", i, energies[i]);

        for (int col=0; col<data.length; col++) {
            String head = spotNames[col];
            int colIndex = -1;  //will become the column index in the table
            for (int i=0; i<data[col].length; i++) {
                double value = data[col][i];
                if (zeroReplacement > 0 && Math.abs(value)<zeroReplacement)
                    value = Math.copySign(zeroReplacement, value+Double.MIN_VALUE);

                if (colIndex < 0) {
                    rt.setValue(head, i, value);
                    colIndex = rt.getColumnIndex(head);
                } else
                    rt.setValue(colIndex, i, value); //faster not to look up the column index each time
            }
        }
        return rt;
    } */

}
