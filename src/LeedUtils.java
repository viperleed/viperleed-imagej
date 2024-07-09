import ij.*;
import ij.ImagePlus;
import ij.process.*;
import ij.gui.Overlay;
import ij.util.Tools;
import ij.measure.ResultsTable;
import ij.gui.GenericDialog;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Arrays;
import java.util.Formatter;
import java.util.Locale;
import java.awt.Choice;
import java.awt.Frame;
import java.awt.event.*;
import java.io.File;
import java.io.FileReader;


/**
 *  This class contains various static utility methods
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

public class LeedUtils {
    /** The name of the command key on Macs, control otherwise */
    public static final String CTRL = IJ.isMacintosh() ? "\u2318" : "CTRL";

    /** Removes the extension from a filename */
    public static String removeExtension(String filename) {
        int i = filename.lastIndexOf('.');
        if (i>0 && i>filename.length()-5)
            filename = filename.substring(0, i);
        return filename;
    }

    /** Returns the current date & time in a given format, e.g. "yyyy-MM-dd HH:mm:ss" */
    public static String getDateFormatted(String format) {
        SimpleDateFormat formatter = new SimpleDateFormat(format);
        Date date = new Date();
        return formatter.format(date);
    }

    /** Converts a string to a format suitable for CSV:
     *  Replaces quote characters " with double "" and any sequence of return
     *  or line feed characters with backslash-n, then adds surrounding quotes */
    public static String toCsvString(String str) {
        str = str.replace("\"", "\"\"");
        str = str.replaceAll("[\\r\\n]+","\\\\n");
        return '"' + str + '"';
    }

    /** Takes a string as read from CSV (without quotes! these are already removed
     *  when reading a csv file into a ResultsTable) and converts it to a normal string
     *  Replaces inner double quote characters "" with two single ones " and the
     *  backslash-n sequence with a linefeed. */
    public static String fromCsvString(String str) {
        str = str.replace("\"\"", "\"");
        str = str.replace("\\n","\n");
        return str;
    }

    /** Returns how often a given character occurs in a String */
    public static int countInString(String str, char c) {
        int count = 0;
        for (int i=0; i<str.length(); i++)
            if (str.charAt(i) == c)
                count++;
        return count;
    }

    /** Returns the characters that str1, str2 have in common and at the same position
     *  Non-common inside sections are replaced by a single 'X'. */
    public static String getCommon(String str1, String str2) {
        StringBuilder sb = new StringBuilder();
        for (int i=0; i<Math.min(str1.length(), str2.length()); i++)
            if (str1.charAt(i) == str2.charAt(i))
                sb.append(str1.charAt(i));
            else
                sb.append((char)0);
        String str = sb.toString().trim();
        str = str.replaceAll("\\000+","X"); //replace multiple non-common characters with a single 'X'
        
        return str;
    }

    /** Returns the indices of a spot designation, with or without round brackets,
     *  with or without additional group number in square brackets. Allowed formats are, e.g.
     *  (1|0) [1] or (1/2;3/2) or 0.3333,-1.333[5].
     *  Returns null on error. */
    public static double[] spotIndices(String spotName) {
        String[] hAndK = spotName.replaceAll("[()]", "").replaceAll("\\[.*\\]","").trim().split("[,;_|\\s]+");
        if (hAndK.length != 2) return null;
        double h = LeedUtils.toNumber(hAndK[0]);
        double k = LeedUtils.toNumber(hAndK[1]);
        if (Double.isNaN(h) || Double.isNaN(k)) return null;
        else return new double[] {h, k};
    }

    /** Converts a string with a fraction like 2/3 or a decimal number to a double.
     *  Returns NaN in case of failure. */
    public static double toNumber(String str) {
        str = str.trim();
        if (str.indexOf('/') > 0) {
            String[] strings = Tools.split(str,"/");
            if (strings.length != 2)
                return Double.NaN;
            else
                return Tools.parseDouble(strings[0])/Tools.parseDouble(strings[1]);
        } else
            return Tools.parseDouble(str);
    }

    /** Converts a string of comma-delimited numbers to a double[] array */
    public static double[] numbersFromString(String str) {
        String[] parts = str.split(","); //(no trailing empty Strings)
        double[] out = new double[parts.length];
        for (int i=0; i<parts.length; i++)
            out[i] = Tools.parseDouble(parts[i]);
        return out;
    }

    /** Converts a string String like "80-200.5" to a two-element double[] array.
     *  Both numbers must be non-negative and the second number not smaller than the first one.
     *  Returns null on error */
    public static double[] rangeFromString(String str) {
        String[] parts = str.split("[-\u2013]", 3); // minus or en-dash (trailing empty Strings are kept)
        if (parts.length != 2) return null;
        double[] out = new double[2];
        for (int i=0; i<parts.length; i++) {
            out[i] = Tools.parseDouble(parts[i]);
            if (!(out[i] >= 0)) return null;
        }
        if (!(out[1] >= out[0])) return null;
        return out;
    }

    /** Converts a String with fractions or decimal numbers to a numeric array. Delimiters may be
     *  any sequence of spaces, tabs, or any of these characters: ,;|
     *  Leading or trailing delimiters are ignored.
     *  Returns an array with 0 elements if there is any invalid number or NaN */
    public static double[] toNumberArray(String str) {
        str = str.replaceAll("[\\s,;|]", " ").trim();
        String[] parts = Tools.split(str);
        double[] v = new double[parts.length];
        for (int i=0; i<parts.length; i++) {
            v[i] = toNumber(parts[i]);
            if (Double.isNaN(v[i])) return new double[0];
        }
        return v;
    }

    /** Reads the number (integer or floating with decimal point, not scientific format)
     *  immediately following one of the keywords.
     *  The keywords are tried in the sequence given, the first match is returned.
     *  Returns NaN if str is null, none of the keywords is found or the first keyword
     *  of the list encountered in the String is not followed by a valid number. */
    public static double getNumberFromString(String str, String[] keywords) {
        if (str == null) return Double.NaN;
        for (String key : keywords) {
            int i = str.indexOf(key);
            if (i >= 0) {
                int offset = i+key.length();
                if (offset >= str.length()) return Double.NaN;  //key at the end
                int nPoint=0;
                int p = offset;
                char c = str.charAt(p);
                if (!(isDigit(c) || c=='-')) return Double.NaN;  //does not start like a number
                while (++p < str.length()) {
                    c = str.charAt(p);
                    if (c == '.') nPoint++;
                    else if (!isDigit(c)) break;
                    if (nPoint > 1) break;
                }
                return Tools.parseDouble(str.substring(offset, p));
            }
        }
        return Double.NaN;
    }

    /** Locates and returns an integer in square brackets in the string, or Integer.MIN_VALUE if none */
    public static int getIntInBrackets(String str) {
        if (str == null) return Integer.MIN_VALUE;
        int braIndex = str.indexOf('[');
        int ketIndex = str.indexOf(']');
        if (braIndex >= 0 && ketIndex > braIndex+1) {
            str = str.substring(braIndex+1, ketIndex).trim();
            try {
            return Integer.parseInt(str);
            } catch (NumberFormatException e) {
            return Integer.MIN_VALUE;
            }
        } else
            return Integer.MIN_VALUE;
    }

    /** Returns a String representation of the number with a given number of significant digits.
     *  Uses scientific notation for large or small numbers.
     *  Integers are show without decimals unless scientific notation is used.
     *  When not using scientifc notation and the number rounded to an integer
     *  has more than dgigits than 'significantDigits', no further rounding is
     *  done, e.g. 123456.78 is shown as 123457, not 123000 if 'significantDigits'=3.
     *  The maximum String length is 'significantDigits'+5 for positive and
     *  'significantDigits'+7 for negative numbers with magnitude up to 9e+99 */
    public static String d2s(double x, int significantDigits) {
        boolean isInteger = x==(int)x;
        int digits = significantDigits - (int)Math.floor(Math.log10(Math.abs(x))+1e-10) - 1;
        if (digits < -5 || digits > significantDigits + 2 || digits > 9)
            return (new Formatter(null, Locale.US)).format("%."+significantDigits+"g", x).toString();
        else 
            return IJ.d2s(x, digits < 0 || isInteger ? 0 : digits);
    }

    /** Returns a String representation of the number with the accuracy of 32-bit floats */
    public static String d2s(double x) {
        if (x==(int)x)
            return IJ.d2s(x, 0);
        else
            return Float.toString((float)x);
    }

    /** Converts an array to a String of comma-delimited numbers (no comma at the end)
     *  with the accuracy of 32-bit floats */
    public static String toString(double[] a) {
        StringBuilder sb = new StringBuilder(a.length*20);
        for (int i=0; i<a.length; i++) {
            sb.append(d2s(a[i]));
            if (i < a.length-1)
                sb.append(",");
        }
        return sb.toString();
    }

    /** Reads an array from a String of comma-delimited values */
    public static double[] arrayFromString(String str) {
        String[] parts = Tools.split(str, ",");
        int length = parts.length;
        if (parts[length-1].trim().length() == 0)
            length--;       //nothing after last comma? then ignore it
        double[] a = new double[length];
        for (int i=0; i<length; i++)
            a[i] = Tools.parseDouble(parts[i]);
        return a;
    }

    /** Returns whether a line is a comment line starting with any of these characters: #%!*
     *  Null or empty lines also qualify as comment lines. */
    public static boolean isCommentLine(String line) {
        if (line == null || line.trim().length() == 0) return true;
        char c = line.charAt(0);
        if (c=='%' || c=='!' || c=='#' || c=='*') return true;
        return false;
    }

    static boolean isDigit(char c) {
        return c>='0' && c<='9';
    }

    /** Returns for how many charcters the begin of the two strings is equal */
    public static int nMatchingBeginCharacters(String str1, String str2) {
        if (str1 == null || str2 == null) return 0;
        int commonLength = Math.min(str1.length(), str2.length());
        for (int i=0; i<commonLength; i++)
            if (str1.charAt(i) != str2.charAt(i)) return i;
        return commonLength;
    }

    /** Returns for how many charcters the end of the two strings is equal */
    public static int nMatchingEndCharacters(String str1, String str2) {
        if (str1 == null || str2 == null) return 0;
        int commonLength = Math.min(str1.length(), str2.length());
        for (int i=0; i<commonLength; i++)
            if (str1.charAt(str1.length()-i-1) != str2.charAt(str2.length()-i-1)) return i;
        return commonLength;
    }

    /** Returns whether a number is integer within 1e-10 (absolute) numerical tolerance */
    public static boolean isInteger(double x) {
        return Math.abs(x - Math.round(x)) <= 1e-10;
    }

    /** Returns the last element of a double[] array */
    public static double lastElement(double[] a) {
        return a[a.length-1];
    }

    /** Returns the last element (last subarray) of a double[][] array */
    public static double[] lastElement(double[][] a) {
        return a[a.length-1];
    }

    /** Returns the first index of a given value in the array,
     *  or -1 if that value is not in the array */
    public static int arrayIndexOf(int[] a, int value) {
        for (int i=0; i<a.length; i++)
            if (a[i] == value)
                return i;
        return -1;
    }

    /** Returns the first index of a given object in the array,
     *  or -1 if that object is not in the array or null.
     *  Compares the object reference with 'equals', not '=='.
     *  The array may contain null values; these are ignored. */
    public static int arrayIndexOf(Object[] a, Object theObject) {
        for (int i=0; i<a.length; i++)
            if (a[i] != null && a[i].equals(theObject))
                return i;
        return -1;
    }

    /** Returns the index of the array element nearest to the given value.
     *  If 'minusIfOutside' is true and the value is further outside
     *  than 0.5 times the nearest step, returns -1.
     *  It the array is null or empty, also returns -1. */
    public static int getNearestIndex(double value, double[] a, boolean minusIfOutside) {
        if (a==null) return -1;
        double leastDistance = Double.MAX_VALUE;
        int iOfLeastDistance = -1;
        double nearestStep = 0;
        for (int i=0; i<a.length; i++) {
            if (Math.abs(value-a[i]) < leastDistance) {
                leastDistance = Math.abs(value-a[i]);
                iOfLeastDistance = i;
                if (i>0)
                    nearestStep = Math.abs(a[i]-a[i-1]);
                else if (i+1<a.length)
                    nearestStep = Math.abs(a[i+1]-a[i]);
            }
        }
        if (minusIfOutside && leastDistance > 0.5*nearestStep)
            iOfLeastDistance = -1;
        return iOfLeastDistance;
    }

    /** Returns the minimum of an integer array */
    public static int getMin(int[] array) {
        int min = Integer.MAX_VALUE;
        for (int i : array)
            if (i < min) min = i;
        return min;
    }

    /** Returns the maximum of an integer array */
    public static int getMax(int[] array) {
        int max = Integer.MIN_VALUE;
        for (int i : array)
            if (i > max) max = i;
        return max;
    }

    /** Returns the number of 'true' elements in a boolean array */
    public static int countTrue(boolean[] array) {
        int count = 0;
        for (boolean b : array)
            if (b) count++;
        return count;
    }

    /** For a boolean array, examines the elements with an index given by the int array, and
     *  returns the number of true booleans that the int array points to */
    public static int countTrue(boolean[] bArray, int[] indices) {
        int count = 0;
        for (int i : indices) {
            boolean b = bArray[i];
            if (b) count++;
        }
        return count;
    }

    /** Returns the number of non-NaN elements in the array */
    public static int countNonNaN(double[] array) {
        int count = 0;
        for (double x : array)
            if (!Double.isNaN(x)) count++;
        return count;
    }

    /** Returns the number of non-positive and infinite elements in the array */
    public static int countNonPositives(double[] array) {
        int count = 0;
        for (double x : array)
            if (!(x > 0) || Double.isInfinite(x)) count++;
        return count;
    }

    /** Returns the number of non-null entries in the array */
    public static int countNonNull(Object[] array) {
        int count = 0;
        for (Object o : array)
            if (o != null) count++;
        return count;
    }

    /** Returns the number of non-positive entries in the array */
    public static int countNonPositives(int[] array) {
        int count = 0;
        for (int i : array)
            if (i < 0) count++;
        return count;
    }

    /** Returns whether a number has the given flag set (strictly speaking,
     *  whether any of the bits given by the bitmask 'flag' is set) */
    public static boolean flagSet(int number, int flag) {
        return (number&flag) != 0;
    }

    /** Returns the sum over all array values */
    public static double getArraySum(double[] a) {
        double sum = 0;
        for (double x : a)
            sum += x;
        return sum;
    }

    /** Returns the maximum of the array values in the given range.
     *  Returns NaN if there are no usable data in the range */
    public static double getArrayMax(double[] a, int start, int end) {
        double max = -Double.MAX_VALUE;
        for (int i=start; i<end; i++)
            if (a[i] > max) max = a[i];
        return max == -Double.MAX_VALUE ? Double.NaN : max;
    }

    /** Returns the sum over the array values from 'beginIndex' (inclusive) to 'endIndex' (exclusive) */
    /* currently unused
    public static double getArraySum(double[] a, int beginIndex, int endIndex) {
        double sum = 0;
        for (int i=beginIndex; i<endIndex; i++)
            sum += a[i];
        return sum;
    }*/

    /** Returns whether an array is sorted in ascending sequence */
    public static boolean isSorted(int[] a) {
        for (int i=1; i<a.length; i++)
            if (a[i] < a[i-1]) return false;
        return true;
    }

    /** Returns the indices of the array in ascending sequence of the values */
    public static int[] sortedIndices(int[] a) {
        long[] sortMe = new long[a.length];
        for (int i=0; i<a.length; i++)
            sortMe[i] = ((long)a[i]) << 32 | i;
        Arrays.sort(sortMe);
        int[] out = new int[a.length];
        for (int i=0; i<a.length; i++)
            out[i] = (int)(sortMe[i] & 0xffffffff);
        return out;
    }

    /** Returns the 'indices' array sorted by the ascending values in the second array.
     *  Only the first 'length' values are read.
     *  Both input arrays must have length >= 'length'.
     *  The input arrays remain unchanged.
     *  If values are the same, sorts by ascending indices. */
    public static int[] sortedIndices(int[] indices, int[] values, int length) {
        long[] sortMe = new long[length];
        for (int i=0; i<length; i++)
            sortMe[i] = ((long)values[i]) << 32 | indices[i];
        Arrays.sort(sortMe);
        int[] out = new int[length];
        for (int i=0; i<length; i++)
            out[i] = (int)(sortMe[i] & 0xffffffff);
        return out;
    }

    /** Converts an array of Double objects to an array of the corresponding values */
    public static double[] doubleFromDouble(Object[] in) {
        double[] out = new double[in.length];
        for (int i=0; i<in.length; i++)
            out[i] = ((Double)in[i]).doubleValue();
        return out;
    }

    /** Returns the input array with the data outside the range set to NaN.
     *  Returns a new array if the input would be modified */
    public static double[] restrictRange(double[] data, int iEStart, int iEEnd) {
        boolean modificationRequired = false;
        for (int i=0; i<iEStart; i++)
            if (!Double.isNaN(data[i])) {
                modificationRequired = true;
                break;
            }
        if (!modificationRequired && iEEnd >= 0)
            for (int i=iEEnd; i<data.length; i++)
            if (!Double.isNaN(data[i])) {
                modificationRequired = true;
                break;
            }
        if (!modificationRequired)
            return data;
        data = (double[])data.clone(); //copy, don't modify the original
        for (int i=0; i<iEStart; i++)
            data[i] = Double.NaN;
        if (iEEnd >= 0)
            for (int i=iEEnd; i<data.length; i++)
                data[i] = Double.NaN;
        return data;
    }

    /** Returns the limits of valid data as an array {start, end, start, end, ...}
     *  where the valid (non-NaN) data of each range run from start to end-1.
     *  Returns a two-element array with -1 as start limit if zero range.
     *  Call with iEStart <=0 and iEEnd=-1 for the full range */
    public static int[] getRangeLimits(double[] data, int iEStart, int iEEnd) {
        if (iEStart < 0) iEStart = 0;
        if (iEEnd < 0) iEEnd = data.length;
        LeedIntegerArray limits = new LeedIntegerArray(2); //in most cases, there are only two limits, then we need not copy the array
        int rangeStart = -1;
        for (int i=iEStart; i<=iEEnd; i++) {
            boolean isNaN = i==iEEnd || Double.isNaN(data[i]);
            if (rangeStart < 0 && !isNaN)
                rangeStart = i;
            else if (rangeStart>=0 && isNaN) { //end of range of valid (non-NaN) data
                limits.add(rangeStart);
                limits.add(i); //range end (exclusive)
                rangeStart = -1;
            }
        }
        if (limits.size() > 0)
            return limits.getArray();
        else
            return new int[] {-1, 0};
    }

    /** For debug, converts the output of getRangeLimits to a String */
    public static String rangeLimitsString(int[] rangeLimits) {
        if (rangeLimits[0] < 0) return "Empty range";
        String str = "";
        for (int i=0; i<rangeLimits.length/2; i++)
            str += rangeLimits[2*i]+"-"+rangeLimits[2*i+1]+" ";
        return str;
    }

    /** Returns as a two-element array the first and last+1 index
     *  where the array is not NaN.
     *  If there are only NaN values, returns [-1, -1].
     *  The output array can be supplied; then the limits in it are extended to reflect the range
     *  where any data, the previous and the current ones, are not NaN. */
    public static int[] getStartEndNum(double[] array) {
        int[] startEnd = new int[] {-1, -1};
        for (int i=0; i<array.length; i++)
            if (!Double.isNaN(array[i])) {
                startEnd[0] = i;
                break;
            }
        if (startEnd[0] < 0) return startEnd;
        int checkLastFrom = Math.max(startEnd[0], startEnd[1]);
        for (int i=array.length-1; i>checkLastFrom; i--)
            if (!Double.isNaN(array[i])) {
                startEnd[1] = i+1;
                break;
            }
        return startEnd;
    }

    /** Creates an array with the sequence from 0 to length-1 */
    public static double[] createSequence(int length) {
        double[] out = new double[length];
        for (int i=0; i<length; i++)
            out[i] = i;
        return out;
    }

    /** Creates a copy of a 2D double[][] array */
    public static double[][] duplicate(double[][] a) {
        double[][] out = new double[a.length][];
        for (int i=0; i<a.length; i++)
            out[i] = a[i].clone();
        return out;
    }

    /** Converts a byte[] or short[] array to float[] */
    public static float[] toFloat(Object array) {
        if (array instanceof byte[]) {
            byte[] a = (byte[])array;
            float[] out = new float[a.length];
            for (int i=0; i<a.length; i++)
                out[i] = a[i] & 0xff;
            return out;
        } else {
            short[] a = (short[])array;
            float[] out = new float[a.length];
            for (int i=0; i<a.length; i++)
                out[i] = a[i] & 0xffff;
            return out;
        }
    }

    /** Adds 'x' to the values of array 'input' and returns the values
     *  in the 'output' array if supplied, or a new array if output=null.
     *  When the 'input' and 'output' are identical, modifies the 'input' array. */
     //UNUSED
/*    public static double[] addToArray(double[] input, double x, double[] output) {
        if (output == null)
            output = new double[input.length];
        for (int i=0; i<Math.min(input.length, output.length); i++)
            output[i] = input[i] + x;
        return output;
    } */

    /** Returns the index of the heading String that is the energy,
     *  or -1 of none found */
    public static int getEnergyColumn(String[] headings) {
        int energyCol = -1;
        for (int i=0; i<headings.length; i++) {
            String head = headings[i].trim().toLowerCase();
            if (head.startsWith("energy") || head.equals("ev") || head.equals("e") || head.matches("e[\\s(\\[].*")) {
                energyCol = i;
                break;
            }
        }
        return energyCol;
    }

    /** Returns the absolute value of the (energy) step for an (energy) array,
     *  assuming the values are evenly spaced */
    public static double getEnergyStep(double[] a) {
        double step = (a[a.length-1] - a[0])/(a.length - 1);
        return Math.abs(step);
    }

    /** Returns the number of digits for writing the energy in the file */
    public static int getEnergyDigits(double[] energies) {
        int energyDigits = 0;
        double energyStep = getEnergyStep(energies);
        if (LeedUtils.isInteger(energyStep))
            return 0;
        if (LeedUtils.isInteger(energyStep*10))
            return 1;
        return (int)Math.max(0, 2.5 - Math.log10(energyStep));
    }


    /** Returns magnitude of the smallest step between tow successive values in an array */
    public static double getSmallestStep(double[] a) {
        double minStep = 0;
        for (int i=1; i<a.length; i++)
            if (Math.abs(a[i]-a[i-1]) > minStep)
                minStep = Math.abs(a[i]-a[i-1]);
        return minStep;
    }

    /** For an array with linearly function, returns first value and
     *  increment (step height) as a two-element array.
     *  Returns NaN, NaN if it is not a linear function
     *  (nonlinearity > 0.5 average step heights) */
    public static double[] getFirstAndIncrement(double[] a) {
        double step = (a[a.length-1] - a[0])/(a.length - 1);
        double yLin=a[0]+step;
        for (int i=1; i<a.length; i++, yLin+=step)
            if (Math.abs(a[i] - yLin) > 0.5*Math.abs(step))
                return new double[] {Double.NaN, Double.NaN};
        return new double[] {a[0], step};
    }

    /** Creates a ResultsTable for energies and beam data, e.g. I(V) curves.
     *  The heading of the energy column is written as "E"; the headings of
     *  the other columns according to 'headings'.
     *  If 'zeroReplacement' is nonzero, values with an absolute value less
     *  than this number are replaced by this number (with the correct sign,
     *  if nonzero). This should avoid reading values of zero as undefined
     *  by some LEED programs. */
    public static ResultsTable getResultsTable(double[] energies, String[] headings, double[][] data, double zeroReplacement) {
        if (headings.length != data.length)
            throw new RuntimeException("Number of headings different from number of data arrays: "+headings.length+" != "+data.length);
        ResultsTable rt = new ResultsTable(energies.length);

        for (int i=0; i<energies.length; i++)
            rt.setValue("E", i, energies[i]);

        for (int col=0; col<data.length; col++) {
            String head = headings[col];
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
    }

    /** With toArray=true, returns a moving-average filtered input array 'a'.
     *  Otherwise, returns a four-element array with [0] the minimum and
     *  [1] the maximum of the moving-average filtered array 'a', and
     *  [2,3] the repective indices. 
     *  Filtering is done by averaging over (2*radius + 1) values.
     *  Regions within 'radius' from the array boundaries or 'NaN' are filled
     *  by the nearest average of at least (2*radius + 1) values.
     *  In non-NaN regions shorter than (2*radius + 1), the output is NaN. */
    public static double[] movingAverage(double[] a, int radius, boolean toArray) {
        double[] out = new double[toArray ? a.length : 4];
        Arrays.fill(out, Double.NaN);
        double min = Double.MAX_VALUE, max = -Double.MAX_VALUE;
        int iOfMin = -1, iOfMax = -1;
        int[] rangeLimits = getRangeLimits(a, 0, a.length);
        for (int range=0; range<rangeLimits.length/2; range++) {
            int rStart = rangeLimits[2*range];
            int rEnd = rangeLimits[2*range+1];
            if (rStart>=0 && rEnd-rStart > 2*radius) {
                double rSum = 0;
                for (int i=rStart; i<rStart+2*radius; i++)
                    rSum += a[i];
                for (int i=rStart+2*radius; i<rEnd; i++) {
                    rSum += a[i];
                    if (toArray)
                        out[i-radius] = rSum*(1.0/(2*radius+1));
                    else {
                        if (min > rSum) { min = rSum; iOfMin = i-radius; }
                        if (max < rSum) { max = rSum; iOfMax = i-radius; }
                    }
                    rSum -= a[i-2*radius];
                }
                if (toArray) {
                    for (int i=rStart; i<rStart+radius; i++)    //near-end values: take the same as first/last smoothed value
                        out[i] = out[rStart+radius];
                    for (int i=rEnd-radius; i<rEnd; i++)
                        out[i] = out[rEnd-radius-1];
                } else {
                    out[0] = min/(2*radius+1);
                    out[1] = max/(2*radius+1);
                    out[2] = iOfMin;
                    out[3] = iOfMax;
                }
            }
        }
        return out;
    }

    /* DEBUG  static{
    double[] data = new double[] {0, 1, -2, 3, -4, 5, -6, 7, -8, 9, 10, 6, 3, 1, 0};
    double[] out = minimumFilter(data,2);
    IJ.log(Arrays.toString(out));}/**/

    /** A moving minimum filter, reports the local minimum within an environment of
     *  2*radius+1 length. Output is NaN for indices where the 2*radius+1 environemnt
     *  is not fully defined (contains NaN values or reaches beyond the boudaries of
     *  the input array) */
    public static double[] minimumFilter(double[] a, int radius) {
        double[] out = new double[a.length];
        Arrays.fill(out, Double.NaN);
        int[] rangeLimits = getRangeLimits(a, 0, a.length);
        for (int range=0; range<rangeLimits.length/2; range++) {
            int rStart = rangeLimits[2*range];
            int rEnd = rangeLimits[2*range+1];
            if (rStart>=0 && rEnd-rStart > 2*radius) {
                double lastMin = Double.NaN;
                for (int i=rStart+radius; i<rEnd-radius; i++) {
                    double min = a[i+radius];
                    if (!(min <= lastMin)) {                    //new value is not smaller?
                        if (!Double.isNaN(lastMin) && a[i-radius-1] > lastMin) {
                            min = lastMin;                      //minimum unchanged
                        } else {                                //otherwise determine the minimum the slow way
                            for (int j=i-radius; j<i+radius; j++)
                                if (a[j] < min)
                                    min = a[j];
                        }
                    }
                    out[i] = min;
                    lastMin = min;
                }
            }
        }
        return out;
    }

    /** Inverts a symmetric n*n matrix; the matrix is replaced by its inverse.
     *  Returns false on error (if not a singular matrix).
     *  This is essentially Gauss-Jordan elimination with the bottom left of the
     *  matrix used as storage for the output. It needs a reasonably well-conditioned
     *  matrix (more complicated algorithms may be more robust). */
    public static boolean invertSymmetricMatrix(double[][] matrix) {
        int n = matrix.length;
        if (matrix[0].length != n)
            throw new RuntimeException("invertSymmetricMatrix - not a square matrix");
        //DEBUG for (int r=0; r<n; r++) {String s=""; for (int c=0; c<n; c++) s+=IJ.d2s(matrix[r][c], -3)+", "; IJ.log(s);}

        boolean[] processed = new boolean[n];   //which rows we have inverted already
        double[] tmp1 = new double[n];
        double[] tmp2 = new double[n];
        double pivot0 = 1e-200;                 //this will become the absolute value of the first (largest) pivot
        for (int i=0; i<n; i++) {               //loop: process all matrix lines&columns
            double max = pivot0 * 1e-14;        //pivot must be bigger than this
            int iPivot = -1;
            for (int j=0; j<n; j++)             //find pivot
                if (!processed[j]) {
                    double abs = Math.abs(matrix[j][j]);
                    if (abs > max)  {
                        max = abs;
                        iPivot = j;
                    }
            }
            if (iPivot < 0) {                   //no pivot large enough, looks singular
            if (IJ.debugMode) IJ.log("singular "+n+"x"+n+" matrix after "+i+" steps");
                return false;
            }
            if (i==0)
                pivot0 = max;
            tmp2[iPivot] = 1.;                  //invert row 'iPivot' now:
            double norm = 1./matrix[iPivot][iPivot]; //multiplication with 1/pivot is faster than division by pivot
            tmp1[iPivot] = norm;
            matrix[iPivot][iPivot] = 0.;
            for (int j=0; j<iPivot; j++) {      //read the full row from the top-right, 1st part as column, 2nd part as row
                tmp2[j] = matrix[j][iPivot];
                tmp1[j] = matrix[j][iPivot]*(processed[j] ? norm:-norm);
                matrix[j][iPivot] = 0.;
            }
            for (int j=iPivot+1; j<n; j++) {
                tmp2[j] = processed[j] ? -matrix[iPivot][j] : matrix[iPivot][j];
                tmp1[j] = -matrix[iPivot][j]*norm;
                matrix[iPivot][j] = 0.;         //this row has been done
            }
            for (int j=0; j<n; j++)
                for (int k=j; k<n; k++)
                    matrix[j][k] += tmp2[j]*tmp1[k];
            processed[iPivot] = true;
        } //for i
        for (int i=1; i<n; i++)                 //fill bottom left to make it symmetric again
            for (int j=0; j<i; j++)
                matrix[i][j] = matrix[j][i];
        return true;
    }

    private static long lastBeepTimeMillis;
    /** Shows an error message in the 'Log' window, puts the Log window to the foreground,
     *  and beeps */
    public static void logError(String text) {
        IJ.log(text);
        IJ.selectWindow("Log");
        long time = System.currentTimeMillis();
        if (time - lastBeepTimeMillis > 2000)
            IJ.beep();      //don't beep more often than every 2 sec
        lastBeepTimeMillis = time;
    }

    /** Sets the items of a choice according to the specified list. Does nothing if no change.
     *  Returns whether there was a change */
    public static boolean setChoiceItems(Choice choice, String[] items) {
        if (!choiceItemsEqual(choice, items)) {
            choice.removeAll();
            for (String s : items)
                choice.add(s);
            choice.setSize(choice.getPreferredSize());
            return true;
        }
        return false;
    }

    /** Tests whether the choice items are equal to the list supplied */
    static boolean choiceItemsEqual(Choice choice, String[] items) {
        if (choice.getItemCount() != items.length)
            return false;
        for (int i=0; i<items.length; i++)
            if (!choice.getItem(i).equals(items[i]))
                return false;
        return true;
    }

    /** Displays a dialog with 'yes' and 'cancel' and returns true unless 'cancel' or ESC was pressed */
    public static boolean yesCancelDialog(String title, Frame parent, String message, String yesLabel) {
        GenericDialog gd = new GenericDialog(title, parent == null ? IJ.getInstance() : parent);
        gd.addMessage(message);
        gd.setOKLabel(yesLabel);
        gd.showDialog();
        return !gd.wasCanceled();
    }

    /** Duplicates a stack slice; like Duplicator.run but without IJ.showStatus() */
    public static ImagePlus duplicateSlice(ImagePlus stackImp, int slice) {
        ImageStack stack = stackImp.getStack();
        ImageProcessor ip = (stack == null || stack.size() == 0) ? null : stack.getProcessor(slice);
        if (ip == null) return null;
        ImagePlus imp = stackImp.createImagePlus();
        imp.setProcessor(ip);
        imp.setTitle(stackImp.getTitle());
        Overlay overlay = stackImp.getOverlay();
		if (overlay!=null && !stackImp.getHideOverlay()) {
			Overlay overlay2 = overlay.duplicate();
			overlay2.crop(slice, slice);
			imp.setOverlay(overlay2);
        }
        return imp;
    }

    /** Returns a mask from an ImagePlus with an 8-bit binary image as an ImageProcessor.
     *  The mask pixels are 0 for background (white) and 255 for foreground (black) pixels,
     *  irrespective of whether the image has an inverted LUT or not.
     *  There is no checking whether it is a binary image. Retunrs null if makImp is null. */
    public static ByteProcessor getMaskIp(ImagePlus maskImp) {
        ImageProcessor maskIp = maskImp==null ? null : maskImp.getProcessor();
        if (!(maskIp instanceof ByteProcessor)) return null;
        if (maskIp.isInvertedLut()) {
            return (ByteProcessor)maskIp;
        } else {                    //not an inverted LUT? create a copy with inverted pixel values
            ByteProcessor maskIp2 = (ByteProcessor)maskIp.createProcessor(maskIp.getWidth(), maskIp.getHeight());
            byte[] maskPixels  = (byte[])maskIp.getPixels();
            byte[] maskPixels2 = (byte[])maskIp2.getPixels();
            for (int i=0; i<maskPixels.length; i++)
                maskPixels2[i] = (byte)(maskPixels[i]^0xff);
            return maskIp2;
        }
    }

    /** Returns a mask from an ImagePlus with an 8-bit binary image as a pixels array.
     *  The mask pixels are 0 for background (white) and 255 for foreground (black) pixels,
     *  irrespective of whether the image has an inverted LUT or not.
     *  There is no checking whether it is a binary image. Returns null if maskImp is null. */
    public static byte[] getMaskPixels(ImagePlus maskImp) {
        ByteProcessor maskIp = getMaskIp(maskImp);
        return maskIp == null ? null : (byte[])maskIp.getPixels();
    }

    /** Returns the hash code of the pixels array of an image (must not be null).
     *  If a stack, the current slice is sude. */
    public static int getHashCode(ImagePlus imp) {
        ImageProcessor ip = imp.getProcessor();
        Object pixels = ip.getPixels();
        if (pixels instanceof byte[])
            return Arrays.hashCode((byte[])pixels);
        else if (pixels instanceof short[])
            return Arrays.hashCode((short[])pixels);
        else if (pixels instanceof float[])
            return Arrays.hashCode((float[])pixels);
        //else if (pixels instanceof int[])       //RGB images not used in ViPErLEED
        //    return Arrays.hashCode((int[])pixels);
        else
            throw new IllegalArgumentException("Unsupported image data type: "+getSimpleClassName(pixels));
    }

    /** Human-readable class name of an object. Returns "null" if the 'obj' argument is null */
    static String getSimpleClassName(Object obj) {
        if (obj == null) return "null";
        String s = obj.getClass().getSimpleName();
        if (s.length() == 0) s = "<anonymous>";
        return s;
    }

    /** Returns whether a file with a given path exists as a normal file (not a directory) and is readable */
    public static boolean fileOk(String path) {
        if (path == null || path.length() == 0) return false;
        File file = new File(path);
        return file.isFile() && file.canRead();
    }

    /** Returns whether a file starts with a given magic code (ascii only) */
    public static boolean hasMagic(File file, String magic) {
        if (!file.canRead()) return false;
        FileReader reader = null;
        try {
            reader = new FileReader(file);
            int length = magic.length();
            char[] cbuf = new char[length];
            reader.read(cbuf);
            reader.close();
            return Arrays.equals(cbuf, magic.toCharArray());
        } catch (Exception e) {
            try {
                reader.close();
            } catch (Exception e2) {}
            return false;
        }
    }
    /** Returns whether a mouse event should trigger a popup menu. Returns false in case of a null argument*/
    public static boolean isPopupTrigger(MouseEvent e) {
        if (e == null) return false;
        int flags = e.getModifiers();
        return e.isPopupTrigger()||(!IJ.isMacintosh()&&(flags&java.awt.Event.META_MASK)!=0);
    }
}
