import ij.*;
import ij.IJ;
import ij.gui.*;
import ij.util.Tools;
import ij.io.OpenDialog;
import java.awt.*;
import java.io.File;
import java.util.Arrays;
import java.util.ArrayList;


/** This class contains a list of spots, each with designation, h & k indices,
 *  and a spot group (symmetry-equivalent spots are in the same group).
 *  These are read from a Spot Pattern Definition File.
 *
 *  Format of the Spot Pattern Definition File:
 *
 *  The spot pattern definition file is based on the csv (comma-separated values) format.
 *  All leading and trailing spaces and empty lines are ignored
 *  Comment lines can appear anywhere and start with any of these
 *  characters: #%!*
 *  The first non-comment line is the header (ignored when reading).
 *  Then for each spot, there is a line with comma-delimited items:
 *  indices ( h k ) as fractions (opionally in brackets) separated by whitespace, underscore, semicolon, or vertical bar,
 *  h(float), and k(float),
 *  gx, and gy in cartesian coordinates,
 *  beam group index: the same number for symmetry-equivalent beams. Spots of the same group must be in adjacent lines.
 *  further items may be present (currently igonred)
 *
 *  Example:
(  h     k  ),    h   ,     k  ,     gx  ,    gy   , group
(  1/2   0  ), 0.50000, 0.00000,  1.50000,  0.00000,     0
( -1/2   1/2),-0.50000, 0.50000, -0.75000,  1.29903,     0
(  0    -1/2), 0.00000,-0.50000, -0.75000, -1.29903,     0
(  1     0  ), 1.00000, 0.00000,  3.00000,  0.00000,     1
( -1    -1  ),-1.00000, 1.00000, -1.50000,  2.59808,     1
(  0    -1  ), 0.00000,-1.00000, -1.50000, -2.59808,     1
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

public class LeedSpotPattern {
    public static final int UNKNOWN_GROUP = -9999;    //group number if spot group is unknown (not too negative since this may become first group; then loops would run very long)
    static final int N_NEAREST = 20;    //maximum number of nearest spots to list for each spot
    int size;
    String[]  spotNames;      // spot name, comma-separated h,k, no brackets, e.g. "1/2,3"
    String[]  spotNamesForOvly; // spot name, but with ',' replaced by '\n' for superstructure (for overlay)
    double[]  spotHs;         // h index of spot
    double[]  spotKs;         // k index of spot
    double[]  spotKx;         // cartesian reciprocal-space coordinate kx
    double[]  spotKy;         // cartesian reciprocal-space coordinate ky
    int[]     spotGroups;     // group of symmetry-equivalent spots that the spot belongs to
    int[][]   nearestSpots;   // for each spot, the list of indices of nearest spots, sorted by distance (ascending)
    boolean[] isSuperstr;     // for each spot, whether it is a superstructure spot
    boolean   missingGroup;   // group information missing
    double    minDistanceSqr; // distance^2 of the two spots closest to each other in reciprocal space
    String  path;
    String  errorMessage;
    boolean hasSuperstructure;
    boolean hasEquivalentBeams;

    /** Creates a LeedSpotPattern by reading from a file.
     *  After the constructor was called, <code>getErrorMessage()</code> returns null if successful,
     *  or an error String on failure */
    public LeedSpotPattern(String path) {
        if (path == null || path.length() == 0) {
            errorMessage = "No spot pattern file given";
            return;
        }
        String str = IJ.openAsString(path);
        if (str == null || str.startsWith("Error:")) {      // failure
            errorMessage = "Cannot open spot pattern file";
            return;
        }
        if (str.indexOf("SUPERLATTICE")>=0 && str.indexOf("bulkGroup")>=0) {
            errorMessage = "Wrong spot pattern file format (looks like input for GUI Pattern Simulator).";
            return;
        }
        this.path = path;
        String[] lines = Tools.split(str, "\n");
        int n = -1;             // minus one for the header
        for (int l=0; l<lines.length; l++) {
            lines[l] = lines[l].trim();
            if (LeedUtils.isCommentLine(lines[l]))
                lines[l] = null;
            else
                n++;
        }
        makeArrays(n, true);
        spotNames = new String[n];
        spotNamesForOvly = new String[n];
        spotHs = new double[n];
        spotKs = new double[n];
        spotKx = new double[n];
        spotKy = new double[n];
        spotGroups = new int[n];
        isSuperstr = new boolean[n];
        boolean waitForHeader = true;  //the first non-comment line is header
        int i = 0;
        int lastGroup = Integer.MIN_VALUE;
        for (int l=0; l<lines.length; l++) {
            String line = lines[l];
            if (line != null) {
                if (waitForHeader) {   //skip the header line
                    waitForHeader = false;
                    continue;
                }
                String[] items = line.split(",");
                if (items.length < 6) {
                    errorMessage = "Spot pattern file: Line "+(l+1)+" too short:\n"+line;
                    return;
                }
                String withoutBrackets = items[0].trim().replaceFirst("\\(\\s*", "").replaceFirst("\\s*\\)", "");
                spotNames[i] = withoutBrackets.trim().replaceFirst("[\\s;_|]+",",");
                spotHs[i] = Tools.parseDouble(items[1]);
                spotKs[i] = Tools.parseDouble(items[2]);
                spotKx[i] = Tools.parseDouble(items[3]);
                spotKy[i] = Tools.parseDouble(items[4]);
                double spotGroup = (Tools.parseDouble(items[5]));
                if (Double.isNaN(spotHs[i]) || Double.isNaN(spotKs[i]) || !((int)spotGroup == spotGroup)) {
                    errorMessage = "Spot pattern file: Invalid number in line"+(l+1)+":\n"+line;
                    return;
                }
                spotGroups[i] = (int)spotGroup;
                lastGroup = spotGroups[i];
                isSuperstr[i] = !LeedUtils.isInteger(spotHs[i]) || !LeedUtils.isInteger(spotKs[i]);
                spotNamesForOvly[i] = isSuperstr[i] ? spotNames[i].replace(',','\n') : spotNames[i];
                if (isSuperstr[i])
                    hasSuperstructure = true;
                i++;
            }
        }
        size = n;
        if (size < 2 && errorMessage == null)
            errorMessage = "Spot pattern file: Only "+size+" spot(s) in file, this is not enough";
        makeNearestTable();
        checkForEquivalentBeams();
    }

    /** Creates a LeedSpotPattern from an array of table headings of an I(V) file in the form:
     *    E, 1/2|0 [1], 0|1/2 [1], ...
     *  where the numbers in brackets are the spot groups. An energy column may be present or not.
     *  There must be no duplicate spots.
     *  The spot pattern created this way has no Cartesian coordinates Kx, Ky!
     *  Calls that rely on Cartesian coordinates such as getNearestSpotIndices and getMinDistance
     *  will cause an exception! See isFullSpotPattern().
     *  The spot pattern has only h and k indices and groups.
     *  If the headings contain no groups in square brackets and 'groupsRequired' is true,
     *  an empty LeedSpotPattern will be created, i.e. size() will return 0, and other calls
     *  will result in an exception. The same is true in case of other errors.
     *  After the constructor was called, <code>getErrorMessage()</code> returns null if successful,
     *  or an error String on failure */
    public LeedSpotPattern(String[] headings, boolean groupsRequired) {
        int energyColumn = LeedUtils.getEnergyColumn(headings);
        int n = energyColumn >= 0 ? headings.length-1 : headings.length;
        makeArrays(n, false);
        
        for (int iHead=0, iSpot=0; iHead<headings.length; iHead++) {
            if (iHead != energyColumn) {
                int braPos = headings[iHead].indexOf('[');
                int ketPos = headings[iHead].indexOf(']');
                if (braPos < 0 || ketPos < braPos) {
                    if (groupsRequired) {
                        errorMessage = "Missing group in '"+headings[iHead].trim()+"'";
                        return;
                    }
                }
                double group = braPos>=0 && ketPos>braPos ?
                        Tools.parseDouble(headings[iHead].substring(braPos+1, ketPos)) : UNKNOWN_GROUP;
                if (!(group ==(int)group))  //NaN or non-integer
                    group = UNKNOWN_GROUP;
                if (groupsRequired && group == UNKNOWN_GROUP) {
                    errorMessage = "Invalid group in '"+headings[iHead].trim()+"'";
                    return;
                }

                int hAndKEndPos = braPos>=0 ? braPos : headings[iHead].length() ;
                String[] hAndK = headings[iHead].substring(0, hAndKEndPos).replaceAll("[()]", "").trim().split("[,;_|\\s]+");
                if (hAndK.length != 2) {
                    errorMessage = "Not a pair of h, k indices in '"+headings[iHead].trim()+"'";
                    return;
                }
                double h = LeedUtils.toNumber(hAndK[0]);
                double k = LeedUtils.toNumber(hAndK[1]);
                if (Double.isNaN(h) || Double.isNaN(k)) {
                    errorMessage = "Non-numeric or missing h, k index in '"+headings[iHead].trim()+"'";
                    return;
                }
                int alreadyPresentIndex = getSpotIndex(h, k, iSpot);
                if (alreadyPresentIndex >= 0) {        //a duplicate spot
                    errorMessage = "Duplicate spot (inconsistent name/group?) '"+headings[iHead].trim()+"'";
                    return;
                }
                spotNames[iSpot] = hAndK[0].trim() + ',' + hAndK[1].trim();
                spotNamesForOvly[iSpot] = isSuperstr[iSpot] ? spotNames[iSpot].replace(',','\n') : spotNames[iSpot];
                spotHs[iSpot] = h;
                spotKs[iSpot] = k;
                spotGroups[iSpot] = (int)group;
                isSuperstr[iSpot] = !LeedUtils.isInteger(spotHs[iSpot]) || !LeedUtils.isInteger(spotKs[iSpot]);
                if (isSuperstr[iSpot])
                    hasSuperstructure = true;
                if (group == UNKNOWN_GROUP)
                    missingGroup = true;
                iSpot++;
            }
        }
        size = n;
        checkForEquivalentBeams();
    }

    /** Creates a LeedSpotPattern by merging two LeedSpotPatterns.
     *  When groupsRequired is true, both LeedSpotPatterns must have spot groups.
     *  If this is not the case or in case of discordant group information,
     *  getErrorMessage() will return a non-null error text.
     *  The output is sorted in sequence of increasing group number. */
    public LeedSpotPattern(LeedSpotPattern spotPattern0, LeedSpotPattern spotPattern1, boolean groupsRequired) {
        int n = spotPattern0.size() + spotPattern1.size();
        boolean withKxKy = spotPattern0.spotKx != null && spotPattern1.spotKx != null;
        makeArrays(n, withKxKy);
        int[] loHiG0 = spotPattern0.getLowestAndHighestGroup();
        int[] loHiG1 = spotPattern1.getLowestAndHighestGroup();
        if (groupsRequired && loHiG1 == null)
            errorMessage = "Information on equivalent spots missing in second file";
        if (groupsRequired && loHiG0 == null)
            errorMessage = "Information on equivalent spots missing in first file";
        if (errorMessage != null) return;
        if (loHiG0 != null && loHiG1 != null) {     // merge with groups, sorted by group
            int lowestGroup  = Math.min(loHiG0[0], loHiG1[0]);
            int highestGroup = Math.max(loHiG0[1], loHiG1[1]);
            size = 0;
            for (int g=lowestGroup; g<=highestGroup; g++) {
                int[] spots0 = spotPattern0.getAllSpotsForGroup(g);
                int[] spots1 = spotPattern1.getAllSpotsForGroup(g);
                for (int i0 : spots0) {
                    int i1 = spotPattern1.getSpotIndex(spotPattern0.spotHs[i0], spotPattern0.spotKs[i0]);
                    if (i1 >= 0) {
                        if (spotPattern1.spotGroups[i1] != g) {
                            errorMessage = "Groups of symmetry-equivalent beams disagree:\n"+
                                    spotPattern0.getNameWithGroup(i0, true) +" vs. "+spotPattern1.getNameWithGroup(i1, true);
                            size = 0;
                            return;
                        }
                        int index1 = LeedUtils.arrayIndexOf(spots1, i1);
                        if (index1 < 0) throw new RuntimeException("Internal error merging spot patterns");
                        spots1[index1] = -1;        //don't handle this spot twice
                    }
                    addSpot(spotPattern0, i0);
                }
                for (int i1 : spots1) {
                    if (i1 < 0) continue;
                    addSpot(spotPattern1, i1);
                }
            }
        } else {                                    //merge without group information
            ArrayList<String> spotNames = new ArrayList<String>(n);
            for (int i0=0, i1=0; i0<spotPattern0.size() || i1<spotPattern0.size(); ) {
                String name0 = i0<spotPattern0.size() ? spotPattern0.getName(i0) : null;
                String name1 = i1<spotPattern1.size() ? spotPattern1.getName(i1) : null;
                if (name0 != null && name0.equals(name1)) { //common spot, add it
                    addSpot(spotPattern0, i0);
                    i0++; i1++;
                } else if (name0 == null ||         //#0 is exhausted or comes later, add #1
                    spotPattern1.getIndex(name0) > i1) {
                    if (getIndex(name1)<0)          //add spot if not present already
                        addSpot(spotPattern1, i1);
                    i1++;
                } else if (name1 == null ||         //#1 is exhausted or comes later, add #0
                    spotPattern0.getIndex(name1) > i0) {
                    if (getIndex(name0)<0)          //add spot if not present already
                        addSpot(spotPattern0, i0);
                    i0++;
                } else {                            //otherwise add both
                    if (getIndex(name0)<0)
                        addSpot(spotPattern0, i0);
                    if (getIndex(name1)<0)
                        addSpot(spotPattern1, i1);
                    i0++; i1++;
                }
            }
        }
        if (withKxKy) makeNearestTable();
        checkForEquivalentBeams();
    }

    /** Creates the arrays with a size of 'n'.
     *  The Kx, Ky arrays (cartesian coordinates) are created if withKxKy=true. */
    private void makeArrays(int n, boolean withKxKy) {
        spotNames = new String[n];
        spotNamesForOvly = new String[n];
        spotHs = new double[n];
        spotKs = new double[n];
        spotGroups = new int[n];
        isSuperstr = new boolean[n];
        if (withKxKy) {
            spotKx = new double[n];
            spotKy = new double[n];
        }
    }
    /** Adds a new spot from a different spot pattern, where i is the index in that spot pattern.
     *  The arrays must have sufficient size. */
    private void addSpot(LeedSpotPattern spotPattern, int i) {
        spotNames[size] = spotPattern.spotNames[i];
        spotNamesForOvly[size] = spotPattern.spotNamesForOvly[i];
        spotHs[size] = spotPattern.spotHs[i];
        spotKs[size] = spotPattern.spotKs[i];
        spotGroups[size] = spotPattern.spotGroups[i];
        isSuperstr[size] = spotPattern.isSuperstr[i];
        if (spotKx != null) {
            spotKx[size] = spotPattern.spotKx[i];
            spotKy[size] = spotPattern.spotKy[i];
        }
        if (isSuperstr[size]) hasSuperstructure = true;
        if (spotGroups[size] == UNKNOWN_GROUP) missingGroup = true;
        size++;
    }

    /** Sets the 'hasEquivalentBeams' variable, to true if there is
     *  at least one group with at least two symmetry-equivalent spots.
     *  Sets 'hasEquivalentBeams' to false if group information is not available. */
    private void checkForEquivalentBeams() {
        if (spotGroups == null) {
            hasEquivalentBeams = false;
            return;
        }
        int[] loHiGr = getLowestAndHighestGroup();
        if (loHiGr == null) {
            hasEquivalentBeams = false;
            return;
        }
        for (int g=loHiGr[0]; g<=loHiGr[1]; g++) {
            int i0 = getFirstSpotForGroup(g);
            for (int i=i0+1; i<size; i++) {
                if (spotGroups[i] == g) {
                    hasEquivalentBeams = true;
                    return;
                }
            }
        }
    }

    /** Returns the number of spots; 0 if not successfully read */
    public int size() {
        return size;
    }

    /** Returns whether this is a complete spotPattern read from a file
     *  (with Cartesian values in the reciprocal space).
     *  Returns false for spot patterns read 'on the fly' from the headings of a data file */
    public boolean isFullSpotPattern() {
        return spotKx != null;
    }

    /** Returns the error message during reading or null if ok */
    public String getErrorMessage() {
        return errorMessage;
    }

    /** Returns the file path where this was read, or null */
    public String getPath() {
        return path;
    }

    /** Returns the filename where this was read, or null */
    public String getFileName() {
        return (new File(path)).getName();
    }

    /** Returns the name (without brackets, e.g. "1/2,2/5") of the i-th spot */
    public String getName(int i) {
        return spotNames[i];
    }

    /** Returns the name for the overlay (without brackets, and with
     *  comma replaced by line feed for superstructure spots), for the i-th spot */
    public String getNameForOvly(int i) {
        return spotNamesForOvly[i];
    }

    /** Returns the name with brackets and group number, e.g. "(1/2,2/5) [5]" of the i-th spot.
     *  The square brackets with group number are omitted for unknown groups.
     *  If 'replaceComma' is true, replaces the comma by a vertical bar. */
    public String getNameWithGroup(int i, boolean replaceComma) {
        String str = '('+spotNames[i];
        if (replaceComma) str = str.replace(',', '|');
        str += ")";
        if (spotGroups[i] != UNKNOWN_GROUP) str+= " [" + spotGroups[i] + ']';
        return str;
    }

    /* unused * Returns the Miller index 'h' of the i-th spot */
    /*public double getH(int i) {
        return spotHs[i];
    }*/

    /* unused * Returns the Miller index 'k' of the i-th spot */
    /*public double getK(int i) {
        return spotKs[i];
    }*/

    /** Returns the reciprocal-space x-coordinate of the i-th spot */
    public double getKx(int i) {
        return spotKx[i];
    }

    /** Returns the reciprocal-space x-coordinate of the i-th spot */
    public double getKy(int i) {
        return spotKy[i];
    }

    /** Returns the spot group of the i-th spot */
    public int getGroup(int i) {
        return spotGroups[i];
    }

    /** Returns the first spot index for a given group, or -1 if none */
    public int getFirstSpotForGroup(int group) {
        return LeedUtils.arrayIndexOf(spotGroups, group);
    }

    /** Returns the indices of all spots for a given group, or an empty array if none */
    public int[] getAllSpotsForGroup(int group) {
        LeedIntegerArray array = new LeedIntegerArray(12);
        for (int i=0; i<size; i++)
            if (spotGroups[i] == group)
                array.add(i);
        return array.toArray();
    }

    /** Returns the group numbers of the lowest and highest group as a two-element array.
     *  Returns null if there are no group names. */
    public int[] getLowestAndHighestGroup() {
        int lowestGroup = Integer.MAX_VALUE;
        int highestGroup = Integer.MIN_VALUE;
        for (int i=0; i<size; i++) {
            int group = spotGroups[i];
            if (group != UNKNOWN_GROUP && group <  lowestGroup)  lowestGroup = group;
            if (group != UNKNOWN_GROUP && group > highestGroup) highestGroup = group;
        }
        if (lowestGroup > highestGroup)
            return null;
        else
            return new int[] {lowestGroup, highestGroup};
    }

    /** Returns whether this pattern is a superstructure */
    public boolean isSuperstructure() {
        return hasSuperstructure;
    }

    /** Returns whether this pattern has symmetry-equivalent beams,
     *  i.e. beams belonging to the same group. */
    public boolean hasEquivalentBeams() {
        return hasEquivalentBeams;
    }

    /** Returns whether this pattern has information on goups of symmetry-equivalent beams */
    public boolean hasGroups() {
        return !missingGroup;
    }

    /** Returns whether the given spot is a superstructure spot */
    public boolean isSuperstructure(int i) {
        return isSuperstr[i];
    }

    /** Returns the index of the spot with h, k matching within 0.01 or -1 if none */
    public int getSpotIndex(double h, double k) {
        return getSpotIndex(h, k, size);
    }

    /** Returns the index of the spot with h, k matching within 0.01 or -1 if none, for indices < 'end' */
    private int getSpotIndex(double h, double k, int end) {
        for (int i=0; i<end; i++)
            if (Math.abs(h - spotHs[i]) < 0.01 && Math.abs(k - spotKs[i]) < 0.01) return i;
        return -1;
    }

    /** For a given spot index, returns an array with up to N_NEAREST indices of
     *  the nearest neighbor spots as determined by their cartesian coordinates.
     *  The array is sorted by ascending distance.
     *  Do not modify this array! */
    public int[] getNearestSpotIndices(int i) {
        return nearestSpots[i];
    }

    /** Returns the smallest distance in reciprocal (cartesian) space between any two spots
     *  of this pattern (i.e., considering all possible pairs of two spots) */
    public double getMinDistance() {
        return Math.sqrt(minDistanceSqr);
    }

    /** Returns an array of all spot names */
    public String[] getAllSpotNames() {
        return Arrays.copyOf(spotNames, size);
    }

    /** Returns an array of all spot names with groups (where available),
     *  with '|' (not comma) between h and k */
    public String[] getAllSpotNamesWithGroup() {
        String[] out = new String[size];
        for (int i=0; i<size; i++)
            out[i] = getNameWithGroup(i, /*replaceComma=*/true);
        return out;
    }

    /** Returns the indices corresponding to the spot names given as a list.
     *  If there are indices that have not been found or have an invalid syntax,
     *  getErrorMessage() returns at text on it (otherwise it returns null after this call),
     *  and the corresponding indices are set to -1 */
    public int[] getIndices(String[] names) {
        int[] indices = new int[names.length];
        errorMessage = null;
        String error = "";
        int nErr = 0;
        final int maxNErr = 5;
        for (int i=0; i<names.length; i++) {
            int group = LeedUtils.getIntInBrackets(names[i]);
            String[] hAndK = names[i].replaceAll("[()]", "").replaceAll("\\[.*\\]","").trim().split("[,;_|\\s]+");
            double h = Double.NaN, k = Double.NaN;
            if (hAndK.length == 2) {
                h = LeedUtils.toNumber(hAndK[0]);
                k = LeedUtils.toNumber(hAndK[1]);
            }
            if (Double.isNaN(h) || Double.isNaN(k)) {
                if (nErr < maxNErr) {
                    if (nErr > 0) error += "; ";
                    error += "invalid spot: "+names[i];
                }
                nErr++;
                indices[i] = -1;
            } else {
                int index = getSpotIndex(h, k);
                indices[i] = index;
                if (index < 0) {
                    if (nErr < maxNErr) {
                        if (nErr > 0) error += "\n";
                        error += "spot not found: "+names[i];
                    }
                    nErr++;
                } else if (group > Integer.MIN_VALUE && group != getGroup(index)) {
                    if (nErr < maxNErr) {
                        if (nErr > 0) error += "\n";
                        error += getName(index)+"["+getGroup(index)+"]: group different from spotPattern: "+group;
                    }
                    nErr++;
                }
            }
        //DEBUG IJ.log("("+i+") "+names[i]+" -> "+(i>=0 ? ("("+indices[i]+") "+spotNames[indices[i]]) : "NONE")+" g="+spotGroups[indices[i]]);
        }
        if (nErr >= maxNErr)
            error += "\n...";
        if (nErr > 0) errorMessage = error;
        return indices;
    }

    /** Returns the index corresponding to the given spot name, or -1 if not found.
     *  The name my include round brackets around the Miller indices (h,k) or not.
     *  The h and k values may be separated by any of comma, semicolon, underscore,
     *  vertical bar or whitespace, as well a combinations thereof.
     *  The h and k values may be integer or decimal numbers or fractions like -1/2.
     *  The name may include a spot group in square brackets; this is ignored */
    public int getIndex(String name) {
        String[] hAndK = name.replaceAll("[()]", "").replaceAll("\\[.*\\]","").trim().split("[,;_|\\s]+");
        double h = Double.NaN, k = Double.NaN;
        if (hAndK.length == 2) {
            h = LeedUtils.toNumber(hAndK[0]);
            k = LeedUtils.toNumber(hAndK[1]);
        }
        if (Double.isNaN(h) || Double.isNaN(k)) return -1;
        else return getSpotIndex(h, k);
    }

    /** Creates the table of sorted nearest spots for each spot and finds the minimum distance between any two spots */
    void makeNearestTable() {
        if (size<2) return;
        int nN = Math.min(size-1, N_NEAREST);   //number of nearest spots to use
        nearestSpots = new int[size][nN];
        double[] nDistSqr = new double[nN];
        minDistanceSqr = Double.MAX_VALUE;
        for (int i=0; i<size; i++) {
            int nEntries = 0;
            for (int j=0; j<size; j++) if (j != i) {
                double distSqr = sqr(spotKx[i] - spotKx[j]) + sqr(spotKy[i] - spotKy[j]);
                if (distSqr < minDistanceSqr)
                    minDistanceSqr = distSqr;
                int n = nEntries-1;
                for (; n>=0; n--)
                    if (nDistSqr[n] < distSqr)
                        break;
                n++;    //the previous one was still further than the current one
                //IJ.log(i+" enter "+j+" as #"+n);
                if (n < nN) {
                    //IJ.log(i+" shift up "+n+" to "+Math.min(nEntries-1, nN-2)+" New entry at#"+n+" d="+(float)distSqr);
                    for (int k=Math.min(nEntries-1, nN-2); k>=n; k--) {   //shift up all higher entries
                        nDistSqr[k+1] = nDistSqr[k];
                        nearestSpots[i][k+1] = nearestSpots[i][k];
                    }
                    nDistSqr[n] = distSqr;
                    nearestSpots[i][n] = j;
                    if (nEntries < nN) nEntries++;
                }
            }
            //DEBUG String s = i+" nearest: ";for (int k=0; k<nN; k++) {s+=nearestSpots[i][k];s+=",";}IJ.log(s);
        }

    }

    /** Creates a string including the first 10 items for debug */
    public String toString() {
        StringBuilder sb = new StringBuilder(100);
        sb.append("LeedSpotPattern["+size+"]: ");
        for (int i=1; i<Math.min(10, size); i++)
            sb.append(spotNames[i]+"="+IJ.d2s(spotHs[i],2)+','+IJ.d2s(spotKs[i],2)+"; ");
        if (size>10) sb.append("...");
        return sb.toString();
    }

    /** Returns the reciprocal-space components of the vectors (1,0) and (0,1)
     *  as a four-element array kx10, ky10, kx01, ky01.
     *  The current version requires at least one spot (H, 0) and
     *  one spot (0, K) to be present.
     *  Returns null on failure. */
    /* currently unused 
    public static double[] getKxy10and01() {
        if (spotKx == null)     //can't do we have no cartesian reciprocal-space coordinates
            return null;
        double[] out = new double[4];
        boolean got10 = false;
        boolean got01 = false;
        for (int i=0; i<size; i++) {
            if (!got10 && spotHs[i] != 0 spotKs[i] == 0) {    //(H, 0) spot found
                out[2] = spotKx[i]/spotHs[i];
                out[3] = spotKy[i]/spotHs[i];
                got10 = true;
            }
            if (!got01 && spotHs[i] == 0 spotKs[i] != 0) {    //(0, K) spot found
                out[0] = spotKx[i]/spotKs[i];
                out[1] = spotKy[i]/spotKs[i];
                got01 = true;
            }
            if (got10 && got01)
                return out;
        }
        return null;
    } */

    /** Shows a file open dialog and returns the spot pattern selected, or null when cancelled.
     *  In case of an error reading the file, displays an error message and returns null. */
    public static LeedSpotPattern openWithDialog(String defaultPath) {
        OpenDialog od;
        String lastDirectory = OpenDialog.getLastDirectory();
        if (defaultPath != null && defaultPath.length() > 0) {
            File file = new File(defaultPath);
            od = new OpenDialog("Select Spot Pattern .csv File", file.getParent(), file.getName());
        } else
            od = new OpenDialog("Select Spot Pattern .csv File");
        od.setLastDirectory(lastDirectory);
        String patternFilePath = od.getPath();
        if (patternFilePath == null) return null;
        LeedSpotPattern spotPattern = new LeedSpotPattern(patternFilePath);
        String error = spotPattern.getErrorMessage();
        if (error != null) {
            IJ.error("Error reading LEED spot pattern", error);
            return null;
        } else
            return spotPattern;
    }

    double sqr(double x) {return x*x;}
}
