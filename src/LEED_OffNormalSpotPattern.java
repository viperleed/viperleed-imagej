import ij.*;
import ij.IJ;
import ij.gui.*;
import ij.io.OpenDialog;
import ij.io.SaveDialog;
import ij.process.*;
import ij.measure.Measurements;
import ij.plugin.PlugIn;
import ij.util.Tools;
import java.lang.ref.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
/**
 *  This class modifies a spot pattern to correctly handle the reduced symmetry
 *  for off-normal incidence.
 *  The groups of non-equivalent beams are modified accordingly.
 */

/** This ImageJ plugin is part of the ViPErLEED package for LEED I(V) analysis.
 *  Licensed under GNU General Public License v3.0 or later (GPL-3.0-or-later),
 *  https://www.gnu.org/licenses/gpl-3.0.html
 *  This means that the collection of ViPErLEED ImageJ plugins is free software:
 *  you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation,
 *  either version 3 of the License, or (at your option) any later version.
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  The collection of ViPErLEED ImageJ plugins is distributed
 *  in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  The authors may decide later to put part of the auxiliary code in this work into
 *  the public domain, to allow incorporation into ImageJ if desired.
 *  (ImageJ is in the public domain.)
 *
 *  The help text is licensed under the Creative Commons Attribution 4.0 (CC BY 4.0) license.
 *
 *  When using and/or modifying this program for scientific work, please cite
 *  the paper describing it:
 *  M. Schmid, F. Kraushofer, A. M. Imre, T. Kißlinger, L. Hammer, U. Diebold, and M. Riva,
 *  ViPErLEED package II: Spot tracking, extraction and processing of I(V) curves,
 *  Phys. Rev. Research 7, 013006 (2025).
 *  @author Michael Schmid, IAP/TU Wien, 2025
 */



public class LEED_OffNormalSpotPattern implements PlugIn {
    static final double MAX_TOLERANCE = 1;      // maximum tolerance value
    static final double MIN_TOLERANCE = 1e-4;   // minimum tolerance value (avoids non-equivalent beams due to rounding errors)
    static final int MAX_SPOTS_PER_GROUP = 48;  // should make sure we have no ArrayIndexOutOfBoundsException
    static String spotPatternPath;              // remembers last spot Pattern used as input
    static String referenceSpot = "1,0";        // remembers the spot for the azimuth orientation
    static double angleVsRefSpot = 0 ;          // remembers the angle of incidence pane vs. referenceSpot
    static double incidenceAngle = 1;           // remembers incidence angle (degrees)
    static double tolerance = 0.01;             // remembers tolerance, how far spots can be from normal to incidence plane

    private static final String HELP_STRING =
            "<html>"+
            "<h1>Create Spot Pattern for Off-Normal Incidence</h1>"+
            "<p>This command takes a spot pattern file for perpendicular incidence and derives from it "+
            "a new spot pattern file with reduced symmetry, to account for off-normal incidence.</p>"+
            "<h2>Dialog parameters</h2>"+
            "<ul>"+
            "<li><b>Spot pattern file (input)</b>: The spot pattern file (.csv file) assuming normal incidence. "+
            "Usually, this file is exported from the pattern simulator of the ViPErLEED GUI.</li>"+
            "<li><b>Incidence angle</b>: How far the incidence is off-normal. "+
            "Together with the tolerance, this determines which spots are considered sufficiently close to equivalent.</li>"+
            "<li><b>Incidence plane is at an angle of</b>: Enter the azimuthal angle of the incidence plane "+
            "with respect to the direction between (0,0) and the reference spot. "+
            "Angles are counted positive in counterclockwise direction.</li>"+
            "<li>...<b>with respect to this spot</b>: A spot present in the spot pattern file, "+
            "used as reference for the angle above.<br>"+
            "If you know two spots that are equivalent at the given incidence angle and are exactly opposite "+
            "[i.e., (h,k) and (-h,-k)], you can enter one of these spots here and enter an angle of 90&deg; in the field above.</li>"+
            "<li><b>Tolerance</b>: Spots that are linked by mirror symmetry "+
            "with respect to the plane of incidence are considered equivalent. "+
            "The tolerance describes how accurately this mirror symmetry must be fulfilled. "+
            "When averaging experimental AND calculated beams that are nearly symmetry equivalent, you can choose "+
            "a tolerance somewhat below the best R factor that you expect. Otherwise use, e.g., 0.01. "+
            "Avoid ridiculously low values to prevent the influence of rounding errors.</li>"+
            "</ul>"+
            "<h2><a name='ivEditorLicense'>License</a></h2>"+
            "<p>The code is licensed under <a href='http://www.gnu.org/licenses/gpl-3.0.html'>GNU General Public License v3.0</a> "+
            "or later (GPL-3.0-or-later).</p>"+
            "<p>&nbsp;&nbsp;&nbsp;&nbsp;The ViPErLEED ImageJ plugin collection is free software: you can redistribute it and/or modify it "+
            "under the terms of the GNU General Public License as published by the Free Software Foundation, "+
            "either version 3 of the License, or (at your option) any later version.<br>&nbsp;&nbsp;&nbsp;&nbsp;"+
            "The ViPErLEED ImageJ plugin collection is distributed in the hope that it will be useful, "+
            "but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. "+
            "See the GNU General Public License for more details.<br>&nbsp;&nbsp;&nbsp;&nbsp;"+
            "You should have received a copy of the GNU General Public License along with these plugins. "+
            "If not, see <a href='https://www.gnu.org/licenses/'>https://www.gnu.org/licenses/</a>.</p>"+
            "<p>The authors may decide later to put part of the auxiliary code in this work into the public domain, "+
            "to allow incorporation into ImageJ if desired (ImageJ is in the public domain).</p>"+
            "<p>The documentation, including the help texts, is licensed under the "+
            "<a href='http://creativecommons.org/licenses/by/4.0/'>Creative Commons Attribution 4.0</a> "+
            "(CC BY 4.0) license.</p>"+
            "<p>When using this program (in its original or modified form) for scientific work, "+
            "please cite the paper describing the program [<a href='#paper'>1</a>].</p>"+
            "<p>A copy of these license texts and the source code is included in the jar archive holding this plugin "+
            "(use an unzip utility to view its contents).</p>"+
            "<p><a name='paper'>[1]</a> M. Schmid, F. Kraushofer, A. M. Imre, T. Ki&szlig;linger, L. Hammer, U. Diebold, and M. Riva<br>"+
            "<i>ViPErLEED package II: Spot tracking, extraction and processing of <i>I</i>(<i>V</i>) curves</i><br>"+
            "<a href='https://doi.org/10.1103/PhysRevResearch.7.013006'>Phys. Rev. Research 7, 013006 (2025)</a>.</p>"+
            "</html>";

    /** Asks for the input SpotPattern file, presents a dialog to set the , undistorts the image with the given slice or the stack and shows the result.
     *  Note that 'slice' is 1-based, i.e., 1 <= slice < n.
     *  Currently, only floating-point stacks are supported.
     *  If an ivAnalyzer is supplied (non-null) and has a fit to the energy-dependent overall
     *  x/y offset, scale and angle, this correction is taken into account. */
    public void run(String arg) {
        if (spotPatternPath == null)
            spotPatternPath = LeedParams.getString(LeedParams.PATTERNFILE);
        if (spotPatternPath == null)
            spotPatternPath = "";
        String lastDirectory = OpenDialog.getLastDirectory();
        String spotPatternDir = (new File(spotPatternPath)).getParent();
        if (spotPatternDir != null && (new File(spotPatternDir)).exists())
            OpenDialog.setLastDirectory(spotPatternDir);

        // D I A L O G
        GenericDialog gd = new GenericDialog("Create Spot Pattern for Off-Normal Incidence...");
        gd.addFileField("Spot Pattern File (input)", spotPatternPath, 65);
        gd.addNumericField("Incidence angle vs. normal:", incidenceAngle, 0, 6, "\u00b0");
        gd.addNumericField("Incidence plane is at angle of:", angleVsRefSpot, 0, 6, "\u00b0 (ccw) ...");
        gd.addStringField("with respect to this spot:", referenceSpot);
        gd.addNumericField("Tolerance:", tolerance, 3, 6, "(typ. 0.01 to 0.1)");
        gd.setInsets(20,2,0);   //top, left, bottom
        gd.addMessage("\n");    //linefeed to have a multiline label with two lines
        final MultiLineLabel statusLabel = (MultiLineLabel)gd.getMessage();
        DialogListener dialogListener = new DialogListener() {
            public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
                String errorMessage = null;
                String spotPatternPath = gd.getNextString();
                if (!(new File(spotPatternPath).exists()))
                    errorMessage = "No such spot pattern file";
                LeedSpotPattern spotPattern = new LeedSpotPattern(spotPatternPath, /*keepSourceLines=*/true);
                if (errorMessage == null)
                    errorMessage = spotPattern.getErrorMessage();
                if (errorMessage == null && !spotPattern.hasEquivalentBeams())
                    errorMessage = "Error: Spot pattern has no equivalent beams, symmetry can't be reduced";

                double inAngle = gd.getNextNumber();
                if (errorMessage == null && !(inAngle >= 0 && inAngle < 90))
                    errorMessage = "Invalid incidence angle";
                double azAngle = gd.getNextNumber();
                if (errorMessage == null && !(azAngle>-360 && azAngle<360))
                    errorMessage = "Invalid incidence plane angle";
                String refSpotName = gd.getNextString();
                int refSpotIndex = spotPattern.getIndex(refSpotName);
                if (refSpotIndex < 0 || sqr(spotPattern.getKx(refSpotIndex)) + sqr(spotPattern.getKy(refSpotIndex)) < 1e-6)
                    errorMessage = "Error: Spot invalid or not in Spot pattern file";
                double tolerance = gd.getNextNumber();
                if (errorMessage == null && !(tolerance >=MIN_TOLERANCE && tolerance<=MAX_TOLERANCE))
                    errorMessage = "Invalid tolerance value";
                if (errorMessage != null) {
                    statusLabel.setText("\n"+errorMessage);
                    return false;
                }
                String status = makeGroups(spotPattern, /*applyToSpotPattern=*/false,
                        inAngle, azAngle, refSpotName, tolerance);
                statusLabel.setText(status);
                return true;
            }
        };
        gd.addDialogListener(dialogListener);
        dialogListener.dialogItemChanged(gd, null);
        gd.addHelp(HELP_STRING);
        gd.showDialog();

        OpenDialog.setLastDirectory(lastDirectory);
        if (gd.wasCanceled()) return;

        spotPatternPath = gd.getNextString();
        LeedSpotPattern spotPattern = new LeedSpotPattern(spotPatternPath, /*keepSourceLines=*/true);
        String errorMessage = spotPattern.getErrorMessage();
        if (errorMessage != null) {
            IJ.error("Spot pattern file error:\n"+errorMessage);
            return;
        }
        incidenceAngle = gd.getNextNumber();
        angleVsRefSpot = gd.getNextNumber();
        referenceSpot = gd.getNextString();
        tolerance = gd.getNextNumber();
        
        String status = makeGroups(spotPattern, /*applyToSpotPattern=*/true,
            incidenceAngle, angleVsRefSpot, referenceSpot, tolerance);
        if (status.contains("rror")) {
            IJ.error(status);
            return;
        }

        // O U T P U T   F I L E  D I A L O G
        lastDirectory = OpenDialog.getLastDirectory();
        File inputFile = new File(spotPatternPath);
        String inputFilename = inputFile.getName();
        String inputDir = inputFile.getParent();
        String outName = LeedUtils.removeExtension(inputFilename)+"_offNormal.csv";
        IJ.showStatus("Save off-normal spot pattern file as?");
        SaveDialog sd = new SaveDialog("Save off-normal spot pattern file as", inputDir, outName, ".csv");
        String outFileName = sd.getFileName();
        if (outFileName == null)
            return;
        IJ.showStatus("");
        String outDirectory = sd.getDirectory();
        OpenDialog.setLastDirectory(lastDirectory);

        // W R I T E   S P O T   P A T T E R N   F I L E
        ArrayList<String> initialCommentLines = spotPattern.getInitialCommentLines();
        
        initialCommentLines.add("# Modified for "+incidenceAngle+" deg off-normal incidence:");
        initialCommentLines.add("# Azimuth "+angleVsRefSpot+" deg w.r.t. "+referenceSpot+" tolerance="+tolerance);
        initialCommentLines.add("#");


        // S A V E   O U T P U T
        spotPattern.save(outDirectory+outFileName);
        if (spotPattern.getErrorMessage() != null)
            IJ.error(spotPattern.getErrorMessage());
    }

    /**
     * Creates the new list of groups, adds the appropriate output lines to the StringBuilder
     * and returns a string starting with "Error" in case of error, or a message on the
     * old/new number of groups.
     * @param spotPattern The original spot pattern. Must include the original lines of the file.
     * @param applyToSPotPattern When true, the spot pattern is modified with the new group numbers
     * @param incidenceAngle (vs. normal, in degrees)
     * @param angleVsRefSpot azimuthal angle of incidence plane, vs 'referenceSpot'
     * @param referenceSpot Direction to a spot as reference ofr azimuthal angle 'angleVsRefSpot'
     * @param tolerance determines how exactly spots must be symmetry-equivalent
     * @return error message or two-line diagnostic message on success
     */
    String makeGroups(LeedSpotPattern spotPattern, boolean applyToSpotPattern,
            double incidenceAngle, double angleVsRefSpot, String referenceSpot, double tolerance) {

        int refSpotIndex = spotPattern.getIndex(referenceSpot);
        if (refSpotIndex < 0 || sqr(spotPattern.getKx(refSpotIndex)) + sqr(spotPattern.getKy(refSpotIndex)) < 1e-6)
            return "Error: Spot "+referenceSpot+" is invalid or not in spot pattern";

        double refSpotAngle = Math.atan2(spotPattern.getKy(refSpotIndex), spotPattern.getKx(refSpotIndex));
        double incidenceAzmuth = refSpotAngle + Math.toRadians(angleVsRefSpot);
        double cosAngle = Math.cos(incidenceAzmuth);
        double sinAngle = Math.sin(incidenceAzmuth);

        // M A K E   S U B G R O U P S   F O R   N O N _ E Q U I V A L E N T   S P O T S
        int[] subGroupNumbers = new int[spotPattern.size()];
        boolean[] spotHandled = new boolean[spotPattern.size()];
        int[] groupLimits = spotPattern.getLowestAndHighestGroup();
        int[] nSubGroupsMinus = new int[-groupLimits[0]+1];
        int[] nSubGroupsPlus = new int[groupLimits[1]+1];
        for (int g=groupLimits[0]; g<=groupLimits[1]; g++) {
            int[] spotsForGroup = spotPattern.getAllSpotsForGroup(g);
            if (spotsForGroup.length == 0)
                continue;
            if (spotsForGroup.length == 1) {
                incrementNSubGroups(g, nSubGroupsMinus, nSubGroupsPlus);
                continue;
            }
            int currentSubGroup = 0;
            for (int i=0; i<spotsForGroup.length; i++) {
                int spotI = spotsForGroup[i];
                if (spotHandled[spotI])
                    continue;
                subGroupNumbers[spotI] = currentSubGroup;
                spotHandled[spotI] = true;
                double minProjection = getProjection(spotPattern, spotI, cosAngle, sinAngle);
                double maxProjection = minProjection;
                for (int j=0; j<spotsForGroup.length; j++) {    //search for others that may be equivalent
                    int spotJ = spotsForGroup[j];
                    if (spotHandled[spotJ])
                        continue;
                    double projection = getProjection(spotPattern, spotJ, cosAngle, sinAngle);
                    //if (g<4)IJ.log(spotPattern.getName(spotJ)+" vs. "+spotPattern.getName(spotI)+": proj="+(float)projection+" vs. "+(float)minProjection+"-"+(float)maxProjection+" incAngle="+incidenceAngle);
                    if ((projection - minProjection)*(incidenceAngle*incidenceAngle) <= tolerance
                            && (maxProjection - projection)*(incidenceAngle*incidenceAngle) <= tolerance) {
                        subGroupNumbers[spotJ] = currentSubGroup;
                        spotHandled[spotJ] = true;
                        if (projection < minProjection) {
                            minProjection = projection;
                            j = 0;
                        }
                        if (projection > maxProjection) {
                            maxProjection = projection;
                            j = 0;
                        }
                    }
                }
                currentSubGroup++;
                incrementNSubGroups(g, nSubGroupsMinus, nSubGroupsPlus);
            }
        }

        // N E W   G R O U P   N U M B E R S  & S T A T I S T I C S
        int[][] groupSizeHistograms = new int[2][MAX_SPOTS_PER_GROUP];    //number of spots per group size before and after symmetry reduction
        int[] newGroups = new int[spotPattern.size()];
        int nMinusGroups = 0;
        for (int i=0; i<nSubGroupsMinus.length; i++)
            nMinusGroups += nSubGroupsMinus[i];
        int currentNewGroup = -nMinusGroups;
        for (int g=groupLimits[0]; g<=groupLimits[1]; g++) {
            int[] spotsForGroup = spotPattern.getAllSpotsForGroup(g);
            groupSizeHistograms[0][spotsForGroup.length]++;               //statistics for original
            for (int sub=0; sub<getNSubGroups(g, nSubGroupsMinus, nSubGroupsPlus); sub++) {
                int nSpotsPerSubGroup = 0;
                for (int i=0; i<spotsForGroup.length; i++) {
                    int spotI = spotsForGroup[i];
                    if (subGroupNumbers[spotI] == sub) {
                        newGroups[spotI] = currentNewGroup;
                        nSpotsPerSubGroup++;
                    }
                }
                groupSizeHistograms[1][nSpotsPerSubGroup]++;              //statistics for new spotpattern
                currentNewGroup++;
            }
        }
        // A P P L Y
        if (applyToSpotPattern) {
            for (int i=0; i<spotPattern.size(); i++) {
                String origLine = spotPattern.getSourceLine(i);
                String[] items = origLine.split(",");
                if (items.length < 6)
                    return "Internal spot pattern file error, line too short:\n"+origLine;
                        int commaPos = 0;
                for (int j=0; j<5; j++)     //find the 5th comma, which precedes the group index
                    commaPos = origLine.indexOf(',', commaPos+1);
                String out = origLine.substring(0, commaPos+1); //everything before the group number
                out += LeedUtils.pad(newGroups[i], 5);          //the new group number
                commaPos = origLine.indexOf(',', commaPos+1);
                if (commaPos > 0)
                    out += origLine.substring(commaPos);        //everything behind the group number
                spotPattern.setSourceLine(i, out);
            }
        }
        return "Input: "+writeStatistics(groupSizeHistograms[0])
                +"\nOutput: "+writeStatistics(groupSizeHistograms[1]);
    }


    /** Increments the number of subgroups for the 'g'-th group (g may be negative or postitve):
     *  If g is negative, the array element at -g of nSubGroupsMinus is incremented,
     *  otherwise, the array element at +g of nSubGroupsPlus is incremented. */
    private static void incrementNSubGroups(int g, int[] nSubGroupsMinus, int[] nSubGroupsPlus) {
        if (g < 0)
            nSubGroupsMinus[-g]++;
        else
            nSubGroupsPlus[g]++;
    }

    /** Returns the number of subgroups for the 'g'-th group (g may be negative or postitve) */
    private static int getNSubGroups(int g, int[] nSubGroupsMinus, int[] nSubGroupsPlus) {
        if (g < 0)
            return nSubGroupsMinus[-g];
        else
            return nSubGroupsPlus[g];
    }

    /** Returns the projection of the Cartesian reciprocal-lattice vector on the direction given by the
     *  cosine and sine of an angle */
    private static double getProjection(LeedSpotPattern spotPattern, int spotI, double cosAngle, double sinAngle) {
        double kx = spotPattern.getKx(spotI);
        double ky = spotPattern.getKy(spotI);
        return kx*cosAngle + ky*sinAngle;
    }

    /** Returns a string describing the histogram of group sizes */
    private String writeStatistics(int[] groupSizeHistogram) {
        int nGroups=0;
        int smallestGrp=0, largestGrp=0;
        for (int i=2; i<groupSizeHistogram.length; i++) {
            nGroups += groupSizeHistogram[i];
            if (groupSizeHistogram[i]>0) {
                largestGrp = i;
                if (smallestGrp==0) smallestGrp=i;
            }
        }
        String s = groupSizeHistogram[1]+" single spot(s), "+nGroups+" Groups of size ";
        if (smallestGrp==largestGrp) s += smallestGrp;
        else {
            s += smallestGrp;
            if (groupSizeHistogram[smallestGrp]+groupSizeHistogram[largestGrp] == nGroups)
                s += "&";
            else
                s += "-";
            s += largestGrp;
        }
        return s;
    }

    static double sqr(double x) { return x*x; }
}
