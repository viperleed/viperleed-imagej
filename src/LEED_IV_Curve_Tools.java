import ij.*;
import ij.gui.*;
import ij.io.*;
import ij.macro.*;
import ij.plugin.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.File;

/** Various LEED I(V) Curve Utility Functions
 *  - data correction (as function of energy)
 *  - limit energy range
 *  - extract integer/superstructure/beams or beams (not) in other file
 *  - set group according to SpotPattern or delete group information
 *  - modify basis: matrix for h, k
 */


/** This ImageJ plugin is part of the ViPErLEED package for LEED I(V) analysis.
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


public class LEED_IV_Curve_Tools implements PlugIn, DialogListener {
    static final String PLUGIN_NAME = "LEED I(V) Curve Tools";
    /* main tasks */
    static final int N_TASKS = 5;
    static final int CORRECT_I=0, E_RANGE=1, EXTRACT=2, MODIFY_BASIS=3, CHANGE_GROUPS=4;
    static final String[] TASK_NAMES = new String[]{"Correct intensity values", "Limit energy range",
            "Extract beams", "Modify basis", "Modify groups of equivalent beams"};
    /* sub-tasks */
    static final int EXTRACT_INT=0, EXTRACT_SUP=1, EXTRACT_FILE_YES=2, EXTRACT_FILE_NO=3;
    static final String[] EXTRACT_NAMES = new String[]{"Integer beams", "Superstructure beams", "Beams present in file:", "Beams NOT in file:"};
    static final int GROUPS_FROM_PATTERN=0, GROUPS_INEQUIVALENT=1, GROUPS_DELETED=2;
    static final String[] CHANGE_GROUPS_NAMES = new String[]{"Set group numbers from spot pattern file:", "Make all beams inequivalent", "Delete group numbers"};
    /* The following should be remembered until ImageJ is closed */
    static String filePath;
    static boolean[] taskEnabled = new boolean[N_TASKS];
    static String macroCode = "Econv=80; I = I/(1-exp(-(E+Econv)/230));";
    static String energyRangeStr = "80-200";
    static int extractWhat;
    static String otherFilePath;
    static int changeGroupTask;
    static String spotPatternPath = LeedParams.getString(LeedParams.PATTERNFILE);
    static String spot10new = "0,1";
    static String spot01new = "1,0";
    /* GUI-related */
    ArrayList[] taskComponents = new ArrayList[N_TASKS];    //Components that are en/disabled if a task is en/disabled
    ArrayList[] taskComponents2 = new ArrayList[N_TASKS];   //Components that are conditionally en/disabled if a task is en/disabled
    private static final String IMAGEJ_MACRO_LANGUAGE_URL = "http://imagej.net/ij/developer/macro/macros.html";
    private static final String HELP_STRING =
            "<html>"+
            "<h1>"+PLUGIN_NAME+"</h1>"+
            "<p>This plugin provides utility functions for modifying <i>I</i>(<i>V</i>) data files (.csv format):</p>"+
            "<ul>"+
            "<li><b>"+TASK_NAMES[CORRECT_I]+"</b>: Modifies the intensity values. "+
            "The function for modifiying the intensities should be specified in the "+
            "<a href='"+IMAGEJ_MACRO_LANGUAGE_URL+"'>ImageJ macro language</a>. "+
            "Variables 'I' (current intensity value) and 'E' (energy) are provided, and the result should be assigned to the variable 'I'. "+
            "For simple expressions that do not contain '=' signs, 'if' clauses, etc., assignment to 'I' can be omitted, "+
            "e.g. '<tt>I/sqrt(E)</tt>' is equivalent to '<tt>I=I/sqrt(E)</tt>' and divides all intensities by the square root of the energy.<br />"+
            "An energy-dependent correction can be used, e.g., to account for energy-dependent detection efficiency of a "+
            "micro-channel plate in an MCP-LEED system. "+
            "A typical correction like this might be '<tt>Econv=80; I = I/(1-exp(-(E+Econv)/230));</tt>' "+
            "where '<tt>Econv=80</tt>' is the conversion voltage (the bias at the entrance side of the MCP). "+
            "Of course, you can also put the value of Econv directly into the equation, e.g. '<tt>I/(1-exp(-(E+80)/230));</tt>'."+
            "</li>"+
            "<li><b>"+TASK_NAMES[E_RANGE]+"</b>: Extracts the given energy range. "+
            "Specify the energy range as, e.g., &quot;<tt>80-200</tt>&quot;. "+
            "Both limits are included."+
            "</li>"+
            "<li><b>"+TASK_NAMES[EXTRACT]+"</b>: Extracts only certain beams and omits the others. "+
            "One can extract only integer beams, only superstructure beams, or the beams present (or not present) in a given file "+
            "[some other .csv file with <i>I</i>(<i>V</i>) curves]."+
            "</li>"+
            "<li><b>"+TASK_NAMES[MODIFY_BASIS]+"</b>: Applies a transformation matrix to the spot indices. "+
            "You can select the new indices of the (1,0) and (0,1) spots. "+
            "You can use this function, e.g., to swap the <i>h</i> and <i>k</i> indices: Change (1,0) to (0,1) and (0,1) to (1,0). "+
            "Swapping <i>h</i> and <i>k</i> is often needed for crystals with threefold symmetry, where (1,0) and (0,1) are not equivalent, "+
            "but is is initially not known which is which.<br>"+
            "If the spot labels have included group numbers in square brackets (to indicate the symmetry-equivalent beams, "+
            "see below), these group numbers will be deleted."+
            "</li>"+
            "<li><b>"+TASK_NAMES[CHANGE_GROUPS]+"</b>: An <i>I</i>(<i>V</i>) curve file typically contains the <i>group number</i> "+
            "for each beam. The group number is the number in square brackets in the name of the beam given in the heading. "+
            "Beams with the same group number are considered symmetry-equivalent.<br />"+
            "These group numbers can be set to the values of a given spot pattern file,\u00b9 "+
            "changed to a running index such that all beams are considered symmetry-inequivalent, "+
            "or, alternatively, the group numbers can be deleted. (In the latter case, "+
            "any future operations that require knowledge on symmetry-equivalent beams "+
            "will ask for a spot pattern file.\u00b9)"+
            "</li>"+
            "</ul>"+
            "<p>You can select one or more of the functions in the above list.</p>"+
            "<p>\u00b9A <i>spot pattern file</i> is a .csv file containing a list of the beams with additional information on them, "+
            "such as the group number, which defines groups of symmetry-equivalent beams. "+
            "A spot pattern file can be created with the ViPErLEED gui (select 'export').<br />"+
            "Note that replacing the group numbers with those from a given spot pattern file may also "+
            "change the beam names to those of the spot pattern file. E.g., if the beam is currently "+
            "named '0.5;1.0', its name may change to '1/2| 1 [8]' "+
            "(if it is given like this in the spot pattern file), where '[8]' is the group number. "+
            "(Beam indices may be given in decimal notation or as fractions; beams are consider identical "+
            "if both indices agree within 0.01.)</p>"+
            "<h2><a name='ivEditorLicense'>License</a></h2>"+
            "<p>The code is licensed under <a href='http://www.gnu.org/licenses/gpl-3.0.html'>GNU General Public License v3.0</a> "+
            "or later (GPL-3.0-or-later). "+
            "The authors may decide later to put part of the auxiliary code in this work into the public domain, "+
            "to allow incorporation into ImageJ if desired (ImageJ is in the public domain).</p>"+
            "<p>This documentation is licensed under the <a href='http://creativecommons.org/licenses/by/4.0/'>Creative Commons Attribution 4.0</a> "+
            "(CC BY 4.0) license.</p>"+
            "<p>When using this program (in its original or modified form) for scientific work, "+
            "please cite the paper describing the program [<a href='#paper'>1</a>].</p>"+
            "<p>You should find a copy of these license texts and the source code in the zip/jar archive holding this plugin "+
            "(use an unzip utility to view its contents).</p>"+
            "<h2>References</h2>"+
            "<p><a name='paper'>[1]</a> M. Schmid, F. Kraushofer, A. M. Imre, T. Kißlinger, L. Hammer, U. Diebold, and M. Riva, "+
            "<i>ViPErLEED package II: Spot tracking, extraction and processing of I(V) curves</i>, "+
            "Phys. Rev. Research, 2024. <a href='https://arxiv.org/abs/2406.18413/'>arXiv:2406.18413</a></p>"+
            "</html>";

    @SuppressWarnings("unchecked")
    public void run(String arg) {
        if (IJ.versionLessThan("1.52u")) return;
        for (int iTask=0; iTask<N_TASKS; iTask++) {
            taskComponents[iTask] = new ArrayList();
            taskComponents2[iTask] = new ArrayList();
        }
        // I N I T I A L   D I A L O G
        if (spotPatternPath==null)
            spotPatternPath = LeedParams.getString(LeedParams.PATTERNFILE);
        GenericDialog gd = new GenericDialog(PLUGIN_NAME);
        gd.addFileField("I(V) Curve File", filePath, 50);
        gd.setInsets(12, 0, 0);            //(top, left, bottom): at the very left 
        gd.addMessage("Action(s) to do:");
        for (int iTask=0; iTask<N_TASKS; iTask++) {
            gd.setInsets(5, 0, 0);            //(top, left, bottom): at the very left 
            gd.addCheckbox(TASK_NAMES[iTask], taskEnabled[iTask]);
            switch(iTask) {
                case CORRECT_I:
                    gd.addStringField("Macro code", macroCode, 50);
                    taskComponents[iTask].add(gd.getStringFields().lastElement());
                    taskComponents[iTask].add(gd.getLabel());
                    break;
                case E_RANGE:
                    gd.addStringField("Energy range", energyRangeStr, 50);
                    taskComponents[iTask].add(gd.getStringFields().lastElement());
                    taskComponents[iTask].add(gd.getLabel());
                    break;
                case EXTRACT:
                    gd.addChoice("Extract", EXTRACT_NAMES, EXTRACT_NAMES[extractWhat]);
                    taskComponents[iTask].add((Component)gd.getChoices().lastElement());
                    taskComponents[iTask].add((Component)gd.getLabel());
                    gd.addFileField("File defining the beams", otherFilePath, 50);
                    taskComponents2[iTask].add((Component)gd.getStringFields().lastElement());
                    taskComponents2[iTask].add((Component)gd.getLabel());
                    break;
                case MODIFY_BASIS:
                    gd.addStringField("Change (1,0) to", spot10new, 20);
                    taskComponents[iTask].add((Component)gd.getStringFields().lastElement());
                    taskComponents[iTask].add((Component)gd.getLabel());
                    gd.addStringField("Change (0,1) to", spot01new, 20);
                    taskComponents[iTask].add((Component)gd.getStringFields().lastElement());
                    taskComponents[iTask].add((Component)gd.getLabel());
                    break;
                case CHANGE_GROUPS:
                    gd.addChoice("Action", CHANGE_GROUPS_NAMES, CHANGE_GROUPS_NAMES[changeGroupTask]);
                    taskComponents[iTask].add((Component)gd.getChoices().lastElement());
                    taskComponents[iTask].add((Component)gd.getLabel());
                    gd.addFileField("Spot pattern file", spotPatternPath, 50);
                    taskComponents2[iTask].add((Component)gd.getStringFields().lastElement());
                    taskComponents2[iTask].add((Component)gd.getLabel());
                    break;
            }
        }
        gd.addMessage("You can drag and drop (.csv-)files into the 'file' fields.");

        gd.addDialogListener(this);
        Button okButton = gd.getButtons()[0];
        okButton.setEnabled(dialogItemChanged(gd, null));
        gd.addHelp(HELP_STRING);
        gd.showDialog();
        if (gd.wasCanceled()) return;

        filePath = gd.getNextString();
        double[] energyRange = null;
        double[][] basisMatrix = null;
        for (int iTask=0; iTask<N_TASKS; iTask++) {
            taskEnabled[iTask] = gd.getNextBoolean();
            switch(iTask) {
                case CORRECT_I:
                    macroCode = gd.getNextString();
                    break;
                case E_RANGE:
                    energyRangeStr = gd.getNextString();
                    energyRange = LeedUtils.rangeFromString(energyRangeStr);
                    break;
                case EXTRACT:
                    extractWhat = gd.getNextChoiceIndex();
                    String otherFilePathTmp = gd.getNextString();
                    if (extractNeedsFile(extractWhat))
                        otherFilePath = otherFilePathTmp;
                    break;
                case MODIFY_BASIS:
                    spot10new = gd.getNextString();
                    spot01new = gd.getNextString();
                    basisMatrix = makeBasisMatrix(spot10new, spot01new);
                    if (basisMatrix == null) {
                        IJ.error(PLUGIN_NAME, "Error: Invalid basis transformation");
                        return;
                    }
                    break;
                case CHANGE_GROUPS:
                    changeGroupTask = gd.getNextChoiceIndex();
                    String spotPatternPathTmp = gd.getNextString();
                    if (changeGroupTask == GROUPS_FROM_PATTERN)
                        spotPatternPath = spotPatternPathTmp;
                    break;
             }
        }
        // R E A D   D A T A
        LeedIVData ivData = LeedIVData.fromFile(filePath, /*spotPattern=*/null, /*flags=*/0);
        if (ivData == null) return;
        // P R O C E S S   D A T A
        boolean ok = true;
        for (int iTask=0; iTask<N_TASKS; iTask++) {
            if (taskEnabled[iTask]) {
                switch(iTask) {
                    case CORRECT_I:
                        correctIntensities(ivData, macroCode);
                        break;
                    case E_RANGE:
                        boolean success = ivData.cropEnergyRange(energyRange);
                        if (!success) {
                            IJ.error(PLUGIN_NAME, "Error: Empty data set when extracting energies "+
                                    IJ.d2s(energyRange[0])+"\u2013"+IJ.d2s(energyRange[1])+" from "+
                                    IJ.d2s(ivData.energies[0])+"\u2013"+IJ.d2s(LeedUtils.lastElement(ivData.energies)));
                            return;
                        }
                        break;
                    case EXTRACT:
                        LeedSpotPattern spotPattern = new LeedSpotPattern(ivData.spotNames, /*groupsRequired=*/false);
                        ivData.setSpotPattern(spotPattern, /*showWarning=*/true); //no warning expected anyhow
                        if (extractNeedsFile(extractWhat))
                            ok = extractFileBased(ivData, spotPattern, otherFilePath, extractWhat);
                        else
                            extractBeams(ivData, spotPattern, extractWhat);
                        if (ivData.data.length == 0) {
                            IJ.error(PLUGIN_NAME, "Error: no data extracted (empty data set)");
                            return;
                        }
                        break;
                    case MODIFY_BASIS:
                        ok = modifyBasis(ivData, basisMatrix);
                        break;
                    case CHANGE_GROUPS:
                        ok = changeGroups(ivData, spotPatternPath, changeGroupTask);
                        break;
                }
                if (!ok) return;
            }
        }
        // O U T P U T
        String fileName = (new File(filePath)).getName();
        String outName = LeedUtils.removeExtension(fileName)+"_.csv";
        IJ.showStatus("Save I(V) file as?");
        SaveDialog sd = new SaveDialog("Save I(V) file as", outName, ".csv");
        IJ.showStatus("");
        String outDirectory = sd.getDirectory();
        String outFileName = sd.getFileName();
        if (outFileName != null)
            ivData.save(outDirectory+outFileName);
    }

    /** Checks the dialog input for validity and returns false if invalid.
     *  Also enables/disables components */
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        boolean ok = true;                              //return value. (Don't return before components are enabled/disabled)
        String filePath = gd.getNextString();
        if (!LeedUtils.fileOk(filePath)) ok = false;
        boolean anyTaskSelected = false;
        for (int iTask=0; iTask<N_TASKS; iTask++) {
            boolean enabled = gd.getNextBoolean();      //whether this task is enabled
            boolean enabled2 = false;                   //whether to enable subtask-specific fields of this task
            if (enabled) anyTaskSelected = true;
            switch(iTask) {
                case CORRECT_I:
                    String macroCode = gd.getNextString();
                    if (enabled && !checkMacro(macroCode)) ok = false;
                    break;
                case E_RANGE:
                    String energyRangeStr = gd.getNextString();
                    double[] energyRange = LeedUtils.rangeFromString(energyRangeStr);
                    if (enabled && energyRange==null) ok = false;
                    break;
                case EXTRACT:
                    int extractWhat = gd.getNextChoiceIndex();
                    enabled2 = enabled && extractNeedsFile(extractWhat);
                    String otherFilePath = gd.getNextString();
                    if (enabled2 && !LeedUtils.fileOk(otherFilePath)) ok = false;
                    break;
                case CHANGE_GROUPS:
                    int changeGroupTask = gd.getNextChoiceIndex();
                    enabled2 = enabled && changeGroupTask == GROUPS_FROM_PATTERN;
                    String spotPatternPath = gd.getNextString();
                    if (enabled2 && !LeedUtils.fileOk(spotPatternPath)) ok = false;
                    break;
                case MODIFY_BASIS:
                    String spot10new = gd.getNextString();
                    String spot01new = gd.getNextString();
                    double[][] basisMatrix = makeBasisMatrix(spot10new, spot01new);
                    if (basisMatrix == null) ok = false;
            }
            for (Object component : taskComponents[iTask]) {    //enable/disable fields
                ((Component)component).setEnabled(enabled);
                if (component instanceof TextField) 
                    ((Component)component).setForeground(enabled ? Color.BLACK : Color.GRAY);
            }
            for (Object component : taskComponents2[iTask]) {
                ((Component)component).setEnabled(enabled2);
                if (component instanceof TextField) 
                    ((Component)component).setForeground(enabled2 ? Color.BLACK : Color.GRAY);
            }
        }
        return ok && anyTaskSelected;
    }

    /** Checks whether the 'extract' needs an 'other' file*/
    boolean extractNeedsFile(int extractWhat) {
        return extractWhat == EXTRACT_FILE_NO || extractWhat == EXTRACT_FILE_YES;
    }

    /** Checks whether the macro code looks reasonable at first glance */
    boolean checkMacro(String macroCode) {
        Program pgm = (new Tokenizer()).tokenize(macroCode);
        if (!pgm.hasWord("I")) {
            IJ.showStatus("Error: Macro code must contain intensity 'I'");
            return false;
        }
        Interpreter interp = new Interpreter();
        interp.setIgnoreErrors(true);       //don't show dialog box on errors
        macroCode = toFullMacroCode(macroCode);
        interp.run(macroCode, null);
        String error = interp.getErrorMessage();
		if (error != null) {
			IJ.showStatus("Macro code error: "+error);
            return false;
		} else {
            IJ.showStatus("");
            return true;
        }
    }

    /** Creates the full macro code needed by the Interpreter */
    String toFullMacroCode(String macroCode) {
        if (!macroCode.contains("="))
            macroCode = "I="+macroCode;     //allow also simplified syntax like "2.0*I" for "I=2.0*I"
        macroCode =
			"var E,I;\n"+
			"function dummy() {}\n"+        //11 tokens up to here. When changing, modify PCStart in correctIntensities, below.
			macroCode+";\n";
        
        return macroCode;
    }

    /** Change intensities based on macro code. Throws an exception in case of macro error. */
    void correctIntensities(LeedIVData ivData, String macroCode) {
        Interpreter interp = new Interpreter();
        macroCode = toFullMacroCode(macroCode);
		interp.run(macroCode, null);
        int PCStart = 11;                   
        for (int ic=0; ic<ivData.data.length; ic++) { //for all columns (spots)
            for (int i=0; i<ivData.energies.length; i++) {
                interp.setVariable("E", ivData.energies[i]);
                interp.setVariable("I", ivData.data[ic][i]);
                interp.run(PCStart);
                ivData.data[ic][i] = interp.getVariable("I");
            }
        }
    }

    /** Keeps only the beams present or absent in a given file */
    boolean extractFileBased(LeedIVData ivData, LeedSpotPattern spotPattern, String otherFilePath, int extractWhat) {
        LeedIVData ivData2 = LeedIVData.fromFile(otherFilePath, /*spotPattern=*/null, /*flags=*/0);
        if (ivData2 == null) return false;
        LeedSpotPattern spotPattern2 = new LeedSpotPattern(ivData2.spotNames, false);
        for (int ic=0; ic<ivData.spotNames.length; ic++) {
            String spotName = ivData.spotNames[ic];
            boolean isPresent = spotPattern2.getIndex(spotName) >= 0;
            if ((extractWhat == EXTRACT_FILE_NO && isPresent) ||
                    extractWhat == EXTRACT_FILE_YES && !isPresent)
                ivData.spotIndices[ic] = Integer.MIN_VALUE; //mark as invalid
        }
        ivData.deleteIfSpotIndex(Integer.MIN_VALUE);
        return true;
    }

    /** Keeps only the integer or superstructure beams */
    void extractBeams(LeedIVData ivData, LeedSpotPattern spotPattern, int extractWhat) {
        for (int ic=0; ic<ivData.spotIndices.length; ic++) {
            boolean isSuperstr = spotPattern.isSuperstructure(ivData.spotIndices[ic]);
            if ((extractWhat == EXTRACT_INT && isSuperstr) ||
                    extractWhat == EXTRACT_SUP && !isSuperstr)
                ivData.spotIndices[ic] = Integer.MIN_VALUE; //mark as invalid
        }
        ivData.deleteIfSpotIndex(Integer.MIN_VALUE);
    }

    /** Changes the group numbers in the column headings */
    boolean changeGroups(LeedIVData ivData, String spotPatternPath, int changeGroupTask) {
        switch (changeGroupTask) {
            case GROUPS_FROM_PATTERN:
                LeedSpotPattern spotPattern = new LeedSpotPattern(spotPatternPath);
                if (spotPattern.size() < 2) {
                    IJ.error("Error reading Spot Pattern", spotPattern.getErrorMessage());
                    return false;
                }
                ivData.setSpotPattern(spotPattern, /*showWarning=*/false);
                if (ivData.spotIndices == null) return false;
                for (int ic=0; ic<ivData.spotNames.length; ic++)
                    ivData.spotNames[ic] = spotPattern.getNameWithGroup(ivData.spotIndices[ic], /*replaceComma=*/true);
                break;
            case GROUPS_INEQUIVALENT:
            case GROUPS_DELETED:
                for (int ic=0; ic<ivData.spotNames.length; ic++) {
                    String spotName = ivData.spotNames[ic];
                    spotName = spotName.replaceAll("\\[.*\\]","").trim(); //remove spot group in []
                    if (changeGroupTask == GROUPS_INEQUIVALENT)
                        spotName += " ["+ic+"]";
                    ivData.spotNames[ic] = spotName;
                }
                break;
        }
        return true;
    }

    /** Modifies the basis by applying a matrix. */
    boolean modifyBasis(LeedIVData ivData, double[][] basisMatrix) {
        String[] spotNames = ivData.spotNames;
        for (int n=0; n<spotNames.length; n++) {
            double[] spotIndices = LeedUtils.spotIndices(spotNames[n]);
            double[] newIndices = new double[2];
            for (int i=0; i<2; i++)
                for (int j=0; j<2; j++)
                    newIndices[j] += basisMatrix[i][j]*spotIndices[i];
            spotNames[n] = toString(newIndices[0])+"|"+toString(newIndices[1]);
        }
        ivData.setSpotPattern(null, /*showWarning=*/false);    //remove previous spot pattern
        return true;
    }

    /** Returns the matrix needed to convert (1,0) into spot10new and (0,1) into spot01new.
     *  Returns null on error (also if the two vectors are collinear) */
    double[][] makeBasisMatrix(String spot10new, String spot01new) {
        double[][] matrix = new double[2][];
        for (int i=0; i<2; i++) {
            String newName = i==0 ? spot10new : spot01new;
            double[] spotIndices = LeedUtils.spotIndices(newName);
            if (spotIndices == null) return null;
            matrix[i] = spotIndices;
        }
        if (Math.abs(matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]) > 1e-6)
            return matrix;
        else
            return null;    //degenerate matrix or NaN
    }

    /** Converts a spot index to a String. For non-integers, the output is formatted as a fraction
     *  if it can be expressed as such with a denominator <= 99. */
    String toString(double x) {
        if (LeedUtils.isInteger(x))
            return Integer.toString((int)Math.round(x));
        for (int denom=2; denom<=99; denom++)
            if (LeedUtils.isInteger(x*denom))
                return (int)Math.round(x*denom)+"/"+denom;
        return IJ.d2s(x, 4);
    }

}
