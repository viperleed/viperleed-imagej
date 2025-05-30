import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.io.*;
import ij.plugin.*;
import ij.util.Tools;
import ij.measure.ResultsTable;
import ij.macro.*;
import ij.plugin.frame.Recorder;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.File;
import java.util.concurrent.atomic.*;

/** ViPErLEED spot tracker for extracting I(V) curves from LEED movies.
 *  See the LeedSpotTrackerHelp class for more info.
 *
 *  This is a singleton class (only one instance must run at a time).
 *  Also several classes called by the Spot Tracker rely on this being a singleton
 *  (e.g. when using the current values of the parameters in LeedParams like
 *  global variables) */

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
 *  @author Michael Schmid, IAP/TU Wien, 2019-2024
 * 
 *  The collection of ViPErLEED ImageJ plugins is free software:
 *  you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  The collection of ViPErLEED ImageJ plugins is distributed
 *  in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.

 */

public class LEED_Spot_Tracker implements PlugIn, DialogListener, ImageListener, ActionListener, WindowListener, MacroExtension, KeyListener {
    static final public String VERSION = "1.02";
    static final String PLUGIN_NAME = "LEED Spot Tracker";              //for dialog && error message titles
    static final String LOC_KEY_D = "leedSpotTracker.location";         //position of dialog
    static final String LOC_KEY_S = "leedSpotTracker.stackloc";         //position of SpotTracking stack
    static final String SCREENFIT_KEY = "leed.screenfit";               //in log.txt file
    // for input images / image stacks
    static final int N_IMP_FIELDS = 5;                  //number of image input Choices: input, dark, flat, dark2, mask
    static final int MAIN=0, DARK=1, FLAT=2, DARK2=3, MASK=4;          //image Choices in this sequence
    static final String[] IMP_NAMES = new String[] {"Main Input Stack",//names of input images
            "Dark Frame(s)", "Flat Field(s)", "Dark Frame for Flat F.", "Mask"};
    static final String[] IMP_SHORT_NAMES = new String[] {"iv", "dark", "flat", "dark2", "mask"};//short nemes for macro recording
    static final String[] IMP_SEARCH_NAMES = new String[] {"iv", "dark", "flat", "dark", "mask"};//when in image title, use it
    static final String NO_INPUT_STRING = "---";                        //used to select no input image in choice
    static final String STACK_SUFFIX = "_SpotTracking";                 //added to main input title, for the stack we work on
    // gui: buttons and labels next to them
    static final int N_BUTTONS = 8; //number of further buttons
    static final int FLAT_BUTTON=0, E_BUTTON=1, PATTERN_BUTTON=2, INDICES_BUTTON=3, RADIUS_BUTTON=4, TRACK_BUTTON=5, SAVE_BUTTON=6, MORE_BUTTON=7;
    static final String[] BUTTON_TITLES  = new String[] {"Dark&Flat Processing", "Set Energies, I0, t", "Set Pattern File", "Set Indices",
            "Set Integration Radius", "Track Spots", "Save Data", "More>>"};

    // Status flags, describe what has been achieved so far. These are bit masks for bitwise OR
    static final int MAIN_OK=1,    //we have a main LEED I(V) input image stack
            MASK_OK=2,              //we have a mask and maskRoi
            ENERGY_OK=4,            //energy range is known
            PATTERN_OK=8,           //we have a spot pattern file
            INDICES_OK=16,          //user has entered indices for spots
            TRACK_OK=32,            //spot tracking is done
            SAVE_OK=64;             //data are saved

    static final String[] MORE_MENU_TEXTS = new String[] {"Open Images/Movies...", "Select stack slice for...", "Create Mask...", //More>> popup menu entries
            "Highlight beams...", "Delete highlighted beams...", "Un-highlight beams", "Beam statistics",
            "Plot...", "Close plots...", "Undistort image/stack...",
            "List parameters", "Read parameters from log file...",
            "Compare parameters with log file...", "Metadata keywords...", "Reset all options & close"};
    static final int OPEN_FILE=0, SELECT_STACK_SLICE=1, CREATE_MASK=2,  //More>> menu entries
            HIGHLIGHT_SPOTS=3, DELETE_HIGHLIGHTED=4, UNHIGHLIGHT_SPOTS=5, SPOT_STATISTICS=6,
            PLOT_BEAM=7, CLOSE_PLOTS=8, UNDISTORT=9,
            LIST_PARAMETERS=10, READ_PARAMETERS=11,
            COMPARE_PARAMETERS=12, DATA_KEYWORDS=13, RESET_ALL_PARAMS=14;
    static final int[] POPUP_REQUIRED_STATUS = new int[] {               //required status for More>> popup items
            0, MAIN_OK|ENERGY_OK, MAIN_OK,
            TRACK_OK, TRACK_OK, TRACK_OK, TRACK_OK,
            TRACK_OK, 0, INDICES_OK,
            0, 0,
            0, 0, 0};

    static final String BADNESS_WARNING_TITLE = "WARNING: Beam(s) uncertain (highlighted)"; //title of table of bad spots
    //macro recording commands
    static final String MACRO_COMMAND_SETINPUT   = "LEEDsetInput";      //set input image or stack
    static final String MACRO_COMMAND_UPDATE     = "LEEDupdate";        //update 'stack' or 'meta'(data)
    static final String MACRO_COMMAND_SETVALUE   = "LEEDsetValue";      //set value for parameter
    static final String MACRO_COMMAND_SETARRAY   = "LEEDsetArray";      //set value for energiesEtc sub-array
    static final String MACRO_COMMAND_SETPATTERN = "LEEDsetPattern";    //set pattern file
    static final String MACRO_COMMAND_SELECTSLICE= "LEEDselectSliceForE"; //select Stack slice for this energy
    static final String MACRO_COMMAND_SETINDICES = "LEEDsetIndices";    //set ScreenFitter
    static final String MACRO_COMMAND_TRACKSPOTS = "LEEDtrackSpots";    //track spots
    static final String MACRO_COMMAND_SAVEDATA   = "LEEDsaveData";      //save data
    static final String MACRO_COMMAND_CLOSE      = "LEEDclose";         //close (what)
    static final String MACRO_COMMAND_GETSTATUS  = "LEEDgetStatus";     //get status(what)
    static final String MACRO_COMMAND_MAKETABLE  = "LEEDmakeTable";     //make a table of 'meta', 'index', 'spots'

    // gui-related
    static LEED_Spot_Tracker instance;
    NonBlockingGenericDialog gd;
    Choice[] impFieldChoices;       //Choices for input images/stacks
    Label[]  impFieldLabels     = new Label[N_IMP_FIELDS];     //labels for the input images/stacks
    ImagePlus[] impFieldImps    = new ImagePlus[N_IMP_FIELDS]; //ImagePlus selected for each image input field
    ImagePlus[] oldImpFieldImps = new ImagePlus[N_IMP_FIELDS]; //previous ImagePluses, to check for changes
    int[] impPreviousStackSize  = new int[N_IMP_FIELDS];       //the previous stack size for these images/stacks
    int[] impPreviousSlice      = new int[N_IMP_FIELDS];       //the slice previously displayed
    int[] impPreviousHashCode   = new int[N_IMP_FIELDS];       //previous hashcode, to check for changes
    int stackImpPreviousSlice   = -100;
    Button[] buttons = new Button[N_BUTTONS];
    Label[] buttonLabels = new Label[N_BUTTONS];
    //TextField[] textFields;
    PopupMenu popupMenu = new PopupMenu();
    Button okButton, cancelButton;
    Dialog currentFrontDialog;
    boolean componentsEnabled = true;                           //false when gui is blocked during lengthy operations
    static long lastRecorderHintsShownMillis = 0;               //for macro recording, when the comments have been shown last

    // describing the work
    // Metadata. When adding items here or changing the sequence, also modify LeedEnergyI0Selector.keyParams!
    static final int N_DATA_E_ETC = 10;    //size of energiesEtc array (number of data: ENERGY, TIME, CURRENT_I0, ...)
                                    //NOTE: when increasing this size, also modify LeedEnergyI0Selector.keyParams array!
    static final int ENERGY=0, TIME=1, CURRENTI0=2, CURRENTI00=3, TEMPERATURE=4,
            AUX=5, E_LEEM=6, PROCESSEDI0=7, BACKGROUND=8, BG_OVER_I0=9;  //indices in energiesEtc array
    public static final String[] ENERGY_ETC_NAMES = new String[] {"energy", "time", "I0", "I00", "temperature",
            "aux", "energy", "processedI0", "background intens.", "background/I0"};
    public static final String[] ENERGY_ETC_SHORTNAMES = new String[] {"E", "t", "I0", "I00", "T",
            "aux", "LEEM", "I0pr", "bgd", "bg/I0"};

    double[][] energiesEtc = new double[N_DATA_E_ETC][];    //energies, beam currents+offsets, times. Each starts with index 0 for stack slice 1
    int xAxisVariable = ENERGY;     //can be time, temperature, etc. if energy = constant
    double energyMin=LeedRadiusSelector.UNKNOWN_ENERGY, energyMax;  //only used for integration radius, not the energies in LEEM mode
    boolean anyImageChanged;        //makes it rebuild all input file choices
    int status;                     //flags, how far we are with all the input
    volatile boolean spotTrackerClosed;     //tells updateStackThread to stop
    ImagePlus stackImp;             //the (dark- & flat-processed) 'Spot tracking' stack for analysis (the stack that we work on)
    Thread updateStackThread;       //updates the stack for analysis in a separate thread
    AtomicInteger stackUpdateNeeded =
            new AtomicInteger();    //>0 means update needed; 2 tells updateStackThread that input images (except mask) have changed
    int maskCenterX, maskCenterY;   //center of active area in mask
    Roi maskRoi;                    //mask as a Roi
    LeedSpotPattern spotPattern;    //from spot pattern file, defines which spots we can have and their reciprocal-sapce coordinates
    LeedScreenFitter screenFitter;  //converts kx, ky to screen coordinates
    static double[] screenFitterArray; //remembers conversion between kx, ky and screen coordinates
    private Thread overlayThread;   //paints the overlay on the 'Spot tracking' stack
    private Object overlaySynchronizer = new Object(); //avoids concurrent call of overlay painting
    int xOfIndexInput;              //stack slice-1 where indices have been defined; spot tracking starts here
    LeedIVAnalyzer ivAnalyzer;      //gets data and holds it
    ArrayList<String> logLines;     //parameters for log file in human-readable form
    ArrayList<String> logParamLines;  //parameters for log file in machine-raeadable form
    String postTrackingChanges;     //non-null if changes occurred between tracking and 'save'
    ImagePlus indexInputImp;        //the stack slice used for index input, gets saved with the data

    // The following is for ImageJ Compile&Run only, to ensure everything is updated if it is not directly called here.
    // (we do not check for classes called by the LeedCurvePlotEditor, assuming these are up to date)
    LeedSpotAnalyzer dummy; LeedPlotter dummy2; LeedIndexSelector dummy3; LeedIntegerArray dummy4;
    LEED_Data_Quality_Statistics dummy5; LeedSmoother dummy6; LeedRFactor dummy7; Leed2DRegression dummy8;
    LeedScreenFitter dummy9; LeedAtomicFloatArray dummy10; LeedLinearRegression dummy11; LeedUtils dummy12;
    LeedOverlay dummy13; LeedDifferenceStack dummy14; Averaging_LEED_Stack dummy15; Linear_Fit_LEED_Stack dummy16;
    LeedFlatFitter dummy17; LeedDarkFlatVirtualStack dummy18; Open_LEED_Movie dummy19; Open_Aida_LEED_Video dummy20;
    LeedLabeledField dummy21; LeedImageSmoother dummy22; LeedEllipseEdgeFinder dummy23;

    //DEBUG static long t0;

// ----------------  T H E   G U I  ----------------

    /** This method is called by ImageJ when selecting "LEED Spot Tracker" from the Plugins menu.
     *  The method creates the GUI panel (unless we have it already on the screen). */
    @SuppressWarnings("unchecked")
    public void run(String arg) {
        if (IJ.versionLessThan("1.54k")) return;
        if (instance != null) {         //single-instance listener (otherwise params could be a mess)
            instance.gd.toFront();
            return;
        }
        instance = this;
        LeedParams.initialize();
        //DEBUG t0 = System.currentTimeMillis();
        LeedDarkFlatVirtualStack.prefetchAllowed = true;
        if (Prefs.blackBackground) {
            Prefs.blackBackground = false;
            String warning = "NOTE: "+PLUGIN_NAME+" has disabled 'Black Background' in Process>Binary>Options.";
            IJ.showStatus(warning);
            LeedUtils.logError(warning);
            LeedUtils.logError("note that the mask should be black in the (usable) screen area (pixel value 255),\n"+
                    "white outside (pixel value 0).");
        }
        if (Recorder.record) {
            recordMacro(null, null);    //show macro help
            Recorder.disableCommandRecording(); //don't record next GenericDialog
        }

        IJ.setForegroundColor(0,0,0);   //set fill/clear colors, as required for editing a mask
        IJ.setBackgroundColor(255,255,255);


        gd = new NonBlockingGenericDialog("LEED Spot Tracker");  // C R E A T E   G U I

        // Choices for input images
        String[] dummyStrings = new String[] {NO_INPUT_STRING};
        for (int i=0; i<N_IMP_FIELDS; i++) {
            gd.addChoice(IMP_NAMES[i], dummyStrings, null);
            impFieldLabels[i] = gd.getLabel();
        }
        Vector choices = gd.getChoices();
        impFieldChoices = (Choice[])(choices.toArray(new Choice[0]));
        for (int i=0; i<N_BUTTONS; i++)
            addButtonAndLabel(i);
        buttons[MORE_BUTTON].add(popupMenu);
        GUI.scalePopupMenu(popupMenu);
        gd.addHelp(LeedSpotTrackerHelp.getTrackerHelp());
        gd.setOKLabel("Done");
        Button[] okCancel = gd.getButtons();
        okButton = okCancel[0];
        cancelButton = okCancel[1];

        refreshImageChoices(null, /*autofill=*/true);
        if (impFieldImps[MAIN] != null)
            refreshImageChoices(null, true); //re-run if input image selected

        setPatternFile(LeedParams.getString(LeedParams.PATTERNFILE)); //must not be earlier because it can enable buttons
        enableAndHighlightComponents(true);

        gd.addDialogListener(this);
        ImagePlus.addImageListener(this);
        for (Button b : buttons)
            if (b != null)
                b.addActionListener(this);
        gd.addWindowListener(this);         //while sub-dialogs are active, we intercept the close box in windowClosing
        gd.removeWindowListener(gd);        //avoids disposing (we intercept closing while subdialog is open)
        Point loc = Prefs.getLocation(LOC_KEY_D);
        if (loc!=null) gd.setLocation(loc);
        ArrayList components = new ArrayList(gd.getChoices());
        Collections.addAll(components, buttons);
        Collections.addAll(components, gd.getButtons());
        for (Object c : components)         //keylistener to intercept <esc>
            if (c != null) {
                ((Component)c).removeKeyListener(gd);
                ((Component)c).addKeyListener(this);
            }
        EventQueue.invokeLater(new Runnable() {public void run() {
                    IJ.showStatus("");      //remove 'wait...' message when the stack is shown
                }});

        gd.showDialog();   // S H O W   D I A L O G   (everything happens during this phase)

        cleanup();
        if (flagSet(status, TRACK_OK) && !flagSet(status, SAVE_OK))
            saveData(null, null, true);
    }

    /** Returns the instance of the SpotTracker, or null if none or it does not have a stack or the dialog is closed */
    public static LEED_Spot_Tracker getInstance() {
        LEED_Spot_Tracker theInstance = instance;
        GenericDialog gd = theInstance == null ? null : theInstance.gd;
        return theInstance != null && theInstance.stackImp != null && gd != null && gd.isVisible() ?
                theInstance : null;
    }

    /** Things to do when the main dialog gets closed */
    void cleanup() {
        spotTrackerClosed = true;
        Prefs.saveLocation(LOC_KEY_D, gd.getLocation());
        LeedParams.saveToPrefs();
        ImagePlus.removeImageListener(this);
        LeedDarkFlatVirtualStack.prefetchAllowed = false;    //no  prefetch for any stacks that might be created while closing
        stopUpdateStackThread();
        LeedDarkFlatVirtualStack stack = getStack();
        if (stack != null)
            stack.terminatePrefetch();
        instance = null;
        ImageWindow stackWin = getStackWindow();
        if (stackWin != null && stackWin.isVisible()) {
            stackWin.removeWindowListener(this);
            stackWin.addWindowListener(stackWin); //revert: ImageJ becomes responsible for the close box
            stackWin.toFront();
        }
        LeedCloseDialog.showDialogAndClose(this, stackImp);
        closeBadnessWarning();
    }

    void addButtonAndLabel(int buttonIndex) {
        buttons[buttonIndex] = new Button(BUTTON_TITLES[buttonIndex]);
        buttonLabels[buttonIndex] = new Label("    ");
        Panel panel = new Panel();
        panel.setLayout(new FlowLayout(FlowLayout.LEFT, 0, 0));
        panel.add(buttons[buttonIndex]);
        panel.add(buttonLabels[buttonIndex]);
        gd.addPanel(panel);
    }

    /** Sets the text of a button label and updates the width.
     *  (in Swing would need to be in the EventQueue; here in awt not scrictly necessary, but does not hurt) */
    void setButtonLabelText(final int buttonIndex, final String text) {
        EventQueue.invokeLater(new Runnable() {public void run() {
                Label label = buttonLabels[buttonIndex];
                label.setText(text);
                label.setSize(label.getPreferredSize());
                gd.pack();
        }});
    }

// ----------------  I N P U T   I M A G E   S E L E C T I O N  ----------------

    /** Sets up the Choice menus for selecting the images/image stacks.
     *  If a closed or invalid image is supplied, also sets the corresponding 'impFieldImps' to null
     *  Automatically selects input images if 'autofill' is true. */
    @SuppressWarnings("unchecked")
    private void refreshImageChoices(ImagePlus closedImage, boolean autofill) {
        //IJ.log("refreshChoices auto="+autofill+" imp[0]="+impFieldImps[0]);
        int setToUnknown = -1;     //if an image has been closed, this entry must be changed (to unknown)
        if (closedImage != null) for (int i=0; i<N_IMP_FIELDS; i++) {
            if (closedImage == impFieldImps[i]) {
                setToUnknown = i;
                break;
            }
        }
        int[] imageIDs = WindowManager.getIDList();
        if (imageIDs == null) imageIDs = new int[0];
        ArrayList[] imps = new ArrayList[N_IMP_FIELDS]; //NB generices not allowed in arrays, no ArrayList<ImagePlus>[]
        for (int i=0; i<N_IMP_FIELDS; i++)
            imps[i] = new ArrayList<ImagePlus>(imageIDs.length+1);
        for (int id : imageIDs) {           // examine all open images whether they can be used
            ImagePlus imp = WindowManager.getImage(id);
            for (int type=0; type<N_IMP_FIELDS; type++)
                if (isGoodAsInput(imp, type, impFieldImps[MAIN]))
                    imps[type].add(imp);
        }
        for (int i=0; i<N_IMP_FIELDS; i++) {
            Choice choice = impFieldChoices[i];
            int nExtraOptions = (imps[i].size() == 0 || i != MAIN || i == setToUnknown) ?
                    1 : 0;  //in these cases we have to add a "---" = "nothing" option
            String[] titles = new String[imps[i].size() + nExtraOptions];
            if (nExtraOptions > 0) titles[0] = NO_INPUT_STRING;
            int newOptionSelected = -1;
            for (int jImp=0; jImp<imps[i].size(); jImp++) {             //check all listed images to see whether we still have the previous choice
                ImagePlus imp = (ImagePlus)(imps[i].get(jImp));
                String title = imp.getTitle();                          //and enter their names into the list for the choice
                titles[jImp+nExtraOptions] = title;
                if (newOptionSelected < 0 && imp == impFieldImps[i] && imp != closedImage)
                    newOptionSelected = jImp+nExtraOptions;             //keep previous if we have it
            }
            if (impFieldImps[i] == null && autofill) {                  //input unknown, autofill
                if (i == DARK2) {
                    if (impFieldImps[FLAT] == null)                     //no dark2 (for flat) if we have no flat
                        continue;
                    if (impFieldImps[FLAT].getProp(LeedDarkFlatVirtualStack.PROCESSED_FLAT_INFOKEY) != null)
                        continue;                                       //no autofill for dark2 if the flat has already the dark2 in it
                }
                int bestMatchQuality = -1;
                for (int jImp=0; jImp<imps[i].size(); jImp++) {
                    ImagePlus imp = (ImagePlus)(imps[i].get(jImp));
                    String title = imp.getTitle();

                    int guessedType = -1;  //negative: unknown
                    int matchQuality = 0;  //the higher the better
                    for (int i2=0; i2<IMP_SEARCH_NAMES.length; i2++) {
                        if (title.toLowerCase().indexOf(IMP_SEARCH_NAMES[i2])>=0) {
                            if (guessedType == -1 || (i != MAIN && guessedType == MAIN))
                                guessedType = i2;                       //e.g. "dark" in the name, "iv_dark" also counts as "dark" ("iv"= key for main)
                            else if (guessedType == DARK && i2 == FLAT) guessedType = -2; //"dark" AND "flat" in the name, uncertain type
                        }
                    }
                    if (guessedType == MAIN) {
                        matchQuality = 2;
                    } else if (guessedType < 0) { //unknown is ok for main
                        guessedType = MAIN;
                        matchQuality = 1;
                    } else if (impFieldImps[MAIN] != null)      //not main (flat, dark, mask): prefer if name starts like the name of main
                        matchQuality = LeedUtils.nMatchingBeginCharacters(impFieldImps[MAIN].getTitle(), title);
                    int wantedType = i==DARK2 ? DARK : i;
                    if (guessedType==wantedType && matchQuality > bestMatchQuality) {
                        newOptionSelected = jImp+nExtraOptions;
                        bestMatchQuality = matchQuality;
                    }
                }
                if (newOptionSelected >= nExtraOptions)
                    impFieldImps[i] = (ImagePlus)imps[i].get(newOptionSelected - nExtraOptions);
            }
            if (i == setToUnknown || newOptionSelected < 0) {
                newOptionSelected = 0;
                ImagePlus imp = (nExtraOptions > 0) ? null : (ImagePlus)imps[i].get(0);
                impFieldImps[i] = imp;
            }
            if (LeedUtils.setChoiceItems(choice, titles))   //update options for images
                gd.pack(); //if a change
            if (newOptionSelected >= 0)
                choice.select(newOptionSelected);           //select previous image or the one supplied by autofill
            if (choice.getSelectedIndex() == 0 && nExtraOptions > 0)
                impFieldImps[i] = null;                     //currently 'nothing' entry "---" is selected
            //IJ.log("refreshImageChoices "+IMP_SEARCH_NAMES[i]+" newSel="+newOptionSelected+" imp="+impFieldImps[i]);
        }
        checkForImpChanges();
    }

    /** Returns whether an image is useful as a given input type (MAIN, DARK... MASK), with a given main image */
    boolean isGoodAsInput(ImagePlus imp, int type, ImagePlus mainImp) {
        if (imp == null) return false;
        if (imp.getProcessor() == null) return false;       // in case the image has been closed asynchronously
        if (type != MAIN && mainImp == null) return false;
        if (imp == stackImp) return false;                  // the spot tracking stack (dark-&flat-corrected input) is not good as an input
        int bitDepth = imp.getBitDepth();
        if (bitDepth == 24) return false;                   // RGB images not supported
        int nSlices = imp.getImageStackSize();
        if (type == MAIN) return nSlices > 1;               // every other stack is good as main input
        if (type != MAIN && (imp.getWidth() != mainImp.getWidth() || imp.getHeight() != mainImp.getHeight()))
            return false; //For other images (dark, flat, mask) they must match in size
        if (type == MASK)
            return nSlices == 1 && bitDepth == 8 && imp.getProcessor().isBinary();  //mask must be binary
        else
            return true;    //dark, flat, dark2: no other requirements than what has been checked above.
    }

    /** Analyzes whether the mask is valid and gets its boundingbox */
    void analyzeMask() {
        ImagePlus maskImp = impFieldImps[MASK];
        ImagePlus mainImp = impFieldImps[MAIN];
        ImageProcessor maskIp = maskImp == null ? null : maskImp.getProcessor();
        //DEBUG Object maskPixels = maskIp == null ? null : maskIp.getPixels();
        //DEBUG if ((maskPixels instanceof byte[]) && impPreviousSlice[MASK] > 0 && Arrays.hashCode((byte[])maskPixels) == 
        //DEBUG       return;
        //DEBUG IJ.log("analyzeMask");
        if (maskImp == null || mainImp == null ||
                maskImp.getWidth() != mainImp.getWidth() || maskImp.getHeight() != mainImp.getHeight() ||
                !maskIp.isBinary() || maskImp.getImageStackSize() != 1) {
            killMask();
            return;
        }
        int foregroundValue = maskIp.isInvertedLut() ? 255 : 0;
        maskIp.setThreshold(foregroundValue, foregroundValue, ImageProcessor.NO_LUT_UPDATE);
        if (nForegroundCorners((ByteProcessor)maskIp, foregroundValue) > 1) {
            String warning = "Warning: Usually a mask should have background (white) pixels at the corners.";
            IJ.showStatus(warning);
            LeedUtils.logError(warning+"\nThis one has doesn't: "+maskImp.getTitle());
        }
        maskRoi = ij.plugin.filter.ThresholdToSelection.run(maskImp);
        if (maskRoi == null) {
            IJ.error("Error: Mask is empty");
            killMask();
            return;
        }
        Rectangle bounds = maskRoi.getBounds();
        maskCenterX = bounds.x + bounds.width/2;
        maskCenterY = bounds.y + bounds.height/2;
        changeStatus(MASK_OK, TRACK_OK | SAVE_OK);
        LeedDarkFlatVirtualStack stack = getStack();
        if (stack != null)
            triggerUpdateStack(false);
    }

    /** Deselects the mask (if any) */
    void killMask() {
        maskRoi = null;
        impPreviousSlice[MASK] = -1;
        LeedDarkFlatVirtualStack stack = getStack();
        if (stack != null)
            triggerUpdateStack(false);

        changeStatus(0, MASK_OK | INDICES_OK | TRACK_OK | SAVE_OK);
        refreshImageChoices(impFieldImps[MASK], false);    //mask is not usable (any more)
        impFieldImps[MASK] = null;
    }

    /** Counts the number of foreground pixels of the mask.
     *  Requires that maskIp is binary, i.e., all pixel values must be either 0 or 255. */
    int nForegroundCorners(ByteProcessor maskIp, int foregroundValue) {
        int xmax = maskIp.getWidth() - 1;
        int ymax = maskIp.getHeight() - 1;
        int n255corners = (maskIp.getPixel(0, 0) + maskIp.getPixel(0, ymax) +
                maskIp.getPixel(xmax, 0) + maskIp.getPixel(xmax, ymax))/255;
        return foregroundValue == 255 ? n255corners : 4 - n255corners;
    }
// ----------------  T H E   S T A C K  ----------------

    /** Creates a stack or updates the SpotTracking stack (i.e., the dark&flat-corrected input)
     *  if the input changes.
     *  Contains a synchronized section to avoid creating two SpotTracking stacks */
    void createOrUpdateStack() {
        ImagePlus mainImp = impFieldImps[MAIN];
        if (mainImp == null) return;
        for (int i=1; i<N_IMP_FIELDS; i++) {
            ImagePlus imp = impFieldImps[i];
            if (!isGoodAsInput(imp, i, mainImp)) {
                impFieldImps[i] = null;
                impFieldChoices[i].select(NO_INPUT_STRING);
            }
        }
        LeedDarkFlatVirtualStack stack = getStack();
        synchronized(this) {
            boolean makeNewStack = stack == null || stackImp == null || stack.getMainImp() != mainImp;
            if (makeNewStack)
                createNewStack();
            else
                triggerUpdateStack(true);
        }

        Object inputInfo = mainImp == null ? null :
            mainImp.getProperty("Info");    //transfer 'info' property
        String info = inputInfo == null ? "" : (String)inputInfo;
        if (info.length() > 0 && !info.endsWith("\n")) info += "\n";
        if (impFieldImps[DARK] != null)
            info += "Dark Frame: "+impFieldImps[DARK].getTitle()+"\n";
        if (impFieldImps[FLAT] != null) {
            info += "Flat Field: "+impFieldImps[FLAT].getTitle()+"\n";
            if (impFieldImps[DARK2] != null)
                info += "Dark for Flat: "+impFieldImps[DARK2].getTitle()+"\n";
        }
        ImagePlus stackImp = this.stackImp;
        if (info.length() > 0)
            stackImp.setProperty("Info", info);

        readEnergiesEtc(stackImp.getStack());
        drawOverlays();
        setButtonLabelText(FLAT_BUTTON, getFlatProcessingString());
        for (int i=0; i<N_IMP_FIELDS; i++)
            updateImpHashEtc(i);
    }

    /** Creates a new stack for analysis from main input, dark, flat, and sets the output image */
    void createNewStack() {
        ImagePlus mainImp = impFieldImps[MAIN];
        if (mainImp == null) return;
        ImagePlus stackImp = this.stackImp;
        if (stackImp != null) {     //close the old stack
            Thread updateStackThread = this.updateStackThread;
            if (updateStackThread != null) {
                updateStackThread.interrupt();
                this.updateStackThread = null;
            }
            ImageStack oldStack = stackImp.getStack();
            if (oldStack != null && (oldStack instanceof LeedDarkFlatVirtualStack))
                ((LeedDarkFlatVirtualStack)oldStack).terminatePrefetch();
            ImageWindow win = stackImp.getWindow();
            if (win != null) {
                win.removeWindowListener(this);
                Point loc = win.getLocation();
                Prefs.saveLocation(LOC_KEY_S, loc);
            }
            stackImp.changes = false;
            stackImp.close();
            this.stackImp = null;           //avoids re-creating the stack when it is closed
        }

        LeedDarkFlatVirtualStack stack = new LeedDarkFlatVirtualStack(mainImp);
        stackUpdateNeeded.set(2);
        startUpdateStackThread();

        String title = LeedUtils.removeExtension(impFieldImps[MAIN].getTitle())+STACK_SUFFIX;
        title = WindowManager.makeUniqueName(title);

        Recorder.suspendRecording();        //avoids recording 'selectWindow'
        stackImp = stack.getImagePlus(title);
        stackImp.show();
        this.stackImp = stackImp;
        ImageWindow win = stackImp.getWindow();
        //DEBUG IJ.log("createNewStack win="+win);
        //DEBUG try{throw new RuntimeException("");}catch(Exception e){IJ.log(" from\n"+e.getStackTrace()[1]+'\n'+e.getStackTrace()[2]+(e.getStackTrace().length>3?"\n"+e.getStackTrace()[3]:""));}

        if (win != null) {
            Point loc = Prefs.getLocation(LOC_KEY_S);
            if (loc != null) win.setLocation(loc);
            win.addWindowListener(this);
            win.removeWindowListener(win);  //disable close box; the Spot Tacker takes over as windowListener
        }
        if (Recorder.record)
            stackImp.waitTillActivated();
        Recorder.resumeRecording();
    }

    /** Triggers updating the stack for analysis in case input images/stacks have changed
     *  (newInput=true) or the processing has changed (newInput = false)
     *  This function does not wait until the stack for analysis is ready.
     *  Therefore, functions that work on the stack must call waitForStackUpdate. */
    void triggerUpdateStack(boolean newInput) {
        //DEBUG IJ.log("triggerUpdateStack new="+newInput+" flat="+impFieldImps[FLAT]);
        //DEBUG try{throw new RuntimeException("");}catch(Exception e){IJ.log(" from\n"+e.getStackTrace()[1]+'\n'+e.getStackTrace()[2]+(e.getStackTrace().length>3?"\n"+e.getStackTrace()[3]:""));}

        ImagePlus stackImp = this.stackImp;
        if (stackImp == null) return;
        if (!spotTrackerClosed && (updateStackThread == null || !updateStackThread.isAlive()))
            startUpdateStackThread();
        synchronized(updateStackThread) {
            if (newInput)
                stackUpdateNeeded.set(2);
            else
                stackUpdateNeeded.compareAndSet(0, 1); //won't change 2 to 1
            LeedDarkFlatVirtualStack stack = getStack();
            if (stack == null)
                return;
            stack.isChanging = true;
            updateStackThread.interrupt();
        }
    }

    /** Starts the background thread for updating the stack for analysis in case
     *  input images/stacks have changed (stackUpdateNeeded = 2) or
     *  the processing has changed (stackUpdateNeeded = 1) */
    void startUpdateStackThread() {
        updateStackThread = new Thread(new Runnable() {
            public void run() {
                try {
                    while (true) {
                        boolean newInput;
                        //DEBUG IJ.log("updateStackThread waiting t="+(System.currentTimeMillis()-t0));
                        synchronized(Thread.currentThread()) {
                            Thread.currentThread().interrupted();   //clear interrupted status
                            if (stackUpdateNeeded.get() == 0) {
                                LeedDarkFlatVirtualStack stack = getStack();
                                if (stack != null) stack.isChanging = false;
                                try {
                                    Thread.currentThread().wait();
                                } catch (InterruptedException e) {};
                            }
                            newInput = stackUpdateNeeded.get() > 1;
                            stackUpdateNeeded.set(0);
                        }
                        if (spotTrackerClosed) return;
                        //DEBUG IJ.log("updateStackThread notified, newInput="+newInput+" t="+(System.currentTimeMillis()-t0));
                        updateStack(newInput); //will reset stackUpdateNeeded
                        //DEBUG IJ.log("updateStackThread updateStack done t="+(System.currentTimeMillis()-t0));
                    }
                }
                catch (Exception e) {IJ.handleException(e);}
            }
        }, "update LEED_IV_stack");
        updateStackThread.start();
        //DEBUG IJ.log("updateThread started "+updateStackThread+" t="+(System.currentTimeMillis()-t0));
    }

    /** Stops the background thread for updating the stack for analysis */
    void stopUpdateStackThread() {
        Thread updateStackThread = this.updateStackThread;
        if (updateStackThread != null) {
            updateStackThread.interrupt();
            this.updateStackThread = null;
        }
    }

    /** Called upon changes of any of the input images/stacks (new/other/no image/stack)
     *  or the processing method.
     *  Called from the updateStackThread.
     *  The relevant methods of LeedDarkFlatVirtualStack are synchronized to avoid problems
     *  with non-synchronous execution */
    void updateStack(boolean newInput) {
        //IJ.log("updateStack "+stackImp);
        LeedDarkFlatVirtualStack stack = getStack();
        if (stack == null) return;
        ImagePlus mainImp = impFieldImps[MAIN];
        if (mainImp == null) return;
        stack.isChanging = true;
        setButtonLabelText(FLAT_BUTTON, "wait...");
        ImagePlus theStackImp = stackImp;
        //DEBUG IJ.log("updateStack new="+newInput+" flat="+impFieldImps[FLAT]);
        if (newInput) {
            stack.setInputDarkFlatDark2Mask(mainImp, impFieldImps[DARK], impFieldImps[FLAT], impFieldImps[DARK2],
                    impFieldImps[MASK],
                    (int)LeedParams.get(LeedParams.DARKAVERAGING),
                    (int)LeedParams.get(LeedParams.FLATAVERAGING),
                    (int)LeedParams.get(LeedParams.DARK2AVERAGING),
                    (int)LeedParams.get(LeedParams.FLATFITORDER), theStackImp);
        } else {
            stack.setAveragingAndPoly((int)LeedParams.get(LeedParams.DARKAVERAGING),
                    (int)LeedParams.get(LeedParams.FLATAVERAGING),
                    (int)LeedParams.get(LeedParams.DARK2AVERAGING),
                    (int)LeedParams.get(LeedParams.FLATFITORDER),
                    impFieldImps[MASK], theStackImp);
        } //these reset stack.isChanging()
        if (Thread.currentThread().isInterrupted()) return;
        stack.isChanging = false;
        drawOverlays();
        if (Thread.currentThread().isInterrupted()) return;
        setButtonLabelText(FLAT_BUTTON, getFlatProcessingString()); //update the gui in the EventQueue (not really needed in awt, only in swing)
        if (theStackImp != null)
            updateContrast(theStackImp, impFieldImps[MASK]);
    }

    /** Asks the user for patience while the stack is updated,
     *  otherwise returns true immediately.
     *  Returns false if the user presses 'cancel' or no foreground stack exists.
     *  Must NOT be called from the EventQueue */
    boolean waitForStackUpdate() {
        LeedDarkFlatVirtualStack stack = getStack();
        if (stack == null) return false;
        if (stack.isChanging)
            return LeedWaitForStackDialog.wait(this, stack);
        else
            return true;
    }

    /** Returns the stack, or null if not available, in a multithreading-safe manner */
    LeedDarkFlatVirtualStack getStack() {
        ImagePlus stackImp = this.stackImp;
        if (stackImp == null)
            return null;
        ImageStack stack = stackImp.getStack();
        if (!(stack instanceof LeedDarkFlatVirtualStack) || stack.size()==0)
            return null;                            //(after disposal, imp.getStack would return a newly created stack)
        return (LeedDarkFlatVirtualStack)stack;
    }

    /** Returns the window of the Spot tracking stack, or null, in a multithreading-safe manner */
    ImageWindow getStackWindow() {
        ImagePlus stackImp = this.stackImp;
        if (stackImp == null)
            return null;
        return stackImp.getWindow();
    }

    /** Changes brightness&contrast of the main stack, if it is too bad.
     *  The main image stack imp must be floating point type and the mask a binary
     *  byte image (background 0, foreground 255) */
    static void updateContrast(ImagePlus imp, ImagePlus maskImp) {
        byte[] maskPixels = LeedUtils.getMaskPixels(maskImp);
        if (maskPixels == null) return;
        FloatProcessor fp = imp == null ? null : (FloatProcessor)imp.getProcessor();
        float[] pixels = fp == null ? null : (float[])fp.getPixels();
        if (pixels == null) return;
        double sum=0;
        int count=0;
        float min = Float.MAX_VALUE, max = Float.MIN_VALUE;
        for (int p=0; p<pixels.length; p++) {
            if (maskPixels[p] == 0) continue;
            if (pixels[p] > max) max = pixels[p];
            if (pixels[p] < min) min = pixels[p];
            sum += pixels[p];
            count++;
        }
        double avg = sum/count;
        double currentMin = fp.getMin(), currentMax = fp.getMax();
        double currentDiff = currentMax - currentMin;

        boolean updateMin = currentMin < 1.2*min - 0.2*avg || currentMin > 0.7*avg + 0.3*max;
        boolean updateMax = currentMax < 0.5*min+0.5*avg || currentMax > 1.5*max - 0.5*avg;
        if (updateMin) currentMin = 1.1*min - 0.1*avg;
        if (currentMax - currentMin < currentDiff) currentMax = currentMin + currentDiff;
        if (updateMax && currentMax < currentMin + currentDiff) currentMax = currentMin + currentDiff;
        if (currentMax < 0.5*min+0.5*avg) currentMax = 0.5*min+0.5*avg;
        if (currentMax > 1.5*max - 0.5*avg) currentMax = 1.1*max - 0.1*avg;
        if (updateMin || updateMax) {
            fp.setMinAndMax(currentMin, currentMax);
            //IJ.log("updateMinMax min,avg,max="+(float)min+","+(float)avg+","+(float)max+" new:"+(float)currentMin+"-"+(float)currentMax);
            imp.updateAndDraw();
        }
    }

    /** Creates the overlay with the spot names */
    void drawOverlays() {
        synchronized (overlaySynchronizer) {
            final ImagePlus stackImp = this.stackImp;
            final Thread previousOverlayThread = overlayThread;
            if (stackImp == null) return;
            overlayThread = new Thread(
                new Runnable() {    //after (accidentally) closing and re-opening the stack, we have to wait until it is there
                    final public void run() {
                        if (previousOverlayThread != null && previousOverlayThread.isAlive()) {
                            previousOverlayThread.interrupt();
                            try {
                                previousOverlayThread.join();
                            } catch (InterruptedException e) {return;}
                        }
                        try {
                            Thread.sleep(100);
                        } catch (InterruptedException e) {return;}
                        LeedIVAnalyzer ivAnalyzerUsed = ivAnalyzer;
                        if (flagSet(status, TRACK_OK) && ivAnalyzerUsed != null) {
                            ivAnalyzerUsed.showOverlay(stackImp);   //also sets energy&mask overlays
                        } else {
                            drawBasicOverlay();
                            stackImp.draw();
                        }
                    }
                }, PLUGIN_NAME+"_makeOverlays");
            overlayThread.start();
        }
    }

    /** Creates the basic overlay with mask, energy and spot shape */
    void drawBasicOverlay() {
        if (flagSet(status, ENERGY_OK) || flagSet(status, MASK_OK))
            LeedOverlay.add(stackImp, flagSet(status, ENERGY_OK) ? energiesEtc[xAxisVariable] : null, xAxisVariable==ENERGY,
                    flagSet(status, MASK_OK) ? maskRoi : null);
    }        

    void readEnergiesEtc(final ImageStack stack) {
        ImageStack darkStack = impFieldImps[DARK] != null ? impFieldImps[DARK].getStack() : null;
        double[][] tempData = LeedEnergyI0Selector.decode(this, stack, darkStack);
        LeedEnergyI0Selector.tempDataToFinal(tempData, energiesEtc, stack.size(), null);   //tempSource=null means 'from File'
        updateEnergiesEtc(false, false);  //more checking, sets energyMin,energymax, sets status, shows it next to button
    }

    /** Sets the independent variable, ENERGY, TIME etc.
     *  Call updateEnergiesEtc thereafter. */
    void setXAxisVariable(int xAxisVariable) {
        this.xAxisVariable = xAxisVariable;
    }

    // ----------------  A C T I O N   B U T T O N   &   M E N U   R E S P O N S E  ----------------

    void askForFlatProcessing() {
        final LeedDarkFlatVirtualStack stack = getStack();
        final ImagePlus stackImp = this.stackImp;
        if (stack == null) return;  //should never happen
        int status = stack.doProcessingDialog();
        if (status > 0) {           //the dialog may have led to changes
            IJ.showProgress(0.0);
            triggerUpdateStack(false);
            if (status > 1) {       //show processed flat field
                IJ.showStatus("Showing flat field...");
                new Thread(new Runnable() {
                    public void run() {
                        try{
                            IJ.wait(100);
                            if (!waitForStackUpdate()) return;
                            ImagePlus imp = stack.getProcessedFlat();
                            if (imp == null)
                                return;
                            if (stackImp != null)
                                imp.setCalibration(stackImp.getCalibration());
                            imp.show();
                            IJ.showProgress(1.0);
                            IJ.showStatus("");
                        } catch (Exception e) { IJ.handleException(e); };
                    }
                }, "Show processed flat field").start();
            }
        }
        if (Recorder.record)
            recordMacro(MACRO_COMMAND_UPDATE, "'stack'");
    }

    /** The dialog for energies, I0, etc. */
    void askForEnergies() {
        final LeedDarkFlatVirtualStack stack = getStack();
        if (stack == null) return;
        final LEED_Spot_Tracker spotTracker = this;
        new Thread(
            new Runnable() {    //separate thread, so this is not on the EventQueue and we can have a non-blocking dialog
                final public void run() {
                    spotTracker.enableAndHighlightComponents(false);
                    ImageStack darkStack = impFieldImps[DARK] == null ? null : impFieldImps[DARK].getStack();
                    (new LeedEnergyI0Selector()).runSelectionDialog(spotTracker, stack, darkStack, energiesEtc);
                }
            }, PLUGIN_NAME+"_setEnergiesI0etc"
        ).start();  //when done, it calls updateEnergiesEtc and enableAndHighlightComponents(true);
    }

    /** Asks the user for a file with the spot pattern */
    void askForPatternFile() {
        LeedSpotPattern spotPattern = LeedSpotPattern.openWithDialog(LeedParams.getString(LeedParams.PATTERNFILE));
        if (spotPattern == null) return;
        if (spotPattern.size() < 3) {
            IJ.error(PLUGIN_NAME, "Error: LEED Pattern has <3 beams");
            return;
        }
        LeedParams.setString(LeedParams.PATTERNFILE, spotPattern.getPath());
        this.spotPattern = spotPattern;
        updateSpotPattern();
        if (Recorder.record)
            recordMacro(MACRO_COMMAND_SETPATTERN, '"'+spotPattern.getPath()+'"');
    }

    /** Asks the user to selects spot indices to determine the basis
     *  and saves it in the parameters */
    void askForIndices() {
        LeedIndexSelector indexSelector = new LeedIndexSelector(this, stackImp, impFieldImps[MASK], maskRoi,
                energiesEtc, xAxisVariable, spotPattern, screenFitterArray);
        enableAndHighlightComponents(false); // block all buttons&choices
        new Thread(indexSelector).start();   // separate thread, because it is nonmodal and requires user interaction
    }

    void askForRadius() {
        new LeedRadiusSelector(this, stackImp, energyMin, energyMax, spotPattern.isSuperstructure()).run();
        showRadii();
        drawOverlays();
    }

    /** Tracks the spots, analyzes the spots over all slices and creates output tables including I(V) curves*/
    void analyzeIV(boolean interactive) {
        closeBadnessWarning();
        energiesEtc[PROCESSEDI0] = LeedEnergyI0Selector.getProcessedI0(energiesEtc);
        if (energiesEtc[PROCESSEDI0] != null) {
            int nNonPositive = LeedUtils.countNonPositives(energiesEtc[PROCESSEDI0]);
            if (nNonPositive > 0) {
                IJ.error(PLUGIN_NAME, "Error: (processed) I0 has "+nNonPositive+" non positive values");
                changeStatus(0, ENERGY_OK);
                return;
            }
        }
        ivAnalyzer = new LeedIVAnalyzer(this, spotPattern, stackImp, impFieldImps[MASK], maskRoi,
                energiesEtc[ENERGY], energiesEtc[PROCESSEDI0], energiesEtc[xAxisVariable], ENERGY_ETC_NAMES[xAxisVariable],
                screenFitter, xOfIndexInput, interactive);
        enableAndHighlightComponents(false);    // block all buttons&choices
        if (interactive)
            new Thread(ivAnalyzer).start();     // separate thread, because it is slow (keep ImageJ responsive)
        else
            ivAnalyzer.run();
    }

    /** Saves the data. In interactive mode, 'dir' = null.
     *  If 'askToDiscard' is true, asks the user whether to discard the data
     *  until saving was successful or the user has pressed 'cancel'
     *  Returns null if successful, otherwise an error text starting with ERROR: */
    String saveData(String dir, String prefix, boolean askToDiscard) {
        boolean interactive = dir==null;
        if (dir == null || dir.length() ==0)
            dir = LeedParams.getString(LeedParams.SAVEDIRECTORY);
        if (!interactive && !(new File(dir)).isDirectory())
            return "ERROR: No such directory: "+dir;
        if (interactive && !askToDiscard && postTrackingChanges != null) {
            YesNoCancelDialog yncDialog = new YesNoCancelDialog(getStackWindow(), PLUGIN_NAME,
                    "Warning: "+postTrackingChanges+"\nchanged after spot tracking.",
                    /*yesLabel*/"Track again", /*noLabel*/"Save old tracking results");
            if (yncDialog.cancelPressed()) {
                return null;
            } else if (yncDialog.yesPressed()) {
                analyzeIV(/*interactive*/true);
                return null;
            }
        }
        if (prefix == null || prefix.length() == 0)
            prefix = LeedUtils.removeExtension(impFieldImps[MAIN].getTitle());
        LeedDataSaver dataSaver = new LeedDataSaver(this, dir, prefix, ivAnalyzer,
            energiesEtc[xAxisVariable], ENERGY_ETC_NAMES[xAxisVariable], spotPattern, stackImp, indexInputImp, askToDiscard);
        if (interactive) {
            Thread thread = new Thread(dataSaver);          // separate thread, saving stack is slow.
            thread.start();
        } else
            dataSaver.runInMacro();
        return null;
    }

    /** Dialog to open images or movied as input; when appropriate, these will be entered it into the relevant field(s) */
    void openFiles() {
        LeedFileOpenDialog fileOpenDialog = new LeedFileOpenDialog(this);
        enableAndHighlightComponents(false);    // block all buttons&choices
        new Thread(fileOpenDialog).start();     // separate thread (otherwise event queue, no progress bar, no IJ.showStatus)
    }

    /** Called by the LeedFileOpenDialog upon completion.
     *  imps is the list of ImagePlus objects (main, dark, ... mask) or null. */
    void setOpenFiles(ImagePlus[] imps) {
        if (imps == null) return;
        boolean changes = false;
        String error = "";
        for (int type=0; type<imps.length; type++) {
            if (imps[type] != null) {
                if (isGoodAsInput(imps[type], type, impFieldImps[MAIN])) {
                    impFieldImps[type] = imps[type];
                    changes = true;
                    if (Recorder.record)
                        recordMacro(MACRO_COMMAND_SETINPUT, '"'+IMP_SHORT_NAMES[type]+"\", \""+imps[type].getTitle()+'"');
                } else
                    error += "\n"+imps[type].getTitle()+"\nis not suitable as "+IMP_NAMES[type];
            }
        }
        if (changes) {
            checkForImpChanges();
            for (int type=0; type<N_IMP_FIELDS; type++) {
                ImagePlus imp = impFieldImps[type];
                if (imp == null)                // update what the choices show
                    impFieldChoices[type].select(0);
                else
                    impFieldChoices[type].select(imp.getTitle());
            }
        }
        if (error.length() > 0)
            IJ.error(PLUGIN_NAME, "ERROR"+error);
        IJ.showStatus("");
    }

    /** Asks for the energy (or x axis variable) and sets the stack slice */
    void selectStackSlice() {
        int sliceNumber = stackImp.getCurrentSlice();
        double currentValue = energiesEtc[xAxisVariable][sliceNumber-1];
        String prompt = "Select stack slice for "+ENERGY_ETC_NAMES[xAxisVariable]+":";
        double value = IJ.getNumber(prompt, currentValue);
        if (value == IJ.CANCELED) return;
        sliceNumber = energyToSlice(value, energiesEtc[xAxisVariable], false);
        if (sliceNumber > 0) {
            stackImp.setSlice(sliceNumber);
            if (getEnergyArray() != null && Recorder.record)
                recordMacro(MACRO_COMMAND_SELECTSLICE, LeedUtils.d2s(sliceToEnergy(sliceNumber)));
        } else
            IJ.showStatus("ERROR: No suitable "+ENERGY_ETC_NAMES[xAxisVariable]+" found.");
    }

    /** Creates a mask, with user interaction */
    void createMask() {
        LeedMaskCreator leedMaskCreator = new LeedMaskCreator(this, impFieldImps[MAIN], stackImp);
        (new Thread(leedMaskCreator, "Create LEED Mask")).start();
    }

    /** Reads the parameters of a previous session from a _log.txt file and applies them */
    void readParametersFromFile() {
        Object dataOrError = LeedParams.readFromFile(null, false);
        if (dataOrError == null) return;
        if (dataOrError instanceof String) {
            IJ.error(PLUGIN_NAME, (String)dataOrError);
            return;
        }
        setPatternFile(LeedParams.getString(LeedParams.PATTERNFILE));
        if (allFlagsSet(status, MAIN_OK|MASK_OK))
            triggerUpdateStack(false);      //dark & flat processing may have changed
        if (flagSet(status, MAIN_OK)) {     //metadata sources may have changed
            ImageStack darkStack = impFieldImps[DARK] == null ? null : impFieldImps[DARK].getStack();
            LeedEnergyI0Selector.updateFromParams(this, getStack(), darkStack);
            updateEnergiesEtc(false, false);
        }
        changeStatus(0, INDICES_OK|TRACK_OK|SAVE_OK);
        if (dataOrError != null) {
            try {
                double[] screenFitterArray = (double[])dataOrError;
                LeedScreenFitter screenFitter = new LeedScreenFitter(screenFitterArray); //exception on failure (array too short)
                this.screenFitter = screenFitter;
                this.screenFitterArray = screenFitterArray;
            } catch (Exception e) {
                IJ.error(e.getMessage());
            }
        }
    }

    /** Compares the numeric parameters of a previous session from a _log.txt file with the current ones */
    void compareParametersWithFile() {
        Object dataOrError = LeedParams.readFromFile(null, /*compare=*/true);
        if (dataOrError == null) return;
        if (dataOrError instanceof String)
            IJ.error(PLUGIN_NAME, (String)dataOrError);
    }

    /** Resets all parameters and exits. */
    void resetAllParams() {
        GenericDialog gd = new GenericDialog(PLUGIN_NAME);
        gd.addMessage("Close this Spot Tracker session, discard everything\nand reset all parameters to defaults?");
        gd.enableYesNoCancel("Cancel", "Discard & Reset"); //default button (usually OK') is 'cancel'
        Button[] buttons = gd.getButtons();
        buttons[1].setVisible(false);      //hide 'Cancel' button, we have the OK for it.
        gd.showDialog();
        if (gd.wasCanceled() || gd.wasOKed()) return;
        LeedParams.reset();
        this.gd.dispose();
        cleanup();
    }

// ----------------  U P D A T E   G U I   &   D A T A  ----------------

    /** Called upon new input stack and after user changes the energy I0 etc in LeedEnergyI0Selector.
     *  With checkProcessedI0, checks processed I0 for negative values and shows the corresponding plot.
     *  With showPlot, plots I0, I00 (if present) and the processed I0. */
    void updateEnergiesEtc(boolean checkProcessedI0, boolean showPlot) {
        int energyAxis = energiesEtc[ENERGY]==null && energiesEtc[E_LEEM]!=null && LeedParams.get(LeedParams.LEEMMODE)!=0 ?
            E_LEEM : ENERGY;
        if (energyAxis == E_LEEM && xAxisVariable == ENERGY)
            xAxisVariable = E_LEEM;
        else if (energyAxis == ENERGY && xAxisVariable == E_LEEM)
            xAxisVariable = ENERGY;
        if (energiesEtc[energyAxis] != null) {      //get energyMin, energyMax for the radiusSelector and check them:
            double[] minMax = Tools.getMinMax(energiesEtc[energyAxis]);
            energyMin = minMax[0]; energyMax = minMax[1];
            if (xAxisVariable == ENERGY &&
                    (!LeedEnergyI0Selector.energyOK(energyMin) || !LeedEnergyI0Selector.energyOK(energyMax)))
                energiesEtc[ENERGY] = null;         //not evenly spaced energies in LEED I(V) mode
        }
        if (energiesEtc[ENERGY]==null) {            //in LEEM mode or if we have no energy (e.g. time-dependent data):
            energyMin = LeedRadiusSelector.UNKNOWN_ENERGY;
            energyMax = Double.NaN;
        }
        double[] processedI0 = checkProcessedI0 ? LeedEnergyI0Selector.getProcessedI0(energiesEtc) : null;
        int nI0nonPositive = processedI0 == null ? 0 : LeedUtils.countNonPositives(processedI0);
        if (xAxisVariable>=0 && xAxisVariable < energiesEtc.length && energiesEtc[xAxisVariable] != null && nI0nonPositive == 0)
            changeStatus(ENERGY_OK, 0);
        else
            changeStatus(0, ENERGY_OK | INDICES_OK | TRACK_OK | SAVE_OK);
        setButtonLabelText(E_BUTTON, LeedEnergyI0Selector.getStatusText(energiesEtc, xAxisVariable));
        drawOverlays();
        showRadii();
        if (showPlot || nI0nonPositive > 0) {
            String plotTitle = nI0nonPositive > 0 ?
                    "ERROR: "+nI0nonPositive+" Non-Positive (processed) I0 Values" :
                    impFieldImps[MAIN].getTitle()+" I0";
            Plot plot = LeedEnergyI0Selector.getProcessedI0Plot(plotTitle, energiesEtc, xAxisVariable);
            if (plot != null) {
                plot.show();
                LeedIVAnalyzer.addToPlotImpList(plot.getImagePlus(), 'i');
            }
        }
    }

    /** Reads the spot pattern file and sets status if successful */
    boolean setPatternFile(String patternFilePath) {
        //IJ.log("setPatternFile:"+patternFilePath);
        if (!LeedUtils.fileOk(patternFilePath))
            return false;
        LeedSpotPattern spotPattern = new LeedSpotPattern(patternFilePath);
        if (spotPattern.size() < 3) {
            IJ.error(PLUGIN_NAME,"LEED pattern in "+patternFilePath+"\ninvalid or < 3 spots\n"+spotPattern.getErrorMessage());
            return false;
        }
        this.spotPattern = spotPattern;
        updateSpotPattern();
        return true;
    }

    /** After setting or changing the spot pattern, sets the status and text next to the button */
    void updateSpotPattern() {
        changeStatus(PATTERN_OK, INDICES_OK|TRACK_OK|SAVE_OK);
        setButtonLabelText(PATTERN_BUTTON, spotPattern.getFileName());
        showRadii();
    }

    /** Sets the conversion from k space to screen coordinates.
     *  Called by LeedIndexSelector when the "Set spot indices" is finished
     *  and a valid LeedScreenFitter was obtained (must be non-null). */
    void setScreenFitter(LeedScreenFitter screenFitter, int xOfIndexInput) {
        this.screenFitter = screenFitter;
        this.xOfIndexInput = xOfIndexInput;
        this.screenFitterArray = screenFitter.toArray();
        setButtonLabelText(INDICES_BUTTON,screenFitter.getStatusText()+" spots");
        changeStatus(INDICES_OK, TRACK_OK | SAVE_OK);
        if (Recorder.record) {
            recordMacro(MACRO_COMMAND_SELECTSLICE, LeedUtils.d2s(sliceToEnergy(xOfIndexInput+1)));
            recordMacro(MACRO_COMMAND_SETINDICES, "newArray("+LeedUtils.toString(this.screenFitterArray)+")");
        }
    }

    /** Called after successfully saving the data */
    void setDataSaved(String directory, String prefix) {
        changeStatus(SAVE_OK, 0);
        if (Recorder.record)
            recordMacro(MACRO_COMMAND_SAVEDATA, "\""+directory+"\", \""+prefix+"\"");
    }

    //Shows the integration radii next to the Radius button
    void showRadii() {
        //try{throw new RuntimeException("");}catch(Exception e){IJ.log("showRadii spot="+spotPattern+" energ="+energiesEtc[ENERGY]+" from\n"+e.getStackTrace()[1]+'\n'+e.getStackTrace()[2]);}
        String str = (spotPattern == null || !flagSet(status, ENERGY_OK)) ?
                "" : LeedRadiusSelector.getStatusText(energyMin, energyMax, spotPattern);
        setButtonLabelText(RADIUS_BUTTON, str);
    }

    /** Called when we have IV curves */
    void setIVAnalyzer(LeedIVAnalyzer ivAnalyzer, int nGoodSpots) {
        this.ivAnalyzer = ivAnalyzer;
        setButtonLabelText(TRACK_BUTTON, ivAnalyzer.getStatusText(/*verbose=*/false));
        setButtonLabelText(SAVE_BUTTON, "");
        if (nGoodSpots > 0) {
            changeStatus(TRACK_OK, SAVE_OK);
            ResultsTable badnessTable = ivAnalyzer.makeSpotTable(false);
            if (badnessTable != null)
                badnessTable.show(BADNESS_WARNING_TITLE);
            energiesEtc[BACKGROUND] = ivAnalyzer.getBackgroundIntensities();
        }
        if (energiesEtc[BACKGROUND] != null && energiesEtc[CURRENTI0] != null) {    //calculate background/I0 (in case we want to plot it)
            int n = energiesEtc[BACKGROUND].length;
            energiesEtc[BG_OVER_I0] = new double[n];
            for (int i=0; i<n; i++) {
                double i0 = energiesEtc[CURRENTI0][i];
                if (energiesEtc[CURRENTI00] != null) i0 -= energiesEtc[CURRENTI00][i];
                energiesEtc[BG_OVER_I0][i] = energiesEtc[BACKGROUND][i]/i0;
            }
        }
        makeLogStrings();
        postTrackingChanges = null;
        indexInputImp = LeedUtils.duplicateSlice(stackImp, xOfIndexInput+1);
        if (Recorder.record)
            recordMacro(MACRO_COMMAND_TRACKSPOTS, "");
    }

    String getFlatProcessingString() {
        LeedDarkFlatVirtualStack stack = getStack();
        if (stack == null)
            return "";
        else
            return stack.getFlatProcessingString();
    }

    /** Creates the info written into the log */
    void makeLogStrings() {
        logLines = new ArrayList<String>(50);
        logLines.add(LEED_Spot_Tracker.PLUGIN_NAME+" LEED I(V) Analysis ");
        logLines.add(LeedUtils.getDateFormatted("yyyy-MM-dd HH:mm"));
        logLines.add("");
        logLines.add("Input Images & Image Stacks:");
        for (int i=0; i<N_IMP_FIELDS; i++)
            if (impFieldImps[i] != null)
                logLines.add(IMP_NAMES[i] + ": " + impFieldImps[i].getTitle());
        String info = impFieldImps[MAIN] != null ? impFieldImps[MAIN].getInfoProperty() : null;
        if (info != null) {
            logLines.add("Input information:");
            logLines.add(info);
        }
        if (impFieldImps[FLAT] != null) {
            logLines.add("");
            logLines.add("Dark&flat processing: "+getFlatProcessingString());
        }
        logLines.add("");
        logLines.add("Pattern file: "+LeedParams.getString(LeedParams.PATTERNFILE));
        logLines.add("");
        logLines.add("Energy range, I0 etc.: " + LeedEnergyI0Selector.getStatusText(energiesEtc, xAxisVariable));
        logLines.add("");
        logLines.add("Integration Radius: "+LeedRadiusSelector.getStatusText(energyMin, energyMax, spotPattern));
        logLines.add("");
        logLines.add("Indices defined at " + (xAxisVariable == ENERGY ?
                LeedUtils.d2s(sliceToEnergy(energiesEtc[ENERGY], xOfIndexInput+1))+" eV" :
                "#"+xOfIndexInput+" ("+ENERGY_ETC_NAMES[xAxisVariable]+"="+LeedUtils.d2s(sliceToEnergy(energiesEtc[xAxisVariable], xOfIndexInput+1))
                ));
        logLines.add(screenFitter.getStatusText()+" spots");
        logLines.add("");
        logLines.add(ivAnalyzer.getStatusText(/*verbose=*/true));
        logLines.add("");
        logParamLines = new ArrayList<String>(10);
        logParamLines.add("Parameters (in machine-readable form):");
        for (String str : LeedParams.getParameterLines())
            logParamLines.add(str);
        logParamLines.add(SCREENFIT_KEY+"="+LeedUtils.toString(screenFitter.toArray()));
    }

    public void addToLog(String str) {
        logLines.add(str);
    }

    public ArrayList<String> getLogLines() {
        ArrayList<String> allLines = new ArrayList<String>();
        allLines.addAll(logLines);
        allLines.addAll(logParamLines);
        return allLines;
    }

    /** closes the "WARNING: Beam(s) uncertain (highlighted)" table */
    void closeBadnessWarning() {
        Window win = WindowManager.getWindow(BADNESS_WARNING_TITLE);
        if (win != null) {
            win.dispose();
            WindowManager.removeWindow(win);
        }
    }

    /** Called by the IV Curve Editor to highlight the spots of the group selected there.
     *  May be called with argument null to de-highlight the spots */
    public void highlightSpotWithGroup(String spotName) {
        ImagePlus stackImp = this.stackImp;
        LeedSpotPattern spotPattern = this.spotPattern;
        if (stackImp == null || spotPattern == null || !allFlagsSet(status, TRACK_OK))
            return;
        HashSet<String> spotNames = null;
        if (spotName != null) {
            int iSpot = spotPattern.getIndex(spotName);
            if (iSpot < 0)
                return;
            int group = spotPattern.getGroup(iSpot);
            int[] spots = spotPattern.getAllSpotsForGroup(group);
            spotNames = new HashSet<String>();
            for (int is : spots)
                spotNames.add(spotPattern.getNameForOvly(is));
        }
        LeedOverlay.highlightCircles(stackImp, spotNames, /*strength=*/0);
        stackImp.draw();
    }

// ----------------  S T A T U S   C H A N G E  ----------------

    /** Sets the status and does the next processing step if possible.
     *  'plusFlag' is the new status, 'minusFlag' is a bitwise or
     *  of status flags that should be reset.
     *  Notifies the user by enabling/disabling buttons and painting
     *  the required button red, unless the gui is disabled. */
    void changeStatus(int plusFlag, int minusFlag) {
        int oldStatus = status;
        status &= ~minusFlag;
        status |= plusFlag;
        if (flagSet(minusFlag, MAIN_OK|MASK_OK|ENERGY_OK|PATTERN_OK|INDICES_OK|TRACK_OK|SAVE_OK))
            setButtonLabelText(SAVE_BUTTON, "");
        if (flagSet(minusFlag, MAIN_OK|MASK_OK|ENERGY_OK|PATTERN_OK|INDICES_OK|TRACK_OK))
            setButtonLabelText(TRACK_BUTTON, "");
        if (flagSet(minusFlag, MAIN_OK|MASK_OK|ENERGY_OK|PATTERN_OK|INDICES_OK))
            setButtonLabelText(INDICES_BUTTON, "");
        if (flagSet(plusFlag, MAIN_OK|MASK_OK))
            createOrUpdateStack();
        if (IJ.debugMode) IJ.log("changeStatus +"+plusFlag+" -"+minusFlag+" status="+status+" inputOk="+flagSet(status, MAIN_OK));
        if (componentsEnabled)
            enableAndHighlightComponents(true);
    }

    /** Sets the enable & highlight status of dialog components (buttons, choices, etc.)
     *  according to the current status.
     *  Disables all components if 'enable' = false; */
    void enableAndHighlightComponents(boolean enable) {
        componentsEnabled = enable;
        boolean[] enableButtons = new boolean[N_BUTTONS];
        Arrays.fill(enableButtons, enable);
        // enable/disable
        LeedDarkFlatVirtualStack stack = getStack();
        //IJ.log("enableAndHighlightComponents status="+status+" enable="+enable+" stack="+stack);
        enableButtons[FLAT_BUTTON] = enable && allFlagsSet(status, MAIN_OK|MASK_OK) && stack!= null && stack.hasProcessingOptions();
        enableButtons[E_BUTTON] = enable && allFlagsSet(status, MAIN_OK);
        enableButtons[INDICES_BUTTON] = enable && allFlagsSet(status, MAIN_OK|MASK_OK|PATTERN_OK);
        enableButtons[RADIUS_BUTTON] = enable && allFlagsSet(status, MAIN_OK|PATTERN_OK|ENERGY_OK);
        enableButtons[TRACK_BUTTON] = enable && allFlagsSet(status, MAIN_OK|MASK_OK|PATTERN_OK|ENERGY_OK|INDICES_OK);
        enableButtons[SAVE_BUTTON] = enable && allFlagsSet(status, MAIN_OK|MASK_OK|PATTERN_OK|ENERGY_OK|INDICES_OK|TRACK_OK);
        for (int i=0; i<N_BUTTONS; i++)
            buttons[i].setEnabled(enableButtons[i]);
        boolean unsavedData = flagSet(status, TRACK_OK) && !flagSet(status, SAVE_OK);
        for (Choice c : impFieldChoices)
            c.setEnabled(enable);
        //red text for the next input required
        boolean trackAgain = postTrackingChanges != null && flagSet(status, TRACK_OK) && !flagSet(status, SAVE_OK);
        highlightImpField(MAIN, !flagSet(status, MAIN_OK));
        highlightImpField(MASK, flagSet(status, MAIN_OK) && !flagSet(status, MASK_OK));
        highlightImpField(FLAT, NO_INPUT_STRING.equals(impFieldChoices[FLAT].getSelectedItem()) &&
                !NO_INPUT_STRING.equals(impFieldChoices[DARK2].getSelectedItem())); //no flat but dark2
        highlightComponent(buttons[E_BUTTON], flagSet(status, MAIN_OK) && !flagSet(status, ENERGY_OK));
        highlightComponent(buttons[PATTERN_BUTTON], !flagSet(status, PATTERN_OK));
        highlightComponent(buttons[INDICES_BUTTON], flagSet(status, ENERGY_OK) && !flagSet(status, INDICES_OK));
        highlightComponent(buttons[TRACK_BUTTON], !flagSet(status, TRACK_OK) || trackAgain);
        highlightComponent(buttons[SAVE_BUTTON], !flagSet(status, SAVE_OK) && !trackAgain);
        highlightComponent(buttons[MORE_BUTTON], !flagSet(status, MAIN_OK) || !flagSet(status, MASK_OK));

        boolean allDone = allFlagsSet(status, MAIN_OK|MASK_OK|PATTERN_OK|ENERGY_OK|INDICES_OK|TRACK_OK|SAVE_OK);
        okButton.setEnabled(enable && allDone);
        cancelButton.setEnabled(enable);
        cancelButton.setVisible(!allDone);
        currentFrontDialog = null;
    }

    /** Sets the fields (Label and Choice) for an input image/stack red or black.
     *  On MacOS, changing the color of a Choice will be ignored. */
    void highlightImpField(int fieldIndex, boolean makeRed) {
        highlightComponent(impFieldChoices[fieldIndex], makeRed);
        highlightComponent(impFieldLabels[fieldIndex], makeRed);
    }

    void highlightComponent(Component c, boolean makeRed) {
        if (!c.isEnabled()) makeRed = false;
        c.setForeground(makeRed ? Color.RED : Color.BLACK);
    }

    static boolean flagSet(int number, int flag) {
        return (number&flag) != 0;
    }

    /** With a bitwise or of the flags, returns whether all flags are set in 'number'*/
    static boolean allFlagsSet(int number, int flags) {
        return (number&flags) == flags;
    }

    /** Defines a dialog that should be brought to the front when the main dialog is brought to the front.
     *  Also avoids closing the main dialog while that other dialog is open. */
    void setCurrentFrontDialog(Dialog d) {
        this.currentFrontDialog = d;
    }

    /** Checks whether new images have been selected and reacts accordingly */
    void checkForImpChanges() {
        boolean updateStack = false;
        for (int i=0; i<N_IMP_FIELDS; i++) {
            if (oldImpFieldImps[i] != impFieldImps[i]) {
                //IJ.log("checkForImpChanges: changed= "+IMP_NAMES[i]);
                if (impFieldImps[i] != null)
                    impPreviousSlice[i] = impFieldImps[i].getCurrentSlice();
                if (i == MAIN) {
                    if (impFieldImps[MAIN] == null)
                        changeStatus(0, MAIN_OK | MASK_OK | INDICES_OK | TRACK_OK | SAVE_OK);
                    else
                        changeStatus(MAIN_OK, INDICES_OK | TRACK_OK | SAVE_OK);
                } else if (i == MASK) {
                    if (impFieldImps[MASK] == null)
                        changeStatus(0, MASK_OK | INDICES_OK | TRACK_OK | SAVE_OK);
                    else {
                        analyzeMask(); //changes the status, this will update the stack
                        updateStack = false;
                    }
                } else { //dark, flat: updates stack
                    updateStack = true;
                    if (i == FLAT && impFieldImps[FLAT] != null &&
                            impFieldImps[FLAT].getProp(LeedDarkFlatVirtualStack.PROCESSED_FLAT_INFOKEY) != null) {
                        impFieldChoices[DARK2].select(0);       //we got a processed flat field, it should have no dark frame
                        impFieldImps[DARK2] = null;
                        LeedParams.set(LeedParams.FLATFITORDER, Math.min(LeedParams.get(LeedParams.FLATFITORDER), 0)); //no polynomial fit
                    }
                }
                updateImpHashEtc(i);
            }
        }
        //IJ.log("checkForImpChanges: updateStack= "+updateStack);
        if (updateStack) {
            createOrUpdateStack();
            changeStatus(0, TRACK_OK | SAVE_OK);
        }
        oldImpFieldImps = (ImagePlus[])impFieldImps.clone();
    }

    /** Shows the popup menu of the More>>button */
    void showPopupMenu(ActionEvent e) {
        popupMenu.removeAll();
        for (int i=0; i<MORE_MENU_TEXTS.length; i++) {
            MenuItem mi = new MenuItem(MORE_MENU_TEXTS[i]);
            boolean enable = allFlagsSet(status, POPUP_REQUIRED_STATUS[i]);
            if ((i == DELETE_HIGHLIGHTED || i == UNHIGHLIGHT_SPOTS) && LeedOverlay.getHighlighted(stackImp) == null)
                enable = false;
            if (i == CLOSE_PLOTS)
                enable = LeedIVAnalyzer.anyOpenPlots((char)0);
            mi.setEnabled(enable);
            mi.addActionListener(this);
            popupMenu.add(mi);
        }
        Component src = (Component)e.getSource();
        popupMenu.show(src, src.getWidth()/2, src.getHeight()/2);
    }

    /** Called if the user changes parameters: Energies, I0... or radii
     *  Changes between spot tracking & saving lead to a question whether to really save  */
    void setPostTrackingChanges(String what) {
        if (!flagSet(status, TRACK_OK)) return;
        if (postTrackingChanges==null)
            postTrackingChanges = what;
        else if (!postTrackingChanges.contains(what))
            postTrackingChanges += ", "+what;
        enableAndHighlightComponents(true); //update dialog: 'Track', not 'Save' should be red.
    }
// ----------------  C A L L B A C K    ( O F   L I S T E N E R S ) ----------------

    /** This callback method is called when the user changes choices (for the input files). */
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        try {
            if (e == null) return true;         //e=null on final call, when dialog is closed
            for (int i=0; i<N_IMP_FIELDS; i++)
                if (e.getSource() == impFieldChoices[i]) {
                    String title = impFieldChoices[i].getSelectedItem();
                    impFieldImps[i] = WindowManager.getImage(title);
                    if (i == MAIN)
                        refreshImageChoices(null, false);
                    checkForImpChanges();
                    if (Recorder.record)
                        recordMacro(MACRO_COMMAND_SETINPUT, '"'+IMP_SHORT_NAMES[i]+"\", \""+title+'"');
                }
        } catch (Exception ex) { IJ.handleException(ex); return false; }
        return okButton.isEnabled();            //return value is used for okButton.setEnabled; we have done this already
    }

    /** Callback routine handling buttons */
    public void actionPerformed(ActionEvent e) {
        try {
            if (e.getSource() instanceof Button) {
                int srcIndex = LeedUtils.arrayIndexOf(buttons, e.getSource());
                switch (srcIndex) {
                    case FLAT_BUTTON:
                        askForFlatProcessing(); return;
                    case E_BUTTON:
                        askForEnergies(); return;
                    case INDICES_BUTTON:
                        askForIndices(); return;
                    case PATTERN_BUTTON:
                        askForPatternFile(); return;
                    case RADIUS_BUTTON:
                        askForRadius(); return;
                    case TRACK_BUTTON:
                        analyzeIV(true); return;
                    case SAVE_BUTTON:
                        saveData(null, null, false); return;
                    case MORE_BUTTON:
                        showPopupMenu(e); return;
                    default:
                }
            } else {    // MenuItems of PopupMenu
                int srcIndex = LeedUtils.arrayIndexOf(MORE_MENU_TEXTS, e.getActionCommand());
                switch (srcIndex) {
                    case OPEN_FILE:
                        openFiles(); return;
                    case SELECT_STACK_SLICE:
                        selectStackSlice(); return;
                    case CREATE_MASK:
                        createMask(); return;
                    case HIGHLIGHT_SPOTS:
                        LeedHighlightDialogs.highlightSpots(stackImp, spotPattern); return;
                    case DELETE_HIGHLIGHTED:
                        String logStr = LeedHighlightDialogs.deleteHighlightedSpots(this, stackImp, spotPattern, ivAnalyzer);
                        if (logStr != null) {
                            addToLog(logStr);
                            changeStatus(0, SAVE_OK);
                        }
                        return;
                    case UNHIGHLIGHT_SPOTS:
                        LeedOverlay.highlightCircles(stackImp, null, 1); stackImp.draw(); return;
                    case SPOT_STATISTICS:
                        ivAnalyzer.makeSpotTable(true).show("Beam statistics"); return;
                    case PLOT_BEAM:
                        LeedBeamPlotter.doDialogAndPlot(spotPattern, ivAnalyzer, screenFitter, energiesEtc); return;
                    case CLOSE_PLOTS:
                        LeedCloseDialog.showDialogAndClose(this, null); return;
                    case UNDISTORT:
                        ImageStack stack = getStack();
                        if (stack == null || stack.size() == 0) return;
                        int stackSlice = stackImp.getCurrentSlice();
                        (new LeedUndistorter()).undistort(LeedUtils.removeExtension(impFieldImps[MAIN].getTitle()),
                                this, stack, maskRoi, stackSlice, energiesEtc, screenFitter, ivAnalyzer);
                        return;
                    case LIST_PARAMETERS:
                        LeedParams.makeTable(PLUGIN_NAME+" Parameters"); return;
                    case READ_PARAMETERS:
                        readParametersFromFile(); return;
                    case COMPARE_PARAMETERS:
                        compareParametersWithFile(); return;
                    case DATA_KEYWORDS:
                        LeedParams.showParamsDialog(); return;
                    case RESET_ALL_PARAMS:
                        resetAllParams(); return;
                    default:
                }
            }
            LeedUtils.logError("Internal Error: Event handling missing for "+e.getActionCommand());
        } catch (Exception ex) { IJ.handleException(ex); }
    }

    /** When an image gets opened, the choice menus for the images have to be updated */
    public void imageOpened(ImagePlus imp) {
        refreshImageChoices(null, false);
    }

    /** Called when an image got closed; we have to react if this is an input.
     *  Also the choice menus for the inputs must be updated.
     *  If it is the I(V) measurement stack window which should not be closed, 
     *  recreates the window */
    public void imageClosed(ImagePlus imp) {
        for (ImagePlus imp2 : impFieldImps)
            if (imp == imp2)
                refreshImageChoices(imp, false);
        if (imp == stackImp && gd != null && gd.isVisible() && flagSet(status, MAIN_OK)) {
            createNewStack();
            drawOverlays();
            gd.toFront();
        }
    }

    /** Changes of image data (or stack slice, which is not relevant and would just create
     *  a lot of unnecessary activity) */
    public void imageUpdated(ImagePlus imp) {
        try {
            if (imp == stackImp) {
                int currentSlice = stackImp.getCurrentSlice();
                if (currentSlice != stackImpPreviousSlice) {
                    updateContrast(stackImp, impFieldImps[MASK]);
                    stackImpPreviousSlice = currentSlice;
                    updateCurveEditors();
                }
                return;
            }
            for (int i=0; i<N_IMP_FIELDS; i++) {
                ImagePlus fieldImp = impFieldImps[i];
                if (fieldImp == imp) {
                    if (imp.getImageStackSize() != impPreviousStackSize[i]) {   //stack size changed
                        if (isGoodAsInput(imp, i, impFieldImps[MAIN])) {
                            oldImpFieldImps[i] = null;      //changing slices is like selecting a different input
                        } else {
                            impFieldImps[i] = null;         //not suitable as this input type any more
                            anyImageChanged = true;         //we will update Choices for input images when the dialog comes to front
                        }
                        checkForImpChanges();
                        return;
                    }
                    
                    if (imp.getCurrentSlice() == impPreviousSlice[i]) { //slice unchanged, has contents changed?
                        int hashCode = LeedUtils.getHashCode(imp);
                        if (hashCode == impPreviousHashCode[i])
                            return;
                        impPreviousHashCode[i] = hashCode;
                        if (imp == impFieldImps[MASK])
                            analyzeMask();
                        else {
                            LeedDarkFlatVirtualStack stack = getStack();
                            if (stack != null)
                                stack.update(stackImp);
                        }
                    } else {                                //slice changed, probably no contents changed
                        impPreviousSlice[i] = imp.getCurrentSlice();
                        updateImpHashEtc(i);
                    }
                    return;
                }
            }
            anyImageChanged = true;        //other image updated, we will update Choices for input images when the dialog comes to front
        } catch (Exception e) {IJ.handleException(e);}
    }

    /** Updates the saved stack size, stack slice, and pixel data hashcode */
    private void updateImpHashEtc(int impType) {
        ImagePlus imp = impFieldImps[impType];
        if (imp == null) {
            impPreviousSlice[impType] = -1;
            return;
        }
        impPreviousStackSize[impType] = imp.getImageStackSize();
        impPreviousSlice[impType] = imp.getCurrentSlice();
        impPreviousHashCode[impType] = LeedUtils.getHashCode(imp);
    }

    /** Synchronized the vertical line of the Curve Editors with the current energy */
    void updateCurveEditors() {
        double energy = getEnergy();
        if (Double.isNaN(energy)) return;
        ArrayList<LeedCurvePlotEditor> curveEditors = LeedCurvePlotEditor.getAllCurvePlotEditors();
        //DEBUG String str="";for (LeedCurvePlotEditor ce:curveEditors) str+="\n"+ce.pathPrefix;IJ.log("E="+energy+" curve editors: "+str);
        if (curveEditors != null)
            for (LeedCurvePlotEditor curveEditor : curveEditors)
                curveEditor.setSpotTrackerEnergy(energy);
    }

    /** Returns the current energy, or NaN if unknown */
    public double getEnergy() {
        int sliceNumber = stackImp == null ? 0 : stackImp.getCurrentSlice();
        double[] energies = energiesEtc[ENERGY] != null ? energiesEtc[ENERGY] : energiesEtc[E_LEEM];
        if (energies == null || sliceNumber <= 0 || sliceNumber > energies.length)
            return Double.NaN;
        return energies[sliceNumber-1];
    }

    /** Puts the dialog to the front */
    public void toFront() {
        gd.setVisible(true);
        gd.toFront();
    }

    /** WindowListener callback, called when the window becomes the
     *  active window. When the main SpotTracker dialog becomes active,
     *  and a daughter dialog is open, puts that daugher dialog to the front.
     *  Otherwise, when the main SpotTracker dialog becomes active and
     *  images have changed in the meanwhile, updates the choices for
     *  the input images. */
    public void windowActivated(WindowEvent e) {
        Window win = e.getWindow();
        if (win == gd) {
            if (anyImageChanged)
                refreshImageChoices(null, false);
            if (currentFrontDialog != null) {
                currentFrontDialog.toFront();
                return;
            }
        }
        if (win instanceof WindowListener) //The WindowListener of ImageJ needs to know 
            ((WindowListener)win).windowActivated(e);
    }

    /** WindowListener callback, called when the user clicks on
     *  the close box of the spot tracking stack or the dialog.
     *  [We cannot avoid closing the spot tracking stack via the
     *  ImageJ menu or shortcut; in that case it is re-created]
     *  Also called when the user tries to close the main dialog
     *  while a non-modal sub-dialog is active. */
    public void windowClosing(WindowEvent e) {
        Window win = e.getWindow();
        if (win == getStackWindow() && flagSet(status, MAIN_OK) && gd.isVisible()) {
            IJ.beep();
            gd.toFront();
            IJ.showStatus("Stack is protected; close the "+PLUGIN_NAME+" dialog");
        } else if (win == gd && currentFrontDialog != null && currentFrontDialog.isVisible()) {
            IJ.beep();
            currentFrontDialog.toFront();
            IJ.showStatus("Close '"+currentFrontDialog.getTitle()+"' first");
        } else if (win instanceof WindowListener)
            ((WindowListener)win).windowClosing(e); //ImageJ cares about closing
    }

    //other WindowListener callbacks
    public void windowClosed(WindowEvent e) {}
    public void windowOpened(WindowEvent e) {}
    public void windowIconified(WindowEvent e) {}
    public void windowDeiconified(WindowEvent e) {}
    public void windowDeactivated(WindowEvent e) {}

    /** Intercept 'esc' for closing and show the spot tracking stack on <return> */
    public void keyPressed(KeyEvent e) {
        int keyCode = e.getKeyCode();
        if (keyCode == KeyEvent.VK_ENTER) {     //Shift-Enter takes spot tracking window to the top, enter ImageJ
            ImagePlus stackImp = this.stackImp;
            Window win = stackImp == null || (e.getModifiers() & KeyEvent.SHIFT_MASK) == 0 ?
                    IJ.getInstance() : stackImp.getWindow();
            if (win != null) win.toFront();
            e.setKeyCode(0);                    //spoof the key code: don't close the GenericDialog
        } else if (keyCode == KeyEvent.VK_ESCAPE) {
            if (flagSet(status, MAIN_OK))       //ESC closes only if blank
                e.setKeyCode(0);                //spoof the key code: don't abort the GenericDialog
            else
                gd.dispose();
        }
    }

    //other KeyListner callbacks
    public void keyReleased(KeyEvent e) {}
    public void keyTyped(KeyEvent e) {}

// ----------------  M A C R O   I N T E R F A C E  ----------------

    public static void startMacro() {
        if (instance == null) {
            Thread thread = new Thread(new Runnable() {public void run() {
                        new LEED_Spot_Tracker().run("");
                    }}, PLUGIN_NAME);
            thread.start();
            long time = System.currentTimeMillis();
            while (instance == null) {
                IJ.wait(20);
                if (System.currentTimeMillis() - time > 2000)
                    IJ.error("Macro Timeout - Could not start "+PLUGIN_NAME);
            }
            while (instance.gd==null || !instance.gd.isVisible()) {
                IJ.wait(100);
                if (System.currentTimeMillis() - time > 50000) //50 s timeout; initial processing of stacks can take some time
                    IJ.error("Macro Timeout waiting for "+PLUGIN_NAME+" dialog window");
            }
        }
        Functions.registerExtensions(instance);             // register macro extensions
    }

    /** ImageJ Macro Commands that can be called with 'Ext.<command>(args) must be defined here.
     *  The return value is always a String, so it has to be converted to a number with parseFloat or parseInt if required!
     *  Note that you have to start with a macro command 'call("LEED_Spot_Tracker.startmacro"); */
    private ExtensionDescriptor[] extensions = {
		ExtensionDescriptor.newDescriptor(MACRO_COMMAND_SETINPUT, this, ARG_STRING, ARG_STRING),
		ExtensionDescriptor.newDescriptor(MACRO_COMMAND_UPDATE, this, ARG_STRING),
		ExtensionDescriptor.newDescriptor(MACRO_COMMAND_SETVALUE, this, ARG_STRING, ARG_NUMBER),
		ExtensionDescriptor.newDescriptor(MACRO_COMMAND_SETARRAY, this, ARG_STRING, ARG_ARRAY),
        ExtensionDescriptor.newDescriptor(MACRO_COMMAND_SETPATTERN, this, ARG_STRING),
        ExtensionDescriptor.newDescriptor(MACRO_COMMAND_SELECTSLICE, this, ARG_NUMBER),
        ExtensionDescriptor.newDescriptor(MACRO_COMMAND_SETINDICES, this, ARG_ARRAY),
        ExtensionDescriptor.newDescriptor(MACRO_COMMAND_TRACKSPOTS, this),
        ExtensionDescriptor.newDescriptor(MACRO_COMMAND_SAVEDATA, this, ARG_STRING, ARG_STRING),
        ExtensionDescriptor.newDescriptor(MACRO_COMMAND_CLOSE, this, ARG_STRING),
        ExtensionDescriptor.newDescriptor(MACRO_COMMAND_MAKETABLE, this, ARG_STRING),
        ExtensionDescriptor.newDescriptor(MACRO_COMMAND_GETSTATUS, this, ARG_STRING)
	};

    /** This method is called when operating through ImageJ macros
     *  (macro extensions, Ext.* macro calls) */
    public String handleExtension(final String name, final Object[] args) {
        String macroResult = null;
        //DEBUG IJ.log("macro call "+name+" @"+(System.currentTimeMillis()-t0));
        try {
            if (name.equals(MACRO_COMMAND_SETINPUT)) {
                macroResult = setInputImageViaMacro((String)args[0], (String)args[1]);
            } else if (name.equals(MACRO_COMMAND_UPDATE)) {
                String what = ((String)args[0]).toLowerCase();
                if (what.contains("stack")) {
                    if (allFlagsSet(status, MAIN_OK|MASK_OK)) {
                        createOrUpdateStack();
                        if (!waitForStackUpdate()) {
                            Macro.abort();
                            return null;
                        }
                    }
                } else if (what.contains("meta")) {
                    LeedDarkFlatVirtualStack stack = getStack();
                    if (stack == null)
                        macroResult = "ERROR: Must specify input stacks before update metadata";
                    else {
                        ImageStack darkStack = impFieldImps[DARK] == null ? null : impFieldImps[DARK].getStack();
                        macroResult = LeedEnergyI0Selector.updateFromParams(this, stack, darkStack);
                        updateEnergiesEtc(false, false);
                    }
                }
            } else if (name.equals(MACRO_COMMAND_SETVALUE)) {
                macroResult = LeedParams.setValue((String)args[0], ((Double)args[1]).doubleValue());
                setButtonLabelText(FLAT_BUTTON, getFlatProcessingString());
                setButtonLabelText(E_BUTTON, LeedEnergyI0Selector.getStatusText(energiesEtc, xAxisVariable));
                showRadii();
            } else if (name.equals(MACRO_COMMAND_SETARRAY)) {
                macroResult = setArrayFromMacro((String)args[0], (Object[])args[1]);
            } else if (name.equals(MACRO_COMMAND_SETPATTERN)) {
                boolean patternOk = setPatternFile((String)args[0]);
                if (!patternOk) macroResult = "ERROR: No or invalid Spot Pattern File";
            } else if (name.equals(MACRO_COMMAND_SELECTSLICE)) {
                if (!waitForStackUpdate()) {
                    Macro.abort();
                    return null;
                }
                double energy = (Double)args[0];
                int sliceNumber = energyToSlice(energy, getEnergyArray(), true);
                if (stackImp==null || sliceNumber<1 || sliceNumber>stackImp.getImageStackSize())
                    macroResult = "ERROR: Energy "+LeedUtils.d2s(energy)+" not in stack";
                else
                    stackImp.setSlice(sliceNumber);
            } else if (name.equals(MACRO_COMMAND_SETINDICES)) {
                if (!waitForStackUpdate()) {
                    Macro.abort();
                    return null;
                }
                Object[] array = (Object[])args[0];
                double[] screenFitterArray = LeedUtils.doubleFromDouble(array);
                macroResult = setScreenFitterFromArray(screenFitterArray);
            } else if (name.equals(MACRO_COMMAND_TRACKSPOTS)) {
                if (allFlagsSet(status, MAIN_OK|MASK_OK|PATTERN_OK|ENERGY_OK|INDICES_OK)) {
                    if (!waitForStackUpdate()) {
                        Macro.abort();
                        return null;
                    }
                    analyzeIV(false);
                } else
                    macroResult = "ERROR: Cannot track spots yet";
            } else if (name.equals(MACRO_COMMAND_SAVEDATA)) {
                if (allFlagsSet(status, MAIN_OK|MASK_OK|PATTERN_OK|ENERGY_OK|INDICES_OK|TRACK_OK))
                    saveData((String)args[0], (String)args[1], false);
                else
                    macroResult = "ERROR: No data to save yet";
            } else if (name.equals(MACRO_COMMAND_CLOSE)) {
                String arg = (String)args[0];
                Macro.setOptions(arg);      //injects the argument as GenericDialog options for the LeedCloseDialog
                Thread.currentThread().setName("Run$_"+Thread.currentThread().getName()); //makes the GenericDialog accept macro options
                boolean closeTracker = arg.toLowerCase().contains("tracker");
                if (closeTracker) {
                    this.gd.dispose();
                    cleanup();
                }
                LeedCloseDialog.showDialogAndClose(this, closeTracker ? stackImp : null);
                Thread.currentThread().setName(Thread.currentThread().getName().substring("Run$_".length()));
            } else if (name.equals(MACRO_COMMAND_MAKETABLE)) {
                String arg = (String)args[0];
                macroResult = makeTable(arg);
            } else if (name.equals(MACRO_COMMAND_GETSTATUS)) {
                return getStatus((String)args[0]);
            }
        } catch (Exception e) {
            IJ.handleException(e);
            Macro.abort();
        }
        if (macroResult != null && macroResult.startsWith("ERROR")) {
            Interpreter interpreter = Interpreter.getInstance();
            if (interpreter != null)
                macroResult += "\nin macro line number "+interpreter.getLineNumber();
            IJ.error(PLUGIN_NAME, macroResult);
        }
        return macroResult;
    }

    /** Sets an array of energiesEtc from a macro. An array length of 0 sets the array to null */
    String setArrayFromMacro(String name, Object[] array) {
        int inputType = LeedUtils.arrayIndexOf(ENERGY_ETC_NAMES, name);
        if (inputType < 0)
            return "ERROR: Unknown array type '"+name+"'";
        if (array.length == 0) {
            energiesEtc[inputType] = null;
            if (inputType < LeedEnergyI0Selector.N_DIALOG_DATA)
                LeedEnergyI0Selector.setDataSource(inputType, LeedEnergyI0Selector.NONE);
        } else {
            ImageStack stack = getStack();
            if (stack != null && stack.getSize() != array.length)
                return "ERROR: Array length ("+array.length+") does not fit stack size ("+stack.getSize()+")";
            double[] data = LeedUtils.doubleFromDouble(array);
            if (inputType < LeedEnergyI0Selector.N_DIALOG_DATA)
                LeedEnergyI0Selector.setDataSource(inputType, LeedEnergyI0Selector.PREV);
            if (inputType == ENERGY && LeedParams.get(LeedParams.LEEMMODE) != 0 ) {
                inputType = E_LEEM;
                energiesEtc[ENERGY] = null;
            }
            energiesEtc[inputType] = data;
        }
        updateEnergiesEtc(false, false);
        return null;
    }

    /** Sets the input image for the type as defined in IMP_SHORT_NAMES to the given image.
     *  Returns null if ok, error message if failed */
    String setInputImageViaMacro(String inputTypeStr, String inputName) {
        int inputType = LeedUtils.arrayIndexOf(IMP_SHORT_NAMES, inputTypeStr);
        if (inputType < 0)
            return "ERROR: Invalid input type '"+inputTypeStr+"'";
        boolean setToNull = NO_INPUT_STRING.equals(inputName) || inputName.length()==0;
        ImagePlus imp = setToNull ? null : WindowManager.getImage(inputName);
        if (!setToNull) {
            if (imp == null)
                return "ERROR: No open image '"+inputName+"'";
            if (!isGoodAsInput(imp, inputType, impFieldImps[MAIN])) {
                String err = "ERROR: Image '"+inputName+"' is not suitable as "+IMP_NAMES[inputType];
                if (inputType != MAIN)
                    err += "\n or not same size as main input";
                return err;
            }
        }
        impFieldImps[inputType] = imp;
        if (inputType == MAIN)
            refreshImageChoices(null, false);     //when main is selected, also sets non-fitting aux images to null
        impFieldChoices[inputType].select(inputName);   //display update, does not trigger an event
        checkForImpChanges();                           //updates status and triggers updateStack
        return null;
    }

    /** Creates a table of data; for use in macros (see MACRO_COMMAND_MAKETABLE).
     *  'title' must contain one of "meta" (metadata: energies, I0...),
     *  "indices" (array for "+MACRO_COMMAND_SETINDICES+"), "beam" (statistics) */
    String makeTable(String title) {
        ResultsTable rt = null;
        String titleLowerCase = title.toLowerCase();
        if (titleLowerCase.contains("meta")) {
            rt = new ResultsTable();
            for (int c=0; c<energiesEtc.length; c++)
                if (energiesEtc[c] != null)
                    rt.setValues(ENERGY_ETC_NAMES[c], energiesEtc[c]);
        } else if (titleLowerCase.contains("indices")) {
            if (!flagSet(status, INDICES_OK))
                return "ERROR: Indices not set yet yet";
            rt = new ResultsTable();
            rt.setValues("IndicesData",screenFitter.toArray());
        } else if (titleLowerCase.contains("beam")) {
            if (!flagSet(status, TRACK_OK)) return "ERROR: No spot tracking data";
            rt = ivAnalyzer.makeSpotTable(true);
        } else return "ERROR: Underfined table type '"+title+"'";
        rt.show(title);
        return null;
    }

    /** Sets the screenFitter, i.e. the conversion from k space to screen coordinates
     *  Returns null if ok, error message if failed */
    String setScreenFitterFromArray(double[] screenFitterArray) {
        if (!allFlagsSet(status, MAIN_OK|MASK_OK|PATTERN_OK))
            return "ERROR: Cannot set reciprocal-space indices yet";
        int stackSlice = stackImp.getCurrentSlice();
        double energy = sliceToEnergy(energiesEtc[ENERGY], stackSlice);
        int spotBackgrShape = (int)LeedParams.get(LeedParams.BACKGROUNDTYPE);
        double radius = LeedRadiusSelector.radius(energy, false);
        double azBlurRadians = Math.toRadians(LeedParams.get(LeedParams.AZIMUTHBLURANGLE));
        LeedIndexSelector indexSelector = new LeedIndexSelector(this, stackImp, impFieldImps[MASK],
                maskRoi, energiesEtc, xAxisVariable, spotPattern, screenFitterArray);
        double[][] xyMax = LeedIndexSliceSelector.findSpotMaxima(stackImp, impFieldImps[MASK],
            maskRoi, spotBackgrShape, radius, azBlurRadians, LeedParams.get(LeedParams.MINSIGNIFICANCEINDEX));
        if (xyMax[0].length < LeedIndexSliceSelector.MIN_N_MAXIMA)
            return "ERROR: Not enough maxima in image (only "+xyMax[0].length+")";
        LeedIndexSliceSelector.showMaxOverlay(stackImp, stackSlice, xyMax, radius);
        LeedScreenFitter screenFitter = indexSelector.useScreenFitterArray(xyMax[0], xyMax[1], stackSlice, energy,
                /*searchRadius default*/Double.NaN, /*showLabels=*/!Interpreter.isBatchMode());
        if (screenFitter == null)
            return("ERROR: Array in "+MACRO_COMMAND_SETINDICES+" does not fit");
        else
            setScreenFitter(screenFitter, stackSlice-1);
        return null;
    }

    /** Returns the status string on the SpotTracker panel, where 'what' can be:
     *  'iv', 'dark', 'mask', ... for the input image/stack or the following for the text next to button:
     *  'darkflat', 'energies', 'pattern', 'indices', 'radius', 'track', 'save';
     *  or 'all' for all of these */
    String getStatus(String what) {
        what = what.toLowerCase();
        if (what.equals("all")) {
            StringBuilder sb = new StringBuilder(500);
            for (String item : Tools.split("iv,dark,flat,dark2,mask,dark&flat_proc,energies,pattern,indices,radius,track,save",",")) {
                sb.append(item);
                sb.append(": ");
                sb.append(getStatus(item));
                sb.append('\n');
            }
            return sb.toString();                
        } else if (what.contains("dark") && what.contains("flat")) {
            LeedDarkFlatVirtualStack stack = getStack();
            return stack == null ? "" : stack.getFlatProcessingString();
        } else if (what.contains("energ")) {
            if (flagSet(status, ENERGY_OK)) return LeedEnergyI0Selector.getStatusText(energiesEtc, xAxisVariable);
            else return "";
        } else if (what.contains("pattern")) {
            if (spotPattern == null) return "";
            else return spotPattern.getFileName();
        } else if (what.contains("indices")) {
            if (flagSet(status, INDICES_OK)) return screenFitter.getStatusText()+" spots";
            else return "";
        } else if (what.contains("radi")) {
            if (spotPattern == null || !flagSet(status, ENERGY_OK)) return "";
            else return LeedRadiusSelector.getStatusText(energyMin, energyMax, spotPattern);
        } else if (what.contains("track")) {
            if (flagSet(status, TRACK_OK)) return ivAnalyzer.getStatusText(false);
            else return "";
        } else if (what.contains("save")) {
            if (flagSet(status, SAVE_OK)) return ivAnalyzer.getStatusText(false);
            else return "";
        }

        for (int i=N_IMP_FIELDS-1; i>=0; i--) { //is 'what' an image/stack? reverse: first search for flat2, then for flat
            if (what.contains(IMP_SHORT_NAMES[i])) {
                ImagePlus imp = impFieldImps[i];
                return imp==null ? "" : imp.getTitle();
            }
        }
        return "ERROR: Undefined status item: '+what'";
    }

    /** Returns to ImageJ the list of macro commands supported */
    public ExtensionDescriptor[] getExtensionFunctions() {
        return extensions;
    }

    /** Records a macro String with a given name and argument */
    static void recordMacro(String name, String argument) {
        boolean isMacro = Thread.currentThread().getName().endsWith("Macro$")  && !Recorder.recordInMacros;
        if (isMacro) return;
        long millis = System.currentTimeMillis();
        if (millis - lastRecorderHintsShownMillis > 1000*60*5) {
            if (Recorder.scriptMode()) {
                Recorder.recordString("//Error: "+PLUGIN_NAME+" command recording works only in 'Macro' mode\n");
            } else {
                Recorder.recordString("// Start a macro for the "+PLUGIN_NAME+" with:\n"+
                        "     call(\"LEED_Spot_Tracker.startMacro\");\n"+
                        "// Commands:\n"+
                        "//   Ext."+MACRO_COMMAND_SETINPUT+"(type, name);  // selects the named input image/stack ('iv', 'dark', 'mask', ...)\n"+
                        "//   Ext."+MACRO_COMMAND_UPDATE+"(what);          // updates 'stack' (after changing dark&flat processing) or 'meta' (matadata: energy, I0...)\n"+
                        "//   Ext."+MACRO_COMMAND_SETVALUE+"(name, value); // sets a numeric parameter *\n"+
                        "//   Ext."+MACRO_COMMAND_SETARRAY+"(name, array); // sets a metadata array (energy, I0, ...) for the stack slices\n"+
                        "//   Ext."+MACRO_COMMAND_SELECTSLICE+"(value); // sets the current stack slice to the desired energy value\n"+
                        "//   Ext."+MACRO_COMMAND_SETPATTERN+"(pathToPatternfile); // sets spot pattern file\n"+
                        "//   Ext."+MACRO_COMMAND_SETINDICES+"(array);     // sets the spot indices (array with fit data for screen coordinates vs. k)\n"+
                        "//   Ext."+MACRO_COMMAND_TRACKSPOTS+"();          // runs spot tracking\n"+
                        "//   Ext."+MACRO_COMMAND_SAVEDATA+"(dir, prefix); // saves the results (dir must exist)\n"+
                        "//   Ext."+MACRO_COMMAND_CLOSE+"(what);           // closes 'what', e.g. 'tracker stack'\n"+
                        "//   Ext."+MACRO_COMMAND_MAKETABLE+"(title);      // creates a table; 'title' must contain one of \"meta\" (metadata: energies, I0...), \"indices\" (array for "+MACRO_COMMAND_SETINDICES+"), or \"beam\" (statistics)\n"+
                        "//   Ext."+MACRO_COMMAND_GETSTATUS+"(what);       // returns image/stack name ('what' = \"iv\", \"dark\", \"flat\", \"dark2\", \"mask\") or text next to button ('what' = \"darkflat\", \"energies\", \"pattern\", \"indices\", \"radius\", \"track\", \"save\"), or \"all\"\n"+
                        "// *NOTE: Call "+MACRO_COMMAND_UPDATE+"('stack') after "+MACRO_COMMAND_SETVALUE+" commands affecting dark&flat processing and\n"+
                        "//       "+MACRO_COMMAND_UPDATE+"('meta') after "+MACRO_COMMAND_SETVALUE+" commands affecting stack metadata (I0, I00 etc.)! \n\n"+
                        "// Use this Macro Recorder to find out the command arguments. \n\n"+
                        "// For recording all parameter values, use More>>List parameters while the Macro Recorder is open. \n\n");
                }
            lastRecorderHintsShownMillis = millis;
        }
        if (!Recorder.scriptMode() && name!=null)
            Recorder.recordString("Ext."+name+ "("+argument+");\n");
    }

// ----------------  U T I L I T I E S    ----------------
    /** Returns the energy array or null */
    public double[] getEnergyArray() {
        double[] energies = energiesEtc[ENERGY];
        if (energies == null && LeedParams.get(LeedParams.LEEMMODE)!=0) energies = energiesEtc[E_LEEM];
        return energies;
    }

    /** Returns the energy corresponding to a given stack slice.
     *  Note that stack slices start with 1 and run to stack size */
    public double sliceToEnergy(int stackSlice) {
        return sliceToEnergy(getEnergyArray(), stackSlice);
    }

    /** Returns the energy corresponding to a given stack slice.
     *  Note that stack slices start with 1 and run to stack size */
    public static double sliceToEnergy(double[] energies, int stackSlice) {
        double energy = energies == null ? LeedRadiusSelector.UNKNOWN_ENERGY : energies[stackSlice-1];
        return energy;
    }

    /** Returns the stack slice nearest to the given energy in the energies array
     *  (or another array)
     *  Returns -1 if outside the range and 'minusIfOutside' is true.
     *  Also returns -1 if the energies array is null or has only equal values. */
    public static int energyToSlice(double energy, double[] energies, boolean minusIfOutside) {
        if (energies == null)
            return -1;
        int sliceNumber = LeedUtils.getNearestIndex(energy, energies, minusIfOutside);
        if (sliceNumber>=0) sliceNumber++; //stack slices start with 1
        return sliceNumber;
    }
}
