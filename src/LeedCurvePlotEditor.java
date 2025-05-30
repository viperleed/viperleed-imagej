import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.io.*;
import ij.measure.*;
import ij.text.*;
import ij.plugin.*;
import ij.util.Tools;
import ij.util.FloatArray;
import java.awt.*;
import java.awt.geom.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;

/** The main part of the interactive I(V) Curve Editor, called from LEED_Curve_Editor,
 *  This class shows a plot of the I(V) curves and Y functions, based on the ImageJ
 *  Plot class, and buttons for the gui. It also implements the gui functions.
 * 
 *  Implementation note: In contrast to the Spot Tracker, which is a singleton,
 *  there can be multiple instances of the I(V) Curve Editor (this class).
 *  Therefore, we must not use values from the ImageJ Prefs or LeedParams
 *  but one has to pass the parameters to the constructor of this class.
 */

/** This code is part of the ViPErLEED package for LEED I(V) analysis.
 *  Licensed under GNU General Public License v3.0 or later (GPL-3.0-or-later),
 *  https://www.gnu.org/licenses/gpl-3.0.html
 *  The authors may decide later to put part of the auxiliary code in this work into the public domain,
 *  to allow incorporation into ImageJ if desired (ImageJ is in the public domain).
 *  The help text is licensed under the Creative Commons Attribution 4.0 (CC BY 4.0) license.
 *  When using and/or modifying this program for scientific work, please cite
 *  the paper describing it:
 *  M. Schmid, F. Kraushofer, A. M. Imre, T. Kißlinger, L. Hammer, U. Diebold, and M. Riva,
 *  ViPErLEED package II: Spot tracking, extraction and processing of I(V) curves,
 *  Phys. Rev. Research, 2024. 
 *  @author Michael Schmid, IAP/TU Wien, 2020-2024
 */

public class LeedCurvePlotEditor implements Runnable, MouseListener, MouseMotionListener,
        ActionListener, ImageListener, KeyListener, MouseWheelListener {

    public static final double MIN_RANGE_OVER_V0I = 4;                  //we should get contiguous data over at least this length
    public static final double BLEND_OVER_V0I = 4;                      //blending over this energy range
    public static final double INTERPOLATE_OVER_V0I = 0.5;              //for holes less than this, we interpolate in the main LEED_CURVE_EDITOR
    static final double NEGATIVES_TOLERANCE = 1.0;                      //how much the R factor may change (comprared to noisLimit) when the curve shifts up due to neg. intensities
    static final double MAX_EQUIV_I_RATIO = 3.0;                        //if intensity ratio is larger, beams averaging does not take ratio into account
    static final double NOISE_SMOOTH_EXPONENT = 0.15;                   //Noise dependent smoothing scales with this power of avg noise
    static final int    E_DECIMALS = 1;                                 //energies & smooth parameters are shown with 1 decimal
    static final double NOISE_Y_RANGE = 0.25;                           //fixed y axis maximum for noise function

    public static final String PLUGIN_NAME = "LEED I(V) Curve Editor";

    private static final double GUI_SCALE = Prefs.getGuiScale();
    private static final int    SNAP_RANGE = (int)(5*GUI_SCALE);        //snap roi to limits within this (in pixels)
    // buttons
    private static final int N_BUTTONS=9;
    private static final int HELP=0, SAVE=1, PREV=2, NEXT=3, ON_OFF=4, SMOOTH_LESS=5, SMOOTH_MORE=6, ROI=7, OPTIONS=8;
    //private static final String[] BUTTON_NAMES = new String[]
    //        {"Save", "Previous beams(s)", "Next beam(s)", "Beam on/off", "Smooth less", "Smooth more", "Select E range", "Options"};
    private static final int[] BUTTON_X = new int[] {-1, 0, -1, 0, 0, -1, 0, -1, 0};       //button x positions, 0 is rightmost, -1 left of it
    private static final int[] BUTTON_Y = new int[] {0, 0, 1, 1, 2, 3, 3, 4, 4};          //button y positions, 0 is topmost
    private static final Color[] BUTTON_COLOR = new Color[] {                       //button symbol colors (active)
        Color.BLUE, Color.BLUE, Color.BLACK, Color.BLACK, Color.YELLOW, Color.BLUE, Color.BLUE, Color.YELLOW, Color.BLACK
    };
    private static final float[][] BUTTON_SYMBOL_X = new float[][] {                //symbol polygon x coords in 0-1 range
        new float[] {0.65f, 0.925f, 0.775f, 0.775f},                                //'T' next to '?' for help
        new float[] {.75f, .875f, .875f, .125f, .125f, .725f,
                .725f, .15f, .15f, .78f, .78f},                                     //save
        new float[] {.75f, .25f, .75f}, new float[] {.25f, .75f, .25f},             //prev, next
        //new float[] {.25f, .75f, .5f, .75f, .25f}, new float[] {.125f, .375f, .875f}, //no, yes (x & ok hook)
        new float[] {0.5f, 0.5f},                                                   //on/off
        new float[] {.125f, .25f, .375f, .5f, .625f, 0.75f, 0.875f},                //smooth less
        new float[] {.125f, .25f, .375f, .5f, .625f, 0.75f, 0.875f},                //smooth more
        new float[] {.25f, .75f, .75f, .25f, .25f}, null                            //roi, tools ('tools created later)
    };
    private static final float[][] BUTTON_SYMBOL_Y = new float[][] {                //symbol polygon y coords in 0-1 range (0 top)
        new float[] {0.65f, 0.65f, 0.65f, 0.925f},                                  //'T' next to '?' for help
        new float[] {0.125f, .25f, .875f, .875f, .125f, .125f,
                .45f, .45f, .53f, .53f, .2f},//save
        new float[] {.25f, .5f, .75f}, new float[] {.25f, .5f, .75f},               //prev, next
        //new float[] {.25f, .75f, .5f, .25f, .75f}, new float[] {.625f, .875f, .125f}, //no, yes
        new float[] {0.3f, 0.7f},                                                   //on/off
        new float[] {.75f, .875f, .375f, .5f, .1125f, 0.75f, 0.625f},               //smooth less
        new float[] {.875f, .75f, .5f, .25f, .125f, .25f, .5f},                     //smooth more
        new float[] {.25f, .25f, .75f, .75f, .25f}, null                            //roi, tools ('tools created later)
    };
    private static final int BUTTON_SIZE = (int)Math.round(GUI_SCALE*24);
    private static final int BUTTON_SPACE = BUTTON_SIZE/3;
    private static final int MAX_BUTTONS = N_BUTTONS + 12;    //assume no more than 12 equivalent beams
    private static final Font BUTTON_FONT = new Font(Font.DIALOG, Font.BOLD, BUTTON_SIZE*3/4);

    private static final Color C_BG_ACTIVE = Color.GRAY;    //button background colors
    private static final Color C_BG_INACTIVE = Color.LIGHT_GRAY;

    private static final Font smallFont = new Font("SansSerif", Font.PLAIN, (int)Math.round(GUI_SCALE*12)); //Font for beam names
    private static final Font normFont = new Font("SansSerif", Font.PLAIN, (int)Math.round(GUI_SCALE*16));    //Font for group name

    private static final String SET_COMMENT   = "Set comment for group";      //context menu entries
    private static final String DELETE_COMMENT = "Delete comment for group";
    private static final String LAST_COMMENT  = "Set to last comment entered/deleted";
    private static final String LIST_COMMENTS = "List Comments";
    private static final String SAVE_NO_SMOOTH = "Save without smoothing";
    private static final String GO_FIRST      = "Go to first group";
    private static final String GO_LAST       = "Go to highest group";
    private static final String GO_PREV_COMMENT   = "Go to previous group with a comment";
    private static final String GO_NEXT_COMMENT  = "Go to next group with a comment";
    private static final String GO_PREV_BAD   = "Go to previous 'bad beam agreement' ["+LeedUtils.CTRL+" <]";
    private static final String GO_NEXT_BAD   = "Go to next 'bad beam agreement' ["+LeedUtils.CTRL+" >]";
    private static final String GO_HIGHEST    = "Go to highest selected group";
    private static final String GO_INDEX      = "Go to beam index/beam group...";
    private static final String SELECT_ALL    = "Select all groups";
    private static final String DESELECT_ALL  = "Deselect all groups";
    private static final String SET_SMOOTHING = "Set smoothing (default/for all)...";
    private static final String BEST_SMOOTHING = "Set best smoothing [vs. simulated I(V)]";
    private static final String SET_E_RANGE_NO_SNAP = "Set energy range to current ROI (no snap)";
    private static final String DEFAULT_E_RANGE = "Default energy range";
    private static final String FULL_E_RANGE  = "Full energy range";
    private static final String SET_E_LIMITS  = "Set energy limits (default/for all)...";
    private static final String AUTO_SELECT   = "Auto Select...";
    private static final String SET_MIN_E_SPAN = "Set minimum energy span...";
    private static final String OPTION_DIALOG = "Options...";
    private static final String SET_BAD_LIMIT = "Set 'bad-agreement' limit...";
    private static final String LIST_BAD =      "List 'bad-agreement' beams";
    private static final String EXCLUDE_RANGE = "Exclude selected energies";
    private static final String NO_EXCLUDE_BEAM = "Cancel excluded energies";
    private static final String[][] CONTEXT_MENU_ITEMS = new String[][] {   //context menu for each button
            {SET_COMMENT, LAST_COMMENT, DELETE_COMMENT, LIST_COMMENTS},     //HELP (menus for comments)
            {SAVE_NO_SMOOTH},                                               //SAVE
            {GO_FIRST, GO_PREV_BAD, GO_PREV_COMMENT, GO_INDEX},             //PREV
            {GO_LAST, GO_HIGHEST, GO_NEXT_BAD, GO_NEXT_COMMENT, GO_INDEX},  //NEXT
            {SELECT_ALL, DESELECT_ALL, AUTO_SELECT},                        //ON_OFF,
            {SET_SMOOTHING, BEST_SMOOTHING},                                //SMOOTH_LESS
            {SET_SMOOTHING, BEST_SMOOTHING},                                //SMOOTH_MORE
            {SET_E_RANGE_NO_SNAP, DEFAULT_E_RANGE, FULL_E_RANGE, SET_E_LIMITS, SET_MIN_E_SPAN, AUTO_SELECT}, //ROI
            {OPTION_DIALOG, SET_BAD_LIMIT, LIST_BAD},                       //OPTIONS
            {EXCLUDE_RANGE, NO_EXCLUDE_BEAM}                                //beam on/off buttons
    };
    private static final int ROI_REQUIRED = 1;                              //flag: this context menu item needs a roi
    private static final int EXCLUDED_REQUIRED = 2;                         //flag: this context menu item needs an excluded range
    private static final int COMMENT_REQUIRED = 4;                          //flag: this context menu item needs a group with comment
    private static final int ANY_COMMENT_REQUIRED = 8;                      //flag: this context menu item needs at least one comment
    private static final int LAST_COMMENT_REQUIRED = 16;                    //flag: this context menu item needs a previous comment
    private static final int CURRENT_ON_REQUIRED = 32;                      //flag: this context menu item needs the current item selected
    private static final int[][] CONTEXT_MENU_FLAGS = new int[][] {         //for each of the main buttons
            {0, LAST_COMMENT_REQUIRED, COMMENT_REQUIRED, ANY_COMMENT_REQUIRED}, //HELP (menus for comments)
            null,                                                           //SAVE
            {0, 0, ANY_COMMENT_REQUIRED , 0},                               //PREV
            {0, 0, 0, ANY_COMMENT_REQUIRED , 0},                            //NEXT
            null,                                                           //ON_OFF
            {0, CURRENT_ON_REQUIRED}, {0, CURRENT_ON_REQUIRED},             //SMOOTH_LESS, SMOOTH_MORE
            {ROI_REQUIRED, 0, 0, 0, 0, 0},                                  //SELECT
            null,                                                           //OPTIONS
            {ROI_REQUIRED, EXCLUDED_REQUIRED}                               //beam on/off buttons
    };

    private static final int MIN_SMOOTH = LeedSmoother.MIN_HALFWIDTH; //smooth kernel halfwidth (everything below is no smoothing
    private static final int MIN_PLOT_E_POINTS = 200;       //at least this number of points on the x axis

    private static final int N_PLOT_DATA = 4;                //number of data sets in processedData (for plotting)

    private static final int AVG = 0, AVG_SMOOTH = 1, Y_AVG = 2, Y_SMOOTH = 3;    //data sets in processedData
    private static final String[] PLOT_DATA_NAMES = new String[] {"Avg", "Smooth", "Y", "Y_Smooth"};
    private static final Color[] PLOT_COLOR = new Color[] { //colors for additional data sets
        Color.RED, Color.BLACK, Color.RED, Color.BLACK
    };

    private static final String EDIT_FILE_COLUMNS = "spot,useIn,useOut,smooth,start,end,noInStart,noInEnd,comment";

    private static final String IMP_PROPERTY_KEY = "LEED_curveEditor";  //ImagePlus has the reference to the LeedCurvePlotEditor under this key

    ImagePlus  plotImp;

    private Rectangle[] buttonRects = new Rectangle[MAX_BUTTONS];
    private boolean[] buttonEnabled = new boolean[MAX_BUTTONS];
    private boolean[] buttonActive = new boolean[MAX_BUTTONS];
    private int mouseOverButton = -1;
    static int smallestMouseWheelStep = 100; //the smallest mouse wheel increment seen so far (increments are machine-specific)

    private Thread guiBlockingThread;   //while this thread is alive, we keep the gui blocked (used for findBestSmoothing)
    private Thread spotTrackerHighlightThread;  //highlighting in the Spot tracker can be slow, it runs in a separate thread
    private String spotTrackerHighlightName;    //name of the beam that should be highlighted in the SpotTracker (together with its group)

    int previousWidth, previousHeight;    //previous size of plot, to detect changed size
    double plotXmin, plotXmax, plotYmin, plotYmax;             //current y range, to react to changed y range (for plotting 'Y' function)
    double plotXdefMin, plotXdefMax, plotYdefMin, plotYdefMax; //default plot limits for current group, remember to get 'R' of plot working

    private String pathPrefix;          //directory and name of input file excluding "_Int.csv"
    private String editFilePath;        //directory and name of file with previous edit parameters, or null
    private boolean changesDone;        //whether data have changed and we should write a new Edit file
    //input data
    LeedSpotPattern spotPattern;
    double[]   energies;        //energies (same for all columns = beams)
    double     eStep;           //energy step
    double[][] intData;         //intensity data for each column (=each beam)
    int[]      spotIndices;     //spot (=beam) number in spotPattern for each column of the data
    //parameters (we must not read these from prefs since we can have several instances of this class)
    double     v0i, smoothEv, minSpanEv, curveStartEv, curveEndEv, noiseLimit;
    boolean    autoSelect;      //start with auto selection of groups and ranges
    int        maxSmooth;       //no more smoothing like this (kernel halfwidth)
    boolean    plotNoise;       //option: whether to plot the noise function
    boolean    showSpotTrackerEnergy = true; //whether to show the energy of the SpotTracker
    boolean    highlightInSpotTracker = true; //whether to highlight beams in the SpotTracker (may be slow)
    boolean    fixedXrange;     //option: whether to fix the x range
    boolean    noiseDependentSmooth=true;  //remembers last setting of 'noise-dependent smoothing'
    String     lastComment;     //last comment entered or deleted
    double     badnessLimit;    //limit for check of bad R factors of inequivalent beams
    String     theobeamsPath;   //theobeams.csv file to find best smoothing
    LeedIntegerArray badBeamColumns;    //column indices of beams not fitting the average
    FloatArray badEnergies;     //energies of the worst agreement between equivalent curves
    FloatArray badnessValues;   //value of the worst agreement between equivalent curves

    int rFactorType = LeedRFactor.R_PENDRY;

    //setting for individual beams or groups (these are in the edit file)
    boolean[]  useInput;        //for the beams (=data columns), whether the beam is selected as input
    boolean[]  useOutput;       //for the first beam of each group whether the group is selected as output
    int[]      smoothR;         //smooth radius (relevant only for the first beam of each group)
    int[][]    outputRange;     //for the first beam of each group, start (inclusive) & end (exclusive) indices for output. Not null.
    int[][]    noInputRange;    //for each input beam, start, end of excluded range. null=fullrange
    int[]      selectedRangeLengths;   //length of actually used range for each group (uses column indices)
    String[]   comment;         //optional comment for each beam group

    int[] defaultOutputRange;   //indices in energy array for start and (exclusive) end. Clone before using as output range for individual groups.

    int group = -1;             //currently selected group
    boolean groupChanged = true;
    int[] currentColumns;       //data columns of currently selected group
    Color[] currentColors;
    double spotTrackerEnergy = Double.NaN; //energy of spotTracker to mark
    int lowestGroup, highestGroup;
    int totalGroups, selectedGroups;  //for statistics: all groups in file, number of groups switched on
    double selectedRange;       //cumulated for all groups
    double [][] processedData = new double[N_PLOT_DATA][];  //data to plot in addition to individual curves: averaged and smoothed, Y function

    private static Dialog helpDialog;
    private static final String HELP_STRING =
            "<html>"+
            "<h1>ViPErLEED I(V) Curve Editor</h1>"+
            "<p><b>This ImageJ plugin is used to select usable data from LEED <i>I</i>(<i>V</i>) curves and smooth them.</b></p>"+
            "<p>Version: "+LEED_Spot_Tracker.VERSION+"</p>"+
            "<p>This documentation contains a description of the <a href='#buttons'>GUI elements (buttons)</a>, "+
            "the typical <a href='#workflow'>work flow</a>, <a href='#ivEditorHints'>useful hints</a> "+
            "and <a href='#ivEditorLicense'>license information</a>.</p>"+
            "<h2><a name='mainPlot'>Main Plot</a></h2>"+
            "<ul>"+
            "<li>The <b>main part </b> of the plot area shows the <i>I</i>(<i>V</i>) curve(s) of a group of symmetry-equivalent beams, "+
            "their average (red) and the smoothed average in the selected range (black).</li>"+
            "<li>The <b>top part</b> of the plot area shows the <i>Y</i> function "+
            "(modified logarithmic derivative used in Pendry's <i>R</i> factor), "+
            "for the averaged <i>I</i>(<i>V</i>) curve (red) and the smoothed <i>I</i>(<i>V</i>) curve (black). "+
            "This is the function that is compared to the corresponding function of the calculated <i>I</i>(<i>V</i>) curves.</li>"+
            "<li>Small, colored <b>ticks</b> below the <i>Y</i> function indicate the limits of the data ranges of the individual <i>I</i>(<i>V</i>) curves. "+
            "Thicker and longer ticks indicate an excluded energy range of an input beam.</li>"+
            "</ul>"+
            "<h2><a name='buttons'>Buttons</a></h2>"+
            "<ul>"+
            "<li><b>Help</b> ('?' symbol) Shows this help window. Right-click to set or modify a <b>comment</b> (annotation) for the current group "+
            "of symmetry-equivalent beams. (This is indicated by the small 'T' symbol.) "+
            "Comments are saved in the <a href='#editFile'>edit file</a>.</li>"+
            "<li><b>Save</b> saves the csv file with the final curves. Right-click to save the data without smoothing.</li>"+
            "<li>The <b>&lt;</b> and <b>&gt;</b> buttons switch to the previous/next group of symmetry-equivalent beams. "+
            "Keyboard shortcut: '&lt;' or ',' and '&gt;' or '.' (on non-English keyboards, ',' and '.' keys work only "+
            "if this is the main function of the key, without pressing the shift key). The &lt;Page Up&gt;, &lt;Page Down&gt;, &lt;Home&gt;, "+
            "and &lt;End&gt; keys also work. Alternatively, you can use the mouse wheel.<br>"+
            "Right-click for more navigation options: First, last, or last selected beam group. "+
            "You can also jump to the previous or next beam group with a comment, or "+
            "previous or next group with bad agreement of equivalent beams (for the latter, only selected groups and energy ranges count). "+
            "You can also enter the beam indices or group number for jumping to a group."+
            "</li>"+
            "<li><b>Group On/Off</b> (on/off symbol) selects whether the curves of the currently "+
            "shown group of symmetry-equivalent beams should be used at all.<br>"+
            "Keyboard shortcut: SPACE bar.<br>"+
            "Right-click to select or deselect all groups. "+
            "(Note that 'Select all groups' does not select groups where the energy span is "+
            "less than the <a href='#energyRange'>minimum span</a>.) "+
            "The '<a href='#autoSelect'>Auto Select</a>' function available via right-clicking does not only switch groups on/off "+
            "but also selects the energy ranges, maximizing a figure of merit that depends on the noise.</li>"+
            "<li><b>Smooth More/Less</b> selects how much the curve (or average of curves) will be smoothed. "+
            "You can modify smoothing also with the &lt;ALT&gt; key down and the mouse wheel. Press &lt;SHIFT&gt; for a faster change.<br>"+
            "Right-click to change the (default) smoothing for all deselected (switched-off) groups "+
            "or to apply a given smoothing value to all selected groups "+
            "(or the current group and all higher groups). When applying the smoothing value, you can select 'noise-dependent smoothing'. "+
            "This option sets slightly stronger smoothing for curves that are more noisy than the currently selected one, "+
            "and slightly weaker smoothing for less noisy curves.<br>"+
            "The 'Set best smoothing [vs. simulated I(V)]' option asks for a file of simulated <i>I</i>(<i>V</i>) curves (<tt>THEOBEAMS.csv</tt> file). "+
            "It then tries different smoothing parameters to find the smoothing that results in the lowest <i>R</i> factor "+
            "between smoothed experimental and simulated curves, and applies this smoothing. "+
            "This should be done only as the very last step, after refinement of the structure and inner potential.</li>"+
            "<li><a name='energyRange'><b>Set Energy Range</b></a> (yellow rectangle): "+
            "Use the ImageJ Rectangle tool from the main ImageJ Toolbar to select the energy range and then press this button.<br>"+
            "Keyboard shortcut: '#' or '0'.<br>"+
            "Only the left and right sides of the selection rectangle are taken into account. "+
            "When the left or right side is close to the limit of one of the input curves, "+
            "and a slight shift of the boundary would avoid the necessity to smoothly fade in or fade out of that curve, "+
            "the range snaps to this point. "+
            "You can override the snap with a right-click on the 'Set Energy Range' button and "+
            "'Set energy range to current ROI (no snap)'.<br>"+
            "For a given beam group, only one energy range can be selected. In other words, you can't have a gap in the output data. "+
            "(This restriction is required by TensErLEED, which is used as a backend for "+
            "LEED-<i>I</i>(<i>V</i>) calculations in <tt>viperleed.calc</tt>.)<br>"+
            "When there is no selection and 'Set Energy Range' button is pressed, the default energy range is used. "+
            "This is defined by the limits from the dialog window shown when opening the Curve Editor, "+
            "unless it has been modified by 'Set energy limits' from the context menu of the 'Set Energy Range' button. "+
            "The 'Set Energy Range' button also selects a curve if it was deselected.<br>"+
            "By right-clicking the 'Set Energy Range' button, one can also set the full energy range for the curve, "+
            "or the default energy range (as given by the limits). "+
            "You can also modify these limits (the default energy range) or the minimum energy span "+
            "(previously specified in the starting dialog), "+
            "and optionally restrict the energy ranges of all groups to obey the new limits. "+
            "<b><a href='#autoSelect'>Automatic selection</a></b> of beam groups and energy ranges (depending on the noise) "+
            "is also available from this context menu.<br>"+
            "<i>Hints for manual selection</i> - "+
            "Select regions where the <i>Y</i> function at the top of the plot is not too noisy. "+
            "If there is no sufficiently long low-noise range, deselect the whole curve. "+
            "Rule of thumb: For 0.5 eV energy steps, in the ideal case, the noise (peak\u2013peak) "+
            "should be on average less than ~20\u201330% of the full vertical range of the <i>Y</i> function. "+
            "Select the high-energy limit such that high-energy regions with an average noise larger than this noise are excluded.<br>"+
            "If symmetry-equivalent beams differ substantially in the low-energy region, exclude these low energies. "+
            "(Deviations between symmetry-equivalent beams at low energies only are typically due to residual magnetic fields; "+
            "deviations at all energies are more likely due deviations from normal incidence. "+
            "Very large deviations may indicate that you have the wrong spot pattern file, i.e., the symmetry is wrong.)<br>"+
            "</li>"+
            "<li><a name='ivEditorOptions'><b>Options</b></a> (gearwheel symbol): Use the 'Options...' dialog to change the V0i value, "+
            "use a fixed <i>x</i> axis with the full energy range for all curves, "+
            "or to show also a noise estimate on the plot (blue). When the noise is plotted, "+
            "the <i>y</i> scale for this curve is always fixed (range 0&ndash;"+(float)NOISE_Y_RANGE+") "+
            "and the blue horizontal line is the noise limit for automatic selection.<br>"+
            "The Options... dialog also lets you <b>synchronize</b> the current I(V) Curve Editor with the <b>'SpotTracking'</b> display of "+
            "the Spot Tracker: This shows the current energy of the Spot Tracker as a black vertical line in the I(V) Curve Editor. "+
            "The beams of the current beam group in the I(V) Curve Editor can be highlighted in the Spot Tracker "+
            "(unless the highlighted beams of the Spot Tracker still mark the dubious beams of spot tracking; these are not modified).<br>"+
            "By right-clicking the 'Options' button, you can select the threshold for classifying "+
            "'<b>bad agreement</b>' of symmetry-equivalent beams and "+
            "create a list of these 'bad-agreement' beams. "+
            "A typical value for this threshold is between the threshold value for "+
            "<a href='#autoSelect'>automatic selection</a> and twice that value. "+
            "'Bad agreement' is based on a 'local <i>R</i> factor' "+
            "(within a small energy window) between the beam and the average of all selected beams. Note that Pendry's <i>R</i> factor is "+
            "very sensitive at the minima, much less sensitive at the maxima. You can see the reason for the disagreement "+
            "when moving the mouse over the button of the beam with 'bad agreement' and comparing its <i>Y</i> function to that "+
            "of the average.<br>"+
            "You can also jump to the previous or next group with bad agreement by right-clicking the '&lt;' or '&gt;' button, respectively; "+
            "energies of bad agreement are then flagged by thick, vertical lines.<br>"+
            "'Bad agreement' is analyzed only for selected beams and energy ranges.</li>"+
            "<li><b>Individual curve buttons</b> (colored '<b>-</b>' symbols): Press a button to select/deselect the curve "+
            "of the respective beam (if more than one). "+
            "You can highlight a curve by moving the mouse over the button. "+
            "Then, the status line at the bottom shows the <i>R</i> factor of this beam vs. the average of the currently selected beams.<br>"+
            "If a beam has bad data in a given energy range (e.g., due to a dust grain on the LEED screen), "+
            "select that range with the ImageJ rectangle tool (only the left and right side of the rectangle are taken into account). "+
            "Then right-click the button for this beam and and exclude this energy range with 'Exclude selected energies'. "+
            " In the plot, the excluded energy range for a single beam is indicated by thick, long vertical lines "+
            "below the <i>Y</i> function; then there is a 'Cancel excluded energies' command "+
            "available with right-clicking the button for the respective curve.</li>"+
            "</ul>"+
            "<h2><a name='workflow'>Typical workflow</a></h2>"+
            "<p>If you have taken more than one <i>I</i>(<i>V</i>) movie for averaging, "+
            "average their <i>I</i>(<i>V</i>) curves before using the I(V) Curve Editor.</p>"+
            "<p>When opening the (averaged) <i>I</i>(<i>V</i>) curve for the first time, "+
            "start with <b><a name='autoSelect'>automatic selection</a></b> (also available later with right-clicking the on/off or 'Set Energy Range' button) "+
            "to select the energy ranges for all curves. Automatic selection is based on an "+
            "estimate of the noise of the curve and its impact on the <i>R</i> factor. "+
            "Selection is done such that the average noise is less than the given noise limit and a figure of merit "+
            "is optimized. The figure of merit increases with selected range and decreases with noise. "+
            "Typical values for the noise limit are around 0.05; lower values are more selective, higher values "+
            "lead to a larger total energy range (larger 'database'), at the cost of increased noise.</p>"+
            "<p>The next step is setting the <b>smoothing</b> value. Select a curve that starts at a medium energy and has medium noise. "+
            "Then adjust the smoothing (e.g., using the mouse wheel with the &lt;ALT&gt; key pressed) such that 'real' features of the <i>I</i>(<i>V</i>) "+
            "curves are preserved as best as possible while suppressing the noise. "+
            "(Noise is often visible as wiggles in the smoothed <i>Y</i> curve where the <i>I</i>(<i>V</i>) curves have no apparent wiggles. "+
            "These wiggles in the <i>Y</i> curve should be strongly attenuated.)<br>"+
            "Then right-click on a smoothing button and apply the current smoothing to all curves, with "+
            "'noise-dependent smoothing' enabled. Check that the smoothing is not too strong for the first beam groups "+
            "where low energies are included; if so, adjust smoothing for these groups. "+
            "(Very sharp peaks at low energies mainly occur for 5d metals or if very low energies, below 50&nbsp;eV, are included.)</p>"+
            "<p>Finally, check for <b>bad agreement</b> of symmetry-equivalent beams ("+LeedUtils.CTRL+"-previous group, "+
            LeedUtils.CTRL+"-next group). The energy of the worst agreement of a given beam with the average "+
            "will be marked by a vertical line. You may move the mouse over the button for this beam to compare its "+
            "<i>Y</i> function with that of the average. Sometimes, the bad agreement is caused by a defect of the LEED screen "+
            "(if only one out of many beams is an outlier) and a short section of the 'bad' beam can be excluded: "+
            "Select the 'bad' energy range and exclude it for the beam by right-clicking on its button. "+
            "In other cases, especially if such a 'bad' region is close to the ends of the energy range, "+
            "it makes sense to set the energy range for the beam group such that the region of poor agreement is avoided. "+
            "If many groups show large disagreement of the symmetry-equivalent beams, "+
            "consider repeating the measurement with better adjustment of perpendicular incidence and/or "+
            "better compensation of residual magnetic fields.</p>"+
            "<p>Finally, save the data. You can now use them as input for structure optimization (file <tt>EXPBEAMS.csv</tt>).</p>"+
            "<p>Towards the end of a LEED <i>I</i>(<i>V</i>) study, when you have a good best-fit model from the simulation calculations, "+
            "you may want to try the 'Set best smoothing [vs. simulated I(V)]' option (right-click one of the smoothing buttons) "+
            "for the last tweak (typically with 'noise-dependent smoothing' enabled). Usually, the improvement of the <i>R</i> factor "+
            "will be marginal.</p>"+
            "<h2><a name='ivEditorHints'>Good to know</a></h2>"+
            "<p><a name='editFile'>Even if you do not save the final result, the Curve Editor keeps track of what you did "+
            "and saves an 'Edit File' in the same directory as the input. "+
            "If you close the window and open the editor again, it will propose to use that Edit File to continue.</a></p>"+
            "<p>When the mouse is above a button, the status line at the bottom displays more information.</p>"+
            "<p>When the selected curves are averaged, they are first normalized and the slow trends of the intensity "+
            "between the beginning and end of the overlap region are equalized. "+
            "Additionally, if a curve starts or ends inside the selected energy range, smooth fade-in or fade-out is used to avoid jumps.</p>"+
            "<p>For zooming into the plot, select a rectangle and press the '+' key.<br>"+
            "If the mouse pointer is to the left of the plot area or below it, gray symbols allow you to modify the plot limits (axis ranges). "+
            "The 'R' at the bottom left resets the plot to the original limits, "+
            "and arrow-like triangles near the ends of the axes allow you to shrink or extend the range. "+
            "The 'F' ('full range, fit all') gets everything into the plot window, "+
            "including energies where no data are available for the current curve(s). "+
            "You can set a fixed energy range corresponding to this full range in the <a href='#ivEditorOptions'>Options</a>.</p>"+
            "<p>If you want a different color for the selection rectangle, you can set it in ImageJ: Edit&gt;&gt;Options&gt;&gt;Colors...</p>"+
            "<p>You can have more than one I(V) Curve Editor window at the same time, to compare different sets of <i>I</i>(<i>V</i>) curves. "+
            "These windows are synchronized, i.e., they show the same beam group (if available) and the same energy range on the <i>x</i> axis.</p>"+
            "<p>Smoothing is performed using a 4th-degree modified-sinc smoother [<a href='#msSmooth'>2</a>]. "+
            "The best smoothing parameter is typically 0.6&nbsp;V0i "+
            "(very good data or energies &lt; 50 eV) to 1.3&nbsp;V0i (noisy data). "+
            "The smoothing parameter value (in electronvolts) corresponds to the 'window length' of "+
            "a moving-average filter with the same noise suppression for white noise. (The smoothing algorithm used preserves the "+
            "shape of the <i>I</i>(<i>V</i>) curves much better than a moving-average filter with the same noise suppression.)</p>"+
            "<p>You can copy the contents of any help window like this one ("+LeedUtils.CTRL+"-a to select all and "+LeedUtils.CTRL+"-c). "+
            "Paste it into a text editor for printing, annotating, etc.</p>"+
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
            "<h2>References</h2>"+
            "<p><a name='paper'>[1]</a> M. Schmid, F. Kraushofer, A. M. Imre, T. Kißlinger, L. Hammer, U. Diebold, and M. Riva, "+
            "<i>ViPErLEED package II: Spot tracking, extraction and processing of I(V) curves</i>, "+
            "Phys. Rev. Research, 2024. <a href='https://arxiv.org/abs/2406.18413/'>arXiv:2406.18413</a></p>"+
            "<p><a name='msSmooth'>[2]</a> M. Schmid, D. Rath, and U. Diebold, "+
            "<i>Why and how Savitzky-Golay Filters should be replaced</i>, "+
            "<a href='https://doi.org/10.1021/acsmeasuresciau.1c00054'>ACS Measurement Science Au <b>2</b>, 185 (2022)</a>.</p>"+
            "</html>";

    public LeedCurvePlotEditor(String pathPrefix, LeedSpotPattern spotPattern, LeedIVData ivData, String editFilePath,
            double v0i, double smoothEv, double minSpanEv, double curveStartEv, double curveEndEv, double noiseLimit, boolean autoSelect) {
        this.pathPrefix = pathPrefix;
        this.spotPattern = spotPattern;
        this.energies = ivData.energies;
        this.intData = ivData.data;
        this.spotIndices = ivData.spotIndices; //spot number in spotPattern for each column of the data
        this.editFilePath = editFilePath;
        this.v0i = v0i;
        this.smoothEv = smoothEv;
        this.minSpanEv = minSpanEv;
        this.curveStartEv = curveStartEv;
        this.curveEndEv = curveEndEv;
        this.noiseLimit = noiseLimit;
        this.autoSelect = autoSelect;
    }

    /** The Spot Editor is started by this, usually in a separate thread */
    public void run() {
        try {
            IJ.showStatus("Preparing...");
            useInput = new boolean[spotIndices.length];
            useOutput = new boolean[spotIndices.length];
            noInputRange = new int[spotIndices.length][];
            outputRange = new int[spotIndices.length][];
            smoothR = new int[spotIndices.length];
            selectedRangeLengths = new int[spotIndices.length];
            Arrays.fill(selectedRangeLengths, -1);
            Arrays.fill(useInput, true);
            comment = new String[spotIndices.length];

            eStep = LeedUtils.getEnergyStep(energies);  //energies and energy range defaults...
            if (!(eStep > 0)) {
                IJ.error(PLUGIN_NAME, "ERROR: Energies not evenly ascending");
                return;
            }
            if (v0i < eStep) v0i = eStep;

            defaultOutputRange = energyRangeToIndices(curveStartEv, curveEndEv);
            for (int i=0; i<spotIndices.length; i++)
                outputRange[i] = (int[])defaultOutputRange.clone();
            maxSmooth = LeedSmoother.invNoiseGainSquareToM(2.0*v0i/eStep) + 2;
            int defaultSmooth = LeedSmoother.invNoiseGainSquareToM(smoothEv/eStep); // default for smoothing (kernel halfwidth 'm')
            Arrays.fill(smoothR, defaultSmooth);

            if (editFilePath != null && editFilePath.length() > 0)
                readEdit(editFilePath);                 //read parameters from previous edit file

            int[] lowestHighest = spotPattern.getLowestAndHighestGroup();
            if (lowestHighest == null) {
                IJ.error(PLUGIN_NAME, "ERROR: No data found or beam groups missing");
                return;
            }
            lowestGroup = lowestHighest[0];             //lowest & highest group in spot pattern
            highestGroup = lowestHighest[1];
            group = highestGroup;
            if (getColumns(highestGroup) == null)       //highest & lowest group actually present in data
                highestGroup = getPreviousGroup();
            group = lowestGroup;
            if (getColumns(lowestGroup) == null)
                lowestGroup = getNextGroup();
            if (getColumns(lowestGroup) == null) {
                IJ.error(PLUGIN_NAME, "ERROR: No data found for any group");
                return;
            }
            if (autoSelect)
                for (int g = lowestGroup; g <= highestGroup; g++)
                    autoSelectGroup(g, noiseLimit);

            totalGroups = getNGroups(true);             //some statistics for the status line
            selectedGroups = getNGroups(false);
            selectedRange = getSelectedEnergyRange();   //cumultative range selected

            Arrays.fill(buttonEnabled, true);

            setGroup(lowestGroup, false);               //creates the plot and buttons
            IJ.showStatus("");
            IJ.setTool(Toolbar.RECTANGLE);
        } catch (Exception e) {
            IJ.handleException(e);
        }
    }

    /** Returns an ArrayList of all LeedCurvePlotEditors currently open */
    public static ArrayList<LeedCurvePlotEditor> getAllCurvePlotEditors() {
        ArrayList<LeedCurvePlotEditor> list = new ArrayList<LeedCurvePlotEditor>();
        int[] imageIDs = WindowManager.getIDList();
        for (int imageID : imageIDs) {
            ImagePlus imp = WindowManager.getImage(imageID);
            if (imp == null) continue;
            Object curvePlotEditor = imp.getProperty(IMP_PROPERTY_KEY);
            if (curvePlotEditor instanceof LeedCurvePlotEditor)
                list.add((LeedCurvePlotEditor)curvePlotEditor);
        }
        return list;
    }

    /** Shows the current group in a new Plot and PlotWindow; creates the PlotWindow if required.
     *  'calculate' determines whether the curves should be (re-)calculated */
    void makeAndShowPlot(boolean calculate) {
        int group = this.group;
        int[] columns = this.currentColumns;
        if (LeedUtils.countTrue(useInput, columns)==0)    //if no input columns selected (yet), enable them all
            for (int c : columns)
                useInput[c] = true;

        if (calculate) calculatePlotData(columns);
        Color[] colors = currentColors;
        final Plot plot = makePlot(columns, colors, processedData);

        PlotWindow pw = getPlotWindow();
        if (pw == null) {
            Rectangle max = GUI.getMaxWindowBounds(new Point(10,10));
            if (max != null)
                plot.setSize(max.width*3/4, max.height*3/4);
            pw = plot.show();
            plotImp = pw.getImagePlus();
            ImagePlus.removeImageListener(this);    //make sure we don't register twice; removing too often does not hurt
            ImagePlus.removeImageListener(this);    //especially ImagePlus is critical since it is static
            ImagePlus.addImageListener(this);
            plotImp.setProperty(IMP_PROPERTY_KEY, this);
        } else {
            Plot oldPlot = (Plot)(plotImp.getProperty(Plot.PROPERTY_KEY));
            if (oldPlot != null) {
                int flags = Plot.COPY_SIZE | Plot.COPY_LABELS | Plot.COPY_AXIS_STYLE;
                if (!groupChanged) flags |= Plot.X_RANGE | Plot.Y_RANGE;
                plot.useTemplate(oldPlot, flags);
                double[] limits = plot.getLimits();
                if (limits[0] != plotXdefMin || limits[1] != plotXdefMax || limits[2] != plotYdefMin || limits[3] != plotYdefMax) {
                    plot.setLimits(plotXdefMin, plotXdefMax, plotYdefMin, plotYdefMax); //sets the default limits (a hack to get the 'R' working)
                    plot.draw();                    //The initial draw() always uses the default limits
                    plot.setLimits(limits);         //This does not modify the default limits. It works only after draw().
                    plot.update();
                }
            }
            pw.drawPlot(plot);
            plotImp.draw();
        }
        pw.removeMouseWheelListener(pw);            //avoid vertical panning with mouse wheel
        plotSizeChanged(plot);                      //remember plot size & y range
        if (groupChanged) plotRangeChanged(plot);
        groupChanged = false;

        ImageCanvas canvas = getCanvas();
        if (LeedUtils.arrayIndexOf(canvas.getMouseListeners(), this) < 0) {
            canvas.removeMouseListener(canvas);        //mouse events not handled by this class will be forwarded to the canvas
            canvas.addMouseListener(this);
            canvas.addKeyListener(this);
            canvas.addMouseMotionListener(this);
            canvas.addMouseWheelListener(this);
        }
    }

    /** Checks whether the plot size has been changed and remembers the size */
    boolean plotSizeChanged(Plot plot) {
        boolean changed = false;
        ImageProcessor ip = plot.getProcessor();
        if (previousWidth != ip.getWidth())  changed = true;
        previousWidth = ip.getWidth();
        if (previousHeight != ip.getHeight())  changed = true;
        previousHeight = ip.getHeight();
        return changed;
    }

    /** Checks whether the y range of the plot has changed and remembers the y range */
    boolean plotRangeChanged(Plot plot) {
        boolean changed = false;
        double[] limits = plot.getLimits();         // xMin, xMax, yMin, yMax
        if (plotYmin != limits[2] || plotYmax != limits[3])  changed = true; //we have to create new Y function
        plotXmin = limits[0];
        plotXmax = limits[1];
        plotYmin = limits[2];
        plotYmax = limits[3];
        return changed;
    }

    /** Creates the plot. The input curves to be plotted are given by 'columns' and
     *  for each of them a corresponding color should be provided.
     *  The averaged (and smoothed) curves must be provided. */
    Plot makePlot(int[] columns, Color[] colors, double[][] processedData) {
        String title = PLUGIN_NAME+ " - "+(new File(pathPrefix)).getName();
        Plot plot = new Plot(title, "Energy", "Intensity");
        PlotWindow pw = getPlotWindow();

        float[] floatEnergies = Tools.toFloat(energies);
        double min = 0, max = 0;
        for (int ic=0; ic<columns.length; ic++) {      // plot the input IV curves
            double[] data = intData[columns[ic]];
            double[] minMax = Tools.getMinMax(data);
            if (minMax[0] < min) min = minMax[0];
            if (minMax[1] > max) max = minMax[1];
            plot.setColor(colors[ic]);
            plot.setLineWidth(1);
            plot.addPoints(floatEnergies, Tools.toFloat(data), null, Plot.LINE, spotPattern.getNameWithGroup(spotIndices[columns[ic]], true));
        }
        if (groupChanged) {                         // no previous range, we have to make a new plot
            int firstI = 0, lastI = energies.length;
            if (!fixedXrange) {
                int[] firstLast = getStartEndIndex(columns);
                firstI = firstLast[0];
                lastI = firstLast[1];
                if (firstI > energies.length) firstI = 0;   //may be Integer.MAX_VALUE if empty data column
                if (lastI - firstI < MIN_PLOT_E_POINTS) {        //don't stretch IV curves too much
                    lastI = Math.min(firstI + MIN_PLOT_E_POINTS, energies.length);
                    firstI = Math.max(lastI - MIN_PLOT_E_POINTS, 0);
                }
            }
            double firstE = energies[firstI] - 0.005*(energies[lastI-1] - energies[firstI]);
            double lastE = energies[firstI] + 1.05*(energies[lastI-1] - energies[firstI]);
            max = min+1.33*(max-min);               //leave space for 'Y' curves
            plot.setLimits(firstE, lastE, min, max);
            plotXdefMin = firstE; plotXdefMax = lastE; plotYdefMin = min; plotYdefMax = max; //remember default range for 'R' of plot
        } else {                                    // use range as currently displayed
            min = plotYmin;
            max = plotYmax;
        }
        double minY=Double.MAX_VALUE, maxY=-Double.MAX_VALUE;
        for (int d=Y_AVG; d<N_PLOT_DATA; d++) {     //data range for Y_AVG, Y_SMOOTH
            double[] pData = processedData[d];
            if (pData == null) continue;
            double[] minMax = Tools.getMinMax(pData);
            if (minY > minMax[0]) minY = minMax[0];
            if (maxY < minMax[1]) maxY = minMax[1];
        }
        for (int d=0; d<N_PLOT_DATA; d++) {         //plot additional curves: average, smoothed, 'Y'
            double[] pData = processedData[d];
            if (pData == null) continue;
            if (d >= Y_AVG) {        //scaling of 'Y' data to make them appear in the upper quarter
                for (int i=0; i<pData.length; i++)
                    pData[i] = min + (max-min) * (0.76 + (pData[i]-minY)*(0.235/(maxY-minY)));
            }
            Color color = PLOT_COLOR[d];
            int lineWidth = d == Y_AVG ? 1 : 2;
            if (d == Y_SMOOTH && !useOutput[columns[0]]) {
                color = Color.GRAY;                 //for deselected curves, smoothed Y in thin gray
                lineWidth = 1;
            }
            plot.setColor(color);
            plot.setLineWidth(lineWidth);
            plot.addPoints(floatEnergies, Tools.toFloat(pData), null, Plot.LINE, PLOT_DATA_NAMES[d]);
        }
        if (plotNoise) {
            double[] ivCurve = calculateAverage(columns, /*iEStart=*/0, /*iEEnd=*/energies.length);
            if (ivCurve != null) {
                double[] noise = calculateNoise(ivCurve);
                float[] floatNoise = new float[noise.length];
                double noiseScale = NOISE_Y_RANGE*max/noiseLimit;
                for (int i=0; i<noise.length; i++)
                    floatNoise[i] = (float)(noise[i]*noiseScale);
                plot.setColor(Color.BLUE);
                plot.setLineWidth(1);
                    plot.addPoints(floatEnergies, floatNoise, null, Plot.LINE, "noise (x"+IJ.d2s(noiseScale,2,5)+")");
                plot.drawLine(energies[0], noiseLimit*noiseScale, energies[energies.length-1], noiseLimit*noiseScale);
            }
        }

        double yTick0 = min+0.75*(max-min);         //plot ticks where the curves end
        double yTick1 = min+0.76*(max-min);
        double yTick2 = min+0.6*(max-min);          //long and thick ticks for excluded range
        plot.setLineWidth(1);
        for (int ic=0; ic<columns.length; ic++) {
            int[] rangeLimits = LeedUtils.getRangeLimits(intData[columns[ic]], 0, -1);
            if (rangeLimits[0] < 0) continue;
            for (int j=0; j<rangeLimits.length/2; j++) {
                plot.setLineWidth(1);
                int startIE = rangeLimits[2*j];
                plot.setColor(startIE==0 ? Color.BLACK : colors[ic]);
                plot.drawLine(energies[startIE], yTick0, energies[startIE], yTick1);
                int endIE = rangeLimits[2*j+1];
                plot.setColor(endIE==energies.length ? Color.BLACK : colors[ic]);
                plot.drawLine(energies[endIE-1], yTick0, energies[endIE-1], yTick1);
            }
            int[] excludedRange = noInputRange[columns[ic]];
            if (excludedRange != null) {
                plot.setColor(colors[ic]);
                plot.setLineWidth(3);
                if (excludedRange[0] > 0)
                    plot.drawLine(energies[excludedRange[0]-1], yTick0, energies[excludedRange[0]-1], yTick2);
                if (excludedRange[1] > 0 && excludedRange[1] < energies.length)
                    plot.drawLine(energies[excludedRange[1]], yTick1, energies[excludedRange[1]], yTick2);
            }
        }
        if (showSpotTrackerEnergy && !Double.isNaN(spotTrackerEnergy)) {
            plot.setColor(Color.BLACK);                 //mark Spot Tracker energy
            plot.setLineWidth(1);
            plot.drawLine(spotTrackerEnergy, -0.5*max, spotTrackerEnergy, 2*max);
        }

        int highlightedI = mouseOverButton-N_BUTTONS;   //index of highlighted curve; >= 0 if any

        if (badBeamColumns != null) {                   //mark the bad beams by vertical lines
            for (int ic=0; ic<columns.length; ic++) {
                int indexInBad = badBeamColumns.indexOf(columns[ic]);
                if (indexInBad >= 0) {
                    plot.setLineWidth(highlightedI >= 0 ? 1 : 3);
                    plot.setColor(colors[ic]);
                    double badEnergy = badEnergies.get(indexInBad);
                    plot.drawLine(badEnergy, -0.5*max, badEnergy, 2*max);
                }
            }
        }

        if (highlightedI >= 0) {                        //highlight 'mouseOver' curve
            int highlightedCol = columns[highlightedI];
            double[] data = getIntData(highlightedCol);
            plot.setColor(colors[highlightedI]);
            plot.setLineWidth(3);
            plot.addPoints(energies, data, Plot.LINE);  //thick curve
            plot.setLineWidth(2);
            double[] smoothedHighlighted = smooth(data, smoothR[columns[0]]);
            smoothedHighlighted = LeedUtils.restrictRange(smoothedHighlighted, outputRange[columns[0]][0], outputRange[columns[0]][1]);
            double[] yOfHighlighted = LeedRFactor.getYcurve(smoothedHighlighted, rFactorType, v0i/eStep, null);

            double[] minMax = Tools.getMinMax(yOfHighlighted);
            for (int i=0; i<yOfHighlighted.length; i++)
                    yOfHighlighted[i] = min + (max-min) * (0.76 + 0.235*(yOfHighlighted[i]-minMax[0])/(minMax[1]-minMax[0]));
            plot.addPoints(energies, yOfHighlighted, Plot.LINE); //Y function of highlighted
        }

        return plot;
    }

    /** Returns the data columns of spots in a given group, or null if none.
     *  The columns returned are in ascending sequence */
    int[] getColumns(int group) {
        LeedIntegerArray cols = new LeedIntegerArray(12);
        for (int col=0; col<spotIndices.length; col++) {
            if (spotIndices[col] < 0 || intData[col] == null) continue;
            if (spotPattern.getGroup(spotIndices[col]) == group)
                cols.add(col);
        }
        return (cols.size() > 0) ? cols.toArray() : null;
    }

    /** Returns the index of the next group, or the current one, if there are no more */
    int getNextGroup() {
        for (int g = group+1; g <= highestGroup; g++)
            if (getColumns(g) != null) return g;
        return group;
    }

    /** Returns the index of the previous group, or the current one, if there are none below */
    int getPreviousGroup() {
        for (int g = group-1; g >= lowestGroup; g--)
            if (getColumns(g) != null) return g;
        return group;
    }

    /** Returns group number of the highest enabled group */
    int getHighestEnabledGroup() {
        int group = lowestGroup;
        for (int g = lowestGroup; g <= highestGroup; g++) {
            int[] columns = getColumns(g);
            if (columns!= null && useOutput[columns[0]])
                group = g;
        }
        return group;
    }

    /** Returns the group number for a given spot index or group number encoded as a string,
     *  Integer.MIN_VALUE if not found */
    int getGroup(String str) {
        int group = Integer.MIN_VALUE;
        int spotIndex = spotPattern.getIndex(str);
        if (spotIndex >= 0)
            group = spotPattern.getGroup(spotIndex);
        else {
            double d = Tools.parseDouble(str.trim());
            if (!Double.isNaN(d)) group = (int)d;
        }
        if (group != Integer.MIN_VALUE && getColumns(group) == null)
            group = Integer.MIN_VALUE;    //a group with no data columns can't be used
        return group;
    }


    /** Counts all or the selected groups.
     *  Requires that the data columns for each group are adjacent. */
    int getNGroups(boolean getTotal) {
        int nGroups = 0;
        int lastGroup = Integer.MIN_VALUE;
        for (int col=0; col<spotIndices.length; col++) {
            int spotIndex = spotIndices[col];
            if (spotIndex < 0) continue; //can happen if we have the wrong pattern file
            int group = spotPattern.getGroup(spotIndex);
            if (group != lastGroup) {
                lastGroup = group;
                if (getTotal || useOutput[col])
                    nGroups++;
            }
        }
        return nGroups;
    }

    /** Determines the cumulated selected energy range of all beams */
    double getSelectedEnergyRange() {
        int range = 0;
        for (int g = lowestGroup; g <= highestGroup; g++) {
            int[] columns = getColumns(g);
            if (columns != null && useOutput[columns[0]]) {
                int colRange = selectedRangeLengths[columns[0]];
                if (colRange < 0)  //range not calculated yet
                    selectedRangeLengths[columns[0]] = (int)Math.round(getSpan(columns) / eStep);
                range += selectedRangeLengths[columns[0]];
            }
        }
        return range*eStep;
    }

    /** Calculates average, smoothed, and 'Y' functions for plotting.
     *  Smoothed functions are displayed only if the beam is on. */
    void calculatePlotData(int[] columns) {
        boolean inUse = useOutput[columns[0]];
        int iEStart = inUse ? outputRange[columns[0]][0] : 0;  // show full range if deselected
        int iEEnd = inUse ? outputRange[columns[0]][1] : energies.length;
        processedData[AVG] = calculateAverage(columns, iEStart, iEEnd);
        if (processedData[AVG] == null) {
            Arrays.fill(processedData, null);
        } else {
            if (inUse)
                processedData[AVG] = avoidHoles(processedData[AVG], iEStart, iEEnd);
            processedData[Y_AVG] = LeedRFactor.getYcurve(processedData[AVG], rFactorType, v0i/eStep, processedData[Y_AVG]);
            processedData[AVG_SMOOTH] = smooth(processedData[AVG], smoothR[columns[0]]);
            if (useOutput[columns[0]]) {
                processedData[AVG_SMOOTH] = LeedUtils.restrictRange(processedData[AVG_SMOOTH], iEStart, iEEnd);
                processedData[AVG_SMOOTH] = avoidNegatives(processedData[AVG_SMOOTH], iEStart, iEEnd);
                processedData[Y_SMOOTH] = LeedRFactor.getYcurve(processedData[AVG_SMOOTH], rFactorType, v0i/eStep, processedData[Y_SMOOTH]);
                int[] rangeLimits = LeedUtils.getRangeLimits(processedData[AVG_SMOOTH], 0, -1); //also update range
                selectedRangeLengths[columns[0]] = rangeLimits[0] >= 0 ?
                        rangeLimits[1] - rangeLimits[0] : 0;
            } else {
                processedData[Y_SMOOTH] = LeedRFactor.getYcurve(processedData[AVG_SMOOTH], rFactorType, v0i/eStep, processedData[Y_SMOOTH]);
                processedData[AVG_SMOOTH] = null;                   //don't show smoothed intensities if not selected
            }
        }
    }

    /** Returns the intensities for the given column, with the data in its exclusion range replaced by NaNs.
     *  If there is no exclusion range or it does not overlap with the input, returns the original intData
     *  array for the column */
    double[] getIntData(int column) {
        return getIntData(column, noInputRange[column]);
    }

    /** Returns the intensities for the given column, with the data in the exclusion range replaced by NaNs.
     *  The exclusionRange should contain two indices, start (first point where data are excluded) and end
     *  (last point where data are excluded + 1).
     *  If the exclusionRange is null or it does not overlap with the input, returns the original intData
     *  array for the column */
    double[] getIntData(int column, int[] exclusionRange) {
        double[] data = intData[column];
        if (exclusionRange == null || exclusionRange[0] < 0 || exclusionRange[1] < 0)
            return data;
        boolean cloned = false;
        for (int i=exclusionRange[0]; i<Math.min(exclusionRange[1], energies.length); i++) {
            if (!cloned && !Double.isNaN(data[i])) {
                data = (double[])data.clone();
                cloned = true;
            }
            data[i] = Double.NaN;
        }
        return data;
    }

    /** Returns the spot name for a given column */
    String getSpotName(int column) {
        return spotPattern.getName(spotIndices[column]);
    }

    /** Returns the first and last energy (& data) index for plotting the data
     *  as a two-element integer array.
     *  Returns [Integer.MAX_VALUE, -1] if there are no valid data. */
    int[] getStartEndIndex(int[] columns) {
        int start = Integer.MAX_VALUE, end = -1;
        for (int ic=0; ic<columns.length; ic++) {
            if (!useInput[columns[ic]] || intData[columns[ic]]==null)
                continue;
            int[] startEnd = LeedUtils.getStartEndNum(intData[columns[ic]]);
            if (startEnd[0] < 0)
                continue;
            if (startEnd[0] < start)
                start = startEnd[0];
            if (startEnd[1] > end)
                end  = startEnd[1];
        }
        return new int[] {start, end};
    }

    /** Checks the (local) R factors between equivalent beams of the selected beams and ranges.
     *  Returns true if the given group has a local R factor >= ?? and sets the badEnergies,
     *  badBeams arrays.
     *  Note that the R factor is a measure of squared deviations. Thus, for uncorrelated deviations,
     *  averaging over n curves reduces the impact on the average to 1/n. Since the deviations are
     *  not necessarily uncorrelated we only reduce the result with a factor of 1/sqrt(n).
     *  */
    boolean badEquivalentBeams(int group, double badnessLimit) {
        int[] columns = getColumns(group);
        if (columns == null || !useOutput[columns[0]]) return false;
        int nCurves = 0;
        for (int ic=0; ic<columns.length; ic++)
            if (useInput[columns[ic]])
                nCurves ++;
        if (nCurves < 2) return false;                          //we need at least two columns to get an R factor between beams and the average

        badBeamColumns = new LeedIntegerArray();
        badEnergies = new FloatArray();
        badnessValues = new FloatArray();
        FloatArray badRFactors = new FloatArray();
        int smoothR = (int)Math.round(LeedSmoother.invNoiseGainSquareToM(v0i/eStep)); //we compare curves smoothed with 1.0*V0i
        int movingMinimumR = (int)Math.round(0.5*v0i/eStep);    //radius of moving minimum filter to remove spikes of squared Y difference
        int movingAverageR = (int)Math.round(3*v0i/eStep);      //squared Y difference will be smoothed over 3*movingAverageR+1 values
        int oneV0i = (int)Math.round(1.0*v0i/eStep);            //to check whether curves are defined in the vicinity of a point, must be < movingAverageR

        double[] avg = calculateAverage(columns, outputRange[columns[0]][0], outputRange[columns[0]][1]);
        if (avg == null) return false;
        double[] avgSmoothed = smooth(avg, smoothR);
        avgSmoothed = LeedUtils.restrictRange(avgSmoothed, outputRange[columns[0]][0], outputRange[columns[0]][1]);
        double[] avgYFunction = LeedRFactor.getYcurve(avgSmoothed, rFactorType, v0i/eStep, null); //restrict to avoid shift by out-of-range negatives
        for (int ic=0; ic<columns.length; ic++) {
            if (useInput[columns[ic]]) {
                double[] cSmoothed = smooth(getIntData(columns[ic]), smoothR);
                cSmoothed = LeedUtils.restrictRange(cSmoothed, outputRange[columns[0]][0], outputRange[columns[0]][1]);
                double[] cYFunction = LeedRFactor.getYcurve(cSmoothed, rFactorType, v0i/eStep, null);
                double[] yDiffSqr = calculateYDiffSqr(cYFunction, avgYFunction);
                /* We have to eliminate peaks of the Y function difference caused by
                 * slight shifts of the minima of the I(V) curves (the R factor is extremly
                 * sensitive to these, but this does not necessarily indicate poor agreement
                 * between symmetry-equivalent beams). We use a running minimum filter for this */
                double[] yDiffSqrFiltd = LeedUtils.minimumFilter(yDiffSqr, movingMinimumR);
                double[] diffAnalysis = LeedUtils.movingAverage(yDiffSqrFiltd, movingAverageR, false);
                double worstRFactor = diffAnalysis[1];          //maximum value
                if (worstRFactor > badnessLimit) {
                    int iOfWorstRFactor = (int)diffAnalysis[3]; //index of maximum
                    int nCurvesAtBad = 0;
                    for (int ic2=0; ic2<columns.length; ic2++)
                        if (!Double.isNaN(intData[columns[ic2]][iOfWorstRFactor]) &&
                                !Double.isNaN(intData[columns[ic2]][iOfWorstRFactor-oneV0i]) &&
                                !Double.isNaN(intData[columns[ic2]][iOfWorstRFactor+oneV0i]))
                            nCurvesAtBad++;                     //count curves that are well-defined around 'bad' spot
                    worstRFactor /= Math.sqrt(nCurvesAtBad-1);  //if there are many curves, deviations count less
                    if (worstRFactor > badnessLimit) {
                        badBeamColumns.add(columns[ic]);
                        badEnergies.add((float)energies[iOfWorstRFactor]);
                        badnessValues.add((float)worstRFactor);
                    }
                }
            }
        }
        if (badBeamColumns.size() == 0) {
            badBeamColumns = null;
            return false;
        } else
            return true;
    }

    /** After 'go to prev/next bad group', returns a String on the badness of
     *  the group, or null if the badness is less than the limit */
    String getBadnessString() {
        if (badBeamColumns == null) return null;
        int[] columns = currentColumns;
        double worstValue = 0;
        int indexOfWorst = -1;
        for (int ic=0; ic<columns.length; ic++) {
            int indexInBad = badBeamColumns.indexOf(columns[ic]);
            if (indexInBad >= 0) {
                if (badnessValues.get(indexInBad) > worstValue) {
                    worstValue = badnessValues.get(indexInBad);
                    indexOfWorst = indexInBad;
                }
            }
        }
        if (indexOfWorst >= 0) {
            return "Badness "+IJ.d2s(worstValue)+" at "+IJ.d2s(badEnergies.get(indexOfWorst),0)+" eV";
        } else
            return null;
    }


    /** Depending on the noise, calculates and sets automatic limits for a given spot group.
     *  Also selects the group if previously deselected and low noise, deselects the group
     *  if noise is too high. Deselected individual beams are not touched.
     *  Negative group numbers (symmetry-forbidden) are always deselected.
     *  Does nothing if the group does not exist. */
    void autoSelectGroup(int group, double noiseLimit) {
        int[] columns = getColumns(group);
        if (columns == null) return;
        if (LeedUtils.countTrue(useInput, columns)==0)          //if no input columns selected (yet), enable them all
            for (int c : columns)
                useInput[c] = true;

        int[] autoSelectLimits = getAutoSelectLimits(group, noiseLimit);
        if (autoSelectLimits == null) return;

        if (autoSelectLimits.length > 0) {
            if (outputRange[columns[0]][0] != autoSelectLimits[0] ||
                    outputRange[columns[0]][1] != autoSelectLimits[1]) 
                changesDone = true;
            outputRange[columns[0]][0] = autoSelectLimits[0];
            outputRange[columns[0]][1] = autoSelectLimits[1];
        }
        if (autoSelectLimits.length > 0 != useOutput[columns[0]])
            changesDone = true;
        useOutput[columns[0]] = autoSelectLimits.length > 0 && group >= 0;
    }

    /** Returns the indices [start, end[ of "autoselect" for the given group.
     *  Returns null if the group does not exist and an empty array if the group
     *  should be deselected. */
    int[] getAutoSelectLimits(int group, double noiseLimit) {
        int[] columns = getColumns(group);
        if (columns == null) return null;
        double[] ivCurve = calculateAverage(columns, /*iEStart=*/0, /*iEEnd=*/energies.length);
        if (ivCurve == null) return new int[0];             //averaging results in no data, switch group off
        double[] noise = calculateNoise(ivCurve);

        int minSpan = (int)Math.round(minSpanEv/eStep);
        /* Check for negatives in a curve smoothed by 1.0*V0i */
        int smoothMoreR = (int)Math.round(LeedSmoother.invNoiseGainSquareToM(v0i/eStep));
        double[] smoothedMore = smooth(ivCurve, smoothMoreR);

        int[] rangeLimits = LeedUtils.getRangeLimits(noise, defaultOutputRange[0], defaultOutputRange[1]);
        for (int iRange=0; iRange<rangeLimits.length/2; iRange++) {
            int rStart = rangeLimits[2*iRange];
            int rEnd = rangeLimits[2*iRange+1];
            //if(rStart>=0)IJ.log("["+group+"] Check Minima in Range "+iRange+"/"+(rangeLimits.length/2)+": "+energies[rStart]+"-"+energies[rEnd-1]);
            if (rEnd - rStart < minSpan) continue;

            final int MIN=0, START=1, END=2;                //indices in ArrayList arrays
            ArrayList<int[]>minimaList = new ArrayList<int[]>();
            double minValue = 0;
            int nNegatives = 0;
            int iAtMinimum = -1;
            boolean atStart = true;
            for (int i=rStart; i<=rEnd; i++) {              //find all minima with values < 0
                if (i<rEnd && smoothedMore[i] < 0) {
                    nNegatives++;
                    if (smoothedMore[i] < minValue) {
                        minValue = smoothedMore[i];
                        iAtMinimum = i;
                    }
                } else {
                    if (nNegatives > 0) {         //new minimum, but ignore very short negative sections (except at boundaries)
                        if (nNegatives*eStep >= 0.4*v0i || atStart || i==rEnd)
                            minimaList.add(new int[]{iAtMinimum, 0, 0}); 
                        nNegatives = 0;
                        minValue = 0;
                    }
                    atStart = false;
                }
            }
            /* Assuming that a negative value poses a problem, find the range where the values surrounding the
             * minimum are not trustworthy ("bad"): less positive than 1.5*|minimum value|.
             * Note that all intervals are in the form [start, end[:
             * E.g., 'badEnd' is the index of the first non-bad point. */
            for (int iMin=0; iMin<minimaList.size(); iMin++) {
                int[] minAndRange = minimaList.get(iMin);
                int i0 = minAndRange[MIN];
                minValue = smoothedMore[i0];
                int badStart = i0;                          //find the "bad" range towards the left
                while (badStart > rangeLimits[2*iRange] && smoothedMore[badStart-1] < -1.5*minValue) badStart--;
                minAndRange[START] = badStart;
                int badEnd = i0+1;                          //find the "bad" range towards the left
                while (badEnd < rangeLimits[2*iRange+1] && smoothedMore[badEnd] < -1.5*minValue) badEnd++;
                minAndRange[END] = badEnd;
            }
            for (int iMin0=0; iMin0<minimaList.size()-1; iMin0++) { //merge minima with overlapping ranges
                int[] minAndRange0 = minimaList.get(iMin0);
                for (int iMin1=iMin0+1; iMin1<minimaList.size(); iMin1++) {
                    int[] minAndRange1 = minimaList.get(iMin1);
                    if (minAndRange0[END] >= minAndRange1[START]) {
                        minAndRange0[START] = Math.min(minAndRange0[START],minAndRange1[START]);
                        minAndRange0[END]   = Math.max(minAndRange0[END],  minAndRange1[END]);
                        minAndRange0[MIN] = smoothedMore[minAndRange0[MIN]] < smoothedMore[minAndRange0[MIN]] ?
                                minAndRange0[MIN] :  minAndRange1[MIN];
                        minimaList.remove(iMin1);
                        iMin0 = -1;                             //after merging, restart everything (maybe also an earlier minimum can be merged in)
                        break;                                  //(leave inner loop, restart outer loop with iMin0 = 0)
                    }
                }
            }

            for (int iMin=0; iMin<minimaList.size(); iMin++) {  //close to the start of the data range?
                int[] minAndRange = minimaList.get(iMin);
                if ((minAndRange[START] - rangeLimits[2*iRange])*eStep < v0i) { 
                    rangeLimits[2*iRange] = minAndRange[END];   //then skip up to the first good point
                    minimaList.remove(iMin);
                    iMin--;
                }
            }
            for (int iMin=minimaList.size()-1; iMin>=0; iMin--) {   //close to the end of the data range?
                int[] minAndRange = minimaList.get(iMin);
                if ((rangeLimits[2*iRange+1] - minAndRange[END])*eStep < v0i) {
                    rangeLimits[2*iRange+1] = minAndRange[START];   //then skip bad points to the end
                    minimaList.remove(iMin);
                }
            }
            if (rangeLimits[2*iRange+1] - rangeLimits[2*iRange] < minSpan) continue;    //range got too short
    
            if (minimaList.size() == 0) continue;
            double[] minValues = new double[minimaList.size()];
            for (int iMin=0; iMin<minimaList.size(); iMin++)
                minValues[iMin] = smoothedMore[minimaList.get(iMin)[MIN]];
            int[] minRanks = Tools.rank(minValues);             //the worst (deepest) minimum
            int lowestMinI = minRanks[0];
            double lowestMinValue = minValues[lowestMinI];
            int leftGoodEnd = rangeLimits[2*iRange];            //find ranges where we can compare the R factor with and without shifting up
            while (smoothedMore[leftGoodEnd] > 0.2*lowestMinValue) //(no range check: should stop before lowestMinValue)
                leftGoodEnd++;       
            while (leftGoodEnd>rangeLimits[2*iRange] && smoothedMore[leftGoodEnd-1]<=0)
                leftGoodEnd--;                                  //don't accept even slight negatives at the end (only in the interior)
            int rightGoodStart = rangeLimits[2*iRange+1];
            while (smoothedMore[rightGoodStart-1] > 0.2*lowestMinValue)
                rightGoodStart--;
            while (rightGoodStart < rangeLimits[2*iRange+1]-1 && smoothedMore[rightGoodStart]<=0)
                rightGoodStart++;
            boolean leftIsGood = (leftGoodEnd - rangeLimits[2*iRange])*eStep > 4*v0i;
            boolean rightIsGood = (rangeLimits[2*iRange+1] - rightGoodStart)*eStep > 4*v0i;

            if (leftIsGood || rightIsGood) {                    //we have enough data range to check whether shifting would hurt the R factor
                double[] smoothedMoreShifted = (double[])smoothedMore.clone();
                for (int i=0; i<smoothedMoreShifted.length; i++)
                    smoothedMoreShifted[i] += -lowestMinValue;
                double[] leftRFactorData=null, rightRFactorData=null;
                if (leftIsGood)
                    leftRFactorData = LeedRFactor.getRFactor(smoothedMore, smoothedMoreShifted, rangeLimits[2*iRange], leftGoodEnd, 0, v0i/eStep, null);
                if (rightIsGood)
                    rightRFactorData = LeedRFactor.getRFactor(smoothedMore, smoothedMoreShifted, rightGoodStart, rangeLimits[2*iRange+1], 0, v0i/eStep, null);
                double leftRFactor = leftIsGood ? leftRFactorData[rFactorType] : Double.NaN;
                double rightRFactor = rightIsGood ? rightRFactorData[rFactorType] : Double.NaN;
                double rFactor = Double.NaN;
                if (leftIsGood && rightIsGood) {
                    final int N_OVERLAP = LeedRFactor.N_OVERLAP;
                    rFactor = (leftRFactor*leftRFactorData[N_OVERLAP] + rightRFactor*rightRFactorData[N_OVERLAP])/
                            (leftRFactorData[N_OVERLAP] + rightRFactorData[N_OVERLAP]); //average over whole range
                } else
                    rFactor = leftIsGood ? leftRFactor : rightRFactor;
                    
                if (rFactor > NEGATIVES_TOLERANCE*noiseLimit) {
                    int[] minAndRange = minimaList.get(lowestMinI);
                    if (minAndRange[START] - rangeLimits[2*iRange] < minSpan) {
                        rangeLimits[2*iRange] = minAndRange[END];       //if not enough space to the left, shrink range
                        iRange--; continue;                             //and (re-)analyze this smaller range
                    } else if (rangeLimits[2*iRange+1] - minAndRange[END] < minSpan) {
                        rangeLimits[2*iRange+1] = minAndRange[START];   //if not enough space to the right, shrink range
                        iRange--; continue;                             //and (re-)analyze this smaller range
                    } else {                                            //if space at left&right, split the range and redo
                        rangeLimits = splitRangeLimits(rangeLimits, iRange, minAndRange[START], minAndRange[END]);
                        iRange--; continue;                             //and (re-)analyze the two new ranges
                    }
                } else
                    continue;                               //upshift to remove the negative value is harmless for the R factor
            } else {
                rangeLimits[2*iRange] = -1;                 //if the 'good' ranges are too short at the left&right, don't use the data at all
                rangeLimits[2*iRange+1] = -1;
            }
        } //for range

        /* We define a figure of merit (FoM) as length of the data range/(average noise + offset)
         * where offset = 0.5*noiseLimit. The energy ranges selected should maximize the FoM.
         * Note that the FoM increases if we extend the range by one point and noise+offset for
         * the new point is less than 2 * the average of noise+offset over the current range */
        final double offset = 0.5*noiseLimit;
        //IJ.log("noiseLimit="+noiseLimit+" offset="+(float)offset+", "+(rangeLimits.length/2)+" range(s):"+energies[rangeLimits[0]]+"..."+energies[rangeLimits[rangeLimits.length-1]-1]);
        double bestFoM = -Double.MAX_VALUE;
        int bestRangeStart = -1, bestRangeEnd = -1;             //the range with ther best FoM found so far
        double[] rSumNoise = new double[noise.length+1];        //running sum over noise+offset
        for (int iRange=0; iRange<rangeLimits.length/2; iRange++) {
            int rStart = rangeLimits[2*iRange];
            int rEnd = rangeLimits[2*iRange+1];
            if (rEnd - rStart < minSpan) continue;
            double rSum = 0;                                    //running sum over noise+offset, to accelerate calculation of averages
            for (int i=rStart; i<rEnd; i++) {
                rSum += noise[i] + offset;
                rSumNoise[i+1] = rSum;
            }
            int iStart = rStart;
            while (iStart < rEnd - minSpan) {
                while (iStart < rEnd && noise[iStart] > noiseLimit)
                    iStart++;                                   //skip noisy start (if any)
                int nextIStart = -1;
                boolean lastNoisy = false;                      //(actually true, but 'false' prevents setting nextIStart already now)
                for (int i=iStart+1; i<=rEnd; i++) {
                    boolean noisy = i == rEnd || noise[i] > noiseLimit;
                    if (noisy && !lastNoisy || i == rEnd) {                  //low-noise section ends
                        //IJ.log(energies[iStart]+"-"+energies[i-1]+" (n="+(i-iStart)+"): try to extend. avgNoise="+IJ.d2s((rSumNoise[i] - rSumNoise[iStart] - (i-iStart)*offset)/(i-iStart),4));
                        int iFirst = iStart, iStop = i;
                        boolean extendLeft = false, extendRight = false;
                        do {                                    //can we extend the range to the left and/or right?
                            extendLeft = iFirst > rStart &&
                                    (noise[iFirst-1]+offset)*(iStop-iFirst+1) < 2*(rSumNoise[iStop] - rSumNoise[iFirst-1]) &&   //FoM must not decrease
                                    (rSumNoise[iStop] - rSumNoise[iFirst-1]) < (noiseLimit +  offset)*(iStop-iFirst+1);         //avg noise below limit
                            extendRight = iStop < rEnd-1 &&
                                    (noise[iStop]+offset)*(iStop-iFirst+1) < 2*(rSumNoise[iStop+1] - rSumNoise[iFirst]) &&
                                    (rSumNoise[iStop+1] - rSumNoise[iFirst]) < (noiseLimit +  offset)*(iStop-iFirst+1);
                            if (extendLeft && extendRight) {    //if we can extend to both directions, extend in the direction of lower noise
                                if (noise[iFirst-1] < noise[iStop])
                                    extendRight = false;
                                else
                                    extendLeft = false;
                            }
                            if (extendLeft)  iFirst--;          //we can extend the range to the left without decreasing the FoM
                            if (extendRight) iStop++;           //we can extend the range to the right
                        } while (extendLeft || extendRight);
                        //IJ.log(energies[iFirst]+"-"+energies[iStop-1]+" pre-snap");
                        int[] range = new int[]{iFirst, iStop};
                        int[] snappedRange = snapRange(range, columns, (int)Math.round(0.5*v0i/eStep));
                        iFirst = snappedRange[0];               //new range may be slightly lower if this avoids fade in/fade out
                        iStop = snappedRange[1];
                        int n = iStop-iFirst;
                        if (n >= minSpan) {
                            double figureOfMerit = n*n/(rSumNoise[iStop] - rSumNoise[iFirst]);
                            //IJ.log(energies[iFirst]+"-"+energies[iStop-1]+" (n="+n+"): lastNoisy="+lastNoisy+" avgNoise="+IJ.d2s((rSumNoise[iStop] - rSumNoise[iFirst] - n*offset)/n,4)+"m FoM="+IJ.d2s(figureOfMerit)+" bestFoM="+IJ.d2s(bestFoM));
                            if (figureOfMerit > bestFoM && (rSumNoise[iStop] - rSumNoise[iFirst]) < (noiseLimit + offset)*n) {
                                bestFoM = figureOfMerit;        //low-noise section and better than the previous
                                bestRangeStart = iFirst;
                                bestRangeEnd = iStop;
                            }
                        }
                    } else if (lastNoisy && !noisy && nextIStart < 0)
                        nextIStart = i;
                    lastNoisy = noisy;
                }
                if (nextIStart >= 0)
                    iStart = nextIStart;
                else
                    break;
            }
        }

        if (bestRangeStart >= 0)
            return new int[] {bestRangeStart, bestRangeEnd};
        else
            return new int[0];
    }

    /** Retuns a new array of start, end, start, end... where the 'iRange'-th pair is split
     *  and a hole is inserted */
    int[] splitRangeLimits(int[] rangeLimits, int iRange, int holeStart, int holeEnd) {
        int[] newRangeLimits = new int[rangeLimits.length+2];
        System.arraycopy(rangeLimits, 0, newRangeLimits, 0, 2*iRange+1);
        newRangeLimits[2*iRange+1] = holeStart;
        newRangeLimits[2*iRange+2] = holeEnd;
        System.arraycopy(rangeLimits, 2*iRange+1, newRangeLimits, 2*iRange+3, rangeLimits.length-(2*iRange+1));
        return newRangeLimits;
    }

    /** Calculates the average over equivalent beams, with a soft transition
     *  (fade in/fade out) where data for beams end.
     *  Columns with useInput[col]=false and excluded input data are ignored.
     *  Does not apply output limits. */
    double[] calculateAverage(int[] columns, int iEStart, int iEEnd) {
        ArrayList<double[]> dataList = new ArrayList<double[]>(columns.length);
        for (int ic=0; ic<columns.length; ic++)
            if (useInput[columns[ic]])
                dataList.add(getIntData(columns[ic]));
        return LeedCurveAverager.calculateAverage(dataList, iEStart, iEEnd, /*minOverlap=*/2, v0i/eStep);
    }

    /** Returns smoothed data, or the original array if radius is <= MIN_SMOOTH. */
    double[] smooth(double[] data, int radius) {
        if (radius <= MIN_SMOOTH)
            return data;
        LeedSmoother smoother = new LeedSmoother(radius);
        double[] result = smoother.smooth(data);
        return result;
    }

    /** Avoids negative values in the range [iEStart, iEEnd[ by adding a constant,
     *  i.e. shifting up the curve up if necessary. Does not overwrite the data */
    double[] avoidNegatives(double[] data, int iEStart, int iEEnd) {
        double min=0;
        for (int i=iEStart; i<iEEnd; i++)
            if (data[i] < min) min = data[i];
        if (min < 0) {                              //negative values?
            data = (double[])data.clone();          //don't overwrite the original
            for (int i=0; i<data.length; i++)
                data[i] -= min;
        }
        return data;
    }

    /** Returns an estimate of the average noise impact on the R factor in the range given
     *  by the 'start' (inclusive) and 'end' (exclusive) indices for the given group. */
    double estimateNoise(int group, int start, int end) {
        int[] columns = getColumns(group);
        double[] ivCurve = calculateAverage(columns, start,end);
        return estimateNoise(ivCurve, start, end);
    }

    /** Returns an estimate of the average noise impact on the R factor in the range
     *  given by the 'start' (inclusive) and 'end' (exclusive) indices, for the given
     *  (unsmoothed) ivCurve. */
    double estimateNoise(double[] ivCurve, int start, int end) {
        if (ivCurve == null) return Double.NaN;
        double[] noise = calculateNoise(ivCurve);
        int n = 0;
        double sum = 0;
        for (int i=start; i<end; i++)
            if (!Double.isNaN(noise[i])) {
                sum += noise[i];
                n++;
            }
        return sum/n;
    }

    /** Returns an array with a (smoothed) noise estimate for the 'Y' function of the R factor,
     *  assuming that smoothing will be done with a smoothing parameter of 1.0*V0i.
     *  The output is calculated in a manner similar to an energy-dependent R factor
     *  between the raw and smoothed curve. It is calculated from the squared difference
     *  of the Y functions for the raw and smoothed input, with a smoothing parameter
     *  of 0.55*V0i, which should have almost no effect on noise-free I(V) curves.
     *  Since the R factor is based on square deviations of the Y function,
     *  the impact of white noise should be approximately proportional to the bandwidth.
     *  Taking 1/v0i as a bandwidth of unity, the (noise-relevant) bandwidths of the
     *  two input curves for 'calculateNoise' are v0i/eStep (raw) and 
     *  (1/invNoiseGainSqr)*v0i/eStep. Thus, assuming that the final I(V) curve will
     *  be smoothed with a smoothing parameter of 1.0*v0i (bandwidth=1 in thise units)
     *  and white noise, the squared difference mentioned above (raw 'noise' function)
     *  is higher than the final R factor by a factor of (1-1/invNoiseGainSqr)*v0i/eStep.*/
    double[] calculateNoise(double[] ivCurve) {
        double invNoiseGainSqr = 0.55*v0i/eStep;                        //weak smoothing that should not affect noise-free I(V) curves
        int smoothR = (int)Math.round(LeedSmoother.invNoiseGainSquareToM(invNoiseGainSqr));
        invNoiseGainSqr = LeedSmoother.mToInvNoiseGainSquare(smoothR);  //actual value (after rounding smoothR)
        double[] rawY = LeedRFactor.getYcurve(ivCurve, rFactorType, v0i/eStep, null);
        double[] smoothed = smooth(ivCurve, smoothR);
        double[] smoothY = LeedRFactor.getYcurve(smoothed, rFactorType, v0i/eStep, null);
        double[] yDifference = calculateYDiffSqr(rawY, smoothY);
        int movingAvgRadius = (int)Math.round(v0i/eStep);
        double[] noise = LeedUtils.movingAverage(yDifference, movingAvgRadius, true);
        for (int i=0; i<noise.length; i++)
            noise[i] *= eStep/(v0i*(1-1/invNoiseGainSqr));              //correct for bandwidths assuming white noise
        int[] rangeLimits = LeedUtils.getRangeLimits(noise, 0, noise.length);
        for (int j=0; j<rangeLimits.length/2; j++) {
            if (rangeLimits[2*j] > 0)
                noise[rangeLimits[2*j]-1] = noise[rangeLimits[2*j]];    //extend to previous point since the Y curve misses one point at each end
            if (rangeLimits[0] >= 0 && rangeLimits[2*j+1] < noise.length)
                noise[rangeLimits[2*j+1]] = noise[rangeLimits[2*j+1]-1];//extend to next point
        }
        return noise;
    }

    /** Returns an array of normalized squared deviations between the Y functions
     *  of the raw and smoothed curve; this is a measure of the noise.
     *  Normalization is similar to the Pendry R factor (divided by the 2*average square of smoothY).
     *  Thus, the is roughly like an "energy-dependent R factor" between the smoothed and
     *  unsmoothed curve. */
    double[] calculateYDiffSqr(double[] rawY, double[] smoothY) {
        int n = smoothY.length;
        double sumY2 = 0;
        int nData = 0;
        for (int i=0; i<n; i++)
            if (!Double.isNaN(smoothY[i])) {
                sumY2 += smoothY[i]*smoothY[i];
                nData++;
            }
        double norm = 0.5*nData/sumY2;     //reciprocal of average 2*sumY2/smoothY.length;
        double[] deviationSqr = new double[n];
        for (int i=0; i<n; i++)
            deviationSqr[i] = sqr(rawY[i] - smoothY[i])*norm;
        return deviationSqr;
    }

    /** If there are holes in the data, keeps only the longest valid range
     *  and sets the data in the shorter ranges to NaN. This is required
     *  because TensErLEED can have only one contiguous range per beam. */
    double[] avoidHoles(double[] data, int iEStart, int iEEnd) {
        int[] rangeLimits = LeedUtils.getRangeLimits(data, iEStart, iEEnd);
        int nRanges = rangeLimits.length/2;
        if (nRanges <= 1) return data; //only one (or no) range, nothing to do
        data = (double[])data.clone(); //copy, don't modify the original
        int maxLength = 0, jOfMax = 0;
        for (int j=0; j<nRanges; j++) {
            int length = rangeLimits[2*j+1] - rangeLimits[2*j];
            if (length > maxLength) {
                jOfMax = j;
                maxLength = length;
            }
        }
        for (int j=0; j<nRanges; j++) {
            if (j != jOfMax)
                Arrays.fill(data, rangeLimits[2*j], rangeLimits[2*j+1], Double.NaN);
        }
        return data;
    }

    /** Returns the R factor of a given curve (with a given index for the current group) vs the average */
    double getRvsAverage(int curveIndex) {
        int[] columns = currentColumns;
        int nOtherCurves = 0;
        for (int ic=0; ic<columns.length; ic++)
            if (ic != curveIndex && useInput[columns[ic]])
                nOtherCurves++;
        if (nOtherCurves == 0)
            return Double.NaN;   //nothing to compare with

        int iEStart = outputRange[columns[0]][0];
        int iEEnd = outputRange[columns[0]][1];
        int[] rangeLimits = LeedUtils.getRangeLimits(intData[columns[curveIndex]], iEStart, iEEnd);
        if (rangeLimits[0] < 0)
            return Double.NaN;   //no data in range
        if (iEStart < rangeLimits[0]) iEStart = rangeLimits[0];
        if (iEEnd > rangeLimits[1]) iEEnd = rangeLimits[1];
        double[] smoothed1 = smooth(processedData[AVG], smoothR[columns[0]]);
        smoothed1 = LeedUtils.restrictRange(smoothed1, iEStart, iEEnd);
        double[] yCurve1 = LeedRFactor.getYcurve(smoothed1, rFactorType, v0i/eStep, null);
        double[] smoothed2 = smooth(intData[columns[curveIndex]], smoothR[columns[0]]);
        smoothed2 = LeedUtils.restrictRange(smoothed2, iEStart, iEEnd);
        double[] yCurve2 = LeedRFactor.getYcurve(smoothed2, rFactorType, v0i/eStep, null);
        return LeedRFactor.getRFactor(yCurve1, yCurve2);
    }

    /** Finds the smoothing value with the lowest R factor against the calculated
     *  I(V) curve 'theobeams'...
     *  Since this may take a bit and we want a progress indicator, we  */
    void findBestSmoothing(final String theobeamsPath, final boolean noiseDependentSmooth, final int firstGroup, final int lastGroup) {
        IJ.showProgress(0.0);
        showStatus("Find Best Smoothing...  WAIT");
        guiBlockingThread = new Thread(new Runnable() {public void run() {
                findBestSmoothing1(theobeamsPath, noiseDependentSmooth, firstGroup, lastGroup);
            }}, "findBestSmoothing");
        guiBlockingThread.start();
    }

    /** Finds the smoothing value with the lowest R factor against the calculated
     *  I(V) curve 'theobeams' and applies this value (a solid shift to account for V0r
     *  is allowed).
     *  The range of groups, firstGroup--lastGroup (inclusive) determines for which
     *  groups the R factor will be calculated and the smoothing will be applied. */
    void findBestSmoothing1(final String theobeamsPath, final boolean noiseDependentSmooth, final int firstGroup, final int lastGroup) {
        /* get ivData and prepare them for R factor measurement */
        LeedIVData theoData = LeedIVData.fromFile(theobeamsPath, spotPattern,
                LeedIVData.E_ASCENDING|LeedIVData.ZERO_IS_NAN);     //IJ.error on failure
        if (theoData == null) {
            IJ.showProgress(1.0);
            return;
        }
        LeedIVData unsmoothed = makeIVData(null, /*smooth=*/false, firstGroup, lastGroup);
        Object commonData = LEED_R_Factor_Between_Datasets.getCommonIVData(
                new LeedIVData[] {theoData, unsmoothed}, spotPattern);
        if (commonData instanceof String) {
            IJ.error(LEED_Curve_Editor.PLUGIN_NAME, (String)commonData);
            IJ.showProgress(1.0);
            return;
        }
        LeedIVData theoDataForR = ((LeedIVData[])commonData)[0];
        LeedIVData unsmoothedForR = ((LeedIVData[])commonData)[1];
        
        /* preparation of (noise-dependent) smoothing */
        int currentSmoothR = smoothR[currentColumns[0]];
        double eVcurrentSmooth = LeedSmoother.mToInvNoiseGainSquare(currentSmoothR)*eStep;
        double currentSmoothE = LeedSmoother.mToInvNoiseGainSquare(currentSmoothR);
        double maxSmoothE = LeedSmoother.mToInvNoiseGainSquare(maxSmooth);
        double stepForR = LeedUtils.getEnergyStep(unsmoothedForR.energies); //energy step of data for R factor calculation
        int maxSmoothForR = LeedSmoother.invNoiseGainSquareToM(maxSmoothE/stepForR);
        double currentGroupNoise = Double.NaN;
        double[] noises=null;
        if (noiseDependentSmooth) {
            currentGroupNoise = estimateNoise(group, outputRange[currentColumns[0]][0], outputRange[currentColumns[0]][1]);
            if (!(currentGroupNoise > 0)) {
                IJ.error(LEED_Curve_Editor.PLUGIN_NAME, "Error: Cannot determine noise for current beam group");
                IJ.showProgress(1.0);
                return;
            }
            noises = new double[unsmoothedForR.data.length];
            for (int ic=0; ic<noises.length; ic++) {
                int spotIndex = unsmoothedForR.spotIndices[ic];
                int icOrig = LeedUtils.arrayIndexOf(unsmoothed.spotIndices, spotIndex); //column number in original data
                if (icOrig >= 0) {                             
                    noises[ic] = estimateNoise(unsmoothed.data[icOrig],
                            0, unsmoothed.data[icOrig].length);     //data are limited to range already, no need to limit it here
                    if (!(noises[ic] > 0))
                        LeedUtils.logError("Best smoothing internal error: noise="+LeedUtils.d2s(noises[ic])+" for "+unsmoothedForR.spotNames[ic]);
                } else
                    LeedUtils.logError("Best smoothing internal error: missing data for "+unsmoothedForR.spotNames[ic]);
            }
        }
        /* Try all smoothing parameters */
        double[] eVSmooth = new double[maxSmoothForR+1];
        double[] rFactor  = new double[maxSmoothForR+1];
        double bestRFactor = Double.MAX_VALUE;
        int bestRadius = 0;
        for (int smoothRadius = 0; smoothRadius <= maxSmoothForR; smoothRadius++) {
            IJ.showProgress(smoothRadius, maxSmoothForR);
            eVSmooth[smoothRadius] = LeedSmoother.mToInvNoiseGainSquare(smoothRadius)*eStep;
            LeedIVData smoothedForR = unsmoothedForR;
            if (smoothRadius > 0) {
                smoothedForR = smoothedForR.duplicateData();
                for (int ic=0; ic<smoothedForR.data.length; ic++) {
                    int smoothForR = noiseDependentSmooth && noises[ic] > 0 ?
                            smoothRadiusWithNoise(eVSmooth[smoothRadius], currentGroupNoise, noises[ic]) :
                            LeedSmoother.invNoiseGainSquareToM(eVSmooth[smoothRadius]/stepForR);
                    if (smoothForR > maxSmoothForR)
                        smoothForR = maxSmoothForR;
                    smoothedForR.data[ic] = (new LeedSmoother(smoothForR)).smooth(smoothedForR.data[ic]);
                }
            }
            rFactor[smoothRadius] = LEED_R_Factor_Between_Datasets.getRFactor(
                    theoDataForR, smoothedForR, v0i, rFactorType, /*allowShift=*/true, /*showError=*/true);
            if (Double.isNaN(rFactor[smoothRadius])) return;
            if (IJ.debugMode) IJ.log("R="+rFactor[smoothRadius]+" @ smoothing "+eVSmooth[smoothRadius]);
            if (rFactor[smoothRadius] < bestRFactor) {
                bestRFactor = rFactor[smoothRadius];
                bestRadius = smoothRadius;
            }
        }
        /* It is a shallow minimum. We take the lowest smoothing where the R factor is not more than 1.01 times the best one. */
        double eVBestSmooth = 0;
        for (int i=0; i<rFactor.length; i++)
            if (rFactor[i] < 1.01*bestRFactor) {
                eVBestSmooth = eVSmooth[i];
                bestRFactor = rFactor[i];
                break;
            }
        IJ.showProgress(1.0);
        /* Apply the smoothing */
        for (int g=firstGroup; g<=lastGroup; g++) {
            int[] columns = getColumns(g);
            if (columns == null) continue;
            if (useOutput[columns[0]])
                smoothR[columns[0]] = noiseDependentSmooth && eVBestSmooth>0 ?
                        smoothRadiusWithNoise(eVBestSmooth, currentGroupNoise, estimateNoise(g, outputRange[columns[0]][0], outputRange[columns[0]][1])) :
                        LeedSmoother.invNoiseGainSquareToM(eVBestSmooth);
        }
        makeAndShowPlot(true); //update the plot
        showStatus("Best Smoothing: "+IJ.d2s(eVBestSmooth, 1)+", R="+IJ.d2s(bestRFactor));
    }

    /** Returns the noise-dependent smooth radius for a group with newNoise,
     *  if the reference group (usually the current group) with noise 'currentNoise'
     *  has a smooth paarmeter of 'eVcurrentSmooth'. */
    int smoothRadiusWithNoise(double eVcurrentSmooth, double currentNoise, double newNoise) {
        if (!(eVcurrentSmooth > 0)) return 0;
        double eVNewSmooth = eVcurrentSmooth*Math.pow(newNoise/currentNoise, NOISE_SMOOTH_EXPONENT);
        return LeedSmoother.invNoiseGainSquareToM(eVNewSmooth/eStep);
    }

    /** Increases (+1) or decreases (-1) the smooth level, larger increase with shift down.
     *  Sets useOutput to true (switches on the curve) if currently false. */
    void changeSmooth(int delta, InputEvent e) {
        if (e != null && (e.getModifiers() & Event.SHIFT_MASK)!=0) delta *=5;
        int col = currentColumns[0];
        smoothR[col] += delta;
        if (smoothR[col] < MIN_SMOOTH-1) smoothR[col] = MIN_SMOOTH-1;
        if (smoothR[col] > maxSmooth) smoothR[col] = maxSmooth;
        if (!useOutput[col])
            useOutput[col] = true;
        makeAndShowPlot(true);
        showButtonsAndComment();
        showStatus(delta > 0 ? SMOOTH_MORE : SMOOTH_LESS);
        changesDone=true;
    }

    //sets useInput or useOutput
    void setUse(boolean[] useWhat, int col, boolean b) {
        int[] columns = currentColumns;
        if (useWhat == useInput && b && LeedUtils.countTrue(useInput, columns)==0)
            useOutput[columns[0]] = true;         //turning on one curve turns the beam group on
        useWhat[col] = b;
        if (useWhat == useInput && !b && LeedUtils.countTrue(useInput, columns)==0)
            useOutput[columns[0]] = false;        //turning off all curves turns the beam group off
        if (useWhat == useOutput && b && LeedUtils.countTrue(useInput, columns)==0)
            for (int c : columns)
                useInput[c] = true;
        makeAndShowPlot(true);
        showButtonsAndComment();
        selectedGroups = getNGroups(false);
        selectedRange = getSelectedEnergyRange();
        changesDone = true;
    }

    /** Selects or deselects all groups. Asks for confirmation unless currently everything or nothing is selected. */
    void selectAllGroups(boolean select) {
        if (selectedGroups >= 1 && selectedGroups < totalGroups) {
            String what = select ? "Select all" : "Deselect all";
            boolean ok = LeedUtils.yesCancelDialog(PLUGIN_NAME, getPlotWindow(),
                    "Currently "+selectedGroups+"/"+totalGroups+" selected.\n"+what+"?\n"+
                    "You cannot undo this!",
                    what);
            if (!ok) return;
        }
        for (int g = lowestGroup; g <= highestGroup; g++) {
            int[] columns = getColumns(g);
            if (columns == null) continue;
            if (select) {
                double span = getSpan(columns);
                if (span >= minSpanEv - 0.1*eStep)
                    useOutput[columns[0]] = true;   //select only if sufficient energy span
            } else
                useOutput[columns[0]] = false;      //deselect
        }
        makeAndShowPlot(true);
        showButtonsAndComment();
        selectedGroups = getNGroups(false);
        selectedRange = getSelectedEnergyRange();
    }

    /** Sets the output energy range for the given column (first of group)
     *  according to the current roi boundaries (only left&right boundaries are relevant). */
    void selectRangeFromRoi(int column, boolean snapToLimits) {
        selectEnergyRange(column, getRoiEnergyRange(snapToLimits));
    }

    /** Sets the output energy range for the given column (first of group)
     *  according to the energy range.
     *  Call with null for the default range for output (full range for input).
     *  If the group was not selected yet, selects the group (if successful). */
    void selectEnergyRange(int column, double[] eMinMax) {
        int[] oldRange = outputRange[column];
        if (eMinMax == null)
            outputRange[column] = (int[])defaultOutputRange.clone();
        else
            outputRange[column] = energyRangeToIndices(eMinMax[0], eMinMax[1]);
        double span = getSpan(currentColumns);
        if (!askUseSpan(span, "Set range anyhow?", "Set this range"))
            outputRange[column] = oldRange;
        else
            if (!useOutput[currentColumns[0]]) setUse(useOutput, currentColumns[0], true);

        makeAndShowPlot(true);
        changesDone=true;
        selectedRange = getSelectedEnergyRange();
    }

    /** Excludes the energy range of a single beam (out of the group of symmetry-equivalent beams) */
    void excludeEnergyRangeFromRoi(int column) {
        double[] eMinMax = getRoiEnergyRange(/*snapToLimits=*/false);
        int[] exclusionRange = energyRangeToIndices(eMinMax[0], eMinMax[1]);
        int[] oldExclusionRange = noInputRange[column];
        double[] dataWithExclusion = getIntData(column, exclusionRange);
        if (Arrays.equals(intData[column], dataWithExclusion)) {
            IJ.error(LEED_Curve_Editor.PLUGIN_NAME, "Error: No energy range of beam ("+getSpotName(column)+") selected");
            return;
        } else if (LeedUtils.countNonNaN(dataWithExclusion) < 2*v0i/eStep) {
            IJ.error(LEED_Curve_Editor.PLUGIN_NAME, "Remaining energy range of beam ("+getSpotName(column)+
                    ") is too short or zero.\nBeam deselected.");
            setUse(useInput, column, false);
            double span = getSpan(currentColumns);
            if (span < minSpanEv)
                setUse(useOutput, currentColumns[0], false);
            return;
        } else {
            noInputRange[column] = exclusionRange;
            double span = getSpan(currentColumns);
            if (!askUseSpan(span, "Exclude range anyhow?", "Exclude"))
                noInputRange[column] = oldExclusionRange;
            setUse(useInput, column, true); //also updates plot and statistics
            changesDone=true;
        }
    }

    /** Returns the energy index range corresponding to the energy range [eMin, eMax[.
     *  The resulting indices are always within 0 <= index0 < index1 <= energies.length */
    int[] energyRangeToIndices(double eMin, double eMax) {
        double eStart = energies[0];
        double eLast = energies[energies.length-1];
        int iStart = (int)((eMin-eStart)/(eLast-eStart)*(energies.length-1)+0.5);
        int iEnd   = (int)((eMax-eStart)/(eLast-eStart)*(energies.length-1)+0.5)+1;
        iStart = limitInt(iStart, 0, energies.length-1);
        iEnd = limitInt(iEnd, iStart+1, energies.length);
        return new int[] {iStart, iEnd};
    }

    /** Limits an integer to the range >=start and >=0 and <=end */
    int limitInt(int i, int start, int end) {
        if (i < start) i = start;
        if (i > end) i = end;
        if (i < 0) i = 0; //just  to make sure, in case start or end are < 0
        return i;
    }

    /** Returns the energy range corresponding to the current roi, or null if no area roi or error.
     *  If 'snapToLimits' is true, restricts the range to the range limits and snaps to nearby
     *  boundaries of input data to avoid fade in/ fade out. */
    double[] getRoiEnergyRange(boolean snapToLimits) {
        Roi roi = plotImp.getRoi();
        if (roi==null || !roi.isArea()) return null;
        Rectangle bounds = roi.getBounds();
        Calibration cal = plotImp.getCalibration();
        double[] eMinMax = new double[] {
                cal.getX(bounds.x), cal.getX(bounds.x+bounds.width)};
        if (snapToLimits) {
            int[] columns = this.currentColumns;
            double eVsnapRange = SNAP_RANGE*cal.pixelWidth;  //in eV
            if (eVsnapRange < 0.5*v0i) eVsnapRange = 0.5*v0i;
            
            int[] range = energyRangeToIndices(eMinMax[0], eMinMax[1]);
            int[] newRange = snapRange(range, columns, (int)Math.round(eVsnapRange/eStep));
            boolean hasSnapped = false;
            if (newRange[0] != range[0]) {
                eMinMax[0] = energies[newRange[0]];
                hasSnapped = true;
            }
            if (newRange[1] != range[1]) {
                eMinMax[1] = energies[newRange[1]-1];
                hasSnapped = true;
            }
            if (hasSnapped) {       //update Roi
                Rectangle2D.Double rect = roi.getFloatBounds();
                rect.x = Math.round(cal.getRawX(eMinMax[0]));
                rect.width = Math.round(cal.getRawX(eMinMax[1])) - rect.x;
                roi.setBounds(rect);
                plotImp.draw();
            }
        }
        return eMinMax;
    }

    /** Returns the input range or a modified copy if one can slightly (by 'snapRange' or less)
     *  shrink the range to fit one of the limits of the input data for a column.
     *  Shrinking to such a limit avoids fading in or fading out the data for that column when averaging. */
    int[] snapRange(int[] range, int[] columns, int snapRange) {
        LeedIntegerArray[] allLimits = new LeedIntegerArray[] {new LeedIntegerArray(), new LeedIntegerArray()};
        for (int ic=0; ic<columns.length; ic++) {
            if (!useInput[columns[ic]]) continue;
            int[] rangeLimits = LeedUtils.getRangeLimits(intData[columns[ic]], 0, -1);
            for (int j=0; j<rangeLimits.length/2; j++) {
                if (rangeLimits[2*j] < 0) continue;
                if (rangeLimits[2*j] != 0)
                    allLimits[0].add(rangeLimits[2*j]);
                if (rangeLimits[2*j+1] != energies.length)
                    allLimits[1].add(rangeLimits[2*j+1]);
            }
        }
        int[] newRange = range;
        boolean hasSnapped = false;
        for (int n=0; n<2; n++) {   //for lower (left) and upper (right) limit
            int bestLimit = n==0 ? -1 : energies.length;
            double bestDeltaE = Double.MAX_VALUE;
            for (int i=0; i<allLimits[n].size(); i++) {
                int iLimit = allLimits[n].get(i);
                if (n==0) {         //left limits: only snap to right
                    if (range[n]<iLimit && iLimit - range[n] <= snapRange &&
                            iLimit > bestLimit) //for left limits, take the rightmost, to avoid fade-in
                        bestLimit = iLimit;
                } else {            //if n==1, right limits: only snap to left
                    if (range[n]>iLimit && range[n] - iLimit <= snapRange &&
                            iLimit < bestLimit) //for right limits, take the leftmost, to avoid fade-out
                        bestLimit = iLimit;
                }
            }
            if (bestLimit >= 0 && bestLimit < energies.length) {
                if (!hasSnapped) newRange = (int[])range.clone();
                newRange[n] = bestLimit;
                hasSnapped = true;
            }
        }
        return newRange;
    }

    /** Returns the energy span that would be selected when switching on the group
     *  with the given columns. */
    double getSpan(int[] columns) {
        int[] currentRange = outputRange[columns[0]];
        if (currentRange == null) currentRange = new int[] {0, energies.length};
        double[] avg = calculateAverage(columns, currentRange[0], currentRange[1]);
        if (avg == null) return 0;

        avg = LeedUtils.restrictRange(avg, currentRange[0], currentRange[1]);
        avg = avoidHoles(avg, currentRange[0], currentRange[1]);
        int[] rangeLimits = LeedUtils.getRangeLimits(avg, currentRange[0], currentRange[1]);
        if (rangeLimits[0] < 0) return 0;
        double span = energies[rangeLimits[1] - 1] - energies[rangeLimits[0]];
        return span;
    }

    /** Saves the I(V) curves */
    void saveCurves(boolean smooth) {
        saveEdit();
        if (selectedGroups < 1) {
            IJ.error(PLUGIN_NAME, "Error: No curves (Groups) selected, nothing to save");
            return;
        }
        File prefixFile = new File(pathPrefix);
        String oldDir = OpenDialog.getDefaultDirectory();
        String postfix = smooth ? "_edited.csv" : "_ed_nosmooth.csv";
        SaveDialog sd = new SaveDialog("Save Edited I(V) Curves", prefixFile.getParent(), prefixFile.getName()+postfix, ".csv");
        String directory = sd.getDirectory();
        String filename = sd.getFileName();
        OpenDialog.setDefaultDirectory(oldDir);
        if (filename == null) return;
        LeedIVData outData = makeIVData(filename, smooth, lowestGroup, highestGroup);
        if (outData == null)
            IJ.error(LEED_Curve_Editor.PLUGIN_NAME, "Error: No valid data, cannot save");
        else
            outData.save(directory+filename);
    }

    /** Returns a LeedIVData object with the smoothed or unsmoothed output
     *  averaged for the group as selected in the current edit session.
     *  The range of groups is restricted to firstGroup-lastGroup (inclusive).
     *  Returns null if there are no valid data. */
    LeedIVData makeIVData(String title, boolean smooth, int firstGroup, int lastGroup) {
        ArrayList<String> groupNames = new ArrayList<String>();
        ArrayList<double[]> groupData = new ArrayList<double[]>();
        LeedIntegerArray groupSpotIndices = new LeedIntegerArray();
        for (int g=firstGroup; g<=lastGroup; g++) {
            int[] columns = getColumns(g);
            if (columns == null || !useOutput[columns[0]]) continue;
            int iEStart = outputRange[columns[0]][0];
            int iEEnd = outputRange[columns[0]][1];
            double[] data = calculateAverage(columns, iEStart, iEEnd);
            if (data == null) continue;
            if (smooth) {
                data = smooth(data, smoothR[columns[0]]);
                data = avoidNegatives(data, iEStart, iEEnd);
            }
            data = LeedUtils.restrictRange(data, iEStart, iEEnd);
            data = avoidHoles(data, iEStart, iEEnd);
            if (LeedUtils.countNonNaN(data) == 0) continue;
            int spotIndex = spotIndices[columns[0]];
            groupNames.add(spotPattern.getNameWithGroup(spotIndex, true));
            groupData.add(data);
            groupSpotIndices.add(spotIndex);
        }
        int n = groupData.size();
        LeedIVData outData = new LeedIVData(title, energies, groupData.toArray(new double[n][]),
                groupNames.toArray(new String[n]), groupSpotIndices.toArray());
        boolean ok = outData.trimEnergies();
        if (!ok) return null;
        outData.trimData();
        return outData;
    }

    /** Saves the edit parameters (data ranges, smooth points, ...) */
    void saveEdit() {
        if (!changesDone) return;
        File dir = (new File(pathPrefix)).getParentFile();
        while (!dir.isDirectory() || !dir.canWrite()) {
            boolean askDir = LeedUtils.yesCancelDialog(PLUGIN_NAME, getPlotWindow(),
                    "Directory for edit file does not exist or is write-protected:\n"+dir.getPath()+
                    "Select new directory?\n"+
                    "If you click 'Cancel', changes will be lost.",
                    "Select...");
            if (!askDir)
                return;         //user pressed cancel
            DirectoryChooser dc = new DirectoryChooser("Select directory for edit file");
            String newDir = dc.getDirectory();
            if (newDir == null)
                return;         //user pressed cancel
            String namePrefix = (new File(pathPrefix)).getName();
            pathPrefix = newDir + namePrefix;
        }
        String path = pathPrefix + "_Edit_" + LeedUtils.getDateFormatted("yyyyMMdd_HHmm")+".csv";
        PrintWriter pw = null;
        try {
            FileOutputStream fos = new FileOutputStream(path);
            BufferedOutputStream bos = new BufferedOutputStream(fos);
            pw = new PrintWriter(bos);

            pw.println(EDIT_FILE_COLUMNS); // spot,useIn,useOut,smooth,start,end,inStart,inEnd,comment
            for (int col=0; col<spotIndices.length; col++) {
                int spot = spotIndices[col];
                pw.print("(");
                pw.print(spotPattern.getName(spot).replace(',', '|'));
                pw.print("),");
                pw.print(useInput[col]  ? "1," : "0,");
                pw.print(useOutput[col] ? "1," : "0,");
                pw.print(smoothR[col]);
                pw.print(',');
                if (outputRange[col] != null) {
                    pw.print(outputRange[col][0]);
                    pw.print(',');
                    pw.print(outputRange[col][1]);
                } else
                    pw.print("-1,-1");
                pw.print(',');
                if (noInputRange[col] != null) {
                    pw.print(noInputRange[col][0]);
                    pw.print(',');
                    pw.print(noInputRange[col][1]);
                } else
                    pw.print("-1,-1");
                pw.print(',');
                if (comment[col] != null && comment[col].trim().length() >= 0)
                    pw.print(LeedUtils.toCsvString(comment[col]));
                pw.println();
            }
            pw.close();
            changesDone = false;
        } catch (Exception e) {
            IJ.error("Error writing edit log\n"+path+"\n"+e);
            if (pw != null)
                try {pw.close();} catch (Exception e2) {}
        }
    }

    /** Reads edit parameters from a csv file with columns EDIT_FILE_COLUMNS:
     *  spot,useIn,useOut,smooth,start,end,inStart,inEnd,comment */
    @SuppressWarnings("deprecation")
    void readEdit(String path) {
        ResultsTable rt = ResultsTable.open2(path); //IJ.error on failure
        if (rt == null) return;
        String[] colHeadings = EDIT_FILE_COLUMNS.split(",");
        int[] rtColumns = new int[colHeadings.length];
        for (int i=0; i<rtColumns.length; i++) {
            rtColumns[i] = rt.getColumnIndex(colHeadings[i]);
            if (i < 6 && rtColumns[i] < 0) {
                IJ.error("Error in Edit log file", "Column '"+colHeadings[i]+"' missing in\n"+path);
                return;
            }
        }
        for (int i=0; i<rt.size(); i++) {
            String name = rt.getStringValue(rtColumns[0], i);
            int spot = spotPattern.getIndex(name);
            if (spot < 0) continue;
            int column = LeedUtils.arrayIndexOf(spotIndices, spot);
            if (column < 0) continue; //we have read parameters for a beam that is currently not present
            useInput[column]  = rt.getValue(rtColumns[1], i) > 0;
            useOutput[column] = rt.getValue(rtColumns[2], i) > 0;
            smoothR[column] = (int)rt.getValue(rtColumns[3], i);
            int start = (int)rt.getValue(rtColumns[4], i);
            int end =   (int)rt.getValue(rtColumns[5], i);
            if (end > energies.length) end = energies.length;
            outputRange[column] = start>=0 ? new int[]{start, end} : (int[])defaultOutputRange.clone();
            if (rtColumns[6] >= 0) {
                start = (int)rt.getValue(rtColumns[6], i);
                end =   (int)rt.getValue(rtColumns[7], i);
                if (end > energies.length) end = energies.length;
                noInputRange[column] = start>=0 ? new int[]{start, end} : null;
            }
            if (rtColumns[8] >= 0) {
                String str = rt.getStringValue(rtColumns[8], i);
                if ("NaN".equals(str)) str = null;
                if (str != null && str.length() > 0)
                    comment[column] = LeedUtils.fromCsvString(str);
            }
        }
    }

    /** Dialog for setting smoothing 'S', best smoothing 'B',  energy min&max 'E',
     *  min energy span 'M', automatic Selection 'A' */
    void showDialog(final char task) {
        double eVcurrentSmooth = LeedSmoother.mToInvNoiseGainSquare(smoothR[currentColumns[0]])*eStep;
        final boolean isFullSpotPattern = spotPattern.isFullSpotPattern();

        String title = "Set (Default) Smoothing";
        if (task == 'B')
            title = "Set Best Smoothing [vs. Simulated I(V)]";
        if (task == 'E')
            title = "Set (Default) Lowest & Highest Energy";
        else if (task == 'M')
            title = "Set Minimum Energy Span per Curve";
        else if (task == 'A')
            title = "Automatic Selection";
        GenericDialog gd = new GenericDialog(title);
        if (task == 'S') {
            gd.addNumericField("Energy range to smooth", eVcurrentSmooth, 1, 6,
                    "eV (typ. "+IJ.d2s(0.6*v0i, E_DECIMALS)+" to "+IJ.d2s(1.3*v0i, E_DECIMALS)+
                    "; current default "+IJ.d2s(smoothEv, E_DECIMALS)+")");
        }
        if (task == 'S' || task == 'B') {
            gd.addCheckbox("Noise-dependent smoothing *", noiseDependentSmooth);
        }
        if (task == 'B') {
            gd.addFileField("Calculated I(V) csv file (theobeams.csv)", theobeamsPath, 35);
            if (!isFullSpotPattern)     //R factor comparison for "best smoothing" needs a full-fledged SpotPattern file
                gd.addFileField("Spot Pattern file", LeedParams.getString(LeedParams.PATTERNFILE), 35);
        } else if (task == 'E') {
            gd.addNumericField("Lowest energy to use", curveStartEv, 0, 6, "eV *");
            gd.addNumericField("Highest energy to use", curveEndEv, 0, 6, "eV *");
        } else if (task == 'M') {
            gd.addNumericField("Min. data span per curve",
                    minSpanEv, 0, 6, "eV (min. "+IJ.d2s(2*v0i,E_DECIMALS)+")");
        } else if (task == 'A') {
            gd.addNumericField("Noise limit",
                    noiseLimit, 3, 6, "(typ. 0.02-0.08)");
        }
        if (task == 'S' || task == 'E') {
            String[] items = new String[]
                    {" Set as default (for groups not selected)", "Apply to this and all earlier selected groups",
                            "Apply to this and all higher selected groups", "Apply to all selected groups"};
            gd.addRadioButtonGroup("Apply", items, items.length, 1, items[0]);
        }
        if (task == 'E' || task == 'M') {
            String deselectStr = "Deselect all selected groups with less than the minimum span";
            if (task == 'E')
                deselectStr += " ("+IJ.d2s(minSpanEv, E_DECIMALS)+")";
            gd.addCheckbox(deselectStr, false);
        }
        final Checkbox deselectCbx = (task == 'E' || task == 'M') ?
                (Checkbox)(gd.getCheckboxes().lastElement()) : null;

        final String[] autoItems = new String[]{"this group", " this and all earlier groups",
                " this and all higher groups", "all"};
        if (task == 'A' || task == 'B')
            gd.addRadioButtonGroup("Apply for", autoItems, autoItems.length, 1, autoItems[0]);

        if (task == 'E')
            gd.addMessage("*  must allow energy span of "+IJ.d2s(minSpanEv, E_DECIMALS)+
                    " in "+IJ.d2s(energies[0], E_DECIMALS)+"-"+
                    IJ.d2s(energies[energies.length-1], E_DECIMALS)+" range\n"+
                    " 'Apply' restricts the range, does not widen it");
        if (task == 'S')
            gd.addMessage("*  Stronger smoothing for curves with higher noise (with 'Apply' to several groups).\n"+
                    "   Best select smoothing for a curve with medium start energy and medium noise, apply for all,\n"+
                    "   then correct smoothing of low-energy data if required (esp. for 5d elements).\n"+
                    "   Apply noise-dependent smoothing AFTER energy range selection ('Auto select')." );

        DialogListener dialogListener = new DialogListener() {
            public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
                if (task == 'S') {
                    double smoothEv = gd.getNextNumber();
                    if (smoothEv < 0) return false;
                } else if (task == 'B') {
                    String theobeamsPath = gd.getNextString();
                    File theobeamsFile = new File(theobeamsPath);
                    if (!theobeamsFile.isFile() || !theobeamsFile.canRead()) return false;
                    if (!isFullSpotPattern) {
                        String patternPath = gd.getNextString();
                        File patternFile = new File(patternPath);
                        if (!patternFile.isFile() || !patternFile.canRead()) return false;
                    }
                }
                if (task == 'S' || task == 'B') {
                    boolean noiseDependentSmooth = gd.getNextBoolean();
                    String radioString = gd.getNextRadioButton();
                    boolean enableNoiseDependent =      //noiseDependentSmooth checkbox current selected and 'apply' only
                            useOutput[currentColumns[0]] &&
                            (radioString.contains("Apply") || LeedUtils.arrayIndexOf(autoItems, radioString) > 0);
                    if (!useOutput[currentColumns[0]]) enableNoiseDependent = false;
                    ((Checkbox)(gd.getCheckboxes().get(0))).setEnabled(enableNoiseDependent);
                } else if (task == 'E') {
                    double curveStartEv = gd.getNextNumber();
                    double curveEndEv = gd.getNextNumber();
                    if (!(curveStartEv > 0 && curveStartEv < energies[energies.length-1] - minSpanEv &&
                            curveEndEv >= curveStartEv+minSpanEv)) return false;
                    String applyHow = gd.getNextRadioButton();
                    deselectCbx.setEnabled(applyHow.startsWith("Apply")); //"deselect < span" makes sense only if apply to all
                } else if (task == 'M') {
                    double minSpanEv = gd.getNextNumber();
                    if (!(minSpanEv >= 2*v0i)) return false;
                } else if (task == 'A') {
                    double noiseLimit = gd.getNextNumber();
                    return (noiseLimit >= 0.001 && noiseLimit <= 0.3);
                }
                return true;
            }
        };
        Button okButton = gd.getButtons()[0];
        okButton.setEnabled(dialogListener.dialogItemChanged(gd, null));
        gd.addDialogListener(dialogListener);

        gd.showDialog();
        if (gd.wasCanceled()) return;

        if (task == 'S') {
            eVcurrentSmooth = gd.getNextNumber();
        } else if (task == 'B') {
            theobeamsPath = gd.getNextString();
            if (!isFullSpotPattern) {                       // try to use the spot pattern from the file
                String patternPath = gd.getNextString();
                LeedSpotPattern sp = new LeedSpotPattern(patternPath);
                String error = sp.getErrorMessage();
                if (error != null) {
                    IJ.error(LEED_Curve_Editor.PLUGIN_NAME+" Spot Pattern", error);
                    return;
                }
                String[] currentSpotPatternNames = this.spotPattern.getAllSpotNames();
                int[] spIndices = sp.getIndices(currentSpotPatternNames);
                error = sp.getErrorMessage();
                if (error != null) {
                    IJ.error(LEED_Curve_Editor.PLUGIN_NAME+" Spot Pattern does not fit", error);
                    return;
                }
                for (int i=0; i<spIndices.length; i++) {
                    if (sp.getGroup(spIndices[i]) != this.spotPattern.getGroup(i)) {
                        IJ.error(LEED_Curve_Editor.PLUGIN_NAME, "Error: Group number definitions incompatible\n"+
                                "in current data and Spot Pattern file:\n"+
                                this.spotPattern.getNameWithGroup(i, false)+
                                " vs. "+sp.getNameWithGroup(spIndices[i], false));
                        return;
                    }
                }
                for (int i=0; i<this.spotIndices.length; i++)   //translate spotIndices of columns to new spotPattern
                    this.spotIndices[i] = spIndices[this.spotIndices[i]];
                this.spotPattern = sp;                          //now we have the the new spotPattern
            }
        }
        if (task == 'S' || task == 'B') {
            noiseDependentSmooth = gd.getNextBoolean();
        } else if (task == 'E') {
            curveStartEv = gd.getNextNumber();
            curveEndEv = gd.getNextNumber();
            defaultOutputRange = energyRangeToIndices(curveStartEv, curveEndEv);
        } else if (task == 'M') {
            minSpanEv = gd.getNextNumber();
        } else if (task == 'A') {
            noiseLimit = gd.getNextNumber();
        }
        boolean apply = false;
        boolean onlyHighGroups = false;
        boolean onlyLowGroups = false;
        if (task == 'S' || task == 'E') {
            String radioString = gd.getNextRadioButton();
            apply = radioString.contains("Apply");
            onlyLowGroups = radioString.contains("earlier");
            onlyHighGroups = radioString.contains("higher");
        }
        boolean deselect = false;
        if (task == 'E' || task == 'M')
            deselect = gd.getNextBoolean();

        if (task == 'M' && !deselect) return; //just setting default, nothing to do

        if (task == 'A'  || task == 'B') {
            String radioString = gd.getNextRadioButton();
            apply = !radioString.equals(autoItems[0]);
            onlyLowGroups = radioString.contains("earlier");
            onlyHighGroups = radioString.contains("higher");
        }

        double currentGroupNoise = Double.NaN;
        if (task == 'S' && apply && noiseDependentSmooth) {
            currentGroupNoise = estimateNoise(group, outputRange[currentColumns[0]][0], outputRange[currentColumns[0]][1]);
            if (!(currentGroupNoise > 0)) {
                IJ.error(LEED_Curve_Editor.PLUGIN_NAME, "Error: Cannot determine noise for current beam group");
                return;
            }
        }

        if (task == 'S' && !apply)
            smoothEv = eVcurrentSmooth;

        int startGroup = lowestGroup;
        int lastGroup = highestGroup;
        if (onlyLowGroups || ((task == 'A' || task == 'B') && !apply))
            lastGroup = group;
        if (onlyHighGroups || ((task == 'A' || task == 'B') && !apply))
            startGroup = group;
        //IJ.log("apply="+apply+" groups="+startGroup+"-"+lastGroup);
        if (task == 'B') {
            findBestSmoothing(theobeamsPath, noiseDependentSmooth, startGroup, lastGroup);
        } else {
            for (int g = startGroup; g <= lastGroup; g++) {  //apply new settings
                int[] columns = getColumns(g);
                if (columns == null) continue;
                if (task == 'S') {
                    if ((apply && useOutput[columns[0]]) || !useOutput[columns[0]]) {
                        smoothR[columns[0]] = noiseDependentSmooth && currentGroupNoise>0 && eVcurrentSmooth>0 ?
                                smoothRadiusWithNoise(eVcurrentSmooth, currentGroupNoise, estimateNoise(g, outputRange[columns[0]][0], outputRange[columns[0]][1])) :
                                LeedSmoother.invNoiseGainSquareToM(eVcurrentSmooth);
                    } else if (!useOutput[columns[0]])
                        smoothR[columns[0]] = LeedSmoother.invNoiseGainSquareToM(eVcurrentSmooth);
                } else if (task == 'E') {
                    if (apply && useOutput[columns[0]] || !useOutput[columns[0]]) {
                        int[] range = outputRange[columns[0]];
                        if (range[0] < defaultOutputRange[0]) range[0] = defaultOutputRange[0];
                        if (range[1] > defaultOutputRange[1]) range[1] = defaultOutputRange[1];
                    } else if (!useOutput[columns[0]])
                        outputRange[columns[0]] = (int[])(defaultOutputRange.clone());
                } else if (task == 'A')
                        autoSelectGroup(g, noiseLimit);
                if (deselect && useOutput[columns[0]]) {
                    double span = getSpan(columns);
                    if (span < minSpanEv)
                    useOutput[columns[0]] = false;
                }
                calculatePlotData(columns);             //updates selectedRangeLengths
            }
        }
        changesDone = true;
        selectedGroups = getNGroups(false);
        selectedRange = getSelectedEnergyRange();
        makeAndShowPlot(true); //update plot
    }

    /** Shows the options dialog */
    void showOptionsDialog() {
        GenericDialog gd = new GenericDialog(LEED_Curve_Editor.PLUGIN_NAME);
        gd.addNumericField("V0i", this.v0i, 1, 6, " eV (min. "+IJ.d2s(eStep+0.0499999,1)+")");
        gd.addCheckbox("Fixed x axis (full energy range)", this.fixedXrange);
        gd.addCheckbox("Show Spot Tracker energy as vertical line", showSpotTrackerEnergy);
        gd.addCheckbox("Highlight current beams in Spot Tracker", highlightInSpotTracker);
        gd.addCheckbox("Plot estimated noise (blue curve)", plotNoise);
        DialogListener dialogListener = new DialogListener() {  //only checks lower limit of numeric inputs
            public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
                double v0i = Math.abs(gd.getNextNumber());
                if (!(v0i >= eStep)) return false;
                return true;
            }
        };
        Button okButton = gd.getButtons()[0];
        okButton.setEnabled(dialogListener.dialogItemChanged(gd, null));
        gd.addDialogListener(dialogListener);
        gd.showDialog();

        if (gd.wasCanceled()) return;
        v0i = Math.abs(gd.getNextNumber());
        if (v0i < eStep) v0i = eStep;
        if (minSpanEv < 2*v0i) {
            if (IJ.showMessageWithCancel(LEED_Curve_Editor.PLUGIN_NAME,
                    "Minimum data span per curve will be set to 2*V0i = "+IJ.d2s(2*v0i, 1)+" eV")) {
                minSpanEv = 2*v0i;
                this.v0i = v0i;
            } else
                return;
        }
        boolean fixedXrange = gd.getNextBoolean();
        groupChanged = fixedXrange!=this.fixedXrange;    //enforces new x axis
        showSpotTrackerEnergy = gd.getNextBoolean();
        if (showSpotTrackerEnergy) {
            LEED_Spot_Tracker spotTracker = LEED_Spot_Tracker.getInstance();
            if (spotTracker != null)
                spotTrackerEnergy = spotTracker.getEnergy();
        } else
            spotTrackerEnergy = Double.NaN;
        highlightInSpotTracker = gd.getNextBoolean();
        if (highlightInSpotTracker)
            showHighlightInSpotTracker(getSpotName(currentColumns[0]));
        else
            showHighlightInSpotTracker(null);
        this.fixedXrange = fixedXrange;
        plotNoise = gd.getNextBoolean();
        makeAndShowPlot(true);  //update plot
    }

    /** Shows a dialog for the comment */
    void showCommentDialog() {
        String text = comment[currentColumns[0]];
        if (text == null) text = "";
        text = text.trim();
        GenericDialog gd = new GenericDialog(LEED_Curve_Editor.PLUGIN_NAME);
        gd.addMessage("Comment for beam group ["+spotPattern.getGroup(spotIndices[currentColumns[0]])+"]:");
        gd.addTextAreas(text, null, 2, 50);
        gd.showDialog();

        if (gd.wasCanceled()) return;
        String newText = gd.getNextText().trim();
        if (!newText.equals(text))
            changesDone = true;
        comment[currentColumns[0]] = newText.length()==0 ? null : newText;
        if (comment[currentColumns[0]] != null)
            lastComment = newText;
        showButtonsAndComment();
    }

    /** Asks the user for the badness limit for equivalent beams
     *  (a measure of bad agreement to the average).
     *  Returns false if canceled or the user has entered a non-positive number. */
    boolean showBadnessDialog() {
        if (!(badnessLimit > 0))
            badnessLimit = 2*noiseLimit;
        badnessLimit = IJ.getNumber("Find inequivalent beams with a 'local' R factor larger then:", badnessLimit);
        if (badnessLimit > 0)
            badEquivalentBeams(group, badnessLimit);    //examine current group
        makeAndShowPlot(false);
        return badnessLimit > 0;
    }

     /** If the span would be too small, asks the user whether to use and returns false if the user clicks 'no'.
      *  Question should be e.g. "Select anyhow?" and the label of the OK button "Use this Group".
      *  If the span is less than one step, displays an error message and returns false. */
    boolean askUseSpan(double span, String question, String okLabel) {
        if (span < 1.1*eStep) {
            IJ.error(PLUGIN_NAME, "Error: No data in selected energy range");
            return false;
        } else if (span < minSpanEv - 0.1*eStep) {
            return LeedUtils.yesCancelDialog(PLUGIN_NAME, getPlotWindow(),
                    "Energy span ("+IJ.d2s(span,E_DECIMALS)+") is less than minimum ("+
                    IJ.d2s(minSpanEv,E_DECIMALS)+")\n"+question,
                    okLabel);
        }
        return true;
    }

    /** Lists all comments */
    void listComments() {
        ResultsTable rt = new ResultsTable();
        for (int ic=0; ic<comment.length; ic++) {
            if (comment[ic] != null) {
                rt.incrementCounter();
                rt.addLabel(spotPattern.getNameWithGroup(spotIndices[ic], false));
                rt.addValue("Comment", comment[ic]);
            }
        }
        if (rt.size() > 0)
            rt.show("Comments - "+(new File(pathPrefix)).getName());
        else
            IJ.error(LEED_Curve_Editor.PLUGIN_NAME, "No Comments");
    }

    /** Lists beams with bad agreement (selected beams only) */
    void listBadBeams() {
        if (LeedUtils.countTrue(useOutput) == 0) {
            IJ.error(LEED_Curve_Editor.PLUGIN_NAME, "No beams selected (only selected beams will be listed)");
            return;
        }
        if (!(badnessLimit > 0)) {                       //we have to ask for the badness value
            if (!showBadnessDialog()) return;
        }
        ResultsTable rt = new ResultsTable();
        for (int g = lowestGroup; g <= highestGroup; g++) {
            if (badEquivalentBeams(g, badnessLimit)) {
                for (int i=0; i<badBeamColumns.size(); i++) {
                    rt.incrementCounter();
                    rt.addValue("Group", g);
                    rt.addValue("Beam", getSpotName(badBeamColumns.get(i)));
                    rt.addValue("Energy", IJ.d2s(badEnergies.get(i),0));
                    rt.addValue("Badness", IJ.d2s(badnessValues.get(i)));
                }
            }
        }
        if (rt.size() > 0)
            rt.show("Bad agreement (>"+IJ.d2s(badnessLimit)+") - "+(new File(pathPrefix)).getName());
        else
            IJ.error(LEED_Curve_Editor.PLUGIN_NAME, "No equivalent beams with bad agreement selected");
    }

   /** Returns the current PlotWindow */
    PlotWindow getPlotWindow() {
        ImagePlus imp = plotImp;
        return imp == null || !imp.isVisible() ? null : (PlotWindow)imp.getWindow();
    }

    ImageCanvas getCanvas() {
        PlotWindow win = getPlotWindow();
        return win == null ? null : win.getCanvas();
    }

    /** Shows the buttons and the comment as overlays */
    void showButtonsAndComment() {
        PlotWindow pw = getPlotWindow();
        int group = this.group;
        int[] columns = currentColumns;
        Color[] dataColors = currentColors;
        buttonEnabled[PREV] = group > lowestGroup;
        buttonEnabled[NEXT] = group < highestGroup;
        buttonEnabled[SMOOTH_LESS] = smoothR[columns[0]] > MIN_SMOOTH-1;
        buttonEnabled[SMOOTH_MORE] = smoothR[columns[0]] < maxSmooth;
        buttonActive[ON_OFF] = useOutput[columns[0]];
        Roi.setDefaultStrokeWidth(1);
        for (int ic=0; ic<columns.length; ic++)
            buttonActive[N_BUTTONS+ic] = useInput[columns[ic]];

        int width = plotImp.getWidth();
        Overlay ovly = new Overlay();
        for (int b=0; b<N_BUTTONS; b++)                 //function buttons
            addButton(b, ovly, width, null, spotIndices[columns[0]]);
        for (int c=0; c<columns.length; c++)            //data curve buttons
            addButton(N_BUTTONS+c, ovly, width, dataColors[c], spotIndices[columns[c]]);
        for (int i=N_BUTTONS+columns.length; i<buttonRects.length; i++)
            buttonRects[i] = null;
        String text = comment[columns[0]];
        if (text != null) {
            Plot plot = pw.getPlot();
            Rectangle plotFrame = plot.getDrawingFrame();
            int yText = (int)Math.round(plotFrame.y + 0.25*plotFrame.height + BUTTON_FONT.getSize()*3/2);
            TextRoi textRoi = new TextRoi(text, width - 2*BUTTON_SIZE -5*BUTTON_SPACE, yText, BUTTON_FONT);
            textRoi.setFontSize(14);
            textRoi.setStrokeColor(Color.BLACK);
            textRoi.setJustification(TextRoi.RIGHT);
            ovly.add(textRoi);
        }
        plotImp.setOverlay(ovly);
    }

    /** Adds a function button, or a data button if b >= N_BUTTONS */
    void addButton(int b, Overlay ovly, int width, Color dataColor, int spotIndex) {
        int ix = b < N_BUTTONS ? BUTTON_X[b] : -1;      //grid position of button
        int iy = b < N_BUTTONS ? BUTTON_Y[b] : BUTTON_Y[N_BUTTONS-1] + 1 + (b - N_BUTTONS);
        int x = width + (ix - 1) * (BUTTON_SIZE + BUTTON_SPACE);
        int y = BUTTON_SPACE + iy * (BUTTON_SIZE + BUTTON_SPACE);

        Color background = buttonActive[b] && buttonEnabled[b] ? C_BG_ACTIVE : C_BG_INACTIVE;
        Roi backg = new Roi(x, y, BUTTON_SIZE, BUTTON_SIZE);
        backg.setFillColor(background);
        ovly.add(backg);
        Roi frame = new Roi(x, y, BUTTON_SIZE, BUTTON_SIZE);
        frame.setStrokeWidth(1);
        frame.setStrokeColor(Color.BLACK);
        ovly.add(frame);
        buttonRects[b] = frame.getBounds();
        Roi symbol = null, symbol2 = null;


        Color symbolColor = dataColor;
        if (b == ON_OFF)
            symbolColor = buttonActive[b] ? BUTTON_COLOR[b] : Color.BLACK;
        else if (b == OPTIONS)
            symbolColor = Color.DARK_GRAY;
        else if (b < N_BUTTONS)
            symbolColor = buttonEnabled[b] ? BUTTON_COLOR[b] : Color.GRAY;
        if (b==ON_OFF) {                    // ON/OFF is a special case, circle
            symbol2 = new OvalRoi(x+BUTTON_SIZE/6, y+BUTTON_SIZE/6, BUTTON_SIZE*2/3+1, BUTTON_SIZE*2/3+1);
        } else if (b == HELP) {             // HELP button: question mark and small 'T' as polygon
            symbol2 = new TextRoi("?", x+BUTTON_SIZE/6, y+BUTTON_SIZE*3/4, BUTTON_FONT);
        } else if (b==OPTIONS) {            // gearwheel symbol for options dialog
            if (BUTTON_SYMBOL_X[b] == null) {
                BUTTON_SYMBOL_X[b] = new float[3*8+1];
                BUTTON_SYMBOL_Y[b] = new float[3*8+1];
                for (int i=0; i<3*8+1; i++) {
                    double angle = (2*Math.PI/(3*8))*i;
                    double r = (i%3==0) ? 0.28 : 0.38;
                    BUTTON_SYMBOL_X[b][i] = 0.5f + (float)(r*Math.cos(angle));
                    BUTTON_SYMBOL_Y[b][i] = 0.5f + (float)(r*Math.sin(angle));
                }
            }
            symbol2 = new OvalRoi(x+BUTTON_SIZE*3/8, y+BUTTON_SIZE*3/8, BUTTON_SIZE/4, BUTTON_SIZE/4);
        }
        if (b < N_BUTTONS) {
            if (BUTTON_SYMBOL_X[b] != null)
                symbol = new PolygonRoi(toSymbolInt(x, BUTTON_SYMBOL_X[b]), toSymbolInt(y, BUTTON_SYMBOL_Y[b]),
                        BUTTON_SYMBOL_X[b].length, Roi.POLYLINE);
        } else {                            // curve color line
            symbol = new Roi(x+BUTTON_SIZE/8, y+BUTTON_SIZE*5/12, BUTTON_SIZE*3/4+1, BUTTON_SIZE/6+1);
        }
        int strokeWidth = (int)Math.round(BUTTON_SIZE*(buttonActive[b] && b!=ON_OFF && b!=OPTIONS ? 1/6f : 1/12f));
        if (strokeWidth < 1) strokeWidth = 1;
        symbol.setStrokeWidth(strokeWidth);
        symbol.setStrokeColor(symbolColor);
        ovly.add(symbol);
        if (symbol2 != null) {
            symbol2.setStrokeWidth(2);
            symbol2.setStrokeColor(symbolColor);
            ovly.add(symbol2);
        }
        if (b >= N_BUTTONS && buttonActive[b]) {   //add white outline to line
            Rectangle r = symbol.getBounds();
            r.x--; r.y--; r.width++; r.height++;
            Roi whiteRect = new Roi(r);
            whiteRect.setStrokeColor(Color.WHITE);
            whiteRect.setStrokeWidth(BUTTON_SIZE*3/48);
            ovly.add(whiteRect);
        }
        if (buttonActive[b]) {  //for pressed button, add black shadow & white "anti-shadow"
            Roi blackPoly = new PolygonRoi(
                    new int[]{x+1, x+1, x+BUTTON_SIZE-2}, new int[]{y+BUTTON_SIZE-2, y+1, y+1},
                    3, Roi.POLYLINE);
            blackPoly.setStrokeColor(Color.BLACK);
            blackPoly.setStrokeWidth(BUTTON_SIZE*3/48);
            ovly.add(blackPoly);
            Roi whitePoly = new PolygonRoi(
                    new int[]{x+BUTTON_SIZE-1, x+BUTTON_SIZE-1, x+2}, new int[]{y+2, y+BUTTON_SIZE-1, y+BUTTON_SIZE-1},
                    3, Roi.POLYLINE);
            whitePoly.setStrokeColor(Color.WHITE);
            whitePoly.setStrokeWidth(BUTTON_SIZE*3/48);
            ovly.add(whitePoly);
        }
        if (b >= N_BUTTONS) {    //add label to curve buttons
            String name = spotPattern.getName(spotIndex);
            boolean twoLine = name.indexOf('/') > 0; //fractional spots are displayed with two lines
            if (!twoLine)
                y += smallFont.getSize()/2;
            String[] lines = twoLine ? name.split(",") : new String[] {name};
            for (int i=0; i<lines.length; i++) {
                y += smallFont.getSize();
                Roi text = new TextRoi(lines[i], x + BUTTON_SIZE + (int)(4*GUI_SCALE), y, smallFont);
                text.setFillColor(Color.WHITE);
                text.setStrokeColor(Color.BLACK);
                ovly.add(text);
            }
        }
        if (b == ON_OFF) {
            TextRoi text = new TextRoi(x-BUTTON_SIZE/6, y-BUTTON_SIZE/4,
                    "("+spotPattern.getName(spotIndex)+")", normFont);
            text.setJustification(TextRoi.RIGHT);
            text.setStrokeColor(Color.BLACK);
            ovly.add(text);
            text = new TextRoi(x-BUTTON_SIZE/6, y+BUTTON_SIZE/2,
                    "["+spotPattern.getGroup(spotIndex)+"]", normFont);
            text.setJustification(TextRoi.RIGHT);
            text.setStrokeColor(Color.BLACK);
            ovly.add(text);
        }
    }

    /** Converts floating-point symbol polygon coordinates to integer for the button size */
    int[] toSymbolInt(int offset, float[] coord) {
        int[] iCoord = new int[coord.length];
        for (int i=0; i<coord.length; i++)
            iCoord[i] = offset + (int)Math.round(coord[i] * BUTTON_SIZE);
        return iCoord;
    }

    /** Returns the number of the button at pixel position x, y or -1 if none */
    int getButton(int x, int y) {
        for (int i=0; i<buttonRects.length; i++) {
            if (buttonRects[i] == null) break;
            if (buttonRects[i].contains(x, y)) return i;
        }
        return -1;
    }

    /** Goes to the given group and displays the plot. Does nothing if group = Integer.MIN_VALUE.
     *  When interactive, also sets the corresponding group of other open windows */
    void setGroup(int g, boolean interactive) {
        if (g == Integer.MIN_VALUE) return;
        int[] columns = getColumns(g);
        if (columns == null) return;
        if (group != g) groupChanged = true;
        group = g;
        currentColumns = columns;
        int[] range = outputRange[columns[0]];
        if (range ==  null || range[0] < 0 || range[1] <= range[0]) {
            outputRange[columns[0]] = (int[])defaultOutputRange.clone();
            range = outputRange[columns[0]];
        }
        //IJ.log("setGrp="+g+" col="+spotPattern.getName(spotIndices[currentColumns[0]])+" etc, #="+currentColumns.length);
        currentColors = LeedPlotter.getColors(columns.length+1); //orange would be at the end, reserved for avg
        makeAndShowPlot(true);
        showButtonsAndComment();
        mouseMoved(null);               //update status line
        if (interactive)
            synchronizeAllCurvePlotEditors(group);

        IJ.setTool(Toolbar.RECTANGLE);
        if (Arrays.equals(range, defaultOutputRange)) {
            plotImp.killRoi();
        } else {                        //set roi according to current range (not if default range)
            Calibration cal = plotImp.getCalibration();
            int left  = (int)Math.round(cal.getRawX(energies[Math.min(range[0], energies.length-1)]));
            int right = (int)Math.round(cal.getRawX(energies[Math.min(range[1], energies.length) - 1]));
            if (left < 2) left = 2;
            if (right > plotImp.getWidth()) right = plotImp.getWidth() - 8;
            Roi roi = new Roi(left, 2, right-left, plotImp.getHeight()-4);
            plotImp.setRoi(roi);
        }
        if (highlightInSpotTracker)
            showHighlightInSpotTracker(getSpotName(columns[0]));
    }

    Object spotTrackerHighlightSynchronizer = new Object();
    /** Highlights the given spot and its group in the Spot Tracker.
     *  This is done in a separate thread since it can be slow if the overlay of the SpotTracker contains
     *  many items (e.g. 1000 slices and 2000 spots) */
    void showHighlightInSpotTracker(final String spotName) {
        if (LEED_Spot_Tracker.getInstance() == null) return;
        spotTrackerHighlightName = spotName;
        synchronized(spotTrackerHighlightSynchronizer) {
            if (spotTrackerHighlightThread == null || !spotTrackerHighlightThread.isAlive()) {
                spotTrackerHighlightThread = new Thread(new Runnable() {public void run() {
                        String lastHighlightName = null;
                        while (true) {
                            LEED_Spot_Tracker spotTracker = LEED_Spot_Tracker.getInstance();
                            PlotWindow pw = getPlotWindow();
                            if (spotTracker == null || pw == null || !pw.isVisible()) {
                                synchronized(spotTrackerHighlightSynchronizer) {
                                    spotTrackerHighlightThread = null;
                                    return;
                                }
                            }
                            synchronized (spotTrackerHighlightSynchronizer) {
                                if (spotTrackerHighlightName == lastHighlightName)
                                    try {
                                        spotTrackerHighlightSynchronizer.wait(5000);     //timeout, so we can eventually quit when everything is closed
                                    } catch (InterruptedException e) {}
                            }
                            String name = spotTrackerHighlightName;
                            if (name != lastHighlightName) {
                                spotTracker = LEED_Spot_Tracker.getInstance();
                                if (spotTracker != null)
                                    spotTracker.highlightSpotWithGroup(name);
                                lastHighlightName = name;
                            }
                        }
                    }
                }, "setSpotTrackerHighlight"); //name of 'new Thread'
                spotTrackerHighlightThread.start();
            }
        }
        synchronized(spotTrackerHighlightSynchronizer) {
            spotTrackerHighlightSynchronizer.notify();
        }
    }

    /** Starting from the current group (exclusive), selects the next or previous group with a comment.
     *  'direction' should be +1 for forward and -1 for backward search. */
    void goCommentedGroup(int direction) {
        for (int g=group+direction; g>=lowestGroup && g<= highestGroup; g+=direction) {
            int[] columns = getColumns(g);
            if (columns != null && comment[columns[0]] != null) {
                setGroup(g, true);
                return;
            }
        }
        IJ.beep();
        showStatus((direction > 0 ? "No more" : "No earlier") + " group with a comment");
    }

    /** Starting from the current group (exclusive), selects the next or previous selected group with
     *  bad agreement of equivalent beams.
     *  Asks for the R factor limit when called the first time or with alt key down.
     *  'direction' should be +1 for forward and -1 for backward search. */
    void goBadGroup(int direction) {
        int firstStep = direction;
        if (!(badnessLimit > 0) || IJ.altKeyDown()) {   //we have to ask for the badness value
            if (IJ.altKeyDown())
                IJ.setKeyUp(KeyEvent.VK_ALT);           //don't ask again if the user has forgotten to leave the alt key
            if (!showBadnessDialog()) return;
            firstStep = 0;                              //after changing the limit, also examine the current group
        }
        for (int g=group+firstStep; g>=lowestGroup && g<= highestGroup; g+=direction) {
            int[] columns = getColumns(g);
            if (columns != null && useOutput[columns[0]] && badEquivalentBeams(g, badnessLimit)) {
                setGroup(g, true);
                return;
            }
        }
        showStatus((direction > 0 ? "No more" : "No earlier") + " group with bad agreement (>"+IJ.d2s(badnessLimit)+") of inequivalent beams");
    }

    /** Goes to the beam or group given by the group. Returns false if the group is not found */
    boolean gotoBeam(String str, boolean interactive) {
        int group = getGroup(str);
        if (group == Integer.MIN_VALUE) return false;
        setGroup(group, interactive);
        return true;
    }

    /** When several LeedCurvePlotEditor windows with the same spot pattern are open,
     *  synchronizes the group shown and x range of them */
    void synchronizeAllCurvePlotEditors(int group) {
        ArrayList<LeedCurvePlotEditor> list = getAllCurvePlotEditors();
        if (list.size() <= 1) return;
        list.remove(this);
        int currentBeam = spotIndices[currentColumns[0]];
        String currentBeamName = spotPattern.getName(currentBeam);
        double[] plotLimits = getPlotLimits();
        for (int i=list.size()-1; i>=0; i--) { //(the list will be modified by removing elements)
            LeedCurvePlotEditor curvePlotEditor = list.get(i);
            boolean success = curvePlotEditor.gotoBeam(currentBeamName, false);
            if (!success) {         //ignore editor windows that don't have the current beam
                list.remove(curvePlotEditor);
                continue;
            }
            double[] limits = getPlotLimits();
            if (limits != null) {   //find the plot energy range that fits all
                if (limits[0] < plotLimits[0]) plotLimits[0] = limits[0];
                if (limits[1] > plotLimits[1]) plotLimits[1] = limits[1];
            }
        }
        for (LeedCurvePlotEditor curvePlotEditor : list)    //set plot energy range for all
            curvePlotEditor.setPlotEnergyRange(plotLimits[0], plotLimits[1]);
    }

    /** The Spot Tracker sets its current energy with this */
    public void setSpotTrackerEnergy(double energy) {
        spotTrackerEnergy = energy;
        if (showSpotTrackerEnergy)
            makeAndShowPlot(false);
    }

    /** Returns the plot range Emin, Emax, Imin, Imax or null if no plot */
    double[] getPlotLimits() {
        ImagePlus imp = plotImp;        //plotImp might become null asynchronously
        Plot plot = imp == null ? null : (Plot)imp.getProperty(Plot.PROPERTY_KEY);
        if (plot == null) return null;
        return plot.getLimits();
    }

    /** Sets the energy range of the plot; leaves the y range unchanged.
     *  Used when curves get synchronized. */
    void setPlotEnergyRange(double minE, double maxE) {
        ImagePlus imp = plotImp;        //plotImp might become null asynchronously
        Plot plot = imp == null ? null : (Plot)imp.getProperty(Plot.PROPERTY_KEY);
        if (plot == null) return;
        double[] limits = plot.getLimits();
        limits[0] = minE;
        limits[1] = maxE;
        plot.setLimits(limits);
        PlotWindow pw = getPlotWindow();
        if (pw != null) pw.getImagePlus().killRoi(); //Roi does not fit after rescaling the plot
        plot.updateImage();
    }

    /** Returns whether the gui should be blocked*/
    boolean guiBlocked() {
        return guiBlockingThread != null && guiBlockingThread.isAlive();
    }

    /** React to buttons. Also called from keyPressed, with e=null */
    void handleButtonPressed(int button, MouseEvent e) {
        if (guiBlocked()) return;
        boolean isRightClick = LeedUtils.isPopupTrigger(e);
        int modifiers = e==null ? 0 : e.getModifiersEx();
        boolean ctrlOrMeta = (modifiers & (InputEvent.CTRL_DOWN_MASK | InputEvent.META_DOWN_MASK)) != 0;
        if (isRightClick || button == OPTIONS) {
            if (CONTEXT_MENU_ITEMS[Math.min(button, N_BUTTONS)] != null)
                showMenu(button, getCanvas(), e);
        } else { //normal (not right click) or no context menu for this button
            if (!buttonEnabled[button]) return;
            switch(button) {
                case HELP:
                    showHelp();
                    break;
                case SAVE:
                    saveCurves(true);
                    break;
                case PREV:
                    if (ctrlOrMeta)
                        goBadGroup(-1);
                    else
                        setGroup(getPreviousGroup(), true);
                    break;
                case NEXT:
                    if (ctrlOrMeta)
                        goBadGroup( 1);
                    else
                        setGroup(getNextGroup(), true);
                    break;
                case ON_OFF:
                    boolean useOut = useOutput[currentColumns[0]];
                    if (useOut || askUseSpan(getSpan(currentColumns), "Select anyhow?", "Use this Group"))
                        setUse(useOutput, currentColumns[0], !useOut);
                    changesDone = true;
                    break;
                case SMOOTH_LESS:
                    changeSmooth(-1, e);
                    break;
                case SMOOTH_MORE:
                    changeSmooth( 1, e);
                    break;
                case ROI:
                    selectRangeFromRoi(currentColumns[0], true);
                    break;
                //case OPTIONS: already handled above
                default:    //buttons for single curves
                    int i = button - N_BUTTONS;
                    boolean useIn = useInput[currentColumns[i]];
                    setUse(useInput, currentColumns[i], !useIn);    //toggle useInput
                break;
            }
        }
        if (e!=null)    //on mouseover or click
            showStatus(button);
        else            //after button activation via keyboard shortcut
            getPlotWindow().showStatus(selectedGroups+"/"+totalGroups+" groups on, E range total = "+Math.round(selectedRange));
            
    }

    /** Shows a popup menu (context menu) for a button */
    void showMenu(int button, Component parent, MouseEvent e) {
        String[] contextMenuItems = CONTEXT_MENU_ITEMS[Math.min(button, N_BUTTONS)];
        int[] contextMenuFlags = CONTEXT_MENU_FLAGS[Math.min(button, N_BUTTONS)];
        PopupMenu menu = new PopupMenu();
        for (int i=0; i<contextMenuItems.length; i++) {
            String itemStr = contextMenuItems[i];
            NumberedMenuItem item = new NumberedMenuItem(itemStr, button);
            item.addActionListener(this);
            if (contextMenuFlags != null) {
                int flags = contextMenuFlags[i];
                if (LeedUtils.flagSet(flags, ROI_REQUIRED) && plotImp.getRoi() == null)
                    item.setEnabled(false);
                if (LeedUtils.flagSet(flags, EXCLUDED_REQUIRED) && noInputRange[currentColumns[button-N_BUTTONS]] == null)
                    item.setEnabled(false);
                if (LeedUtils.flagSet(flags, COMMENT_REQUIRED) && comment[currentColumns[0]] == null)
                    item.setEnabled(false);
                if (LeedUtils.flagSet(flags, ANY_COMMENT_REQUIRED) && LeedUtils.countNonNull(comment) == 0)
                    item.setEnabled(false);
                if (LeedUtils.flagSet(flags, LAST_COMMENT_REQUIRED) && lastComment == null)
                    item.setEnabled(false);
                if (LeedUtils.flagSet(flags, CURRENT_ON_REQUIRED) && !useOutput[currentColumns[0]])
                    item.setEnabled(false);
            }
            menu.add(item);
        }
        parent.add(menu);
        menu.show(parent, e.getX(), e.getY());
    }

    /** Shows information on a button in the status line */
    void showStatus(int button) {
        String str = null;
        switch(button) {
            case HELP:
                str = "Show Help... Right-click for beam group comments (group annotations)";
                break;
            case SAVE:
                str = "Save I(V) curves... (Right-click to save unsmoothed)";
                break;
            case PREV:
                int g = getPreviousGroup();
                int spot = spotIndices[getColumns(g)[0]];
                str = g==group ? "No previous group" : ("Previous group "+spotPattern.getNameWithGroup(spot, false));
                str += ". Right-click for more navigation.";
                break;
            case NEXT: getNextGroup();
                g = getNextGroup();
                spot = spotIndices[getColumns(g)[0]];
                str = g==group ? "No next group" : ("Next group "+spotPattern.getNameWithGroup(spot, false));
                str += ". Right-click for more navigation.";
                break;
            case ON_OFF:
                str = useOutput[currentColumns[0]] ? "Deselect" : "Select";
                if (currentColumns.length > 1) str += " group";
                str += " " + spotPattern.getNameWithGroup(currentColumns[0], false);
                str += ". Right-click to (de)select all.";
                str += " Shortcut: <space>";
                break;
            case SMOOTH_LESS: case SMOOTH_MORE:
                int sm = smoothR[currentColumns[0]];
                str = "Smooth " + (button==SMOOTH_LESS ? "less" : "more") + " (currently ";
                str += sm <= MIN_SMOOTH ?
                        "off" :
                        (IJ.d2s(LeedSmoother.mToInvNoiseGainSquare(sm)*eStep, E_DECIMALS)+" eV");
                str += "). Right-click for more. Shortcut: <ALT+mouse wheel>";
                break;
            case ROI:
                double[] eMinMax = getRoiEnergyRange(false);
                if (eMinMax==null) {
                    str = "Use rectangle tool to set energy range, or press now for default range.";
                    int[] currentRange = outputRange[currentColumns[0]];
                    if (currentRange != null)
                        str += " Current range: "+IJ.d2s(energies[currentRange[0]], E_DECIMALS)+
                                "-"+energies[Math.max(currentRange[1]-1,0)];
                } else
                    str ="Set Energy Range "+IJ.d2s(Math.max(eMinMax[0],curveStartEv), E_DECIMALS)+
                            "-"+IJ.d2s(Math.min(eMinMax[1], curveEndEv), E_DECIMALS);
                str += ". Right-click for more. Shortcut: # or 0";
                break;
            case OPTIONS:
                str = "Options (V0i, fix x axis, ...) and items related to 'bad beam agreement'";
                break;
            default:
                int curveIndex = button - N_BUTTONS;
                if (curveIndex >= 0) {
                    str = "("+getSpotName(currentColumns[curveIndex])+")";
                    double rFact = getRvsAverage(curveIndex);
                    if (!Double.isNaN(rFact))
                        str += ": R factor vs. average = "+IJ.d2s(rFact, 3);
                    str += " Click to ";
                    str += useInput[currentColumns[curveIndex]] ? "Deselect" : "Select";
                    if (currentColumns.length > 1)
                        str += ", right-click to set exclusion range for this beam";
                }
        }
        if (str == null) return;
        showStatus(str);
    }

    /** Shows a status at the bottom and in the ImageJ Toolbar's status line */
    void showStatus(String str) {
        PlotWindow pw = getPlotWindow();
        pw.showStatus(str);
        IJ.showStatus(str);
    }

    public void showHelp() {
        if (helpDialog != null && helpDialog.isVisible())
            helpDialog.toFront();
        else
            helpDialog = new HTMLDialog(PLUGIN_NAME+" Help", HELP_STRING, false); // non blocking
    }


    /** ImageListener */
    public void imageOpened(ImagePlus imp) {
    }

    public void imageClosed(ImagePlus imp) {
        if (imp == plotImp)
            ImagePlus.removeImageListener(this);
            saveEdit();
    }

    /** Checks for changed size and y scale. Needs ImageJ 1.52u37 or later,
     *  where ImageListener callbacks are executed in the EventQueue */
    public void imageUpdated(ImagePlus imp) {
        if (imp != plotImp) return;
        Plot plot = (Plot)(plotImp.getProperty(Plot.PROPERTY_KEY));
        if (plot == null) return; //should never happen
        if (plotSizeChanged(plot))
            showButtonsAndComment();
        if (plotRangeChanged(plot))
            makeAndShowPlot(true);
    }

    /** MouseListener */
    public void mouseClicked(MouseEvent e) {
        ImageCanvas canvas = getCanvas();
        if (canvas != null) canvas.mouseClicked(e);
    }

    public void mousePressed(MouseEvent e) {
        try {
            ImageCanvas canvas = getCanvas();
            int button = getButton(e.getX(), e.getY());
            if (button >= 0)
                handleButtonPressed(button, e);
            else if (canvas != null)
                canvas.mousePressed(e);
        } catch (Exception ex) {
            IJ.handleException(ex);
        }
    }

    public void mouseReleased(MouseEvent e) {
        ImageCanvas canvas = getCanvas();
        if (canvas != null) canvas.mouseReleased(e);
    }
    public void mouseEntered(MouseEvent e) {
        ImageCanvas canvas = getCanvas();
        if (canvas != null) canvas.mouseEntered(e);
    }

    public void mouseExited(MouseEvent e) {
        ImageCanvas canvas = getCanvas();
        if (canvas != null) canvas.mouseExited(e);
    }

    /** KeyListener for shortcuts */
    public void keyPressed(KeyEvent e) {
        if (guiBlocked()) return;
        try {
            int keyCode = e.getKeyCode();
            char keyChar = e.getKeyChar();
            int modifiers = e.getModifiersEx();
            boolean ctrlOrMeta = (modifiers & (InputEvent.CTRL_DOWN_MASK | InputEvent.META_DOWN_MASK)) != 0;
            if (keyCode == KeyEvent.VK_SPACE) { //for deselected curve, use roi if not previously set otherwise
                if (!useOutput[currentColumns[0]] && plotImp.getRoi() != null && Arrays.equals(outputRange[currentColumns[0]], defaultOutputRange))
                    handleButtonPressed(ROI, null);
                else
                    handleButtonPressed(ON_OFF, null);
            } else if (keyCode == KeyEvent.VK_COMMA  || keyCode == KeyEvent.VK_PAGE_UP   || keyChar == '<') {
                if (ctrlOrMeta)
                    goBadGroup(-1);             //prev/next shortcuts with ctrl or cmd down go to group with bad agreement
                else
                    handleButtonPressed(PREV, null);
            } else if (keyCode == KeyEvent.VK_PERIOD || keyCode == KeyEvent.VK_PAGE_DOWN || keyChar == '>') {
                if (ctrlOrMeta)
                    goBadGroup( 1);
                else
                    handleButtonPressed(NEXT, null);
            } else if (keyCode == KeyEvent.VK_HOME)
                setGroup(lowestGroup, true);
            else if (keyCode == KeyEvent.VK_END)
                setGroup(highestGroup, true);
            else if (keyCode == KeyEvent.VK_NUMBER_SIGN || keyChar == '#' || keyCode == '0')
                handleButtonPressed(ROI, null);
        } catch (Exception ex) {
            IJ.handleException(ex);
        }
    }

    public void keyReleased(KeyEvent e) {}
    public void keyTyped(KeyEvent e) {}

    /** MouseMotionListener: highlighting curves, showing hints in status line
     *  May be also called with e=null; then only the number of groups, energy range
     *  (and badness after 'go to prev/next bad') is shown. */
    public void mouseMoved(MouseEvent e) {
        if (guiBlocked()) return;
        try {
            ImageCanvas canvas = getCanvas();
            int button = e==null ? -1 : getButton(e.getX(), e.getY());
            boolean updatePlot = mouseOverButton != button && (button >= N_BUTTONS || mouseOverButton >= N_BUTTONS);
            mouseOverButton = button;
            if (mouseOverButton >= 0) {
                canvas.setShowCursorStatus(false);
                showStatus(mouseOverButton);
            } else {
                String status = "";
                if (e != null) {
                    Calibration cal = plotImp.getCalibration();
                    Plot plot = (Plot)(plotImp.getProperty(Plot.PROPERTY_KEY));
                    int ox = canvas.offScreenX(e.getX());
                    double eV = plot.descaleX(ox);
                    status = "E = "+IJ.d2s(eV, E_DECIMALS)+".  ";
                }
                status += selectedGroups+"/"+totalGroups+" groups, E range total = "+Math.round(selectedRange);
                String badnessStr = getBadnessString();
                if (badnessStr != null)
                    status +=". "+badnessStr;
                getPlotWindow().showStatus(status);
            }
            if (updatePlot)
                makeAndShowPlot(false);    //highlight 'mouseOver' line
        } catch (Exception ex) {
            IJ.handleException(ex);
        }
    }

    public void mouseDragged(MouseEvent e) {};

    /** ActionListener: For context menu */
    public void actionPerformed(ActionEvent e) {
        try {
            Object src = e.getSource();
            int buttonNumber = src instanceof NumberedMenuItem ?
                    ((NumberedMenuItem)src).getNumber() : -1;
            int buttonColumn = buttonNumber >= N_BUTTONS ?
                    currentColumns[buttonNumber-N_BUTTONS] : -1;

            String cmd = e.getActionCommand();
            if (cmd.equals(SET_COMMENT))          showCommentDialog();
            else if (cmd.equals(LIST_COMMENTS))   listComments();
            else if (cmd.equals(GO_FIRST))        setGroup(lowestGroup, true);
            else if (cmd.equals(GO_LAST))         setGroup(highestGroup, true);
            else if (cmd.equals(GO_HIGHEST))      setGroup(getHighestEnabledGroup(), true);
            else if (cmd.equals(GO_PREV_BAD))     goBadGroup(-1);
            else if (cmd.equals(GO_NEXT_BAD))     goBadGroup( 1);
            else if (cmd.equals(GO_PREV_COMMENT)) goCommentedGroup(-1);
            else if (cmd.equals(GO_NEXT_COMMENT)) goCommentedGroup( 1);
            else if (cmd.equals(SAVE_NO_SMOOTH))  saveCurves(false);
            else if (cmd.equals(SET_SMOOTHING))   showDialog('S');
            else if (cmd.equals(BEST_SMOOTHING))  showDialog('B');
            else if (cmd.equals(SELECT_ALL))      selectAllGroups(true);
            else if (cmd.equals(DESELECT_ALL))    selectAllGroups(false);
            else if (cmd.equals(DEFAULT_E_RANGE)) selectEnergyRange(currentColumns[0], null);
            else if (cmd.equals(OPTION_DIALOG))   showOptionsDialog();
            else if (cmd.equals(LIST_BAD))        listBadBeams();
            else if (cmd.equals(SET_BAD_LIMIT))   showBadnessDialog();
            else if (cmd.equals(FULL_E_RANGE))
                selectEnergyRange(currentColumns[0], new double[]{energies[0], energies[energies.length-1]});
            else if (cmd.equals(SET_E_RANGE_NO_SNAP))
                selectRangeFromRoi(currentColumns[0], false);
            else if (cmd.equals(SET_E_LIMITS))    showDialog('E');
            else if (cmd.equals(SET_MIN_E_SPAN))  showDialog('M');
            else if (cmd.equals(AUTO_SELECT))     showDialog('A');
            else if (cmd.equals(EXCLUDE_RANGE))   excludeEnergyRangeFromRoi(buttonColumn);
            else if (cmd.equals(LAST_COMMENT)) {
                comment[currentColumns[0]] = lastComment;
                showButtonsAndComment();
                changesDone = true;
            } else if (cmd.equals(DELETE_COMMENT)) {
                lastComment = comment[currentColumns[0]];
                comment[currentColumns[0]] = null;
                showButtonsAndComment();
                changesDone = true;
            } else if (cmd.equals(NO_EXCLUDE_BEAM)) {
                noInputRange[buttonColumn] = null;
                setUse(useInput, buttonColumn, true);   //also updates everything
            } else if (cmd.equals(GO_INDEX)) {
                String str = IJ.getString("Beam or beam group number to go to", "");
                if (str==null || str.length()==0) return;
                boolean success = gotoBeam(str, true);
                if (!success) {
                    showStatus("No beam or group '"+str+"'");
                    IJ.beep();
                }
            }

        } catch (Exception ex) {
            IJ.handleException(ex);
        }
    }

    /** MouseWheelListener: for smoothing more/less (with alt key) or navigating in the stack */
    public void mouseWheelMoved(MouseWheelEvent e) {
        if (guiBlocked()) return;
        try {
            int scrollAmount = 0;
            if (e.getScrollType() == MouseWheelEvent.WHEEL_UNIT_SCROLL)
                scrollAmount = e.getUnitsToScroll();
            else         //scroll more (by pages), but we still do the same
                scrollAmount = e.getWheelRotation();
            if (scrollAmount == 0) return;
            if (smallestMouseWheelStep > Math.abs(scrollAmount))
                smallestMouseWheelStep = Math.abs(scrollAmount);
            int delta = (int)Math.round(scrollAmount/(double)smallestMouseWheelStep);
            if ((e.getModifiers() & (InputEvent.ALT_DOWN_MASK|InputEvent.ALT_MASK)) != 0)
                changeSmooth(delta, e);
            else
                handleButtonPressed(scrollAmount < 0 ? PREV : NEXT, null);
        } catch (Exception ex) {
            IJ.handleException(ex);
        }

    }

    static double sqr(double x) { return x*x; }

    /** A MenuItem that also remembers a number (used for the number of the button triggering the PopupMenu) */
    class NumberedMenuItem extends MenuItem {
        int number;
        public NumberedMenuItem(String label, int number) {
            super(label);
            this.number = number;
        }
        public int getNumber() {
            return number;
        }
    }
}
