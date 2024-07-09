import ij.*;
import ij.IJ;
import ij.util.Tools;
import ij.gui.GenericDialog;
import ij.plugin.frame.Recorder;
import ij.measure.ResultsTable;
import ij.io.OpenDialog;
import java.io.File;


/**
 * The numeric and String parameters are kept as static variables in this class
 * (which is possible because the Spot Tracker has a single-instance listener)
 * The class also contains static utility methods relate to this class.
 *
 * Parameters are kept in the ImageJ IJ_Prefs.txt file and saved in the log file
 * created together with the other Spot tacker output files.
 * One Prefs key for all numbers, one Prefs each key for the Strings.
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


public class LeedParams {
    static final int CURRENTVERSION = 6;            // in case of renumbering existing items, increase the version number!

    /** Number of numeric parameters */
    public static final int N_PARAM = 39;           //must be increased when adding numeric parameters
    // KEYS of the individual numeric parameters
    public static final int PARAMVERSION=0,         //version of LeedParams; earlier versions need translation table
            BACKGROUNDTYPE=1,                       //circle (0), oval(1), or azimuthalblur(2) see LeedSpotAnalyzer
            RADIUSINFTYSQR=2,                       //spot integration radius squared at very high energies
            RADIUS1EVSQR=3,                         //increase of integration radius^2 at 1 eV for integer spots
            RADIUS1EVSQRS=4,                        //increase of integration radius^2 at 1 eV for superstructure spots
            AZIMUTHBLURANGLE=5,                     //half-angle of azimuthal spot blurring (radians)
            POSITIONAVERAGINGEV=6,                  //energy range for smoothing positions when tracking spots
            MINSIGNIFICANCETRACK=7,                 //minimum significance to accept spot position during tracking
            SEARCHAGAINEV=8,                        //search for spot again after unseen for ... eV
            NEIGHBORBACKGROUND=9,                   //whether to subtract 1/r^2 background from bright neighboring spots
            I0FROMBACKGR=10,                        //true for I0 high-frequency components from background intensity
            MINSIGNIFICANCEINDEX=11,                //minimum significance to accept spot position during index selection;
            MINRANGEEV=12,                          //min energy range for beam; 0 to try measuring also spots below significance threshold
            RFACTORTYPE=13,                         //R factor type, 0 = RPendry, more will be defined in the future
            V0I=14,                                 //V0i, imaginary part of inner potential for R factor, for R factor
            SMOOTHEV=15,                            //default smoothing energy range for R factors
            SPLITEV=16,                             //split energy ranges for R-factor analysis to roughly this length (eV)
            DARKAVERAGING=17,                       //averaging of dark-frame slices (if a stack)
            FLATAVERAGING=18,                       //averaging of flat-field slices (if a stack)
            DARK2AVERAGING=19,                      //averaging of dark-frame slices (if a stack)
            FLATFITORDER=20,                        //polynomial order for log(flat) fitting, 0 for normalization only, or -1 for none
            METADATASOURCES=21,                     //metadata sources, least-significant digit for energy, then t, I0, I00; 0-4 for none/stack/table/man/prev
            // The following 8 entries must be two each for E, t, I0, I00; this the same sequence as in the energiesEtc array of the spotTracker
            EMANUALFIRST=22,                        //if energy is set manually, first value.
            EMANUALLAST=23,                         //if energy is set manually, last value
            TMANUALFIRST=24,                        //if time is set manually, first value.
            TMANUALLAST=25,                         //if time is set manually, last value
            I0OFFSET=26,                            //Beam current at energy=0 (if linear)
            I0SLOPE=27,                             //Slope of beam current (if linear) w.r.t energy
            I00OFFSET=28,                           //Beam-off current at energy=0
            I00SLOPE=29,                            //Slope of beam-off-current w.r.t energy
            SMOOTHI0POINTS=30,                      //points for smoothing I0
            XAXISVARIABLE=31,                       //x axis variable if not energy
            LEEMMODE=32,                            //In LEEM mode, no 1/sqrt(E) spot movement, uses E_LEEM instead of energy
            SAVEPLOTS=33,                           //whether to save spottracker plots
            SAVESTACK=34,                           //whether to save the spottracker image stack
            OPENSTACKSASVIRTUAL=35,                 //open stacks as vitual stacks (read on the fly)
            UNDISTORTSIZE=36,                       //ouput size of 'undistort' (pixels)
            UNDISTORTTASK=37,                       //ouput type of 'undistort'. 0=current image, 1=stack, 2=stack fixed k-scale
            UNDISTORTNORMALIZE=38;                  //normalize undistorted stack by I0

    /** Number of String parameters */
    public static final int N_STR_PARAM = 8; //max 10 (single-digit, #0-9), otherwise readFromFile must be altered!!
    // KEYS of the individual String parameters in the slice labels of the stack
    public static final int PATTERNFILE=0,SAVEDIRECTORY=1, ENERGYKEYS=2,TIMEKEYS=3,CURRENTI0KEYS=4,CURRENTI00KEYS=5,
            TEMPERATUREKEYS=6, AUXKEYS=7, ELEEMKEYS=8; //Warning: not same as LEED_Spot_Tracker.ENERGY ... AUX: PROCESSEDI0 is missing!

    //names of parameters for macro LEEDsetValue
    static final String[] NUMERIC_PARAM_NAMES = new String[] {
            "PARAM_VERSION",                        // PARAMVERSION
            "backgroundType",                       // BACKGROUNDTYPE
            "rInftySqr",                            // RADIUSINFTYSQR
            "r1eVintSqr",                           // RADIUS1EVSQR
            "r1eVsupSqr",                           // RADIUS1EVSQRS
            "AzimuthBlurAngle",                     // AZIMUTHBLURANGLE
            "positionAveraging_eV",                 // POSITIONAVERAGINGEV
            "minSignificanceTrack",                 // MINSIGNIFICANCETRACK
            "searchAgain_eV",                       // SEARCHAGAINEV
            "neighborBackground",                   // NEIGHBORBACKGROUND
            "I0fromBackgr",                         // I0FROMBACKGR
            "minSignificanceSelecting",             // MINSIGNIFICANCEINDEX
            "minRange_eV",                          // MINRANGEEV
            "rFactorType",                          // RFACTORTYPE
            "V0i",                                  // V0I
            "smoothing_eV",                         // SMOOTHEV
            "splitCurves_eV",                       // SPLITEV
            "darkAveraging",                        // DARKAVERAGING (0=off, 1=average all, 2 = linear, >2 = smooth)
            "flatAveraging",                        // FLATAVERAGING (0=off, 1=average all, >2 = smooth)
            "dark2Averaging",                       // DARK2AVERAGING (0=off, 1=average all, 2 = linear, >2 = smooth)
            "flatFitOrder",                         // FLATFITORDER
            "metaDataSources",                      // METADATASOURCES (0-4 for none/stack/table/man/prev)
            "energyManualFirst",                    // EMANUALFIRST
            "energyManualLast",                     // EMANUALLAST
            "timeManualFirst",                      // TMANUALFIRST
            "timeManualLast",                       // TMANUALLAST
            "I0offset",                             // I0OFFSET
            "I0slope",                              // I0SLOPE
            "I00offset",                            // I00OFFSET
            "I00slope",                             // I00SLOPE
            "smoothI0points",                       // SMOOTHI0POINTS
            "xAxisVariable",                        // XAXISVARIABLE
            "LEEMmode",                             // LEEMMODE
            "savePlots",                            // SAVEPLOTS
            "saveStack",                            // SAVESTACK
            "OpenStacksAsVirtual",                  // OPENSTACKSASVIRTUAL
            "UndistortSize",                        // UNDISTORTSIZE
            "UndistortType",                        // UNDISTORTTASK
            "UndistortNormalize"                    // UNDISTORTNORMALIZE
            };

    //short help for numeric parameters
    static final String[] NUMERIC_PARAM_HELP = new String[] {
            "Version of parameter list (for compatibility check)",          // PARAMVERSION
            "Integration/background type, 0=concentric circles, 1=oval bg, 2=azimuth blur", // BACKGROUNDTYPE
            "Integration radius squared, for energy -> infinity",           // RADIUSINFTYSQR
            "Integration radius squared, extrapolated to energy 1 eV",      // RADIUS1EVSQR
            "Superstructure integration radius squared, extrapolated to energy 1 eV",       // RADIUS1EVSQRS
            "Azimuthally blurred spots: half-angle (\u00B0)",               // AZIMUTHBLURANGLE
            "Position averaging range during spot tracking (in eV if x-axis is energy)",    // POSITIONAVERAGINGEV
            "Minimum spot significance for position tracking",              // MINSIGNIFICANCETRACK
            "Search again if spot is unseen that long during tracking (eV if x-axis is energy)", // SEARCHAGAINEV
            "Whether to subtract 1/r² background of bright neighbor spots", // NEIGHBORBACKGROUND
            "Whether to take I00 (beam-off I0) from the dark frame (1=true)",               // I0FROMBACKGR
            "Minimum spot significance when selecting indices",             // MINSIGNIFICANCEINDEX
            "I(V) curves shorter than this are discarded (in eV if x-axis is energy)",       // MINRANGEEV
            "R factor type, currently only 0=R_Pendry",                     // RFACTORTYPE
            "Imaginary part of the inner potential (absolute value), for R factor calculation",     // V0I
            "Smoothing (in eV), for R factor calculation",                  // SMOOTHEV
            "For R factor statistics of equivalent beams, approx. length of sections (eV)", // SPLITEV
            "Averaging of dark frames: 0=off, 1=average all, 2=linear, >2 for smoothing",   // DARKAVERAGING (0=off, 1=average all, 2 = linear, >2 = smooth)
            "Averaging of flat field: 0=off, 1=average all, 2=linear, >2 for smoothing",    // FLATAVERAGING (0=off, 1=average all, >2 = smooth)
            "Averaging of dark fr. for flat field: 0=off, 1=average all, 2 = linear, >2 = smooth",  // DARK2AVERAGING (0=off, 1=average all, 2 = linear, >2 = smooth)
            "2D polynomial fit order for flat field (-1 = no fit, no normalization)",       // FLATFITORDER
            "Sources for metadata, digits for I00, I0, t, E; 0-4 for none/stack/table/manual/previous", // METADATASOURCES
            "Energy of first stack slice when set manually",                // EMANUALFIRST
            "Energy of last stack slice when set manually",                 // EMANUALLAST
            "Time (s) of first stack slice when set manually",              // TMANUALFIRST
            "Time (s) of last stack slice when set manually",               // TMANUALLAST
            "Offset of I0 when set manually (i.e. value at E = 0)",         // I0OFFSET
            "Slope (per eV) of I0 when set manually",                       // I0SLOPE
            "Offset of I00 (beam-off I0) when set manually (i.e. value at E = 0)",          // I00OFFSET
            "Slope (per eV) of I00 (beam-off I0) when set manually",        // I00SLOPE
            "Smoothing of I0, in points, 0=off",                            // SMOOTHI0POINTS
            "The x axis variable, 0=energy, 1=time, 4=temperature, 5=aux, 6=LEEM mode",     // XAXISVARIABLE
            "Whether LEEM mode is on [then spots don't move radially with 1/sqrt(energy)]", // LEEMMODE
            "Whether to save the plots when saving the result of spot tracking",            // SAVEPLOTS
            "Whether to save the spot tracking image stack when saving the result of spot tracking", // SAVESTACK
            "Open LEED movies as virtual stacks (read on demand)",          // OPENSTACKSASVIRTUAL
            "Size of undistorted output (pixels)",                          // UNDISTORTSIZE
            "Undistort output type: 0=image, 1=stack, 2=stack fixed k-scale",               // UNDISTORTTASK
            "Normalize undistorted output by I0 (0=off, 1=on)"              // UNDISTORTNORMALIZE
            };

    // DEFAULT VALUES
    static final double[] DEFAULT_NUM_VALUES = new double[] {
            CURRENTVERSION,                         // PARAMVERSION
            LeedSpotAnalyzer.CIRCLE,                // BACKGROUNDTYPE
            64,                                     // RADIUSINFTYSQR
            1000,                                   // RADIUS1EVSQR
            1000,                                   // RADIUS1EVSQRS
            3.0,                                    // AZIMUTHBLURANGLE (half-angle in degrees)
            30,                                     // POSITIONAVERAGINGEV
            2.5,                                    // MINSIGNIFICANCETRACK
            30.,                                    // SEARCHAGAINEV
            1,                                      // NEIGHBORBACKGROUND
            0,                                      // I0FROMBACKGR
            2.0,                                    // MINSIGNIFICANCEINDEX
            30,                                     // MINRANGEEV
            LeedRFactor.R_PENDRY,                   // RFACTORTYPE
            5.0,                                    // V0I
            5.,                                     // SMOOTHEV
            100,                                    // SPLITEV
            2,                                      // DARKAVERAGING (0=off, 1=average all, 2 = linear, >2 = smooth)
            50,                                     // FLATAVERAGING (0=off, 1=average all, >2 = smooth)
            2,                                      // DARK2AVERAGING (0=off, 1=average all, 2 = linear, >2 = smooth)
            4,                                      // FLATFITORDER
            1111,                                   // METADATASOURCES ('1' digits for reading from stack: 0-4 for none/stack/table/man/prev)
            Double.NaN,                             // EMANUALFIRST
            Double.NaN,                             // EMANUALLAST
            0,                                      // TMANUALFIRST
            Double.NaN,                             // TMANUALLAST
            1,                                      // I0OFFSET
            0,                                      // I0SLOPE
            0,                                      // I00OFFSET
            0,                                      // I00SLOPE
            20,                                     // SMOOTHI0POINTS
            0,                                      // XAXISVARIABLE
            0,                                      // LEEMMODE
            1,                                      // SAVEPLOTS
            0,                                      // SAVESTACK
            1,                                      // OPENSTACKSASVIRTUAL
            512,                                    // UNDISTORTSIZE
            2,                                      // UNDISTORTTASK
            1                                       // UNDISTORTNORMALIZE
            };

    static final String[] defaultStrs = new String[] {
            "",                                     //PATTERFILE
            "",                                     //SAVEDIRECTORY
            "energy=|E=|V=",                        //ENERGYKEYS  Keys without '=', ':' (if any) must be last
            "time=|Clock=|t=",                      //TIMEKEYS
            "current=|i=|i0=|I0=|AD_0=",            //CURRENTKEYS
            "I00=",                                 //CURRENTI00KEYS
            "T=|Temp=|temp=|temperature=",          //TEMPERATUREKEYS
            "aux=|Aux=|AUX=|AD_2="                  //AUXKEYS
            };

    // keys for ImageJ Prefs (file IJ_Prefs.txt)
    public static final String PREFS_KEY = "leed.tracker";
    static final String PREFS_KEY_NUM = PREFS_KEY+"_n";   //key in Prefs file for the String with all numbers
    static final String PREFS_KEY_STR = PREFS_KEY+"_s";   //key prefix in Prefs file for the Strings, gets "leed.iv.s0", "leed.iv.s1" etc.

    // for keywords dialog; null means don't ask
    static final String[] OPTIONS_STRINGLABELS = new String[] {
            null, null,                             //PATTERFILE, SAVEDIRECTORY
            "Keywords for energy in files",         //ENERGYKEYS
            "Keywords for time in files",           //TIMEKEYS
            "Keywords for beam current I0 in files",//CURRENTKEYS
            "Keywords for no-beam current I00 in files", //CURRENTI00KEYS
            "Keywords for temperature in files",    //TEMPERATUREKEYS
            "Keywords for AUX channel in files"     //AUXKEYS
    };

    // Translation table mapping parameters for previous versions to keep compatibility.
    // Values -1 for parameters with no equivalence in current version
    static final int[][] PREVIOUS_VERSION_PARAMS =
            new int[][] {null, null, null, null,    //versions 0-3 not implemented
                new int[] {                     //version 4 till 20220110
                    PARAMVERSION, BACKGROUNDTYPE, RADIUSINFTYSQR, RADIUS1EVSQR, RADIUS1EVSQRS,      //0-4
                    POSITIONAVERAGINGEV, MINSIGNIFICANCETRACK, SEARCHAGAINEV, MINSIGNIFICANCEINDEX, V0I, //5-9
                    SMOOTHEV, SPLITEV, DARKAVERAGING, FLATAVERAGING, DARK2AVERAGING,                //10-14
                    FLATFITORDER, SMOOTHI0POINTS, I0OFFSET, I0SLOPE, I00OFFSET,                     //15-19
                    I00SLOPE, -1, XAXISVARIABLE, -1, -1, -1, -1,                                    //20-24
                    LEEMMODE, NEIGHBORBACKGROUND, SAVEPLOTS, SAVESTACK},                            //25-28
                new int[] {                     //version 5 till 20240610 (still including IV Editor)
                    PARAMVERSION, BACKGROUNDTYPE, RADIUSINFTYSQR, RADIUS1EVSQR, RADIUS1EVSQRS,      //0-4
                    POSITIONAVERAGINGEV, MINSIGNIFICANCETRACK, SEARCHAGAINEV, NEIGHBORBACKGROUND, -1,   //5-9
                    MINSIGNIFICANCEINDEX, MINRANGEEV, RFACTORTYPE, V0I, SMOOTHEV,                    //10-14
                    SPLITEV, DARKAVERAGING, FLATAVERAGING, DARK2AVERAGING, FLATFITORDER,            //15-19
                    METADATASOURCES, EMANUALFIRST, EMANUALLAST, TMANUALFIRST, TMANUALLAST,          //20-24
                    I0OFFSET, I0SLOPE, I00OFFSET, I00SLOPE, SMOOTHI0POINTS,                         //25-29
                    I0FROMBACKGR, XAXISVARIABLE, LEEMMODE, SAVEPLOTS, SAVESTACK,                    //30-34
                    -1, -1, -1, -1, AZIMUTHBLURANGLE,                                               //35-39
                    OPENSTACKSASVIRTUAL, UNDISTORTSIZE, UNDISTORTTASK, UNDISTORTNORMALIZE}          //40-43
            };
            
    // parameters are stored here
    static double[] params = new double[N_PARAM];
    static String[] strParams = new String[N_STR_PARAM];

    static boolean initialized;

    /** When called the first time, reads the parameters from the ImageJ Prefs; uses the defaults if not in the Prefs.
     *  If not done till then, will be called with the first get or set operation. */
    public static void initialize() {
        if (initialized) return;
        boolean ok = false;
        System.arraycopy(DEFAULT_NUM_VALUES, 0, params, 0, DEFAULT_NUM_VALUES.length);
        String numStr = Prefs.get(PREFS_KEY_NUM, (String)null);
        if (numStr != null)
            ok = readNumParamFromString(numStr, params);
        if (!ok) { //we have no prefs yet or read a version that does not fit; it can't be used
            reset();
            initialized=true;
            return;
        }
        // read the strings
        for (int i=0; i<N_STR_PARAM; i++)
            strParams[i] = Prefs.get(PREFS_KEY_STR+i, i<defaultStrs.length ? defaultStrs[i] : (String)null);
        initialized = true;
    }

    /** Saves the parameters in the Prefs */
    public static void saveToPrefs() {
        if (!initialized) initialize();
        // write all numbers in one line
        Prefs.set(PREFS_KEY_NUM, getNumbersLine());
        // write the strings
        for (int i=0; i<N_STR_PARAM; i++)
            if (strParams[i] != null)
                Prefs.set(PREFS_KEY_STR+i, strParams[i]);
        initialized = false;
    }

    /** Returns a String with all numeric parameters in one line */
    static String getNumbersLine() {
        return LeedUtils.toString(params);
    }

    /** Resets all parameters to the defaults (if they have defaults) */
    public static void reset() {
        System.arraycopy(DEFAULT_NUM_VALUES, 0, params, 0, DEFAULT_NUM_VALUES.length);
        System.arraycopy(defaultStrs, 0, strParams, 0, defaultStrs.length);
        initialized = true;
        saveToPrefs();
    }

    /** Returns the numeric parameter for the given key */
    public static double get(int key) {
        if (!initialized) initialize();
        return params[key];
    }

    /** Returns the parameter for the given key as boolean. Values != 0 are considered true. */
    public static boolean getBoolean(int key) {
        return get(key) != 0;
    }

    /** Returns an array of n numeric parameters starting with the given key */
    public static double[] getArray(int key, int n) {
        if (!initialized) initialize();
        double[] array = new double[n];
        System.arraycopy(params, key, array, 0, n);
        return array;
    }

    /** Returns the default numeric parameter for the given key */
    public static double getDefaultValue(int key) {
        return DEFAULT_NUM_VALUES[key];
    }


    /** Returns the String parameter for the given key */
    public static String getString(int key) {
        if (!initialized) initialize();
        return strParams[key];
    }

    /** Returns the String parameter array for the given key.
     *  The String array is stored as one String with '|' characters as separators */
    public static String[] getStringArray(int key) {
        if (!initialized) initialize();
        return Tools.split(strParams[key], "|");
    }

    /** Sets the numeric parameter for the given key */
    public static void set(int key, double value) {
        if (!initialized) initialize();
        if (key == PARAMVERSION) return;
        params[key] = value;
        Prefs.set(PREFS_KEY_NUM, getNumbersLine());
        recordParameter(key, value);
    }

    /** Sets the numeric parameter for the given key
     *  Boolean true and false are written as 1 and 0, respectively. */
    public static void set(int key, boolean b) {
        set(key, b ? 1 : 0);
    }

    /** Sets the numeric parameter for the given macro parameter name. Returns null if ok, error message if failed */
    public static String setValue(String name, double value) {
        if (!initialized) initialize();
        int key = LeedUtils.arrayIndexOf(NUMERIC_PARAM_NAMES, name);
        if (name == null || key <= PARAMVERSION)
            return "ERROR: no parameter '"+name+"'";
        params[key] = value;
        return null;
    }

    /** Sets the numeric parameters starting with the given key with the values from the array.
     *  All array elements are copied. */
    public static void setFromArray(int key, double[] array) {
        if (!initialized) initialize();
        System.arraycopy(array, 0, params, key, array.length);
        Prefs.set(PREFS_KEY_NUM, getNumbersLine());
    }

    /** Sets the String parameter for the given key */
    public static void setString(int key, String str) {
        if (!initialized) initialize();
        strParams[key] = str;
        Prefs.set(PREFS_KEY_STR+key, strParams[key]);
    }

    /** Returns all parameters as String array of lines as they are in the ImageJ Prefs file */
    public static String[] getParameterLines() {
        String[] out = new String[N_STR_PARAM+1];
        out[0] = PREFS_KEY_NUM+"="+getNumbersLine();
        for (int i=0; i<N_STR_PARAM; i++)
            out[i+1] = (PREFS_KEY_STR+i)+"="+ (strParams[i] == null ? "" : strParams[i]);
        return out;
    }

    /** Shows a table with the numeric parameters (mainly for debugging).
     *  When macro recording is on, also macro-records all numeric parameters */
    public static void makeTable(String title) {
        ResultsTable rt = new ResultsTable();
        for (int i=0; i<N_PARAM; i++) {
            rt.incrementCounter();
            rt.addLabel(NUMERIC_PARAM_NAMES[i]);
            rt.addValue("Value", params[i]);
            rt.addValue("Default", DEFAULT_NUM_VALUES[i]);
            rt.addValue("Description", NUMERIC_PARAM_HELP[i]);
            recordParameter(i, params[i]);
        }
        rt.setPrecision(1); //One digit behind decimal point
        rt.show(title);
    }

    /** Macro-records a numeric parameter if recording is on */
    public static void recordParameter(int key, double value) {
        if (Recorder.record && key != PARAMVERSION)
            LEED_Spot_Tracker.recordMacro(LEED_Spot_Tracker.MACRO_COMMAND_SETVALUE, '"'+NUMERIC_PARAM_NAMES[key]+"\", "+LeedUtils.d2s(value));
    }

    /** Dialog for parameters ('Rare Options') that are not accessible through other dialogs */
    public static void showParamsDialog() {
        if (!initialized) initialize();
        boolean reset = false;
        java.awt.Point loc = null;
        do {
            GenericDialog gd = new GenericDialog(LEED_Spot_Tracker.PLUGIN_NAME+" Metadata Keywords");
            for (int s=0; s<N_STR_PARAM; s++) {
                String label = OPTIONS_STRINGLABELS[s];
                if (label == null) continue;
                gd.addStringField(label, strParams[s], 32);
            }
            gd.addMessage("Keywords: separate multiple ones by '|'; include '=' or ':' if present\n"+
                    "Uppercase/lowercase matters\n"+
                    "Deselect or close+open files to apply");
            gd.enableYesNoCancel("OK", "Reset to defaults");
            if (loc != null) gd.setLocation(loc.x, loc.y);
            gd.showDialog();
            if (gd.wasCanceled()) return;
            reset = !gd.wasOKed();
            loc = gd.getLocation();
            for (int s=0; s<N_STR_PARAM; s++) {
                if (OPTIONS_STRINGLABELS[s] == null) continue;
                if (reset)
                    strParams[s] = defaultStrs[s];
                else
                    strParams[s] = gd.getNextString();
            }
        } while (reset); //redisplay dialog after reset
    }

    /** Reads the parameters from a (_log.txt) file.
     *  Returns an array of ScreenFitter parameters or null on success, otherwise a String starting with "ERROR:"
     *  When 'compare is true, creates a list with a comparison of the numeric parameters. */
    public static Object readFromFile(String path, boolean compare) {
        String str = IJ.openAsString(path);
        if (str == null) return null;
        if (str.startsWith("Error:"))
            return "ERROR: Can't read file,"+str.substring("Error:".length());
        String filename = (path==null || path.length()==0) ? OpenDialog.getLastName() : (new File(path)).getName();
        if (filename.endsWith(".txt")) filename = LeedUtils.removeExtension(filename);
        boolean anythingRead = false;
        double[] screenFitterArray = null;
        String[] lines = Tools.split(str, "\n");
        for (int i=0; i<lines.length; i++) {
            String line = lines[i];
            if (line.startsWith(PREFS_KEY_NUM+"=")) {
                double[] output = compare ? new double[N_PARAM] : params;
                boolean ok = readNumParamFromString(line.substring(PREFS_KEY_NUM.length()+1), output);
                if (!ok) return "ERROR: File is from an incompatible old version.";
                anythingRead = true;
                if (compare)
                    tabulateDifferences(filename, output, params);
            } else if (line.startsWith(LEED_Spot_Tracker.SCREENFIT_KEY+"=") && !compare) {
                screenFitterArray = LeedUtils.numbersFromString(line.substring(LEED_Spot_Tracker.SCREENFIT_KEY.length()+1));
                anythingRead = true;
            } else if (line.startsWith(PREFS_KEY_STR) && line.length() > PREFS_KEY_STR.length()+2 && !compare) {
                int nStr = line.charAt(PREFS_KEY_STR.length()) - '0';   //assumes single digit
                if (nStr >= 0 && nStr < N_STR_PARAM) {
                    strParams[nStr] = line.substring(PREFS_KEY_STR.length()+2);  //assumes single digit
                    anythingRead = true;
                }
            }
        }
        if (anythingRead)
            return screenFitterArray;
        else
            return "ERROR: No parameters found in file";
    }

    /** Displays a table of numeric parameter differences */
    private static void tabulateDifferences(String filename, double[] fileParams, double[] params) {
        ResultsTable rt = new ResultsTable();
        for (int i=0; i<N_PARAM; i++) {
            if (fileParams[i] != params[i]) {
                rt.incrementCounter();
                rt.addLabel(NUMERIC_PARAM_NAMES[i]);
                rt.addValue(filename, fileParams[i]);
                rt.addValue("Current", params[i]);
                rt.addValue("Description", NUMERIC_PARAM_HELP[i]);
            }
        }
        if (rt.size() == 0) {
            rt.incrementCounter();
            rt.addLabel("Parameters in "+filename);
            rt.incrementCounter();
            rt.addLabel("are the same as the current parameters");
        }
        rt.setPrecision(2); //Two digits behind decimal point (needed for EditorNoiseLimit)
        rt.show("Numeric Parameter Comparison "+LeedUtils.getDateFormatted(" [HH'h'mm]"));
    }

    /** Reads the numeric parameters from a String of comma-delimited values into the 'output' array
     *  Returns false in case of a version mismatch with no translation table
     *  Leaves the parameter untouched if the value read in is NaN */
    private static boolean readNumParamFromString(String numStr, double[] output) {
        String[] parts = Tools.split(numStr, ",");
        int version = 0;
        int[] translationTable = null;
        for (int i=0; i<Math.min(parts.length, N_PARAM); i++) {
            double v = Tools.parseDouble(parts[i]);
            if (i == PARAMVERSION) {
                version = (int)v;
                if (version > CURRENTVERSION)
                    return false;
                if (version != CURRENTVERSION) {
                    translationTable = PREVIOUS_VERSION_PARAMS[version];
                    if (translationTable == null)
                        return false;
                }
            } else if (!Double.isNaN(v)) {
                if (translationTable == null)
                    output[i] = v;
                else if (i<translationTable.length && translationTable[i] > 0)
                    output[translationTable[i]] = v;
            }
        }
        return true;
    }
}
