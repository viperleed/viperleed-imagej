import ij.*;
import ij.IJ;
import ij.gui.*;
import ij.process.*;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.frame.Recorder;
import ij.util.Tools;
import ij.util.ThreadUtil;
import ij.util.IJMath;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.lang.ref.*;

/** This class implements the functionality of the 'Track Spots' button of the Spot Tracker.
 *  It shows the dialog, does the spot tracking and I(V) measurement on a stack
 *  and then creates plots and tables with the results. */

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


public class LeedIVAnalyzer implements Runnable {
    /** Names of plot files for the plot types */
    public static final String[] PLOT_NAMES = new String[] {"IVcurves", "spotXYdelta", "SpotRadii", "QualityStat", "overallXYdelta"};

    // for shorter notation: Items in LeedSpotAnalyzer output array
    static final int X = LeedSpotAnalyzer.X, Y = LeedSpotAnalyzer.Y,
            INTEGRAL = LeedSpotAnalyzer.INTEGRAL, PSIGNIFICANCE = LeedSpotAnalyzer.PSIGNIFICANCE;
    //  number of extra data columns not created by the LeedSpotAnalyzer: DELTAX, DELTAY, DELATAXS, DELTAYS, INT_I0CORRECTED
    static final int N_EXTRA_COLUMNS = 5;
    //  deviations of initial spot tracking from fit are written into the DELTAX, DELTAY data columns.
    //  DELTAX, DELTAY has three sources:
    //  - detected spots; then data[PSIGNIFICANCE] > minSignificance (and retainedSignificance during findSpots > NO_SIGNIFICANCE)
    //  - retained values at energies close to a detection, then significance during findSpots is in retainedSignificance array
    //  - positions inferred from the neighborsm, then data[PSIGNIFICANCE] = LOW_SIGNIFICANCE
    static final int DELTAX = LeedSpotAnalyzer.OUTPUT_SIZE, DELTAY = DELTAX + 1;
    // data columns for deviations after smoothing
    static final int DELTAXS = LeedSpotAnalyzer.OUTPUT_SIZE+2, DELTAYS = DELTAXS + 1;
    // data columns for I0-corrected intensities
    static final int INT_I0CORRECTED = DELTAYS + 1;
    // names for the extra data columns
    static final String[] DATA_NAMES = new String[] {"dx_raw", "dy_raw", "dx_smooth","dy_smooth", "IVcurves"};
    // Additional columns in regressionResults2D for overall deviations: delta x, delta y fits, scale fit
    static final int REG2D_DELTAX = Leed2DRegression.OUTPUT_SIZE,
            REG2D_DELTAY = Leed2DRegression.OUTPUT_SIZE+1,
            REG2D_SCALE_FIT = Leed2DRegression.OUTPUT_SIZE+2,
            REG2D_ANGLE_FIT = Leed2DRegression.OUTPUT_SIZE+3;
    // the spotAnalyzers NO_SIGNIFICANCE; data[P_SIGNIFICANCE] from the spotAnalyzer is >= this or 0.
    static final double NO_SIGNIFICANCE = LeedSpotAnalyzer.NO_SIGNIFICANCE;
    // if the position can be only inferred from the neighbors, they are set to this significance
    // (must be lower than 0.5* NO_SIGNIFICANCE)
    static final double LOW_SIGNIFICANCE = 0.05;
    // the weight corresponding to LOW_SIGNIFICANCE
    static final double LOW_WEIGHT = 0.005;
    // maximum deviation from calculated position when based on deviation of neighbors, as a multiple of the radius of theintegration disk
    static final double MAX_DEVIATION_RADII = 2;
    // spots with less integral (compared to the highest one) than this don't get 1/r^2 background subtracted
    static final double MIN_REL_INTEGRAL_FOR_BG = 1e-2;

    LEED_Spot_Tracker spotTracker;
    LeedSpotPattern spotPattern;
    ImagePlus stackImp;
    ByteProcessor maskIp;
    Roi maskRoi;
    Rectangle maskRoiRect;
    double[] energies;
    double[] processedI0;
    double[] xAxis;
    String xAxisLabel;
    double[] radiiInt, radiiSup;
    LeedScreenFitter screenFitter;
    int sliceOfIndexInput;                  //slice where ScreenFitter has been applied, minus one
    boolean interactive;                    //whether to show the dialog

    int spotBackgrShape;
    double azBlurRadians;
    double maskCenterX, maskCenterY;        //center of screen for spotBackgrShape=OVAL, or (0,0) spot for AZIMUTH_BLUR
    double minSignificance;
    int minPointsPerBeam;
    boolean debug = false;
    int debugSpot = -1;                     //spot index for debug messages on this spot
    double[][][] data;                      //[type][spotIndex][energyIndex]
    int[] badness;                          // >0 in case of questionable tracking quality
    boolean[] hasSpot;                      //true for spots with valid intensity data
    double[] backgroundIntensities;         //background intensities for I0 correction
    double[][] regressionResults2D;         //x, y offset, scale vs x axis of overall deviation from polynomial model
    double xOffsetInfty = Double.NaN, yOffsetInfty = Double.NaN;  //extrapolated (0,0) position at E -> infinity
    double xOffset100eV = Double.NaN, yOffset100eV = Double.NaN;  //magnetic-field influence on position at 100 eV
    double invScaleInfty = Double.NaN, deltaPhi = Double.NaN;     //1/scale factor at E->infty, estimated work function diff
    double highestIntensity;
    double energyOfHighestI;
    int spotOfHighestI;
    int nGoodSpots;
    int nTooClose;                          //number of spots too close for a circular background
    double progress, progressInc, nextProgressToShow=-1;
    static boolean quickOutput;
    static WeakHashMap<ImagePlus,String> plotImpList = new WeakHashMap<ImagePlus,String>();
    static int currentRunNumber = 1;        //each run gets a higher number, for numbering plots in plotImpList

    /** Creator. If this is not a standard LEED I(V) experiment with spots moving with 1/sqrt(E),
     *  'energies' should be null. */
    public LeedIVAnalyzer(LEED_Spot_Tracker spotTracker, LeedSpotPattern spotPattern,
            ImagePlus stackImp, ImagePlus maskImp, Roi maskRoi,
            double[] energies, double[] processedI0, double[] xAxis, String xAxisLabel,
            LeedScreenFitter screenFitter, int sliceOfIndexInput, boolean interactive) {
        this.spotTracker = spotTracker;
        this.spotPattern = spotPattern;
        this.stackImp = stackImp;
        this.maskIp = LeedUtils.getMaskIp(maskImp);
        this.maskRoi = maskRoi;
        maskRoiRect = maskRoi.getBounds();
        this.energies = energies;
        this.processedI0 = processedI0;
        this.xAxis = xAxis;
        this.xAxisLabel = xAxisLabel;
        this.screenFitter = screenFitter;
        this.sliceOfIndexInput = sliceOfIndexInput;
        this.interactive = interactive;

        spotBackgrShape = (int)LeedParams.get(LeedParams.BACKGROUNDTYPE);
        azBlurRadians = Math.toRadians(LeedParams.get(LeedParams.AZIMUTHBLURANGLE));
        Rectangle maskRect = maskRoi.getBounds();
        if (spotBackgrShape == LeedSpotAnalyzer.AZIMUTH_BLUR) { // azimuthally blurred spots: center for spotAnalyzer is (0,0) spot
            double[] screenXY = screenFitter.screenCoordinates(0.0, 0.0,
                    1./Math.sqrt(LeedRadiusSelector.UNKNOWN_ENERGY), new double[2]);
            maskCenterX = screenXY[0];
            maskCenterY = screenXY[1];
        } else {  //(relevant for spotBackgrShape==OVAL only)
            maskCenterX = maskRect.x + 0.5*maskRect.width;      // center for spotAnalyzer is screen (mask) center
            maskCenterY = maskRect.y + 0.5*maskRect.height;
        }

        boolean hasSuperstructure = spotPattern.isSuperstructure();
        radiiInt = new double[xAxis.length];
        radiiSup = hasSuperstructure ? new double[xAxis.length] : radiiInt;
        for (int i=0; i<xAxis.length; i++) {
            radiiInt[i] = LeedRadiusSelector.radius(getEnergy(i), false);
            if (hasSuperstructure)
                radiiSup[i] = LeedRadiusSelector.radius(getEnergy(i), true);
        }
    }

    /**
     * Displays a dialog, then analyzes the stack and populates the data array
     */
    public void run() {
        try {
            minSignificance = LeedParams.get(LeedParams.MINSIGNIFICANCETRACK);
            double searchAgain = LeedParams.get(LeedParams.SEARCHAGAINEV);
            double positionAveraging = LeedParams.get(LeedParams.POSITIONAVERAGINGEV);
            double minEnergyRange = LeedParams.get(LeedParams.MINRANGEEV);
            boolean subtractNeighborBackground = LeedParams.getBoolean(LeedParams.NEIGHBORBACKGROUND);
            boolean i0FromBackgr = LeedParams.getBoolean(LeedParams.I0FROMBACKGR);

            double v0i = LeedParams.get(LeedParams.V0I);
            double eVsmooth = LeedParams.get(LeedParams.SMOOTHEV);
            double eVseparate = LeedParams.get(LeedParams.SPLITEV);

            boolean useEnergies = xAxis == energies;
            String eVstr = useEnergies ? "eV " : "";
            double step = 1.;
            if (useEnergies) {
                step = LeedUtils.getEnergyStep(energies);
                if (!(step > 0))
                    throw new RuntimeException("Error: Invalid energy step: "+LeedUtils.d2s(step)+" eV");
            }
            if (!(eVseparate < 1e6 && eVseparate > 10*step)) eVseparate = 0;

            boolean canCorrectI0fromBackgr = processedI0 != null &&
                    LeedParams.get(LeedParams.SMOOTHI0POINTS) >= LeedEnergyI0Selector.MIN_SMOOTHING;

            if (hasBadFlat(true)) {
                prepareExit(false);
                return;
            }
            String errStr = processedI0error();
            if (errStr != null) {
                IJ.error(LEED_Spot_Tracker.PLUGIN_NAME, errStr);
                prepareExit(false);
                return;
            }

            if (interactive) {
                GenericDialog gd = new GenericDialog("Spot Tracking Options");
                int v0iIndex = -1;
                gd.addNumericField("Noise rejection", minSignificance, 1, 6, "(1-10; typ. 2.0-3.0)");
                gd.addNumericField("Search beam unseen for", searchAgain, 0, 6, eVstr+"(typ. 30)");
                gd.addNumericField((useEnergies ? "Energy" : "Data")+" range for x, y smoothing", positionAveraging, 0, 6, eVstr+"(typ. 30)");
                gd.addNumericField("Min. "+(useEnergies ? "energy" : "data")+" range per beam", minEnergyRange, 0, 6, eVstr+"(0 to measure sub-threshold)");
                if (spotBackgrShape == LeedSpotAnalyzer.CIRCLE || spotBackgrShape == LeedSpotAnalyzer.OVAL)
                    gd.addCheckbox("Subtract background of bright neighbor spots", subtractNeighborBackground);
                if (canCorrectI0fromBackgr)
                    gd.addCheckbox("Apply fast background changes to I0", i0FromBackgr);
                if (useEnergies && spotPattern.hasEquivalentBeams()) {
                    gd.addMessage("Options for R factor analysis of equivalent beams:");
                    gd.addNumericField("V0i", v0i, 1, 6, "eV");
                    v0iIndex = gd.getNumericFields().size() - 1;
                    gd.addNumericField("Energy range for smoothing", eVsmooth, 1, 6, "eV  *");
                    gd.addNumericField("Split energy range into", eVseparate, 0, 6, "eV sections **");
                    gd.addCheckbox("Quick Plots only ***", quickOutput);
                    gd.addMessage("*    Typically, use V0i. Enter 0 or leave empty for no smoothing.\n"+
                            "**   Enter 0 or leave empty to always analyze the full range\nper pair.\n"+
                            "*** One I/V curve plot of equivalent beams and quality statistics.");
                }
                gd.addCheckbox("Create debug info", false);
                gd.addStringField("Detailed debug info for spot", debugSpot < 0 ? "" : spotPattern.getName(debugSpot));
                gd.addHelp(LeedSpotTrackerHelp.getTrackSpotsHelp());

                final int N_NUMFIELDS = useEnergies && spotPattern.hasEquivalentBeams() ? 7 : 4;
                final double[] minDialogValues = new double[] {NO_SIGNIFICANCE, 0, 2*step, 0, //min values for minSignificance ... minEnergyRange
                        1e-10, 0, 10*step}; //min values for R factor analysis of equivalent beams
                final int v0iIndexF = v0iIndex;
                DialogListener dialogListener = new DialogListener() {  //only checks lower limit of numeric inputs
                    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
                        for (int i=0; i<N_NUMFIELDS; i++) {
                            double x = gd.getNextNumber();
                            if (i == v0iIndexF) x = Math.abs(x); //strictly speaking, V0i is negative. We take the absolute value
                            if (!(x >= minDialogValues[i])) return false;
                        }
                        return true;
                    }
                };
                Button okButton = gd.getButtons()[0];
                okButton.setEnabled(dialogListener.dialogItemChanged(gd, null));
                gd.addDialogListener(dialogListener);
                gd.showDialog();

                if (gd.wasCanceled()) {
                    prepareExit(false);
                    return;
                }

                minSignificance = gd.getNextNumber();       //read dialog fields
                if (!(minSignificance >=1 && minSignificance<=10))
                    minSignificance = LeedParams.get(LeedParams.MINSIGNIFICANCETRACK);
                searchAgain = gd.getNextNumber();
                positionAveraging = gd.getNextNumber();
                minEnergyRange = gd.getNextNumber();
                if (spotBackgrShape == LeedSpotAnalyzer.CIRCLE || spotBackgrShape == LeedSpotAnalyzer.OVAL)
                    subtractNeighborBackground = gd.getNextBoolean();
                if (canCorrectI0fromBackgr)
                    i0FromBackgr = gd.getNextBoolean();
                else
                    i0FromBackgr = false;
                if (useEnergies && spotPattern.hasEquivalentBeams()) {
                    v0i = Math.abs(gd.getNextNumber());
                    eVsmooth = gd.getNextNumber();
                    if (Double.isNaN(eVsmooth)) eVsmooth = 0;
                    eVseparate = gd.getNextNumber();
                    if (!(eVseparate > 10)) eVseparate = 1e7;
                    quickOutput = gd.getNextBoolean();
                }
                debug = gd.getNextBoolean();
                String debugSpotName = gd.getNextString();
                debugSpot = spotPattern.getIndex(debugSpotName);

                LeedParams.set(LeedParams.MINSIGNIFICANCETRACK, minSignificance);
                LeedParams.set(LeedParams.SEARCHAGAINEV, searchAgain);
                LeedParams.set(LeedParams.POSITIONAVERAGINGEV, positionAveraging);
                LeedParams.set(LeedParams.MINRANGEEV, minEnergyRange);
                LeedParams.set(LeedParams.NEIGHBORBACKGROUND, subtractNeighborBackground);
                if (canCorrectI0fromBackgr)
                    LeedParams.set(LeedParams.I0FROMBACKGR, i0FromBackgr);
                if (useEnergies && spotPattern.hasEquivalentBeams()) {
                    LeedParams.set(LeedParams.V0I, v0i);
                    LeedParams.set(LeedParams.SMOOTHEV, eVsmooth);
                    LeedParams.set(LeedParams.SPLITEV, eVseparate);
                }
            }
            if (!spotTracker.waitForStackUpdate()) {
                prepareExit(false);
                return;
            }

            long t0 = System.currentTimeMillis();
            IJ.getInstance().toFront();            //ensure the progress bar can be seen
            int stackSize = stackImp.getNSlices();
            if (stackSize != xAxis.length)
                throw new RuntimeException("Stack size and number of x points (energies) do not match!");

            // the work starts here...
            spotTracker.setButtonLabelText(LEED_Spot_Tracker.TRACK_BUTTON, "tracking...");
            IJ.showStatus("tracking spots...");
            LeedDarkFlatVirtualStack ldfvStack = spotTracker.getStack();
            if (ldfvStack != null) ldfvStack.setShowProgress(false);
            IJ.showProgress(0.0);
            progressInc = 1./(3.5*stackSize);   //two passes to find positions, one to measure, and time for plots

            if (useEnergies) {                  //eV to number of stack slices
                searchAgain /= step;
                positionAveraging /= step;
                minEnergyRange /= step;
            }
            data = new double[LeedSpotAnalyzer.OUTPUT_SIZE+N_EXTRA_COLUMNS][spotPattern.size()][stackSize];
            badness = new int[spotPattern.size()];
            hasSpot = new boolean[spotPattern.size()];
            for (int i=0; i<data.length; i++)
                for (int j=0; j<data[0].length; j++)
                    Arrays.fill(data[i][j], Double.NaN);
            double[] lastDeltaX = new double[spotPattern.size()];        //remembers offsets from calculated position where the
            double[] lastDeltaY = new double[spotPattern.size()];        //   spot has been detected the last time
            int[]   lastSeen = new int[spotPattern.size()];             //remembers stack index when spot has been seen the last time
            Arrays.fill(lastSeen, Integer.MIN_VALUE);
            Arrays.fill(lastDeltaX, Float.NaN);
            backgroundIntensities = new double[xAxis.length];
            double[][] retainedSignificance = new double[spotPattern.size()][xAxis.length]; // of positions not measured, using retained DELTAX,Y
            double[] decayingSignificance = new double[spotPattern.size()];  // significance of last measured position
            double significanceDecay = Math.pow(NO_SIGNIFICANCE/(2*minSignificance),
                    1./searchAgain);                                    // decayingSignificance is multiplied with this each time

            if (debug) {IJ.wait(200); IJ.log("SEARCH: ascend");}
            for (int i=sliceOfIndexInput; i<stackSize; i++)             //from start upwards in energy>; note that i is from 0 to nSlices-1
                findSpots(i, Math.max(i-1, sliceOfIndexInput), lastDeltaX, lastDeltaY, lastSeen, searchAgain,
                        retainedSignificance, decayingSignificance, significanceDecay);
            Arrays.fill(decayingSignificance, 0);
            if (debug) {IJ.wait(200); IJ.log("SEARCH: descend");}
            for (int i=stackSize-2; i>=0; i--)                          //all the way downwards in energy
                findSpots(i, Math.min(i+1, stackSize-1), lastDeltaX, lastDeltaY, lastSeen, searchAgain,
                        retainedSignificance, decayingSignificance, significanceDecay);
            Arrays.fill(decayingSignificance, 0);
            boolean dataAdded = false;
            if (debug) {IJ.wait(200); IJ.log("SEARCH: ascend again");}
            for (int i=0; i<sliceOfIndexInput; i++)                     //from zero upwards in energy
                if (findSpots(i, Math.max(i-1,0), lastDeltaX, lastDeltaY, lastSeen, searchAgain,
                        retainedSignificance, decayingSignificance, significanceDecay));
                    dataAdded = true;
            if (dataAdded || minEnergyRange==0) {                       //only if the last pass has added spots or we want x&y from up&down pass
                if (debug) {IJ.wait(200); IJ.log("SEARCH: ascend beyond start");}
                progressInc *= 0.1;
                for (int i=sliceOfIndexInput; i<stackSize; i++)         //from start up to the end in energy
                    findSpots(i, Math.max(i-1,0), lastDeltaX, lastDeltaY, lastSeen, searchAgain,
                            retainedSignificance, decayingSignificance, significanceDecay);
                progressInc *= 10;
            }
            if (hasBadFlat(true)) {
                prepareExit(false);
                return;
            }
            long t1 = System.currentTimeMillis();
            if (debug) IJ.log("tracking: "+IJ.d2s(0.001*(t1-t0))+" s");

            nGoodSpots = 0;
            for (int spot=0; spot<spotPattern.size(); spot++) {         //which spots have we found?
                int nOk = 0;
                for (int i=0; i<xAxis.length; i++)
                    if (!Double.isNaN(data[X][spot][i])) {
                        nOk++;
                        if (nOk >= Math.min(minEnergyRange, xAxis.length)) {
                            hasSpot[spot] = true;
                            nGoodSpots++;
                            break;
                        }
                    } else if (minEnergyRange==0 && !Double.isNaN(data[DELTAX][spot][i]) && !Double.isNaN(data[INTEGRAL][spot][i]))
                        hasSpot[spot] = true;                       //measure all spots, also sub-threshold
            }
            if (debug) IJ.log("After tracking: nGoodSpots="+nGoodSpots);

            smoothSpotPositions(positionAveraging, minEnergyRange);
            if (minEnergyRange > 0)
                limitExtrapolation((int)positionAveraging);

            // Eliminate spots with very bad tracking or questionable tracking and never reaching
            // reasonable significance (before they might force neigbors to disappear)
            double significanceThresholdForDeletion = minEnergyRange == 0 ? 0 : 1.5*minSignificance;
            makeDeviationStatistics(significanceThresholdForDeletion);

            double noise = Double.NaN;
            if (subtractNeighborBackground)
                noise = getNoise();
            //creation of a background stack might go here (if it will be implemented)

            IJ.showStatus("analyzing intensities...");
            significanceThresholdForDeletion = minEnergyRange == 0 ? 0 : minSignificance; //sub-threshold measurement: no min significance
            measureIV(noise, subtractNeighborBackground, i0FromBackgr, (int)Math.round(minEnergyRange), significanceThresholdForDeletion);
            nGoodSpots = 0;
            for (int spot=0; spot<spotPattern.size(); spot++)
                if (hasSpot[spot])
                    nGoodSpots++;
            if (nGoodSpots > 0) {
                long t3 = System.currentTimeMillis();
                if (debug) IJ.log("measure: "+IJ.d2s(0.001*(t3-t1))+" s");
                if (nGoodSpots >= 3) {
                    do2DRegression();
                    long t3b = System.currentTimeMillis();
                    if (debug) IJ.log("2D regression: "+IJ.d2s(0.001*(t3-t1))+" s");
                    t3 = t3b;
                }
                currentRunNumber++;                                 //newly created plots get a new number
                //position for plots and plot
                Point loc = Prefs.getLocation(ImageWindow.LOC_KEY);
                if (loc == null) loc = new Point(50,50);
                Rectangle bounds = GUI.getMaxWindowBounds(loc);
                int xSpace = PlotWindow.plotWidth + 200;
                int ySpace = PlotWindow.plotHeight + 150;           //ensure we have enough space for plots
                Point newLoc = (Point)(loc.clone());                //at the default window opening position
                if (loc.x + xSpace > bounds.x + bounds.width)
                    newLoc.x = bounds.x + 50;
                if (loc.y + ySpace > bounds.y + bounds.height)
                    newLoc.y = bounds.y + 50;
                if (!newLoc.equals(loc))
                    Prefs.saveLocation(ImageWindow.LOC_KEY, newLoc);
                Recorder.suspendRecording();
                boolean plotAll = !quickOutput;
                if (!makeQuickPlot(quickOutput)) plotAll = true;    //if quick plot fails, show all
                progressInc = (1 - progress)/PLOT_NAMES.length;
                for (int i=0; i<PLOT_NAMES.length; i++) {
                    if (plotAll || PLOT_NAMES[i].equals("QualityStat")) {
                        ImagePlus imp = getPlotImp(i);
                        if (imp==null) continue;
                        imp.show();
                        addProgress(progressInc);
                        if (Recorder.record)
                            imp.waitTillActivated();
                    }
                }
                Recorder.resumeRecording();
                long t4 = System.currentTimeMillis();
                if (debug) IJ.log("plots: "+IJ.d2s(0.001*(t4-t3))+" s");
                showOverlay(stackImp);

                long t5 = System.currentTimeMillis();
                if (debug) IJ.log("overlay: "+IJ.d2s(0.001*(t5-t4))+" s");
                IJ.showStatus("Spot tracking " + IJ.d2s(0.001*(t1-t0),1) + " s, total " + IJ.d2s(0.001*(t5-t0),1));
            } else
                IJ.error(LEED_Spot_Tracker.PLUGIN_NAME, "Sorry, cannot measure any beams.\nSpots too close? (check radius)");
            IJ.showProgress(1.0);
            boolean success = nGoodSpots > 0;
            if (success)
                spotTracker.setIVAnalyzer(this, nGoodSpots);
            prepareExit(success);
        } catch (Exception e) {
            IJ.handleException(e);
            Macro.abort();
            prepareExit(false);
        }
    }

    /** Prepare main panel before we exit */
    void prepareExit(boolean ok) {
        spotTracker.enableAndHighlightComponents(true);
        LeedDarkFlatVirtualStack ldfvStack = spotTracker.getStack();
        if (ldfvStack != null) ldfvStack.setShowProgress(true);
        if (! ok) {
            spotTracker.setButtonLabelText(LEED_Spot_Tracker.TRACK_BUTTON, "");
            spotTracker.changeStatus(0, LEED_Spot_Tracker.TRACK_OK | LEED_Spot_Tracker.SAVE_OK);
        }

    }

    boolean dataAdded;
    /** Tries to find all spots for a given slice index (0...nSlices-1).
     *  When found, writes the position, deviation form the fit, and significance (roughly SNR)
     *  to the respective fields of the data.
     *  Returns whether data have been added. */
    boolean findSpots(final int sliceI, final int lastSliceI, final double[] lastDeltaX, final double[] lastDeltaY,
            final int[] lastSeen, final double searchAgain,
            final double[][] retainedSignificance, final double[] decayingSignificance, final double significanceDecay) {
        final FloatProcessor stackIp = (FloatProcessor)stackImp.getStack().getProcessor(sliceI+1);
        final int nSpots = spotPattern.size();
        final double energy = getEnergy(sliceI);
        final double radiusInt = radiiInt[sliceI];
        final double radiusSup = radiiSup[sliceI];
        final double invSqrtEnergy = 1./Math.sqrt(energy);
        final double dlnk = energies==null ? 0 :
                Math.sqrt(energies[lastSliceI]/energies[sliceI])-1;     //change of ln(k) since last step
        final double[] stepSizeSqr = new double[nSpots];                //difference between new and previous position
        dataAdded = false;

        final AtomicInteger spotNumber = new AtomicInteger(0);          //prepare parallelization: spots are measured in parallel
        int nCPUs = Runtime.getRuntime().availableProcessors();
        if (nCPUs >=4) nCPUs = nCPUs/2;                                 //avoid too many threads (not much gain with hyperthreading, but added overhead)
        int nThreads = Math.min(nCPUs, spotPattern.size());
        final Callable[] callables = new Callable[nThreads];
        for (int t = 0; t < nThreads; t++ ) {
            callables[t] = new Callable() {
                final public Boolean call() {
                    try {
                        double[] xy = new double[4], xy2 = new double[2];
                        double[] spotData = new double[LeedSpotAnalyzer.OUTPUT_SIZE];
                        Leed2DRegression reg2D = new Leed2DRegression();
                        double[] reg2DOut = new double[2];

                        while (true) {                                  //loop over all spots:
                            int spot = spotNumber.getAndIncrement();    //find the next spot that has not been measured
                            if (spot >= nSpots) break;
                            stepSizeSqr[spot] = Double.NaN;
                            spotData[INTEGRAL] = Double.NaN;
                            spotData[PSIGNIFICANCE] = 0;
                            xy = screenFitter.screenCoordinates(spotPattern, spot, invSqrtEnergy, xy);
                            double dxCalc = xy[2] * dlnk;               //calculated movement per step
                            double dyCalc = xy[3] * dlnk;

                            if (!(data[PSIGNIFICANCE][spot][sliceI]>minSignificance)) {     //spot not found at this energy in a previous pass?
                                double radius = spotPattern.isSuperstructure(spot) ? radiusSup : radiusInt;
                                double xs = xy[0], ys = xy[1];
                                if (!Double.isNaN(lastDeltaX[spot])) {                       //have we seen this spot previously?
                                    double deltaX = lastDeltaX[spot];                        //  then correct search coordinates by last offset
                                    double deltaY = lastDeltaY[spot];
                                    double xExpected = xs + deltaX;
                                    double yExpected = ys + deltaY;
                                    if (maskRoi.contains((int)xExpected, (int)yExpected)) {  // ...and try to find it at the corrected position
                                        spotData = LeedSpotAnalyzer.centerAndAnalyzeSpot(xExpected, yExpected, maskCenterX, maskCenterY,
                                                spotBackgrShape, radius, azBlurRadians, minSignificance, stackIp, maskIp, spotData);
                                        if (spot==debugSpot) IJ.log(IJ.d2s(xAxis[sliceI],1)+" Expected "+spotPattern.getName(debugSpot)+
                                                " at "+IJ.d2s(xExpected,1)+","+IJ.d2s(yExpected,1)+". Got signif="+IJ.d2s(spotData[PSIGNIFICANCE],1)+
                                                " at "+IJ.d2s(spotData[X],1)+","+IJ.d2s(spotData[Y],1));
                                    }
                                }
                                if (!(spotData[PSIGNIFICANCE]>minSignificance) &&           // spot not found at this or at nearby energy or not at all
                                        Math.abs(sliceI - lastSeen[spot]) >= searchAgain && retainedSignificance[spot][sliceI] == 0) {
                                    boolean retryWithNoDelta = true;
                                    int nNeighbors = 0;                                     // build a local coordinate system from neighbors
                                    reg2D.clear();
                                    double sumInvDistancePxlSqr = 0;
                                    for (int spot2 : spotPattern.getNearestSpotIndices(spot)) {
                                        double significance2 = retainedSignificance[spot2][lastSliceI]; // >minSignificance for detected spots, decaying when DELTAX,Y is retained
                                        if (significance2 > NO_SIGNIFICANCE && !Double.isNaN(data[DELTAX][spot2][lastSliceI])) {
                                            nNeighbors++;
                                            xy2 = screenFitter.screenCoordinates(spotPattern, spot2, invSqrtEnergy, xy2);
                                            double invDistancePxlSqr = 1./(sqr(xy2[0]-xy[0]) + sqr(xy2[1]-xy[1])); // 1/distance^2 in pixels
                                            sumInvDistancePxlSqr += invDistancePxlSqr;
                                            double weight = significance2*invDistancePxlSqr; //weight strongly decreases with distance
                                            reg2D.addPoint(data[DELTAX][spot2][lastSliceI], data[DELTAY][spot2][lastSliceI],
                                                    spotPattern.getKx(spot2), spotPattern.getKy(spot2), weight);
                                        }
                                    }
                                    if (nNeighbors >= 3) {                                  //get deltaX, deltaY from neighbors
                                        reg2DOut = reg2D.getFit(spotPattern.getKx(spot), spotPattern.getKy(spot), reg2DOut);
                                        if (!(Double.isNaN(reg2DOut[0])) && (               //fit ok and reasonable values?
                                                sqr(reg2DOut[0]) + sqr(reg2DOut[1]) < sqr(radius*MAX_DEVIATION_RADII) ||
                                                sqr(reg2DOut[0]-lastDeltaX[spot]) + sqr(reg2DOut[1]-lastDeltaY[spot]) < sqr(radius*MAX_DEVIATION_RADII))) {
                                            double xFit = xs + reg2DOut[0];
                                            double yFit = ys + reg2DOut[1];
                                            if (maskRoi.contains((int)(xFit), (int)(yFit))) {
                                                spotData = LeedSpotAnalyzer.centerAndAnalyzeSpot(xFit, yFit, maskCenterX, maskCenterY,
                                                        spotBackgrShape, radius, azBlurRadians, minSignificance, stackIp, maskIp, spotData);
                                                if (spot==debugSpot) IJ.log(IJ.d2s(xAxis[sliceI],1)+" "+spotPattern.getName(spot)+" From "+nNeighbors+" neighbors: dx,y="+
                                                        IJ.d2s(reg2DOut[0],1)+","+IJ.d2s(reg2DOut[1],1)+" x,y="+IJ.d2s(xFit,1)+","+IJ.d2s(yFit,1)+". Got signif="+IJ.d2s(spotData[PSIGNIFICANCE],1)+
                                                        " at "+IJ.d2s(spotData[X],1)+","+IJ.d2s(spotData[Y],1)+" Integral="+LeedUtils.d2s(spotData[INTEGRAL],5));
                                                if (!(spotData[PSIGNIFICANCE]>minSignificance) &&
                                                        !Double.isNaN(spotData[INTEGRAL]) && !(decayingSignificance[spot]>=NO_SIGNIFICANCE) &&
                                                        retainedSignificance[spot][sliceI] == 0) { //if not found but inside and we have a position from neighbors:
                                                data[DELTAX][spot][sliceI] = reg2DOut[0];
                                                data[DELTAY][spot][sliceI] = reg2DOut[1];   // remember the positions from neighbors, better than nothing
                                                data[PSIGNIFICANCE][spot][sliceI] = LOW_SIGNIFICANCE; //but low significance

                                                if (sumInvDistancePxlSqr * maskRoiRect.width*maskRoiRect.height > //this is ~nNeighbors /sqr(typ. distance as fraction of screen)
                                                        3*sqr(4))                           //result ok and neighbors not too far (e.g. 3 neighbors <1/4 screen)
                                                    retryWithNoDelta = false;               //estimate from neighbors should be ok, don't retry without this correction
                                                }
                                            } //if (maskRoi.contains(xFit, yFit))
                                        } //if reg2D reasonable
                                    } //if (nNeighbors >= 3)
                                    if (!(spotData[PSIGNIFICANCE]>minSignificance) && retryWithNoDelta &&
                                            maskRoi.contains((int)xs, (int)ys)) {           //not found with delta: also try uncorrected model
                                        spotData = LeedSpotAnalyzer.centerAndAnalyzeSpot(xs, ys, maskCenterX, maskCenterY,
                                            spotBackgrShape, radius, azBlurRadians, minSignificance, stackIp, maskIp, spotData);
                                        if (spot==debugSpot) IJ.log(IJ.d2s(xAxis[sliceI],1)+" Success at uncorrected x,y="+
                                                IJ.d2s(xs,1)+","+IJ.d2s(ys,1)+": Signif="+IJ.d2s(spotData[PSIGNIFICANCE],2)+
                                                " at "+IJ.d2s(spotData[X],1)+","+IJ.d2s(spotData[Y],1));
                                    }
                                } //if spot not found for a while
                                if (spotData[PSIGNIFICANCE]>minSignificance) {              //spot found
                                    if (debug && Math.abs(sliceI - lastSeen[spot]) >= searchAgain)  //spot found first or after a long while
                                        IJ.log("Found "+spotPattern.getName(spot)+" at "+IJ.d2s(xAxis[sliceI],1)+(sliceI>lastSliceI ? "+" : "-")+
                                                " x,y="+IJ.d2s(spotData[X],1)+","+IJ.d2s(spotData[Y],1)+
                                                " signif="+IJ.d2s(spotData[PSIGNIFICANCE])+" integr="+(float)spotData[INTEGRAL]+
                                                (lastSeen[spot] <0 ? "" : " (last@"+IJ.d2s(xAxis[lastSeen[spot]])+")"));
                                    data[LeedSpotAnalyzer.PSIGNIFICANCE][spot][sliceI] = spotData[LeedSpotAnalyzer.PSIGNIFICANCE];
                                    double deltaX = spotData[X] - xy[0];                    //deviation from prediction
                                    double deltaY = spotData[Y] - xy[1];
                                    stepSizeSqr[spot] = sqr(deltaX - lastDeltaX[spot]) + sqr(deltaY - lastDeltaY[spot]); //NaN if new spot
                                    if (lastSeen[spot] < 0 || Math.abs(sliceI - lastSeen[spot]) >= searchAgain ||
                                            !(stepSizeSqr[spot] > sqr(radius))) {           //accept unless a large jump (within 'searchAgain' eV)
                                        dataAdded = true;
                                        hasSpot[spot] = true;
                                        lastDeltaX[spot] = deltaX;
                                        lastDeltaY[spot] = deltaY;
                                        lastSeen[spot] = sliceI;
                                        data[DELTAX][spot][sliceI] = deltaX;
                                        data[DELTAY][spot][sliceI] = deltaY;
                                        if (lastSliceI != sliceI && xAxis==energies) {      //going to check whether it is a stationary spot
                                            
                                            double drCalc = Math.sqrt(sqr(dxCalc) + sqr(dyCalc));   //calculated movement per step
                                            int stepsPerR = (int)Math.round((radius+2)/drCalc);     //steps to move one integration radius+2pxl; we look back that far
                                            if (stepsPerR < 1) stepsPerR = 1;
                                            int oldSlice = sliceI > lastSliceI ? sliceI - stepsPerR : sliceI + stepsPerR;
                                            if (oldSlice >= 0 && oldSlice < xAxis.length && !(data[PSIGNIFICANCE][spot][oldSlice]>minSignificance)) { //if we should have moved enough to check
                                                double dx = spotData[X] - data[X][spot][oldSlice];  //how much the positions move
                                                double dy = spotData[Y] - data[Y][spot][oldSlice];
                                                double ddx = data[DELTAX][spot][sliceI] - data[DELTAX][spot][oldSlice]; //how much the deviations change
                                                double ddy = data[DELTAY][spot][sliceI] - data[DELTAY][spot][oldSlice];

                                                if (sqr(ddx) + sqr(ddy) > 2*(sqr(dx) + sqr(dy))) {  //spot moves less than 70% of deltaX,Y change
                                                    if (debug || spot==debugSpot) IJ.log(spotPattern.getName(spot)+" got stuck at "+
                                                            IJ.d2s(xAxis[sliceI],1)+(sliceI>lastSliceI ? "+" : "-")+
                                                            " x,y="+IJ.d2s(xy[0]+deltaX,1)+","+IJ.d2s(xy[1]+deltaY,1)+
                                                            ". Keep only high significance back to "+IJ.d2s(xAxis[oldSlice],1)+
                                                            " (there x,y="+IJ.d2s(data[X][spot][oldSlice],1)+","+IJ.d2s(data[Y][spot][oldSlice],1)+")");
                                                    lastSeen[spot] = Integer.MIN_VALUE;
                                                    decayingSignificance[spot] = 0;
                                                    int step = sliceI > oldSlice ? 1 : -1;
                                                    for (int i=oldSlice; i!=sliceI+step; i+=step) { //go back and...
                                                        if (data[PSIGNIFICANCE][spot][i]>3*minSignificance) {   //...only keep high-significance data
                                                            lastDeltaX[spot] = data[DELTAX][spot][i];
                                                            lastDeltaY[spot] = data[DELTAY][spot][i];
                                                            lastSeen[spot] = i;
                                                        } else { //invalidate low-significance data
                                                            decayingSignificance[spot] *= significanceDecay;
                                                            data[X][spot][i] = Double.NaN;
                                                            data[Y][spot][i] = Double.NaN;
                                                            if (decayingSignificance[spot] > NO_SIGNIFICANCE) {
                                                                data[DELTAX][spot][i] = lastDeltaX[spot];
                                                                data[DELTAY][spot][i] = lastDeltaY[spot];
                                                                data[PSIGNIFICANCE][spot][i] = LOW_SIGNIFICANCE;
                                                            } else {
                                                                data[DELTAX][spot][i] = Double.NaN;
                                                                data[DELTAY][spot][i] = Double.NaN;
                                                                data[PSIGNIFICANCE][spot][i] = 0;
                                                            }
                                                        }
                                                        retainedSignificance[spot][i] = decayingSignificance[spot];
                                                    }
                                                } //endif stationary spot detected
                                            } //endif check for stationary spot
                                        }
                                    } //endif not a sudden jump
                                } else
                                    stepSizeSqr[spot] = 0;                              //not a jump: we keep the old deltaX, deltaY

                                if (!Double.isNaN(spotData[INTEGRAL])) {                //spot found or not, but expected spot position fully inside mask
                                    data[LeedSpotAnalyzer.INTEGRAL][spot][sliceI] = spotData[LeedSpotAnalyzer.INTEGRAL];
                                    if (!Double.isNaN(lastDeltaX[spot])) {               //use (new or previous) position offset for detection of spots jumping to the same site
                                        data[X][spot][sliceI] = xy[0] + lastDeltaX[spot];
                                        data[Y][spot][sliceI] = xy[1] + lastDeltaY[spot];
                                    }
                                }
                            } else { //if (data[PSIGNIFICANCE][spot][sliceI]>minSignificance)   // spot seen at this energy in previous pass, keep it
                                stepSizeSqr[spot] = sqr(lastDeltaX[spot] - data[DELTAX][spot][sliceI]) +
                                        sqr(lastDeltaY[spot] - data[DELTAY][spot][sliceI]);
                                lastDeltaX[spot] = data[DELTAX][spot][sliceI];
                                lastDeltaY[spot] = data[DELTAY][spot][sliceI];
                                lastSeen[spot] = sliceI;
                            }
                            double significance = Math.min(spotData[LeedSpotAnalyzer.PSIGNIFICANCE], 2*minSignificance);
                            if (significance <= minSignificance) significance = 0;
                            if (significance >= decayingSignificance[spot]) {               //Spot found and good significance
                                decayingSignificance[spot] = significance;
                                retainedSignificance[spot][sliceI] = significance;
                            } else {
                                decayingSignificance[spot] *= significanceDecay;
                                if (data[PSIGNIFICANCE][spot][sliceI] > NO_SIGNIFICANCE &&  //if spot not detected
                                        decayingSignificance[spot] > NO_SIGNIFICANCE &&     //and a previous detection is not too long ago
                                        decayingSignificance[spot] > retainedSignificance[spot][sliceI]) { //and better than what we got from other previous positions
                                    data[X][spot][sliceI] = xy[0] + lastDeltaX[spot];
                                    data[Y][spot][sliceI] = xy[1] + lastDeltaY[spot];
                                    retainedSignificance[spot][sliceI] = decayingSignificance[spot];
                                }
                            }                            
                            if (spot==debugSpot)
                                IJ.log((float)xAxis[sliceI]+" "+spotPattern.getName(spot)+" Result: Integral="+(float)spotData[INTEGRAL]+
                                        " sig="+IJ.d2s(data[PSIGNIFICANCE][spot][sliceI])+
                                        " delta x,y="+IJ.d2s(data[DELTAX][spot][sliceI])+","+IJ.d2s(data[DELTAY][spot][sliceI])+
                                        " x,y="+IJ.d2s(data[X][spot][sliceI])+","+IJ.d2s(data[Y][spot][sliceI]));
                        } //loop over all spots
                    } catch(Exception ex) {IJ.handleException(ex);}
                    return null;
                }
            }; //new Callable()
        }
        ThreadUtil.startAndJoin(callables);
        //check whether two spots are at the same position; if so, delete the one that has jumped here or the newer one
        //it is enough to check pairs that are neighbors in the spotPattern.
        for (int i=0; i<nSpots; i++) {
            if (Double.isNaN(data[X][i][sliceI])) continue;
            int[] neighbors = spotPattern.getNearestSpotIndices(i);
            for (int iN=0; iN<neighbors.length; iN++) {
                int j = neighbors[iN];
                if (j<i) continue;    //check each pair only once, when j>i
                if (Double.isNaN(data[X][j][sliceI])) continue;
                double dSqr = sqr(data[X][i][sliceI] - data[X][j][sliceI]) +
                        sqr(data[Y][i][sliceI] - data[Y][j][sliceI]);
                if (dSqr < sqr(radiusInt)) {            //if two spots at (almost) the same position
                    if (Double.isNaN(data[DELTAX][i][sliceI]) && Double.isNaN(data[DELTAX][j][sliceI])) continue; //both spots invisible? ignore
                    if ((i==debugSpot ||j==debugSpot)) IJ.log(xAxis[sliceI]+" Collision: "+spotPattern.getName(i)+"@("+IJ.d2s(data[X][i][sliceI],1)+","+IJ.d2s(data[Y][i][sliceI],1)+"). "+
                            spotPattern.getName(j)+"@("+IJ.d2s(data[X][j][sliceI])+","+IJ.d2s(data[Y][j][sliceI])+") d="+IJ.d2s(Math.sqrt(dSqr)));
                    boolean iJumpedOrNew = !(stepSizeSqr[i] < sqr(0.5)*sqr(radiusInt)); //new: stepSizeSqr = NaN
                    boolean jJumpedOrNew = !(stepSizeSqr[j] < sqr(0.5)*sqr(radiusInt));
                    int spotToDelete = -1;
                    if (iJumpedOrNew != jJumpedOrNew) {    //delete the one that is new or has jumped > 0.5*integration radius
                        spotToDelete = iJumpedOrNew ? i : j;
                    } else {
                        double deltaIsqr = sqr(data[DELTAX][i][sliceI]) + sqr(data[DELTAY][i][sliceI]);
                        double deltaJsqr = sqr(data[DELTAX][j][sliceI]) + sqr(data[DELTAY][j][sliceI]);
                        boolean deleteI = deltaIsqr > deltaJsqr; //delete the one further form the predicted position
                        spotToDelete = deleteI ? i : j;
                    }
                    data[X][spotToDelete][sliceI] = Double.NaN;
                    data[DELTAX][spotToDelete][sliceI] = Double.NaN;
                    lastDeltaX[spotToDelete] = Float.NaN;
                    if (debug || (i==debugSpot ||j==debugSpot)) IJ.log(IJ.d2s(xAxis[sliceI],1)+": delete spot "+spotPattern.getName(spotToDelete)
                            +" colliding with "+spotPattern.getName(i==spotToDelete ? j : i));
                }
            }
        }
        addProgress(progressInc);
        return dataAdded;
    }

    /** Smooths the spot positions by fitting a linear function (of 1/sqrt(energy)) over a range of the data,
     *  with weights increasing with the significance as determined by the LeedSpotAnalyzer.
     *  There are two passes; the first pass it also bridges gaps where we have no positions */
    void smoothSpotPositions(double positionAveraging, double minEnergyRange) {
        double[] invSqrtEnergy = new double[xAxis.length];
        boolean useEnergies = xAxis == energies;
        for (int i=0; i<xAxis.length; i++)
            invSqrtEnergy[i] = useEnergies ? 1./Math.sqrt(energies[i]) : i;   //we fit spot position deviations over 1/sqrt(E) if energy is the x axis

        double[] weight1  = new double[xAxis.length];   //weights for 1st linear regression pass
        double[] input2x  = new double[xAxis.length];   //output of 1st and input for 2nd linear regression pass, x
        double[] input2y  = new double[xAxis.length];   //output of 1st and input for 2nd linear regression pass, y
        double[] weight2x = new double[xAxis.length];   //weights for 2nd linear regression pass, x
        double[] weight2y = new double[xAxis.length];   //weights for 2nd linear regression pass, y
        LeedLinearRegression xLine = new LeedLinearRegression();
        LeedLinearRegression yLine = new LeedLinearRegression();
        int pointsPerSide = (int)Math.abs(0.5*positionAveraging);
        if (pointsPerSide < 2) pointsPerSide = 2;       //how many fit points we need at each side of the current point (does not work reasonably with 1)
        final double no_weight = significanceToWeight(NO_SIGNIFICANCE); //measured positions have at least this weight
        final double e_infty_weight = 0.1;              //weight of extra point at energy->infinity
        boolean useSubThreshold = minEnergyRange==0;

        for (int spot=0; spot<spotPattern.size(); spot++) {
            if (!hasSpot[spot]) continue;
            // first smoothing pass
            int lastEnergyI = 0;                        //index of last energy with significance>0
            for (int i=0; i<xAxis.length; i++) {
                weight1[i] = Double.isNaN(data[DELTAX][spot][i]) || Double.isNaN(data[INTEGRAL][spot][i]) ? 0 :
                        significanceToWeight(data[LeedSpotAnalyzer.PSIGNIFICANCE][spot][i]);
                if (!useSubThreshold && weight1[i] == LOW_WEIGHT)
                    weight1[i] = 0;                     //don't use positions from neighbors unless sub-threshold is enabled
                if (weight1[i] != 0) lastEnergyI = i;
            }
            xLine.clear();
            yLine.clear();
            if (useEnergies) {
                
                xLine.addPointW(0, 0, e_infty_weight);  //add a point at energy->infty to avoid large slopes
                yLine.addPointW(0, 0, e_infty_weight);
            }
            int nSpotsBefore = 0, nSpotsAtOrAfter = 0;
            int first = -1, last = -1;                  //first and last point currently in the fit
            int nDataPoints = 0;
            for (int i=0; i<xAxis.length; i++) {
if (spot==debugSpot)IJ.log("i="+i+" weight1="+weight1[i]+" atOrAfter="+nSpotsAtOrAfter+" lastEnergyI="+lastEnergyI);
                if (nSpotsAtOrAfter > 0 && weight1[i] != 0) {
                    if (weight1[i] > LOW_WEIGHT) {      //also for sub-threshold, we count only measured points (i.e., > LOW_WEIGHT)
                        nSpotsBefore++;                 //current point becomes a 'before' point
                        nSpotsAtOrAfter--;
                    }
                    while (nSpotsBefore > pointsPerSide) { //remove old point(s) from the fit
                        xLine.addPointW(invSqrtEnergy[first], data[DELTAX][spot][first], -weight1[first]);
                        yLine.addPointW(invSqrtEnergy[first], data[DELTAY][spot][first], -weight1[first]);
                        nDataPoints--;
                        if (weight1[first] > LOW_WEIGHT)
                            nSpotsBefore--;
                        for (int j=first+1; j<i; j++) { //find the next 'old' point
                            if (weight1[j] != 0) {
                                first = j;
                                break;
                            }
                        }
                    }
                }
                int neededSpotsAtOrAfter = pointsPerSide;
                if (weight1[i] != 0) neededSpotsAtOrAfter++;
                for (int j=last+1; (nSpotsAtOrAfter < neededSpotsAtOrAfter) && j<=lastEnergyI; j++) {
                    if (weight1[j] != 0) {              //add this point to the fit
                        xLine.addPointW(invSqrtEnergy[j], data[DELTAX][spot][j], weight1[j]);
                        yLine.addPointW(invSqrtEnergy[j], data[DELTAY][spot][j], weight1[j]);
                        nDataPoints++;
                        last = j;
                        if (weight1[j] > LOW_WEIGHT) {  //also for sub-threshold, we count only measured points (i.e., > LOW_WEIGHT)
                            if (j < i)
                                nSpotsBefore++;
                            else
                                nSpotsAtOrAfter++;
                        }
                        if (first < 0) first = j;
                    }
                }
                if (useEnergies) {
                    input2x[i] = xLine.getFitValue(invSqrtEnergy[i]);
                    input2y[i] = yLine.getFitValue(invSqrtEnergy[i]);
                } else {                                    //if we have no energies as x axis: avoid large values via linear extrapolation, take the lowest reasonable slope
                    input2x[i] = xLine.getFitValueWithMinSlope(invSqrtEnergy[i], /*yErrSqr=*/ 4, nDataPoints);
                    input2y[i] = yLine.getFitValueWithMinSlope(invSqrtEnergy[i], /*yErrSqr=*/ 4, nDataPoints);
                }
                boolean enoughData = xLine.counter >= 2*no_weight || useSubThreshold;   //require at least two measured data points
                weight2x[i] = enoughData ? xLine.getFitWeight(invSqrtEnergy[i]) : 1e-30;
                weight2y[i] = enoughData ? yLine.getFitWeight(invSqrtEnergy[i]) : 1e-30;
                if (spot==debugSpot && (i==xAxis.length/8 || i==(int)Math.round(xAxis.length*0.67) || i==1701))
                    IJ.log(spotPattern.getName(spot)+" 1stSmoothing@"+IJ.d2s(xAxis[i],1)+" Before,After="+nSpotsBefore+","+nSpotsAtOrAfter+
                            "/"+pointsPerSide+" rangeUsed="+IJ.d2s(xAxis[Math.max(first,0)],1)+","+IJ.d2s(xAxis[Math.min(last,xAxis.length-1)],1)+
                            "\n ->x,y="+IJ.d2s(input2x[i])+","+IJ.d2s(input2y[i])+" weight2x,y="+LeedUtils.d2s(weight2x[i],3)+","+LeedUtils.d2s(weight2y[i],3)+
                            " fitSumWx,y="+(float)xLine.counter+","+(float)yLine.counter);
            }

            // 2nd smoothing pass and calculate corrected coordinates
            xLine.clear();
            yLine.clear();
            double[] xy = new double[2];

            for (int iAdd=0; iAdd<Math.min(pointsPerSide, xAxis.length); iAdd++) {
                xLine.addPointW(invSqrtEnergy[iAdd], input2x[iAdd], weight2x[iAdd]);
                yLine.addPointW(invSqrtEnergy[iAdd], input2y[iAdd], weight2y[iAdd]);
            }
            for (int i=0; i<xAxis.length; i++) {
                int iAdd = i + pointsPerSide;
                if (iAdd < xAxis.length) {
                    xLine.addPointW(invSqrtEnergy[iAdd], input2x[iAdd], weight2x[iAdd]);
                    yLine.addPointW(invSqrtEnergy[iAdd], input2y[iAdd], weight2y[iAdd]);
                }
                int iRemove = i - pointsPerSide - 1;
                if (iRemove >= 0) {
                    xLine.addPointW(invSqrtEnergy[iRemove], input2x[iRemove], -weight2x[iRemove]);
                    yLine.addPointW(invSqrtEnergy[iRemove], input2y[iRemove], -weight2y[iRemove]);
                }
                data[DELTAXS][spot][i] = xLine.getFitValue(invSqrtEnergy[i]);
                data[DELTAYS][spot][i] = yLine.getFitValue(invSqrtEnergy[i]);

                xy = screenFitter.screenCoordinates(spotPattern, spot, useEnergies ? invSqrtEnergy[i] : 1./Math.sqrt(getEnergy(i)), xy);
                data[X][spot][i] = xy[0] + data[DELTAXS][spot][i];
                data[Y][spot][i] = xy[1] + data[DELTAYS][spot][i];
                if (spot==debugSpot && (i==xAxis.length/8 || i==xAxis.length*2/3 || i==1701))
                    IJ.log(spotPattern.getName(spot)+" 2ndSmoothing@"+IJ.d2s(xAxis[i],1)+
                        " rangeUsed="+IJ.d2s(xAxis[Math.max(iRemove,0)],1)+","+IJ.d2s(xAxis[Math.min(iAdd,xAxis.length-1)],1)+
                        " dx,dy="+IJ.d2s(data[DELTAXS][spot][i],1)+","+IJ.d2s(data[DELTAYS][spot][i],1));
            }
        }
    }

    /** Weights for the first pass of smoothing/interpolating spot position deviations */
    double significanceToWeight(double significance) {
        final double noSignificance = 0.8*minSignificance;
        final double highSignificance = 5; //weight levels off at high significance, weight is 50% of max at noSig.+hiSig.
        if (Double.isNaN(significance)) return 0;
        if (significance == LOW_SIGNIFICANCE) return LOW_WEIGHT; //weight for spot positions from neighbors
        significance -= noSignificance;
        if (significance <=0) return 0;
        return LOW_WEIGHT + significance/(significance + highSignificance);
    }

    /** Avoids extrapolating positions further than 'maxExtraolate' from the first/last point
     *  where a spot has been detected (i.e., with >minSignificance).
     *  Sets spot screen positions (for later intensity analysis) and the smoothed position
     *  deviations (for plotting) to NaN if the point is outside the allowed extrapolation
     *  range. Sets hasSpot=false for fully insignificant beams. */
    void limitExtrapolation(int maxExtrapolate) {
        for (int spot=0; spot<spotPattern.size(); spot++) {
            int firstSeen = -1;
            int lastSeen = 0;
            for (int i=0; i<xAxis.length; i++) {
                if (data[PSIGNIFICANCE][spot][i] > minSignificance) {
                    if (firstSeen < 0) firstSeen = i;
                    lastSeen = i;
                }
            }
            if (firstSeen < 0) {    //spot never detected
                hasSpot[spot] = false;
                if (debug || spot == debugSpot) IJ.log(spotPattern.getName(spot)+" never above significance limit, deleted");
            } else {
                if (spot == debugSpot) IJ.log(spotPattern.getName(spot)+" above noise "+IJ.d2s(xAxis[firstSeen])+"-"+IJ.d2s(xAxis[lastSeen]));
                if (firstSeen-maxExtrapolate > 0) {
                    for (int i=0; i<firstSeen-maxExtrapolate; i++) {
                        data[X][spot][i] = Double.NaN;
                        data[Y][spot][i] = Double.NaN;
                        data[DELTAXS][spot][i] = Double.NaN;
                        data[DELTAYS][spot][i] = Double.NaN;
                    }
                    if (spot == debugSpot) IJ.log(spotPattern.getName(spot)+" delete extrapolation <= "+IJ.d2s(xAxis[firstSeen-maxExtrapolate-1],1));
                }
                if (lastSeen < xAxis.length - (maxExtrapolate+1)) {
                    for (int i=lastSeen+(maxExtrapolate+1); i<xAxis.length; i++) {
                        data[X][spot][i] = Double.NaN;
                        data[Y][spot][i] = Double.NaN;
                        data[DELTAXS][spot][i] = Double.NaN;
                        data[DELTAYS][spot][i] = Double.NaN;
                    }
                if (spot == debugSpot) IJ.log(spotPattern.getName(spot)+" delete extrapolation >= "+IJ.d2s(xAxis[lastSeen+(maxExtrapolate+1)],1));
                }
            }
        } //for spot
    }

    /** Returns an estimate of the image noise, taken from the first frame */
    double getNoise() {
        int iEnergyForNoise = 0;  //we estimate the background noise for the first image
        return getNoise((FloatProcessor)stackImp.getStack().getProcessor(iEnergyForNoise+1), maskIp, iEnergyForNoise);
    }

    /** Gets an estimate of the average noise of the background.
     *  This is done by high-pass filtering, i.e. convolution with
     *  the kernel (-1, 2, -1) along a line */
    double getNoise(FloatProcessor fp, ByteProcessor maskIp, int iEnergy) {
        ByteProcessor maskMinusSpots = makeMaskMinusSpots(maskIp, iEnergy, null, /*radiusFactor=*/2.0);
        int width = maskIp.getWidth(), height = maskIp.getHeight();
        byte[] bPixels = (byte[])maskMinusSpots.getPixels();
        float[] fPixels = (float[])fp.getPixels();
        double sumSqr = 0;
        int n=0;
        for (int y=0, p=0; y<height; y++) {
            float last = Float.NaN, preLast = Float.NaN;
            for (int x=0; x<width; x++,p++) {
                float newValue = bPixels[p] != 0 ? fPixels[p] : Float.NaN;
                float filtered = 2*last-(preLast+newValue);
                if (!Float.isNaN(filtered)) {
                    sumSqr += filtered*filtered;
                    n++;
                }
                preLast = last; last = newValue;
            }
        }
        double noise = Math.sqrt(sumSqr/(6*n)); //'6' because of kernel: 1^2 + 2^2 +1^2
        //IJ.log("noise="+(float)noise);
        return noise;
    }

    /** Gets the actual I(V) curves in parallel threads (handling different energies).
     *  If desired, also subtracts a 1/r^2 background of bright spots and measures the overall
     *  background intensity for I0 correction
     *  Parallelization: slices are processed independently */
    void measureIV(final double noise, final boolean subtractNeighborBackground,
            final boolean i0FromBackgr, final int minEnergyRange, double significanceThresholdForDeletion) {
        int[] nSpotDataPoints = new int[spotPattern.size()];
        final LeedAtomicFloatArray maxEnergy = new LeedAtomicFloatArray(spotPattern.size(), Float.MAX_VALUE);
        final boolean[] hasCollision = new boolean[spotPattern.size()];
        final ImageStack bgsubDebugStack = debug && subtractNeighborBackground ?
                new ImageStack(stackImp.getWidth(), stackImp.getHeight(), xAxis.length) : null;
        final AtomicInteger energyI = new AtomicInteger(0);                //prepare parallelization
        int nThreads = Math.min(Runtime.getRuntime().availableProcessors(), xAxis.length);
        final Callable[] callables = new Callable[nThreads];

        for (int t = 0; t < nThreads; t++ ) {
            callables[t] = new Callable() {
                final public Boolean call() {
                    try {
                        double[] spotData = new double[LeedSpotAnalyzer.OUTPUT_SIZE];
                        while (true) {  //loop over different energies within a thread
                            int sliceI = energyI.getAndIncrement(); //starts with 0 (i.e., slice number = sliceI+1)
                            if (sliceI >= xAxis.length) break;
                            double energy = getEnergy(sliceI);
                            double radiusInt = radiiInt[sliceI];
                            double radiusSup = radiiSup[sliceI];
                            FloatProcessor stackIp = (FloatProcessor)stackImp.getStack().getProcessor(sliceI+1);
                            ByteProcessor maskMinusSpots = null;
                            for (int spot=0; spot<spotPattern.size(); spot++) {
                                if (!hasSpot[spot]) continue;
                                double xs = data[X][spot][sliceI];
                                double ys = data[Y][spot][sliceI];
                                spotData = LeedSpotAnalyzer.analyzeSpot(xs, ys, maskCenterX, maskCenterY,
                                        spotBackgrShape, spotPattern.isSuperstructure(spot) ? radiusSup : radiusInt, azBlurRadians,
                                        minSignificance, stackIp, maskIp, spotData);
                                for (int d=Y+1; d<spotData.length; d++)
                                    data[d][spot][sliceI] = spotData[d];
                                if (Double.isNaN(spotData[INTEGRAL])) {    //spot outside mask?
                                    data[DELTAXS][spot][sliceI] = Double.NaN;                    //then don't show smoothed positions
                                    data[DELTAYS][spot][sliceI] = Double.NaN;
                                }
                            }
                            // Check for nearby spots. Cutoff is 2.5 integration radii
                            // With circular background, no overlap between background ring and core of neighbor spot
                            // (limit would be 1+sqrt2 integration radii)
                            // For angle blur, the integration area has a semimajor axis of 2 radii, so the integration
                            // areas overlap at a distance of 4 radii in tangential direction.
                            double r0 = spotPattern.isSuperstructure() ? Math.max(radiusSup, radiusInt) : radiusInt;
                            double minDistSqr = sqr(r0*2.5);
                            for (int spot=0; spot<spotPattern.size(); spot++) {
                                if (!hasSpot[spot]) continue;
                                if (Double.isNaN(data[INTEGRAL][spot][sliceI])) continue; //ignore spots if outside
                                int[] neighborSpots = spotPattern.getNearestSpotIndices(spot);
                                for (int j=0; j<Math.min(neighborSpots.length, 12); j++) {
                                    int spot2 = neighborSpots[j];
                                    if (!hasSpot[spot2]) continue;
                                    if (spot2 > spot) continue; //need to check each pair only once: if spot2 < spot (then spot2 has been measured already)
                                    if (Double.isNaN(data[INTEGRAL][spot2][sliceI])) continue;
                                    double distSqr = sqr(data[X][spot][sliceI] - data[X][spot2][sliceI]) +
                                            sqr(data[Y][spot][sliceI] - data[Y][spot2][sliceI]);
                                     if (distSqr < minDistSqr) {
                                        if (!maskRoi.contains((int)data[X][spot][sliceI], (int)data[Y][spot][sliceI]) ||
                                                !maskRoi.contains((int)data[X][spot2][sliceI], (int)data[Y][spot2][sliceI]))
                                            continue;         //ignore nearby spot out of screen (position my be very inaccurate, e.g. high polynomial order)
                                        if (spot==debugSpot||spot2==debugSpot)
                                            IJ.log(IJ.d2s(xAxis[sliceI])+": Collision: "+spotPattern.getNameWithGroup(spot,false)+"&"+spotPattern.getNameWithGroup(spot2,false)+" r="+IJ.d2s(Math.sqrt(distSqr))+"<"+IJ.d2s(Math.sqrt(minDistSqr)));
                                        if (data[PSIGNIFICANCE][spot2][sliceI] > minSignificance ||              //if the nearby neighbor is above the threshold
                                                data[INTEGRAL][spot2][sliceI] > 0.5*data[INTEGRAL][spot][sliceI]) { //or >50% of the current spot
                                            if (xAxis == energies) {
                                                maxEnergy.setMinimumValue(spot, (float)energy);             //...intensity data are invalid
                                            } else {
                                                data[INTEGRAL][spot][sliceI] = Double.NaN;
                                                hasCollision[spot] = true;
                                            }
                                        }
                                        if (data[PSIGNIFICANCE][spot][sliceI] > minSignificance ||
                                                data[INTEGRAL][spot][sliceI] > 0.5*data[INTEGRAL][spot2][sliceI]) {
                                            if (xAxis == energies) {
                                                maxEnergy.setMinimumValue(spot2, (float)energy);
                                            } else {
                                                data[INTEGRAL][spot2][sliceI] = Double.NaN; //spot2 has been measured already, so it remains NaN
                                                hasCollision[spot2] = true;
                                            }
                                        }
                                    }
                                }
                            }
                            FloatProcessor corrFp = stackIp;        //the same or after 1/r^2 background subtraction
                            //subtract 1/r^2 background emanating from bright spots within ring r, 2r
                            if (subtractNeighborBackground) {
                                boolean[] isBrightSpot = new boolean[spotPattern.size()];
                                corrFp = (FloatProcessor)stackIp.duplicate();
                                float[] fPixels = (float[])corrFp.getPixels();
                                int width = stackIp.getWidth(), height = stackIp.getHeight();
                                int nBrightSpots = 0;
                                double maxIntegral = 0;
                                for (int spot=0; spot<spotPattern.size(); spot++) {         //make list of bright spots
                                    if (!hasSpot[spot]) continue;
                                    double radius = spotPattern.isSuperstructure(spot) ? radiusSup : radiusInt;
                                    if (!(data[INTEGRAL][spot][sliceI] > Math.PI*10*noise*radius*radius)) continue; //bright spots: average at least 10 stddev above noise
                                    if (!(data[LeedSpotAnalyzer.PSIGNIFICANCE][spot][sliceI] > minSignificance)) continue; //bright spots must have good significance
                                    isBrightSpot[spot] = true;
                                    if (data[INTEGRAL][spot][sliceI] > maxIntegral) maxIntegral = data[INTEGRAL][spot][sliceI];
                                    nBrightSpots++;
                                }
                                maskMinusSpots = makeMaskMinusSpots(maskIp, sliceI, isBrightSpot, 1.0); //mask with spot circles=0
                                byte[] maskPixels = (byte[])maskMinusSpots.getPixels();
                                int[] spotIndices = new int[nBrightSpots];
                                int[] spotIntegrals = new int[nBrightSpots];
                                int n=0;
                                for (int spot=0; spot<spotPattern.size(); spot++)
                                    if (isBrightSpot[spot]) {
                                        if (!hasSpot[spot]) continue;
                                        if (!(data[INTEGRAL][spot][sliceI] > 0)) continue;
                                        if (data[INTEGRAL][spot][sliceI] < MIN_REL_INTEGRAL_FOR_BG*maxIntegral) continue;
                                        spotIndices[n] = spot;
                                        spotIntegrals[n] = (int)(data[INTEGRAL][spot][sliceI]*(Integer.MAX_VALUE/maxIntegral));
                                        n++;
                                    }
                                int[] sortedIndices = LeedUtils.sortedIndices(spotIndices, spotIntegrals, n);
                                LeedLinearRegression regression = new LeedLinearRegression();
                                for (n=sortedIndices.length-1; n>=0; n--) {                 //starting from brightest spot
                                    int spot = sortedIndices[n];
                                    double xC = data[X][spot][sliceI];
                                    double yC = data[Y][spot][sliceI];
                                    if (n!=sortedIndices.length-1) {                        //measure spot again after background of others has been subtracted (except brightest)
                                        spotData = LeedSpotAnalyzer.analyzeSpot(xC, yC, maskCenterX, maskCenterY,
                                                spotBackgrShape, spotPattern.isSuperstructure(spot) ? radiusSup : radiusInt, azBlurRadians,
                                                minSignificance, corrFp, maskIp, spotData);
                                        for (int d=Y+1; d<spotData.length; d++)
                                            data[d][spot][sliceI] = spotData[d];
                                        xC = spotData[X];                  //use refined positions from last measurement
                                        yC = spotData[Y];
                                    }
                                    regression.clear();                                     //prepare 1/r^2 intensity fit
                                    double radius = spotPattern.isSuperstructure(spot) ? radiusSup : radiusInt;
                                    for (int y=Math.max(0, (int)Math.ceil(yC-2*radius));    //within a ring of outer radius 2r
                                            y<=Math.min(height-1, (int)Math.floor(yC+2*radius)); y++) {
                                        double dx = Math.sqrt(sqr(2*radius) - sqr(y-yC));
                                        for (int x=Math.max(0, (int)Math.round(xC-dx)), p=y*width+x;
                                                x<=Math.min(width-1, (int)Math.round(xC+dx)); x++,p++) {
                                            if (maskPixels[p] == 0) continue;               //only within mask and not within a spot circle
                                            double rSquare = sqr(x-xC) + sqr(y-yC);
                                            regression.addPoint(1./rSquare, fPixels[p]);    //fit: constant+1/r^2
                                        }
                                    }

                                    if (regression.getCounter() < 8 || regression.getMeanDX2() < 0.018/sqr(sqr(radius)))
                                        continue;       //full ring would give us MeanDX2 = 0.037/r^4
                                    double slope = regression.getSlope();
                                    Rectangle rect = maskRoi.getBounds();
                                    double base = 0.1*noise;                                //we stop where a/r^2 would be base
                                    double rBase = Math.sqrt(slope/base);                   //radius where a/r^2 = base
                                    int ystart = Math.max(rect.y, (int)(yC-rBase));
                                    int yend = Math.min(rect.y+rect.height, (int)(Math.ceil(yC+rBase))+1);
                                    for (int y=ystart; y<yend; y++) {
                                        double dx = Math.sqrt(sqr(rBase) - sqr(y-yC));
                                        int xstart = Math.max(rect.x, (int)(xC-rBase));
                                        int xend = Math.min(rect.x+rect.width, (int)Math.ceil(xC+rBase)+1);
                                        for (int x=xstart, p=y*width+x; x<xend; x++,p++) {
                                            double rSqr = sqr(x-xC) + sqr(y-yC);
                                            double subtractMe = slope/rSqr - base;
                                            if (subtractMe > 0 && rSqr > 0.25*radius*radius)
                                                fPixels[p] -= (float)subtractMe;            //subtract the fit except for the very center where it diverges
                                            }
                                    }
                                }
                                for (int spot=0; spot<spotPattern.size(); spot++) {         //re-measure the remaining (not bright) spots
                                    if (!hasSpot[spot] || isBrightSpot[spot]) continue;
                                    if (Double.isNaN(data[INTEGRAL][spot][sliceI])) continue;    //ignore if outside or deleted due to collision
                                    double xs = data[X][spot][sliceI];
                                    double ys = data[Y][spot][sliceI];
                                    spotData = LeedSpotAnalyzer.analyzeSpot(xs, ys, maskCenterX, maskCenterY,
                                            spotBackgrShape, spotPattern.isSuperstructure(spot) ? radiusSup : radiusInt, azBlurRadians,
                                            minSignificance, corrFp, maskIp, spotData);
                                    for (int d=Y+1; d<spotData.length; d++)
                                        data[d][spot][sliceI] = spotData[d];
                                }
                                if (bgsubDebugStack != null) bgsubDebugStack.setProcessor(corrFp, sliceI+1);

                            } //if (subtractNeighborBackground)
                            backgroundIntensities[sliceI] = getBackgroundIntensity(stackIp, maskMinusSpots, maskIp, sliceI);
                            addProgress(progressInc);
                        }
                    } catch(Exception ex) {IJ.handleException(ex);}
                    return null;
                }
            }; //new Callable()
        }
        ThreadUtil.startAndJoin(callables);
        if (bgsubDebugStack != null) new ImagePlus("Background subtracted", bgsubDebugStack).show();

        // Eliminate data above maxEnergy (where spots are too close)
        nTooClose = 0;
        if (xAxis == energies) {
            for (int spot=0; spot<spotPattern.size(); spot++) {
                if (!hasSpot[spot]) continue;
                double maxE = maxEnergy.get(spot);
                if (maxE == Float.MAX_VALUE) continue;
                nTooClose++;
                if (debug) IJ.log(spotPattern.getName(spot)+" too close to neighbor above "+IJ.d2s(maxE, 1));
                for (int i=energies.length-1; i>=0; i--) {
                    if (energies[i] < maxE) break;
                    data[INTEGRAL][spot][i] = Double.NaN;
                }
            }
        } else { // count the colliding spots if the image doesn't scale with 1/sqrt(E)
            for (int spot=0; spot < spotPattern.size(); spot++)
                if (hasCollision[spot])
                    nTooClose++;
        }
        if (debug && nTooClose>0) IJ.log(nTooClose+" beams with neighbors too close at high E (deleted at high E)");

        // I0-related normalization factor
        double[] inverseI0 = null;
        if (processedI0 != null) {
            if (i0FromBackgr)
                LeedEnergyI0Selector.applyI0BackgroundVariations(processedI0, backgroundIntensities);
            inverseI0 = new double[xAxis.length];
            for (int i=0; i<xAxis.length; i++)
                inverseI0[i] = 1.0/processedI0[i];
        }
        //eliminate if too few adjacent points or significance too low. Also find overall maximum intensity
        int minPointsPerBeam = Math.min(minEnergyRange, xAxis.length/2);
        if (minPointsPerBeam <= 0) minPointsPerBeam = 1;
        this.minPointsPerBeam = minPointsPerBeam;
        double maxIntensity = 0;
        int nTooFewPoints = 0;
        for (int spot=0; spot<spotPattern.size(); spot++) {
            if (!hasSpot[spot]) continue;
            int nDataPoints = 0;        // number of points in largest contiguous non-NaN section
            int currentDataPoints = 0;
            double currentMax = 0;
            double maxSignificance = 0;
            for (int i=0; i<xAxis.length; i++) {
                double integral = data[INTEGRAL][spot][i];
                if (!Double.isNaN(integral)) {
                    currentDataPoints++;
                    if (integral > highestIntensity) {
                        highestIntensity = integral;
                        energyOfHighestI = xAxis[i];
                        spotOfHighestI = spot;
                    }
                    if (data[LeedSpotAnalyzer.PSIGNIFICANCE][spot][i] > maxSignificance)
                        maxSignificance = data[LeedSpotAnalyzer.PSIGNIFICANCE][spot][i];
                    double value = integral;
                    if (processedI0 != null)
                        value *= inverseI0[i];
                    if (value > currentMax)
                        currentMax = value;
                } else { //NaN encountered
                    if (currentDataPoints > nDataPoints) nDataPoints = currentDataPoints;
                    currentDataPoints = 0;
                }
            }
            if (currentDataPoints > nDataPoints) nDataPoints = currentDataPoints;
            if (nDataPoints < minPointsPerBeam || maxSignificance < significanceThresholdForDeletion) {
                hasSpot[spot] = false;
                nTooFewPoints++;
                if (spot == debugSpot) IJ.log(spotPattern.getName(spot)+" not enough contiguous data points: "+nDataPoints);
            } else
                if (currentMax > maxIntensity) maxIntensity = currentMax;
            if (spot==debugSpot) IJ.log(spotPattern.getName(spot)+" contiguous data points: "+nDataPoints+", max Significance: "+IJ.d2s(maxSignificance)+(hasSpot[spot] ? " -> kept" : " deleted"));
        }
        if (debug && nTooFewPoints>0) IJ.log(nTooFewPoints+" beams with too few points (<"+minPointsPerBeam+")");
        //create I0-corrected intensities. These are normalized to 1000 for the highest maximum
        if (inverseI0 != null)
            for (int i=0; i<xAxis.length; i++)
                inverseI0[i] *= 1000./maxIntensity;

        for (int spot=0; spot<spotPattern.size(); spot++) {
            if (!hasSpot[spot]) continue;
            for (int i=0; i<xAxis.length; i++)
                data[INT_I0CORRECTED][spot][i] = data[INTEGRAL][spot][i] *
                        (inverseI0 == null ? 1000./maxIntensity : inverseI0[i]);
        }
    }

    /** Returns a binary image equivalent to the mask, but with pixel values of 0 inside
     *  the a circle around the spots, where the radius of the circle is given by the
     *  integration radius multiplied with 'radiusFactor'.
     *  If isBrightSpots is non-null, uses only the spots where the respective array
     *  element of 'isBrightSpots' is true. */
    ByteProcessor makeMaskMinusSpots(ByteProcessor maskIp, int iEnergy, boolean[] isBrightSpot, double radiusFactor) {
        int width = maskIp.getWidth(), height = maskIp.getHeight();
        ByteProcessor mask2 = (ByteProcessor)maskIp.duplicate();
        byte[] pixels = (byte[])mask2.getPixels();
        for (int spot=0; spot<spotPattern.size(); spot++) {
            if (!hasSpot[spot] || (isBrightSpot != null && !isBrightSpot[spot])) continue;
            double radius = spotPattern.isSuperstructure(spot) ?
                    radiusFactor*radiiSup[iEnergy] : radiusFactor*radiiInt[iEnergy];
            double xC = data[X][spot][iEnergy], yC = data[Y][spot][iEnergy];
            for (int y=Math.max(0, (int)Math.ceil(yC-radius));  //within a circle of radius r
                    y<=Math.min(height-1, (int)Math.floor(yC+radius)); y++) {
                double dx = Math.sqrt(sqr(radius) - sqr(y-yC));
                for (int x=Math.max(0, (int)Math.round(xC-dx)), p=y*width+x;
                        x<=Math.min(width-1, (int)Math.round(xC+dx)); x++,p++)
                    pixels[p] = (byte)0;
            }
        }
        return mask2;
    }

    /** Gets the background intensity.
     *  We take the average intensity of the darkest 40% of the pixels, with the spot integration areas
     *  excluded (unless spots are very dense, this does not make much difference).
     *  Note that taking a lower percentage would lead leads to slightly larger variations since
     *  the background is usually darker near the outer edge, and the average intensity there
     *  can be influenced more strongly by a few bright spots than the average over a larger area.
     */
    double getBackgroundIntensity(FloatProcessor stackIp, ByteProcessor maskMinusSpots, ByteProcessor maskIp, int sliceI) {
        if (maskMinusSpots == null)
            maskMinusSpots = makeMaskMinusSpots(maskIp, sliceI, /*isBrightSpot=*/null, 1.0); //mask with spot circles=0
        byte[] maskPixels = (byte[])maskMinusSpots.getPixels();
        float[] fPixels = (float[])stackIp.getPixels();
        final double darkFraction = 0.4;
        int count=0;
        double sum=0, sumSqr=0;
        for (int p=0; p<fPixels.length; p++) {
            if (maskPixels[p] != 0) {
                sum += fPixels[p];
                sumSqr += sqr(fPixels[p]);
                count++;
            }
        }
        double avg = sum*(1./count);
        double stddev = Math.sqrt(sumSqr*(1./count) - avg*avg);
        double histoMin = avg - 4*stddev;
        double histoMax = avg + stddev;
        final int histoSize = 256;          //create a histogram of intensities
        int[] histo = new int[histoSize];
        int nBelow = 0;
        for (int p=0; p<fPixels.length; p++) {
            if (maskPixels[p] != 0) {
                int bin = (int)Math.round((fPixels[p] - histoMin)*(histoSize/(histoMax-histoMin)));
                if (bin < 0)
                    nBelow++;
                else if (bin < histoSize)
                    histo[bin]++;
            }
        }
        int runningCount = nBelow;          //analyze the histogram
        double sumInt = nBelow*histoMin;
        for (int i=0; i<histoSize; i++) {
            int oldCount = runningCount;
            runningCount += histo[i];
            double intensity = histoMin + i*((histoMax-histoMin)/histoSize);
            if (runningCount > darkFraction*count) {
                double pixelsToUse = (darkFraction*count - oldCount);
                sumInt += pixelsToUse*intensity;
                break;
            }
            sumInt += histo[i]*intensity;
        }
        return sumInt/(darkFraction*count);
    }

    /** Calculates statistics of the x, y deviations of the positions from the
     *  expected (ScreenFitter) position, with respect to the positions of the
     *  nearest spots (min 2, max 8, but not too far away).
     *  Outliers in these statistics may have been tracked incorrectly.
     *  These spots are immediately deleted if very bad or if their
     *  maximum significance is below the given threshold.
     *  Otherwise, the spots with suspected bad tracking are entered into the
     *  'badness' table. (If the table is not empty, at the end of tracking,
     *  these spots will be highlighted and a table with them will be shown). */
    void makeDeviationStatistics(double significanceThresholdForDeletion) {
        int[] nValues = new int[xAxis.length];          //for statistics at each energy: avg and stddev of deviation
        double[] sumDev2 = new double[xAxis.length];    //sum of squared deviations at each energy, then normalized
        double[][] deviations2 = new double[spotPattern.size()][xAxis.length];
        double[] maxSignificance = new double[spotPattern.size()];
        boolean hasBadSpot = false;

        for (int pass=0; pass<2; pass++) {              //we may have to redo this to exclude badly tracked neighbors
            for (int spot=0; spot<spotPattern.size(); spot++) {
                if (!hasSpot[spot]) continue;
                int[] nearestSpots = spotPattern.getNearestSpotIndices(spot);
                if (pass != 0 && (badness[spot]==0 || !hasBadNeighbor(spot))) //in the 2nd pass, do only affected spots
                    continue;
                for (int i=0; i<xAxis.length; i++) {
                    if (Double.isNaN(data[INTEGRAL][spot][i]) ||
                            Double.isNaN(data[DELTAXS][spot][i]) ||
                            Double.isNaN(data[DELTAYS][spot][i]))
                        continue;
                    double dx = data[DELTAXS][spot][i];
                    double dy = data[DELTAYS][spot][i];
                    double sumDx=0, sumDy=0;
                    double firstDeltaKsqr=0;
                    int count=0;
                    for (int isp2=0; isp2 < nearestSpots.length && count < 8; isp2++) {
                        int spot2 = nearestSpots[isp2];
                        if (Double.isNaN(data[INTEGRAL][spot2][i]) ||
                                Double.isNaN(data[DELTAXS][spot2][i]) ||
                                Double.isNaN(data[DELTAYS][spot2][i]) ||
                                (pass != 0 && badness[spot2] != 0)) //don't determine deviation w.r.t. a badly tracked neighbor
                            continue;
                        double deltaKsqr = sqr(spotPattern.getKx(spot) - spotPattern.getKx(spot2)) +
                                sqr(spotPattern.getKy(spot) - spotPattern.getKy(spot2));
                        if (count==0)
                            firstDeltaKsqr = deltaKsqr;
                        else if (count > 3 && deltaKsqr > 4.1* firstDeltaKsqr)
                            break;
                        dx += data[DELTAXS][spot][i] - data[DELTAXS][spot2][i];
                        dy += data[DELTAYS][spot][i] - data[DELTAYS][spot2][i];
                        count++;
                    }
                    if (count < 2) continue;    //ignore if only one neighbor
                    double sqrDev = (dx*dx + dy*dy)/count;
                    if (Double.isNaN(sqrDev)) continue;
                    deviations2[spot][i] = sqrDev;
                    if (pass == 0) {
                        sumDev2[i] += sqrDev;
                        nValues[i]++;
                    }
                }
            }
            if (pass == 0)
                for (int i=0; i<xAxis.length; i++)           //at all energies, to mean sqr deviation
                    sumDev2[i] /= nValues[i];

            for (int spot=0; spot<spotPattern.size(); spot++) {
                if (!hasSpot[spot]) continue;
                maxSignificance[spot] = 0;
                double sumBadness2 = 0;
                int nEnergies = 0;
                for (int i=0; i<xAxis.length; i++) {
                    if (Double.isNaN(data[INTEGRAL][spot][i])) continue;
                    if (nValues[i] < 3) continue;           //not enough statistics
                    if (data[LeedSpotAnalyzer.PSIGNIFICANCE][spot][i] > maxSignificance[spot])
                        maxSignificance[spot] = data[LeedSpotAnalyzer.PSIGNIFICANCE][spot][i];
                    double maxDeviation2 = 2*sumDev2[i] + sqr(radiiInt[i]);
                    double badness2 = deviations2[spot][i]/maxDeviation2;
                    //if (spot==debugSpot)IJ.log(xAxis[i]+": deviation="+IJ.d2s(Math.sqrt(deviations2[spot][i]))+", good if ~ "+IJ.d2s(Math.sqrt(maxDeviation2)));
                    sumBadness2 += badness2;
                    nEnergies++;
                }
                double badness2 = sumBadness2/nEnergies;
                int badnessValue = (int)(100*(badness2 - 2.0)); //a more user friedly number; > 400 is awful
                if (badnessValue < 40)
                    badnessValue = 0;          //if >= 40, a questionable spot, list&highlight, user can delete
                else
                    hasBadSpot = true;
                if (spot==debugSpot)
                    IJ.log("DeviationStat pass "+(pass+1)+": "+spotPattern.getName(spot)+" badness:"+badnessValue+" (from "+nEnergies+" points) maxSignif.="+IJ.d2s(maxSignificance[spot]));
                badness[spot] = badnessValue;
            } // for spot
            if (!hasBadSpot)
                return;
        } // for pass
        for (int spot=0; spot<spotPattern.size(); spot++) { //delete the worst spots immediately
            if (badness[spot] >= 400 || (badness[spot] > 0 && maxSignificance[spot] < significanceThresholdForDeletion)) {
                hasSpot[spot] = false;           //either grossly off or always weak and questionable: delete
                if (debug || spot==debugSpot) IJ.log("delete "+spotPattern.getNameWithGroup(spot,false)+" badness="+badness[spot]+ ", maxSignif.="+IJ.d2s(maxSignificance[spot])+" (thresh="+IJ.d2s(significanceThresholdForDeletion)+")");
            }
        }

    }

    /** Calculates 2D regression for measured vs. calculated x&y for the measured spots,
     *  to determine the energy-dependent x&y offset, scale factor and angle.
     *  The energy dependence of the offset is an indication of residual electric & magnetic
     *  fields or misalignment of the e-gun */
    void do2DRegression() {
        regressionResults2D = new double[Leed2DRegression.OUTPUT_SIZE+4][xAxis.length];
        for (int j=0; j<regressionResults2D.length; j++)
            Arrays.fill(regressionResults2D[j], Double.NaN);
        double[] invEnergy = new double[xAxis.length];
        int nData = 0;
        double[] xyCalc = new double[2];
        xyCalc = screenFitter.screenCoordinates(/*kx=*/0, /*ky=*/0, /*invSqrtEnergy=*/100, xyCalc);
        double x00 = xyCalc[0];
        double y00 = xyCalc[1];
        Leed2DRegression regXY00 = new Leed2DRegression();  //for the position of the x&y spot
        Leed2DRegression regScale = new Leed2DRegression(); //for the image scale
        Leed2DRegression check = new Leed2DRegression();    //checks for degenerate set of points (of undistorted pattern)
        for (int i=0; i<xAxis.length; i++) {
            invEnergy[i] = 1./getEnergy(i);
            double invSqrtEnergy = Math.sqrt(invEnergy[i]);
            regXY00.clear();
            regScale.clear();
            check.clear();
            double maxKsqr = -1;
            double minKsqr = Double.MAX_VALUE;
            for (int spot=0; spot<spotPattern.size(); spot++) {
                if (hasSpot[spot] && data[LeedSpotAnalyzer.PSIGNIFICANCE][spot][i] > NO_SIGNIFICANCE && data[LeedSpotAnalyzer.X][spot][i] > 0) {
                    double kx = spotPattern.getKx(spot);
                    double ky = spotPattern.getKy(spot);
                    double kSquare = kx*kx + ky*ky;
                    if (maxKsqr < kSquare)
                        maxKsqr = kSquare;
                    if (minKsqr > kSquare && kSquare > 0)
                        minKsqr = kSquare;
                }
            }
            for (int spot=0; spot<spotPattern.size(); spot++) {
                double weight = significanceToWeight(data[LeedSpotAnalyzer.PSIGNIFICANCE][spot][i]);
                if (hasSpot[spot] && weight > 0 && data[LeedSpotAnalyzer.X][spot][i] > 0) {
                    double xExp = data[LeedSpotAnalyzer.X][spot][i];
                    double yExp = data[LeedSpotAnalyzer.Y][spot][i];
                    xyCalc = screenFitter.screenCoordinates(spotPattern, spot, invSqrtEnergy, xyCalc);
                    double kx = spotPattern.getKx(spot);
                    double ky = spotPattern.getKy(spot);
                    double kSquare = kx*kx + ky*ky;
                    regScale.addPoint(xExp - x00, yExp - y00, xyCalc[0] - x00, xyCalc[1] - y00, weight);
                    weight *= Math.sqrt(minKsqr/(minKsqr+kSquare)); //x,y of (0,0): lower weight for spots further out
                    regXY00.addPoint(xExp - x00, yExp - y00, xyCalc[0] - x00, xyCalc[1] - y00, weight);
                    check.addPoint(0, 0, spotPattern.getKx(spot), spotPattern.getKy(spot), 1);
                }
            }
            double[] resultOfCheck = check.getFit(Double.NaN, Double.NaN, null);  //used to check for degenerate set of points
            if (resultOfCheck[Leed2DRegression.WEIGHT]/resultOfCheck[Leed2DRegression.SUM_DATA_WEIGHTS] > 1e-3*Math.sqrt(maxKsqr)) {
                double[] regScaleOut = regScale.getFit(Double.NaN, Double.NaN, null);
                double[] regXY00Out  = regXY00.getFit(Double.NaN, Double.NaN, null);
                regressionResults2D[Leed2DRegression.X_OFFSET][i] = regXY00Out[Leed2DRegression.X_OFFSET];
                regressionResults2D[Leed2DRegression.Y_OFFSET][i] = regXY00Out[Leed2DRegression.Y_OFFSET];
                regressionResults2D[Leed2DRegression.DETERMINANT][i] = regScaleOut[Leed2DRegression.DETERMINANT];
                regressionResults2D[Leed2DRegression.WEIGHT][i]   = regScaleOut[Leed2DRegression.WEIGHT];
                regressionResults2D[Leed2DRegression.ANGLE][i]    = regXY00Out[Leed2DRegression.ANGLE];
                nData++;
            }
        }
        if (energies != null && nData >= 2) {                       //fit x & y offset vs. 1/sqrt(E) and 1/scale^2 vs. 1/E
            LeedLinearRegression xFit = new LeedLinearRegression();
            LeedLinearRegression yFit = new LeedLinearRegression();
            LeedLinearRegression angleFit = new LeedLinearRegression();
            LeedLinearRegression scaleFit = new LeedLinearRegression();

            for (int i=0; i<xAxis.length; i++) {
                double weight = regressionResults2D[Leed2DRegression.WEIGHT][i];
                if (weight > 0) {
                    double invSqrtEnergy = Math.sqrt(invEnergy[i]);
                    xFit.addPointW(invSqrtEnergy, regressionResults2D[Leed2DRegression.X_OFFSET][i], weight);
                    yFit.addPointW(invSqrtEnergy, regressionResults2D[Leed2DRegression.Y_OFFSET][i], weight);
                    scaleFit.addPointW(invEnergy[i], 1./regressionResults2D[Leed2DRegression.DETERMINANT][i], weight);
                    angleFit.addPointW(invSqrtEnergy, regressionResults2D[Leed2DRegression.ANGLE][i], weight);
                }
            }
            xOffsetInfty = xFit.getOffset();
            yOffsetInfty = yFit.getOffset();
            xOffset100eV = xFit.getSlope()*0.1;
            yOffset100eV = yFit.getSlope()*0.1;
            invScaleInfty = Math.sqrt(1./scaleFit.getOffset());
            deltaPhi   = -scaleFit.getSlope()*(1./scaleFit.getOffset());
            for (int i=0; i<xAxis.length; i++) {
                double invSqrtEnergy = Math.sqrt(invEnergy[i]);
                regressionResults2D[REG2D_DELTAX][i] = xFit.getFitValue(invSqrtEnergy);
                regressionResults2D[REG2D_DELTAY][i] = yFit.getFitValue(invSqrtEnergy);
                regressionResults2D[REG2D_SCALE_FIT][i] =           //convert to percent linear deviation
                        100*(invScaleInfty*Math.sqrt(1./(1 - deltaPhi*invEnergy[i])) - 1.0);
                regressionResults2D[REG2D_ANGLE_FIT][i] = angleFit.getFitValue(invSqrtEnergy);
            }
            double scaleOffsetPercent = scaleFit.getOffset();
            double workFunctionEstimate = (2./100.)*scaleFit.getSlope();
        }
        for (int i=0; i<xAxis.length; i++)
            regressionResults2D[Leed2DRegression.DETERMINANT][i] =  //convert to percent linear deviation
                    100*(Math.sqrt(regressionResults2D[Leed2DRegression.DETERMINANT][i]) - 1.0);
    }

    /** Returns whether a spot has a neighbor with badness != 0*/
    boolean hasBadNeighbor(int spot) {
        int[] nearestSpots = spotPattern.getNearestSpotIndices(spot);
        for (int spot2 : nearestSpots)
            if (badness[spot2] != 0)
                return true;
        return false;
    }

    /** A rough (~1%) approximation to the inverse error function erfc */
   /*unused  double approxInvErfc(double x) {
        double u = -1./Math.log(x);
        return 1./(u*(0.29+0.3*Math.sqrt(u))+0.98*Math.sqrt(u));
    }*/

    /** Returns whether any data are available for a given spot */
    public boolean hasSpot(int spotIndex) {
        return hasSpot[spotIndex];
    }

    /** Returns the data for a given type (see LeedSpotAnalyzer) and spot index */
    public double[] getData(int type, int spotIndex) {
        return data[type][spotIndex];
    }

    /** Returns the number of data types */
    public int getNDataTypes() {
        return LeedSpotAnalyzer.OUTPUT_SIZE + N_EXTRA_COLUMNS;
    }

    /** Returns the name for the given data type */
    public String getDataName(int i) {
        if (i < LeedSpotAnalyzer.OUTPUT_SIZE)
            return LeedSpotAnalyzer.DATA_NAMES[i];
        else
            return DATA_NAMES[i-LeedSpotAnalyzer.OUTPUT_SIZE];
    }

    /** Returns the background intensities averaged over the mask except spots (null if not measured) */
    public double[] getBackgroundIntensities() {
        return backgroundIntensities;
    }

    /** Returns the arrays with the overall deviations from the 2D model for each energy (stack slice).
     *  Subarray indices X_OFFSET=0, Y_OFFSET=1, DETERMINANT=2, ANGLE=3, WEIGHT=4, SUM_DATA_WEIGHTS=5
     *  are defined in Leed2DRegression; Fits REG2D_DELTAX, REG2D_DELTAY etc are defined in this class. */
    public double[][] getRegressionResults2D() {
        return regressionResults2D;
    }

    /** Returns a status string */
    public String getStatusText(boolean verbose) {
        nGoodSpots = 0;
        for (int spot=0; spot<spotPattern.size(); spot++)
            if (hasSpot!=null && hasSpot[spot])
                nGoodSpots++;
        String str = verbose ? "Spot tracking:\n" : "";
        str += nGoodSpots + " beams";
        if (nTooClose > 0) str += ", "+nTooClose+" too close (deleted)";
        if (verbose) str += "; highest raw intensity=" + LeedUtils.d2s(highestIntensity, 4) +
                " at "+LeedUtils.d2s(energyOfHighestI, 4)+", beam="+spotPattern.getName(spotOfHighestI);
        return str;
    }

    /** Returns a plot with results; index 'which' for different plots, corresponding to the PLOT_NAMES indices
     *  Returns null if no more plots available */
    public ImagePlus getPlotImp(int which) {
        LeedPlotter p = null;
        char plotType = 0;
        switch (which) {
            case 0:
                p = new LeedPlotter("I(V) Curves", data, xAxis, hasSpot, spotPattern, xAxisLabel);
                p.plotStack(new int[] {INT_I0CORRECTED},
                        new String[] {"intensity/I0"},
                        null,     //don't multiply data, we have already normalized them to 1000 for highest peak
                        new double[] {Double.NaN, Double.NaN, Double.NaN, Double.NaN});
                plotType = 'I';
                break;

            case 1:
                p = new LeedPlotter("Spot X/Y Deviations (pxl)", data, xAxis, hasSpot, spotPattern, xAxisLabel);
                p.plotStack(new int[] {DELTAX, DELTAXS, DELTAY, DELTAYS, LeedSpotAnalyzer.PSIGNIFICANCE},
                        new String[] {"x(spot) - fit", "x(smooth) - fit", "y(spot) - fit",  "y(smooth) - fit", "significance"},
                        null,     //don't multiply data
                        new double[] {Double.NaN, Double.NaN, -40, 40});
                plotType = 'D';
                break;

            case 2:
                p = new LeedPlotter("Spot Radii", data, xAxis, hasSpot, spotPattern, xAxisLabel);
                p.plotRadii(radiiInt, radiiSup, xAxis==energies, spotBackgrShape);
                plotType = 'r';
                break;
            case 3:
                if (xAxis != energies && LeedParams.get(LeedParams.LEEMMODE) == 0) return null;
                p = new LeedPlotter("I(V) Quality Statistics", data, xAxis, hasSpot, spotPattern, xAxisLabel);
                double v0i = LeedParams.get(LeedParams.V0I);
                double eVsmooth = LeedParams.get(LeedParams.SMOOTHEV);
                double eVseparate = LeedParams.get(LeedParams.SPLITEV);
                p.plotRFactorOfEquivalent(v0i, eVsmooth, eVseparate);
                plotType = 'R';
                break;
            case 4:
                p = new LeedPlotter("Overall X/Y Deviations",  data, xAxis, hasSpot, spotPattern, xAxisLabel);
                boolean showFit = !Double.isNaN(xOffsetInfty);
                double[] xy0 = screenFitter.screenCoordinates(0, 0, /*invSqrtEnergy=arbitrary*/0.1, new double[2]);
                ArrayList<double[]> yDataL = new ArrayList<double[]>();
                String yLabels = "";
                yDataL.add(regressionResults2D[Leed2DRegression.X_OFFSET]);
                yLabels += "x (0,0) - " + IJ.d2s(xy0[0],1) + " (pxl)\n";
                if (showFit) {
                    yDataL.add(regressionResults2D[REG2D_DELTAX]);
                    yLabels += "x (0,0) fit\n";
                }
                yDataL.add(regressionResults2D[Leed2DRegression.Y_OFFSET]);
                yLabels += "y (0,0) - " + IJ.d2s(xy0[1],1) + " (pxl)\n";
                if (showFit) {
                    yDataL.add(regressionResults2D[REG2D_DELTAY]);
                    yLabels += "y (0,0) fit\n";
                }
                yDataL.add(regressionResults2D[Leed2DRegression.DETERMINANT]);
                yLabels += "scale deviation (%)\n";
                String vectorText = showFit ?
                        "(0,0) @ 100 eV: "+IJ.d2s(xOffset100eV,1)+", "+IJ.d2s(yOffset100eV,1)+" pxl" : "";
                if (showFit) {
                    yDataL.add(regressionResults2D[REG2D_SCALE_FIT]);
                    yLabels += "scale deviation fit (\u0394WF="+IJ.d2s(deltaPhi)+")\n";
                }
                yDataL.add(regressionResults2D[Leed2DRegression.ANGLE]);
                yLabels += "rotation angle deviation (\u00b0)\n";
                if (showFit) {
                    yDataL.add(regressionResults2D[REG2D_ANGLE_FIT]);
                    yLabels += "rotation angle dev. fit (\u00b0)\n";
                }
                double vectorYScale = 2./Math.max(maskRoiRect.width, maskRoiRect.height);
                p.plot2Dregression(yDataL.toArray(new double[yDataL.size()][]), yLabels,
                        xOffset100eV, yOffset100eV, vectorYScale, vectorText);
                plotType = 'O';
                break;
            default:
                return null;
        }

        ImagePlus imp = p.getImagePlus();
        if (imp != null) {
            String impTitle = stackImp.getShortTitle();
            if (impTitle.endsWith(LEED_Spot_Tracker.STACK_SUFFIX))
                impTitle = impTitle.substring(0, impTitle.length() - LEED_Spot_Tracker.STACK_SUFFIX.length());
            imp.setTitle(imp.getTitle() + " of "+impTitle+LeedUtils.getDateFormatted(" [HH'h'mm]"));
            addToPlotImpList(imp, plotType);
            return imp;
        }
        else return null;
    }

    /** Creates and shows a quick output plot with the "best" I(V) curves of equivalent spots for aligning,
     *  Returns false is no suitable beam group was found */
    boolean makeQuickPlot(boolean errorIfNoEquivalent) {
        int nSpots = data[0].length;
        LeedIntegerArray nBeamsPerGroup = new LeedIntegerArray();
        for (int spot=0; spot<nSpots; spot++)
            if (hasSpot[spot]) {
                int group = spotPattern.getGroup(spot);
                if (group >= 0)                    //ignore symmetry-forbidden groups
                    nBeamsPerGroup.increment(group);
            }

        int maxPerGroup = 0;
        for (int g=0; g<nBeamsPerGroup.size(); g++)
            if (nBeamsPerGroup.get(g) > maxPerGroup)
                maxPerGroup = nBeamsPerGroup.get(g);
        if (maxPerGroup == 0) {
            if (errorIfNoEquivalent)
                IJ.error(LEED_Spot_Tracker.PLUGIN_NAME,
                        "No equivalent overlapping beams.\n+Producing full plots instead");
            return false;
        }
        // For the plot of selected symmetry-equivalent beams, get the spot group with the highest figure of merit:
        // overlap * sum of intensities in overlap range
        double bestMerit = 0;
        int[] bestBeams = null;
        int[] beams = new int[maxPerGroup];
        int bestGroup = -1;
        int minPerGroup = maxPerGroup;
        if (maxPerGroup > 4) minPerGroup = 3;    //for 6-fold symmetry, also 3 equivalent are enough for alignment
        for (int g=0; g<nBeamsPerGroup.size(); g++)
            if (nBeamsPerGroup.get(g) >= minPerGroup) {
                int nBeams = 0;
                Arrays.fill(beams, -1);
                for (int i=0; i<spotPattern.size(); i++) {
                    if (!hasSpot[i]) continue;
                    if (spotPattern.getGroup(i) == g)
                        beams[nBeams++] = i;
                }
                int nOverlap = 0;
                double sumI = 0;
                for (int e=0; e<xAxis.length; e++) {
                    double sum = 0;
                    for (int i=0; i<nBeams; i++)
                        sum += data[INT_I0CORRECTED][beams[i]][e];
                    if (!Double.isNaN(sum)) {
                        sumI += sum;
                        nOverlap++;
                    }
                }
                double merit = sumI*nOverlap;    //=avg intensity * overlap^2; more weight for overlap than for intensity
                if (merit > bestMerit) {
                    bestBeams = (int[])beams.clone();
                    bestGroup = g;
                    bestMerit = merit;
                }
            }
        if (bestGroup < 0) {
            if (errorIfNoEquivalent)
                IJ.error(LEED_Spot_Tracker.PLUGIN_NAME,
                        "No equivalent overlapping beams.\n+Producing full plots instead");
            return false;
        }

        int n=0;
        for (int i =0; i<bestBeams.length; i++)
            if (bestBeams[i] >= 0)
                n++;

        Plot plot = new Plot("Selected I(V) Curves (group "+bestGroup+")", xAxisLabel, "intensity");
        Color[] colors = LeedPlotter.getColors(n);
        for (int i =0; i<bestBeams.length; i++) {
            int index = bestBeams[i];
            if (index >= 0) {
                plot.setColor(colors[i]);
                double[] yValues = data[INT_I0CORRECTED][index];
                plot.addPoints(xAxis, yValues, Plot.LINE);
                plot.setLabel(i, spotPattern.getName(index));
            }
        }
        plot.setColor(Color.BLACK);
        plot.setLegend(null, Plot.AUTO_POSITION);
        plot.show();
        ImagePlus imp = plot.getImagePlus();
        addToPlotImpList(imp, 'Q');

        return true;
    }

    /** Adds an ImagePlus with a plot to the plotImpList, so they can be all closed.
     *  Plot types are 'R' (R factor statistics), 'Q' (quick I(V) of equiv. beams),
     *  'I' (stack of all I(V) curves), 'D' (position deviations), 'r' (radii),
     *  'i' (I0), and 'm' (custom from More>>Plot...).
     *  Do not use 'N'; it is reserved for 'all but new' in closeAllPlots. */
    public static void addToPlotImpList(ImagePlus imp, char plotType) {
        plotImpList.put(imp, plotType + Integer.toString(currentRunNumber));  //value is e.g. "Q12" for quick plot of run 12
    }

    /** Closes all plots except the types in the String, e.g. "QR" for quick selected curves and for R factor statistics.
     *  See addToPlotImpList for the other types.
     *  If exceptThese contains 'N', closes only old plots, not those of the last run */
    public static void closeAllPlots(String exceptThese) {
        if (exceptThese == null) exceptThese = "";
        boolean concurrentModificationException;
        do {
            concurrentModificationException = false;
            try {
                for (ImagePlus imp : plotImpList.keySet()) {    //may throw a ConcurrentModificationException
                    if (imp == null) continue;
                    if (!imp.isVisible())
                        plotImpList.remove(imp); //image had been closed already
                    else {
                        String typeAndNumber = plotImpList.get(imp);        //e.g. "Q12" for quick plot of run 12
                        char plotType = typeAndNumber.charAt(0);
                        if (exceptThese.indexOf(plotType) >= 0) continue;  //plotType is an exception
                        int runNumber = (int)Tools.parseDouble(typeAndNumber.substring(1));
                        if (exceptThese.indexOf('N') >= 0 && runNumber ==  currentRunNumber) continue; //current run is an exception
                        imp.close();
                    }
                }
            } catch (ConcurrentModificationException e) { concurrentModificationException = true; }
        } while (concurrentModificationException);
    }

    /** Returns whether there are any open plots, of any type if char=0, or the given type:
     *  'R' Rfactor statistics, 'Q' quick IV of equivalnt curves, 'I' I(V) curves, 'D' position deviations, 'r' radii. */
    public static boolean anyOpenPlots(char whichType) {
        boolean concurrentModificationException;
        do {
            concurrentModificationException = false;
            try {
                for (ImagePlus imp : plotImpList.keySet()) {    //may throw a ConcurrentModificationException
                    if (imp == null) continue;
                    if (!imp.isVisible())
                        plotImpList.remove(imp); //image had been closed already
                    else {
                        if (whichType == 0) return true;
                        String typeAndNumber = plotImpList.get(imp);        //e.g. "Q12" for quick plot of run 12
                        char plotType = typeAndNumber.charAt(0);
                        if (whichType == plotType) return true;
                    }
                }
            } catch (ConcurrentModificationException e) { concurrentModificationException = true; }
        } while (concurrentModificationException);
        return false;
    }

    /** Shows the spot positions as overlay. Can be called in a background thread, interruptible */
    public void showOverlay(ImagePlus imp) {
        IJ.showStatus("creating overlay");
        imp.setOverlay(null);               //if we have many spots marked, deleting by name takes too long
        spotTracker.drawBasicOverlay();     //deleting all and adding energies & mask is faster
        LeedOverlay.add(imp, xAxis, energies!=null, maskRoi);
        double minDistance1eV = Double.NaN;
        if (screenFitter != null) {
            double minDistance = spotPattern.getMinDistance() * screenFitter.kToScreenScale();
            minDistance1eV = minDistance * Math.sqrt(getEnergy(sliceOfIndexInput));
        }
        for (int i=0; i<imp.getNSlices(); i++) {
            if (Thread.currentThread().isInterrupted()) return;
            if (!Double.isNaN(minDistance1eV)) {
                float nameSize = (int)Math.round(0.3*minDistance1eV/Math.sqrt(getEnergy(i)));
                LeedOverlay.setFontSize(nameSize);
            }
            double radiusInt = radiiInt[i];
            double radiusSup = radiiSup[i];
            for (int spot=0; spot<spotPattern.size(); spot++)
                if (hasSpot[spot] && !Double.isNaN(data[INTEGRAL][spot][i]))
                    LeedOverlay.add(imp, i+1, data[X][spot][i],
                            data[Y][spot][i],
                            spotPattern.isSuperstructure(spot) ? radiusSup : radiusInt,
                            spotPattern.getNameForOvly(spot));
        }

        HashSet<String> spotNamesForOvly = new HashSet<String>();
        for (int spot=0; spot<spotPattern.size(); spot++)
            if (hasSpot(spot) && badness[spot] > 0)
                spotNamesForOvly.add(spotPattern.getNameForOvly(spot));
        if (spotNamesForOvly.size() > 0)
            LeedOverlay.highlightCircles(imp, spotNamesForOvly, /*strength=*/2);
        imp.draw();
        IJ.showStatus("");
    }

    /** Gets a table with the spots, their groups, no of data points, slice numbers where the spots appear and disappear,
     *  max significance and (if there are spots with badness>0) their badness.
     *  With allSpots=true, lists all spots with valid intensity data.
     *  With allSpots=false, returns null if there are no bad spots, and only lists the spots with badness > 0.
     */
    ResultsTable makeSpotTable(boolean allSpots) {
        boolean hasBadSpot = false;
        for (int spot=0; spot<spotPattern.size(); spot++)
            if (hasSpot(spot) && badness[spot] > 0) {
                hasBadSpot = true;
                break;
            }
        if (!allSpots && !hasBadSpot) return null;

        ResultsTable rt = new ResultsTable();
        int n = xAxis.length;
        int nBeams=0, nBad=0, nGroups=0, nData=0;
        int lastGroup = Integer.MAX_VALUE;
        int minFirst=Integer.MAX_VALUE, maxLast=-1;
        double maxmaxSig = 0;
        int iOfMaxMaxSig = -1;
        for (int spot=0; spot<spotPattern.size(); spot++) {
            if (hasSpot(spot) && (allSpots || badness[spot] > 0)) {
                int first = -1, last = -1, nGood = 0;
                double maxSignificance = 0;
                int iOfMaxSignificance = -1;
                for (int i=0; i<data[INTEGRAL][spot].length; i++) {
                    if (!Double.isNaN(data[INTEGRAL][spot][i])) {
                        if (first < 0) first = i;
                        last = i;
                        if (data[LeedSpotAnalyzer.PSIGNIFICANCE][spot][i] > maxSignificance) {
                            maxSignificance = data[LeedSpotAnalyzer.PSIGNIFICANCE][spot][i];
                            iOfMaxSignificance = i;
                        }
                        nGood++;
                    }
                }
                rt.incrementCounter();
                nBeams++;
                rt.addLabel(spotPattern.getName(spot).replace(',', '|'));
                int group = spotPattern.getGroup(spot);
                rt.addValue("Group", group);
                if (group != lastGroup) {
                    nGroups++;
                    lastGroup = group;
                }
                if (hasBadSpot) {
                    rt.addValue("Badness", badness[spot]);
                    if (badness[spot] != 0) nBad++;
                }
                rt.addValue("DataPoints", nGood);
                nData += nGood;
                rt.addValue("First", xAxis[first]);
                if (first < minFirst)
                    minFirst = first;
                rt.addValue("Last", xAxis[last]);
                if (last > maxLast)
                    maxLast = last;
                rt.addValue("MaxSignificance", maxSignificance);
                if (iOfMaxSignificance >= 0)
                    rt.addValue("at", xAxis[iOfMaxSignificance]);
                else
                    rt.addValue("at", "");
                if (maxSignificance > maxmaxSig) {
                    maxmaxSig = maxSignificance;
                    iOfMaxMaxSig = iOfMaxSignificance;
                }
                    
            }
        }
        if (allSpots && nBeams > 0) {
            rt.incrementCounter();
            rt.addLabel("Total "+nBeams);
            rt.addValue("Group", nGroups);
            if (hasBadSpot)
                rt.addValue("Badness", nBad);
            rt.addValue("DataPoints", nData);
            rt.addValue("First", minFirst < xAxis.length ? xAxis[minFirst] : Double.NaN);
            rt.addValue("Last", maxLast>= 0  ? xAxis[maxLast] : Double.NaN);
            rt.addValue("MaxSignificance", maxmaxSig);
            if (iOfMaxMaxSig >= 0)
                rt.addValue("at", xAxis[iOfMaxMaxSig]);
            else
                rt.addValue("at", "");
        }
        rt.setDecimalPlaces(rt.getColumnIndex("MaxSignificance"), 1);
        double xStep = LeedUtils.getFirstAndIncrement(xAxis)[1];
        if (!Double.isNaN(xStep)) {
            int nDigits = LeedUtils.isInteger(xStep) ? 0 :
                    Math.max(0, (int)Math.ceil(-Math.log10(xStep)+0.9));
            rt.setDecimalPlaces(rt.getColumnIndex("First"), nDigits);
            rt.setDecimalPlaces(rt.getColumnIndex("Last"), nDigits);
            rt.setDecimalPlaces(rt.getColumnIndex("at"), nDigits);
        }

        return rt;
    }

    void addProgress(double progressInc) {
        progress += progressInc;
        if (progress > nextProgressToShow) {
            IJ.showProgress(progress);
            nextProgressToShow = progress + 0.01;
        }
    }

    /** If we have a LeedDarkFlatVirtualStack, checks whether it has a bad flat field,
     *  i.e. a non-positive flat field value. */
    boolean hasBadFlat(boolean showError) {
        ImageStack stack = stackImp.getStack();
        if (stack instanceof LeedDarkFlatVirtualStack) {
            int badSlice = ((LeedDarkFlatVirtualStack)stack).getSliceWithBadFlat();
            if (badSlice >= 0) {
                if (showError)
                    IJ.error(LEED_Spot_Tracker.PLUGIN_NAME, "Error: Non-positive (processed) flat field value in slice "+badSlice);
                return true;
            }
        }
        return false;
    }

    /** Checks whether processed I0 has too many NaNs or any non-positive values
     *  and if so returns an error String on it.
     *  Returns null if no error. */
    String processedI0error() {
        if (processedI0 == null) return null;
        int nNaN = 0;
        for (int i=0; i<processedI0.length; i++) {
            if (Double.isNaN(processedI0[i])) nNaN++;
            if (processedI0[i] <=0) return "Error: (processed) I0 not positive: "+(float)processedI0[i]+" @"+IJ.d2s(xAxis[i],1);
        }
        if (nNaN >= processedI0.length/2)
            return "Error: (processed) I0 contains "+nNaN+"/"+processedI0.length+" invalid numbers";
        return null;
    }

    /** Returns the energy of slice i+1, or LeedRadiusSelector.UNKNOWN_ENERGY if we have no
     *  energies array (this is possible if the energy is constant, but unknown).
     *  Note that the energy might be different from the x axis variable,
     *  e.g. for time-resolved measurements */
    double getEnergy(int i) {
        if (energies != null)
            return energies[i];
        else
            return LeedRadiusSelector.UNKNOWN_ENERGY;
    }

    static double sqr(double x) {
        return x*x;
    }
}
