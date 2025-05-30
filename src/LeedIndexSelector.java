import ij.*;
import ij.IJ;
import ij.gui.*;
import ij.process.*;
import ij.measure.Measurements;
import ij.plugin.filter.MaximumFinder;
import ij.util.Tools;
import ij.plugin.frame.Recorder;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

/**
 *  This class sets the relation between the reciprocal lattice and the screen coordinates
 *  (i.e., the LeedScreenFitter), based on the maxima in the image.
 *  Detection of maxima and the first dialog to select the energy and maximum detection threshold
 *  are in the LeedIndexSliceSelector class; the dialog to name the maxima and for
 *  identification of the maxima by iterative improvement of the LeedScreenFitter
 *  is in this class.
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

public class LeedIndexSelector  implements Runnable, DialogListener, ActionListener, RoiListener {
    static final double RESIDUALS_CONVERGENCE_LIMIT = 0.8;        // if rms residuals are <= this number (in pixels), don't refine fit
    static int firstSpotSelected = -1;  //index in spotPattern
    LEED_Spot_Tracker spotTracker;
    ImagePlus stackImp;
    ImagePlus maskImp;
    Roi maskRoi;
    double[][] energiesEtc;
    int xAxisVariable;                  //index of x axis in energiesEtc
    LeedSpotPattern spotPattern;
    double[] screenFitterArray;         //previous screen fitter values
    LeedScreenFitter screenFitter;      //the result of using this class
    double maskCenterX, maskCenterY, maskDiameter;
    double radius;                      // integration and search radius
    double[] xMax, yMax;                // pixel coordinates of maxima (spots) found
    int nMaxima;                        // number of maxima in xmax, yMax
    int clickedSpotNumber = -1;         // points to user-selected spot in xMax, yMax
    int[] enteredIndices;               // indices in the spotPattern of spots entered manually
    int[] allIndices;                   // indices in the spotPattern as entered or guessed
    int nEntered = 0;                   // number of spots entered manually
    int stackSlice;
    double energy;
    double minScreenDistance;          //min distance of spots on the screen
    double searchRadius=Double.NaN;    // in pixels, 0 or NaN for auto
    double defaultSearchRadius=Double.NaN;// 0 or NaN for auto

    // GUI-related
    GenericDialog gd;
    MultiLineLabel messageLabel;
    MultiLineLabel footnoteLabel;
    Label statusLabel;
    TextField[] stringFields;
    Button  okButton;
    Button  setIndexButton;
    Label   setIndexMessage;
    boolean debug;
    String footnote;


    public LeedIndexSelector(LEED_Spot_Tracker spotTracker, ImagePlus stackImp, ImagePlus maskImp, Roi maskRoi,
            double[][] energiesEtc, int xAxisVariable, LeedSpotPattern spotPattern, double[] screenFitterArray) {
        this.spotTracker = spotTracker;
        this.stackImp = stackImp;
        this.maskImp = maskImp;
        this.maskRoi = maskRoi;
        this.energiesEtc = energiesEtc;
        this.xAxisVariable = xAxisVariable;
        this.spotPattern = spotPattern;
        this.screenFitterArray = screenFitterArray;
    }

    /** With given position of the maxima (spots), tries the (previous) relation
     *  between indices and screen coordinates as given by screenFitterArray
     *  passed with the constructor.
     *  Returns a refined such relation if successful, otherwise null.
     *  This method is used when processing several LEED I(V) movies with
     *  essentially the same geometry. */
    public LeedScreenFitter useScreenFitterArray(double[] xMax, double[] yMax, int stackSlice, double energy,
            double searchRadius, boolean showLabels) {
        this.xMax = xMax;
        this.yMax = yMax;
        this.nMaxima = xMax.length;
        this.stackSlice = stackSlice;
        this.energy = energy;
        this.searchRadius = searchRadius;
        radius = LeedRadiusSelector.radius(energy, false);
        enteredIndices = new int[nMaxima];
        allIndices = new int[nMaxima];
        Arrays.fill(enteredIndices, -1);
        nEntered = 0;
        LeedScreenFitter initialScreenFitter = new LeedScreenFitter(screenFitterArray);
        this.maskCenterX = initialScreenFitter.getX00();    //we use the (0,0) coordinates for subsequent ScreenFitters
        this.maskCenterY = initialScreenFitter.getY00();

        double initialSearchRadius = searchRadius >= 1 ? searchRadius : radius;
        for (int iSpot=0; iSpot<spotPattern.size(); iSpot++) {
            double[] xy = new double[2];
            xy = initialScreenFitter.screenCoordinates(spotPattern, iSpot, 1./Math.sqrt(energy), xy);
            for (int iMax=0; iMax<nMaxima; iMax++) {
                double distSqr = sqr(xMax[iMax] - xy[0]) + sqr(yMax[iMax] - xy[1]);
                if (distSqr < sqr(initialSearchRadius)) {
                    enteredIndices[iMax] = iSpot;
                    nEntered++;
                }
            }
        }
        //IJ.log("nFound="+nFound+" min="+Math.max(nMaxima/4, 3));
        if (nEntered < Math.max(nMaxima/2, 3))            //very few spots found; looks suspicious
            return null;

        findSpotsAndFit();

        double rmsResiduals = screenFitter == null ? Double.NaN :
                screenFitter.getRmsResiduals(spotPattern, allIndices, xMax, yMax);
        if (Double.isNaN(rmsResiduals))
            screenFitter = null;

        if (showLabels && screenFitter != null) {
            showLabels(false);
            IJ.showStatus(screenFitter.getStatusText()+"/"+nMaxima+" spots");
        }
        return screenFitter;
    }


    /** Asks the user to set the basis and calls LEED_Spot_Tracker.setScreenFitter if successful.
     *  Called in a separate thread. */
    public void run() {
        try {
            if (!spotTracker.waitForStackUpdate()) {
                spotTracker.enableAndHighlightComponents(true);
                return;
            }
            stackImp.setOverlay(null);                  //if we have many spots marked, deleting by name takes tool long
            LeedOverlay.add(stackImp, energiesEtc[xAxisVariable], xAxisVariable==LEED_Spot_Tracker.ENERGY, maskRoi); //rather start from scratch
            stackImp.draw();

            Rectangle maskRect = maskRoi.getBounds();
            maskCenterX = maskRect.x + 0.5*maskRect.width;
            maskCenterY = maskRect.y + 0.5*maskRect.height;
            maskDiameter = Math.sqrt(sqr(maskRect.width) + sqr(maskRect.height));

            //show maxima and ask the user to set energy and significance threshold
            Object maximaOrScreenFitter = (new LeedIndexSliceSelector()).getMaxima(spotTracker,
                    this, stackImp, maskImp, maskRoi, energiesEtc[LEED_Spot_Tracker.ENERGY], screenFitterArray);
            if (maximaOrScreenFitter == null) {         //LeedIndexSliceSelector canceled
                spotTracker.enableAndHighlightComponents(true);
                return;
            }
            if (!stackImp.lock()) {                     //lock to avoid changing slice
                IJ.error("Image locked", "The stack must not be in use by other operations.\n"+
                        "If you erroneously get this message,\n"+
                        "use Plugins>Utilities>Reset to unlock.");
                spotTracker.enableAndHighlightComponents(true);
                return;
            }

            stackSlice = stackImp.getCurrentSlice();
            energy = LEED_Spot_Tracker.sliceToEnergy(energiesEtc[LEED_Spot_Tracker.ENERGY], stackSlice);
            radius = LeedRadiusSelector.radius(energy, false);

            if (maximaOrScreenFitter instanceof double[][]) {   //we have to ask the user for labels for the maxima
                double[][] xymax = (double[][])maximaOrScreenFitter;
                xMax = xymax[0];                            //all maxima found in the image
                yMax = xymax[1];
                nMaxima = xMax.length;


                enteredIndices = new int[nMaxima];
                allIndices = new int[nMaxima];
                Arrays.fill(enteredIndices, -1);
                nEntered = 0;

                stackImp.killRoi();
                LeedOverlay.removeSpotOverlays(stackImp);
                for (int i=0; i<nMaxima; i++)
                    LeedOverlay.add(stackImp, stackSlice, xMax[i], yMax[i], radius, null);
                stackImp.updateAndDraw();
                IJ.setTool("point");
                int pointRoiSize = PointRoi.getDefaultSize();
                PointRoi.setDefaultSize(4); //4=XL, 5=XXL
                boolean isSuperstructure = spotPattern.isSuperstructure();
                String[] selectedHK = new String[] {"1", "0"};  //create suggestion for selected spot
                if (firstSpotSelected >=0 && firstSpotSelected < spotPattern.size()) {
                    String[] previousFirstHK = Tools.split(spotPattern.getName(firstSpotSelected), ",;|");
                    if (previousFirstHK.length == 2)
                        selectedHK = previousFirstHK;           //if we have a precious entry, use this
                }
                gd = new NonBlockingGenericDialog("Select And Label Spots");
                gd.setInsets(5, 0, 0); //top, left, bottom
                gd.addMessage("Click on a spot with known index,\nenter the index and press 'Set'");
                messageLabel = (MultiLineLabel)(gd.getMessage());
                gd.addStringField("Spot index (", selectedHK[0]);
                gd.addToSameRow();
                gd.addStringField(", ",selectedHK[1]);
                gd.addToSameRow();
                gd.addMessage(isSuperstructure ? ") *" : ")");
                Panel panel = new Panel();
                setIndexButton = new Button("Set");
                panel.add(setIndexButton);
                setIndexMessage = new Label("                                       ");
                panel.add(setIndexMessage);
                gd.addPanel(panel);
                gd.addStringField("Search tolerance", "");
                gd.addToSameRow();
                gd.addMessage("pixels **");
                footnote = "** leave empty for auto";
                if (isSuperstructure)
                    footnote = "*  Integer or decimal numbers, or fractions, e.g. -5/2\n" + footnote;
                gd.setInsets(5, 0, 0);      //top, left, bottom
                gd.addMessage(footnote + "\n");    //must always be a MultiLineLabel
                footnoteLabel = (MultiLineLabel)(gd.getMessage());
                gd.setInsets(isSuperstructure ? 5 : 0, 0, 5); //top, left, bottom
                gd.addMessage("");  //status field
                gd.addCheckbox("Log debug info", false);
                gd.addHelp(LeedSpotTrackerHelp.getSetIndices2Help());
                statusLabel = (Label)(gd.getMessage());
                Vector stringFieldV = gd.getStringFields();
                stringFields = (TextField[])(stringFieldV.toArray(new TextField[0]));
                okButton = gd.getButtons()[0];
                okButton.setEnabled(false);
                setIndexButton.setEnabled(false);
                setIndexButton.addActionListener(this);
                gd.addDialogListener(this);
                Roi.addRoiListener(this);
                spotTracker.setCurrentFrontDialog(gd);
                EventQueue.invokeLater(new Runnable() {public void run() {
                        Recorder.suspendRecording();    //don't record selecting the point
                    }});

                gd.showDialog();
                if (gd.wasCanceled())
                    screenFitter = null;

                EventQueue.invokeLater(new Runnable() {public void run() {
                        Recorder.resumeRecording();    //resume recording after selecting the point
                    }});

                PointRoi.setDefaultSize(pointRoiSize);
                Roi.removeRoiListener(this);
            } else      // no need to label spots, we got a LeedScreenFitter
                screenFitter = (LeedScreenFitter)maximaOrScreenFitter;

            stackImp.unlock();
            spotTracker.setCurrentFrontDialog(null);
            spotTracker.toFront();
            IJ.setTool(Toolbar.RECTANGLE);
            if (screenFitter !=  null) {
                spotTracker.setScreenFitter(screenFitter, stackSlice-1);
            } else {
                LeedOverlay.removeSpotOverlays(stackImp);
                LeedOverlay.removeNameOverlays(stackImp);
                stackImp.updateAndDraw();
            }
            spotTracker.enableAndHighlightComponents(true);
        } catch (Exception e) {
            IJ.handleException(e);
            stackImp.unlock();
            Roi.removeRoiListener(this);
            spotTracker.enableAndHighlightComponents(true);
        }
    }

    /** This callback method is called when the user changes choices or text fields the dialog. */
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        try{
            if (e == null) return true;             //ignore call when closing the dialog
            debug = gd.getNextBoolean();
            if (e.getSource()==stringFields[2]) {   //search radius
                double newSearchRadius = Tools.parseDouble(stringFields[2].getText());
                if ((!(newSearchRadius == searchRadius) && newSearchRadius > 0) ||
                        (Double.isNaN(newSearchRadius) && !Double.isNaN(searchRadius))) {
                    searchRadius = newSearchRadius;
                    if (nEntered > 0)
                        findSpotsAndShow();
                }
                return (okButton.isEnabled());
            } else
                enableOrDisableSetButton();
            IJ.setTool("point");        //just to make sure
            return false;
        } catch(Exception ex) {IJ.handleException(ex);return false;}
    }

    /** setIndexButton pressed: Adds spot to list of known spots */
    public void actionPerformed(ActionEvent e) {
        if (e.getSource() == setIndexButton) try {
            String hString = stringFields[0].getText();
            String kString = stringFields[1].getText();
            double h = LeedUtils.toNumber(hString);
            double k = LeedUtils.toNumber(kString);
            int index = spotPattern.getSpotIndex(h, k);
            if (index < 0) return; //should never happen
            enteredIndices[clickedSpotNumber] = index;
            //IJ.log(nEntered+" Enter index "+index+":"+spotPattern.getName(index)+" for spot#"+clickedSpotNumber+"/"+nMaxima+" at "+(int)xMax[clickedSpotNumber]+","+(int)yMax[clickedSpotNumber]);
            if (nEntered == 0)
                firstSpotSelected = index;  //remember first spot entered for the next invocation
            nEntered++;
            clickedSpotNumber = -1;
            setIndexButton.setEnabled(false);
            stackImp.killRoi();
            findSpotsAndShow();
        } catch(Exception ex) {IJ.handleException(ex);}
    }

    /** Handles clicking on the stack image with the point tool, and disables usage of other tools */
    public void roiModified(ImagePlus imp, int id) {
        if (imp != stackImp || id == RoiListener.DELETED) return;
        if (imp.getRoi() == null) return;
        if (id != RoiListener.CREATED || !(imp.getRoi() instanceof PointRoi) || imp.getRoi().getPolygon().npoints != 1) {
            imp.killRoi();
            IJ.setTool("point");
            IJ.beep();
            return;
        }
        Roi roi = imp.getRoi();
        Rectangle r = roi.getBounds();
        clickedSpotNumber = -1;
        for (int i=0; i<nMaxima; i++)
            if (sqr(r.x - xMax[i]) + sqr(r.y - yMax[i]) < sqr(radius)) {
                clickedSpotNumber = i;
                roi.setLocation((int)Math.round(xMax[i]), (int)Math.round(yMax[i]));
                enableOrDisableSetButton();
                gd.toFront();
                break;
            }
            if (clickedSpotNumber == -1)
                stackImp.killRoi();
    }

    /** Does the fit and shows the result. Enables the OK button when successful. */
    void findSpotsAndShow() {
        findSpotsAndFit();
        showLabels(true);
        double rmsResiduals = screenFitter == null ? Double.NaN :
                screenFitter.getRmsResiduals(spotPattern, allIndices, xMax, yMax);
        if (Double.isNaN(rmsResiduals))
            screenFitter = null;
        okButton.setEnabled(screenFitter != null);
        if (screenFitter != null) {
            messageLabel.setText("If all spots have correct labels, press OK.\n"+
                    "Otherwise click on a new spot and enter its index.");
            statusLabel.setText(screenFitter == null ? "" : screenFitter.getStatusText()+"/"+nMaxima+" spots");
            if (!Double.isNaN(defaultSearchRadius))
				footnoteLabel.setText(footnote + " ("+IJ.d2s(defaultSearchRadius, 1)+")");
        } else {
            messageLabel.setText("Assignment of spots still unclear.\n"+
                    "Please click on a new spot and enter its index.");
        }
        stackImp.updateAndDraw();
    }

    /** Checks whether the 'Set' button can be enabled, based on indices in pattern and spot clicked.
     *  Also provides feedback on indices in setIndexMessage and disabled the OK button on unfinished input */
    void enableOrDisableSetButton() {
        boolean ok = true;
        if (clickedSpotNumber < 0) {
            ok = false;
            setIndexMessage.setText("No (new) spot clicked");
        } else if (enteredIndices[clickedSpotNumber]>=0) {
            int spot = enteredIndices[clickedSpotNumber];
            setIndexMessage.setText("RENAME THE " + spotPattern.getName(spot) + " SPOT? ");
        } else {
            String hString = stringFields[0].getText();
            String kString = stringFields[1].getText();
            double h = LeedUtils.toNumber(hString);
            double k = LeedUtils.toNumber(kString);
            if (Double.isNaN(h) || Double.isNaN(k)) {
                ok = false;
                setIndexMessage.setText("Invalid number");
            } else {
                int index = spotPattern.getSpotIndex(h, k);
                if (index < 0) {
                    ok = false;
                    setIndexMessage.setText("Not in Pattern File");
                } else if (spotPattern.getGroup(index) < 0) {
                    setIndexMessage.setText("Symmetry-forbidden spot, SURE?");
                } else if (LeedUtils.arrayIndexOf(enteredIndices, index) >=0) {
                    ok = false;
                    setIndexMessage.setText("New index needed");
                } else
                    setIndexMessage.setText("");
            }
        }
        if (screenFitter != null && !ok) {        //user input because previous basis guess not valid
            removeFitIndices();
            screenFitter = null;
            showLabels(true);
        }
        setIndexButton.setEnabled(ok);
        okButton.setEnabled(ok && screenFitter != null);
    }

    void showLabels(boolean highlightEntered) {
        int fontSize = (int)(minScreenDistance*0.25);
        if (fontSize<5) fontSize = 5;
        if (fontSize>20) fontSize = 20;
        LeedOverlay.setFontSize(fontSize);
        LeedOverlay.removeNameOverlays(stackImp);
        HashSet<String> highlightedSpots = new HashSet<String>();
        for (int i=0; i<allIndices.length; i++) {
            int index = allIndices[i];
            if (index >=0) {
                String name = spotPattern.getNameForOvly(index);
                boolean bold = highlightEntered && LeedUtils.arrayIndexOf(enteredIndices, index) >= 0;
                LeedOverlay.addName(stackImp, stackSlice, xMax[i], yMax[i], radius, name, bold);
            }
        }
        LeedOverlay.highlightCircles(stackImp, highlightedSpots, 1);
        stackImp.draw();
    }

    /** Removes the indices of the spots labelled in the fit */
    void removeFitIndices() {
        System.arraycopy(enteredIndices, 0, allIndices, 0, enteredIndices.length);
    }

	 /** Which spots to use in the current searchWhat of findSpotsAndFit */
    static final int USE_ENTERED_NEAR=0, USE_ENTERED_2=1, USE_NEAR=2, USE_NEAR_2=3, USE_ALL_ADDED=4, USE_FIT=5;

    /** Creates the LeedScreenFitter, i.e., function to translate the spot indices into screen coordinates.
     *  Uses the coordinates known from manual input (enteredIndices) and tries to add the other maxima.
     *  When successful, sets the instance variable screenFitter.
     *
     *  What to search for:
     *  (0) near neighbors of the spots added manually
     *  (1) a bit further neighbors of the spots entered manually (up to 2*smallest distance in spot pattern)
     *  (2) near neighbors of any spot known so far
     *  (3) a bit further neighbors of any spot known so far (up to 2*smallest distance in spot pattern)
     *  (4) all neighbors of any spot known so far
     *  (5) remaining spots based on polynomial fit */
    void findSpotsAndFit() {
        removeFitIndices();
        screenFitter = doTheFit(LeedScreenFitter.LINEAR);
        if (screenFitter == null) return;        //fit failed, we don't have enough information for the basis (e.g. only 0,0 spot)

        boolean[] indexFound = new boolean[spotPattern.size()];
        int nIndices = 0;
        for (int index : allIndices)
            if (index >= 0) {
                indexFound[index] = true;
                nIndices++;
            }

        // Estimate how tolerant we can be with positions, based on the minimum distance between any two spots
        // in k space and the linear component of the fit
        // We use 25% of the closest distance between any two spots, but not less than the radius of the integration disk
        double minKDistance = spotPattern.getMinDistance();
        minScreenDistance = screenFitter.kToScreenScale()*minKDistance;
        defaultSearchRadius = Math.max(0.25*minScreenDistance, radius);
        double searchRadius = (Double.isNaN(this.searchRadius) || this.searchRadius<1) ?
                defaultSearchRadius : this.searchRadius;

        double invSqrtEnergy = 1./Math.sqrt(energy);
        double[] screenXY = new double[2];
        //screenXY = new double[4]; needed for derivative consistency check in ScreenFitter
        nIndices = nEntered;
        boolean spotAdded = false;
        int functionType = LeedScreenFitter.LINEAR;
        LeedScreenFitter screenFitter = this.screenFitter;
        boolean firstForFunction = true;  //did we start a new fit level?
        int searchWhat = USE_ENTERED_NEAR;
        // starting with a linear fit, (affine transformation), firstForFunction add neighbors of the spots entered, then their neighbors,
        // finally all spots in the spot pattern description file.
        while (true) {
            if (debug)
                IJ.log(nIndices+" spots so far, try search "+searchWhat+". Previously "+this.screenFitter.getFunctionName()+
                        ", try fit level="+screenFitter.FUNCTION_NAMES[functionType]);
            spotAdded = false;
            if (searchWhat == USE_FIT) {                            //4th searchWhat pass: add spots based on fit coordinates
                for (int index=0; index<spotPattern.size(); index++) {
                    if (indexFound[index]) continue;
                    if (spotPattern.getGroup(index) < 0) continue; //ignore forbidden spot
                    screenXY = screenFitter.screenCoordinates(spotPattern, index, invSqrtEnergy, screenXY);
                    double x = screenXY[0], y=screenXY[1];
                    if (!maskRoi.contains((int)x,(int)y)) continue;
                    for (int j=0; j<allIndices.length; j++) {    //check all maxima...
                        if (allIndices[j] >= 0) continue;   // ...that have no index yet:
                        if (sqr(xMax[j]-x) + sqr(yMax[j]-y) < sqr(searchRadius)) { //position of this spot?
                            allIndices[j] = index;
                            indexFound[index] = true;
                            nIndices++;
                            spotAdded = true;
                            LeedScreenFitter tmpFitter = doTheFit(functionType);
                            if (tmpFitter != null)          //ignore if fitting fails
                                screenFitter = tmpFitter;
                            break;                          //don't search for this spot any more
                        }
                    }
                }
            } else {                                        //1st-3rd searchWhat pass: add neighbors only
                for (int iKnown = 0; iKnown < allIndices.length; iKnown++) {
                    int knownIndex = searchWhat<=USE_ENTERED_2 ? enteredIndices[iKnown] : allIndices[iKnown];
                    if (knownIndex < 0) continue;
                    double kx0 = spotPattern.getKx(knownIndex);
                    double ky0 = spotPattern.getKy(knownIndex);
                    double maxDistanceSqr = (searchWhat==USE_ENTERED_NEAR || searchWhat==USE_NEAR) ?
							sqr(1.5*minKDistance) : sqr(2.2*minKDistance);
                    if (searchWhat==USE_ALL_ADDED) maxDistanceSqr = Double.MAX_VALUE;
                    for (int neighborIndex : spotPattern.getNearestSpotIndices(knownIndex)) {
                        if (indexFound[neighborIndex]) continue;
                        if (spotPattern.getGroup(neighborIndex) < 0) continue; //ignore forbidden spot
                        if (sqr(spotPattern.getKx(neighborIndex)-kx0) + sqr(spotPattern.getKy(neighborIndex)-ky0) > maxDistanceSqr)
							break;                          //further neighbors are even more distant
                        //IJ.log("check neighbor "+spotPattern.getName(neighborIndex)+" of "+spotPattern.getName(knownIndex));
                        screenXY = screenFitter.screenCoordinates(spotPattern, knownIndex, invSqrtEnergy, screenXY);
                        double xKnownFit = screenXY[0];
                        double yKnownFit = screenXY[1];
                        screenXY = screenFitter.screenCoordinates(spotPattern, neighborIndex, invSqrtEnergy, screenXY);
                        double x = screenXY[0] - xKnownFit + xMax[iKnown];
                        double y = screenXY[1] - yKnownFit + yMax[iKnown];
                        if (debug && nIndices<20)
                            IJ.log("search for neighbor "+spotPattern.getName(neighborIndex)+" at "+(int)x+","+(int)y +
                                    " of "+spotPattern.getName(knownIndex)+" at "+(int)xMax[iKnown]+","+(int)yMax[iKnown]);
                        if (!maskRoi.contains((int)x,(int)y)) continue;
                        for (int j=0; j<allIndices.length; j++) {    //check all maxima...
                            if (allIndices[j] >= 0) continue;   // ...that have no index yet:
                            if (sqr(xMax[j]-x) + sqr(yMax[j]-y) < sqr(searchRadius)) { //position of this spot?
                                allIndices[j] = neighborIndex;
                                if (debug && nIndices<20)
                                    IJ.log("found "+spotPattern.getName(neighborIndex)+" (known="+indexFound[neighborIndex]+") at "+IJ.d2s(xMax[j],2)+","+IJ.d2s(yMax[j],2)+" near "+IJ.d2s(x,2)+","+IJ.d2s(y,2));
                                indexFound[neighborIndex] = true;
                                nIndices++;
                                spotAdded = true;
                                LeedScreenFitter tmpFitter = doTheFit(functionType);
                                if (tmpFitter != null)          //ignore if fitting fails
                                    screenFitter = tmpFitter;
                                break;                          //don't search for this neighbor spot any more
                            }
                        }
                    }
                }
            }
            if (!spotAdded || searchWhat < USE_FIT && functionType < LeedScreenFitter.HIGHEST_STD_FUNCTION &&
					LeedScreenFitter.N_PARAM[functionType+1]*4 < nIndices) {   //nothing found for the current searchWhat or extremely high redundancy
                double residuals = screenFitter.getRmsResiduals(spotPattern, allIndices, xMax, yMax);
                double residualsThreshold = Math.min(0.005*maskDiameter, 0.5*searchRadius);		//0.5% of screen diameter or half search radius
                if (debug)
                    IJ.log((!spotAdded ? "Nothing added" : "High redundancy")+" in search "+searchWhat+" with "+screenFitter.getFunctionName()+
                            "; rms residuals="+IJ.d2s(residuals)+" pxl");
                if (searchWhat < USE_FIT && functionType < LeedScreenFitter.HIGHEST_STD_FUNCTION &&
                        LeedScreenFitter.N_PARAM[functionType+1]*2 < nIndices && residuals > residualsThreshold) {
                    LeedScreenFitter tmpFitter = doTheFit(functionType+1);
                    if (tmpFitter != null) {                    //much redundancy (2-fold) and large residuals: first try higher fit order
                        screenFitter = tmpFitter;
						functionType++;
						if (searchWhat == USE_ALL_ADDED)
							searchWhat--;                       //with higher fit order, maybe we can get nearer neighbors first?
						firstForFunction = true;
						continue;
					}
                }
                if (searchWhat < USE_FIT) {
					if (!spotAdded)
						searchWhat++;                           //low residuals or can't go to higher fit order with neighbor search
                    firstForFunction = false;
                    continue;
                } else {										//final pass 'USE_FIT' based on global fit
                    if (!firstForFunction)
                        this.screenFitter = screenFitter;
                    if (functionType < LeedScreenFitter.HIGHEST_STD_FUNCTION &&
                            LeedScreenFitter.N_PARAM[functionType+1] < nIndices) {
                        functionType++;                         //if possible, try fitting higher-order polynomial to find more spots
                        firstForFunction = true;
                        LeedScreenFitter tmpFitter = doTheFit(functionType);
                        if (tmpFitter != null)
                            screenFitter = tmpFitter;
                    } else
                        break;                                  //otherwise done
                }
            } else
                firstForFunction = false;
        }
        if (this.screenFitter == null)
            this.screenFitter = screenFitter;                   //not needed? just to make sure
        
        // if possible, check whether a higher fit order improves the result
        //IJ.log("nInd="+nIndices+"; need for highest fit type: "+LeedScreenFitter.N_PARAM[LeedScreenFitter.HIGHEST_STD_FUNCTION]*1.20);
        for (functionType = this.screenFitter.getMainFunctionType() + 1;
                functionType <= LeedScreenFitter.HIGHEST_STD_FUNCTION && LeedScreenFitter.N_PARAM[functionType]*1.20 < nIndices;
                functionType++) { //factor 1.20: for trying higher fit order, we require 20% redundancy
            //IJ.log("try "+functionType);
            LeedScreenFitter fitter1 = this.screenFitter;
            double rmsResiduals1 = fitter1.getRmsResiduals(spotPattern, allIndices, xMax, yMax);
            if (rmsResiduals1 < RESIDUALS_CONVERGENCE_LIMIT)
                break;
            LeedScreenFitter fitter2 = doTheFit(functionType);
            double rmsResiduals2 = fitter2 == null ?
                    Double.MAX_VALUE : fitter2.getRmsResiduals(spotPattern, allIndices, xMax, yMax);
            int degFreedom1 = 2*nIndices - LeedScreenFitter.N_PARAM[fitter1.getFunctionType()];
            int degFreedom2 = fitter2 == null ?
                    1 :    2*nIndices - LeedScreenFitter.N_PARAM[fitter2.getFunctionType()];
            if (debug)
                IJ.log("Try higher order: "+fitter1.getFunctionName()+"->"+(fitter2==null ? "null" : fitter2.getFunctionName())+
                        ": rms/sqrt(degF) "+IJ.d2s(rmsResiduals1/Math.sqrt(degFreedom1),4)+"->"+IJ.d2s(rmsResiduals2/Math.sqrt(degFreedom2),4));
            if (sqr(rmsResiduals2)/degFreedom2 < 0.9 * sqr(rmsResiduals1)/degFreedom1)  //at least ~10% improvement of chi^2?
                this.screenFitter = fitter2;
            else
                continue;
        }
        if (debug) {    //show pixel coordinates of (0,0)
            double[] zeroSpotXY = this.screenFitter.screenCoordinates(0, 0, invSqrtEnergy, new double[2]);
            IJ.log("(0,0) spot at pixel x="+IJ.d2s(zeroSpotXY[0])+" y="+IJ.d2s(zeroSpotXY[1]));
        }
    }

    /** Tries to obtain a fit for the screen coordinates using the entered coodinates and the
     *  coordinates of the other spots whose indices have been determined so far.
     *  Returns the LeedScreenFitter if successful, or null on failure.
     *  For function types where an alternative exists, also tries the alternative, and returns the best.
     */
    LeedScreenFitter doTheFit(int functionType) {
        int nSpots = 0;
        for (int index : allIndices)
            if (index >= 0) nSpots++;
        double[] x = new double[nSpots];
        double[] y = new double[nSpots];
        double[] kx = new double[nSpots];
        double[] ky = new double[nSpots];
        for (int i=0, iSpot=0; i<allIndices.length; i++) {
            int index = allIndices[i];
            if (index < 0) continue;
            x[iSpot] = xMax[i];
            y[iSpot] = yMax[i];
            kx[iSpot] = spotPattern.getKx(index);
            ky[iSpot] = spotPattern.getKy(index);
            iSpot++;
        }
        LeedScreenFitter screenFitter = new LeedScreenFitter(functionType, energy, maskCenterX, maskCenterY);
        int alternativeType = LeedScreenFitter.ALTERNATIVE[functionType];
        LeedScreenFitter alternativeScreenFitter = alternativeType >= 0 ?
                new LeedScreenFitter(alternativeType, energy, maskCenterX, maskCenterY) : null;
        boolean resultOk = screenFitter.fitSpots(kx, ky, x, y);
        //IJ.log(nSpots+" fitOk="+resultOk+" resid="+screenFitter.getRmsResiduals(kx, ky, x, y));
        boolean alternativeResultOk = alternativeScreenFitter != null ?
                alternativeScreenFitter.fitSpots(kx, ky, x, y) : false;
        if (resultOk && alternativeResultOk) {        //if we have two possible fit functions, take the better one
            double residuals1 = screenFitter.getRmsResiduals(kx, ky, x, y);
            double residuals2 = alternativeScreenFitter.getRmsResiduals(kx, ky, x, y);
            //if(functionType==LeedScreenFitter.FOURTH_ORDER)IJ.log("residuals 4th, cub+rad^5: "+IJ.d2s(residuals1)+", "+IJ.d2s(residuals2));
            return residuals1 < residuals2 ? screenFitter : alternativeScreenFitter;
        } else if (resultOk)
            return screenFitter;
        else if (alternativeResultOk)
            return alternativeScreenFitter;
        else
            return null;
    }

    static double sqr(double x) {return x*x;}
}
