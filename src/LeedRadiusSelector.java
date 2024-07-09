import ij.*;
import ij.util.*;
import ij.IJ;
import ij.gui.*;
import ij.process.*;
import java.awt.*;
import java.util.*;
/**
 *  This class contains the dialog for selecting the integration radius for the spots,
 *  and the code for calculation of the energy-dependent radius.
 *
 *  The energy-dependent radius is calculated from
 *
 *    r^2 = r_infty^2 + r_1^2 / E
 *
 *  where r_infty is the radius at very high energies, and r_1
 *  indicates how much it should increase at low energies.
 *  (The squares of these two values are saved as parameters;
 *  the latter separately for integer and hasSuperstructure spots).
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


public class LeedRadiusSelector  implements Runnable, DialogListener {
    public static final double UNKNOWN_ENERGY = 99.99999990123456789;      //used for x-axis=time, energy, etc. if the energy is unknown. Should be a reasonable energy for saving in prefs.
    static final String[] shapeChoices = new String[] {"Concentric Circles", "Oval background", "Azimuth blur" };
    static final String[] shapesShort = new String[] {"circ.", "oval bg", "az. blur" };
    static final double MIN_RADIUS = 1.5;
    LEED_Spot_Tracker spotTracker;
    ImagePlus stackImp;
    double energyMin;
    boolean isEnergyRange;
    boolean hasSuperstructure;
    LeedLabeledField azBlurAngleField;      //numeric field with label, for en/disabling
    boolean recordingFinalized;             //to avoid macro recording twice

    public LeedRadiusSelector(LEED_Spot_Tracker spotTracker, ImagePlus stackImp,
            double energyMin, double energyMax, boolean hasSuperstructure) {
        this.spotTracker = spotTracker;
        this.stackImp = stackImp;
        this.energyMin = energyMin;
        this.isEnergyRange = energyMax > energyMin && energyMin != UNKNOWN_ENERGY;
        this.hasSuperstructure = hasSuperstructure;
    }

    /** Asks the user to set the basis and stores the values in the LeedParams. */
    public void run() {
        try{
            double radiusInt0 = radius(energyMin, false);                // integer-beam radius at lowest energy
            double radiusSup0 = radius(energyMin, true);                // superstr-beam radius at lowest energy
            int spotBackgrShape = (int)LeedParams.get(LeedParams.BACKGROUNDTYPE);
            double radiusInftySqr = LeedParams.get(LeedParams.RADIUSINFTYSQR);
            double radius1eVSqr   = LeedParams.get(LeedParams.RADIUS1EVSQR);
            double radius1eVSqrS  = LeedParams.get(LeedParams.RADIUS1EVSQRS);
            double azBlurAngleDegrees  = LeedParams.get(LeedParams.AZIMUTHBLURANGLE);
            double radiusInfty = Math.sqrt(radiusInftySqr);
            if (radiusInt0 < radiusInfty) radiusInt0 = radiusInfty;
            if (radiusSup0 < radiusInfty) radiusSup0 = radiusInfty;
            boolean isSame = radiusInt0 == radiusSup0;
            GenericDialog gd = new GenericDialog("Integration & Background Areas");
            String minEstr = isEnergyRange ? " at E="+IJ.d2s(energyMin,1) : "";
            gd.addChoice("Integration/Backgr. shape", shapeChoices, shapeChoices[spotBackgrShape]);
            String labelRInt = "Radius"+minEstr;
            if (hasSuperstructure)
                labelRInt += " (integer spots)";
            gd.addNumericField(labelRInt, radiusInt0,1);
            if (hasSuperstructure) {
                gd.addNumericField("Radius"+minEstr+" (superstructure spots)", radiusSup0,1);
                gd.addToSameRow();
                gd.addCheckbox("same", isSame);
            }
            if (isEnergyRange)
                gd.addNumericField("Radius at E -> infinity", radiusInfty ,1, 6, "*");
            gd.addNumericField("Azimuthal blur angle", azBlurAngleDegrees, 1, 6, "\u00b0");
            azBlurAngleField = LeedLabeledField.numeric(gd);
            azBlurAngleField.setEnabled(spotBackgrShape == LeedSpotAnalyzer.AZIMUTH_BLUR);
            if (isEnergyRange)
                gd.addMessage("* Must be >="+MIN_RADIUS+", but not larger than radius at E="+IJ.d2s(energyMin,1));
            
            gd.addHelp(LeedSpotTrackerHelp.getSetRadiusHelp());
            gd.addDialogListener(this);
            gd.showDialog();

            if (!gd.wasCanceled())
                dialogItemChanged(gd, null); //to make sure we write to prefs (might be also done by GenericDialog)
        } catch (Exception e) {IJ.handleException(e);}
    }

    /** This callback method is called when the user changes choices or text fields the dialog.
     *  It is called with e=null when the dialog is closed with 'OK' */
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        try {
            int spotBackgrShape = gd.getNextChoiceIndex();
            double radiusInt0 = gd.getNextNumber();
            double radiusSup0 = hasSuperstructure ? gd.getNextNumber() : radiusInt0; //if no hasSuperstructure, always large enough
            double radiusInfty = isEnergyRange ? gd.getNextNumber() : radiusInt0;
            double azBlurAngleDegrees = gd.getNextNumber();
            azBlurAngleField.setEnabled(spotBackgrShape == LeedSpotAnalyzer.AZIMUTH_BLUR);
            if (hasSuperstructure) {
                Vector numFields = gd.getNumericFields();
                Checkbox cbx = (Checkbox)(gd.getCheckboxes().get(0));
                if (e != null && gd.getNextBoolean() && (e.getSource() == numFields.get(0) || e.getSource() == cbx)) {
                    ((TextField)numFields.get(1)).setText(((TextField)numFields.get(0)).getText());
                    radiusSup0 = radiusInt0;
                } else if (e != null && e.getSource() == numFields.get(1) && radiusSup0 != radiusInt0)
                    cbx.setState(false);
            }
            boolean valuesOk = isEnergyRange ?
                    radiusInfty >= MIN_RADIUS && radiusInt0 >= radiusInfty && radiusSup0 >= radiusInfty :
                    radiusInt0 >= MIN_RADIUS && radiusSup0 >= MIN_RADIUS;
            if (spotBackgrShape == LeedSpotAnalyzer.AZIMUTH_BLUR && !(azBlurAngleDegrees > 0 && azBlurAngleDegrees <= 30))
                valuesOk = false;
            if (!valuesOk) return false;

            if (e == null && !recordingFinalized) {                //dialog done successfully, write params
                double radius1eVSqrInt = (sqr(radiusInt0) -  sqr(radiusInfty))*energyMin;
                double radius1eVSqrSup = (sqr(radiusSup0) -  sqr(radiusInfty))*energyMin;
                //IJ.log("r("+(int)energyMin+")="+IJ.d2s(radiusInt0)+" r(infty)="+IJ.d2s(radiusInfty)+ " delta r^2(1eV)="+IJ.d2s(radius1eVSqrInt));
                boolean changes = false;
                if (LeedParams.get(LeedParams.RADIUS1EVSQR) != radius1eVSqrInt) changes = true;
                LeedParams.set(LeedParams.RADIUS1EVSQR, radius1eVSqrInt);
                if (hasSuperstructure) {
                    if (LeedParams.get(LeedParams.RADIUS1EVSQRS) != radius1eVSqrSup) changes = true;
                    LeedParams.set(LeedParams.RADIUS1EVSQRS, radius1eVSqrSup);
                }
                LeedParams.set(LeedParams.AZIMUTHBLURANGLE, azBlurAngleDegrees);

                if (radiusInfty > Math.min(radiusInt0, radiusSup0))  //can happen if !isEnergyRange (but even then we need radiusInfty)
                    radiusInfty = Math.min(radiusInt0, radiusSup0);
                if (LeedParams.get(LeedParams.RADIUSINFTYSQR) != sqr(radiusInfty)) changes = true;
                LeedParams.set(LeedParams.RADIUSINFTYSQR, sqr(radiusInfty));
                if (changes)
                    spotTracker.setPostTrackingChanges("Integration radius");
                if (LeedParams.get(LeedParams.BACKGROUNDTYPE) != spotBackgrShape)
                    spotTracker.setPostTrackingChanges("Integration/background area");
                LeedParams.set(LeedParams.BACKGROUNDTYPE, spotBackgrShape);
                recordingFinalized = true;
            }
            return true;
        } catch(Exception ex){IJ.handleException(ex);return false;}
    }

    /** Calculates the integration radius for at the given energy for integer or
     *  hasSuperstructure spots. */
    public static double radius(double energy, boolean isSuperstructure) {
        double radius1eVSqr = LeedParams.get(isSuperstructure ? LeedParams.RADIUS1EVSQRS : LeedParams.RADIUS1EVSQR);
        double radiusInftySqr = LeedParams.get(LeedParams.RADIUSINFTYSQR);
        double r = Math.sqrt(radiusInftySqr + radius1eVSqr/energy);
        return r;
    }

    /** Returns a text describing the radii */
    public static String getStatusText(double energyMin, double energyMax, LeedSpotPattern spotPattern) {
        boolean isEnergyRange = energyMax > energyMin && energyMin != UNKNOWN_ENERGY;
        int spotBackgrShape = (int)LeedParams.get(LeedParams.BACKGROUNDTYPE);
        boolean superstructureDifferent = spotPattern.isSuperstructure() &&
            LeedParams.get(LeedParams.RADIUS1EVSQRS) != LeedParams.get(LeedParams.RADIUS1EVSQR);
        String str = LeedUtils.d2s(LeedRadiusSelector.radius(energyMin, false), 3);
        if (superstructureDifferent)
            str += " (" + LeedUtils.d2s(LeedRadiusSelector.radius(energyMin, true), 3) +  ")";
        if (isEnergyRange)
            str += " -> " + LeedUtils.d2s(LeedRadiusSelector.radius(energyMax, false), 3);
        str += " " + shapesShort[spotBackgrShape];
        if (spotBackgrShape == LeedSpotAnalyzer.AZIMUTH_BLUR)
            str += " " + LeedUtils.d2s(LeedParams.get(LeedParams.AZIMUTHBLURANGLE), 2) + "\u00B0";
        return str;
    }

    static double sqr(double x) {return x*x;}
}
