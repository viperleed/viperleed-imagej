import ij.*;
import ij.IJ;
import ij.util.Tools;
import ij.gui.GenericDialog;
import ij.gui.DialogListener;
import java.util.*;
import java.awt.*;


/**
 * This class contains dialogs for setting/unsetting/deleting highlighted spots,
 * called from the More>> menu of the Spot Tracker.
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
 *  @author Michael Schmid, IAP/TU Wien, 2020-2024
 */

public class LeedHighlightDialogs {
    static final String[] DELETE_OPTIONS = new String[] {" at all energies", " in energy range"};
    static final int DELETE_ALL=0, DELETE_RANGE=1;
    static final String[] ADD_REMOVE_OPTIONS = new String[] {"New", "add to highlighted", " remove from highlighted"};
    static final int NEW=0, ADD=2, REMOVE=2;
    static String highlightedSpotsStr = "";    //remember dialog entries
    static boolean highlightEquivalent;
    static int lastAddRemove; //0=NEW, etc.
    static int deleteOption;  //0=everywhere, 1=energy range
    static double deleteEmin = 0, deleteEmax = 100;

    /** Shows a menu which spots to highlight and them marks them */
    public static void highlightSpots(ImagePlus stackImp, final LeedSpotPattern spotPattern) {
        HashSet<String> highlighted = LeedOverlay.getHighlighted(stackImp);
        final HashSet<String> currentlyHighlighted =
                highlighted == null || highlighted.size() == 0 ? null :  highlighted;
        GenericDialog gd = new GenericDialog(LEED_Spot_Tracker.PLUGIN_NAME+" - Highlight beams");
        gd.addStringField("Beams(s) to highlight *", highlightedSpotsStr, 12);
        gd.addCheckbox("+ symmetry-equivalent beams", highlightEquivalent);
        gd.addMessage("* separate multiple spots by semicolon\nE.g. '1/2, 0; 2|2'");
        if (currentlyHighlighted != null)
            gd.addRadioButtonGroup("", ADD_REMOVE_OPTIONS, 1, ADD_REMOVE_OPTIONS.length, ADD_REMOVE_OPTIONS[lastAddRemove]);
        gd.addMessage("");
        final Label errorLabel = (Label)gd.getMessage();
        DialogListener dialogListener = new DialogListener() {
            public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
                boolean doRemove = false;
                if (currentlyHighlighted != null) {
                    String optionStr = gd.getNextRadioButton();
                    if (optionStr.equals(ADD_REMOVE_OPTIONS[REMOVE]))
                        doRemove = true;
                }
                highlightedSpotsStr = gd.getNextString();
                String[] spotStrs = Tools.split(highlightedSpotsStr, ";");
                int[] spotIndices = spotPattern.getIndices(spotStrs);
                String error = spotPattern.getErrorMessage();
                if (error != null) {
                    errorLabel.setText(error);
                    return false;
                }
                if (doRemove) {
                    for (int iSpot : spotIndices) {
                        if (!currentlyHighlighted.contains(spotPattern.getNameForOvly(iSpot))) {
                            errorLabel.setText("Not highlighted: "+spotPattern.getName(iSpot));
                            return false;
                        }
                    }
                }
                errorLabel.setText("");
                return spotIndices.length > 0;
            }
        };
        Button okButton = gd.getButtons()[0];
        okButton.setEnabled(dialogListener.dialogItemChanged(gd, null));
        gd.addDialogListener(dialogListener);
        gd.showDialog();

        if (gd.wasCanceled()) return;
        highlightedSpotsStr = gd.getNextString();
        highlightEquivalent = gd.getNextBoolean();
        int newAddRemove = NEW;
        if (currentlyHighlighted != null) {
            String optionStr = gd.getNextRadioButton();
            newAddRemove = LeedUtils.arrayIndexOf(ADD_REMOVE_OPTIONS, optionStr); //can be also -1 if "null" selected
            if (newAddRemove < 0) newAddRemove = NEW;
            lastAddRemove = newAddRemove;
        }
        String[] spotStrs = Tools.split(highlightedSpotsStr, ";");
        int[] spotIndices = spotPattern.getIndices(spotStrs);
        HashSet<String> spotNames = newAddRemove == NEW ? new HashSet<String>() : currentlyHighlighted;
        for (int iSpot : spotIndices) {
            if (iSpot < 0) continue;
            if (highlightEquivalent) {
                int group = spotPattern.getGroup(iSpot);
                int[] spots = spotPattern.getAllSpotsForGroup(group);
                for (int s : spots) {
                    if (newAddRemove == REMOVE)
                        spotNames.remove(spotPattern.getNameForOvly(s));
                    else
                        spotNames.add(spotPattern.getNameForOvly(s));
                }
            } else {
                if (newAddRemove == REMOVE)
                    spotNames.remove(spotPattern.getNameForOvly(iSpot));
                else
                    spotNames.add(spotPattern.getNameForOvly(iSpot));
            }
        }
        LeedOverlay.highlightCircles(stackImp, spotNames, 1);
        stackImp.draw();
    }

    /** Shows a menu whether to delete the highlighted spot(s). The return String, if not null, is intended for the log. */
    public static String deleteHighlightedSpots(LEED_Spot_Tracker spotTracker, ImagePlus stackImp, LeedSpotPattern spotPattern, LeedIVAnalyzer ivAnalyzer) {
        if (stackImp==null || spotPattern==null || ivAnalyzer==null) return null; //should never happen
        HashSet<String> highlightedSpots = LeedOverlay.getHighlighted(stackImp);
        if (highlightedSpots == null) return null; //should never happen
        int[] highlightedIndices = spotPattern.getIndices(highlightedSpots.toArray(new String[0]));
        String error = spotPattern.getErrorMessage();
        if (error != null) {IJ.error("Delete Selected Beams - Internal Error", "Spot name(s) not found\n"+error); return null;}

        GenericDialog gd = new GenericDialog(LEED_Spot_Tracker.PLUGIN_NAME+" - Delete highlighted beams");
        final ArrayList<Component> numFieldsLabels = new ArrayList<Component>(4);
        String whatToDelete = highlightedIndices.length == 1 ?
                "the "+spotPattern.getName(highlightedIndices[0]) + " beam." :
                highlightedIndices.length +" beams.";
        gd.addRadioButtonGroup("Delete "+whatToDelete, DELETE_OPTIONS, 2, 0, DELETE_OPTIONS[deleteOption]);
        gd.addNumericField("From energy", deleteEmin, 1);
        numFieldsLabels.add(gd.getLabel());
        numFieldsLabels.add((Component)gd.getNumericFields().lastElement());
        gd.addNumericField("to", deleteEmax, 1);
        numFieldsLabels.add(gd.getLabel());
        numFieldsLabels.add((Component)gd.getNumericFields().lastElement());
        gd.addMessage("You cannot undo deleting (except by 'Track Spots')");
        DialogListener dialogListener = new DialogListener() {
            public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
                CheckboxGroup radioBG = (CheckboxGroup)gd.getRadioButtonGroups().get(0);
                Checkbox cbx = radioBG.getSelectedCheckbox();
                if (cbx == null) return false;
                boolean deleteRange = cbx.getLabel().equals(DELETE_OPTIONS[DELETE_RANGE]);
                for (Component c : numFieldsLabels)
                    c.setEnabled(deleteRange);
                if (deleteRange) {
                    double deleteEmin = gd.getNextNumber();
                    double deleteEmax = gd.getNextNumber();
                    if (!(deleteEmin >=0 && deleteEmax > deleteEmin)) return false;
                }
                return true;
            }
        };
        Button okButton = gd.getButtons()[0];
        okButton.setEnabled(dialogListener.dialogItemChanged(gd, null));
        gd.addDialogListener(dialogListener);
        gd.showDialog();

        if (gd.wasCanceled()) return null;
        String optionStr = gd.getNextRadioButton();
        deleteOption = LeedUtils.arrayIndexOf(DELETE_OPTIONS, optionStr);
        String deletedS = "Deleted ";
        HashSet<String> fullyDeleted = new HashSet<String>();
        if (deleteOption == DELETE_ALL) {
            for (int index : highlightedIndices) {
                ivAnalyzer.hasSpot[index] = false;
                deletedS += spotPattern.getName(index)+"; ";
                fullyDeleted.add(spotPattern.getNameForOvly(index));
            }
            LeedOverlay.highlightCircles(stackImp, null, 1);
        } else {  //DELETE_RANGE
            deleteEmin = gd.getNextNumber();
            deleteEmax = gd.getNextNumber();
            if (!(deleteEmax > deleteEmin)) {
                IJ.error("Delete Selected Beams - Error", "Invalid energy range "+IJ.d2s(deleteEmin,1)+"-"+IJ.d2s(deleteEmax,1)+
                        "\nNothing was deleted.");
                return null;
            }
            double[] energies = ivAnalyzer.energies;
            int nPartiallyDeleted = 0, nFullyDeleted = 0;
            String deleted2S = "Deleted "+IJ.d2s(deleteEmin,1)+"-"+IJ.d2s(deleteEmax,1)+": ";
            for (int index : highlightedIndices) {
                double[] integrals = ivAnalyzer.data[LeedIVAnalyzer.INTEGRAL][index];
                double[] intCorrI0 = ivAnalyzer.data[LeedIVAnalyzer.INT_I0CORRECTED][index];
                int nGood = 0, nDeleted = 0;
                for (int i=0; i<energies.length; i++) {
                    if (Double.isNaN(integrals[i])) continue;
                    if (energies[i] >= deleteEmin && energies[i] <= deleteEmax) {
                        integrals[i] = Double.NaN;
                        intCorrI0[i] = Double.NaN;
                        nDeleted++;
                    } else
                        nGood++;
                }
                if (nGood < ivAnalyzer.minPointsPerBeam) {
                    ivAnalyzer.hasSpot[index] = false;
                    nFullyDeleted++;
                    deletedS += spotPattern.getName(index)+"; ";
                    fullyDeleted.add(spotPattern.getNameForOvly(index));
                } else {
                    nPartiallyDeleted++;
                    deleted2S += spotPattern.getName(index)+"; ";
                }
            }

            if (nPartiallyDeleted == 0)
                LeedOverlay.highlightCircles(stackImp, null, 1); //nothing remains highlighted
            else
                deletedS = nFullyDeleted == 0 ? deleted2S : deletedS + '\n' + deleted2S;
        }
        if (fullyDeleted.size() > 0) {
            LeedOverlay.delete(stackImp, fullyDeleted);
            stackImp.draw();
        }
        spotTracker.setButtonLabelText(LEED_Spot_Tracker.TRACK_BUTTON, ivAnalyzer.getStatusText(/*verbose=*/false));
        return deletedS;
    }
}
