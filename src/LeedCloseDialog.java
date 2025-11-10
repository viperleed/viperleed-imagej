import ij.*;
import ij.IJ;
import ij.util.Tools;
import ij.gui.GenericDialog;
import ij.gui.DialogListener;
import ij.gui.ImageWindow;
import ij.plugin.frame.Recorder;
import java.util.*;
import java.awt.*;


/**
 *  This class contains dialogs whether to close the plots and the stack for I(V) measurement
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
 *  @author Michael Schmid, IAP/TU Wien, 2020-2025
 */



public class LeedCloseDialog {
    static final String CLOSE_S_KEY = LeedParams.PREFS_KEY+".closeS";  //in Prefs, for "close stack"
    static final String CLOSE_P_KEY = LeedParams.PREFS_KEY+".closeP";  //in Prefs, for "close plots"
    static final String EXCEPT_R_KEY = LeedParams.PREFS_KEY+".keepR";   //in Prefs, for "except R factor statistics"
    static final String EXCEPT_Q_KEY = LeedParams.PREFS_KEY+".keepQ";   //in Prefs, for "except selected I(V) curves of equivalent beams"
    static final String EXCEPT_N_KEY = LeedParams.PREFS_KEY+".keepN";   //in Prefs, for "except new ones"

    /** Shows a dialog what to close and closes accordingly.
     *  When the main Spot Tracker panel closes, stackImp is the
     *  virtual stack that the spor tracker works on.
     *  Otherwise stackImp = null. */
    public static void showDialogAndClose(LEED_Spot_Tracker spotTracker, final ImagePlus stackImp) {
        final boolean anyOpenPlots = LeedIVAnalyzer.anyOpenPlots((char)0);
        final boolean anyOpenPlotsR = LeedIVAnalyzer.anyOpenPlots('R');
        final boolean anyOpenPlotsQ = LeedIVAnalyzer.anyOpenPlots('Q');
        if (!anyOpenPlots && (stackImp==null || stackImp.getWindow()==null))
            return;                                             //nothing to close
        if (Recorder.record)
            Recorder.resetCommandOptions();
        GenericDialog gd = new GenericDialog(LEED_Spot_Tracker.PLUGIN_NAME);
        gd.addMessage("Close the following?");
        if (stackImp != null) {
            gd.addCheckbox("Stack "+ stackImp.getTitle(), !anyOpenPlots || Prefs.get(CLOSE_S_KEY, true));
            if (anyOpenPlots)
                gd.addCheckbox("Plots", Prefs.get(CLOSE_P_KEY, true));
        } else
            gd.addMessage("Plots");
        if (anyOpenPlots) {
            if (anyOpenPlotsR) {
                gd.setInsets(0, 30, 0); //top, left (default 20), bottom
                gd.addCheckbox("Except_quality statistics", Prefs.get(EXCEPT_R_KEY, false));
            }
            if (anyOpenPlotsQ) {
                gd.setInsets(0, 30, 0);
                gd.addCheckbox("Except_selected I(V) curves", Prefs.get(EXCEPT_Q_KEY, false));
            }
            gd.setInsets(0, 30, 0);
            gd.addCheckbox("Except_newest", Prefs.get(EXCEPT_N_KEY, false));
        }
        gd.setOKLabel("Close");
        gd.setCancelLabel("close nothing");
        DialogListener dialogListener = new DialogListener() {
            public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
                int nBooleanRead = 0;
                boolean canClose = stackImp == null;
                if (stackImp != null) {
                    boolean closeStack = gd.getNextBoolean();
                    nBooleanRead++;
                    if (closeStack)
                        canClose = true;
                    if (anyOpenPlots) {
                        boolean closePlots = gd.getNextBoolean();
                        nBooleanRead++;
                        if (closePlots)
                            canClose = true;
                        Vector checkboxes = gd.getCheckboxes();
                        for (int i=2; i<checkboxes.size(); i++)
                            ((Checkbox)checkboxes.get(i)).setEnabled(closePlots);
                    }
                    String oldOkLabel = gd.getButtons()[0].getLabel();
                    gd.getButtons()[0].setLabel(canClose ? "Close" : "OK");
                    if ("OK".equals(oldOkLabel) && canClose) gd.pack(); //'Close' needs more space
                }
                if (Recorder.record)
                    for (int i=nBooleanRead; i<gd.getCheckboxes().size(); i++)
                        gd.getNextBoolean();    //read remaining checkboxes so they are macro recorded
                return true;
            }
        };
        Button okButton = gd.getButtons()[0];
        okButton.setEnabled(dialogListener.dialogItemChanged(gd, null));
        gd.addDialogListener(dialogListener);
        gd.showDialog();

        if (gd.wasCanceled()) return;

        boolean closePlots = true;
        if (stackImp != null) {
            boolean closeStack = gd.getNextBoolean();
            Prefs.set(CLOSE_S_KEY, closeStack);
            if (closeStack) {
                ImageWindow win = stackImp.getWindow();
                if (win != null) {
                    win.removeWindowListener(spotTracker);
                    Point loc = win.getLocation();
                    Prefs.saveLocation(LEED_Spot_Tracker.LOC_KEY_S, loc);
                }
                stackImp.close();
            }
            if (anyOpenPlots) {
                closePlots = gd.getNextBoolean();
                Prefs.set(CLOSE_P_KEY, closePlots);
            }
        }
        if (anyOpenPlots && closePlots) {
            String exceptThese = "";
            boolean exceptR = anyOpenPlotsR ? gd.getNextBoolean() : false;
            if (anyOpenPlotsR) Prefs.set(EXCEPT_R_KEY, exceptR);
            if (exceptR) exceptThese += 'R';
            boolean exceptQ = anyOpenPlotsQ ? gd.getNextBoolean() : false;
            if (anyOpenPlotsQ) Prefs.set(EXCEPT_Q_KEY, exceptQ);
            if (exceptQ) exceptThese += 'Q';
            boolean exceptN = gd.getNextBoolean();
            Prefs.set(EXCEPT_N_KEY, exceptN);
            if (exceptN) exceptThese += 'N';
            LeedIVAnalyzer.closeAllPlots(exceptThese);
        }
        if (Recorder.record) {
            String options = Recorder.getCommandOptions();
            if (stackImp != null) options = "tracker "+options;
            LEED_Spot_Tracker.recordMacro(LEED_Spot_Tracker.MACRO_COMMAND_CLOSE, '"'+options+'"');
        }
    }
}
