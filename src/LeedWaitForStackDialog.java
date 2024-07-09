import ij.*;
import ij.gui.*;
import java.awt.*;


/**
 *  This class displays a dialog asking the user to wait while the stack is updated
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
 *  @author Michael Schmid, IAP/TU Wien, 2024
 */

public class LeedWaitForStackDialog {

    /** While the stack is updating, shows a dialog asking the user for patience.
     *  Returns true when the stack is ready and false if the user has pressed 'cancel'.
     *  To avoid a deadlock, this function must NOT be called in the EventQueue. */
    public static boolean wait(final LEED_Spot_Tracker spotTracker, final LeedDarkFlatVirtualStack stack) {
        if (EventQueue.isDispatchThread()) {
            IJ.error(LEED_Spot_Tracker.PLUGIN_NAME, "Internal Error: WaitForStackDialog must not be called in the Event Queue");
            return false;
        }
        final GenericDialog gd = new GenericDialog(LEED_Spot_Tracker.PLUGIN_NAME+" - Wait...");
        gd.addMessage("Please wait until the calculations to prepare\nthe 'Spot Tracking' stack are finished.");
        gd.getButtons()[0].setVisible(false);       //no OK button, only CANCEL
        EventQueue.invokeLater(new Runnable() {public void run() {
                    gd.showDialog();
                }});
        while (!gd.isVisible())                 //wait until the dialog becomes visible
            IJ.wait(10);
        while (true) {                              //wait until...
            if (gd.wasCanceled())                   //the user cancels
                return false;
            if (!stack.isChanging) {                //or the update is done
                gd.dispose();
                return true;
            }
            IJ.wait(50);
        }
    }
}
