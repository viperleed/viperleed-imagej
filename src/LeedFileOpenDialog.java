import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.io.*;
import ij.plugin.*;

import ij.util.Tools;
import ij.plugin.frame.Recorder;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.File;

/** This class shows a dialog to open multiple image or movie files, opens them
 *  and sets the Spot Tracker input files accordingly. */

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

public class LeedFileOpenDialog implements Runnable {
    /** shorter notation for Spot Tracker constants */
    static final int N_IMP_FIELDS = LEED_Spot_Tracker.N_IMP_FIELDS;
    static final int DARK = LEED_Spot_Tracker.DARK, FLAT = LEED_Spot_Tracker.FLAT, DARK2 = LEED_Spot_Tracker.DARK2;
    static final String[] IMP_NAMES = LEED_Spot_Tracker.IMP_NAMES;
    /** command to open the files */
    static final String OPEN_COMMAND = Menus.getCommands().get("Open LEED Movie") != null ?
            "Open LEED Movie" : "Open LEED Movie...";
    LEED_Spot_Tracker spotTracker;

    public LeedFileOpenDialog(LEED_Spot_Tracker spotTracker) {
        this.spotTracker = spotTracker;
    }

    /** The dialog runs in its own thread, to enable IJ.showStatus, progress bar update */
    public void run() {
        ImagePlus[] imps = openFiles();
        spotTracker.setOpenFiles(imps);
        spotTracker.enableAndHighlightComponents(true);
    }

    /** Asks the user for the files and returns them as an array,
     *  with array elements corresponding to the Soot Tracker file types
     *  "iv", "dark", "flat", "dark2", and "mask".
     *  Returns null if cancelled. */
    static ImagePlus[] openFiles() {
        GenericDialog gd = new GenericDialog(LEED_Spot_Tracker.PLUGIN_NAME+" - Open Files");
        for (int type = 0; type < N_IMP_FIELDS; type++)
            gd.addFileField(IMP_NAMES[type], "", 65);
        gd.addCheckbox("Use "+IMP_NAMES[DARK]+" also as "+IMP_NAMES[DARK2], true);
        gd.addCheckbox("Virtual Stack", LeedParams.getBoolean(LeedParams.OPENSTACKSASVIRTUAL));
        gd.addMessage("You can drag&drop onto the input fields.\nLeave the input empty to keep the current movie or image");

        DialogListener dialogListener = new DialogListener() {  //for updating the 'use dark as dark2' checkbox
            public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
                Vector textFields = gd.getStringFields();
                boolean useDarkAsDark2 = gd.getNextBoolean();
                String[] paths = new String[N_IMP_FIELDS];
                boolean allOk = true;
                for (int type = 0; type < N_IMP_FIELDS; type++) {
                    paths[type] = gd.getNextString();
                    if (paths[type].length() > 0) {
                        boolean fileOk = LeedUtils.fileOk(paths[type]);
                        if (!fileOk) allOk = false;
                        ((TextField)textFields.get(type)).setForeground(fileOk ? Color.BLACK : Color.RED);
                    }
                }
                Checkbox cbx = (Checkbox)(gd.getCheckboxes().get(0));
                if (paths[DARK].length() > 0 && paths[DARK2].length() > 0 && !paths[DARK].equals(paths[DARK2]))
                    cbx.setState(false);                        //dark1 & dark2 given, but different
                else if (paths[DARK].length() > 0 && paths[DARK].equals(paths[DARK2]))
                    cbx.setState(true);                         //dark1 & dark2 given and equal
                boolean enabled = paths[DARK].length() > 0 && paths[DARK2].length() == 0 && paths[FLAT].length() > 0; //dark given, dark2 empty
                cbx.setEnabled(enabled);
                return allOk;
            }
        };
        Button okButton = gd.getButtons()[0];
        okButton.setEnabled(dialogListener.dialogItemChanged(gd, null));
        gd.addDialogListener(dialogListener);
        gd.showDialog();

        if (gd.wasCanceled()) return null;
        IJ.showStatus("Opening, wait...");

        boolean useDarkAsDark2 = gd.getNextBoolean();
        boolean useVirtualStack = gd.getNextBoolean();
        LeedParams.set(LeedParams.OPENSTACKSASVIRTUAL, useVirtualStack);
        String[] paths = new String[N_IMP_FIELDS];
        for (int type = 0; type < N_IMP_FIELDS; type++)
            paths[type] = gd.getNextString();

        //Check whether already open
        ImagePlus[] alreadyOpenImp = new ImagePlus[N_IMP_FIELDS];
        boolean anyOpen = false;
        int[] openImages = WindowManager.getIDList();
        if (openImages != null) {
            for (int imageID : openImages) {
                ImagePlus imp = WindowManager.getImage(imageID);
                if (imp == null) continue;
                FileInfo fi = imp.getOriginalFileInfo();
                if (fi == null) continue;
                String impPath = fi.directory + fi.fileName;
                File impFile = new File(impPath);
                for (int type = 0; type < N_IMP_FIELDS; type++) {
                    if (paths[type].length() > 0 && impFile.equals(new File(paths[type]))) {
                        alreadyOpenImp[type] = imp;
                        anyOpen = true;
                    }
                }
            }
        }
        if (anyOpen) {
            String str = "The following file(s) are already open in ImageJ:";
            for (int type = 0; type < N_IMP_FIELDS; type++) {
                if (alreadyOpenImp[type] != null) {
                    ImagePlus imp = alreadyOpenImp[type];
                    str += "\n";
                    str += imp.getImageStackSize() > 1 ? "Movie: " : "Image: ";
                    str += imp.getTitle();
                    String newFileName = (new File(paths[type])).getName();
                    if (!newFileName.equals(imp.getTitle()))
                        str += "\nFile: "+newFileName;                    
                    str += "\n   (for use as "+LEED_Spot_Tracker.IMP_NAMES[type]+")";
                }
            }
            str += "\nUse the open File(s)?";
            YesNoCancelDialog dialog = new YesNoCancelDialog(IJ.getInstance(), LEED_Spot_Tracker.PLUGIN_NAME,
                str, "Use Already Open File(s)", "Open Again");
            if (dialog.cancelPressed())
                return null;
            if (!dialog.yesPressed())
                Arrays.fill(alreadyOpenImp, null);
        }
        //Open the image(s)/movie(s) (unless open already)
        String badFiles = "";
        ImagePlus[] imps = new ImagePlus[N_IMP_FIELDS];
        for (int type = 0; type < N_IMP_FIELDS; type++) {
            if (alreadyOpenImp[type] != null) {
                imps[type] = alreadyOpenImp[type];
            } else {
                String path = paths[type];
                if (path.length() == 0) continue;
                if (type == DARK2 && useDarkAsDark2 && imps[DARK] != null) continue; //don't read dark2 if the same as dark
                File file = new File(path);
                ImagePlus previousImp = WindowManager.getCurrentImage();
                if (file.isFile() && file.canRead()) {
                    Recorder.suspendRecording();
                    Open_LEED_Movie opener = new Open_LEED_Movie();
                    ImagePlus imp = opener.openImagePlus(path, useVirtualStack);
                    if (imp != null) {
                        imp.show();
                        imp.waitTillActivated();
                        imps[type] = imp;
                    } else
                        badFiles += "\n" + path;
                    Recorder.resumeRecording();
                    if (Recorder.record && !Recorder.scriptMode() && imp != null) {
                        String str = "run(\""+OPEN_COMMAND+"\", \"open=["+path+"]";
                        if (useVirtualStack)
                            str += " virtual";
                        str += "\");\n";
                        Recorder.recordString(str);
                    }
                }
            }
        }

        if (badFiles.length() > 0)
            IJ.error(LEED_Spot_Tracker.PLUGIN_NAME, "Error: Could not open file(s):"+badFiles);
        if (useDarkAsDark2 && imps[DARK2] == null)
            imps[DARK2] = imps[DARK];

        return imps;
    }
}

