import ij.*;
import ij.IJ;
import ij.gui.*;
import ij.io.*;
import ij.util.Tools;
import java.io.*;
import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;

/**
 *  This class saves data and plots of LEED I/V measurements
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

public class LeedDataSaver implements Runnable {
    static final String TABLE_EXTENSION = ".csv";
    static final String STACK_EXTENSION = ".tif.zip";

    LEED_Spot_Tracker spotTracker;
    String directory;
    String prefix;
    LeedIVAnalyzer ivAnalyzer;
    double[] xAxis;
    String xAxisLabel;
    LeedSpotPattern spotPattern;
    ImagePlus stackImp;
    ImagePlus indexInputImp;
    boolean askToDiscard;
    boolean interactive;
    String errorMessage;
    static String lastStatusText; //remembers text next to Spot Tracker button

    /**
     *  Creator
     */
    public LeedDataSaver(LEED_Spot_Tracker spotTracker, String directory, String prefix, LeedIVAnalyzer ivAnalyzer,
            double[] xAxis, String xAxisLabel, LeedSpotPattern spotPattern, ImagePlus stackImp, ImagePlus indexInputImp, boolean askToDiscard) {
        this.spotTracker = spotTracker;
        this.directory = directory;
        this.prefix = prefix;
        this.ivAnalyzer = ivAnalyzer;
        this.xAxis = xAxis;
        this.xAxisLabel = xAxisLabel;
        this.spotPattern = spotPattern;
        this.stackImp = stackImp;
        this.indexInputImp = indexInputImp;
        this.askToDiscard = askToDiscard;    //true when window gets closed with unsave data
    }

    /** Interactive operation, for running in a separate thread */
    public void run() {
        interactive = true;
        boolean saved = false;
        try {
            Thread.currentThread().setName("LeedDataSaver");
            do {
                if (!saved && askToDiscard) {
                    GenericDialog gd = new GenericDialog(LEED_Spot_Tracker.PLUGIN_NAME);
                    gd.addMessage("Unsaved I(V) Curve Data");
                    gd.enableYesNoCancel("Save", "Discard");
                    gd.hideCancelButton();
                    Button okButton = gd.getButtons()[0];
                    // okButton.requestFocusInWindow(); //does not work because it is not visible
                    gd.showDialog();
                    if (!gd.wasOKed())
                        return;
                }
                String previousDir = OpenDialog.getDefaultDirectory();
                DirectoryChooser.setDefaultDirectory(this.directory);
                DirectoryChooser dChooser = new  DirectoryChooser("Select where to save files");
                directory = dChooser.getDirectory();
                DirectoryChooser.setDefaultDirectory(previousDir);

                if (directory != null) {
                    LeedParams.setString(LeedParams.SAVEDIRECTORY, directory);
                    spotTracker.setButtonLabelText(LEED_Spot_Tracker.SAVE_BUTTON, "saving, please wait...");

                    saved = askAndSaveData(directory);
                    spotTracker.enableAndHighlightComponents(true);
                    lastStatusText = saved ? "Saved ("+prefix+")" : "";
                    spotTracker.setButtonLabelText(LEED_Spot_Tracker.SAVE_BUTTON, lastStatusText);
                }
            } while (askToDiscard && !saved);
        } catch (Exception e) {
            IJ.handleException(e);
        }
        if (saved)
            spotTracker.setDataSaved(directory, prefix);
        else
            spotTracker.enableAndHighlightComponents(true);
    }

    /** Non-interactive operation.
     *  Returns null if successful, otherwise an error text starting with ERROR: */
    public String runInMacro() {
        boolean saved = askAndSaveData(directory);
        if (saved)
            spotTracker.setDataSaved(directory, prefix);
        return saved ? null :
                "ERROR: " + (errorMessage== null ? "Could not save data" : errorMessage);
    }

    /** Returns the 'saved' text of the last saving */
    public static String getLastStatusText() {
        return lastStatusText;
    }

    boolean askAndSaveData(String directory) {
        boolean savePlots = LeedParams.get(LeedParams.SAVEPLOTS) != 0;
        boolean saveImageStack = LeedParams.get(LeedParams.SAVESTACK) != 0;
        if (!directory.endsWith(File.separator) && !directory.endsWith("/"))
            directory += "/";

        if (interactive) {
            GenericDialog gd = new GenericDialog("Save Options");
            gd.addStringField("File prefix", prefix, 30);
            gd.addCheckbox("Save Plots", savePlots);
            gd.addCheckbox("Save Image Stack", saveImageStack);
            gd.addHelp(LeedSpotTrackerHelp.getSaveDataHelp());
            gd.showDialog();
            if (gd.wasCanceled()) return false;
            prefix = gd.getNextString();
            savePlots = gd.getNextBoolean();
            saveImageStack = gd.getNextBoolean();
            LeedParams.set(LeedParams.SAVEPLOTS, savePlots);
            LeedParams.set(LeedParams.SAVESTACK, saveImageStack);
        }

        int nDataTypes = ivAnalyzer.getNDataTypes();
        int nFiles = nDataTypes;                            // create list of files
        nFiles++; //log file
        if (savePlots) nFiles += LeedIVAnalyzer.PLOT_NAMES.length;
        nFiles++; //image stack or single image indexInputImp
        String[] filenames = new String[nFiles];
        for (int i=0; i<nDataTypes; i++)
            filenames[i] = prefix+"_"+ivAnalyzer.getDataName(i)+TABLE_EXTENSION;
        filenames[nDataTypes] = prefix+"_log.txt";
        if (savePlots)
            for (int i=0; i<LeedIVAnalyzer.PLOT_NAMES.length; i++) {
                filenames[nDataTypes+i+1] = prefix + "_" + LeedIVAnalyzer.PLOT_NAMES[i];
                if (!filenames[nDataTypes+i+1].endsWith(".zip"))
                    filenames[nDataTypes+i+1] += ".zip";
        }
        String stackOrSlice = saveImageStack ? "_stack" : "_spotIndicesImage";
        filenames[filenames.length-1] = prefix+stackOrSlice+STACK_EXTENSION;
        ArrayList<String> existingFiles = new ArrayList<String>();
        for (String filename : filenames) {
            File file = new File(directory + filename);
            if (file.exists())
                existingFiles.add(filename);
        }
        if (interactive && existingFiles.size() > 0) {      //ask whether to overwrite
            String message = "Existing files in\n"+directory+":";
            for (String filename : existingFiles)
                message += "\n - " + filename;
            message += "\n\nOverwrite?";
            GenericDialog gd1 = new GenericDialog("File(s) Exist");
            gd1.addMessage(message);
            gd1.enableYesNoCancel("Overwrite", null);
            gd1.showDialog();
            if (gd1.wasCanceled())
                return false;
        }
        spotTracker.enableAndHighlightComponents(false);    //saving can take some time (esp. for the stack)
        for (int i=0; i<nDataTypes; i++)                    //save tables
            if (!saveDataFile(i, directory, filenames[i])) {
                spotTracker.enableAndHighlightComponents(true);
                return false;
            }
        if (!saveLogfile(directory, filenames[nDataTypes])) //save log
            return false;
        if (savePlots) {                                    //create and save plots
            for (int i=0; i<LeedIVAnalyzer.PLOT_NAMES.length; i++) {
                ImagePlus imp = ivAnalyzer.getPlotImp(i);
                if (imp == null) continue;
                String path = directory+filenames[nDataTypes+i+1];
                IJ.saveAs(imp, "zip", path);
                File file = new File(path);
                if (!file.exists() || file.length()<100) {    //seems saving has failed
                    error("Saving failed\n"+path+"\n(saved file missing or too small)");
                    return false;
                }
            }
        }
        ImagePlus imageToSave = saveImageStack ? stackImp : indexInputImp; //save image stack or one slice
        if (imageToSave != null && imageToSave.getProcessor() != null) {
            IJ.run(imageToSave, "Enhance Contrast", "saturated=1.0");
            String path = directory+filenames[filenames.length-1];  // (The big one last. In case the disk will be full, we have the rest)
            IJ.showStatus(saveImageStack ? "Saving image stack..." : "Saving sample image...");
            IJ.saveAs(imageToSave, "zip", path);
            File file = new File(path);
            if (!file.exists() || file.length()<1000) {    //seems saving has failed
                error("Saving failed\n"+path);
                return false;
            }
            IJ.showStatus("Saving done");
        } else
            LeedUtils.logError("Error: Cannot save image or stack (input prematurely closed?)");
        return true;
    }

    void error(String message) {
        if (interactive)
            IJ.error(LEED_Spot_Tracker.PLUGIN_NAME, "Error: "+message);
        else
            errorMessage = message;
    }

    boolean saveDataFile(int dataType, String directory, String filename) {
        PrintWriter pw = null;
        try {
            FileOutputStream fos = new FileOutputStream(directory+filename);
            BufferedOutputStream bos = new BufferedOutputStream(fos);
            pw = new PrintWriter(bos);

            //IJ.log("datasaver IVAnalyzer="+ivAnalyzer+" .hassp="+ivAnalyzer.hasSpot);
            pw.print(xAxisLabel);
            for (int spot=0; spot<spotPattern.size(); spot++) {
                if (!ivAnalyzer.hasSpot(spot)) continue;
                pw.print(",");
                pw.print(spotPattern.getNameWithGroup(spot, true));
            }
            pw.println();

            for (int i=0; i<xAxis.length; i++) {
                pw.print(IJ.d2s(xAxis[i]));
                for (int spot=0; spot<spotPattern.size(); spot++) {
                    if (!ivAnalyzer.hasSpot(spot)) continue;
                    pw.print(',');
                    double v = ivAnalyzer.getData(dataType, spot)[i];
                    if (!Double.isNaN(v))
                        pw.print((float)v);
                }
                pw.println();
            }

            pw.close();
            return true;
        } catch (Exception e) {
            IJ.error("Error writing file "+filename+"\n"+e);
            //IJ.handleException(e);
            if (pw != null)
                try {pw.close();} catch (Exception e2) {}
            return false;
        }
    }

    boolean saveLogfile(String directory, String filename) {
        PrintWriter pw = null;
        try {
            double[] eMinMax = Tools.getMinMax(xAxis);
            double energyMin = eMinMax[0], energyMax = eMinMax[1];

            FileOutputStream fos = new FileOutputStream(directory+filename);
            BufferedOutputStream bos = new BufferedOutputStream(fos);
            pw = new PrintWriter(bos);
            for (String str : spotTracker.getLogLines())
                pw.println(str);
            pw.close();
            return true;
        } catch (Exception e) {
            IJ.error("Error writing file "+filename+"\n"+e);
            //IJ.handleException(e);
            if (pw != null)
                try {pw.close();} catch (Exception e2) {}
            return false;
        }
    }

    static double sqr(double x) {return x*x;}
}
