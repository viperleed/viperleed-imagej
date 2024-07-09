import ij.plugin.*;
import ij.*;
import ij.gui.*;
import ij.io.*;
import ij.util.Tools;
import ij.measure.*;
import ij.process.*;
import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.zip.*;
import java.text.SimpleDateFormat;

/**
 * This ImageJ plugin opens a ViPErLEED I(V) zip archive, an
 * 'AIDA' (Automatic Image and Data Acquisition) EE2000/EE2010 LEED movie (*.vid format),
 * or an image type readable by ImageJ (except zip).
 * 
 * Reading a ViPErLEED I(V) zip archive is fully handled by this class.
 * AIDA files are read by a separate class, and other file types by ImageJ.
 * 
 * A ViPErLEED I(V) zip archive must contain the images and a .csv index file.
 * The contents of the csv index file must be a line with non-numeric column headers,
 * the list of tiff image file names in the first column, and further data
 * such as energy and beam current (non-numeric data are also possible, but ignored).
 * The csv file may contain empty lines and comment lines starting with any of these characters: #%!*
 * Non-blank comment lines at the top (excluding the comment character) are added to the image info.
 */

/** This ImageJ plugin is part of the ViPErLEED package for LEED I(V) analysis.
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

public class Open_LEED_Movie extends VirtualStack implements PlugIn {
    final static String PLUGIN_NAME = "Open LEED Movie";
    final static String MAGIC_ZIP = "PK"; //zip file starts with these two ascii characters
    int nFrames;
    int width, height; //private in ImageStack; thus we can't use the superclass variables
    File inFile;
    ZipEntry[] zipEntries; //zip entries for the images
    String[] sliceLabels;
    StringBuilder imageInfo = new StringBuilder();
    //dialog parameters. These are also used by Open_AIDA_Video
    public static boolean showTable = false;
    public static boolean hideDialog = false;

    /** The plugin is invoked by ImageJ using this method.
     *  @param arg   String 'arg' may be used to select the path. If it is null or
     *  not a path to a readable file, a file open dialog is shown.
     */
    public void run(String arg) {
        if (IJ.altKeyDown()) hideDialog = false;
        boolean useVirtualStack = LeedParams.getBoolean(LeedParams.OPENSTACKSASVIRTUAL);
        ImagePlus imp = openImagePlus(arg, useVirtualStack);
        if (imp == null)
            return;
        IJ.run(imp, "Enhance Contrast", "saturated=0.35");
        imp.show();
        if (IJ.isMacro())
            imp.waitTillActivated();
    }

   /** Opens the image or image stack as ImagePlus, without displaying it,
    *  and returns a reference to it.
    *  @param path   Path to the file. If it is null or
    *        not a path to a readable file, a file open dialog is shown.
    *  @return The ImagePlus or null on cancel or error
    */
    public ImagePlus openImagePlus(String path, boolean useVirtualStack) {
        
        if (path.length() == 0) path = null;
        if (path != null) {
            File file = new File(path);
            if (!file.canRead())
                path = null;
        }
        boolean interactive = path == null;
        if (interactive) {
            OpenDialog od = new OpenDialog("Open LEED Movie or Image", "");
            String fileName = od.getFileName();
            if (fileName == null) {
                if (IJ.isMacro()) IJ.error(PLUGIN_NAME, "No file selected");
                return null;
            }
            String directory = od.getDirectory();
            path = directory + fileName;
        }
        inFile = new File(path);
        if (!inFile.canRead()) {
            if (IJ.isMacro()) IJ.error(PLUGIN_NAME, "File unreadable or missing:\n"+path);
            return null;
        }
        if (LeedUtils.hasMagic(inFile, Open_Aida_LEED_Video.FILE_START))
            return new Open_Aida_LEED_Video().openImagePlus(path, useVirtualStack, showTable,
                    hideDialog || !interactive);
        if (!LeedUtils.hasMagic(inFile, MAGIC_ZIP)) // not a zip file
            return IJ.openImage(path);
        Object table = readTable(inFile);           // If a zip file, read the csv table and populate zipEntries list
        if (table instanceof File)                  // it's a zip file, but seemingly an ImageJ zip
            return IJ.openImage(path);
        if (table == null)                          // Reading table failed
            return null;

        ResultsTable rt = (ResultsTable)table;      // We got a LEED zip file
        imageInfo.append("ViPErLEED movie\n");
        ImageProcessor ip = nFrames > 0 ? getProcessor(1) : null;
        if (ip == null) return null;
        setBitDepth(ip.getBitDepth());
        width = ip.getWidth();
        height = ip.getHeight();
        if (width==0 || height==0) {
            IJ.error("Invalid image size: "+width+"x"+height+":\nFile "+path+
                    "\nImage "+zipEntries[0].getName());
            return null;
        }
        String fileName = inFile.getName();

        if (interactive && !hideDialog) {                       // dialog
            GenericDialog gd = new GenericDialog("Open LEED .zip File...");
            double sizeMB = width*height*getBitDepth()/(1024.*1024.*8)*nFrames;
            gd.addMessage(fileName+"\n"+nFrames+" frames, "+width+"x"+height+" ("+IJ.d2s(sizeMB, 1)+" MB)");
            gd.addCheckbox("Virtual Stack", useVirtualStack);
            gd.addCheckbox("Show Parameter Table", showTable);
            gd.addCheckbox("Don't Show This Dialog Again *", hideDialog);
            gd.addMessage("* Invoke with ALT key to display the dialog again");
            gd.showDialog();
            if (gd.wasCanceled()) return null;
            useVirtualStack = gd.getNextBoolean();
            showTable = gd.getNextBoolean();
            hideDialog = gd.getNextBoolean();
            LeedParams.set(LeedParams.OPENSTACKSASVIRTUAL, useVirtualStack);
        }
        if (showTable && interactive) {
            String tableName = LeedUtils.removeExtension(fileName)+"_T";
            rt.show(tableName);
        }
        ImageStack stack = makeStack(inFile, rt, useVirtualStack);
        if (stack == null) return null;
        ImagePlus imp = new ImagePlus(WindowManager.makeUniqueName(fileName), stack);
        imp.setProperty("Info", imageInfo.toString());
        setCalibration(imp, rt);
        FileInfo fi = new FileInfo();                           // parameters for 'info'
        fi.fileName = fileName;
        fi.directory = inFile.getParent();
        if (!(fi.directory.endsWith("/") || fi.directory.endsWith("\\")))
            fi.directory += File.separator;
        fi.fileFormat = FileInfo.UNKNOWN;
        fi.width = width;
        fi.height = height;
        fi.nImages = stack.size();
        imp.setFileInfo(fi);

        return imp;
    }

    /** Tries to read the first .csv from the zip file and returns it as a ResultsTable.
     *  Also populates the zipEntries if we find a table.
     *  Returns null on error, and the file if it looks like an ImageJ zip file */
    private Object readTable(File inFile) {
        ArrayList<ZipEntry> entries = getZipEntries(inFile);
        if (entries == null)
            return null;
        if (entries.size() == 1 && entries.get(0).getName().matches(".*\\.tiff?"))    //looks like an ImageJ zip file
            return inFile;
        ResultsTable rt = new ResultsTable();
        for (ZipEntry entry : entries) {    //find the (first) csv file
            if (entry.getName().endsWith(".csv")) {
                ZipFile zipFile = null;
                int n=0;
                try {
                    zipFile = new ZipFile(inFile, ZipFile.OPEN_READ);
                    BufferedReader reader = new BufferedReader(new InputStreamReader(zipFile.getInputStream(entry)));
                    String[] headings = null;
                    do {                            // read comments (if any) and column headings
                        String line = reader.readLine().trim();
                        if (line == null)
                            throw new RuntimeException("Error: Empty table in "+entry.getName());
                        if (line.trim().length() == 0) continue; //ignore empty lines
                        if (LeedUtils.isCommentLine(line)) {
                            line = line.substring(1).trim();
                            if (line.length() > 0) {
                                imageInfo.append(line);
                                imageInfo.append('\n');
                            }
                            continue;
                        }
                        headings = line.split(","); // not a comment or empty line: it must be the headings
                        if (headings.length < 2)
                            throw new RuntimeException("Error: Not enough columns in "+entry.getName()+":\n"+line);
                        for (int i=0; i<headings.length; i++) {
                            headings[i] = headings[i].trim();
                            if (!Double.isNaN(Tools.parseDouble(headings[i])))
                            throw new RuntimeException("Error: Invalid (numeric) column name in "+entry.getName()+":\n"+line);
                            headings[i] = headings[i].replaceAll("\\(.*\\)",""); //delete anything in parentheses
                            int suffix=1;
                            boolean isDuplicateHeading = false;
                            do {
                                int index = LeedUtils.arrayIndexOf(headings, headings[i]);
                                isDuplicateHeading = index < i && index >= 0;
                                if (isDuplicateHeading) {
                                    if (headings[i].matches(".*_[1-9]"))
                                        headings[i] = headings[i].substring(0, headings[i].length()-2);
                                    headings[i] += "_"+suffix;
                                    suffix++;
                                }
                            }
                            while (isDuplicateHeading);
                        }
                    } while (headings == null);
                    zipEntries = new ZipEntry[entries.size()-1];
                    while (true) {                  // read data
                        String line = reader.readLine();
                        if (line == null) break;
                        if (LeedUtils.isCommentLine(line)) continue;
                        String[] items = line.split(",");
                        if (items.length != headings.length)
                            throw new RuntimeException("Error: wrong number of items in "+entry.getName()+":\n"+line);
                        String imgName = items[0].trim();
                        ZipEntry imgEntry = findEntry(entries, imgName);
                        if (imgEntry == null) {
                            LeedUtils.logError("ERROR: Cannot find image file "+imgName+" in zip");
                            //DEBUG String e20=entries.get(20).getName();
                            //DEBUG IJ.log("len="+imgName.length()+" entry20="+e20+"; len="+e20.length());
                            continue;
                        }
                        //DEBUG else IJ.log(""+imgEntry.getTime());
                        zipEntries[n] = imgEntry;
                        rt.addRow();
                        rt.setLabel(imgName, n);
                        for (int i=1; i<Math.min(headings.length, items.length); i++) {
                            String item = items[i].trim();
                            double value = Tools.parseDouble(item);
                            if (Double.isNaN(value))
                                rt.addValue(headings[i], item);
                            else
                                rt.addValue(headings[i], value);
                        }
                        n++;
                    }
                } catch (Exception e) {
                    IJ.handleException(e);
                    rt = null;
                } finally {
                    try {
                        zipFile.close();
                    } catch (Exception ee) {}
                }
                nFrames = n;
                if (nFrames == 0) {
                    IJ.error("Error: Table '"+entry.getName()+"' in zip file\n"+inFile.getName()+"\nhas no valid data.");
                    return null;
                }
                return rt;
            } //if (entry.getName().endsWith(".csv"))
        }
        IJ.error("Error: cannot find required .csv file in "+inFile.getName());
        return null;
    }

    /** Returns a list of all file entries in top level of the zip archive.
     *  Directories in the archive and their contents are ignored.
     *  Displays an exception and returns null in case of error. */
    private ArrayList<ZipEntry> getZipEntries(File inFile) {
        ArrayList<ZipEntry> list = new ArrayList<ZipEntry>(1000);
        ZipFile zipFile = null;
        try {
            zipFile = new ZipFile(inFile, ZipFile.OPEN_READ);
            Enumeration<? extends ZipEntry> entries = zipFile.entries();
            while (entries.hasMoreElements()) {
                ZipEntry entry = entries.nextElement();
                if (!entry.isDirectory())
                    list.add(entry);
            }
        } catch (Exception e) {
            IJ.handleException(e);
            list = null;
        } finally {
            try {
                zipFile.close();
            } catch (Exception ee) {}
        }
        return list;
    }

    /** Finds a zip entry with a given name in the list */
    private ZipEntry findEntry(ArrayList<ZipEntry> entries, String name) {
        for (ZipEntry entry : entries)
            if (entry.getName().equals(name))
                return entry;
        return null;
    }

    /** Returns an ImagePlus with the stack. Uses instance variables nFrames, zipEntries,
     *  writes info of first image to instance variable imageInfo */
    private ImageStack makeStack(File inFile, ResultsTable rt, boolean useVirtualStack) {
        sliceLabels = makeSliceLabels(rt, nFrames);
        ImageStack stack = null;
        ZipFile zipFile = null;
        try {
            zipFile = new ZipFile(inFile, ZipFile.OPEN_READ);
            IJ.showStatus("Opening "+inFile.getName());
            ImagePlus imp0 = getImagePlus(zipFile, zipEntries[0]);
            Object info = imp0.getProperty("Info");
            if (info != null)
                imageInfo.append((String)info);
            if (useVirtualStack) {
                stack = this;
            } else {
                stack = new ImageStack(width, height);
                stack.addSlice(sliceLabels[0], imp0.getProcessor().getPixels());
                for (int i=1; i<nFrames; i++) {
                    IJ.showProgress(i, nFrames);
                    ImageProcessor ip = getProcessor(zipFile, zipEntries[i]);
                    stack.addSlice(sliceLabels[i], ip.getPixels());
                }
            }
        } catch (Exception e) {
            IJ.handleException(e);
            stack = null;
        } finally {
            try {
                IJ.showProgress(1.0);
                zipFile.close();
            } catch (Exception ee) {}
        }
        return stack;
    }

    /** For a virtual stack, reads the image slice from the file 'on the fly' */
    public ImageProcessor getProcessor(int n) {
        if (n < 1 || n > nFrames)
            throw new IllegalArgumentException("Slice "+n+" beyond stack limit ("+nFrames+")");
        ImageProcessor ip = null;
        ZipFile zipFile = null;
        try {
            zipFile = new ZipFile(inFile, ZipFile.OPEN_READ);
            ZipEntry zipEntry = zipEntries[n-1];
            if (zipEntry == null)
                throw new RuntimeException("Missing image in zip file");
            ip = getProcessor(zipFile, zipEntry);
        } catch (Exception e) {
            LeedUtils.logError("Error opening "+inFile.getName()+", slice="+n+"\n"+e);
            ip = new ShortProcessor(width, height);
        } finally {
            try {
                zipFile.close();
            } catch (Exception ee) {}
        }
        ip.setSliceNumber(n);
        return ip;
    }


	/** Returns the image width of the virtual stack.
     *  Must override the superclass since width, height of Stack are private */
	public int getWidth() {
		return width;
	}

	/** Returns the image height of the virtual stack */
	public int getHeight() {
		return height;
	}

    /** For a virtual stack, does nothing */
    public void addSlice(String fileName) {}

    /** For a virtual stack; can only delete the last slice */
    public void deleteSlice(int n) {
        if (n == nFrames)
            deleteLastSlice();
        else
            throw new RuntimeException("Only the last stack slice can be deleted");
    }

    /** For a virtual stack; deletes the last slice */
    public void deleteLastSlice() {
        if (nFrames > 1) nFrames --;
    }
    /** For a virtual stack: Number of slices */
    public int getSize() {
        return nFrames;
    }

    /** For a virtual stack; returns the label of the Nth image. */
    public String getSliceLabel(int n) {
        return sliceLabels[n-1];
    }

    /** For a virtual stack; there is no directory with the images. */
    public String getDirectory() {
        return null;
    }

    /** For a virtual stack; returns the file name of the specified slice, were 1<=n<=nslices. */
    public String getFileName(int n) {
        ZipEntry zipEntry = zipEntries[n-1];
        return zipEntry == null ? null : zipEntry.getName();
    }

    /** Returns the ImageProcessor for a tiff in a given entry.
     *  Note that zipFile must be open already */
    private ImageProcessor getProcessor(ZipFile zipFile, ZipEntry zipEntry) throws Exception {
        ImagePlus imp = getImagePlus(zipFile, zipEntry);
        return imp.getProcessor();
    }

    /** Returns the ImagePlus for the image described by in a given entry.
     *  Note that zipFile must be open already */
    private ImagePlus getImagePlus(ZipFile zipFile, ZipEntry zipEntry) throws Exception {
        InputStream inputStream = zipFile.getInputStream(zipEntry);
        TiffDecoder decoder = new TiffDecoder(inputStream, "");
        FileInfo[] fiArray = decoder.getTiffInfo();  //one-element array (we don't handle stacks in one file)
        inputStream.close();
        if (fiArray == null || fiArray .length <1)
            LeedUtils.logError("Error: can't decode tiff file "+zipEntry.getName());
        fiArray[0].inputStream = zipFile.getInputStream(zipEntry);
        FileOpener fo = new FileOpener(fiArray[0]);
        ImagePlus imp = fo.openImage();
        if (imp == null)
            throw new IOException("Cannot open as image: "+zipEntry.getName());
        return imp;
    }

    private String[] makeSliceLabels(ResultsTable rt, int nFrames) {
        String[] sliceLabels = new String[nFrames];
        StringBuilder label = new StringBuilder(100);
        String[] headings = rt.getHeadings();
        String[] shortHeadings = headings.clone();
        for (int i=0; i<headings.length; i++) {     //shorter names for stack slice labels (these will be truncated for display)
            if ("energy".equalsIgnoreCase(shortHeadings[i]))
                shortHeadings[i] = "E";
            else if ("time".equalsIgnoreCase(shortHeadings[i]) || "times".equalsIgnoreCase(shortHeadings[i]))
                shortHeadings[i] = "t";
            else if ("temperature".equalsIgnoreCase(shortHeadings[i]) || "temperatures".equalsIgnoreCase(shortHeadings[i]))
                shortHeadings[i] = "T";
            else if ("cold_junction".equalsIgnoreCase(shortHeadings[i]))
                shortHeadings[i] = "Tcj";
            else if ("I_sample".equalsIgnoreCase(shortHeadings[i]))
                shortHeadings[i] = "Is";
            else if ("label".equalsIgnoreCase(shortHeadings[i]) ||
                    "date".equalsIgnoreCase(shortHeadings[i]) ||
                    "clock".equalsIgnoreCase(shortHeadings[i]))
                shortHeadings[i] = null;  // don't write filename, day&time to slice labels
        }

        for (int f=0; f<nFrames; f++) {
            label.setLength(0);
            for (int i=0; i<headings.length; i++) {
                String shortName = shortHeadings[i];
                if (shortName == null) continue;                    // column to ignore
                double value = rt.getValue(headings[i], f);
                if (Double.isNaN(value)) continue;
                int digits = 5;
                if (shortName.equals("E"))
                    digits = 1;
                else if (shortName.equals("t"))
                    digits = 2;
                else {
                    if (value < 0.01) digits++;
                    if (value < 0.001) digits++;
                }
                label.append(shortName);
                label.append('=');
                label.append(IJ.d2s(value, digits));
                label.append(", ");
            }
            if (label.length() > 2) {
                label.delete(label.length()-2, label.length());     //remove last comma and space
                sliceLabels[f] = label.toString();
            }
        }
        return sliceLabels;
    }

    private void setCalibration(ImagePlus imp, ResultsTable rt) {
        String[] headings = rt.getHeadings();
        for (int i=0; i<headings.length; i++) {       //try to find a column for calibration
            if ("label".equalsIgnoreCase(headings[i])) continue;
            double[] data = rt.getColumn(headings[i]);
            double[] firstInc = LeedUtils.getFirstAndIncrement(data);
            if (!Double.isNaN(firstInc[0])) {                     //found data with linear increase
                Calibration cal = imp.getCalibration();
                cal.pixelDepth = firstInc[1];
                cal.zOrigin = -firstInc[0]/cal.pixelDepth;
                cal.setUnit("pxl"); //make it calibrated, even if z increment is 1
                if ("energy".equalsIgnoreCase(headings[i]) || "E".equals(headings[i]))
                    cal.setZUnit("eV");
                if ("time".equalsIgnoreCase(headings[i]) || "t".equals(headings[i]))
                    cal.setZUnit("sec");
                break;
            }
        }
    }
}
