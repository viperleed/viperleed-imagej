import ij.plugin.*;
import ij.*;
import ij.gui.*;
import ij.io.*;
import ij.measure.*;
import java.awt.*;
import java.io.*;
import java.util.*;
import java.text.SimpleDateFormat;


/**
 * ImageJ Plugin for reading 'AIDA' (Automatic Image and Data Acquisition)
 * EE2000/EE2010 LEED movies (*.vid format).
 * For the EE2000/EE2010 data acquisition system, see http://www.ee2000.de.
 *
 * This plugin can also use an AIDA LEED movie as a template and write an ImageJ
 * 16-bit stack currently open as foreground image foreground as an AIDA file,
 * using the metadata of the template, but replacing the image data.
 * In that case, the image size and number of images must be the same for both.
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


public class Open_Aida_LEED_Video implements PlugIn {
    static final String FILE_START = "EEVID";	// "magic" of file
    static final int MAX_STRING_LENGTH = 20;    // parameter designations cannot be longer than this
    static final int DATA_OFFSET = 1068;        // data start at this position
    static final SimpleDateFormat DATE_FORMAT = new SimpleDateFormat("yyyy-MM-dd", Locale.US);
    static final SimpleDateFormat TIME_FORMAT = new SimpleDateFormat("HH:mm:ss", Locale.US);
    static final long START_DATE = (new Date(-1,11,30)).getTime(); //day 0 in 'Date' field, 1899-12-30

    int width, height, nFrames, nParam;
    int nFramesReadable;                        // may be less than nFrames if file is truncated
    String[] paramNames;
    String dateTime;
    RandomAccessFile raFile;
    boolean enableSaveAsVid = false;
    //dialog parameters useVirtualStack, showTable, hideDialog: we use the variables of Open_LEED_Movie

  /** The plugin is invoked by ImageJ using this method.
   *  @param arg   String 'arg' may be used to select the path. If it is an empty string,
   *  a file open dialog is shown.
   */
    public void run (String arg) {
        enableSaveAsVid = true;
        if (IJ.altKeyDown()) Open_LEED_Movie.hideDialog = false;
        boolean useVirtualStack = LeedParams.getBoolean(LeedParams.OPENSTACKSASVIRTUAL);
        boolean showTable = Open_LEED_Movie.showTable;
        boolean hideDialog = Open_LEED_Movie.hideDialog;
        ImagePlus imp = openImagePlus(arg, useVirtualStack, showTable, hideDialog);
        if (imp == null) return;
        IJ.run(imp, "Enhance Contrast", "saturated=0.35");
        imp.show();
    }

  /** Opens a .vid file with a LEED movie, without displaying it,
   *  and returns a reference to the ImagePlus.
   *  @param path   Path to the file. If it is null or
   *        not a path to a readable file, a file open dialog is shown.
   *  @return The ImagePlus or null on error/cancel/no output (template)
   */
    public ImagePlus openImagePlus(String path, boolean useVirtualStack, boolean showTable, boolean hideDialog) {
        if (path.length() == 0) path = null;
        if (path != null) {
            File file = new File(path);
            if (!file.canRead())
                path = null;
        }
        if (path == null) {
            OpenDialog od = new OpenDialog("Open AIDA LEED Movie (.vid)", "");
            String fileName = od.getFileName();
            if (fileName == null) return null;
            String directory = od.getDirectory();
            path = directory + fileName;
        }
        ResultsTable rt;
        try {
            rt = readHeaderAndTable(path);      // read header & params
        } catch (Exception e) {
            IJ.error("Error Opening AIDA Video", e.getMessage());
            if (raFile != null) try {
                raFile.close();
            } catch (Exception e2) {}
            return null;
        }
        File file = new File(path);
        String fileName = file.getName();
        String directory = file.getParent();

        String firstDate = rt.getStringValue("Date", 0);
        String firstClock = rt.getStringValue("Clock", 0);
        if (firstDate != null && firstClock != null)
            dateTime = firstDate + ' ' + firstClock;

        ImagePlus img = WindowManager.getCurrentImage();
        boolean useAsTemplate = enableSaveAsVid && img != null &&
                img.getWidth()==width && img.getHeight()==height && img.getBitDepth()==16;

        if (!hideDialog) {                   // dialog
            GenericDialog gd = new GenericDialog("Open AIDA LEED Movie...");
            double sizeMB = nFrames*width*height*2/(1024.*1024.);
            gd.addMessage(fileName+"\n"+nFrames+" frames, "+width+"x"+height+" ("+IJ.d2s(sizeMB, 1)+" MB)");
            gd.addCheckbox("Virtual Stack", useVirtualStack);
            gd.addCheckbox("Show Parameter Table", showTable);
            gd.addCheckbox("Don't Show This Dialog Again *", hideDialog);
            if (useAsTemplate)
                gd.addCheckbox("Store "+img.getShortTitle()+" in Same Format **", false);
            gd.addMessage("* Invoke with ALT key to display the dialog again");
            if (useAsTemplate)
                gd.addMessage("** Does not open the input file chosen.\n"+
                        "   Rather creates a copy thereof with only the images\n"+
                        "   taken from the current foreground image stack\n"+
                        "   Use only if both movies have the same metadata (energy...)!");
            gd.showDialog();
            if (gd.wasCanceled()) return null;
            useVirtualStack = gd.getNextBoolean();
            showTable = gd.getNextBoolean();
            hideDialog = gd.getNextBoolean();
            if (useAsTemplate)
                useAsTemplate = gd.getNextBoolean();
            if (useAsTemplate) {
                hideDialog = false;
                showTable = false;
            }
            LeedParams.set(LeedParams.OPENSTACKSASVIRTUAL, useVirtualStack);
            Open_LEED_Movie.showTable = showTable;
            Open_LEED_Movie.hideDialog = hideDialog;

        }
        if (rt != null && showTable) {
            int dotPos = fileName.lastIndexOf('.');
            String tableName = dotPos > 0 ? fileName.substring(0, dotPos) : fileName+"_T";
            rt.show(tableName);
        }

        int gapBetweenImages = nParam*4;
        int offset = DATA_OFFSET + gapBetweenImages;
        if (useAsTemplate) {
            useTemplate(directory, fileName, img, Math.min(nFrames, img.getNSlices()), gapBetweenImages, offset);
            return null;
        }
        FileInfo fi = new FileInfo();                           // parameters for opening
        fi.fileName = fileName;
        fi.directory = directory;
        fi.fileType = FileInfo.GRAY16_UNSIGNED;
        fi.width = width;
        fi.height = height;
        fi.nImages = nFramesReadable;
        fi.gapBetweenImages = gapBetweenImages;
        fi.intelByteOrder = true;
        fi.sliceLabels = (rt == null) ? null : makeSliceLabels(rt);
        ImagePlus imp = null;
        if (useVirtualStack) {                                  // open
            FileInfo[] infos = new FileInfo[nFramesReadable];
            for (int f=0; f<nFramesReadable; f++) {
                infos[f] = (FileInfo)(fi.clone());
                infos[f].longOffset = offset + (long)f*(width*height*2+gapBetweenImages);
            }
            FileInfoVirtualStack fiv = new FileInfoVirtualStack(infos);
            imp = new ImagePlus(WindowManager.makeUniqueName(fileName), fiv);
        } else {
            fi.offset = offset;
            imp = (new FileOpener(fi)).open(false);
        }

        if (rt != null && imp.getNSlices() == nFramesReadable && nParam > 0) {
            for (int i=0; i<paramNames.length; i++) {       //try to find a column for calibration
                double[] data = rt.getColumn(paramNames[i]);
                double[] firstInc = LeedUtils.getFirstAndIncrement(data);
                if (!Double.isNaN(firstInc[0])) {                     //found data with linear increase
                    Calibration cal = imp.getCalibration();
                    cal.pixelDepth = firstInc[1];
                    cal.zOrigin = -firstInc[0]/cal.pixelDepth;
                    cal.setUnit("pxl"); //make it calibrated, even if z increment is 1
                    if ("Energy".equals(paramNames[i])) cal.setZUnit("eV");
                    if ("Time".equals(paramNames[i])) cal.setZUnit("sec");
                    break;
                }
            }
        }
        if (dateTime != null)
            imp.setProperty("Info", "date: "+dateTime+'\n');

        return imp;
    }

    /** Open the file with given path and read its header as well as the ResultsTable */
    private ResultsTable readHeaderAndTable(String path) throws IOException {
        long lastTime = System.currentTimeMillis();
        File file = new File(path);                             // o p e n
        raFile = new RandomAccessFile(file, "r");
        long fileSize = raFile.length();
        if (fileSize < DATA_OFFSET)
            throw new RuntimeException("Error: File size only "+fileSize+" bytes");
        byte[] byteBuf = new byte[FILE_START.length()];         // m a g i c   c o d e
        raFile.read(byteBuf);
        String fileStart = new String(byteBuf);
        if (!FILE_START.equals(fileStart))
            throw new RuntimeException("Error: This does not look like an AIDA LEED file, starts with '"+fileStart+"', not '"+FILE_START+"'.");
        raFile.seek(0x08);                                      // f i l e     h e a d e r
        nFrames = readShort(raFile);
        raFile.seek(0x10);
        width = readShort(raFile);
        raFile.seek(0x14);
        height = readShort(raFile);
        raFile.seek(0x1c);
        nParam = readShort(raFile);
        raFile.seek(0x20);
        paramNames = readZeroTerminatedStrings(raFile, nParam);

        int bytesPerFrame = 4*nParam + 2*width*height;          // c h e c k   f i l e   s i z e
        long expectedSize = DATA_OFFSET + nFrames*bytesPerFrame;
        nFramesReadable = nFrames;
        int badFrames = (int)((expectedSize - fileSize + bytesPerFrame - 1) / bytesPerFrame);
        if (badFrames > 0)
            nFramesReadable -= badFrames;
        if (nFramesReadable < 1)
            throw new RuntimeException("Error: File size only "+fileSize+" bytes, needed for 1 frame: "+(DATA_OFFSET+bytesPerFrame));
        else if (badFrames > 0)
            IJ.error("Warning", "Aida Video file "+(new File(path)).getName()+
                    "\nis incomplete, "+badFrames+" frames missing at the end");
        ResultsTable rt = new ResultsTable(nFramesReadable);    // r e a d   p a r a m s
        Date date = new Date();
        for (int f=0; f<nFramesReadable; f++) {
            raFile.seek(DATA_OFFSET + f*bytesPerFrame);
            for (int i=0; i<nParam; i++) {
                float value = readFloat(raFile);
                if ("Date".equals(paramNames[i]))
                    rt.setValue(paramNames[i], f, getDate(value, date));
                else if ("Clock".equals(paramNames[i]))
                    rt.setValue(paramNames[i], f, getClock(value, date));
                else
                    rt.setValue(paramNames[i], f, value);
            }
            if (System.currentTimeMillis() - lastTime > 100) {
                IJ.showProgress(f, nFramesReadable);
                lastTime = System.currentTimeMillis();
            }
        }
        raFile.close();
        IJ.showProgress(1.0);
        return rt;
    }

    private String[] makeSliceLabels(ResultsTable rt) {
        String[] sliceLabels = new String[nFramesReadable];
        StringBuilder label = new StringBuilder(100);
        for (int f=0; f<nFramesReadable; f++) {
            label.setLength(0);
            for (int i=0; i<nParam; i++) {
                String name = paramNames[i];
                if ("Date".equals(name) || "Clock".equals(name))
                        continue; //don't write day&time
                String shortName = name;
                int digits = 5;
                if ("Energy".equals(name)) {
                    shortName = "E";
                    digits = 1;
                } else if ("Time".equals(name)) {
                    shortName = "t";
                    digits = 2;
                }
                label.append(shortName);
                label.append('=');
                label.append(IJ.d2s(rt.getValue(name, f), digits));
                label.append(',');
            }
            if (label.length() > 1) {
                label.deleteCharAt(label.length()-1);   //remove last comma
                sliceLabels[f] = label.toString();
            }
        }
        return sliceLabels;
    }

    /** creates a new file from the current one, replacing only the image data with those of img */
    void useTemplate(String directory, String fileName, ImagePlus img, int nFrames, int gapBetweenImages, int offset) {
        boolean nameOk = false;
        String outName = null, outDirectory = null;
        while (!nameOk) {
            SaveDialog sd = new SaveDialog("Save as AIDA-PC video", img.getTitle(), ".vid");
            outName = sd.getFileName();
            if (outName==null)
                return;
            outDirectory = sd.getDirectory();
            nameOk = !(new File(directory+fileName)).equals(new File(outDirectory+outName));
            if (!nameOk) IJ.error("Output file and input file must be different");
        }

        RandomAccessFile infile = null, outfile = null;
        try {
            infile = new RandomAccessFile(new File(directory+fileName), "r");
            outfile = new RandomAccessFile(new File(outDirectory+outName), "rw");
            byte[] buffer = new byte[offset];
            infile.read(buffer, 0, offset);
            outfile.write(buffer, 0, offset);
            for (int f=0; f<nFrames; f++) {
                IJ.showProgress(f, nFrames);
                short[] pixels = (short[])(img.getStack().getProcessor(f+1).getPixels());
                writeShortArray(outfile, pixels);
                if (f < nFrames - 1) {
                    infile.seek(infile.getFilePointer() + 2*pixels.length);
                    infile.read(buffer, 0, gapBetweenImages);
                    outfile.write(buffer, 0, gapBetweenImages);
                }
            }

        } catch (Exception e) {
            IJ.handleException(e);
        } finally {
            try{ infile.close(); } catch (Exception ex) {}
            try{ outfile.close(); } catch (Exception ex) {}
            IJ.showProgress(1.0);
        }
    }


    /** Read 2-byte signed short with Intel (little-endian) byte order
     * (note: RandomAccessFile.readShort has big-endian byte order) */
    static final short readShort(RandomAccessFile raFile) throws IOException {
        int in = raFile.readShort();
        int b0 = in&0xff;
        int b1 = (in&0xff00)>>8;
        return (short) (b0 << 8 | b1);
    }

    /** Read 4-byte float with Intel (little-endian) byte order
     * (note: RandomAccessFile.readInt has big-endian byte order) */
    static final float readFloat(RandomAccessFile raFile) throws IOException {
        int in = raFile.readInt();
        int out = 0;
        for (int i=0; i<3; i++) {
            out |= in&0xff;
            in = in>>>8;
            out = out<<8;
        }
        out |= in;
        return Float.intBitsToFloat(out);
    }

    /** Reads nStrings 0-terminated Strings */
    static final String[] readZeroTerminatedStrings(RandomAccessFile raFile, int nStrings) throws IOException {
        String[] strings = new String[nStrings];
        byte[] byteBuf = new byte[MAX_STRING_LENGTH];
        for (int n=0; n<nStrings; n++) {
            for (int i=0; i<MAX_STRING_LENGTH; i++) {
                byte b = raFile.readByte();
                if (b == 0) {
                    strings[n] = new String(byteBuf, 0, i);
                    break;
                }
                byteBuf[i] = b;
            }
        }
        return strings;
    }

    /** Write array of 2-byte signed short with Intel (little-endian) byte order
     * (note: RandomAccessFile.writeShort has big-endian byte order) */
    static final void writeShortArray(RandomAccessFile raFile, short[] array) throws IOException {
        byte[] buffer = new byte[2*array.length];  //writing bytes directly into RandomAccessFile is too slow
        for (int i=0, p=0; i<array.length;i++) {
            short v = array[i];
            buffer[p++] = (byte)(v&0xff);
            buffer[p++] = (byte)(v>>8);
        }
        raFile.write(buffer);
    }


    /* converts days since 1999-12-30 to date string. Uses 'date' as workspace */
    static final String getDate(float value, Date date) {
        date.setTime(START_DATE + 86400L*1000L*(long)value);
        return DATE_FORMAT.format(date);
    }
    /* converts fractional day to time string. Uses 'date' as workspace */
    static final String getClock(float value, Date date) {
        date.setTime((long)(value*86400*1000));
        return TIME_FORMAT.format(date);
    }
}
