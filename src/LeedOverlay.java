import ij.*;
import ij.IJ;
import ij.gui.*;
import java.awt.*;
import java.util.*;
/**
 *  This class creates the overlays on the I(V) image stack:
 *  Energy label, mask outline, as well as circles and labels for spots.
 *  For simplicity, this class has static methods only;
 *  it can be used with only one image stack at a time.
 *  (This is not a problem since the LEED_SPot_tracker is a singleton)
 *
 *  Note: to display the overlay, use ImagePlus.draw, NOT ImagePlus.updataAndDraw!
 *  (the latter triggers the ImageListener and may cause deadlocks or endless loops)
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


public class LeedOverlay {
    static final String E_ETC = "Eetc";    //roi name for energy, mask, int&bg indicator (circles & beam labels have spot index as name)
    static final Color CIRCLE_COLOR = new Color(0x00c0ff);              //color for spot circles if there is no selection
    static final Color CIRCLE_WEAK_COLOR = new Color(0x0060c0);         //color for deselected spot circles
    static final Color LABEL_COLOR = Color.GREEN;                       //color for spot labels
    static final Color INT_INDICATOR_COLOR = new Color(0x00a0a0);       //color for integration and background area indicator
    static Font font = new Font("SansSerif", Font.PLAIN, 10);           //Font for beam labels. Size will vary
    static Font boldFont = new Font("SansSerif", Font.BOLD, 10);        //Bold font for beam labels. Size will vary
    static HashSet<String> spotNamesHighlighted;
    static int highlightedStrength;                         //sync from CurveEditor has strength 0, manual 1, bad beams 2

    /** Adds an overlay visible in all slices; the Roi gets name 'name'.
     *  All other Rois with this name are removed. Does not update the display. */
    public static void add(ImagePlus imp, Roi roi, String name) {
        if (imp == null || roi == null) return;
        Overlay ovly = getOverlay(imp);
        removeOverlay(ovly, name);
        ovly.add(roi, name);
    }

    /** Adds the mask outline, energy (or xAxis) labels,
     *  and integration area indicator in the corner.
     *  Note that xAxis[0] corresponds to the first stack slice 1.
     *  For regular LEED I(V), isEnergy should be true; then the integration
     *  area depends on the energy.
     *  Does not update the display; use imp.draw() to show the overlay. */
    public static void add(ImagePlus imp, double[] xAxis, boolean isEnergy, Roi maskRoi) {
        if (imp == null) return;
        Overlay ovly = getOverlay(imp);
        removeOverlay(ovly, E_ETC);
        if (maskRoi != null)
            add(imp, maskRoi, E_ETC);

        int width = imp.getWidth();
        int height = imp.getHeight();
        int fontSize = 10+height/100;               //for drawing energy labels
        final Font font = new Font("Arial", Font.PLAIN, fontSize);
        int digits = LeedUtils.getEnergyDigits(xAxis);
        int xIntIndicator = (int)Math.round(width - 0.08*Math.min(width, height)); //where to draw integration&bg indicator
        int yIntIndicator = (int)Math.round(0.08*Math.min(width, height));
        double dx = 0.25*Math.min(width, height);   //for drawing integration&bg indicator, we assume spot halfway from the center
        double dy = -dx;
        int backgroundType = (int)LeedParams.get(LeedParams.BACKGROUNDTYPE);
        double azBlurRadians = Math.toRadians(LeedParams.get(LeedParams.AZIMUTHBLURANGLE));
        
        for (int i=0; i<imp.getNSlices(); i++) {
            if (Thread.currentThread().isInterrupted()) return;
            if (xAxis != null) {
                Roi roi = new TextRoi(fontSize/2, fontSize/2, IJ.d2s(xAxis[i], digits), font);
                roi.setStrokeColor(Color.YELLOW);
                roi.setPosition(i+1);
                ovly.add(roi, E_ETC);
            }
            double energy = isEnergy && xAxis != null ?
                    xAxis[i] : LeedRadiusSelector.UNKNOWN_ENERGY;
            Roi[] rois = LeedSpotAnalyzer.getShapeRois(xIntIndicator, yIntIndicator, dx, dy,
                    backgroundType, LeedRadiusSelector.radius(energy, /*isSuperstructure=*/false), azBlurRadians);
            rois[0].setStrokeColor(INT_INDICATOR_COLOR);
            rois[1].setStrokeColor(INT_INDICATOR_COLOR);
            rois[0].setPosition(i+1);
            rois[1].setPosition(i+1);
            ovly.add(rois[0], E_ETC);
            ovly.add(rois[1], E_ETC);
        }
    }

    /** Sets the font size for the spot names */
    public static void setFontSize(float size) {
        if (size < 2f) size = 2f;
        if (size > 14f) size = 14f;
        if (Math.abs(size - font.getSize2D())/size > 0.05) { //ignore changes by less than 5%
            font = font.deriveFont(Font.PLAIN, size);
            boldFont = font.deriveFont(Font.BOLD, size);
        }
    }

    /** Adds a circle for the spot and a name (if not null).
     *  Does not update the display; use imp.draw() to show the overlay. */
    public static void add(ImagePlus imp, int slice, double x, double y, double radius, String name) {
        Overlay ovly = getOverlay(imp);
        //int diam = (int)Math.round(2*radius);
        //int x0 = (int)Math.round(x-0.5*diam);
        //int y0 = (int)Math.round(y-0.5*diam);
        Roi roi = new OvalRoi(x-radius+0.5, y-radius+0.5, 2*radius, 2*radius); //circle around pixel (0,0) has to left x=y=0
        roi.setPosition(slice);
        roi.setStrokeColor(CIRCLE_COLOR);
        ovly.add(roi, name);
        if (name != null)
            addName(imp, slice, x, y, radius, name, /*bold=*/false);
    }

    /** Adds a spot label; must be called when the overlay is present already.
     *  To save space, fractional spots (which tend to have long labels) use
     *  to lines. Does not update the display; use imp.draw() to show the overlay. */
    public static void addName(ImagePlus imp, int slice, double x, double y, double radius, String name, boolean bold) {
        Overlay ovly = imp.getOverlay();
        boolean twoLine = name.indexOf('\n') > 0;       //fractional spots are displayed with two lines
        int y0 = (int)(y + radius - (twoLine ? font.getSize2D() : 0.2f*font.getSize2D()));
        Roi roi = new TextRoi(name, (int)(x-radius), y0, bold ? boldFont : font);
        roi.setPosition(slice);
        roi.setStrokeColor(LABEL_COLOR);
        ovly.add(roi, name);
    }

    /** Deletes all spot circle overlays. Does not update the display.
     *  This will be too slow for removing overlays from the whole stack;
     *  use it only while only one or a few slices have overlays.
     *  Does not update the display; use imp.draw() for this. */
    public static void removeSpotOverlays(ImagePlus imp) {
        Overlay ovly = getOverlay(imp);
        removeOverlay(ovly, OvalRoi.class);
    }

    /** Deletes all label overlays. Does not update the display.
     *  This will be too slow for removing overlays from the whole stack;
     *  use it only while only one or a few slices have overlays.
     *  Does not update the display; use imp.draw() for this. */
    public static void removeNameOverlays(ImagePlus imp) {
        Overlay ovly = getOverlay(imp);
        removeOverlay(ovly, TextRoi.class);
    }

    /** Shows the spots in the hashset brighter (all others dim); or reverts if nameHashSet=null.
     *  The spot names must be given exactly in the same form as they were given
     *  in the 'add' method (two-line for fractional spots).
     *  This method does nothing if 'strength' is less than the current strength minus one.
     *  Since spots entered as 'bad' by the IVAnalyzer have strength = 2 and spots from
     *  synchronization with the IV Curve Editor have strength = 0, the Curve Editor does not
     *  erase information on bad spots.
     *  Call with spotNamesToHighlight = null to erase highlighting.
     *  Does not update the display; use imp.draw() thereafter. */
    public static void highlightCircles(ImagePlus imp, HashSet<String> spotNamesToHighlight, int strength) {
        if (strength < highlightedStrength - 1) return;
        Overlay ovly = imp.getOverlay();
        if (ovly == null) return;
        for (Iterator<Roi> iterator = ovly.iterator(); iterator.hasNext();) {
            Roi roi = iterator.next();
            if (roi instanceof OvalRoi) {
                if (spotNamesToHighlight != null && spotNamesToHighlight.contains(roi.getName())) {
                        roi.setStrokeColor(CIRCLE_COLOR);
                        roi.setStrokeWidth(roi.getFloatWidth() > 4 ? 2f : 1f);
                } else {  
                    roi.setStrokeWidth(0f);             // unhighlight previously highlighted
                    roi.setStrokeColor(spotNamesToHighlight == null ? CIRCLE_COLOR : CIRCLE_WEAK_COLOR);
                }
            }
        }
        spotNamesHighlighted = spotNamesToHighlight;    //remember the highlighted ones
        highlightedStrength = strength;
    }

    /** Deletes the spots in the overlay with the given names.
     *  Does not update the display; use imp.draw() thereafter.
     *  Synchronized to avoid two concurrent 'remove' operations */
    public static synchronized void delete(ImagePlus imp, HashSet<String> spotNamesToDelete) {
        Overlay ovly = imp.getOverlay();
        if (ovly == null) return;
        for (int i=ovly.size()-1; i>=0; i--) {
            Roi roi = ovly.get(i);
            if (spotNamesToDelete.contains(roi.getName()))
                ovly.remove(i);
        }
    }


    /** Returns a HashSet of highlighted spots, or null if none.
     *  Note that fractional spots will have '\n' as delimiter between h and k */
    static HashSet<String> getHighlighted(ImagePlus imp) {
        if (imp == null || imp.getOverlay() == null) return null;
        return spotNamesHighlighted;
    }

    /** Returns the overlay of the image; creates an overlay if the image has none. */
    static Overlay getOverlay(ImagePlus imp) {
        Overlay ovly = imp.getOverlay();
        if (ovly == null) {
            ovly = new Overlay();
            imp.setOverlay(ovly);
            spotNamesHighlighted = null;
        }
        return ovly;
    }

    /** Remove all overlay Rois with a given name.
     *  Synchronized to avoid two concurrent 'remove' operations */
    static synchronized void removeOverlay(Overlay ovly, String name) {
        ovly.remove(name);
    }

    /** Remove overlay Rois for beams with a given class.
     *  Energy labels, mask and integration indicators are not affected
     *  Synchronized to avoid two concurrent 'remove' operations */
    static synchronized void removeOverlay(Overlay ovly, Class roiClass) {
        spotNamesHighlighted = null;
        if (ovly == null) return;
        for (int i=ovly.size()-1; i>=0; i--) {
            Roi roi = ovly.get(i);
            if ((roi.getClass().equals(roiClass) && (roiClass == OvalRoi.class || roi.getName() != E_ETC)))
                ovly.remove(i);
        }
    }

}
