import ij.*;
import ij.gui.*;
import ij.util.*;
import ij.measure.Calibration;
import ij.process.*;
import ij.plugin.PlugIn;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.io.File;

/** Creates a simulated LEED image at a single energy or LEED I(V) movie
 *  based on a Spot Patter file and a file with I(V) curves.
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
 *  @author Michael Schmid, IAP/TU Wien, 2019-2025
 */


public class LEED_Movie_Simulation extends VirtualStack implements PlugIn, DialogListener {
    static final String PLUGIN_NAME = "LEED Movie Simulation";
    static final double MAX_FWHM = 20;
    static final int MIN_SIZE = 100, MAX_SIZE = 4000;
    static final double MIN_GAMMA = 1, MAX_GAMMA = 4;
    /* The following should be remembered until ImageJ is closed */
    static String ivFilePathS;
    static String spotPatternPathS;
    static double startES = 0;
    static double endES = 999;
    static double fwhmS = 5;
    static double lorentzFractionS = 0.1;
    static int sizeS = 800;
    static double gammaS = 2;
    static boolean leemModeS;
    /* parameters of the simulation */
    String ivFilePath;
    String spotPatternPath;
    double startE = 0;
    double endE = 999;
    double fwhm;
    double lorentzFraction;
    double fwhmSaved;				//keep old value during preview
    double lorentzFractionSaved;
    int size;
    double gamma;
    double gammaSaved;
    boolean leemMode;
    /* data */
    LeedSpotPattern spotPattern;
    LeedIVData ivData;
    int nSlices;
    int startEIndex;
    double maxRadius;				//maximum radius form (0,0) of spot with nonzero intensity
    boolean addEquivalent;			//if symmetry-equivalent spots are not given, we add them
    double maxIntensity;
    /* GUI-related */
    ArrayList<Component> createOnlyComponents = new ArrayList<Component>(); //components active only when creating an image/stack
    Label errorField;
    Checkbox previewCheckbox;
    int currentSliceWhenCalled;		//for detecting animation
    
    /* For modifying an existing simulation */
    ImagePlus targetImp;
    LEED_Movie_Simulation targetStack;

    private static final String HELP_STRING =
            "<html>"+
            "<h1>"+PLUGIN_NAME+"</h1>"+
            "<p>This plugin creates a simulated a LEED image or LEED I(V) movie</p>"+
            "<p>When called with an existing simulated a LEED image or movie at the currently selected "+
            "(=foreground) image in ImageJ, one can also modify the appearance of the spots in that "+
            "existing image or movie.</p>"+
            "<ul>"+
            "<li><b>I(V) File</b>: "+
            "Determines the intensities of each beam. This can be a theobeams.csv or expbeams.csv "+
            "file of viperleed.calc or a (raw or edited) IVcurves.csv file created by the Spot Tracker."+
            "</li>"+
            "<li><b>Spot Pattern File</b>: "+
            "A spot pattern file is a comma-separated list (.csv file) of beam names, <i>h</i>, <i>k</i>, <i>gx</i>, <i>gy</i> "+
            "(reciprocal lattice vector in Cartesian coordinates), and beam group number (equal for symmetry-equivalent beams). "+
            "Such a file can be created with the pattern simulator of the ViPErLEED GUI."+
            "</li>"+
            "<li><b>Energy Start, End</b>: "+
            "The energy range of the movie (inclusive). If this range is wider than the energy range of the I(V) file, "+
            "the I(V) file restricts the energy range. Therefore, you can enter e.g., 0 and 999 for the full "+
            "energy range if the I(V) file. The energy step is determined by the I(V) file."+
            "</li>"+
            "<li><b>LEEM Mode</b>: In LEEM mode, the spots do not move with energy but remain stationary."+
            "</li>"+
            "<li>"+
            "<b>Image Size</b>: Size of the resulting image or movie (square) in pixels."+
            "</li>"+
            "<li><b>FWHM and Lorentzian Fraction</b>: "+
            "Determines the full width at half maximum of the spots. The spots are Gaussians if Lorentzian fraction is 0, "+
            "and pure Lorentzians (with a soft cutoff) if the Lorentzian fraction is 1. For intermediate fractions, "+
            "a Gaussian+Lorentzian sum is used."+
            "</li>"+
            "<li><b>Display Gamma</b>: When >1, reduces the contrast by taking the gamma-th root of the intensity. "+
            "Most computer screens use a display gamma value near 2; thus this value makes the image appear with similar contrast "+
            "as on the actual LEED screen."+
            "</li>"+
            "<li>If the <b>Modify Current Image/Movie Only</b> checkbox is enabled, the FWHM, Lorentzian fraction, and display gamma "+
            "will be modified in the already existing simulated image or movie that is currently in the foreground."+
            "</li>"+
            "</ul>"+
            "<p>The output is a vitual stack, i.e., the images are calculated only if required (on the fly). "+
            "If you want to modify it by any image processing operations, you have to "+
            "duplicate it (Image&gt;Duplicate...).</p>"+
             "<p>The calculation may be slow for large images with many spots and large FWHM. "+
             "Calculations are faster with a Lorentzian fraction of 0.</p>"+
            "<h2><a name='ivEditorLicense'>License</a></h2>"+
            "<p>The code is licensed under <a href='http://www.gnu.org/licenses/gpl-3.0.html'>GNU General Public License v3.0</a> "+
            "or later (GPL-3.0-or-later).</p>"+
            "<p>&nbsp;&nbsp;&nbsp;&nbsp;The ViPErLEED ImageJ plugin collection is free software: you can redistribute it and/or modify it "+
            "under the terms of the GNU General Public License as published by the Free Software Foundation, "+
            "either version 3 of the License, or (at your option) any later version.<br>&nbsp;&nbsp;&nbsp;&nbsp;"+
            "The ViPErLEED ImageJ plugin collection is distributed in the hope that it will be useful, "+
            "but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. "+
            "See the GNU General Public License for more details.<br>&nbsp;&nbsp;&nbsp;&nbsp;"+
            "You should have received a copy of the GNU General Public License along with these plugins. "+
            "If not, see <a href='https://www.gnu.org/licenses/'>https://www.gnu.org/licenses/</a>.</p>"+
            "<p>The authors may decide later to put part of the auxiliary code in this work into the public domain, "+
            "to allow incorporation into ImageJ if desired (ImageJ is in the public domain).</p>"+
            "<p>The documentation, including the help texts, is licensed under the "+
            "<a href='http://creativecommons.org/licenses/by/4.0/'>Creative Commons Attribution 4.0</a> "+
            "(CC BY 4.0) license.</p>"+
            "<p>When using this program (in its original or modified form) for scientific work, "+
            "please cite the paper describing the program [<a href='#paper'>1</a>].</p>"+
            "<p>A copy of these license texts and the source code is included in the jar archive holding this plugin "+
            "(use an unzip utility to view its contents).</p>"+
            "<h2>References</h2>"+
            "<p><a name='paper'>[1]</a> M. Schmid, F. Kraushofer, A. M. Imre, T. Ki&szlig;linger, L. Hammer, U. Diebold, and M. Riva<br>"+
            "<i>ViPErLEED package II: Spot tracking, extraction and processing of <i>I</i>(<i>V</i>) curves</i><br>"+
            "<a href='https://doi.org/10.1103/PhysRevResearch.7.013006'>Phys. Rev. Research <b>7</b>, 013006 (2025)</a>.</p>"+
            "</html>";

    /** This method is called first
     * @param arg unused 
     */
    public void run(String arg) {
		if (IJ.versionLessThan("1.52u")) return;
		ImagePlus imp = WindowManager.getCurrentImage();
        if (imp != null && (imp.getStack() instanceof LEED_Movie_Simulation)) {
			targetImp = imp;
			currentSliceWhenCalled = targetImp.getSlice();
			targetStack = (LEED_Movie_Simulation)imp.getStack();
			targetStack.saveFwhmEtc();
        }
		// D I A L O G
        spotPatternPath = spotPatternPathS == null ?
            spotPatternPath = LeedParams.getString(LeedParams.PATTERNFILE) : spotPatternPathS;
        GenericDialog gd = new GenericDialog(PLUGIN_NAME);
        gd.addFileField("I(V) Curve File", ivFilePathS, 50);
        createOnlyComponents.add((Component)(gd.getStringFields().lastElement()));
        createOnlyComponents.add(gd.getLabel());
        gd.addFileField("Spot Pattern File", spotPatternPath, 50);
        createOnlyComponents.add((Component)(gd.getStringFields().lastElement()));
        createOnlyComponents.add(gd.getLabel());
        gd.addNumericField("Energy Start", startES, 1);
        createOnlyComponents.add((Component)(gd.getNumericFields().lastElement()));
        createOnlyComponents.add(gd.getLabel());
        gd.addToSameRow();
        gd.addNumericField("End", endES, 1, 6, "(leave empty for single image)");
        createOnlyComponents.add((Component)(gd.getNumericFields().lastElement()));
        createOnlyComponents.add(gd.getLabel());
        gd.addCheckbox("LEEM Mode", leemModeS);
        createOnlyComponents.add((Component)(gd.getCheckboxes().lastElement()));
        gd.addNumericField("Image Size", sizeS, 0, 6, "("+MIN_SIZE+"-"+MAX_SIZE+")");
        createOnlyComponents.add((Component)(gd.getNumericFields().lastElement()));
        createOnlyComponents.add(gd.getLabel());
        gd.addNumericField("FWHM", fwhmS, 1, 6, "(max. "+(int)MAX_FWHM+")");
        gd.addNumericField("Lorentzian Fraction", lorentzFractionS, 2, 6, "(0-1)");
        gd.addNumericField("Display Gamma", gammaS, 0, 6, "("+(int)MIN_GAMMA+"-"+(int)MAX_GAMMA+")");
        if (targetStack != null) {
			String imgOrMovie = targetStack.size() == 1 ? "Image" : "Movie";
			gd.addCheckbox("Modify Current "+imgOrMovie+" Only", false);
			previewCheckbox = (Checkbox)(gd.getCheckboxes().lastElement());
		}
        gd.addMessage("You can drag and drop (.csv-)files into the 'file' fields.");
        gd.addMessage("");
        errorField = (Label)gd.getMessage();

        gd.addDialogListener(this);
        Button okButton = gd.getButtons()[0];
        okButton.setEnabled(dialogItemChanged(gd, null));
        gd.addHelp(HELP_STRING);
        gd.showDialog();
        boolean previewRunning = previewCheckbox != null && previewCheckbox.getState();
        if (gd.wasCanceled()) {
			if (previewRunning)
				revertPreview();
			return;
		}
        if (previewRunning)	//target has been modified already during preview
			return;


        ivFilePath = gd.getNextString();
        spotPatternPath = gd.getNextString();
        startE = gd.getNextNumber();
        endE = gd.getNextNumber();
        leemMode = gd.getNextBoolean();
        size = (int)gd.getNextNumber();
        fwhm = gd.getNextNumber();
        lorentzFraction = gd.getNextNumber();
        gamma = gd.getNextNumber();

		if (!(startE >= 0)) {
			startE = 0;
		} else {
			if (Double.isNaN(endE) || endE < startE)
				endE = startE;
			startES = startE;
			endES = endE;
		}
		if (fwhm >= 0 && fwhm <= MAX_FWHM)
			fwhmS = fwhm;
		else
			fwhm = 4;
		if (lorentzFraction >= 0 && lorentzFraction <= 1)
			lorentzFractionS = lorentzFraction;
		else
			lorentzFraction = 0.1;
		if (!(gamma >= MIN_GAMMA))
			gamma = MIN_GAMMA;
		else if (gamma > MAX_GAMMA)
			gamma = MAX_GAMMA;
		else
			gammaS = gamma;
		if (!(size >= MIN_SIZE))
			size = MIN_SIZE;
		else if (size > MAX_SIZE)
			size = MAX_SIZE;
		else
			sizeS = size;

        // R E A D   D A T A
        spotPattern = new LeedSpotPattern(spotPatternPath);
        String err = spotPattern.getErrorMessage();
        if (err != null) {
			IJ.error(err);
			return;
		}
        ivData = LeedIVData.fromFile(ivFilePath, spotPattern, LeedIVData.E_ASCENDING);
        if (ivData == null) return;
        
        ivFilePathS = ivFilePath;
        spotPatternPathS = spotPatternPath;

        // A N A L Y Z E   D A T A
        // find the first index in ivData and stack size
        double[] energies = ivData.energies;
        double[] minMaxE = Tools.getMinMax(energies);
        if (startE < minMaxE[0]) startE = minMaxE[0];
        if (endE < minMaxE[0]) endE = minMaxE[0];
        if (startE > minMaxE[1]) startE = minMaxE[1];
        if (endE > minMaxE[1]) endE = minMaxE[1];
        startEIndex = -1;
        double lastE = 2*energies[0] - energies[1];
        for (int i=0; i<energies.length; i++) {
			double step = energies[i] - lastE;
			if (startEIndex < 0 && energies[i] + 0.3*step > startE)
				startEIndex = i;
			if (i > startEIndex && energies[i] - 0.3*step < endE)
				nSlices = i - startEIndex + 1;
			lastE = energies[i];
		}
		if (nSlices < 1)
			nSlices = 1;

		maxRadius = getMaxRadius();
		addEquivalent = shouldAddEquivalent();

        // C R E A T E   S T A C K
        String ivName = (new File(ivFilePath)).getName();
        String title = (ivName+" Simulation");
        if (WindowManager.getImage(title) != null)
			title = WindowManager.getUniqueName(title);
        ImagePlus newImp = getImagePlus(title);
	double[] firstAndStep = LeedUtils.getFirstAndIncrement(energies);
	if (firstAndStep[1] > 0) {
	    Calibration cal = newImp.getCalibration();
	    cal.pixelDepth = firstAndStep[1]; //energy step
	    cal.zOrigin = -startE/firstAndStep[1];
	}
	newImp.show();
        return;
    }

    /** Creates an ImagePlus showing this stack */
    ImagePlus getImagePlus(String title) {
        ImagePlus outImp = new ImagePlus(title, this);
        makeOverlay(outImp);
        return outImp;
    }
    /** Creates an overly with the energy in the top-left corner */
    void makeOverlay(ImagePlus imp) {
		Overlay ovly = new Overlay();
		int fontSize = (int)Math.sqrt(0.5*size);
		Font font = new Font("Arial", Font.PLAIN, fontSize);
		for (int n=1; n<=nSlices; n++) {
			Roi roi = new TextRoi(fontSize/2, fontSize/2, IJ.d2s(getEnergy(n), 1), font);
			roi.setStrokeColor(Color.YELLOW);
            roi.setPosition(n);
            ovly.add(roi, "energy");
		}
		imp.setOverlay(ovly);
	}

    /** Checks the dialog input for validity and returns false if invalid.
     *  Provides live update for the disply parameters of a target stack.
     *  Also enables/disables components that are not active for a target stack */
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        boolean previewRunning = previewCheckbox != null && previewCheckbox.getState();
        if (previewRunning && targetImp.getSlice() != currentSliceWhenCalled) {
			previewCheckbox.setState(false);
			return error("Cannot modify a running movie");
		}
        String ivFilePath = gd.getNextString();
        String spotPatternPath = gd.getNextString();
        double startE = gd.getNextNumber();
        double endE = gd.getNextNumber();
        boolean leemMode = gd.getNextBoolean();
        double size = (int)gd.getNextNumber();
        double fwhm = gd.getNextNumber();
        double lorentzFraction = gd.getNextNumber();
        double gamma = (int)gd.getNextNumber();

		if (!(fwhm >= 0 && fwhm <= MAX_FWHM)) return error("invalid FWHM");
		if (!(lorentzFraction >= 0 && lorentzFraction <= 1)) return error("invalid Lorentz fraction");
		if (!(gamma >= MIN_GAMMA && gamma <= MAX_GAMMA)) return error("invalid display gamma");

        if (previewRunning) {		//modify an existing stack with preview
			targetStack.setFwhmEtc(fwhm, lorentzFraction, gamma);
			updateImage(targetImp, targetStack);
		} else if (e != null && e.getSource() == previewCheckbox)
			revertPreview();
		enableCreateOnlyFields(!previewRunning);

        if (!LeedUtils.fileOk(ivFilePath)) return error("cannot open I(V) file");
        LeedSpotPattern spotPattern = new LeedSpotPattern(spotPatternPath);
        String err = spotPattern.getErrorMessage();
        if (err != null) return error(err);
        if (!(startE >= 0)) return error("invalid start energy");
        if (endE > 0 && endE < startE) return error("end energy must be at least start energy");
        if (!(size >= MIN_SIZE && size <= MAX_SIZE)) return error("invalid size");

		double workEstimate = spotPattern.size()*sqr(spotPaintRadius(fwhm, lorentzFraction, gamma));
		if (workEstimate > 5e6) {
			errorField.setText("Warning: Possibly slow. Suggestion: "+
					(lorentzFraction == 0 ? "reduce FWHM" : "set Lorentzian fraction to 0"));
			errorField.setForeground(Color.BLUE);
			return true;
		}
        return error(null); //deletes any error message and return true
    }

	/** Shows the error in the error message line and returns false if non-null */
	boolean error(String err) {
		errorField.setText(err == null ? "" : "ERROR: "+err);
		errorField.setForeground(err == null ? Color.BLACK : Color.RED);
		return err == null;
	}

	/** Enable/disable the fields not used for mofifying a stack */
	void enableCreateOnlyFields(boolean b) {
		for (Component c : createOnlyComponents)
			c.setEnabled(b);
	}

	/** Ends preview of the target image/stack */
	void revertPreview() {
		targetStack.restoreFwhmEtc();
		updateImage(targetImp, targetStack);
	}

    /** Saves the FWHM, Lorentzian fraction, and gamma */
    public void saveFwhmEtc() {
        fwhmSaved = fwhm;
        lorentzFractionSaved = lorentzFraction;
        gammaSaved = gamma;
    }

    /** Restores the previously saved FWHM, Lorentzian fraction, and gamma */
    public void restoreFwhmEtc() {
        fwhm = fwhmSaved;
        lorentzFraction = lorentzFractionSaved;
        gamma = gammaSaved;
    }

	/** Sets fwhm, Lorentz fraction, and gamma */
	public void setFwhmEtc(double fwhm, double lorentzFraction, double gamma) {
		this.fwhm = fwhm;
		this.lorentzFraction = lorentzFraction;
		this.gamma = gamma;
	}

	static void updateImage(ImagePlus targetImp, LEED_Movie_Simulation targetStack) {
		ImageProcessor ip = targetImp.getProcessor();
		int slice = targetImp.getSlice();
		float[] pixels = targetStack.getPixels(slice);
		ip.setPixels(pixels);
		double min = ip.getMin(), max = ip.getMax();
		double currentMax = targetStack.getMaxIntenisity();
		if (currentMax < max)
			ip.setMinAndMax(min*currentMax/max, currentMax);
		else if (currentMax > 10*max)
			ip.setMinAndMax(min, 0.1*currentMax);
		targetImp.updateAndDraw();
	}

	/** Returns an ImageProcessor for slice n, where 1 <= n <= nSlices */
	public ImageProcessor getProcessor(int n) {
		float[] pixels = getPixels(n);
		FloatProcessor ip = new FloatProcessor(size, size, pixels);
		ip.setSliceNumber(n);
		return ip;
	}
		
	/** Creates the image for slice n, where 1 <= n <= nSlices*/
	public float[] getPixels(final int n) {
		final double energy = getEnergy(n);
		final int spotR = spotPaintRadius(fwhm, lorentzFraction, gamma);
		final double mid = 0.5*(size -1);
		final double rRange = size*0.475 - 5;
		double theScale = rRange/maxRadius;
		if (!leemMode) theScale
		 *= 1./Math.sqrt(energy);
		final double scale = theScale;
		//count the spots
		int nSpots = 0;
		for (int i=0; i<ivData.data.length; i++) {
			double intensity = ivData.data[i][startEIndex+n-1];
			if (intensity > 0) {
				int spot = ivData.spotIndices[i];
				if (addEquivalent)
					nSpots += spotPattern.getAllSpotsForGroup(spotPattern.getGroup(ivData.spotIndices[i])).length;
				else
					nSpots++;
			}
		}
		if (nSpots == 0) return new float[size*size];
		//paint the spots
        final AtomicInteger iData = new AtomicInteger(0);          //prepare parallelization: spots are painted in parallel
        int nCPUs = Runtime.getRuntime().availableProcessors();
        int nThreads = Math.min(nCPUs, nSpots/4 + 1);
        final Callable[] callables = new Callable[nThreads];
        final float[][] thPixels = new float[nThreads][size*size];
        for (int t = 0; t < nThreads; t++ ) {
			final int it = t;
            callables[t] = new Callable() {
                final public Boolean call() {
					float[] pxls = new float[size*size];
					thPixels[it] = pxls;
					while (true) {
						int i = iData.getAndIncrement();
						if (i >= ivData.data.length) break;
						double intensity = ivData.data[i][startEIndex+n-1];
						if (intensity > 0) {
							int spot = ivData.spotIndices[i];
							if (addEquivalent) {
								int[] eqSpots = spotPattern.getAllSpotsForGroup(spotPattern.getGroup(spot));
								for (int sp : eqSpots)
									paintSpot(pxls, scale, mid, spotR, spotPattern.getKx(sp), spotPattern.getKy(sp), intensity);
							} else
								paintSpot(pxls, scale, mid, spotR, spotPattern.getKx(spot), spotPattern.getKy(spot), intensity);
						}
					}
					return null;
				}
			};
		}
		ThreadUtil.startAndJoin(callables);

		float[] pixels = nThreads == 1 ? thPixels[0] : new float[size*size];
		float maxValue = 0;
		for (int p=0; p<pixels.length; p++) {
			if (nThreads > 1) {			
				//combine images created by all threads
				float sum = 0f;
				for (int t=0; t<nThreads; t++)
					sum += thPixels[t][p];
				pixels[p] = sum;
			}
			float v = pixels[p];
			if (gamma == 2) {
				v = (float)Math.sqrt(v);
				pixels[p] = v;
			} else if (gamma != 1) {
				v = (float)Math.pow(pixels[p], 1./gamma);
				pixels[p] = v;
			}
			if (v > maxValue)
				v = maxValue;
		}
		maxIntensity = maxValue;
        return pixels;
    }

	/** Paints a spot */
	void paintSpot(float[] pixels, double scale, double mid, int spotR, double kx, double ky, double intensity) {
		double x = mid + kx*scale;
		double y = mid - ky*scale; 	//ImageJ has y top down, we want bottom up
		if (fwhm <= 1) {		//small spot should not hide between pixels
			x = Math.round(x);
			y = Math.round(y);
		}
		int xmin = (int)Math.round(x-spotR);
		if (xmin < 0) xmin = 0;
		int xmax = (int)Math.round(x+spotR);
		if (xmax >= size) xmax = size - 1;
		int ymin = (int)Math.round(y-spotR);
		if (ymin < 0) ymin = 0;
		int ymax = (int)Math.round(y+spotR);
		if (ymax >= size) ymax = size - 1;
		
		for (int iy=ymin; iy<=ymax; iy++) {
			int p = iy*size + xmin;
			for (int ix=xmin; ix<=xmax; ix++, p++) {
				double rsqr = (y-iy)*(y-iy) + (x-ix)*(x-ix);
				if (rsqr > spotR*spotR) continue;
				double v = (1. - lorentzFraction) * Math.exp(rsqr*(-4*Math.log(2)/sqr(fwhm)));
				if (lorentzFraction > 0) {
					v += lorentzFraction * (1./(1 + rsqr*sqr(2./fwhm)) - 1./(1 + sqr(spotR*2./fwhm)));
				}
				v *= intensity;
				pixels[p] += (float)v;
			}
		}
	}

	/** Returns the radius around a spot where painting is needed */
	static int spotPaintRadius(double fwhm, double lorentzFraction, double gamma) {
		return (int)Math.ceil(fwhm*(lorentzFraction == 0 ? 3 : 10*gamma));
	}

	/** Returns the maximum intensity */
	public double getMaxIntenisity() {
		return maxIntensity;
	}

    /** Calculates and returns the maximum radius from the (0,0) spot */
    double getMaxRadius() {
		double maxRSqr = 0;
		for (int n=1; n<=nSlices; n++) {
			double energy = getEnergy(n);
			for (int i=0; i<ivData.data.length; i++) {
				if (ivData.data[i][startEIndex+n-1] > 0) {
					int spot = ivData.spotIndices[i];
					double kx = spotPattern.getKx(spot);
					double ky = spotPattern.getKy(spot);
					double rSqr = kx*kx + ky*ky;
					if (!leemMode) rSqr *= 1./energy;
					if (rSqr > maxRSqr)
						maxRSqr = rSqr;
				}
			}
		}
		return Math.sqrt(maxRSqr);
	}

	/** Returns whether we should add symmetry-equivalent beams */
	boolean shouldAddEquivalent() {
		if (!spotPattern.hasEquivalentBeams())
			return false;
		LeedIntegerArray groups = new LeedIntegerArray(ivData.spotIndices.length);
		for (int i=0; i<ivData.spotIndices.length; i++) {
			int spot = ivData.spotIndices[i];
			int group = spotPattern.getGroup(spot);
			if (groups.indexOf(group) >= 0)
				return false;		//at least one group appears twoce or more often
			else
				groups.add(group);
		}
		return true;
	}

	/** returns the energy for a given stack slice */
	double getEnergy(int n) {
		return ivData.energies[startEIndex+n-1];
	}

    /** Returns the image width of the virtual stack */
    public int getWidth() {
        return size;
    }

    /** Returns the image height of the virtual stack */
    public int getHeight() {
        return size;
    }

    /** Returns the number of images in this virtual stack */
    public int getSize() {
        return nSlices;
    }

    /** Adds a slice to this virtual stack, unsupported */
    public void addSlice(String sliceLabel, Object pixels) { addSlice(""); }
    public void addSlice(String sliceLabel, ImageProcessor ip) { addSlice(""); }
    public void addSlice(String sliceLabel, ImageProcessor ip, int n) { addSlice(""); }
    public void addSlice(String fileName) {
        throw new RuntimeException("Adding Slices NOT Supported");
    }

    /** Deletes the specified image from this virtual stack, unsupported */
    public void deleteSlice(int n) {
        throw new RuntimeException("Deleting Slices NOT Supported");
	}

    /** Returns the label of the specified slice in this virtual stack (if it is one). */
    public String getSliceLabel(int n) {
        return "E="+IJ.d2s(getEnergy(n),1);
    }

	/** square of a number */
	double sqr(double x) { return x*x; }
}
