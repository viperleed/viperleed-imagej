import ij.IJ;

/** This class contains help messages for the ViPErLEED Spot Tracker.
 */

/** This code is part of the ViPErLEED package for LEED I(V) analysis.
 *  Licensed under GNU General Public License v3.0 or later (GPL-3.0-or-later),
 *  https://www.gnu.org/licenses/gpl-3.0.html
 *  The authors may decide later to put part of the auxiliary code in this work into the public domain,
 *  to allow incorporation into ImageJ if desired (ImageJ is in the public domain).
 *  The help text is licensed under the Creative Commons Attribution 4.0 (CC BY 4.0) license.
 *  When using and/or modifying this program for scientific work, please cite
 *  the paper describing it:
 *  M. Schmid, F. Kraushofer, A. M. Imre, T. Kißlinger, L. Hammer, U. Diebold, and M. Riva,
 *  ViPErLEED package II: Spot tracking, extraction and processing of I(V) curves,
 *  Phys. Rev. Research, 2024. 
 *  @author Michael Schmid, IAP/TU Wien, 2020-2024
 */

public class LeedSpotTrackerHelp {

    //The following texts are appended to form the main help dialog. Parts are used for subdialogs.
    private static final String PREHTML =
            "<html>";
    private static final String HELP_TITLE_TRACKER =
            "<h1>ViPErLEED Spot Tracker</h1>"+
            "<p><b>This ImageJ plugin extracts intensity vs. energy curves [a.k.a. <i>I</i>(<i>V</i>) curves] from "+
            "low-energy electron diffraction (LEED) movies.</b></p>"+
            "<p>Version: "+LEED_Spot_Tracker.VERSION+"</p>"+
            "<p>This documentation contains a description of the GUI elements, "+
            "<a href='#spottrackerHints'>useful hints</a> "+
            "and <a href='#spottrackerLicense'>license information</a>.</p>";
    private static final String HELP_INTRO_TRACKER =
            "<h2>Input Images/Stacks</h2>"+
            "<dl>"+
            "<dt>Main Input Stack (required)</dt>"+
            "<dd>The movie with the LEED images (must have evenly spaced electron energies, "+
            "unless intensity vs. time or some other data channel at constant energy should be measured).</dd>"+
            "<dt>Dark Frame and Flat Field (optional)</dt><dd>";
    private static final String HELP_DARKFLATEXPLAIN =
            "These images should be used to correct for the dark current of the camera (and background illumination, if any) "+
            "as well as inhomogeneities of the LEED screen, optics, and camera. The <em>flat field</em> is an image "+
            "(or image stack) of a diffuse diffraction pattern with no spots, only showing the inhomogeneities "+
            "of the optics and camera. The flat field is typically recorded by moving a polycrystalline "+
            "sample holder to the sample position; make sure that the distance from the electron source "+
            "is exactly the same as that of the sample surface. (Since the intensity is spread out over all the screen, "+
            "you may use a higher beam current or longer exposure than for the LEED movie of the sample.) "+
            "A flat field is especially important if the LEED image has a high "+
            "background intensity, e.g., for sample temperatures around or above the Debye temperature of the sample. "+
            "<em>Dark frame</em> images (when present) are subtracted from the input stack and flat field before "+
            "the flat-field correction. The dark frame is best measured with a high negative Wehnelt voltage, "+
            "suppressing the beam, but keeping light emitted by the filament as well as any glow due to (field-emitted) "+
            "electrons reaching the fluorescent screen.<br>"+
            "The formula is<br>"+
            "<pre>            corrected_pixel_value = (input - dark)/(flat - dark2)</pre>"+
            "where <tt>dark2</tt> is the dark frame for the flat field (usually the same as the dark frame "+
            "of the main input, unless the main input and flat field differ in the exposure time or gain). "+
            "If a dark frame is not given, nothing is subtracted there.<br>&nbsp;&nbsp;&nbsp;&nbsp;"+
            "The <tt>(flat - dark2)</tt> values must be positive in the area defined by the mask (see below). "+
            "In case of negative or zero values (e.g., due to excessive noise), attempts to track (and measure) the LEED spots "+
            "will result in an error message.<br>&nbsp;&nbsp;&nbsp;&nbsp;"+
            "For the dark and flat images, one can provide a stack. "+
            "If the stack has the same number of images as the main input, it is used image-by-image for the "+
            "corresponding main input image. There is no check whether the energies are the same as for the main stack; "+
            "this is your responsibility. Preprocessing of these correction image stacks (e.g., averaging over a "+
            "few energies to reduce noise, 2D fitting or at least intensity normalization for the flat field) "+
            "should be selected with <a href='#darkFlatProcessing'>Dark &amp; Flat Processing</a>. "+
            "If the number of images of a dark or flat stack differs from that of the main input stack, they are averaged "+
            "and this average is used for all images of the main input. "+
            "(If the <tt>flat</tt> is a single image or an average, the <tt>dark2</tt> also uses the average over the stack, "+
            "even if it has the same number of images as the main input.) "+
            "The type of dark & amp;flat processing is displayed in abbreviated form next to the "+
            "<a href='#darkFlatProcessing'>Dark &amp; Flat Processing</a> button. "+
            "Some sort of averaging (to reduce the noise) is advisable, especially for the flat field.<br>&nbsp;&nbsp;&nbsp;&nbsp;"+
            "Note that image processing with dark and flat requires a <em>mask</em> (see below). "+
            "The mask defines the area where the correction will be performed. "+
            "Without a proper mask, the correction will result in arbitrary values outside the LEED screen; "+
            "this will prevent proper setting of brightness &amp; contrast and may make the LEED screen appear black or dull.";
    private static final String HELP_INTRO_TRACKER2 =
            "</dd><dt>Mask (required)</dt><dd>An 8-bit image (&quot;binary image&quot;) with white pixels "+
            "for the background and black for the foreground. The foreground area determines the screen area "+
            "where intensities may be measured. The mask image must NOT contain pixel values other than 0 or 255. "+
            "You can use <i>Analyze&gt;Histogram</i> to see which pixel values are present in an image; "+
            "in the histogram panel click on &quot;Log&quot; (logarithmic <i>y</i> axis scale) to see also values that occur "+
            "for a few pixels only. You can use <i>Image&gt;Adjust&gt;Threshold</i> (and &quot;Apply&quot;) "+
            "to ensure that all pixel values are 0 or 255. Images with in-between values cannot be selected as a mask.<br>&nbsp;&nbsp;&nbsp;&nbsp;"+
            "You have to create the Mask image; start with <a href='#createMask'>More&gt;&gt;Create Mask</a>. "+
            "As soon as a mask has been selected, it is shown as an outline on the image stack used for spot tracking. "+
            "Editing of the mask image will be reflected by the outline almost immediately.</dd>"+
            "</dl>"+
            "<p>All these images or image stacks must be open in ImageJ; you can use <a href='#openImagesMovies'>More&gt;&gt;Open Images/Movies...</a> "+
            " to open them (and select them as the input, if appropriate). Alternatively, you can use the "+
            "<i>Plugins&gt;ViPErLEED&gt;Open LEED Movie...</i> command (which can also create a table of metadata when opening a movie) "+
            "or the various ImageJ <i>File&gt;Open</i> and <i>File&gt;Import</i> commands "+
            "(e.g. <i>File&gt;Import&gt;Image Sequence...</i>) "+
            "for combining multiple images into an image stack (i.e., a movie). "+
            "Image stacks may be in memory or (if the open/import command allows) they can be <em>virtual stacks</em>, "+
            "i.e., the active image is read from disk as required, but not always kept in memory (needs less RAM).<br>&nbsp;&nbsp;&nbsp;&nbsp;"+
            "All these images or stacks must have the same width and height (in pixels). "+
            "Do not close these images or stacks while using the ViPErLEED Spot Tracker.</p>";
    private static final String HELP_DARKFLATMENU =
            "<h2><a name='darkFlatProcessing'>Dark &amp; Flat Processing</a></h2>"+
            "<p>If the dark frame and flat field (and the dark2; "+
            "this is the dark frame subtracted from the flat field) are stacks with "+
            "the same number of slices as the main input stack, one can select several options:<br>"+
            "- No averaging: The correction images should be applied slice by slice. In this case, "+
            "any noise of the correction images will increase the noise of the result.<br>"+
            "- Average all: Best noise suppression, but no energy dependence (usually not suitable for the flat field).<br>"+
            "- Linear fit: This provides good noise suppression, but keeps linear drift "+
            "with energy or time. This option can be useful especially for the dark frames.<br>"+
            "- Smoothing: Select this option for averaging over a few energies (usually used for the flat field). "+
            "The <i>Smooth amount</i> determines the noise suppression. "+
            "(For a smooth amount of <i>n</i>, the noise suppression is roughly equivalent to that of averaging over <i>n</i> "+
            "stack slices, but the algorithm is actually different from moving averages, with better suppression "+
            "of rapid fluctuations and less memory use.)</p>"+
            "<p><i>Flat field processing type</i>: Before applying the flat field "+
            "(and after subtraction of its dark frame, if provided) to the main image, "+
            "in most cases further processing of the flat field should be applied: "+
            "If the flat field is a stack with 1:1 correspondence to the input "+
            "(i.e., same number of slices), one should select at least normalization, i.e., keeping the average "+
            "intensity of the stack (across the mask area) constant. "+
            "Furthermore, usually the flat-field intensity decreases "+
            "towards the edge, because diffuse scattering intensity from a polycrystalline sample tends "+
            "to peak at a scattering angle of 180&deg; (backscattering). Then one should correct the (background-corrected)"+
            " flat field by a fit of a 2nd-order or 4th-order polynomial to its logarithm (for 2nd-order, "+
            "this would approximately correspond to division by a best-fit Gaussian; "+
            "the 4th-order fit is usually preferable, but slower to compute). "+
            "Polynomial fits should be only applied if the long-range intensity variations of the flat field "+
            "are due to the scattering process in obtaining the flat field, not if these variations are due "+
            "inhomogeneous efficiency of the LEED screen or vignetting of the camera.<br>&nbsp;&nbsp;&nbsp;&nbsp;"+
            "One can examine the flat field (after all processing steps) by "+
            "<i>Show processed flat field</i>. The output of <i>Show processed flat field</i> "+
            "can be also used for further processing of the flat field by any image processing method of your choice, "+
            "then the result of this processing can be used as a flat field further on. "+
            "In this case, do not forget to deselect the dark frame for the flat field "+
            "and deactivate the polynomial fit (unless you really want it).</p>"+
            "<p>Note that processing steps that require all images of an input stack, such as averaging over all "+
            "images or a linear fit over all images of a stack require some computing time. "+
            "It may take several seconds until the processed result becomes available; this may also "+
            "cause a delay in opening the Spot Tracker main panel when it automatically selects the input "+
            "images or stacks. Expect delays especially if the inputs are supplied as virtual stacks "+
            "(i.e., not kept in memory but read from disk on demand).</p>";
    private static final String HELP_SETENERGIES =
            "<h2><a name='setEnergies'>Set Energies, I0, t</a></h2>"+
            "<p>The energy range of the input stack, the beam current I0 (for normalization), and the time "+
            "(with respect to the first slice) are usually read from the slice labels of the main input stack; they are "+
            "read in when reading the main image stack. If these data were not available, one can set them manually, "+
            "as a linear function or from a table that has been opened by ImageJ. [Note that one can open <tt>.csv</tt> "+
            "and <tt>.tsv</tt> files (comma- or tab-separated values) as tables with ImageJ, if they contain exactly "+
            "one header line.] "+
            "The &quot;Set Energies, I0, t&quot; dialog also provides an entry for I00, the offset of the beam current I0. "+
            "I00 can be also read from the I0 values of the dark frame(s) (if present); "+
            "either using the individual values for each energy or a linear fit "+
            "(if the dark frames are a stack with 1:1 correspondence to the main input stack). "+
            "I00 may be read from the dark frame(s) only if they were obtained without emitted electron current "+
            "(i.e., by setting the Wehnelt voltage to a highly negative value, not by setting the screen voltage to 0). "+
            "If I00 is provided, the actual beam current is calculated as <tt>I0 - I00</tt>. "+
            "Usually, the beam current should be smoothed; the amount of smoothing can be entered (0 for no smoothing). "+
            "The amount of smoothing is given as the number of points in a moving-averages filter "+
            "with the same noise suppression; actually a modified sinc kernel [<a href= '#msSmooth'>2</a>] is used. "+
            "I0, I00, and the smoothed and I00-corrected beam current are plotted after closing this input window "+
            "if 'Plot I0' is selected. The plot does not include corrections based on the background intensity "+
            "measured during spot tracking (after spot tracking, you can use <a href='#morePlot'>More&gt;&gt;Plot...</a> for this).</p>"+
            "<p>This input window also has an option to set LEEM mode. In LEEM mode, it is assumed that "+
            "the spots do not move with energy as they would in normal LEED experiments.</p>"+
            "<p>For experiments with constant electron energy, such as studies of phase transitions "+
            "(over temperature) or evolution over time, one can specify an <i>x</i> axis other than the energy.</p>"+
            "<p>In the Spot Tracker panel, the data present are indicated to the right of the 'Set Energies, I0, t' button; the sources of "+
            "these data are given in square brackets (if not the input stack): <tt>[T]</tt> table, "+
            "<tt>[!]</tt> manual input of a linear relationship, <tt>[?]</tt> kept previous data. "+
            "For I00, <tt>[D]</tt> indicates values from the dark frames (lowercase <tt>[d]</tt> if I00 is a linear fit).</p>";
    private static final String HELP_PATTERNFILE =
            "<h2><a name='patternFile'>Set Pattern File</a></h2>"+
            "<p>This file describes which beams (spots) are expected and their relative positions in a Cartesian "+
            "coordinate system. This is a comma-separated (csv) file, one line for each diffraction maximum (beam), "+
            "with the beam label, beam h and k indices, Cartesian reciprocal-space coordinates "+
            "<i>g</i><sub>x</sub>, <i>g</i><sub>y</sub>, and a beam group "+
            "index (same for equivalent beams). Such a file can be created with the pattern simulator of the ViPErLEED GUI "+
            "(there, use Export from the File Menu). For creating the Spot Pattern file in the ViPErLEED GUI, "+
            "instead of manually entering the symmetry, you can also read the symmetry from the "+
            "<tt>experiment-symmetry.ini</tt> file created when running a LEED <i>I</i>(<i>V</i>) simulation with <tt>viperleed.calc</tt>.</p>";
    private static final String HELP_SETINDICES1 =
            "<h2><a name='setIndices'>Set Indices</a></h2>"+
            "<p>Firstly, an input window appears and asks you to select a stack slice "+
            "[an energy in usual <i>I</i>(<i>V</i>) measurements] where as many spots as possible are clearly visible, "+
            "including spots near the edge of the screen. Use the slider below the '<tt>..._SpotTracking</tt>' image stack "+
            "to select the energy (For energy selection, you cannot use one of the image stacks supplied as the input, "+
            "unless it is coupled with the '<tt>..._SpotTracking</tt>' stack via <i>Analyze&gt;Tools&gt;Synchronize Windows</i>). "+
            "You may use the <tt>Auto</tt> button of the "+
            "<i>Image&gt;Adjust&gt;Brightness&amp;Contrast panel</i> (B&amp;C; keyboard shortcut "+LeedUtils.CTRL+"-SHIFT-C when the stack "+
            "is in the foreground) or the sliders of the B&amp;C "+
            "panel to adjust the contrast. If many spots are not recognized (not marked by "+
            "a circle), or many circles do not mark diffraction maxima, you can adjust the noise rejection. "+
            "Note that also an improper value of the <a href='#setRadius'>integration radius</a> may prevent recognition of the spots. "+
            "Having a few 'false' spots marked does not hurt as long as these are not close to the position where "+
            "a diffraction maximum would be expected.</p>"+
            "<p>When spot indices have been defined previously and fit the maxima, the maxima will be labeled "+
            "and you can click 'Labels are OK' if the identification is correct. Otherwise, you have to manually name "+
            "one or a few spots in the next step.</p>";
    private static final String HELP_SETINDICES2 =
            "<p>For naming the spots, a second ";
    private static final String HELP_SETINDICES3 =
            "input window asks you to select a spot by clicking on it and enter its "+
            "h, k indices (for superstructures, also fractions are allowed). The spot must be listed in the "+
            "<a href='#patternFile'>pattern file</a>. After pressing <i>Set</i>, the program tries "+
            "to name the other spots. Examine the outcome (you may select the image and press '+' or "+LeedUtils.CTRL+"-'+' "+
            "to zoom in at the cursor position). If the spots are named correctly, press OK; otherwise select "+
            "another spot and enter its designation (best, choose a spot that has not been correctly identified). "+
            "For LEED patterns recorded with oblique incidence (&gt; 10&deg; off-normal) "+
            "it will be necessary to enter several spots. In that case, start with naming a few nearby spots. "+
            "In case of a mistake, you can click on a spot and rename it. As soon as the program has tried to identify "+
            "the spot pattern, the input window also shows the polynomial order it uses to convert reciprocal-space "+
            "coordinates to pixel coordinates and the root-mean-square (rms) deviation of the spots from the fit "+
            "in pixels. For small images (e.g., 512 x 512 pixels), the rms deviation should not be much more than about "+
            "two pixels; otherwise carefully check whether the assignment is correct (especially if spots are dense "+
            "in some places, but with larger gaps in between, as it happens for some superstructure cells). "+
            "Especially in case of distorted images, one can also "+
            "try to adjust the search tolerance (larger values allow for more distortions, but are more likely to "+
            "result in incorrect spot assignment).</p>";
    private static final String HELP_SETRADIUS =
            "<h2><a name='setRadius'>Set Integration Radius</a></h2>"+
            "<p>The integration radius should be selected such that the spots are fully enclosed by a circle "+
            "with the given radius. For low energies, where additional broadening due finite domain and "+
            "terrace sizes (as well as the electron source) may occur, "+
            "the integration radius can be given separately for integer and superstructure spots.</p>"+
            "<p>The output after spot tracking contains a <a href='#spotRadii'>plot with the measured spot radii</a> "+
            "and a line with half the integration radius as a function of energy. "+
            "The point cloud should be mostly below the line. "+
            "If there is no danger of the background area (see below) overlapping with neighboring spots, "+
            "and noise is low, set an integration radius about 2.5 times the spot radii measured during spot tracking. "+
            "(Then, the line in the plot should be at ~1.2 times the typical values of the point cloud.) Large integration radii "+
            "lead to lower noise for strong spots, but poorer performance for weak ones (more noise and more errors "+
            "in case of an uneven background). "+
            "One may use up to 3 times the typical spot radii if (i) the noise of the images is very low, and (ii) the background "+
            "due to inelastic scattering is very low and flat, and (iii) all spots are bright. "+
            "Even in case of very close spots and high noise, the radius values should not be less than ~2 times "+
            "the typical spot radii. This means that the line in the "+
            "<a href='#spotRadii'>spot radii plot</a> should still be above the center of the point cloud. "+
            "Usually it is a good idea to try spot tracking with different radii and compare the "+
            "<a href='#qualityPlot'><i>I</i>(<i>V</i>) quality plots</a> to find out which value works best. "+
            "Keep in mind that the radii must be chosen such that the background area of each spot "+
            "does not overlap with any of the neighboring spots.</p>"+
            "<p>Apart from the radius values, there is a selection for the shapes of the integration area and the background "+
            "(see the paper about the Spot Tracker [<a href='#paper'>1</a>]):<br>"+
            "- <em>Concentric Circles</em> means that the integration area is a circular disk and the background "+
            "intensity is averaged over a circular ring with outer radius \u221a"+"2 times the radius of the integration circle. "+
            "This option is best if spots are very close at high energies, and it also provides the best suppression "+
            "of local variations in the diffuse background intensity due to scattering by phonons.<br>"+
            "- <em>Oval background</em> also uses a circular integration disk, but an elliptical outer border"+
            " of the background area, with a semi-major axis twice the radius "+
            "of the integration circle and the semi-minor axis equal to that radius. "+
            "(The background area has the shape of two large ears.) "+
            "This option is best in the case of a background intensity varying with the distance from the screen center "+
            "and it also maximizes the energy range that can be evaluated for most spots, as the background area "+
            "does not touch the inner or outer border when the circular background would already touch it.<br>"+
            "- <em>Azimuth blur</em> is for spots blurred in azimuthal (tangential) direction, due to small variations of "+
            "the azimuthal orientation of the domains or crystallites. In this case, the integration disk is an ellipse "+
            "with the semiminor axis being the radius defined above. The semimajor axis is oriented in tangential direction "+
            "[perpendicular to the direction to the (0,0) spot] and gets enlarged "+
            "to fit the azimuthal blur angle entered. "+
            "For near-circular spots [near the (0,0) spot], the background area has a circular outer border, "+
            "with \u221a"+"2\u00b7<tt>radius</tt>. "+
            "With increasing ellipticity (integration ellipse elongated in tangential direction) "+
            "the outer background radius remains constant as long as the integration ellipse fits inside "+
            "the circular background outline. "+
            "Then the background outline becomes stretched in tangential direction, i.e., it becomes an ellipse. "+
            "The total area remains \u221a"+"2 times the integration area, thus the background area is about 41% of the "+
            "integration area. The background intensity is measured near the borders in the radial direction from the spot. "+
            "Thus, the background intensity is sensitive to nonlinear background variations over the radius, but insensitive to "+
            "neighboring spots in the azimuthal direction. "+
            "(This is the opposite to the oval background, where the background is evaluated in tangential direction from the spot.)</p>"+
            "<p>For both <em>Oval background</em> and <em>Azimuth blur</em>, note that there is no thorough checking "+
            "for whether the background area of a spot overlaps with the integration disk of any of its neighbors. "+
            "Thus, use these modes only if the spots are well-separated up to the highest energy.</p>"+
            "<p>The type and size of the integration and background areas are indicated at the top right of the SpotTracking stack. "+
            "(Not updated while the dialog is open.)</p>";
    private static final String HELP_TRACKSPOTS1 =
            "<h2><a name='trackSpots'>Track Spots</a></h2>"+
            "<p>Starts spot tracking and intensity measurements (may take up to a minute). "+
            "Tracking starts at the energy where the spot indices "+
            "have been entered via <a href='#setIndices'>Set Indices</a>, then the program tries to follow "+
            "the beams and to discover other beams listed in the <a href='#patternFile'>Pattern File</a> "+
            "if they become visible.</p>";
    private static final String HELP_TRACKSPOTS2 =
            "<p>The numeric values of the parameters determining the tracking process can be usually left "+
            "at their defaults. (When not measuring as a function of energy, the numeric parameters given in 'eV' below "+
            "rather refer to the number of frames, i.e., data points of the output.)</p>"+
            "<dl>"+
            "<dt>Noise Rejection (default "+IJ.d2s(LeedParams.getDefaultValue(LeedParams.MINSIGNIFICANCETRACK),1)+"):</dt>"+
            "<dd>Determines how strong a spot must be to be taken as valid. "+
            "The value sets the threshold of the significance (red curve in the "+
            "<a href='#positionDeviations'>Position Deviation</a> plot generated by spot tracking). "+
            "Increase this value if spots get trapped by screen inhomogeneities during tracking "+
            "or if spot tracking finds spots where there are none at any energy. If the noise threshold is too low, "+
            "stray spots might also wander into another (real, but weak) spot, inducing deletion of that spot.</dd>"+
            "<dt>Search beam unseen for (default "+IJ.d2s(LeedParams.getDefaultValue(LeedParams.SEARCHAGAINEV),0)+" eV):</dt>"+
            "<dd>If a beam has not been detected over an energy "+
            "interval of more than this value, information about its deviation from the predicted position is discarded and "+
            "the program tries to find it as if it had not been detected before.</dd>"+
            "<dt><a name='smoothDeviations'>Energy range for x, y smoothing (default "+
            IJ.d2s(LeedParams.getDefaultValue(LeedParams.POSITIONAVERAGINGEV),0)+" eV)</a>:</dt>"+
            "<dd>Before measuring intensities, the deviations of the spot positions from the predicted value are smoothed, "+
            "to obtain a continuous movement of the integration circle. "+
            "A larger value of this parameter results in smoother movement of the integration circle. "+
            "In case of local distortions of the LEED pattern, a large value may result in inaccurate positions, however. "+
            "Measurements down to very low energies "+
            "(&lt;30 eV) may require a smaller value for this parameter. This parameter also determines how far the "+
            "position may be extrapolated from the first/last detection of a spot. (Set 'Min. energy range per beam', below, "+
            "to 0 for unlimited extrapolation.)</dd>"+
            "<dt>Min. energy range per beam (default "+IJ.d2s(LeedParams.getDefaultValue(LeedParams.MINRANGEEV),0)+" eV):</dt>"+
            "<dd>Beams with a total energy range less than this limit are ignored. "+
            "When this parameter is set to 0 there is no such limit, and the program tries to measure all spots "+
            "in the spot pattern file, irrespective of whether they are above the significance threshold. "+
            "(see 'Noise Rejection', above.) "+
            "Then, for spots that are below the noise threshold at all energies, the positions "+
            "are inferred from the neighboring stronger spots. This is useful, e.g., for checking whether minority areas with a "+
            "given superstructure are present. "+
            "In this case, check the positions and make sure you don't mistake small background variations for a real signal.</dd>"+
            "<dt>Subtract background of bright neighbor spots:</dt><dd>Since spots often have a long-range background "+
            "(tails of the spot profile), the background for nearby spots is not uniform. "+
            "Use this option only if the spots show circular symmetry, not for elongated spots. "+
            "This option fits a 1/<i>r</i>&#178; background in the vicinity of bright spots "+
            "and subtracts that background before evaluating the intensity of weaker spots. "+
            "This option is not available if spot shape 'Azimuth blur' is selected. "+
            "Depending on the spot profiles, for good results, this option may need a larger integration radius "+
            "than the lower limit of twice the measured spot radii. "+
            "The background-subtracted stack will be displayed and can be examined if 'Create debug info' "+
            "is selected.</dd>"+
            "<dt>Apply fast background changes to I0:</dt>"+
            "<dd>This option is only present if I0 is available and I0 smoothing is on, see "+
            "<a href='#setEnergies'>Set Energies, I0, t</a>. "+
            "It measures the average background intensity and corrects the (processed) beam currents "+
            "based on fast fluctuations of the background intensity. In this context, 'fast' means fluctuations that "+
            "would be filtered out by smoothing of I0. This option is useful if the background intensity can be "+
            "determined with less noise than I0. Typically, only fast variations of the background "+
            "intensity are due to I0 variations, slower changes of the background are due to energy-dependent scattering. "+
            "This option works well if the background intensity is high (temperature &gt; Debye temperature), I0 smoothing "+
            "is not much more than 5 eV/&Delta;<i>E</i> (where &Delta;<i>E</i> is the energy step), and if dark frames are available, "+
            "i.e., when zero light corresponds to a pixel value of zero. "+
            "The ViPErLEED hardware measures I0 with very low noise; for such data it may be better to disable I0 smoothing "+
            "in case of I0 jumps instead of using this option. "+
            "To examine whether this option works well, plot the processed I0 vs. (raw) I0 with <a href='#morePlot'>More&gt;&gt;Plot...</a> "+
            "after spot tracking and check whether the processed I0 is less noisy than the raw I0. "+
            "You can also plot the background intensity divided by the (unsmoothed, but I00-corrected) I0. "+
            "Use the &quot;Apply fast background changes to I0&quot; option only if background/I0 "+
            "is substantially more noisy (with logarithmic <i>y</i> axis) than the background intensity.</dd>"+
            "</dl>"+
            "<p>The other parameters in this input window only affect the <a href='#qualityPlot'>quality statistics plot</a> "+
            "for analysis of equivalent beams; "+
            "these options are only available for <i>I</i>(<i>V</i>) measurements, not for, e.g., time-dependent measurements. "+
            "These are the V0i parameter of Pendry's <i>R</i> factor  (strictly speaking, its absolute value), "+
            "the smoothing before calculating <i>R</i> factors (default "+IJ.d2s(LeedParams.getDefaultValue(LeedParams.SMOOTHEV),0)+" eV), "+
            "and the energy ranges handled separately (default "+IJ.d2s(LeedParams.getDefaultValue(LeedParams.SPLITEV),1)+" eV): "+
            "<i>I</i>(<i>V</i>) curves are split into regions roughly this size, "+
            "and for statistics, the <i>R</i> factor between symmetry-equivalent beams is evaluated for each of these regions. "+
            "(Splitting into ~100 eV sections is done because the intensities of the low- and high-energy regions "+
            "are often very different; the plot shows mutual <i>R</i> factors vs. intensity.)</p>"+
            "<p>One can also decide to create only two quick plots (<a href='#selectedCurves'>Selected IV Curves</a> and "+
            "<i>I</i>(<i>V</i>) Quality Statistics) instead of all plots. For large datasets, this avoids the time-consuming creation "+
            "of large plot stacks. This is useful if the analysis is only meant to check the alignment of the sample, "+
            "i.e., whether the <i>I</i>(<i>V</i>) curves of equivalent beams agree, or the impact of different parameter values on the quality "+
            "of the result. (When saving the data and 'Save Plots' "+
            "is selected, the full set of plots will be saved irrespective of the 'Quick plot only' setting.)</p>";
    private static final String HELP_TRACKSPOTS3 =
            "<p>Spot tracking sometimes creates a window named <tt>WARNING: Beam(s) uncertain (highlighted)</tt>. "+
            "In this case, check these beams (marked by thick light blue circles in the main spot-tracker image stack; "+
            "you may have to set an appropriate energy to see them) and use "+
            "<a href='#deleteHighlighted'>More&gt;&gt;Delete highlighted beams</a> if appropriate. "+
            "This warning may also come up if spot tracking has worked well, but the choice of the energy for "+
            "<a href='#setIndices'>Set Indices</a> did not allow the program to determine a good model for the "+
            "distortions of the LEED pattern.</p>"+
            "<p>Spot tracking creates several plots:</p>"+
            "<dl><dt><a name='selectedCurves'>Selected I(V) Curves</a></dt><dd>If there are symmetry-equivalent beams, "+
            "a plot of one set of such beams. "+
            "This is mainly for judging the alignment and residual electric and magnetic fields. "+
            " (The selection criteria for this set of curves are a large common energy range, "+
            "a large number of equivalent beams, and high intensity.)</dd>"+
            "<dt>I(V) curves</dt><dd>The stack of all <i>I</i>(<i>V</i>) curves extracted. Note that the beam groups are in square "+
            "brackets; <i>I</i>(<i>V</i>) curves belonging to the same group should agree (symmetry-equivalent beams). "+
            "Negative group numbers denote symmetry-forbidden beams; these may have a finite (true or apparent) intensity, e.g., "+
            "due to deviations from perpendicular incidence or an uneven background.</dd>"+
            "<dt><a name='positionDeviations'>Spot Position X/Y Deviations</a></dt><dd>For each beam, "+
            "a plot with the deviations of the measured spot position "+
            "from the position expected (as a function of energy), and this curve after smoothing as defined by the "+
            "&quot;<a href='#smoothDeviations'>Energy range for x, y smoothing</a>&quot; parameter. "+
            "The smoothed deviations are used to obtain the "+
            "positions of the integration disks for intensity measurement. Especially if the density of spots "+
            "is high or there are defects on the screen, these plots should be examined to check whether spot "+
            "tracking has been successful. A large jump or a position 'running away' from a certain energy on "+
            "will indicate a problem with tracking. "+
            "These plots also include a line with the 'significance' of the spots, i.e., the value used to "+
            "determine whether a spot is above the noise threshold (below the threshold, the significance "+
            "is plotted as 0).</dd>"+
            "<dt><a name='spotRadii'>Spot Radii</a></dt><dd>A scatter plot with all spot sizes measured, as a function of energy. "+
            "In case of a superstructure, the superstructure spots are marked by different colors than the integer ones. "+
            "The spots are assumed to be ellipses; the extension in radial and tangential direction is plotted. "+
            "Here, 'radial' refers to the direction to the screen center, except for 'azimuth blur' mode, where "+
            "the radial direction is defined as the direction to the (0,0) spot. "+
            "The spot size plotted is an approximation to the standard deviation &sigma; of a Gaussian peak. "+
            "Since the peaks are usually not Gaussians, the sizes depend somewhat on the integration radius. "+
            "The plot also shows a line with half the integration radius; it should trace the center of the point cloud "+
            "or be at larger radius values than the center line of the cloud. Ignore the outliers, "+
            "these are usually very weak spots. In &quot;azimuth blur&quot; mode, "+
            "the tangential size of the spots is scaled down, as if the image was distorted to make the integration ellipses circular. "+
            "Nevertheless, in &quot;azimuth blur&quot; mode, this criterion "+
            "is mainly relevant for the size &sigma; in the radial direction. It is less important in the azimuthal direction "+
            "because the background area does not extend from the ellipse in azimuthal direction. "+
            "See <a href='#setRadius'>Set Integration Radius</a>.</dd>"+
            "<dt><a name='qualityPlot'>I(V) Quality Statistics</a></dt><dd>This plot shows the <i>R</i> factors between symmetry-equivalent "+
            "beams in various ways:<br>"+
            "(1) A scatter plot, where the horizontal axis is the average intensity of the <i>I</i>(<i>V</i>) curve "+
            "or section of that curve (normalized to 1000 for the highest intensity occurring anywhere). "+
            "Typically, due to noise, weak beams have higher <i>R</i> factors between symmetry-equivalent beams, "+
            "i.e., the point cloud is &quot;higher&quot; on the left side of the plot. "+
            "If the experiment was done very well (no residual magnetic fields, sample surface exactly normal to the electron beam), "+
            "the points for the highest intensities will have mutual <i>R</i> values of less than 0.1.<br>"+
            "(2) A line that shows the <i>R</i> factor against the cumulative energy range of all curve-section pairs "+
            "with an <i>R</i> factor better than the given number. (If there are many symmetry-equivalent beams, this is more than "+
            "the energy range available up to that <i>R</i> factor because many pairs can be selected from one group of beams.) "+
            "Beams where no symmetry-equivalent curves are available are ignored.<br>"+
            "(3) If there are beams with negative intensities (after smoothing), "+
            "the plot shows red triangles: "+
            "For these points, the <i>x</i> axis gives the absolute value of the most negative intensity, "+
            "and the <i>y</i> axis gives the total energy range (per beam) where the intensity is negative. "+
            "Thus, the red symbols, if any, should be as far to the bottom left as possible.</dd>"+
            "<dt>Overall X/Y Deviations</dt><dd>This plot contains information on systematic deviations of the spot positions "+
            "from the expected position. It is based on 2D linear regression of the measured vs. expected position. "+
            "This fitting procedure yields an estimated position of the (0,0) spot (which is usually hidden by the electron source). "+
            "An energy dependence of the (0,0) position can be caused by, e.g., residual (magnetic) fields "+
            "or an off-axis position of the filament in the electron source. "+
            "The text and arrow at the top right indicate the displacement of the (0,0) spot at 100&nbsp;eV vs. "+
            "the extrapolated position at infinite energies (where the electron beam is not deflected). "+
            "This calculation assumes magnetic fields, which yield a deflection proportional to 1/&radic;<i>E</i>. "+
            "The direction of the arrow marks the direction of the movement of the (0,0) spot in the "+
            "'SpotTracking' movie with decreasing energy. (This is valid unless the plot window is resized, "+
            "which would distort the plot area and, hence, the arrow.)<br>&nbsp;&nbsp;&nbsp;&nbsp;"+
            "The plot also contains a line for the deviation of the scale factor (with respect to the usual 1/&radic;<i>E</i> scale) from 100%. "+
            "Small scale deviations can be due to the difference of the work functions of the filament and "+
            "that of the sample and its surrounding. "+
            "There is also a fit line for the scale factor, which assumes that deviations are due to such work function differences. "+
            "For conventional LEED optics without additional electric fields, "+
            "large scale deviations (more than ~2&ndash;3% at 50&nbsp;eV) "+
            "may indicate surface charging. The deviation of the scale factor is usually much larger for MCP-LEED, "+
            "where fringe-field electrodes are present.</dd></dl>";
    private static final String HELP_SAVEDATA =
            "<h2><a name='saveData'>Save Data</a></h2>"+
            "<p>Saves the results as csv files in a directory chosen by the user. </p>"+
            "<p>The user can choose a prefix (all file names start with this) and also select whether to "+
            "save plots (as created by 'Track Spots') and the processed image stack (usually a big file).</p>"+
            "<p>Files created:</p>"+
            "<ul>"+
            "<li>The main outcome: The intensities, corrected for the (processed) beam current I0, file <tt>&lt;prefix&gt;_IVcurves.csv</tt>. "+
            "This file should be processed with the LEED I(V) Curve Editor to create the experimental <i>I</i>(<i>V</i>) curves used for "+
            "comparison with the calculated ones.</li>"+
            "<li>The uncorrected intensities (not taking I0 into account) are written into <tt>_rawInt.csv</tt>. "+
            "For all beams tracked also the following information is "+
            "written as csv files: The <i>x</i> and <i>y</i> positions (in pixels) as well as the raw "+
            "and smoothed deviations of the positions from the expected ones "+
            "(<tt>_x</tt>, <tt>_y</tt>, <tt>_dx_raw</tt>, <tt>_dy_raw</tt>, <tt>_dx_smooth</tt>, <tt>_dy_smooth</tt>), "+
            "the significance (set to ~0 when it is below the noise threshold; file <tt>_signif</tt>), "+
            "the spot radii in radial and tangential direction (<tt>_rSize</tt>, <tt>_tSize</tt>; "+
            "see <a href='#spotRadii'>Spot Radii</a>), "+
            "as well as the intensity in the background region around the integration area and its standard deviation "+
            "(<tt>_backgr</tt>, <tt>_bsigma</tt>).</li>"+
            "<li>A <tt>.log</tt> file, which contains the parameters used for spot tracking.</li>"+
            "<li>If 'Save Plots' is selected, one file for each of the plots mentioned "+
            "in the <a href='#trackSpots'>'Track Spots' help</a> is saved.</li>"+
            "<li>If 'Save Stack' is selected, the 'spotTracking' image stack is saved, i.e., the input with dark/flat "+
            "corrections applied and labeled integration circles (a large file). "+
            "Otherwise, one image with beam-index designations (<tt>spotIndicesImage.tif.zip</tt>) is saved "+
            "(at the energy of <a href='#setIndices'>Set Indices</a>). "+
            "The image or stack is saved in .tif.zip format. It should be opened in ImageJ without prior unpacking.</li>"+
            "</ul>";
    private static final String HELP_MOREMENU =
            "<h2>More&gt;&gt; Popup Menu</h2>"+
            "<dl>"+
            "<dt><a name='openImagesMovies'>Open Images/Movies...</a></dt><dd>Displays a dialog for opening the input files. Files selected as "+
            "input that do not fit the requirements (e.g., wrong size, or a mask that is not a binary image) "+
            "will be opened, but not selected as the respective input. "+
            "Opening all files using this command is usually more convenient than opening and selecting them individually. "+
            "If your computer does not have enough RAM for all image stacks, "+
            "(including the processed one, after dark-frame and flat-field correction) "+
            "select &quot;Virtual Stack&quot;. Then only the images currently required will be read from disk.</dd>"+
            "<dt>Select stack slice for...</dt><dd>Asks the user for the energy (or other <i>x</i>-axis variable, e.g., 'time'; "+
            "see <a href='#setEnergies'>Set Energies, I0, t</a>) and selects the stack slice (i.e., the image) for this energy.</dd>"+
            "<dt><a name='createMask'>Create mask...</a></dt><dd>Facilitates creating a mask for a stack of LEED images. "+
            "If a flat field is available, the tool is based on the flat field as an input (averaging over the stack slices); "+
            "otherwise it uses the main input stack (the standard deviation of the intensity vs. energy for each pixel). "+
            "There are two sliders: Initially, move the threshold slider to obtain a smooth outline of the screen area "+
            "and of the edge of the electron source (with the arm holding it). "+
            "In the second step, you may shrink or grow that area by a few pixels to refine it. "+
            "In addition, if the program can guess the outline of the LEED screen, you can select "+
            "&quot;Limit to elliptical fit&quot; and you can try the checkbox to correct for the "+
            "radius dependence of the intensity, which is typically present in the input images. "+
            "To evaluate the result, you may use <i>Image&gt;Adjust&gt;Brightness&amp;Contrast</i> "+
            "to better see the border of the LEED image in the &quot;SpotTracking&quot; stack. "+
            "When done, select the mask created in the main Spot Tracker panel. "+
            "If required, use the standard "+
            "<a href='https://imagej.net/ij/docs/guide/146-19.html'>ImageJ selection tools</a> "+
            "(oval, polygon, freehand) "+
            "and the <i>Edit&gt;Fill</i>, <i>Clear</i> and <i>Clear Outside</i> commands to further refine the mask.</dd>"+
            "<dt>Highlight beams...</dt><dd>Shows the selected beam(s) with a thicker circle. Useful in case "+
            "of complex superstructures and to delete these from the output (see &quot;Delete highlighted beams...&quot;, below).</dd>"+
            "<dt><a name='deleteHighlighted'>Delete highlighted beams...</a></dt><dd>Deletes these beams completely or in a given energy range. "+
            "These beams will not appear in the output of a successive <a href='#saveData'>Save Data</a> operation "+
            "(but they remain circled and named in the stack). "+
            "One cannot undo deletion except by running <a href='#trackSpots'>Track Spots</a> again.</dd>"+
            "<dt>Un-highlight beams</dt><dd>Restores normal view after 'Highlight beams'.</dd>"+
            "<dt>Beam statistics</dt><dd>Creates a table with one line for each beam successfully tracked: "+
            "Beam group, number of data points, energy of the first and last data point, "+
            "highest significance, and the energy where the highest significance is achieved. "+
            "There is also a summary line with the total numbers and overall minima/maxima.</dd>"+
            "<dt><a name='morePlot'>Plot...</a></dt><dd>Creates a plot of various data available. Examples include (processed) I0 or "+
            "data for the beam(s) requested, such as intensity. See the list of .csv files in "+
            "<a href='#saveData'>Save Data</a> for the items available.</dd>"+
            "<dt>Close plots...</dt><dd>Closes all or a given subset of plots created by 'Track Spots'.</dd>"+
            "<dt>Undistort image/stack</dt><dd>Uses the &quot;Spot Tracking&quot; stack as an input and "+
            "applies the distortion corrections determined when the <a href='#setIndices'>spot indices were set</a>. "+
            "The undistorted image has the position of the (0,0) spot at the center. "+
            "If the outer border of the LEED screen in the undistorted image is far from elliptical this indicates that "+
            "there have not been enough spots close to the border to correctly determine the distortions. "+
            "The user can select whether to apply the correction only to the current slice "+
            "(the image currently shown) or the whole stack. In the latter case, except for LEEM mode, "+
            "there are two options for the image scale: "+
            "&quot;scale like input&quot; creates a stack where the spots are moving as in the input, while "+
            "&quot;fixed k-space scale&quot; de-magnifies low-energy images such that the spots of a given beam "+
            "are always at the same pixel position. Stack undistortion should be run after spot tracking; then the "+
            "overall deviations caused by residual fields are taken into account. "+
            "The output can be normalized, by dividing by the processed I0 (if available), to eliminate "+
            "brightness variations due to changes of the incident electron beam current. "+
            "The output of stack undistortion is a virtual stack with caching. "+
            "Thus, the images are created only when needed. This means that the output is quickly visible, "+
            "but some image operations using the undistorted output will be initially slow, "+
            "until all slices of that stack have been created. "+
            "This also means that the stack created may become partly or fully unavailable "+
            "when one of the input images/stacks is closed. "+
            "Note that averaging over the slices of an undistorted stack "+
            "(<i>Image&gt;Stacks&gt;Z Project...</i>) "+
            "requires ImageJ version 1.54j or later. "+
            "Further processing of the undistorted stack (e.g., background subtraction, resizing) "+
            "requires that you duplicate this stack.</dd>"+
            "<dt>List parameters</dt><dd>Creates a table with all numeric parameters; mainly for debugging. "+
            "When in macro-recording mode (<i>Plugins&gt;Macros&gt;Record...</i>), "+
            "also writes the macro commands for setting all parameter values to the Recorder window.</dd>"+
            "<dt>Read parameters from log file...</dt><dd>Reads all parameters from a user-selectable file, usually "+
            "a <tt>_log.txt</tt> file saved in a previous spot-tracking session. (Experts might also modify the "+
            "machine-readable part of such a log file to create a setup file; better use a macro for this.)</dd>"+
            "<dt>Compare parameters with log file...</dt><dd>Reads all numeric parameters from a user-selectable file, "+
            "usually a <tt>_log.txt</tt> file saved in a previous spot-tracking session, and compares it with "+
            "the parameters of the current session (shows a table of differences).</dd>"+
            "<dt>Metadata keywords...</dt><dd>Defines the names (case-sensitive) of data columns for energy, "+
            "I0, ... in the input files.  Separate multiple possibilities with a vertical bar. The keywords should be the full prefix, "+
            "as it appears below the window title in the 'SpotTracking' stack, including the '=' sign. Example: "+
            "'<tt>energy=|E=|V=</tt>' allows 'E=123.5' but not 'Energy=123.5'.</dd>"+
            "<dt>Reset all options &amp; close</dt><dd>Resets all options and parameters to their default settings "+
            "and closes the Spot Tracker.</dd>"+
            "</dl>";
    private static final String HELP_GOOD_TO_KNOW =
            "<h2><a name='spottrackerHints'>Good to Know</a></h2>"+
            "<p>On the main ViPErLEED Spot Tracker panel, red color usually indicates where user interaction is "+
            "required next.</p>"+
            "<p>The <tt>Enter</tt> key brings the main ImageJ panel to the foreground. When the Spot Tracker panel is "+
            "in the foreground, <tt>SHIFT-Enter</tt> brings the 'SpotTracking' image stack to the foreground.</p>"+
            "<p>In ImageJ, <i>slices</i> are the single images that an <i>image stack</i> (e.g., a movie) consists of. "+
            "Keyboard shortcuts for moving to the previous/next slice are ',' or '&lt;', and '.' or '&gt;', respectively. "+
            "You can also use the mouse wheel. "+
            "Right-click on the animation (play) symbol at the bottom left to change the speed of playing a movie "+
            "(it may be slower than the speed in frames per second, fps, if calculating the frames takes time).</p>"+
            "<p>If you have a movie with strong brightness variations, you may try setting 'Auto contrast stacks' "+
            "in <i>Edit&gt;Options&gt;Appearance...</i>&nbsp; Bright spots will appear saturated, but you may better see the faint ones.</p>"+
            "<p>Why is there a <tt>(V)</tt> at the end of title of many image stacks? The <tt>(V)</tt> stands for "+
            "<i>virtual stack</i>, i.e., the images are not in memory but read from disk or calculated on the fly if "+
            "required. The stack named <tt>..._SpotTracking</tt> (created by the ViPErLEED Spot Tracker) is a virtual "+
            "stack with caching; its slices are kept in memory as long as there is sufficient memory. "+
            "You can click on the status line at the bottom of the main ImageJ toolbar to see the used and total memory. "+
            "Use <i>Edit&gt;Options&gt;Memory&amp;Threads...<i> to change the maximum amount of "+
            "memory allocated to ImageJ (requires restarting ImageJ thereafter; "+
            "make sure you allocate less than the RAM available in your computer).</p>"+
            "<p>If you get an error message saying that all memory has been used, it may help to add"+
            "<pre>   -XX:SoftRefLRUPolicyMSPerMB=100</pre> "+
            "to the options for the Java Virtual Machine (JVM). These options are in the ImageJ.cfg file in the ImageJ folder. "+
            "On MacOS, these options are in the info.plist inside the ImageJ package: "+
            "Right-click on ImageJ in the Finder and select &quot;Show package contents&quot;. "+
            "The same command with a value of 10 instead of 100 would be even stronger, telling Java "+
            "to be very restrictive with memory spent for caches. "+
            "The downside is that this will usually lead to longer execution times.</p>"+
            "<p>ImageJ can save images (also image stacks and plots) as .zip files. These are compressed .tif files "+
            "that can be directly opened in ImageJ without unpacking (also by drag&amp;drop onto the status line of the main ImageJ panel). "+
            "I(V) movies recorded as ViPErLEED .zip archives can NOT be opened by the ImageJ <i>File&gt;Open</i> command "+
            "(or drag&amp;drop  on the status line); use Open LEED Movie or "+
            "<a href='#openImagesMovies'>More>>Open images/movies</a> for these.<br>"+
            "The file types .tif and .zip are the preferred formats for saving, since they contain the full information "+
            "and image resolution. The mask can be also saved as .png file. "+
            "NEVER save images as .jpg unless they are meant for illustration only; even then mind that .jpg "+
            "has usually very poor quality for sharp borders, as lines, text, etc.</p>"+
            "<p>Single plots (but not stacks of plots) can be modified, e.g., one can zoom into a rectangle selection "+
            "with the '+' key, and reset the range with the gray 'R' appearing under the mouse pointer in the bottom left-corner. "+
            "With the mouse over the left or lower border, there are also other gray symbols; "+
            "their meaning is explained in the help text appearing in the status line of the plot "+
            "when the mouse cursor is above such a symbol. <br>"+
            "You can set the default size of plots (and final size of plot stacks) with <i>Edit&gt;Options&gt;Plots...<i></p>"+
            "<p>Many functions of the "+LEED_Spot_Tracker.PLUGIN_NAME+" can be controlled and automated via "+
            "the ImageJ macro language. "+
            "Record macro commands with <i>Plugins&gt;Macros&gt;Record...</i></p>"+
            "<p>For smoothing in ViPErLEED, the strength of smoothing is given as the inverse square of the noise gain for white noise. "+
            "E.g., '9 points' means that smoothing will reduce white noise to approximately 1/3 of its original value. "+
            "This corresponds to the noise suppression of a moving average over 9 points. "+
            "The actual smoothing algorithm averages over more points with a weight function. "+
            "For I0 data and <i>I</i>(<i>V</i>) curves a modified sinc smoother [<a href='#msSmooth'>2</a>] is used.</p>"+
            "<p>You can copy the contents of any help window like this one ("+LeedUtils.CTRL+"-a to select all and "+LeedUtils.CTRL+"-c). "+
            "Paste it into a text editor for printing, annotating, etc.</p>";
    private static final String HELP_LICENSE =
            "<h2><a name='spottrackerLicense'>License</a></h2>"+
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
            "(use an unzip utility to view its contents).</p>";
    private static final String HELP_REFERENCES_HEADING =
            "<h2>References</h2>";
    private static final String HELP_FOOTNOTE1 =
            "<p><a name='paper'>[1]</a> M. Schmid, F. Kraushofer, A. M. Imre, T. Kißlinger, L. Hammer, U. Diebold, and M. Riva, "+
            "<i>ViPErLEED package II: Spot tracking, extraction and processing of I(V) curves</i>, "+
            "Phys. Rev. Research, 2024. <a href='https://arxiv.org/abs/2406.18413/'>arXiv:2406.18413</a>.</p>";
    private static final String HELP_FOOTNOTE2 =
            "<p><a name='msSmooth'>[2]</a> M. Schmid, D. Rath, and U. Diebold, "+
            "<i>Why and how Savitzky-Golay Filters should be replaced</i>, "+
            "<a href='https://doi.org/10.1021/acsmeasuresciau.1c00054'>ACS Measurement Science Au <b>2</b>, 185 (2022)</a>.</p>";
    private static final String ENDHTML =
            "</html>";

    /** Returns the help string for the main SpotTracker panel */
    public static final String getTrackerHelp() {
        return  PREHTML + HELP_TITLE_TRACKER + HELP_INTRO_TRACKER + HELP_DARKFLATEXPLAIN +
                HELP_INTRO_TRACKER2 + HELP_DARKFLATMENU + HELP_SETENERGIES + HELP_PATTERNFILE +
                HELP_SETINDICES1 + HELP_SETINDICES2 + HELP_SETINDICES3 +
                HELP_SETRADIUS + HELP_TRACKSPOTS1 + HELP_TRACKSPOTS2 + HELP_TRACKSPOTS3 + HELP_SAVEDATA +
                HELP_MOREMENU + HELP_GOOD_TO_KNOW + HELP_LICENSE +
                HELP_REFERENCES_HEADING +HELP_FOOTNOTE1 + HELP_FOOTNOTE2 + ENDHTML;
    }

    /** Returns the help string for the Dark & Flat processing dialog */
    public static final String getDarkFlatProcessingHelp() {
        String str = HELP_DARKFLATMENU +
                "<h2>The Basics of using a Dark Frame &amp; Flat Field</h2><p>" +
                HELP_DARKFLATEXPLAIN + "<p>";
        return makeOutput(str, null, true);
    }

    /** Returns the help string for the "Set Energy, I0, t..."  dialog */
    public static final String getEnergyI0Help() {
        String str = HELP_SETENERGIES;
        return makeOutput(str, HELP_FOOTNOTE2, true);
    }

    /** Returns the help string for the "Select Stack Slice" dialog */
    public static final String getSetIndices1Help() {
        String str = HELP_SETINDICES1;
        return makeOutput(str, null, false);
    }

    /** Returns the help string for the "Select and Label Spots" dialog */
    public static final String getSetIndices2Help() {
        String str = "<h1>Select and Label Spots</h1>" +
                "<p>The " +
                HELP_SETINDICES3;
        return makeOutput(str, null, false);
    }

    /** Returns the help string for the "Set Integration Radius" dialog */
    public static final String getSetRadiusHelp() {
        String str = HELP_SETRADIUS;
        return makeOutput(str, HELP_FOOTNOTE2, true);
    }

    /** Returns the help string for the "Track Spots" dialog */
    public static final String getTrackSpotsHelp() {
        String str = "<h1>Spot Tracking Options</h1>" +
                HELP_TRACKSPOTS2;
        return makeOutput(str, null, true);
    }

    /** Returns the help string for the "Save Data" dialog */
    public static final String getSaveDataHelp() {
        String str = HELP_SAVEDATA;
        return makeOutput(str, null, true);
    }

    /** Creates help dialog of menus; optionally with a footnote.
     *  For modal dialogs, includes a note that the window has to be closed to proceed. */
    private static final String makeOutput(String str, String footnotes, boolean modal) {
        str = str.replace("h2>", "h1>");        //convert subheading of main help to main heading of subtopic help
        str = str.replace("href='#", "name='"); //hides internal links; they would point outside
        str += "<br><hr>";
        if (footnotes != null)
            str += HELP_REFERENCES_HEADING + HELP_FOOTNOTE1 + footnotes;
        if (modal) str = "<small>[Close this help window to proceed]</small>" + str;
        return PREHTML + str + ENDHTML;
    }
}
