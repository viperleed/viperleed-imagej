import ij.plugin.frame.PlugInFrame;
import ij.*;
import ij.gui.*;
import java.awt.*;
import java.awt.event.*;
import java.util.Hashtable;
import java.util.HashMap;


/** This ImageJ plugin is a dispatcher for the ViPErLEED commands.
 *  It is inspired by the ImageJ 'Recent Commands' class
 *  (Plugins>Utilities>Commands...; ij.plugin.frame.Commands).
 */
 
 /** This ImageJ plugin is part of the ViPErLEED package for LEED I(V) analysis.
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

public class LEED_Commands extends PlugInFrame implements MouseListener {
    public static final String LOC_KEY = "leedcommands.loc";	//ImageJ remembers position and width of the panel with these keys
    private static Frame instance;								//This frame is a singleton (one one instance)
    private static final String help = "Help";
	//Commands as they appear in the list, and the corresponding ImageJ class call
    private static final String[] commandsAndClasses = {
        "Open LEED Movie", "Open_LEED_Movie",
        "Spot Tracker", "LEED_Spot_Tracker",
        "I(V) Curve Editor", "LEED_Curve_Editor",
        "----- I(V) Curve Utilities -----", null,
		"  Average I(V) Curves", "LEED_Average_Datasets",
		"  Stitch I(V) Curves", "LEED_Stitch_Datasets",
		"  I(V) Curve Interpolation", "LEED_IV_Curve_Interpolation",
		"  R Factor Between Data Sets", "LEED_R_Factor_Between_Datasets",
		"  I(V) Curve Quality Statistics", "LEED_Data_Quality_Statistics",
		"  I(V) Curve Tools", "LEED_IV_Curve_Tools",
        help, null
    };
    private Color labelBackground, highlightBackground;
    private HashMap<String, String> commands = new HashMap<String, String>(commandsAndClasses.length/2);
    private static Dialog helpDialog;
    private static final String HELP_STRING =
            "<html>"+
            "<h1>ViPErLEED Commands Overview</h1>"+
            "<p>The <b>ViPErLEED Commands</b> window provides a quick way to access the commands. "+
            "These commands are also available via the ImageJ menu Plugins>ViPErLEED or the ImageJ command finder "+
            "(Plugins>Utilities>Find Commands, or type "+LeedUtils.CTRL+"-L)</p>"+
            "<ul>"+
            "<li><b>Open LEED Movie</b> opens a ViPErLEED movie stored as a .zip file or an 'AIDA' LEED movie "+
            "(Automatic Image and Data Acquisition, EE2000/EE2010; <tt>*.vid</tt> format). LEED movies are shown as "+
            "<em>Image Stacks</em> in ImageJ.</li>"+
            "<li>The <b>Spot Tracker</b> is used to extract <i>I</i>(<i>V</i>) curves from a LEED movie.</li>"+
            "<li>The <b>I(V) Curve Editor</b> is used for selection and smoothing of experimental <i>I</i>(<i>V</i>) curves.</li>"+
            "<li><b>Average I(V) Curves</b> averages <i>I</i>(<i>V</i>) data.<br>"+
            "If you have measured several LEED movies (at different distances to the sample), "+
            "use this tool for averaging the <i>I</i>(<i>V</i>) curves, to reduce the noise (including the noise introduced by the grids). "+
            "Thereafter, use the I(V) Curve Editor.</li>"+
            "<li><b>Stitch I(V) Curves</b> joins sets of <i>I</i>(<i>V</i>) curves that have different energy ranges (with a small overlap) "+
            "and possibly different scale factors (e.g., camera gain, exposure time, etc.).</li>"+
            "<li><b>I(V) Curve Interpolation</b> interpolates data to different (usually finer) energy steps.</li>"+
            "<li><b>R Factor Between Data Sets</b> compares two files with <i>I</i>(<i>V</i>) curves and lists the <i>R</i> factor between them "+
            "for each beam and the overall <i>R</i> factor. "+
            "(Values obtained with <tt>viperleed.calc</tt> may be slightly different, due to subtle differences of the algorithms.)</li>"+
            "<li><b>I(V) Curve Quality Statistics</b> provides statistics for assessing the quality of <i>I</i>(<i>V</i>) measurements, "+
            "such as <i>R</i> factors between symmetry-equivalent beams and negative intensities. "+
            "This tool can create a list or a plot; the latter is the same as obtained from the Spot Tracker.</li>"+
            "<li><b>I(V) Curve Tools</b> can modify the intensities in <i>I</i>(<i>V</i>) curve files, "+
            "extract a subset of the beams or an energy range, transform beam indices, "+
            "and/or modify the information on symmetry-equivalent beams in <i>I</i>(<i>V</i>) curve files.</li>"+
            "</ul>"+
            "<h2>Good to know</h2>"+
            "<p>If you want to have the <i>LEED Commands</i> panel present when ImageJ starts, add "+
            "<tt><font color='#cc0000'>run(&quot;LEED Commands&quot;);</font></tt> "+
            "(with the semicolon) to the &quot;AutoRun&quot; macro of your StartupMacros file (In the ImageJ/macros directory) or "+
            "to the list of startup commands by selecting Edit&gt;Options&gt;Startup from the ImageJ Menu.</p>"+
            "<p>You can copy the contents of any help window like this one ("+LeedUtils.CTRL+"-a to select all and "+LeedUtils.CTRL+"-c). "+
            "Paste it into a text editor for printing, annotating, etc.</p>"+
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
            "<p><a name='paper'>[1]</a> M. Schmid, F. Kraushofer, A. M. Imre, T. Kißlinger, L. Hammer, U. Diebold, and M. Riva, "+
            "<i>ViPErLEED package II: Spot tracking, extraction and processing of I(V) curves</i>, "+
            "Phys. Rev. Research, 2024. <a href='https://arxiv.org/abs/2406.18413/'>arXiv:2406.18413</a></p>"+
            "</html>";

    /** The constructor creates a "LEED Commands" panel, or brings the current one (when existing) to the foreground. */
    public LEED_Commands() {
        super("ViPErLEED");
        if (instance!=null) {
            WindowManager.toFront(instance);
            return;
        }
        instance = this;
        WindowManager.addWindow(this);
        /* populate HashMap of Commands*/
        Hashtable ijCommands = Menus.getCommands();
        for (Object ijCommand : ijCommands.keySet()) {
			String ijClassPath = (String)ijCommands.get(ijCommand);
			for (int i=0; i<commandsAndClasses.length; i+=2) {
				if (commandsAndClasses[i+1]!= null && ijClassPath.endsWith(commandsAndClasses[i+1]))
					commands.put(commandsAndClasses[i], (String)ijCommand);
			}
		}
		/* set up colors */
		Color foreground = getForeground();
		Color background = getBackground();
		if (foreground == null && background == null) {
			background = new Color(0xdddddd);
			foreground = Color.BLACK;
		} else if (foreground == null) {
			foreground = brightness(background) > 128 ? Color.BLACK : Color.WHITE;
		} else if (background == null) {
			background = brightness(foreground) > 128 ? Color.BLACK : new Color(0xdddddd);
		}
		boolean lightBackground = brightness(foreground) < brightness(background);
		labelBackground = background;
		highlightBackground = lightBackground ? brighter(labelBackground) : labelBackground.darker();
		if (difference(labelBackground,highlightBackground) < 200)
			highlightBackground = lightBackground ? labelBackground.darker() : brighter(labelBackground);
		/* set up gui */
        GridLayout layout = new GridLayout(commandsAndClasses.length/2, 1);
        setLayout(layout);
        ImageJ ij = IJ.getInstance();
        addKeyListener(ij);
        for (int i=0; i<commandsAndClasses.length; i+=2) {
			boolean isActive = commandsAndClasses[i+1] != null || commandsAndClasses[i].equals(help);
			Label label = new Label(commandsAndClasses[i]);
			label.addMouseListener(this);
			label.addKeyListener(ij);
			label.setBackground(labelBackground);
			//if (!isActive) label.setForeground(Color.BLUE);
			GUI.scale(label);
			add(label);
		}
        pack();
        Point loc = Prefs.getLocation(LOC_KEY);
        if (loc!=null)
            setLocation(loc);
        show();
    }

	/** Shows the help text */
    public void showHelp() {
        if (helpDialog != null && helpDialog.isVisible())
            helpDialog.toFront();
        else
            helpDialog = new HTMLDialog("ViPErLEED Commands", HELP_STRING, false); // non blocking
    }

    /** Overrides PlugInFrame.close(). */
    public void close() {
        super.close();
        instance = null;
        Prefs.saveLocation(LOC_KEY, getLocation());
    }

    /** The MouseListener callback executes the command */
    public void mousePressed(MouseEvent e) {
        Label label = (Label)e.getSource();
		String itemStr = label.getText();
		if (itemStr.equals(help))
			showHelp();
		else {
			String command = commands.get(itemStr);
			if (command != null)
				IJ.doCommand(command);
		}
    }

    public void mouseReleased(MouseEvent e) {}
    public void mouseClicked(MouseEvent e) {}
    public void mouseEntered(MouseEvent e) {
		highlight(e, true);
	}
    public void mouseExited(MouseEvent e) {
		highlight(e, false);
	}

	/** Highlights the label or de-highlights it (when the mouse enters or exits) */
	void highlight(MouseEvent e, boolean high) {
		Label label = (Label)e.getSource();
		String itemStr = label.getText();
		if (itemStr.equals(help) || commands.get(itemStr) != null)
			label.setBackground(high ? highlightBackground : labelBackground);
	}

	/** Brightness of a color in 0-255 range */
	double brightness(Color c) {
		return 0.299*c.getRed() + 0.587*c.getGreen() + 0.114*c.getBlue();
	}

	/** Rough indication of the perceptive difference between two colors */
	double difference(Color c1, Color c2) {
		double diff = 0.299*sqr(c2.getRed() - c1.getRed())
				+ 0.587*sqr(c2.getGreen() - c1.getGreen())
				+ 0.114*sqr(c2.getBlue() - c1.getBlue());
		double sumBright = brightness(c1) + brightness(c2);
		diff *= Math.min(0.05*sumBright, 1);		//very bright or dark colors may be saturated; reduce difference value
		diff *= Math.min(0.05*(520-sumBright), 1);
		return diff;
	}

	/** A brighter color since Color.brighter() leaves very dark colors almost unchanged */
	Color brighter(Color c) {
		Color brighter = c.brighter();
		if (difference(c, brighter) > 500)
			return brighter;
		else
			return new Color(Math.min(c.getRed()+64, 255), Math.min(c.getGreen()+64, 255), Math.min(c.getBlue()+64, 255));
	}

    static double sqr(double x) {return x*x;}
}
