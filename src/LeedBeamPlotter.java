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
 *  This class plots one or more curves for a single beam.
 *  It implements the Spot tracker's More>>Plot... command
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

public class LeedBeamPlotter {
	static final String[] SCREEN_FITTER_NAMES = new String[] {"x_calc (screen)", "y_calc (screen)", "x_calc (vs 0,0)", "y_calc (vs 0,0)", };
	static final int NCOLUMNS = 2; 		//columns in radioButtons & checkboxGroup
	static final String SPACER = "---";	//for padding energiesEtc to two columns
	static final String[] COLUMN_HEADINGS = new String[] {"Y Axis", ""}; //for y axis checkboxgroup, length must be NCOLUMNS
	static final String defaultMessage = "* separate beams by ';', e.g. '1/2,0; 2|2'";
	static String spotsString = "1,0";    //of the spots we plot
	static String xAxisVariable = "energy";
	static String[] yVariables;
    /**
     *  This main function displays the dialog and creates the plot
     */
    public static void doDialogAndPlot(final LeedSpotPattern spotPattern, final LeedIVAnalyzer ivAnalyzer, LeedScreenFitter screenFitter, double[][] energiesEtc) {
		LeedIntegerArray energyEtcColumnsA = new LeedIntegerArray(energiesEtc.length);
		for (int i=0; i<energiesEtc.length; i++) {
			if (energiesEtc[i] != null)
				energyEtcColumnsA.add(i);
		}
		int[] energyEtcColumns = energyEtcColumnsA.toArray();  //which columns are non-null in energiesEtc
		String[] energyEtcNames = new String[energyEtcColumns.length];
		for (int i=0; i<energyEtcColumns.length; i++)
			energyEtcNames[i] = LEED_Spot_Tracker.ENERGY_ETC_NAMES[energyEtcColumns[i]];

		final String[] variableNames = concatenate(energyEtcNames, SCREEN_FITTER_NAMES, LeedSpotAnalyzer.DATA_NAMES, LeedIVAnalyzer.DATA_NAMES);
		final int nVariables = variableNames.length;
		int nRows = (variableNames.length+NCOLUMNS-1)/NCOLUMNS;
		boolean[] yAxisSelected = new boolean[nVariables];
		if (yVariables != null) {
			for (String yVariable : yVariables) {
				int index = LeedUtils.arrayIndexOf(variableNames, yVariable);
				if (index > 0)
					yAxisSelected[index] = true;
			}
		}

		GenericDialog gd = new GenericDialog("Plot for Beam(s)");
		gd.addStringField("Beam(s) *", spotsString, 12);
		gd.setInsets(5, 0, 0); //top left bottom
		gd.addRadioButtonGroup("X Axis", variableNames, nRows, NCOLUMNS, xAxisVariable);
		Label xAxisLabel = (Label)gd.getMessage();
		gd.addCheckboxGroup(nRows, NCOLUMNS, variableNames, yAxisSelected, COLUMN_HEADINGS);
		gd.setInsets(5, 0, 0); //top left bottom
		gd.addMessage(defaultMessage);  //becomes error message
        final Label errorLabel = (Label)gd.getMessage();
        if (xAxisLabel != null) {
			Font font = xAxisLabel.getFont();
			Font boldFont = font == null ? new Font("Dialog", Font.BOLD, 12) : font.deriveFont(Font.BOLD);
			xAxisLabel.setFont(boldFont);
        }
        disableSpacerCheckboxes(gd);
        

		DialogListener dialogListener = new DialogListener() {
            public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
				String spotsString = gd.getNextString();
				String xAxisVariableName = gd.getNextRadioButton();
				boolean[] yAxisSelected = getCheckboxStates(gd, nVariables);

                String[] spotStrs = Tools.split(spotsString, ";");
                int[] spotIndices = spotPattern.getIndices(spotStrs);
                String error = spotPattern.getErrorMessage();
                if (error != null) {
                    showError(errorLabel, error);
                    return false;
                }
                int xIndex = LeedUtils.arrayIndexOf(variableNames, xAxisVariableName);
                if (xIndex < 0) {
                    showError(errorLabel, "No x axis selected");
                    return false;
                }
				int[] yIndices = indicesOf(yAxisSelected);
				if (yIndices.length == 0) {
                    showError(errorLabel, "No y axis selected");
                    return false;
				}
				boolean hasIVanalyzerData = getIVanalyzerIndex(xIndex) >= 0;
				for (int yIndex : yIndices)
					if (getIVanalyzerIndex(yIndex) >= 0)
						hasIVanalyzerData = true;
                if (hasIVanalyzerData)
					for (int spot : spotIndices)
						if (!ivAnalyzer.hasSpot(spot)) {
							showError(errorLabel, "Spot not detected: "+spotPattern.getName(spot));
							return false;
						}
                errorLabel.setText(defaultMessage);
                errorLabel.setForeground(Color.BLACK);
				return true;
            }
        };
        Button okButton = gd.getButtons()[0];
        okButton.setEnabled(dialogListener.dialogItemChanged(gd, null)); //check default input
        gd.addDialogListener(dialogListener);
        gd.showDialog();

        if (gd.wasCanceled()) return;
        spotsString = gd.getNextString();
		xAxisVariable = gd.getNextRadioButton();
		yAxisSelected = getCheckboxStates(gd, nVariables);

		String[] spotStrs = Tools.split(spotsString, ";");
		int[] spotIndices = spotPattern.getIndices(spotStrs);
		int[] yIndices = indicesOf(yAxisSelected);
		yVariables = new String[yIndices.length];
		for (int i=0; i<yIndices.length; i++)
			yVariables[i] = variableNames[yIndices[i]];
		int nCurves = spotIndices.length * yIndices.length;

		Color[] colors = LeedPlotter.getColors(nCurves);
		String title = "Spot Plot "+LeedUtils.getDateFormatted(" [HH'h'mm]");
		String yAxisLabel = yIndices.length == 1 ? variableNames[yIndices[0]] : "value";
		Plot plot = new Plot(title, xAxisVariable, yAxisLabel);
		int iColor = 0;
		String legendLabels = "";
		for (int i=0; i<spotIndices.length; i++) {
			int spotIndex = spotIndices[i];
			double[] xData = getVariableArray(spotPattern, ivAnalyzer, screenFitter,
					energyEtcColumns, energiesEtc, xAxisVariable, spotIndex);
			for (String yVariable : yVariables) {
				double[] yData = getVariableArray(spotPattern, ivAnalyzer, screenFitter,
						energyEtcColumns, energiesEtc, yVariable, spotIndex);
				boolean isSpotDependent = LeedUtils.arrayIndexOf(LEED_Spot_Tracker.ENERGY_ETC_NAMES, yVariable) <0;
				if (i==0 || isSpotDependent) {		//(plot spot-independent data only once)
					plot.setColor(colors[iColor++]);
					plot.addPoints(xData, yData, Plot.LINE);
					if (isSpotDependent)
						legendLabels += spotPattern.getName(spotIndex)+" ";
					legendLabels += yVariable+"\n";
				}
			}
		}
		plot.setColor(Color.BLACK);
		plot.addLegend(legendLabels);
		plot.show();
		plot.setLimitsToFit(true);
		LeedIVAnalyzer.addToPlotImpList(plot.getImagePlus(), 'm');
	}

	static void showError(Label errorLabel, String errorText) {
		errorLabel.setText("ERROR "+errorText);
		errorLabel.setForeground(Color.RED);
	}

	/** Returns the array containing the variable values for a given variable index and spot */
	static double[] getVariableArray(LeedSpotPattern spotPattern, LeedIVAnalyzer ivAnalyzer, LeedScreenFitter screenFitter,
			int[] energyEtcColumns, double[][] energiesEtc, String variableName, int spotIndex) {
		int ivAnalyzerIndex = LeedUtils.arrayIndexOf(LeedSpotAnalyzer.DATA_NAMES, variableName); //one of SpotAnalyzer names?
		if (ivAnalyzerIndex < 0) {
			ivAnalyzerIndex = LeedUtils.arrayIndexOf(LeedIVAnalyzer.DATA_NAMES, variableName);
			if (ivAnalyzerIndex >= 0)
				ivAnalyzerIndex += LeedSpotAnalyzer.DATA_NAMES.length;	//these are added columns of the ivAnalyzer: "dx_raw" etc.
		}
		if (ivAnalyzerIndex >= 0)
			return ivAnalyzer.getData(ivAnalyzerIndex, spotIndex);
		int variableIndex = LeedUtils.arrayIndexOf(LEED_Spot_Tracker.ENERGY_ETC_NAMES, variableName); //energy, I0 etc.?
		if (variableIndex == LEED_Spot_Tracker.ENERGY &&		//both ENERGY & E_LEEM have name "energy"
				energiesEtc[LEED_Spot_Tracker.ENERGY] == null && energiesEtc[LEED_Spot_Tracker.E_LEEM] != null)
			variableIndex = LEED_Spot_Tracker.E_LEEM;
		if (variableIndex >= 0)
			return energiesEtc[variableIndex];
		//otherwise ScreenFitter x_calc, y_calc (screen or vs 00)
		variableIndex = LeedUtils.arrayIndexOf(SCREEN_FITTER_NAMES, variableName);
		if (variableIndex < 0)
			throw new RuntimeException("Internal error: No match for data named "+variableName);
		double[] energies = energiesEtc[LEED_Spot_Tracker.ENERGY];
		double[] out = new double[getNonNullSubarrayLength(energiesEtc)];
		double[] xy = new double[2];
		int direction = variableIndex % 2;		//0=x, 1=y
		for (int i=0; i<out.length; i++) {
			double energy = energies != null ? energies[i] : LeedRadiusSelector.UNKNOWN_ENERGY;
			xy = screenFitter.screenCoordinates(spotPattern, spotIndex, 1./Math.sqrt(energy), xy);
			out[i] = xy[direction];
			if (variableIndex >= 2)
				out[i] -= direction==0 ? screenFitter.getX00() :  screenFitter.getY00();
		}
		return out;
	}

	/** Returns the length of the first non-null array contained in a[][],
	 *  or -1 if no non-null array */
	static int getNonNullSubarrayLength(double[][] a) {
		for (int i=0; i<a.length; i++)
			if (a[i] != null)
				return a[i].length;
		return -1;
	}

	/** For a given index of a variable, returns the index in the IVanalyzer's data array,
	 *  or -1 if not a variable in the IVanalyzer's data */
	static int getIVanalyzerIndex(int variableIndex) {
		int index = variableIndex - SCREEN_FITTER_NAMES.length; //first variables are data from screen fitter
		if (index < 0 || index >= LeedSpotAnalyzer.DATA_NAMES.length + LeedIVAnalyzer.DATA_NAMES.length)
			return -1;
		return index;
	}

	/** Returns the state of the checkboxes as boolean array */
	static boolean[] getCheckboxStates(GenericDialog gd, int nCheckboxes) {
		boolean[] out = new boolean[nCheckboxes];
		for (int i=0; i<nCheckboxes; i++)
			out[i] = gd.getNextBoolean();
		return out;
	}

	/** Disables all checkboxes of this container and its children */
	static void disableSpacerCheckboxes(Container container) {
		for (Component component : container.getComponents()) {
			if ((component instanceof Checkbox) && SPACER.equals(((Checkbox)component).getLabel()))
				((Checkbox)component).setEnabled(false);
			else if (component instanceof Container)
				disableSpacerCheckboxes((Container)component);
		}
	}

	/** Returns an int array with the indices of the boolean elements that are true */
	static int[] indicesOf(boolean[] bArray) {
		LeedIntegerArray iArray = new LeedIntegerArray(5);
		for (int i=0; i<bArray.length; i++)
			if (bArray[i]) iArray.add(i);
		return iArray.toArray();
	}

	/** Concatenates String arrays. The first two must be non-null.
	 *  If the length of the first one is not an integer multiple of NCOLUMNS, itis padded by spacers. */
	static String[] concatenate(String[] a, String[] b, String[] c, String[] d) {
		int aLen = a.length;
		int aLenPadded = ((a.length + NCOLUMNS-1)/NCOLUMNS)*NCOLUMNS;
		int bLen = b.length;
		int cLen = c== null ? 0 : c.length;
		int dLen = d== null ? 0 : d.length;
		String[] out = new String[aLenPadded + bLen + cLen + dLen];
		System.arraycopy(a, 0, out, 0, aLen);
		for (int i=aLen; i<aLenPadded; i++)
			out[i] = SPACER;
		System.arraycopy(b, 0, out, aLenPadded, bLen);
		if (c != null) System.arraycopy(c, 0, out, aLenPadded+bLen, cLen);
		if (d != null) System.arraycopy(d, 0, out, aLenPadded+bLen+cLen, dLen);
		return out;
	}
}
