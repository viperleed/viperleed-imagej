import ij.*;
import ij.IJ;
import ij.gui.*;
import ij.process.*;
import ij.measure.Measurements;
import ij.plugin.filter.MaximumFinder;
import ij.util.Tools;
import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * This class creates most of the plots of the Spot Tracker
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

public class LeedPlotter {
    static final String X_LABEL = "energy";
    static final double MID_SIGNIFICANCE=10;
    double[][][] data;  //[type][spot][energyIndex]
    double[] xAxis;
    boolean[] showSpot;
    LeedSpotPattern spotPattern;
    String xAxisLabel;
    Plot plot;
    String title;
    static final int ENERGY=0, INTEGER=1, SUPERSTRUCTURE=2; //indices for arrays in integrationRadiusData
    float[][] integrationRadiusData;                        //float arrays for plotting integration-halfradius

    /** Creator
     *  'limits' should be x_left, x_right, y_bottom, y_top.
     *  'limits' may be null for auto or contain NaN values for auto */
    public LeedPlotter(String title, double[][][] data, double[] xAxis, boolean[] showSpot, LeedSpotPattern spotPattern, String xAxisLabel) {
        this.title=title;
        this.data = data;
        this.xAxis = xAxis;
        this.showSpot = showSpot;
        this.spotPattern = spotPattern;
        this.xAxisLabel = xAxisLabel;
    }

    /** Creates a plot stack with one image for each spot that has showSpot true.
     *  Each plot can show one or more types of data from the data[][][] array, as given by the columns variable.
     *  When columnFactors is not null, the corresponding factor is applied to the respective column
     *  Use show() to show the plot stack */
    public Plot plotStack(int[] columns, String[] columnLabels, double[] columnFactors, double[] limits) {
        limits = getLimits(limits);
        plot = new Plot(title, xAxisLabel, title);
        plot.setLimits(limits[0], limits[1], limits[2], limits[3]);
        Color[] colors = getColors(columns.length);
        for (int spot=0; spot<data[columns[0]].length; spot++) {
            if (showSpot[spot]) {
                for (int c=0; c<columns.length; c++) {
                    int column = columns[c];
                    plot.setColor(colors[c]);
                    double[] yValues = data[column][spot];
                    if (columnFactors != null && !Double.isNaN(columnFactors[c])) {
                        double[] yNew = new double[yValues.length];
                        for (int i=0; i<yValues.length; i++)
                            yNew[i] = yValues[i]*columnFactors[c];
                        yValues = yNew;
                    }
                    plot.addPoints(xAxis, yValues, Plot.LINE);
                    if (columnLabels[c] != null)
                        plot.setLabel(c, columnLabels[c]);
                }
                plot.setColor(Color.BLACK);
                plot.setFont(Font.BOLD, 18);
                plot.setJustification(Plot.RIGHT);
                String label = spotPattern.getName(spot)+" ["+spotPattern.getGroup(spot)+"]";
                plot.addLabel(0.93, 0.08, label);
                plot.setFont(Font.PLAIN, 14);
                if (columns.length > 1)
                    plot.setLegend(null, Plot.AUTO_POSITION);
                plot.draw();    //workaround to avoid wrong range for 1st
                if (columns.length > 1)
                    plot.setLegend(null, Plot.AUTO_POSITION);
                plot.setLimits(limits);
                plot.addToStack();
            }
        }
        return plot;
    }

    /** Plots the radius as function of energy for all spots, with larger symbols more more significant spots.
     *  If xIsEnergy is true, also plots half of the spot integration radius for comparison */
    public Plot plotRadii(double[] radiiInt, double[] radiiSup, boolean xIsEnergy, int spotBackgrShape) {
        plot = new Plot(title, xAxisLabel, "spot radius");
        double xStep = LeedUtils.getSmallestStep(xAxis);

        //array indices: 0, 2 superstructure; 1, 3 integer; 0, 1 low-significance, 2, 3 high-significance
        //this is the sequence of plotting (more important on top of less important)
        int[] nData = new int[4];   //number of super/integ spot points with [0,1] significance <= MID_SIGNIFICANCE and [1] above
        for (int spot=0; spot<data[0].length; spot++) {
            if (showSpot[spot]) {
                boolean isSuperstructure = spotPattern.isSuperstructure(spot);
                for (int i=0; i<xAxis.length; i++) {
                    int arrayIndex = makeArrayIndex(isSuperstructure, data[LeedSpotAnalyzer.PSIGNIFICANCE][spot][i]);
                    if (arrayIndex >=0) {
                        nData[arrayIndex]++;
                    }
                }
            }
        }
        int nSupData = nData[0] + nData[2];
        int nIntData = nData[1] + nData[3];

        if (nSupData + nIntData == 0)                           //nothing to plot at all?
             return null;
        float[][] xData = new float[nData.length][];            //arrays where we transfer all data
        float[][] rSize = new float[nData.length][];
        float[][] tSize = new float[nData.length][];
        for (int j=0; j<nData.length; j++) {
            xData[j] = new float[nData[j]];
            rSize[j] = new float[nData[j]];
            tSize[j] = new float[nData[j]];
        }
        int[] count = new int[nData.length];
        for (int spot=0; spot<data[0].length; spot++) {
            if (showSpot[spot]) {
                boolean isSuperstructure = spotPattern.isSuperstructure(spot);
                double offset = 0;
                if (nSupData>0 && nIntData>0 && xIsEnergy)      //slightly shift x values to avoid plotting integer & superstructure on top of each other
                    offset = isSuperstructure ? -0.17*xStep : 0.17*xStep;
                for (int i=0; i<xAxis.length; i++) {
                    int arrayIndex = makeArrayIndex(isSuperstructure, data[LeedSpotAnalyzer.PSIGNIFICANCE][spot][i]);
                    if (arrayIndex >=0) {
                        xData[arrayIndex][count[arrayIndex]] = (float)(xAxis[i] + offset);
                        rSize[arrayIndex][count[arrayIndex]]  = (float)data[LeedSpotAnalyzer.R_SIZE][spot][i];
                        tSize[arrayIndex][count[arrayIndex]]  = (float)data[LeedSpotAnalyzer.T_SIZE][spot][i];
                        count[arrayIndex]++;
                    }
                }
            }
        }

        boolean manyIntSpots = nIntData/xAxis.length > 50;      //if we have many spots, we use smaller symbols, and weak colors for low significance
        boolean manySupSpots = nSupData/xAxis.length > 50;
        int nSpotTypes = (nIntData>0 && nSupData>0) ? 2 : 1;    //if we have both, type 0=superstructure, 1=integer
        Color[] goodColors = getColors(2*nSpotTypes);           //for each spot type, avg & major radius (different indices than data!)
        Color[] weakColors = new Color[goodColors.length];      //colors for the low-significance spots
        for (int i=0; i<goodColors.length; i++) {
            boolean isInt = nSupData == 0 || i >= 2;            //it's the color for integer spots
            boolean isFaint = isInt ? manyIntSpots : manySupSpots;
            Color c = goodColors[i];
            weakColors[i] = isFaint ?
                    new Color(makeBrighter(c.getRed()), makeBrighter(c.getGreen()), makeBrighter(c.getBlue())) :
                    c;
        }
        //plot the measured spot sizes
        for (int j=0; j<nData.length; j++) {                    //0=sup&lo, 1=sup&hi, 2=int&lo, 3=int&hi significance
            boolean highSignificance = j >= 2;
            if (!highSignificance && nData[j+2] == 0)
                highSignificance = true;                        //if we have no high-significance data, treat low significance as high
            if (nData[j] == 0)
                continue;                                       //no data? skip plotting
            int colorIndexBaseHalf = nSpotTypes == 2 ? j%2 : 0; //1 for integer beams if we have both int & sup, otherwise 0
            int colorIndexBase = 2*colorIndexBaseHalf;          //2 for integer beams if we have both int & sup, otherwise 0
            Color[] colors = highSignificance ? goodColors : weakColors;
            Color tSizeColor = colors[colorIndexBase];
            Color rSizeColor = colors[colorIndexBase+1];
            boolean isInt = nSupData == 0 || colorIndexBase >= 2;
            boolean manySpots = isInt ? manyIntSpots : manySupSpots;
            String intOrSup = "";
            if (nSpotTypes == 2) intOrSup = isInt ? " (integer)" : " (superstructure)";
            plot.setLineWidth(highSignificance && !manySpots ? 2 : 1);
            plot.setColor(tSizeColor);
            String tSizeLabel = highSignificance ? "Spot tangential size"+intOrSup : null;    //only high-significance markers in the legend
            plot.addPoints(xData[j], tSize[j], /*yErrorBars*/null, Plot.DOT, tSizeLabel);
            plot.setColor(rSizeColor);
            String rSizeLabel = highSignificance ? "Spot radial size"+intOrSup : null;
            plot.addPoints(xData[j], rSize[j], /*yErrorBars*/null, Plot.DOT, rSizeLabel);
        }
        //plot half integration radii
        integrationRadiusData = new float[3][];
        integrationRadiusData[ENERGY] = Tools.toFloat(xAxis);    //we use float arrays so we can change them 'on the fly' and the plot changes, in case we ever make it 'live')
        integrationRadiusData[INTEGER] = nIntData > 0 ? new float[xAxis.length] : null;
        integrationRadiusData[SUPERSTRUCTURE] = nSupData > 0 ? new float[xAxis.length] : null;
        calcHalfIntegrationRadius(radiiInt, radiiSup, integrationRadiusData);
        boolean doInteger = integrationRadiusData[INTEGER] != null;
        boolean doSuperstructure = integrationRadiusData[SUPERSTRUCTURE] != null &&
                (!doInteger || !Arrays.equals(radiiInt, radiiSup)); //don't plot both if they are the same

        plot.setLineWidth(2);
        if (doInteger) {
            String label = "Integrat. radius/2";
            if (doSuperstructure) label += " (integer)";
            plot.setColor(Color.BLUE);
            plot.addPoints(integrationRadiusData[ENERGY], integrationRadiusData[INTEGER], /*yErrorBars*/null, Plot.LINE, label);
        }
        if (doSuperstructure) {
            String label = "Integrat. radius/2";
            if (doInteger) label += " (superstructure)";
            plot.setColor(Color.BLACK);
            plot.addPoints(integrationRadiusData[ENERGY], integrationRadiusData[SUPERSTRUCTURE], /*yErrorBars*/null, Plot.LINE, label);
        }

        plot.setColor(Color.BLACK);
        plot.setLineWidth(1);
        plot.setLegend(null, Plot.AUTO_POSITION);
        double[] limits = getLimits(null);
        double radiusMax = 0;
        if (xIsEnergy && (nIntData > 0 || nSupData > 0)) {
            double[] radiusMinMax = Tools.getMinMax(integrationRadiusData[SUPERSTRUCTURE] != null ?
                integrationRadiusData[SUPERSTRUCTURE] : integrationRadiusData[INTEGER]);
            radiusMax = radiusMinMax[1];
        } else
            radiusMax = LeedRadiusSelector.radius(LeedRadiusSelector.UNKNOWN_ENERGY, doSuperstructure);
        plot.setLimits(limits[0], limits[1], 0, 1.5*radiusMax);
        return plot;
    }

    /** Returns -1 if significance is not > 0, otherwise the array index for plotRadii:
     *  0 or 2 superstructure; 1 or 3 integer; 0 or 1 low-significance, 2 or 3 high-significance */
    int makeArrayIndex(boolean isSuperstructure, double significance) {
        if (significance > 0) {
            int index = isSuperstructure ? 0 : 1;
            if (significance > MID_SIGNIFICANCE) index += 2;
            return index;
        } else
            return -1;
    }
        

    public Plot plotRFactorOfEquivalent(double v0i, double eVsmooth, double eVseparate) {
        //create new data and indices arrays only for spots where we have more than one per group
        LeedIntegerArray beamsPerGroup = new LeedIntegerArray();
        for (int i=0; i<data[LeedIVAnalyzer.INT_I0CORRECTED].length; i++) {
            if (showSpot[i]) {
                int group = spotPattern.getGroup(i);
                if (group >= 0)            //ignore symmetry-forbidden groups
                    beamsPerGroup.increment(group);
            }
        }
        ArrayList<double[]> intDataL = new ArrayList<double[]>();
        LeedIntegerArray indicesA = new LeedIntegerArray();
        for (int i=0; i<data[LeedIVAnalyzer.INT_I0CORRECTED].length; i++) {
            if (showSpot[i]) {
                int group = spotPattern.getGroup(i);
                if (group >= 0 && beamsPerGroup.get(spotPattern.getGroup(i)) > 1) {
                    intDataL.add(data[LeedIVAnalyzer.INT_I0CORRECTED][i]);
                    indicesA.add(i);
                }
            }
        }
        int[] indices = indicesA.toArray();
        double[][] intData = intDataL.toArray(new double[0][0]);

        plot = (Plot)LEED_Data_Quality_Statistics.analyze(title, LEED_Data_Quality_Statistics.PLOT, xAxis, indices, intData,
            spotPattern, v0i, eVsmooth, eVseparate);
        return plot;
    }

    /** Gets plot limits as an array; if the x limits (first two numbers) are not given, replaces them by the energy range.
     *  May be also called with a null array */
    double[] getLimits(double[] limits) {
        if (limits == null)
            limits = new double[] {Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] minmax = Tools.getMinMax(xAxis);
        if (Double.isNaN(limits[0]))
            limits[0] = minmax[0];
        if (Double.isNaN(limits[1]))
            limits[1] = minmax[1];
        return limits;
    }

    /** Extracts from the given data column for the given spot the values with the given significance level,
     *  where 0 is high significance and 1 is low significance
     */
    double[] getWithSignificance(int column, int spot, int wantedLevel, double[] output) {
        double[] input = data[column][spot];
        double[] sigData = data[LeedSpotAnalyzer.PSIGNIFICANCE][spot];
        if (output == null) output = new double[input.length];
        for (int i=0; i<input.length; i++) {
            double sig = sigData[i];
            int level = sig>MID_SIGNIFICANCE ? 0 : 1;
            output[i] = (level == wantedLevel) ? input[i] : Double.NaN;
        }
        return output;
    }

	/** Creates an 8-bit color value closer to 255 */
	private static int makeBrighter(int v) {
		return v> 190 ? 255 : 255 - (int)(0.4*(255-v));
	}

    /** Calculates the integration radii/2 */
    void calcHalfIntegrationRadius(double[] radiiInt, double[] radiiSup, float[][] integrationRadiusData) {
        boolean doInteger = integrationRadiusData[INTEGER] != null;
        boolean doSuperstructure = integrationRadiusData[SUPERSTRUCTURE] != null;
        for (int i=0; i<xAxis.length; i++) {
            if (doInteger)
                integrationRadiusData[INTEGER][i] = 0.5f*(float)radiiInt[i];
            if (doSuperstructure)
                integrationRadiusData[SUPERSTRUCTURE][i] = 0.5f*(float)radiiSup[i];
        }
    }

    /** Creates a plot of the 2D regression, i.e., the position of the (0,0) spot and the scale factor.
     *  yLabels should be linefeed-terminated one item per yData.
     *  Also plots a vector for the displacement of the (0,0) due to stray magnetic fields (unless NaN);
     *  the length in plot frame heights is obtained by multiplying with 'vectorYScale'.
     */
    public Plot plot2Dregression(double[][] yData, String yLabels,
            double vectorX, double vectorY, double vectorYscale, String vectorText) {
        plot = new Plot(title, xAxisLabel, "deviation");
        Color[] colors = getColors(yData.length);
        for (int d=0; d<yData.length; d++) {
            plot.setColor(colors[d]);
            plot.addPoints(xAxis, yData[d], Plot.LINE);
        }
        plot.setColor(Color.BLACK);
        if (!Double.isNaN(vectorX)) {           //draw vector in the top-right corner
            Rectangle plotFrame = plot.getDrawingFrame();
            double vectorXscale = vectorYscale*plotFrame.height/(double)plotFrame.width;
            double vecXsc = vectorX*vectorXscale;
            double vecYsc = vectorY*vectorYscale;
            double xCenter = 0.98 - 0.5*Math.abs(vecXsc);
            double yCenter = 0.06 + 0.5*Math.abs(vecYsc);
            double xEnd = xCenter + 0.5*vecXsc;
            double yEnd = yCenter + 0.5*vecYsc;
            plot.setLineWidth(2);
            plot.drawNormalizedLine(xEnd-vecXsc, yEnd-vecYsc, xEnd, yEnd);
            double norm = 1./Math.sqrt(vecXsc*vecXsc + vecYsc*vecYsc);
            plot.drawNormalizedLine(xEnd, yEnd, xEnd-norm*(0.01*vecXsc+0.005*vecYsc), yEnd-norm*(0.01*vecYsc-0.005*vecXsc));
            plot.drawNormalizedLine(xEnd, yEnd, xEnd-norm*(0.01*vecXsc-0.005*vecYsc), yEnd-norm*(0.01*vecYsc+0.005*vecXsc));
            plot.setJustification(Plot.RIGHT);
            plot.setLimitsToFit(/*updateImg=*/false);
            plot.addLabel(0.99, 0.05, vectorText);
        }
        plot.setLineWidth(1);
        plot.setLegend(yLabels, Plot.AUTO_POSITION);
        return plot;
    }

    /** Returns an ImagePlus containing the plot, but does not show it */
    public ImagePlus getImagePlus() {
		if (plot == null) return null;
        ImagePlus imp = plot.getImagePlus();
        imp.setTitle(title);
        return imp;
    }

    /** Returns the Plot */
    public Plot getPlot() {
        return plot;
    }

    /** Returns an array of n different colors */
    public static Color[] getColors(int n) {
        int brightSteps = n/5+1;           //number of different brightnesses (main distinction is hue)
        double brightnessRange = brightSteps==2 ? 0.5 : 0.7;
        double minBrightness = 1. - 0.3*brightnessRange;
        Color[] colors = new Color[n];
        for (int i=0; i<n; i++) {
            double bright = brightSteps<2 ? 1. :
                    minBrightness + brightnessRange*(i%brightSteps)/(brightSteps-1);
            colors[i] = getColor(i/(double)n, bright);
        }
        return colors;
    }

    /** Returns a color with hue from 0 to 1 cycling 360 degrees; bright=1 is standard brightness */
    private static Color getColor (double hue, double bright) {
        final double rW=0.299, gW=0.587, bW=0.114;     //weight of the 3 colors R,G,B
        //final double sqrt3 = Math.sqrt(3);
        final double vtW = Math.sqrt(rW*rW-bW*bW);
        final double[] weights = new double[] {rW, gW, bW};
        final double[][] basicColors = new double[][]{  //have roughly the same contrast to white
            new double[]{.7,.6,0},      //orange
            new double[]{0,.6,0},       //green
            new double[]{0.05,0.1,1},   //blue
            new double[]{1,0,1},        //violet
            new double[]{1,0,0}         //red
        };
        hue = hue%1;
        int color1i = (int)(hue*basicColors.length);
        int color2i = (color1i+1)%basicColors.length;
        double remainder = hue*basicColors.length - color1i;
        //IJ.log("col1i="+color1i+" remain="+(float)remainder);
        double[] color = sumColor (basicColors[color1i], basicColors[color2i], 1.-remainder, remainder);
        double dist111 = 0;                 //distance of (111) plane in scaled color space * sqrt3
        for (int i=0; i<3; i++) {           //change brightness
            color[i] *= bright;
            dist111 += color[i]*weights[i];
        }
        // get max saturation, but make sure the color is 0 to 1 in unscaled space.
        // the brightness (weighted) remains fixed, i.e. move in a plane normal to (1,1,1) in scaled space
        double[] inPlane = new double[3];   //"hue" vector in scaled color space
        double bestSaturation = Double.MAX_VALUE; //factor for saturation
        for (int i=0; i<3; i++) {
            inPlane[i] = color[i]*weights[i] - dist111*weights[i];
            if (inPlane[i] < -1e-50) {           //avoid a color less than 0
                double maxSaturation = 1./(1. - color[i]/dist111);
                if (bestSaturation > maxSaturation) bestSaturation = maxSaturation;
            } /*else*/ if (inPlane[i] > 1e-50) { //avoid a color greater than 1 (weight in scaled space)
                double maxSaturation = ((1-dist111)*weights[i])/inPlane[i];
                if (bestSaturation > maxSaturation) bestSaturation = maxSaturation;
            }
        }
        for (int i=0; i<3; i++)   //white point + upscaled in-plane vector
            color[i] = (dist111*weights[i] + bestSaturation*inPlane[i])/weights[i];
        return new Color((int)(255*color[0]+.49), (int)(255*color[1]+.49), (int)(255*color[2]+.49));
    }

    /** weighted sum of two colors */
    private static double[] sumColor (double[] c1, double[] c2, double w1, double w2) {
        return new double[] {
            c1[0]*w1 + c2[0]*w2,
            c1[1]*w1 + c2[1]*w2,
            c1[2]*w1 + c2[2]*w2
        };
    }
}
