import ij.*;
import ij.IJ;
import ij.gui.*;
import ij.process.*;
import ij.plugin.PlugIn;
import java.util.Random;
import java.util.Arrays;
/**
 *  Finds an elliptical ROI to the edge of the LEED screen.
 *  Assumptions:
 *  - There is a "groove" of intensity around most of the LEED screen
 *    (after bandpass filtering in radial direction)
 *  - The edge is the steepest descent inside the groove.
 *  - The outline is an eellipse that is roughly circular and centered,
 *    and the outline is at a distance larger than 0.35*min(width, height)
 *    from the image center and mostly inside the image. For a centered
 *    circle, this means that the diameter must be more than 70% of the
 *    smallest image side length.
 * 
 *  The underlying math:
 *
 *  In polar coordinates, an off-axis circle of radius R with the center at (d,0)
 *  can be roughly fit by
 *    ln(r) = c + e*cos(phi)
 *  From the fit parameters c, e we get
 *    ln(R + d) = c + e,  ln(R - d) = c - e,   and, hence
 *    R = ½[exp(c + e) + exp(c - e)],  d = ½[exp(c + e) - exp(c - e)]
 *  We use this fit in the first step, to get an initial equation for the position
 *  of the center (but fitting cos, sin, because the center can be dsiplaced
 *  in x and y).
 *
 *  A centered ellipse is given by (x/a)^2 +(y/b)^2 = 1, hence in polar coordinates
 *  with x = r cos(phi), y = r sin(phi), and dividing both sides by r², we have
 *    1/r² = cos²(phi)/a² + sin²(phi)/b²
 *         = (1/a² + 1/b²)/2 + cos(2phi)*(1/a² - 1/b²)/2 = c + d*cos(2phi)
 *  where c, d are fit components of 1/r². Then the semi-axes a, b are
 *    a = 1/sqrt(c+d), b = 1/sqrt(c-d)
 *  here, the cosine fit component 'd' is negative. When taking the absolute value
 *  (as we do), we have to reverse its sign.
 *  This fit is used in the second step, where we use the center position from the first step.
 *
 *  In the final step, the steepest radial slope in a ring-like zone inside the fit ellipse
 *  determines the size of the ellipse is used (keeping the center, angular orientation and
 *  aspect ration constant).
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

public class LeedEllipseEdgeFinder implements PlugIn {
    final static int MIN_DATA_POINTS = 30;          //minimum number of data points for different angles for fitting (must be > 3*(N_FITPARAMS+1))
    final static double HIGHPASS_CUTOFF = 0.05;     //cutoff of radial high-pass filter before edge detection, as fraction of image size
    final static double LOWPASS_CUTOFF  = 0.005;    //cutoff of radial low-pass filter before finding the minimum, as fraction of image size
    final static double MAX_DEVIATION   = 0.01;     //maximum relative deviation of the radius in a fit; everything else is outliers
    final static int N_TRIES = 1000;                //attemps of random sample consensus (RANSAC), roughly 98% success for 60% outliers
    final static int N_FITPARAMS = 5;
    final static int CONST=0, COS=1, SIN=2, COS2=3, SIN2=4; //indices of fit params

    /** For testing only, this class can be run as ImageJ PlugIn */
    public void run(String arg) {
        ImagePlus imp = IJ.getImage();
        if (imp == null) return;
        Roi ellipseRoi = getEllipse(imp.getProcessor());
        if (ellipseRoi != null)
            imp.setRoi(ellipseRoi);
    }
    /** Returns an ij.gui.EllipseRoi at the edge of the ellipse with a strong negative radial gradient */
    static Roi getEllipse(ImageProcessor ip) {
        ip.setInterpolationMethod(ImageProcessor.BILINEAR);
        int width = ip.getWidth();
        int height = ip.getHeight();
        double xC = 0.5*(width-1);  //center coordinates
        double yC = 0.5*(height-1);
        int rMin = (int)Math.round(0.35*Math.min(width, height));
        int rMax = (int)Math.round(0.55*Math.min(width, height));
        // find the strongest descending step for all angles
        int nAngles = (int)Math.round(2*Math.PI*rMax);
        double[][] fitXValues = new double[nAngles][N_FITPARAMS];
        double[] grooveRadius = new double[nAngles];
        double[] grooveCurv = new double[nAngles];
        double[] radiusProfile = new double[rMax-rMin];
        double cutoffWavelengthH = HIGHPASS_CUTOFF*0.5*(width + height);
        int highpassKernelHalfwidth = LeedSmoother.bandwidthToM(1./cutoffWavelengthH);
        LeedSmoother smootherH = new LeedSmoother(highpassKernelHalfwidth);
        double cutoffWavelengthL = LOWPASS_CUTOFF*0.5*(width + height);
        int lowpassKernelHalfwidth = LeedSmoother.bandwidthToM(1./cutoffWavelengthL);
        LeedSmoother smootherL = new LeedSmoother(lowpassKernelHalfwidth);
        LeedMultiVariateFitter fitter1 = new LeedMultiVariateFitter(N_FITPARAMS);
        LeedMultiVariateFitter fitter2 = new LeedMultiVariateFitter(N_FITPARAMS);
        double[] bestFitParameters = null;
        // two or three runs; the first for finding the center; the final one for accurate ellipse fitting
        int lastTask = 1;
        for (int iTask = 0; iTask <= lastTask; iTask++) {
            //DEBUG Plot plot = new Plot("radial profile filtd", "r-"+rMin, "I");
            for (int iA=0; iA<nAngles; iA++) {
                double cos = Math.cos(iA*(2*Math.PI/nAngles));
                double sin = Math.sin(iA*(2*Math.PI/nAngles));
                double maxCurv = 0;
                double grooveR = -1;
                for (int r=rMin; r<rMax; r++) {
                    double x = xC + r*cos;
                    double y = yC + r*sin;
                    radiusProfile[r-rMin] = (x>= 0 && y>=0 && x<=width-1 && y<=height-1) ?
                            ip.getInterpolatedValue(x, y) : Double.NaN;
                }
                double[] radiusProfileSmooth = smootherH.smooth(radiusProfile);
                for (int r=rMin; r<rMax; r++)
                    radiusProfile[r-rMin] -= radiusProfileSmooth[r-rMin];       //subtract smoothed, to get a high-pass filter
                double[] radiusProfileFiltd = smootherL.smooth(radiusProfile);  //low-pass filter to remove noise
                //find the most pronounced minimum; typically, this is the "groove" around the screen
                for (int r=rMin+1; r<rMax-1; r++) {
                    double left = radiusProfileFiltd[r-rMin-1], mid = radiusProfileFiltd[r-rMin], right=radiusProfileFiltd[r-rMin+1];
                    if (mid < left && mid < right) {
                        double curv = 0.5*(left + right) - mid;     //curvature
                        if (curv > maxCurv) {
                            maxCurv = curv;
                            double b = 0.25*(right-left);
                            grooveR = r - b/curv;      //minimum position with subpixel precision, from a parabola through 3 points
                        }
                    }
                }
                //DEBUG if (iTask==0) plot.setColor(java.awt.Color.BLACK);
                //DEBUG if (iTask==0) plot.addPoints(null, radiusProfileFiltd, Plot.LINE);
                //DEBUG if (iTask==0) plot.setColor(java.awt.Color.RED);
                //DEBUG if (iTask==0) plot.addPoints(new double[]{grooveR-rMin,grooveR-rMin}, new double[]{-20,20}, Plot.LINE);
                //DEBUG if (iTask==0) plot.addToStack();
                fitXValues[iA][CONST] = 1.;
                fitXValues[iA][COS] = cos;
                fitXValues[iA][SIN] = sin;
                fitXValues[iA][COS2] = 2*sqr(cos) - 1;
                fitXValues[iA][SIN2] = 2*cos*sin;
                grooveRadius[iA] = grooveR;
                grooveCurv[iA] = maxCurv;
            }
            //DEBUG if (iTask==0) plot.show();
            // analyze curvatures at the groove (the minimum) and find average
            int grooveCount = 0;
            double sumCurv = 0;
            for (int iA=0; iA<nAngles; iA++) {
                if (grooveCurv[iA] > 0) {
                    sumCurv += grooveCurv[iA];
                    grooveCount++;
                }
            }
            if (grooveCount < 2*MIN_DATA_POINTS)
                return null;                                //we should definitely have more than 40 inital points
            double avgCurvature = sumCurv/grooveCount;
            double curvThreshold = 0.02*avgCurvature;
            int nAboveThreshold = 0;
            for (int iTimeout=0; iTimeout<31; iTimeout++) { //timeout: 0.8^31 ~ 1e-3; we don't expect a threshold that low
                nAboveThreshold = 0;
                for (int iA=0; iA<nAngles; iA++)
                    if (grooveCurv[iA] >= curvThreshold)
                        nAboveThreshold++;
                if (nAboveThreshold > 0.6*grooveCount) break; //we want more than 60% of all points above the threshold
                curvThreshold *= 0.8;                       //otherwise decrease the threshold
            }
            if (nAboveThreshold < MIN_DATA_POINTS)          //not enough points for fitting
                return null;
            // Fit a + b*scos + c*sin into subsets of the data, in the spirit of
            // random sample consensus (RANSAC)
            int nEqualFits = 0;
            double bestFigureOfMerit = 0;
            Random random = new Random(0);

            for (int iTry=0; iTry<N_TRIES; iTry++) {
                fitter1.clear();
                // select N_FITPARAMS+1 random points for the initial fit (RANSAC)
                int count = 0;
                int[] selectedPoints = new int[N_FITPARAMS];
                while (count < N_FITPARAMS+1) {
                    int iA = random.nextInt(grooveCurv.length);
                    if (grooveCurv[iA] >= curvThreshold) {
                        for (int i=0; i<count; i++)
                            if (Math.abs(selectedPoints[i] - iA) < 2)
                                continue;                   //don't select the same point or neighboring points
                        double yValue = iTask != lastTask ? Math.log(grooveRadius[iA]) : 1./sqr(grooveRadius[iA]);
                        fitter1.addPoint(yValue, fitXValues[iA]);
                        if (count < selectedPoints.length-1)
                            selectedPoints[count] = iA;
                        count++;
                    }
                }
                // refine: apply the fit to all inliers, and redo the fit for these (in two cycles)
                double maxDeviation = iTask != lastTask ? MAX_DEVIATION : MAX_DEVIATION*2*fitter1.getFitParameters()[CONST];
                for (int iRefine=0; iRefine <2; iRefine++) {
                    LeedMultiVariateFitter previousFitter = iRefine==0 ? fitter1 : fitter2;
                    LeedMultiVariateFitter currentFitter  = iRefine==0 ? fitter2 : fitter1;
                    currentFitter.clear();
                    for (int iA=0; iA<nAngles; iA++) {
                        double yValue = iTask != lastTask ? Math.log(grooveRadius[iA]) : 1./sqr(grooveRadius[iA]);
                        double fitValue = previousFitter.getFitValue(fitXValues[iA]);
                        double deviationSqr = sqr(fitValue - yValue);
                        boolean isInlier = deviationSqr < sqr(maxDeviation);
                        if (isInlier)
                            currentFitter.addPoint(yValue, fitXValues[iA]);
                    }
                    if (currentFitter.getNPoints() < MIN_DATA_POINTS/3) break;
                }
                // evaluate the fit: how many inliers and how good are these
                LeedMultiVariateFitter previousFitter = fitter1;
                if (previousFitter.getNPoints() < MIN_DATA_POINTS/3) continue;
                int nInliers = 0;
                double sumDeviationsSqr = 0;

                for (int iA=0; iA<nAngles; iA++) {
                    double yValue = iTask != lastTask ? Math.log(grooveRadius[iA]) : 1./sqr(grooveRadius[iA]);
                    double fitValue = previousFitter.getFitValue(fitXValues[iA]);
                    double deviationSqr = sqr(fitValue - yValue);
                    boolean isInlier = deviationSqr < sqr(maxDeviation);
                    if (isInlier) {
                        sumDeviationsSqr += deviationSqr;
                        nInliers++;
                    }
                }
                double meanDeviationSqr = sumDeviationsSqr/nInliers;
                double figureOfMerit = nInliers*nInliers/(meanDeviationSqr + 0.1*sqr(maxDeviation));
                if (figureOfMerit > bestFigureOfMerit) {
                    bestFigureOfMerit = figureOfMerit;
                    bestFitParameters = previousFitter.getFitParameters();
                    nEqualFits = 1;
                } else if (figureOfMerit == bestFigureOfMerit && Arrays.equals(bestFitParameters, previousFitter.getFitParameters()))
                    nEqualFits++;
            } //for iTry
                if (!(bestFigureOfMerit > 0)) return null;

            //DEBUG String s="task="+iTask+" nEqualFits="+nEqualFits; for (int i=0; i<N_FITPARAMS; i++) s+=" "+(char)('a'+i)+"="+(float)bestFitParameters[i]; IJ.log(s);
            double c = bestFitParameters[CONST];
            double e = Math.sqrt(sqr(bestFitParameters[COS]) + sqr(bestFitParameters[SIN]));
            if (e == 0) e = 1;                                  //avoid division by zero if both sin and cos components are zero

            double offset = iTask != lastTask ?                 //calculate distance from the center
                    0.5*(Math.exp(c + e) - Math.exp(c - e)) :   //when fitting ln(r)
                    -0.5*e/(c*Math.sqrt(c));                    //when fitting 1/r^2

            xC += bestFitParameters[COS]*offset/e;
            yC += bestFitParameters[SIN]*offset/e;
            if (offset > 0.03*rMax)
                lastTask = 2;                                   //when the center is rather far off, one more round for refining

        //DEBUG double[]v = new double[nAngles];
        //DEBUG for (int iA=0; iA<nAngles; iA++) {v[iA]=0; for (int j=0; j<N_FITPARAMS; j++) v[iA]+=bestFitParameters[j]*fitXValues[iA][j];}
        //DEBUG double[]y = new double[nAngles];
        //DEBUG for (int iA=0; iA<nAngles; iA++) y[iA] = iTask != lastTask ? Math.log(grooveRadius[iA]) : 1./sqr(grooveRadius[iA]);
        //DEBUG Plot plot1 = new Plot("EllFit", "angle index", "r-"+rMin);
        //DEBUG plot1.addPoints(null, y, Plot.CIRCLE);
        //DEBUG plot1.setColor(java.awt.Color.RED);
        //DEBUG plot1.addPoints(null, v, Plot.LINE);
        //DEBUG plot1.show();

        } //for iTask

        double ellipt = Math.sqrt(sqr(bestFitParameters[COS2]) + sqr(bestFitParameters[SIN2]));
        double semiMajor = 1./Math.sqrt(bestFitParameters[CONST] - ellipt);
        double semiMinor = 1./Math.sqrt(bestFitParameters[CONST] + ellipt);
        double angle = 0.5*Math.atan2(-bestFitParameters[SIN2], -bestFitParameters[COS2]);
        //DEBUG IJ.log("semiMajor="+IJ.d2s(semiMajor)+" semiMinor="+IJ.d2s(semiMinor)+" angle="+IJ.d2s(angle));
        if (Double.isNaN(angle)) angle = 0;
        double cosMajor = Math.cos(angle), sinMajor = Math.sin(angle);
        // From the border inwards, find the steepest radial descent inside the groove (allow up to 5% smaller)
        int nRadii = (int)Math.round(0.05*semiMajor+1);
        double bestSizeFactor = Double.NaN;
        double lastSizeFactor = 1.0;
        double lastIntensity = Double.NaN;
        double lastSlope = -Double.MAX_VALUE;
        double largestSlope = 0;
        for (int iR=0; iR<nRadii; iR++) {
            double sizeFactor = (semiMajor + 1 - iR)/semiMajor;
            Roi roi = new EllipseRoi(
                    xC + semiMajor*sizeFactor*cosMajor, yC + semiMajor*sizeFactor*sinMajor, //one vertex
                    xC - semiMajor*sizeFactor*cosMajor, yC - semiMajor*sizeFactor*sinMajor, //other vertex
                    semiMinor/semiMajor);                                                   //aspect ratio
            double intensity = getIntensityOnRoiLine(ip, roi);
            double slope = intensity - lastIntensity;       //typically positive; we have an increase towards the center
            if (slope > largestSlope) {                     //higher slope than ever before? reset
                largestSlope = slope;
                bestSizeFactor = Double.NaN;
            }
        //DEBUG IJ.log((float)sizeFactor+": I="+(float)intensity+" dI="+(float)slope);
            lastIntensity = intensity;
            if (slope <= lastSlope && Double.isNaN(bestSizeFactor)) //the first decrease of the slope
                bestSizeFactor = lastSizeFactor;            //we are beyond the position of steepest descent, take previous
            lastSlope = slope;
            lastSizeFactor = sizeFactor;
        }
        xC += 0.5;      //0.5 pixel shift: in ImageJ, boundary coordinates are to the top left of pixels
        yC += 0.5;
        Roi roi = new EllipseRoi(
                xC + semiMajor*bestSizeFactor*cosMajor, yC + semiMajor*bestSizeFactor*sinMajor,
                xC - semiMajor*bestSizeFactor*cosMajor, yC - semiMajor*bestSizeFactor*sinMajor,
                semiMinor/semiMajor);
        return roi;
    }

    /** Measures the median pixel value along the roi border line (for area rois only) */
    public static double getIntensityOnRoiLine(ImageProcessor ip, Roi roi) {
        int width = ip.getWidth();
        int height = ip.getHeight();
        FloatPolygon fp = roi.getFloatPolygon();    //contains vertex coordinates of the approximation polygon
        double length = roi.getLength();
        int nPoints = (int)Math.ceil(length);       //number of points that we have to measure
        double inc = length/nPoints;
        int iSegmentEnd = -1;
        double xVertexStart = 0;
        double yVertexStart = 0;
        double xVertexEnd = fp.xpoints[fp.npoints-1];
        double yVertexEnd = fp.ypoints[fp.npoints-1];
        double segmentLength = 0;
        double positionOnSement = 0;
        double[] intensities = new double[nPoints];
        int count = 0;
        for (int i=0; i<nPoints; i++) {
            if (iSegmentEnd < 0 || (positionOnSement > segmentLength && iSegmentEnd < fp.npoints-1)) {
                iSegmentEnd++;                      //we enter a new segment on the polygon, with this index of the end vertex
                positionOnSement -= segmentLength;
                xVertexStart = xVertexEnd;
                yVertexStart = yVertexEnd;
                xVertexEnd = fp.xpoints[iSegmentEnd];
                yVertexEnd = fp.ypoints[iSegmentEnd];
                segmentLength = Math.sqrt(sqr(xVertexEnd - xVertexStart) + sqr(yVertexEnd - yVertexStart));
            }
            double x = xVertexStart + (xVertexEnd - xVertexStart)*positionOnSement/segmentLength;
            double y = yVertexStart + (yVertexEnd - yVertexStart)*positionOnSement/segmentLength;
            if (x>= 0 && y>=0 && x<=width-1 && y<=height-1) {
                double intensity = ip.getInterpolatedValue(x, y);
                if (!Double.isNaN(intensity))
                    intensities[count++] = intensity;
            }
            positionOnSement += inc;
        }
        if (count < 2) return Double.NaN;
        Arrays.sort(intensities, 0, count);
        double median = (count&1) == 0 ?            //for even array length, use average of center two
                0.5*(intensities[count/2-1] + intensities[count/2]) : intensities[count/2];
        return median;
    }


    static double sqr(double x) {
        return x*x;
    }
}
