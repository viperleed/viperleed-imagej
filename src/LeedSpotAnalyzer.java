import ij.IJ;
import ij.util.Tools;
import ij.ImagePlus;
import ij.process.*;
import ij.gui.*;
import java.util.Arrays;

/**
 *  This class analyzes a single spot with subpixel precision (aperture photometry).
 *  For background shapes CIRCLE and OVAL, the spot intensity is integrated over
 *  a circle and the background is a linear fit over an area around that circle.
 *  The outer border of the background area can be either a circle (with radius*sqrt2)
 *  or,  for OVAL, an ellipse with the same semiminor axis as the integration
 *  circle radius and the semimajor axis twice as large. The ellipse is oriented
 *  with the major axis perpendicular to the line between the spot and the screen
 *  center. For spots tangentially blurred due to poor azimuthal alignment (e.g.,
 *  thin films with imprefect azimuthal orientation), use shape=AZIMUTH_BLUR.
 *  Then the spot integration area is an ellipse with 2:1 ratio between major
 *  and minor axis, and the outer border of the background is the smallest circle
 *  encompassing the ellipse. In this cse, 'radius' is the semiminor axis of the
 *  integration ellipse.
 *  In all three cases, the area of the integration circle and the background area
 *  are equal. Since the noise is usually dominated by shot noise and the
 *  inhomogeneity of the grid, which both have a higher impact at higher intensities,
 *  the noise of the output (after background subtraction) will be dominated by
 *  the integration over the spot, not by the noise of the background.
 *  In the worst case, at very low spot intensities, the noise in the integration
 *  and background areas will be equal, thus the background noise increases the
 *  noise of the output by a factor of sqrt(2).
 *
 *  Output items:
 *  X - x coordinate of the center of the integration circle with subpixel resolution
 *  Y - y coordinate of the center of the integration circle with subpixel resolution
 *  INTEGRAL - integrated intensity of the spot, after background subtraction
 *  R_SIZE -  Assuming a 2D Gaussian, the sigma in radial direction with respect to
 *            the center given by x0, y0. This is based on the normalized central
 *            2nd moments.  For spots close to the center (within 'radius'), this rather gives
 *            the semiminor axis of the ellipse given by the normalized central 2nd moments.
 *            A correction is applied to account for small integration radii where sigma
 *            from 2nd moments is overestimated. The corrections assumes a circular Gaussian
 *            spot. (For AZIMUTH_BLUR, the correction is not very accurate and assumes
 *            an elliptical spot in tangetial direction, with an aspect ratio according
 *            to the integration ellipse).
 *  T_SIZE  - Assuming a 2D Gaussian, the sigma in tangential direction with respect to
 *            the center given by x0, y0. This is based on the normalized central
 *            2nd moments.  For spots close to the center (within 'radius'), , this rather gives
 *            the semimajor axis of the ellipse given by the normalized central 2nd moments.
 *            For AZIMUTH_BLUR, the value is scaled by the aspect ratio of the integration ellipse,
 *            so the ratio between the T_SIZE item and the spot radius can be interpreted
 *            as for circular spots.
 *            A correction is applied to account for small integration radii where sigma
 *            from 2nd moments is overestimated. The corrections assumes a circular Gaussian
 *            spot. (For AZIMUTH_BLUR, the correction is not very accurate and assumes
 *            an elliptical spot in tangetial direction, with an aspect ratio according
 *            to the integration ellipse).
 *  PSIGNIFICANCE - by which factor the spot is higher than the standard deviation of the background.
 *            The spot height is taken as the geometric mean of the peak height
 *            (assuming a 2D Gaussian as described above) and the average intensity
 *            in the integration circle.
 *            Values below NO_SIGNIFICANCE = 0.5 are set to 0 (these values are reserved for
 *            significance of positions inferred from neighboring spots).
 *            If the spot is not inside the mask, the value is NaN.
 *  BACKGR  - average background intensity
 *  BSIGMA  - standard deviation of background intensity
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


public class LeedSpotAnalyzer {
    /** Indices of the items in the output array. Note that the coordinates x, y must be first */
    public static final int X=0, Y=1, INTEGRAL=2, R_SIZE=3, T_SIZE=4, PSIGNIFICANCE=5, BACKGR=6, BSIGMA=7;
    /** Names of these data columns */
    public static final String[] DATA_NAMES = new String[] {"x", "y", "rawInt","rSize", "tSize", "signif", "backgr", "bsigma"};
    /** shapes of background area */
    public static final int CIRCLE=0, OVAL=1, AZIMUTH_BLUR=2;
    /** measured significance values below this are set to 0. */
    public static final double NO_SIGNIFICANCE = 0.5;
    static final int OUTPUT_SIZE=8;     //size of output array
    static final int MAX_ITERATIONS = 5;    //maximum number of iterations for iterative center determination
    static final double MIN_XYSQR_CONVERGENCE = 0.3; //INTEGRAL position must be converged better than this (squared deviation from last step)
    static final double MIN_SIGMAS_FOR_XY = 1.5;     //INTEGRAL must be about this many sigmas above background variations to evaluate x,y

    /** Determines the exact position and analyzes a spot with a given starting position for the search.
     *  Before doing the measurement of the spot properties, tries to position the spot accurately.
     *  @param xs, ys The starting position for searching, in pixel coordinates. The spot must be nearby
     *    (within 'radius').
     *  @param x0, y0  The screen center, or for shape=AZIMUTH_BLUR, the position of the (0,0) spot.
     *  @param shape This determines the geometry of the integration and background area:
     *    The integration area is circle with the given radius, except for AZIMUTH_BLUR, where
     *    an ellipse with the major axis in tangetial direction [with respect to the (0,0) spot
     *    taken at the center] is used, The semiminor axis is 'radius' is used, and the length of
     *    the major axis is given by radius and azBlurRadians.
     *    For shapes CIRCLE and OVAL, the background region has twice the area of the integration region;
     *    the outer border of the background is circular for CIRCLE and elliptical for OVAL,
     *    with the long axis in tangential direction. With shape=AZIMUTH_BLUR, as long as the
     *    integration ellipse has an aspect ratio less than sqrt(2), the outer background border
     *    is a circle with r=sqrt(2)*radius (where radius is the semiminor axis of the integration ellipse).
     *    For more elongated integration ellipses, the background outline is an ellipse
     *    obtained by stretching the integration ellipse in radial direction [towards and away
     *    from the (0,0) spot] by a factor of sqrt 2.
     *  @param radius Radius of the integration area (semiminor axis for shape=AZIMUTH_BLUR).
     *    The radius also defines how far the refined position may deviate from the initial position.
     *  @param azBlurRadians For shape=AZIMUTH_BLUR, in the limit of large distances from the (0,0) spot,
     *    the semimajor axis of the integration ellipse over the distance from the (0,0) spot.
     *    this corresponds to the half-angle in radians under which the integration ellipse would be seen
     *    from the position of the (0,0) spots. For finite distances from the (0,0) spot, the
     *    semimajor axis smoothly changes from 'radius' when close to the (0,0) spot (i.e., we have
     *    an integration circle) to the distance from the (0,0) spot times azBlurRadians at very high
     *    distances.
     *  @param minSignificance Threshold for accepting a spot. The value is typically about 2;
     *    larger numbers lead to stronger background suppression
     *  @param image The input image.
     *  @param mask  Mask with values different from 0 defining the allowed area for measurements.
     *    Measurements are taken only if all pixels of the integration and background regions
     *    have mask pixel values > 0. Image and mask must have the same size.
     *  @param output An output array may be supplied to reduce unnecessary creation of objects;
     *    if it is null, a new array is created for the output.
     *  @return  The output array contains the values corresponding go X, Y, ... BSIGMA as explained
     *    in the description of the class (see above).
     *    If the position cannot be determined with sufficient accuracy, the X & Y output values
     *    are set to NaN and PSIGNIFICANCE in the output array is less than minSignificance.
     *    If the position is partly outside the valid area defined by the mask, all return values
     *    including the significance are NaN.
     *   */
    public static double[] centerAndAnalyzeSpot(double xs, double ys, double x0, double y0,
            int shape, double radius, double azBlurRadians, double minSignificance,
            FloatProcessor image, ByteProcessor mask, double[] output) {
        double lastX = xs, lastY = ys;
        double[] firstResult = null;
        for (int iter=0; iter<MAX_ITERATIONS; iter++) {
            if (output == null || output.length < OUTPUT_SIZE)
                output = new double[OUTPUT_SIZE];
            output = analyzeSpot(lastX, lastY, x0, y0, shape, radius, azBlurRadians, minSignificance, image, mask, output);
            if (iter == 0) {
                if (Double.isNaN(output[X]))        //can't refine position, significance too low
                    return output;
                else
                    firstResult = (double[])output.clone();
            } else {                                //check how far the spot has moved from the initial guess
                if (shape == AZIMUTH_BLUR) {
                    if (Double.isNaN(output[X]))    //significance too low, abort position refinement
                        break;
                    double rhoX = xs - x0;          //vector from center to spot, x component
                    double rhoY = ys - y0;          //vector from center to spot, y component
                    double rho = Math.sqrt(rhoX*rhoX + rhoY*rhoY);
                    if (rho >= radius) {            //far enough from the center for really using AZIMUTH_BLUR
                        double cosR = rhoX/rho;
                        double sinR = rhoY/rho;
                        double deltaR = (output[X] - xs)*cosR + (output[Y] - ys)*sinR;  //radial displacement from xs, ys
                        double deltaT = (output[Y] - ys)*cosR - (output[X] - xs)*sinR;  //tangential displacement from xs, ys
                        double radiusT = Math.sqrt(sqr(radius) + sqr(azBlurRadians)*(rhoX*rhoX + rhoY*rhoY));
                        if (!(sqr(deltaR/radius) + sqr(deltaT/radiusT) <= 1))
                            break;                   //too far from xs, ys; abort position refinement
                    } else if (!(sqr(output[X] - xs) + sqr(output[Y] - ys) <= sqr(radius)))
                        break;                       //too far from xs, ys for CIRCLE (AZ_BLUR at the center); abort
                } else if (!(sqr(output[X] - xs) + sqr(output[Y] - ys) <= sqr(radius))) //CIRC or OVAL
                    break;   //too far from original position or NaN (=signif. too low), abort position refinement

                double deltaSqr = sqr(output[X] - lastX) + sqr(output[Y] - lastY);
                if (deltaSqr < MIN_XYSQR_CONVERGENCE)
                    return output;                  //position converged, done
                lastX = output[X];
                lastY = output[Y];
            }
        }
        // not converged or running into the border
        if (firstResult[PSIGNIFICANCE] > minSignificance && Double.isNaN(output[PSIGNIFICANCE]))
            return output;                  //we have had a peak, but ran into the border, return all NaNs
        else {
            firstResult[X] = Double.NaN;
            firstResult[Y] = Double.NaN;
            firstResult[PSIGNIFICANCE] = 0; //probably peak too low
            return firstResult;
        }
    }

    /** Analyzes a spot at a given position xs, ys, with screen center at x0, y0,
     *  using integration circle radius r and background shape CIRCLE or OVAL.
     *  The output array contains the values corresponding go X, Y, ... BSIGMA as explained above.
     *  Image and mask must have the same size; measurements are taken only at positions where the
     *  integration area and background are have mask pixels >0.
     *  The INTEGRAL, PRADIUS ... BACKGR, BSIGMA values are NaN if the integration or background area
     *  is partly or fully outside the image bounds or mask foreground.
     *  If the position cannot be determined with sufficient accuracy, the X and Y values
     *  in the output array are set to NaN, and significance below NO_SIGNIFICANCE.
     *  If the position is partly outside the mask, all return values including the significance are NaN.
     *  An output array may be supplied to reduce unnecessary creation of objects; if it is null, a new
     *  array is created for the output. */
    public static double[] analyzeSpot(double xs, double ys, double x0, double y0,
            int shape, double radius, double azBlurRadians, double minSignificance,
            FloatProcessor image, ByteProcessor mask, double[] output) {
        int width = image.getWidth();
        int height = image.getHeight();
        if (output == null || output.length < OUTPUT_SIZE)
            output = new double[OUTPUT_SIZE];
        Arrays.fill(output, Double.NaN);
        if (xs<radius || xs>width-1-radius || ys<radius || ys>height-1-radius)
            return output;              //clearly outside? all NaN
        double rhoX = xs - x0;          //vector from center to spot, x component
        double rhoY = ys - y0;          //vector from center to spot, y component
        double rho = Math.sqrt(rhoX*rhoX + rhoY*rhoY);
        boolean atCenter = rho < radius;  //at the center, radial and tangential directions are undefined
        double cosR = rhoX/rho;
        double sinR = rhoY/rho;
        if (atCenter) shape = CIRCLE;
        int bgShape = shape;
        double bgRadius = radius*Math.sqrt(2); // radius of background circle (if we want one)

        double xBounding = 0, yBounding = 0;
        double invRadiusT = 0;          //AZIMUTH_BLUR: 1/radius in tangential direction
        if (shape == AZIMUTH_BLUR) {
            double radiusT = Math.sqrt(sqr(radius) + sqr(azBlurRadians)*(rhoX*rhoX + rhoY*rhoY));
            invRadiusT = 1./radiusT;
            if (radiusT < bgRadius) {   //integration ellipse fits inside circular background outline
                bgShape = CIRCLE;
            } else {                    //background area with s*radius in radial direction
                xBounding = Math.sqrt(2*sqr(radius*cosR)+sqr(radiusT*sinR))+0.5;
                yBounding = Math.sqrt(2*sqr(radius*sinR)+sqr(radiusT*cosR))+0.5;
            }
        }
        if (bgShape == CIRCLE) {
            xBounding = bgRadius+0.5;
            yBounding = bgRadius+0.5;
        } else if (bgShape == OVAL) {
            xBounding = radius*Math.sqrt(sqr(cosR)+4*sqr(sinR))+0.5;
            yBounding = radius*Math.sqrt(sqr(sinR)+4*sqr(cosR))+0.5;
        }
        int xmin = Math.max((int)Math.floor(xs-xBounding), 0);  //Rectangle enclosing the pixels that we have to measure
        int ymin = Math.max((int)Math.floor(ys-yBounding), 0);  //(inclusive)
        int xmax = Math.min((int)Math.ceil(xs+xBounding), width-1);
        int ymax = Math.min((int)Math.ceil(ys+yBounding), height-1);

        // find background linear fit
        double sumW=0, sumX=0, sumY=0, sumXX=0, sumXY=0, sumYY=0, sumV=0, sumXV=0, sumYV=0;
        for (int y=ymin; y<=ymax; y++) {
            for (int x=xmin; x<=xmax; x++) {
                double dx = x - xs, dy = y - ys;
                double totalWeight = 0;
                if (bgShape == CIRCLE)
                    totalWeight = getCircleWeight(dx, dy, bgRadius);
                else if (shape==OVAL)
                    totalWeight = getOvalWeight(dx, dy, cosR, sinR, radius);
                else                    //if (bgShape==AZIMUTH_BLUR)
                    totalWeight = getOvalWeight(dx, dy, cosR, sinR, 1./bgRadius, invRadiusT);
                //IJ.log(x+","+y+","+totalWeight);

                if (totalWeight == 0) continue;
                if (mask.get(x,y) == 0) {       //pixel not in mask, but should be measured?
                    if (totalWeight < 0.5) continue;  //ignore out-of-mask if low weight (i.e., mask touches the border)
                    else return output;         //return all NaN because we can't measure
                }
                double peakWeight = shape == AZIMUTH_BLUR?
                        getOvalWeight(dx, dy, cosR, sinR, 1./radius, invRadiusT) :
                        getCircleWeight(dx, dy, radius);    //INTEGRAL weight
                //IJ.log(x+","+y+","+(totalWeight-peakWeight));
                if (peakWeight == 1) continue;  //don't care about INTEGRAL area pixels
                double w = totalWeight - peakWeight; //actual background weight, excluding INTEGRAL area
                double v = image.getf(x, y);    //pixel value;
                sumW += w; sumX += w*dx; sumY += w*dy;
                sumXX += w*dx*dx; sumXY += w*dx*dy; sumYY += w*dy*dy;
                sumV += w*v; sumXV += w*dx*v; sumYV += w*dy*v;
            }
        }

        double oneOverDenominator = 1.0 / (sumW*(sqr(sumXY)- sumXX*sumYY)-2*sumX*sumXY*sumY + sumXX*sqr(sumY)+ sumYY*sqr(sumX));
        double xSlope = (sumW*(sumYV*sumXY - sumXV*sumYY) + sumV*(sumX*sumYY - sumXY*sumY) + sumY*(sumXV*sumY - sumYV*sumX)) * oneOverDenominator;
        double ySlope = (sumW*(sumXV*sumXY - sumYV*sumXX) + sumV*(sumY*sumXX - sumXY*sumX) + sumX*(sumYV*sumX - sumXV*sumY)) * oneOverDenominator;
        if (Double.isNaN(xSlope)) xSlope=0; //should never happen
        if (Double.isNaN(ySlope)) ySlope=0;
        double offset = (sumV*(sumXY*sumXY - sumXX*sumYY) + sumXV*(sumX*sumYY - sumXY*sumY) + sumYV*(sumXX*sumY - sumX*sumXY)) * oneOverDenominator;

        // For the data with the fit subtracted, find the stddev of the backgound area.
        // For the INTEGRAL area, find the 0th, 1st, and 2nd moments (0th moment = integral, 1st moment for position)
        double sumBW=0, sumB=0, sumB2=0;
        double sumPW=0, sumP=0, sumXP=0, sumYP=0, sumXXP=0, sumXYP=0, sumYYP=0;
        for (int y=ymin; y<=ymax; y++) {
            for (int x=xmin; x<=xmax; x++) {
                double dx = x - xs, dy = y - ys;
                double totalWeight = 0;
                if (bgShape == CIRCLE)
                    totalWeight = getCircleWeight(dx, dy, bgRadius);
                else if (shape==OVAL)
                    totalWeight = getOvalWeight(dx, dy, cosR, sinR, radius);
                else                    //if (bgShape==AZIMUTH_BLUR)
                    totalWeight = getOvalWeight(dx, dy, cosR, sinR, 1./bgRadius, invRadiusT);
                if (totalWeight == 0) continue;
                double pw = shape == AZIMUTH_BLUR ?
                        getOvalWeight(dx, dy, cosR, sinR, 1./radius, invRadiusT) :
                        getCircleWeight(dx, dy, radius);    //INTEGRAL weight
                double bw = totalWeight - pw;                    //backgound weight (0 in INTEGRAL area)
                double pixelValue = image.getf(x, y);
                double v = pixelValue - (offset + dx*xSlope + dy*ySlope);
                sumBW += bw; sumB += bw*v; sumB2 += bw*v*v;
                if (pw > 0) {
                    sumPW += pw;
                    sumP += pw*v; sumXP += pw*v*dx; sumYP += pw*v*dy;
                    sumXXP += pw*v*dx*dx; sumXYP += pw*v*dx*dy; sumYYP += pw*v*dy*dy;
                }
            //if(Math.abs(xs-254)+Math.abs(ys-315)<4 && bw>0.5)IJ.log(dx+","+dy); //list relative coordinates of background
            }
        }
        double peakIntegral = sumP;
        output[INTEGRAL] = peakIntegral;                        // the integral is also needed for sub-threshold spots (unless out-of-mask)
        double peakX = sumXP*(1/sumP);                          // coordinates relative to xs, ys: from 1st moment
        double peakY = sumYP*(1/sumP);
        // calculate covariance matrix
        double peakSigmaXsqr = sumXXP*(1/sumP) - peakX*peakX;   // 2nd central moment   sigma_xx = mu_20'
        double peakSigmaYsqr = sumYYP*(1/sumP) - peakY*peakY;   // 2nd central moment   sigma_yy = mu_02'
        double peakXY = sumXYP*(1/sumP) - peakX*peakY;          // mixed central moment sigma_xy = mu_11'
        double peakMinorSqr = 0.5*(peakSigmaXsqr + peakSigmaYsqr) - Math.sqrt
                (sqr(0.5*(peakSigmaXsqr - peakSigmaYsqr)) + sqr(peakXY));   //smaller eigenvalue of moment matrix: shorter extension

        if (peakMinorSqr <= 0) {                                // not a maximum: at least one negative 2nd moment
            output[PSIGNIFICANCE] = 0.0;                        // everything else (except integral) in the output remains NaN
            return output;
        }
        double peakMajorSqr = 0.5*(peakSigmaXsqr + peakSigmaYsqr) + Math.sqrt
                (sqr(0.5*(peakSigmaXsqr - peakSigmaYsqr)) + sqr(peakXY));   //larger eigenvalue of moment matrix: longer extension
        // Rotate coordinates to radial and tangential direction (unless we are at the center)
        // ur = x*cosR + y*sinR, ut = y*cosR - x*sinR; thus:
        // sum(ur^2) = sum(x^2)*cosR^2 + sum(y^2)*sinR^2 + 2*sum(x*y)*cosR*sinR
        // sum(ut^2) = sum(y^2)*cosR^2 + sum(x^2)*sinR^2 - 2*sum(x*y)*cosR*sinR
        // sum(ur*ut) = (sum(y^2) + sum(x^2))*cosR*sinR + sum(x*y)*(cosR^2 - sinR^2)
        double peakSigmaRsqr = atCenter? peakMinorSqr :
                peakSigmaXsqr*sqr(cosR) + peakSigmaYsqr*sqr(sinR) + 2*peakXY*(cosR*sinR);
        double peakSigmaTsqr = atCenter? peakMajorSqr :
                peakSigmaYsqr*sqr(cosR) + peakSigmaXsqr*sqr(sinR) - 2*peakXY*(cosR*sinR);
        if (shape == AZIMUTH_BLUR && !atCenter) {
            double peakMixedRT = (peakSigmaYsqr - peakSigmaXsqr)*(cosR*sinR) + peakXY*(sqr(cosR)-sqr(sinR));
            peakMixedRT *= radius*invRadiusT;
            peakSigmaTsqr *= sqr(radius*invRadiusT);            // scale size down in tangential (blurred) direction
            peakMajorSqr = 0.5*(peakSigmaRsqr + peakSigmaTsqr) + Math.sqrt
                (sqr(0.5*(peakSigmaRsqr - peakSigmaTsqr)) + sqr(peakMixedRT));   //larger eigenvalue after downscaling
        }
        double peakSigmaT = Math.sqrt(peakSigmaTsqr);
        double peakSigmaR = Math.sqrt(peakSigmaRsqr);

        double backIntensity = sumB*(1/sumBW);                  // should be close to zero, but not exact (pixel quantization, mask edge)
        double backIntensitySigma = Math.sqrt(sumB2*(1/sumBW) - sqr(backIntensity)) + 1e-100;   //(avoid divide by zero)
        double peakHeight = peakIntegral *                      // estimate height minus background, assuming a Gaussian
                (1/(2*Math.PI*Math.max(Math.sqrt(peakMinorSqr*peakMajorSqr), 0.01*radius*radius)));

        double heightSignificance = peakHeight/backIntensitySigma;              //how much the maximum is above the background stddev
        double meanSignificance = peakIntegral/(sumPW*backIntensitySigma);      //how much the average in the integration area is above -"-
        double significance = Math.sqrt(heightSignificance*meanSignificance);   //we take the geometric mean of the two
        if (significance > 1e3) significance = 1e3;
        double pRadius = Math.sqrt(peakSigmaR*peakSigmaT);
        double pRadiusRatio = pRadius*(1./radius);
        // Heuristic correction factor for sigma from second moments in case that sigma
        // is not much smaller than the integration radius. With this correction, the error is
        // < 5% for a Gaussian spot with circular symmetry and radius > 0.7 sigma for CIRCLE and OVAL.
        // For AZIMUTH_BLUR, the correction depends on the ratio of major & minor radii. For simplicity,
        // we take the same correction as for 
        double radiusCorrectionFactor = 1;
        if (shape==OVAL)
            radiusCorrectionFactor = 1. + 0.001581/sqr(0.4925 - pRadiusRatio);
        else // CIRCLE or AZIMUTH_BLUR
            radiusCorrectionFactor = 1. + 0.001199/sqr(0.4744 - pRadiusRatio);
        peakSigmaR *= radiusCorrectionFactor;
        peakSigmaT  *= radiusCorrectionFactor;
        output[R_SIZE] = peakSigmaR;
        output[T_SIZE] = peakSigmaT;
        // Conditions for a valid position:
        // - The geometrical mean of INTEGRAL and estimated height must be 'minSignificance'
        //   above background stddev (noise and other background variations, including these
        //   due to the tails of the peak).
        // - Both eigenvalues of the 'sigma ellipse' must be positive (already checked above)
        //   and not too large, i.e., we must have a clear maximum, not a minimum in all directions.
        //   Note that a top-hat profile with the size of a circular integration disk would result in pMajor=0.5*radius.
        double pMajorMax = radius * 0.49;
        if (peakIntegral > 0 && peakMinorSqr > 0 && peakMajorSqr > 0 && peakMajorSqr < pMajorMax*pMajorMax
                && significance > NO_SIGNIFICANCE)
            output[PSIGNIFICANCE] = significance;
        else
            output[PSIGNIFICANCE] = 0.0;
        output[BACKGR] = offset;
        output[BSIGMA] = backIntensitySigma;
        if (significance >= minSignificance) {
            output[X] = xs + peakX;
            output[Y] = ys + peakY;
        }
        /*DEBUG if(Math.abs(xs-322.23)+Math.abs(ys-168.02)<0.2) {
        IJ.log(IJ.d2s(xs)+","+IJ.d2s(ys)+" + "+IJ.d2s(peakX)+","+IJ.d2s(peakY)+" int="+peakIntegral+" offs="+(float)offset);
        IJ.log("  sxx,xy,yy="+IJ.d2s(peakSigmaXsqr)+","+IJ.d2s(peakXY)+","+IJ.d2s(peakSigmaYsqr)+" offsX,Y="+IJ.d2s(peakX,1)+","+IJ.d2s(peakY,1)+" rScale="+IJ.d2s(radius*invRadiusT,4)+" cosR="+IJ.d2s(cosR,4)+" sinR="+IJ.d2s(sinR,4)+" sin¹+cos¹-1="+(float)(sinR*sinR+cosR*cosR-1));
        //IJ.log("raw rr,st,tt="+IJ.d2s(peakSigmaRsqr)+","+IJ.d2s(rtMixedRaw)+","+IJ.d2s(sigmaTsqrRaw)+" xyMajor²="+IJ.d2s(xyMajorSq,5)+" rtMajor²="+IJ.d2s(rtMajorRaw,5)+" rtMinor²="+IJ.d2s(rtMinorRaw,5)+" corrMajor²="+IJ.d2s(peakMajorSqr)+" Minor²"+IJ.d2s(peakMinorSqr,4));
        IJ.log("  rSiz="+IJ.d2s(output[R_SIZE])+" tSiz="+IJ.d2s(output[T_SIZE])+" rint="+IJ.d2s(radius)+" slope="+(float)xSlope+","+(float)ySlope);
        IJ.log("  p/b area="+IJ.d2s(sumPW,1)+","+IJ.d2s(sumBW,1)+" peakHeight="+(float)(peakHeight)+" sigmaBg="+(float)backIntensitySigma);
        IJ.log("  peakIntegral="+(float)peakIntegral+" minorSqr="+IJ.d2s(peakMinorSqr,1)+" majorSqr="+IJ.d2s(peakMajorSqr,1)+" pMajor="+IJ.d2s(Math.sqrt(peakMajorSqr),1)+" pMajorMax="+IJ.d2s(pMajorMax,1));
        IJ.log("  signif="+IJ.d2s(output[PSIGNIFICANCE])+" rawSig="+(float)significance+" heightSig="+IJ.d2s(heightSignificance)+" meanSig="+IJ.d2s(meanSignificance));
        }        /**/
        return output;
    }

    /** Returns the forground+background weight for a pixel dx, dy from the center of the integration circle.
     *  (to get the background weight alone, subtract the weight it would have in the foreground integration circle.
     *  Weights are 1 in the center, 0 outside, with a decay in a 1-pixel transition zone.
     *  The direction from the center of the LEED screen (needed for oval background) is given by cosR, sinR.
     *  The shape is encoded in invRadiusT:0=circle, -1=oval, positive values for AZIMUTH_BLUR */
    /* currently unused
        static double getTotalWeight(double dx, double dy, double cosR, double sinR, int shape, double bgRadius, double radiusT) {
        if (shape==CIRCLE)
            return getCircleWeight(dx, dy, bgRadius);
        else if (shape==OVAL)
            return getOvalWeight(dx, dy, cosR, sinR, bgRadius);
        else                    //if (shape==AZIMUTH_BLUR)
            return getOvalWeight(dx, dy, cosR, sinR, 0.5/bgRadius, 1./radiusT);
    } */

    /** Returns the weight for a pixel dx, dy from the center of an integration circle with given radius.
     *  Weights are 1 in the center, 0 outside, with a decay in a 1-pixel transition zone. */
    static double getCircleWeight(double dx, double dy, double radius) {
        double rSqr = dx*dx+dy*dy;
        if (rSqr >= sqr(radius+0.5))
            return 0;
        else if (rSqr <= sqr(radius-0.5))
            return 1;
        else
            return radius+0.5 - Math.sqrt(rSqr);
    }

    /** Returns the weight for a pixel dx, dy from the center of a background oval in tangential direction,
     *  with given semiminor axis 'radius' (in radial direction); the semimajor axis is twice as large.
     *  This is the outer border of the background area for OVAL background.
     *  Weights are 1 in the center, 0 outside, with a decay in a 1-pixel transition zone.
     *  The direction from the center of the LEED screen is given by cosR, sinR. */
    static double getOvalWeight(double dx, double dy, double cosR, double sinR, double radius) {
        double xp = sinR*dx - cosR*dy;          //point in rotated system
        double yp = cosR*dx + sinR*dy;
        double rSqr = xp*xp + 4*yp*yp;
        if (rSqr > sqr(2*radius+1)) return 0;   //for sure outside
        if (rSqr < sqr(2*radius-1)) return 1;   //for sure inside
        double corr = 0.5 * (1 - Math.abs(xp/(2*radius+0.5))); //correction to make transition zone 1 pxl wide everywhere
        if (corr < 0) corr = 0;
        double r = Math.sqrt(rSqr);
        if (r >= (2*radius+0.5) + corr)
            return 0;
        else if (r <= (2*radius-0.5) - corr)
            return 1;
        else
            return ((2*radius+0.5) + corr - r)/(1 + 2*corr);
    }

    /** Returns the weight for a pixel dx, dy from the center of an oval in tangential direction,
     *  with given semiminor axis (in radial direction) and semimajor axis in tangetial direction.
     *  'invRadiusR' is the reciprocal of the semiminor axis, and 'invRadiusT' the reciprocal of
     *  the semimajor axis. Note that the axis in radial direction MUST be smaller than that in
     *  tangential direction.
     *  Weights are 1 in the center, 0 outside, with a decay in an approx. 1-pixel-wide transition zone.
     *  The direction from the center defining the radial and tangential directions is given by cosR, sinR.
     *  The equation for the width of the transition zone is an approximation for elongated ovals:
     *  For an axis ratio of 3:1, the transition zone is between 0.9 and 1 pixels wide (0.82 to 1 pxl for 4:1).
     *  The accuracy is better for more circular ovals. */
    static double getOvalWeight(double dx, double dy, double cosR, double sinR, double invRadiusR, double invRadiusT) {
        double xp = (sinR*dx - cosR*dy)*invRadiusT;     //point in rotated & scaled system; x is tangential, y is radial
        double yp = (cosR*dx + sinR*dy)*invRadiusR;
        double rSqr = xp*xp + yp*yp;                    // 1 at the border
        if (rSqr > sqr(1 + 0.5*invRadiusR)) return 0;   // for sure outside
        if (rSqr < 1 - invRadiusR) return 1;            // for sure inside
        double xpNormSq = xp*xp/(xp*xp+yp*yp);
        // dr is the width of th etransition zone in scaled coordinates.
        // In the last term, xpNormSq^(3/2) would be more accurate, but slower:
        double dr = invRadiusR * (1 - (1.-invRadiusT/invRadiusR)*xpNormSq);
        if (rSqr >= sqr(1 + 0.5*dr))
            return 0.0;
        else if (rSqr <= sqr(1 - 0.5*dr))
            return 1.0;
        else
            return (0.5 + (1 - Math.sqrt(rSqr))/dr);
    }

    static double sqr(double x) {
        return x*x;
    }

    /** Returns two Rois with the shape of the integration area and background area,
     *  for showing it at the corner of the 'Spot tracking' stack as an overlay.
     *  The Roi is centered at (x, y).
     *  For the shape, it is assumed that the distance from the relevant center
     *  (screen center for OVAL, 0,0 spot for AZIMUTH_BLUR) is dx, dy in pixels;
     *  at least one of these must be nonzero.
     *  N.B. we use EllipseRois also for circles since circles (ovalRoi) have
     *  integer coordinates and would not appear concentric in many cases.
     */
    public static Roi[] getShapeRois(int x, int y, double dx, double dy,
            int shape, double radius, double azBlurRadians) {
        Roi intRoi = new EllipseRoi(x-radius, y, x+radius, y, 1.0);
        Roi bgRoi = null;
        double rho = Math.sqrt(dx*dx + dy*dy);
        double cosR = dx/rho;
        double sinR = dy/rho;
        switch (shape) {
            case CIRCLE:
                bgRoi = new EllipseRoi(x-Math.sqrt(2)*radius, y, x+Math.sqrt(2)*radius, y, 1.0);
                break;
            case OVAL:
                bgRoi = new EllipseRoi(x-2*radius*sinR, y+2*radius*cosR, x+2*radius*sinR, y-2*radius*cosR, /*aspectRatio=*/0.5);
                break;
            case AZIMUTH_BLUR:
                double radiusT = Math.sqrt(sqr(radius) + sqr(azBlurRadians)*(dx*dx + dy*dy));
                double radiusTbg = Math.max(Math.sqrt(2)*radius, radiusT);
                intRoi = new EllipseRoi(x-radiusT*sinR, y+radiusT*cosR, x+radiusT*sinR, y-radiusT*cosR, radius/radiusT);
                bgRoi = new EllipseRoi(x-radiusTbg*sinR, y+radiusTbg*cosR, x+radiusTbg*sinR, y-radiusTbg*cosR, radius*Math.sqrt(2)/radiusTbg);
                break;
        }
        return new Roi[] {intRoi, bgRoi};
    }

    /** Constructor for tests only: Creates an image with the response to
     *  a pixel (delta function, unit pulse) at different positions
     *  (Use with simple Compile&Run in ImageJ).
     *  Comment out when not in use. */ /*
    public LeedSpotAnalyzer() {
        final int SIZE = 50;                    // image size
        final int OUTTYPE = BACKGR;
        final int SHAPE = AZIMUTH_BLUR;//OVAL;//
        final double RADIUS = 10;
        final double AZBLURRADIANS = 0.05;     //az blur angle in rad
        double[] output = new double[OUTPUT_SIZE];
        FloatProcessor fpIn = new FloatProcessor(SIZE, SIZE);
        FloatProcessor fpOut = new FloatProcessor(SIZE, SIZE);
        ByteProcessor mask = new ByteProcessor(SIZE, SIZE);
        mask.setRoi(0, 0, SIZE, SIZE);
        mask.setValue(255);
        mask.fill();
        for (int y=0; y<SIZE; y++) {
            for (int x=0; x<SIZE; x++) {
                // int x=20, y=20;{{ //alternatively: no loop and write values in analyzeSpot: 
                fpIn.putPixelValue(x, y, 1000);
                output = analyzeSpot(25, 25, 421, 75, SHAPE, RADIUS, AZBLURRADIANS, 2., fpIn, mask, output);
                fpIn.putPixelValue(x, y, 0);
                fpOut.putPixelValue(x, y, output[OUTTYPE]);
            }
        }
        (new ImagePlus(DATA_NAMES[OUTTYPE], fpOut)).show();
    }/**/
}
