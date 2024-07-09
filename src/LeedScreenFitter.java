import ij.*;
import ij.IJ;
import ij.gui.*;
import ij.process.*;
/**
 *  Fits a 2D polynomial for the conversion from kx, ky to screen coordinates
 *
 *  Note that there are primary fit types (levels, linear to 5th order, up to
 *  HIGHEST_STD_FUNCTION) and alternative fit types for the levels "square"
 *  (=2nd order) and 4th order. These alternative fit types have fewer
 *  parameters but include higher order terms and may yield a better fit
 *  if the distortions only depend on the radius [e.g., (0,0) spot in the
 *  center and camera close to the screen or sample surface not in the center
 *  of the screen curvature.]
 *
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

public class LeedScreenFitter {
    public static final int LINEAR=0, SQUARE=1, CUBIC=2, FOURTH_ORDER=3, FIFTH_ORDER=4, LIN_EDGE_DISTORT=5, CUBE_EDGE_DISTORT=6; //function types

    public static final int[] N_PARAM = new int[]       // number of parameters per direction for each function type
            {3, 6, 10, 15, 21, 5, 12};
    public static final int[] ALTERNATIVE = new int[]   // alternative fit type using the same number or fewer parameters
            {-1, LIN_EDGE_DISTORT, -1, CUBE_EDGE_DISTORT, -1, -1, -1};
    public static final int HIGHEST_STD_FUNCTION = FIFTH_ORDER; // all higher ones are alternative types
    static final int[] MAIN_TYPE = new int[]            // main fit type for each function type
            {LINEAR, SQUARE, CUBIC, FOURTH_ORDER, FIFTH_ORDER, SQUARE, FOURTH_ORDER};
    static final String[] FUNCTION_NAMES = new String[] // fit names
            {"linear", "square", "cubic", "4th order", "5th order", "linear+r^3", "cubic+r^5"};
    static final double MAX_NONLINEARITY_SQR = 10;      // (x^2+y^2) must not be less than 1/10 of that from linear
    static final int ARRAY_OFFSET = 6;                  // in toArray, fromArray: up to 6 initial values before parameter arrays
    int functionType;                                   // which function we use to translate indices to screen coordinates
    double maskCenterX, maskCenterY;                    // internally, x & y are relative to this
    double[] xParams, yParams;                          // fit parameters for x, y: offest, 2 linear terms in kx, ky, 2nd-order, ...
    double kFactor;                                     // k must be multiplied by this to enter fit function
    double invSqrtEnergy;                               // 1/sqrt of the energy where the fit was done
    int nFitIndices;                                    // how many spot indices we have used for the fit
 	private String statusString;						// e.g. "cubic, delta rms 1.6 pxl, 50 spots"

    /** Constructor, sets the function type (LINEAR, etc., the energy, and the screen center */
    public LeedScreenFitter(int functionType, double energy, double maskCenterX, double maskCenterY) {
        this.functionType = functionType;
        this.invSqrtEnergy = 1.0/Math.sqrt(energy);
        this.maskCenterX = maskCenterX;
        this.maskCenterY = maskCenterY;
    }

    /** Constructor setting all variables and parameters from one array of the following:
     *  functionType, invSqrtEnergy, kFactor, maskCenterX, maskCenterY, <spare>, xParams, yParams */
    public LeedScreenFitter(double[] a) {
        this.functionType = (int)a[0];
        int nParam = N_PARAM[functionType];
        if (a.length < ARRAY_OFFSET+2*nParam)
            throw new RuntimeException(LEED_Spot_Tracker.PLUGIN_NAME+": Error reading screen basis\nArray too short ("+a.length+")");
        this.kFactor = a[1];
        this.invSqrtEnergy = a[2];
        this.maskCenterX = a[3];
        this.maskCenterY = a[4];
        //a[5] spare
        this.xParams = new double[nParam];
        this.yParams = new double[nParam];
        System.arraycopy(a, ARRAY_OFFSET, xParams, 0, nParam);
        System.arraycopy(a, ARRAY_OFFSET+nParam, yParams, 0, nParam);
    }

    /** Creates one array with all variables and parameters:
     *  functionType, invSqrtEnergy, kFactor, maskCenterX, maskCenterY, <spare>, xParams, yParams */
    public double[] toArray() {
        int nParam = N_PARAM[functionType];
        double[] a = new double[ARRAY_OFFSET+2*nParam];
        a[0] = functionType;
        a[1] = kFactor;
        a[2] = invSqrtEnergy;
        a[3] = maskCenterX;
        a[4] = maskCenterY;
        //a[5] spare
        System.arraycopy(xParams, 0, a, ARRAY_OFFSET, nParam);
        System.arraycopy(yParams, 0, a, ARRAY_OFFSET+nParam, nParam);
        return a;
    }

    /** Fits a given function type to a number of spots with energy-corrected reciprocal-space coordinates kx, ky,
     *  and screen coordinates x, y.
     *  For a linear fit, if there are not enough points given or they are collinear, tries to add points to obtain
     *  a fit: The (0,0) spot is assumed to be in the mask center position, and it is assumed that the
     *  screen coordinates are a scaled and rotated version of the reciprocal-space coordinates.
     *  Returns false on error (singular matrix), and saves the result in instance variables xParams, yParams, kFactor.
     *  */
    public boolean fitSpots(double[] kx, double[] ky, double[] x, double[] y) {
        if (functionType > LINEAR && x.length < N_PARAM[functionType])
            return false;                // would be overfitting; we try enhancing the input for linear fits only

        // add to spots if not enough for a linear fit
        int nSpots = kx.length < 3 ? 3 : kx.length;
        boolean allCollinear = allCollinear(kx, ky);
        if (allCollinear && kx.length >= 3)             // at least 3 spots, but not enough?
		nSpots = kx.length + 1;
        if (nSpots > kx.length) {						// add one or two spots to create a basis
            double[] kxN = new double[nSpots];
            double[] kyN = new double[nSpots];
            double[] xN = new double[nSpots];
            double[] yN = new double[nSpots];
            System.arraycopy(kx, 0, kxN, 0, kx.length);
            System.arraycopy(ky, 0, kyN, 0, kx.length);
            System.arraycopy(x, 0,  xN, 0,  kx.length);
            System.arraycopy(y, 0,  yN, 0,  kx.length);
			int i00spot = -1;                       	//was the (0,0) spot entered?
			for (int i=0; i<kx.length; i++) {
				if (kx[i] == 0 && ky[i] == 0) {
					i00spot = i;
					break;
				}
			}
            if (allCollinear) {
                double highestKsqr = 0;
                int iOfHighestK = 0;
                for (int i=0; i<kx.length; i++) {
                    double kSqr = sqr(kx[i]) + sqr(ky[i]);
                    if (kSqr > highestKsqr) {
                        highestKsqr = kSqr;
                        iOfHighestK = i;
                    }
                }
                double xCenter = i00spot<0 ? maskCenterX : x[i00spot]; // if we have a (0,0) spot, use it as the center for getting the 3rd spot
                double yCenter = i00spot<0 ? maskCenterY : y[i00spot];

                int iNewSpot = x.length;
                kxN[iNewSpot] = -ky[iOfHighestK];       // add spot rotated 90deg, assuming that one system is right-handed
                kyN[iNewSpot] =  kx[iOfHighestK];       // and one left-handed (screen coordintes have y down)
                xN[iNewSpot] = xCenter + (y[iOfHighestK] - yCenter);
                yN[iNewSpot] = yCenter - (x[iOfHighestK] - xCenter);
                //IJ.log("center="+(int)maskCenterX+","+(int)maskCenterY+" add kxy="+IJ.d2s(kxN[iNewSpot], 3,5)+","+IJ.d2s(kyN[iNewSpot], 3,5)+" at "+(int)xN[iNewSpot]+","+(int)yN[iNewSpot]);
            }
            int lastI = kxN.length-1;
            if (kxN[lastI] == 0 && kyN[lastI] == 0) {   // last spot not known?
				if (i00spot >= 0) return false;
                xN[lastI] = maskCenterX;				// add suspected spot: (0,0) in the screen center
                yN[lastI] = maskCenterY;
            }
            kx = kxN; ky = kyN;
            x = xN;   y = yN;
        }

        // The fit actually uses kx*kFactor, ky*kFactor; these products are, on average, numbers close to one;
        // this makes the simple numeric fit procedure (using the normal equations) more stable.
        double sumKsqr = 0;
        for (int i = 0; i < nSpots; i++)
            sumKsqr += sqr(kx[i]) + sqr(ky[i]);
        sumKsqr *= 1./(2*nSpots);                       //mean square k value
        kFactor = 1./Math.sqrt(sumKsqr*invSqrtEnergy);  //normalization factor to avoid numeric problems in matrix inversion
        //IJ.log("FIT: type="+functionType+" nSpots="+nSpots+" kFact="+kFactor+" 1/sqrtE="+invSqrtEnergy);

        int nParam = N_PARAM[functionType];
        xParams = new double[nParam];
        yParams = new double[nParam];
        double[][] matrix = new double[nParam][nParam];
        double[] xVec = new double[nParam];
        double[] yVec = new double[nParam];
        double[] functP = new double[nParam];

        for (int i=0; i < kx.length; i++) {
            //IJ.log("point kx,ky="+kx[i]+","+ky[i]+" x,y="+(int)x[i]+","+(int)y[i]);
            for (int iPar=0; iPar<nParam; iPar++) {
                    double fct = functionPart(functionType, iPar, kx[i]*(kFactor*invSqrtEnergy), ky[i]*(kFactor*invSqrtEnergy));
                    functP[iPar] = fct;
                    xVec[iPar] += (x[i] - maskCenterX)*fct;
                    yVec[iPar] += (y[i] - maskCenterY)*fct;
                    matrix[iPar][iPar] += fct*fct;
                    for (int jPar=0; jPar<iPar; jPar++) {
                            double fctProd = fct*functP[jPar];
                            matrix[iPar][jPar] += fctProd;
                            matrix[jPar][iPar] += fctProd;
                    }
            }
        }
        boolean ok = LeedUtils.invertSymmetricMatrix(matrix);
        if (!ok) return false;

        for (int i=0; i<nParam; i++)                    //calculate fit parameters
            for (int j=0; j<nParam; j++) {
                    xParams[i] += matrix[i][j]*xVec[j];
                    yParams[i] += matrix[i][j]*yVec[j];
        }
        //IJ.log("function "+functionType+" residuals="+IJ.d2s(getRmsResiduals(kx, ky, x, y)));
        //for (int i=0; i<N_PARAM[functionType]; i++)IJ.log("xyParam["+i+"]="+IJ.d2s(xParams[i], 3, 5)+",  "+IJ.d2s(yParams[i], 3, 5));
        nFitIndices = nSpots;
        return true;
    }

    /** Returns whether all points are on a line or two points are identical. */
    boolean allCollinear(double[] kx, double[] ky) {
        for (int i1=0; i1<kx.length; i1++)
            for (int i2=i1+1; i2<kx.length; i2++)
				for (int i3=i2+1; i3<kx.length; i3++)
					if (!collinear(kx[i2]-kx[i1], ky[i2]-ky[i1], kx[i3]-kx[i2], ky[i3]-ky[i2]))
						return false;
        return true;
    }

	/** Returns whether two vectors are collinear */
	boolean collinear(double x1, double y1, double x2, double y2) {
		return sqr(x1*y2 - x2*y1) < 1e-8*(sqr(x1) + sqr(y1))*(sqr(x2) + sqr(y2)); //cross product << product of lengths
	}

    /** Returns the rms distance between fit and measured spot positions, for the energy where the fit was done */
    double getRmsResiduals(double[] kx, double[] ky, double[] x, double[] y) {
        double sumResidualSqr = 0;
        double[] xy = new double[2];
        for (int i=0; i < kx.length; i++) {
            xy = screenCoordinates(kx[i], ky[i], invSqrtEnergy, xy);
            sumResidualSqr += sqr(xy[0] - x[i]) + sqr(xy[1] - y[i]);
        }
        double residuals = Math.sqrt(sumResidualSqr/kx.length);
        setStatusString(residuals);
        return residuals;
    }

    /** Returns the rms distance between fit and measured spot positions, using only spots with indices >=0,
     *  for the energy where the fit was done */
    double getRmsResiduals(LeedSpotPattern spotPattern, int[] indices, double[] x, double[] y) {
        double sumResidualSqr = 0;
        double[] xy = new double[2];
        int nSpots=0;
        for (int i=0; i < indices.length; i++) {
            int index = indices[i];
            if (index >= 0) {
                xy = screenCoordinates(spotPattern.getKx(index), spotPattern.getKy(index), invSqrtEnergy, xy);
                sumResidualSqr += sqr(xy[0] - x[i]) + sqr(xy[1] - y[i]);
                nSpots++;
            }
        }
        double residuals = Math.sqrt(sumResidualSqr/nSpots);
        setStatusString(residuals);
        return residuals;
    }

	/** Returns information on the fit, e.g. "cubic, delta rms 1.6 pxl, 50" <add " spots"> */
	public String getStatusText() {
		return statusString;
	}

	private void setStatusString(double residuals) {
		statusString = getFunctionName() + ", \u0394rms="+IJ.d2s(residuals,1) + "pxl, " + getNFitIndices();
	}

    /** Returns the number of indices (spots) used for the fit */
    int getNFitIndices(){
        return this.nFitIndices;
    }

    /** Returns the calculated screen coordinates x, y for a given spot
     *  and 1/square root of the energy. The array xy, if not null, will hold the result.
     *  If the xy array has a length of 4 or more, writes dx/d ln k and dy/d ln k into
     *  xy[2] and xy[3] */
    public double[] screenCoordinates(LeedSpotPattern spotPattern, int spotIndex, double invSqrtEnergy, double[] xy) {
        double kx = spotPattern.getKx(spotIndex);
        double ky = spotPattern.getKy(spotIndex);
        return screenCoordinates(kx, ky, invSqrtEnergy, xy);
    }

    /** Returns the smallest scale factor in x&y of the linear component at the energy of the fit */
    public double kToScreenScale() {
        return kFactor*invSqrtEnergy*Math.sqrt(Math.min(sqr(xParams[1])+sqr(yParams[1]), sqr(xParams[2])+sqr(yParams[2])));
    }

    /** Returns the function type used; for alternative types, returns the main type */
    public int getMainFunctionType() {
        return MAIN_TYPE[functionType];
    }

    /** Returns the function type used; may be an alternative type */
    public int getFunctionType() {
        return functionType;
    }

    /** Returns the name of the fit function*/
    public String getFunctionName() {
        return FUNCTION_NAMES[functionType];
    }

    /** Returns the calculated screen coordinates x, y for given reciprocal-space (cartesian) coordinates
     *  and 1/square root of the energy. The xy array for the output must be provided.
     *  If the xy array has a length of 4 or more, writes dx/d ln k and dy/d ln k into xy[2] and xy[3].
     *  Implementation note: The brackets for grouping products allow the compiler to better optimize
     *  (common subexpressions, but better readability than Horner's method) */
    public double[] screenCoordinates(double kx, double ky, double invSqrtEnergy, double[] xy) {
        kx *= kFactor*invSqrtEnergy;
        ky *= kFactor*invSqrtEnergy;
        double x = kx*xParams[1] + ky*xParams[2];   //offsets will be added at the very end
        double y = kx*yParams[1] + ky*yParams[2];
        double linSqr = x*x + y*y;
        if (functionType == LIN_EDGE_DISTORT) {
            double ksqr = sqr(kx) + sqr(ky);
            x += ksqr*(kx*xParams[3] + ky*xParams[4]);
            y += ksqr*(kx*yParams[3] + ky*yParams[4]);
        } else if (functionType >= SQUARE) {
            x += kx*kx*xParams[3] + kx*ky*xParams[4] + ky*ky*xParams[5];
            y += kx*kx*yParams[3] + kx*ky*yParams[4] + ky*ky*yParams[5];
        }
        if (functionType == CUBIC || functionType == FOURTH_ORDER || functionType == FIFTH_ORDER ||
				functionType == CUBE_EDGE_DISTORT) {
            x += kx*kx*kx*xParams[6] + kx*kx*ky*xParams[7] + kx*(ky*ky)*xParams[8] + ky*(ky*ky)*xParams[9];
            y += kx*kx*kx*yParams[6] + kx*kx*ky*yParams[7] + kx*(ky*ky)*yParams[8] + ky*(ky*ky)*yParams[9];
        }
        if (functionType == FOURTH_ORDER || functionType == FIFTH_ORDER) {
            x += kx*kx*(kx*kx)*xParams[10] + kx*kx*(kx*ky)*xParams[11] + kx*kx*(ky*ky)*xParams[12] +
                    kx*ky*(ky*ky)*xParams[13] + ky*ky*(ky*ky)*xParams[14];
            y += kx*kx*(kx*kx)*yParams[10] + kx*kx*(kx*ky)*yParams[11] + kx*kx*(ky*ky)*yParams[12] +
                    kx*ky*(ky*ky)*yParams[13] + ky*ky*(ky*ky)*yParams[14];
        } else if (functionType == CUBE_EDGE_DISTORT) {
            double ksqr2 = sqr(sqr(kx) + sqr(ky));
            x += ksqr2*(kx*xParams[10] + ky*xParams[11]);
            y += ksqr2*(kx*yParams[10] + ky*yParams[11]);
        }
        if (functionType == FIFTH_ORDER) {
            x += kx*(kx*kx)*(kx*kx)*xParams[15] + (kx*kx)*(kx*kx)*ky*xParams[16] + (kx*kx)*(kx*(ky*ky))*xParams[17] +
                    (kx*kx)*(ky*(ky*ky))*xParams[18] + (kx*(ky*ky))*(ky*ky)*xParams[19] + ky*(ky*ky)*(ky*ky)*xParams[20];
            y += kx*(kx*kx)*(kx*kx)*yParams[15] + (kx*kx)*(kx*kx)*ky*yParams[16] + (kx*kx)*(kx*(ky*ky))*yParams[17] +
                    (kx*kx)*(ky*(ky*ky))*yParams[18] + (kx*(ky*ky))*(ky*ky)*yParams[19] + ky*(ky*ky)*(ky*ky)*yParams[20];
        }
		//consistency check for debugging
		/*  if (true) {
            double xt=0,yt=0, xm=0, ym=0;
		    for (int i=0;i<N_PARAM[functionType];i++) {
			  double xp=functionPart(functionType,i,kx,ky)*xParams[i];
			  xt+=xp;
			  double yp=functionPart(functionType,i,kx,ky)*yParams[i];
			  yt+=yp;
			  if (Math.abs(xp)>xm) xm=Math.abs(xp); if (Math.abs(yp)>ym) ym=Math.abs(yp);
		    }
		    if(Math.abs(x+xParams[0]-xt)>1e-3)IJ.log(FUNCTION_NAMES[functionType]+" x="+(x+xParams[0])+" vs "+xt+" max="+xm);
		    if(Math.abs(y+yParams[0]-yt)>1e-3)IJ.log(FUNCTION_NAMES[functionType]+" y="+(y+yParams[0])+" vs "+yt+" max="+ym);
		  } */

        // Calculation of gradients: used to check that the function is not too far from nonlinear, which would result in
        // putting high-order spots to near the center and does not fold back
        // Another check for folding back is the deviation from the linear terms, see MAX_NONLINEARITY_SQR
        // If required, here we also calculate dx/d ln k, dy/d ln k
        boolean calculate_dxy_dlnk = xy.length >= 4;

        if (kx == 0 && ky == 0) {
            if (calculate_dxy_dlnk) {
                xy[2] = 0; xy[3] = 0;
            }
        } else if (functionType > LINEAR && (x*x+y*y)*MAX_NONLINEARITY_SQR < linSqr) {
            x = Double.NaN;    // x, y much less than what the linear coefficients say? Then it is invalid
            y = Double.NaN;
        } else if (functionType > LINEAR || calculate_dxy_dlnk) {
            double dxdkx = xParams[1], dxdky = xParams[2];   //derivatives dx/dkx etc.
            double dydkx = yParams[1], dydky = yParams[2];
            if (functionType == LIN_EDGE_DISTORT) {
                dxdkx += (kx*kx*3+ky*ky)*xParams[3] + kx*ky*2*xParams[4];   //derivative of (kx^3 + kx*ky^2) & (kx^2*ky + ky^3)
                dxdky += kx*ky*2*xParams[3] + (kx*kx+ky*ky*3)*xParams[4];
                dydkx += (kx*kx*3+ky*ky)*yParams[3] + kx*ky*2*yParams[4];
                dydky += kx*ky*2*yParams[3] + (kx*kx+ky*ky*3)*yParams[4];
            } else if (functionType >= SQUARE) {
                dxdkx += kx*2*xParams[3] + ky*xParams[4];                   //derivative of kx*kx & kx*ky & ky*ky
                dxdky += kx*xParams[4] + ky*2*xParams[5];
                dydkx += kx*2*yParams[3] + ky*yParams[4];
                dydky += kx*yParams[4] + ky*2*yParams[5];
            }
            if (functionType == CUBIC || functionType == FOURTH_ORDER || functionType == CUBE_EDGE_DISTORT || functionType == FIFTH_ORDER) {
                dxdkx += kx*kx*3*xParams[6] + kx*ky*2*xParams[7] + ky*ky*xParams[8]; //derivative of kx^3 & kx^2*ky & kx*ky*^2 & ky^3
                dxdky += kx*kx*xParams[7] + kx*ky*2*xParams[8] + ky*ky*3*xParams[9];
                dydkx += kx*kx*3*yParams[6] + kx*ky*2*yParams[7] + ky*ky*yParams[8];
                dydky += kx*kx*yParams[7] + kx*ky*2*yParams[8] + ky*ky*3*yParams[9];
            }
            double dxdkx3 = dxdkx, dxdky3 = dxdky; // derivatives using fit terms up to 3rd order
            double dydkx3 = dydkx, dydky3 = dydky;
            if (functionType == FOURTH_ORDER || functionType == FIFTH_ORDER) { //derivative of kx^4 & kx^3*ky & kx^2*ky^2 & kx*ky^3 & ky^4
                dxdkx += kx*kx*kx*4*xParams[10] + kx*kx*ky*3*xParams[11] + kx*(ky*ky)*2*xParams[12] + ky*ky*ky*xParams[13];
                dxdky += kx*kx*kx*xParams[11] + kx*kx*ky*2*xParams[12] + kx*(ky*ky)*3*xParams[13] + ky*ky*ky*4*xParams[14];
                dydkx += kx*kx*kx*4*yParams[10] + kx*kx*ky*3*yParams[11] + kx*(ky*ky)*2*yParams[12] + ky*ky*ky*yParams[13];
                dydky += kx*kx*kx*yParams[11] + kx*kx*ky*2*yParams[12] + kx*(ky*ky)*3*yParams[13] + ky*ky*ky*4*yParams[14];
            } else if (functionType == CUBE_EDGE_DISTORT) { //derivative of kx^5 + 2 kx^3*ky^2 + kx*ky^4 & kx^4*ky + 2 kx^2*ky^3 + ky^5
                dxdkx += (kx*kx*(kx*kx)*5 + kx*kx*(ky*ky)*(2*3) + ky*ky*(ky*ky))*xParams[10] +
                        (kx*kx*kx*ky + kx*ky*(ky*ky))*4*xParams[11];
                dxdky += (kx*kx*(kx*ky) + kx*ky*(ky*ky))*4*xParams[10] +
                        (kx*kx*(kx*kx) + kx*kx*(ky*ky)*(2*3) + ky*ky*(ky*ky)*5)*xParams[11];
                dydkx += (kx*kx*(kx*kx)*5 + kx*kx*(ky*ky)*(2*3) + ky*ky*(ky*ky))*yParams[10] +
                        (kx*kx*(kx*ky) + kx*ky*(ky*ky))*4*yParams[11];
                dydky += (kx*kx*(kx*ky) + kx*ky*(ky*ky))*4*yParams[10] +
                        (kx*kx*(kx*kx) + kx*kx*(ky*ky)*(2*3) + ky*ky*(ky*ky)*5)*yParams[11];
            }
            double dxdkx4 = dxdkx, dxdky4 = dxdky; // derivatives using fit terms up to 4th order
            double dydkx4 = dydkx, dydky4 = dydky;
            if (functionType == FIFTH_ORDER) { //derivative of kx^5 & kx^4*ky & kx^3*ky^2 & kx^2*ky^3 & kx*ky^4 & ky^5 [15-20]
				dxdkx += kx*kx*(kx*kx)*5*xParams[15] + kx*kx*(kx*ky)*4*xParams[16] + kx*kx*(ky*ky)*3*xParams[17] + kx*ky*(ky*ky)*2*xParams[18] + ky*ky*(ky*ky)*xParams[19];
                dxdky += kx*kx*(kx*kx)*xParams[16] + kx*kx*(kx*ky)*2*xParams[17] + kx*kx*(ky*ky)*3*xParams[18] + kx*ky*(ky*ky)*4*xParams[19] + ky*ky*(ky*ky)*5*xParams[20];
				dydkx += kx*kx*(kx*kx)*5*yParams[15] + kx*kx*(kx*ky)*4*yParams[16] + kx*kx*(ky*ky)*3*yParams[17] + kx*ky*(ky*ky)*2*yParams[18] + ky*ky*(ky*ky)*yParams[19];
                dydky += kx*kx*(kx*kx)*yParams[16] + kx*kx*(kx*ky)*2*yParams[17] + kx*kx*(ky*ky)*3*yParams[18] + kx*ky*(ky*ky)*4*yParams[19] + ky*ky*(ky*ky)*5*yParams[20];
            }
            /*if (true) {
				//consistency test of gradient for debugging
				double deltaX=-2e-7, deltaY=1e-7;
				double dx=dxdkx*deltaX + dxdky*deltaY;
				double dy=dydkx*deltaX + dydky*deltaY;
				double[] sc2=screenCoordinates((kx+deltaX)/(kFactor*invSqrtEnergy), (ky+deltaY)/(kFactor*invSqrtEnergy), invSqrtEnergy, new double[2]);
				double dxApprox=sc2[0]-maskCenterX-xParams[0]-x;
				double dyApprox=sc2[1]-maskCenterY-yParams[0]-y;
				if (Math.abs(dx-dxApprox)>1e-3*Math.abs(dx) && Math.abs(dx-dxApprox)>1e-14)
					IJ.log(FUNCTION_NAMES[functionType]+" dx="+dx+" approx="+dxApprox);
				if (Math.abs(dy-dyApprox)>1e-3*Math.abs(dy) && Math.abs(dy-dyApprox)>1e-14)
					IJ.log(FUNCTION_NAMES[functionType]+" dy="+dy+" approx="+dyApprox);
			}*/

            if (functionType > LINEAR && !(kx == 0 && ky == 0)) {
                // in the following we care about the movement becoming reversed, which can occur for high kx, ky, where the function folds back into the screen area
                double xLin = kx*xParams[1] + ky*xParams[2];    //direction from the center with linear terms;
                double yLin = kx*yParams[1] + ky*yParams[2];    //  the spot should move close to this direction with increasing |k|
                double xMove = dxdkx*kx + dxdky*ky;             //the spot moves in this direction with increasing |k|
                double yMove = dydkx*kx + dydky*ky;
                if (!collinearWithin30deg(xLin, yLin, xMove, yMove)) {
                    x = Double.NaN;                             //invalid result, position and motion not collinear within 30 degrees
                    y = Double.NaN;
                //DEBUG IJ.log("invald: kx,ky="+IJ.d2s(kx,5)+","+IJ.d2s(ky,5)+" xyLin="+(float)xLin+","+(float)yLin+" move3="+(float)xMove+","+(float)yMove);
                } else if (functionType > CUBIC && functionType != LIN_EDGE_DISTORT) {
                    // In these cases the gradient may be reversed twice, so folding back & forth can happen.
                    // We therefore check for folding back of the 3rd order terms (which fold back where
                    // higher orders can go forth again). This assumes that higher orders are small corrections
                    // and we don't eliminate false positives
                    double xMove3 = dxdkx3*kx + dxdky3*ky;      //based on therms up to 3rd order, the spot moves in this direction with increasing |k|
                    double yMove3 = dydkx3*kx + dydky3*ky;
                    if (!collinearWithin30deg(xLin, yLin, xMove3, yMove3)) { //if (xLin*xMove3 + yLin*yMove3 <= 0) {
                        x = Double.NaN;                         //invalid result, 3rd-order terms
                        y = Double.NaN;
                        //DEBUG IJ.log("invald: kx,ky="+IJ.d2s(kx,5)+","+IJ.d2s(ky,5)+" xyLin="+(float)xLin+","+(float)yLin+" move3="+(float)xMove3+","+(float)yMove3);
                    }
                    if (functionType == FIFTH_ORDER) {          //also check for 4th-order terms folding back
                        double xMove4 = dxdkx4*kx + dxdky4*ky;  //based on therms up to 4th order, the spot moves in this direction with increasing |k|
                        double yMove4 = dydkx4*kx + dydky4*ky;
                        if (!collinearWithin30deg(xLin, yLin, xMove4, yMove4)) { //if (xLin*xMove4 + yLin*yMove4 <=0 ) {
                            x = Double.NaN;                     //invalid result, 4rd-order fit function folds back
                            y = Double.NaN;
                            //DEBUG IJ.log("fold back kx,y="+(float)kx+","+(float)ky+" lin="+(float)xLin+","+(float)yLin+" mov4="+(float)xMove4+","+(float)yMove4);
                        }
                    }
                }
            }
            if (calculate_dxy_dlnk) {
                xy[2] = dxdkx * kx + dxdky * ky;
                xy[3] = dydkx * kx + dydky * ky;
            }
        }
        if (xy == null)
            xy = new double[2];
        xy[0] = x + xParams[0] + maskCenterX;
        xy[1] = y + yParams[0] + maskCenterY;
        return xy;
    }

    /** Returns the fit parameters in x direction. */
    /*  currently unused
    public double[] getXparameters() {
        return xParams;
    }  */

    /** Returns the fit parameters in y direction. */
    /*  currently unused
    public double[] getYparameters() {
        return yParams;
    } */

    /** Returns the linear fit parameters for the given 1/sqrt(E) as a 4-element array: 
     *  dx/dkx, dy/dkx, dx/dky, dy/dky */
    public double[] getLinearMatrix(double invSqrtEnergy) {
        return new double[] {kFactor*invSqrtEnergy*xParams[1], kFactor*invSqrtEnergy*yParams[1],
                kFactor*invSqrtEnergy*xParams[2], kFactor*invSqrtEnergy*yParams[2]};
    }

    /** Returns the x position of the (0,0) spot in screen coordinates */
    public double getX00() {
        return xParams[0] + maskCenterX;
    }

    /** Returns the y position of the (0,0) spot in screen coordinates */
    public double getY00() {
        return yParams[0] + maskCenterY;
    }

    /** Returns the fit function contribution for a given parameter number
     *  and kx*sqrt(energy), ky*sqrt(energy). */
    static double functionPart(int functionType, int iParam, double kx, double ky) {
        switch(iParam) {
            case 0: return 1;
            case 1: return kx;
            case 2: return ky;
            case 3: return functionType == LIN_EDGE_DISTORT ? (sqr(kx) + sqr(ky))*kx :
                    kx*kx;
            case 4: return functionType == LIN_EDGE_DISTORT ? (sqr(kx) + sqr(ky))*ky :
                    kx*ky;
            case 5: return ky*ky;
            case 6: return kx*kx*kx;
            case 7: return kx*kx*ky;
            case 8: return kx*ky*ky;
            case 9: return ky*ky*ky;
            case 10: return functionType == CUBE_EDGE_DISTORT ? sqr(sqr(kx) + sqr(ky))*kx :
                    kx*kx*(kx*kx);
            case 11: return functionType == CUBE_EDGE_DISTORT ? sqr(sqr(kx) + sqr(ky))*ky :
                    kx*kx*(kx*ky);
            case 12: return kx*kx*(ky*ky);
            case 13: return kx*ky*(ky*ky);
            case 14: return ky*ky*(ky*ky);
            case 15: return (kx*kx)*(kx*kx)*kx;
            case 16: return (kx*kx)*(kx*kx)*ky;
            case 17: return (kx*kx)*kx*(ky*ky);
            case 18: return (kx*kx)*ky*(ky*ky);
            case 19: return kx*(ky*ky)*(ky*ky);
            case 20: return ky*(ky*ky)*(ky*ky);
            default: throw new RuntimeException("Invalid parameter number: "+iParam);
        }
    }

    /** Distance squared between two points in 2D space */
    double distanceSqr(double[] xy1, double[] xy2) {
        return sqr(xy1[0] - xy2[0]) + sqr(xy1[1] - xy2[1]);
    }

    /** Returns whether two non-zero vectors are collinear within 30 degrees.
     *  Returns false if one or both vectors have length zero.
     *  Used for checking whether the mapping folds back because of high nonlinearity
     *  far from the center (where the mapping is not valid any more)
     *  Not that even for 20° from perpendicular incidence, the spots still move towards
     *  the (0,0) spot within less than 30 deg deviation **/
    static boolean collinearWithin30deg(double x1, double y1, double x2, double y2) {
        double innerProduct = x1*x2 + y1*y2;
        double len1sqr = x1*x1 + y1*y1;
        double len2sqr = x2*x2 + y2*y2;
        return innerProduct > 0 &&
                sqr(innerProduct) >= 0.75 * len1sqr * len2sqr; //sqr(cos(30)) = 0.75
    }

    static double sqr(double x) {return x*x;}
}
