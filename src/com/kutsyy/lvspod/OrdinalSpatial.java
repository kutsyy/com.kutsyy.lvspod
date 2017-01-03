package com.kutsyy.lvspod;

import com.kutsyy.util.*;
import java.io.*;

/**
 *  Main class of the libarary.<BR>
 *   Created by <A
 *  href="http://www.kutsyy.com">Vadim Kutsyy</A> <BR>
 *
 * This class should be called for all computations<BR>
 *
 *@author     <A href="http://www.kutsyy.com">Vadim Kutsyy</A>
 *@created    December 1, 2000
 */
public final class OrdinalSpatial extends OrdinalSpatialAbstract {
    /**
     *  Description of the Field
     */
    public OrdinalSpatialQL[] TNearestAvg;
    private String Type = "median";
    private OrdinalSpatialAbstract T;
    private OrdinalSpatialMdA TMdA;
    private OrdinalSpatialMnA TMnA;
    private OrdinalSpatialCLTrueXPLE TTrueXPLE;
    private OrdinalSpatialCLTrueX TTrueX;
    /**
     *  Insert the method's description here. Creation date: (9/15/00 11:25:35
     *  AM)
     *
     *@param  XSize  int
     *@param  Theta  double[]
     *@param  Phi    double
     */
    private OrdinalSpatialDataAgumentation TDataAgumentation;
    //private OrdinalSpatialQL[] TpartialAvgSingle;
    //private OrdinalSpatialSingleAveregeQL SingleAveregeQL;
    private final static int[][] TimeSeriesNeighborgood = {
            {
            -1
            }, {
            0
            }
            };
    private final static int[][] HNeighborgood = {
            {
            -1, 1
            }, {
            0, 0
            }
            };
    private final static int[][] VNeighborgood = {
            {
            0, 0
            }, {
            -1, 1
            }
            };


    /**
     *  Empty constructor
    public OrdinalSpatial() {
        super();
    }


    /**
     *  Constructor for the OrdinalSpatial object
     *
     *@param  t  Description of Parameter
     */
    public OrdinalSpatial(OrdinalSpatial t) {
        super(t);
    }


    /**
     *  Insert the method's description here.
     *
     *@param  y  int[]
     */
    public OrdinalSpatial(int[] y) {
        super(y);
    }


    /**
     *  Insert the method's description here.
     *
     *@param  y  int[]
     *@param  z  double[][]
     */
    public OrdinalSpatial(int[] y, double[][] z) {
        super(y, z);
    }


    /**
     *  Insert the method's description here.
     *
     *@param  y             Description of Parameter
     *@param  xLocation     Description of Parameter
     *@param  yLocation     Description of Parameter
     *@param  z             Description of Parameter
     *@param  Neighborgood  Description of Parameter
     */
    public OrdinalSpatial(int[] y, int[] xLocation, int[] yLocation, double[][] z,
            int[][] Neighborgood) {
        super(y, xLocation, yLocation, z, Neighborgood);
    }


    /**
     *  put your documentation comment here
     *
     *@param  XSize  Description of Parameter
     *@param  Theta  Description of Parameter
     *@param  Phi    Description of Parameter
     */
    public OrdinalSpatial(int XSize, double[] Theta, double Phi) {
        this(XSize, Theta, Phi, null, null, TimeSeriesNeighborgood);
    }


    /**
     *  Insert the method's description here. Creation date: (9/15/00 11:25:35
     *  AM)
     *
     *@param  XSize         int
     *@param  Theta         double[]
     *@param  Phi           double
     *@param  Neighborgood  int[][]
     */
    public OrdinalSpatial(int XSize, double[] Theta, double Phi, int[][] Neighborgood) {
        this(XSize, Theta, Phi, null, null, Neighborgood[0]);
    }

    /**
     *  Insert the method's description here. Creation date: (9/15/00 11:25:10
     *  AM)
     *
     *@param  XSize  int
     *@param  Theta  double[]
     *@param  Phi    double
     *@param  Beta   double[]
     *@param  z      double[][][]
     */
    public OrdinalSpatial(int XSize, double[] Theta, double Phi, double[] Beta,
            double[][][] z) {
        this(XSize, 1, Theta, Phi, Beta, z, TimeSeriesNeighborgood);
    }


    /**
     *  Insert the method's description here. Creation date: (9/15/00 11:25:10
     *  AM)
     *
     *@param  XSize         int
     *@param  Theta         double[]
     *@param  Phi           double
     *@param  Beta          double[]
     *@param  z             double[][][]
     *@param  Neighborgood  int[][]
     */
    public OrdinalSpatial(int XSize, double[] Theta, double Phi, double[] Beta,
            double[][][] z, int[][] Neighborgood) {
        this(XSize, 1, Theta, Phi, Beta, z, Neighborgood);
    }


    /**
     *  Insert the method's description here. Creation date: (9/15/00 11:25:10
     *  AM)
     *
     *@param  XSize         int
     *@param  Theta         double[]
     *@param  Phi           double
     *@param  Beta          double[]
     *@param  z             double[][][]
     *@param  Neighborgood  int[]
     */
    public OrdinalSpatial(int XSize, double[] Theta, double Phi, double[] Beta,
            double[][][] z, int[] Neighborgood) {
        this(XSize, 1, Theta, Phi, Beta, z, neighborgood1to2(Neighborgood));
    }


    /**
     *  Insert the method's description here. Creation date: (9/15/00 11:25:35
     *  AM)
     *
     *@param  XSize         int
     *@param  Theta         double[]
     *@param  Phi           double
     *@param  Neighborgood  int[]
     */
    public OrdinalSpatial(int XSize, double[] Theta, double Phi, int[] Neighborgood) {
        this(XSize, Theta, Phi, null, null, Neighborgood);
    }


    /**
     *  Insert the method's description here. Creation date: (9/15/00 12:32:21
     *  PM)
     *
     *@param  Neighborgood  int[]
     *@param  XSize         Description of Parameter
     *@param  YSize         Description of Parameter
     *@param  Theta         Description of Parameter
     *@param  Phi           Description of Parameter
     */
    /**
     *  Insert the method's description here. Creation date: (9/15/00 12:32:21
     *  PM)
     *
     *  Insert the method's description here. Creation date: (9/15/00 12:32:21
     *  PM) Insert the method's description here. Creation date: (9/15/00
     *  12:32:21 PM) Insert the method's description here. Creation date:
     *  (9/15/00 12:32:21 PM) Insert the method's description here. Creation
     *  date: (9/15/00 12:32:21 PM) Insert the method's description here.
     *  Creation date: (9/15/00 12:32:21 PM) Insert the method's description
     *  here. Creation date: (9/15/00 12:32:21 PM) Insert the method's
     *  description here. Creation date: (9/15/00 12:32:21 PM) Insert the
     *  method's description here. Creation date: (9/15/00 12:32:21 PM) Insert
     *  the method's description here. Creation date: (9/15/00 12:32:21 PM)
     *  Insert the method's description here. Creation date: (9/15/00 12:32:21
     *  PM) Insert the method's description here. Creation date: (9/15/00
     *  12:32:21 PM) Insert the method's description here. Creation date:
     *  (9/15/00 12:32:21 PM) Insert the method's description here. Creation
     *  date: (9/15/00 12:32:21 PM) Insert the method's description here.
     *  Creation date: (9/15/00 12:32:21 PM) Insert the method's description
     *  here. Creation date: (9/15/00 12:32:21 PM) Insert the method's
     *  description here. Creation date: (9/15/00 12:32:21 PM) Insert the
     *  method's description here. Creation date: (9/15/00 12:32:21 PM) Insert
     *  the method's description here. Creation date: (9/15/00 12:32:21 PM)
     *  Insert the method's description here. Creation date: (9/15/00 12:32:21
     *  PM) Insert the method's description here.
     *
     *@param  Neighborgood  int[]
     *@param  XSize         Description of Parameter
     *@param  YSize         Description of Parameter
     *@param  Theta         Description of Parameter
     *@param  Phi           Description of Parameter
     */
    public OrdinalSpatial(int XSize, int YSize, double[] Theta, double Phi,
            int[][] Neighborgood) {
        this(XSize, YSize, Theta, Phi, null, null, Neighborgood);
    }


    /**
     *  Insert the method's description here.
     *
     *@param  YSize         int
     *@param  XSize         int
     *@param  Theta         double[]
     *@param  Phi           double
     *@param  Beta          double[] is a true value for covariates coffitints.
     *@param  z             double[][][] is a covariates coeffitient
     *      [Xloc][Yloc][cavariate]
     *@param  Neighborgood  Description of Parameter
     */
    public OrdinalSpatial(int XSize, int YSize, double[] Theta, double Phi,
            double[] Beta, double[][][] z, int[][] Neighborgood) {
        if (XSize < 2 || YSize < 2) {
            TimeSeries = true;
            if (XSize < 2) {
                if (YSize < 2) {
                    new IllegalAccessException("Xsize or YSize has be biger than 1");
                }
                XSize = YSize;
            }
            YSize = 1;
        }
        N = XSize * YSize;
        if (Beta != null && Beta.length > 0) {
            if (Beta.length != z[0][0].length) {
                throw new IllegalArgumentException("z and beta must agree in dimentions");
            }
            if (z.length != XSize || z[0].length != YSize) {
                throw new IllegalArgumentException("z must agree in dimentions with XSize and YSIze");
            }
            K = Beta.length;
            TrueBeta = (double[]) Beta.clone();
            beta = new double[K];
            Z = new double[N][K];
        }
        if (TimeSeries) {
            NeighborDefinition = new int[2][Neighborgood[0].length];
            NeighborDefinition[0] = (int[]) Neighborgood[0].clone();
            for (int i = 0; i < Neighborgood[0].length; i++) {
                NeighborDefinition[1][i] = 0;
            }
        }
        else {
            NeighborDefinition = Neighborgood;
        }
        TruePhi = Phi;
        TrueTheta = (double[]) Theta.clone();
        xLoc = new int[N];
        yLoc = new int[N];
        for (int j = 0; j < YSize; j++) {
            for (int i = 0; i < XSize; i++) {
                int jj = j * YSize + i;
                xLoc[jj] = i;
                yLoc[jj] = j;
                for (int k = 0; k < K; k++) {
                    Z[jj][k] = z[i][j][k];
                }
            }
        }
        createLoc();
        L = Theta.length;
        if (!TimeSeries) {
            createSigma(Phi);
        }
        TrueX = new double[N];
        generateY();
        theta = (double[]) Theta.clone();
        phi = Phi;
    }


    /**
     *  Insert the method's description here.
     */
    /**
     *  Insert the method's description here.
     *
     *  Insert the method's description here. Insert the method's description
     *  here. Insert the method's description here. Insert the method's
     *  description here. Insert the method's description here. Insert the
     *  method's description here. Insert the method's description here. Insert
     *  the method's description here. Insert the method's description here.
     *  Insert the method's description here. Insert the method's description
     *  here. Insert the method's description here. Insert the method's
     *  description here. Insert the method's description here. Insert the
     *  method's description here. Insert the method's description here. Insert
     *  the method's description here. Insert the method's description here.
     *  Insert the method's description here. Insert the method's description
     *  here.
     *
     */
    public void generateY() {
        //TrueX = Rnd.mvnor_sqrtSigma(sqrtSigma);
        generateX();
        if (K > 0) {
            double[] xt = La.times(Z, TrueBeta);
            for (int i = 0; i < K; i++) {
                for (int j = 0; j < N; j++) {
                    TrueX[j] += xt[i];
                }
            }
        }
        Y = new int[N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < TrueTheta.length; j++) {
                if (TrueX[i] > TrueTheta[j]) {
                    Y[i] = j + 1;
                }
            }
        }
        theta = (double[]) TrueTheta;
        phi = TruePhi;
        return;
    }


    /**
     *  Insert the method's description here.
     *
     *@return    double
     */
    public double lA() {
        return lA(Type);
    }


    /**
     *  Insert the method's description here.
     *
     *@param  phi    double
     *@param  theta  double[]
     *@return        double
     */
    public double lA(double phi, double[] theta) {
        return lA(Type, phi, theta);
    }


    /**
     *  Insert the method's description here.
     *
     *@param  type  Description of Parameter
     *@return       double
     */
    public double lA(String type) {
        type = type.toLowerCase();
        if (type.equals("mean")) {
            return lMnA();
        }
        if (type.equals("median")) {
            return lMdA();
        }
        return lA(Type);
    }


    /**
     *  Insert the method's description here.
     *
     *@param  phi    double
     *@param  theta  double[]
     *@param  type   Description of Parameter
     *@return        double
     */
    public double lA(String type, double phi, double[] theta) {
        type = type.toLowerCase();
        if (type.equals("mean")) {
            return lMnA(phi, theta);
        }
        if (type.equals("median")) {
            return lMdA(phi, theta);
        }
        return lA(Type, phi, theta);
    }


    /**
     *  Insert the method's description here.
     *
     *@return    double
     */
    public double lCLTrueX() {
        return lCLTrueX(phi, beta);
    }


    /**
     *  Description of the Method
     *
     *@param  Phi  Description of Parameter
     *@return      Description of the Returned Value
     */
    public double lCLTrueX(double Phi) {
        return lCLTrueX(Phi, beta);
    }


    /**
     *  Description of the Method
     *
     *@param  Phi   Description of Parameter
     *@param  Beta  Description of Parameter
     *@return       Description of the Returned Value
     */
    public double lCLTrueX(double Phi, double[] Beta) {
        T = new OrdinalSpatialCLTrueX(TrueX, Phi, Loc, Beta, Z, NeighborDefinition);
        return T.lL();
    }


    /**
     *  Description of the Method
     *
     *@return    Description of the Returned Value
     */
    public double lCLTrueXPLE() {
        return lCLTrueX(phi, beta);
    }


    /**
     *  Description of the Method
     *
     *@param  Phi  Description of Parameter
     *@return      Description of the Returned Value
     */
    public double lCLTrueXPLE(double Phi) {
        return lCLTrueX(Phi, beta);
    }


    /**
     *  Description of the Method
     *
     *@param  Phi   Description of Parameter
     *@param  Beta  Description of Parameter
     *@return       Description of the Returned Value
     */
    public double lCLTrueXPLE(double Phi, double[] Beta) {
        T = new OrdinalSpatialCLTrueXPLE(TrueX, Phi, Loc, Beta, Z, NeighborDefinition);
        return T.lL();
    }


    /**
     *  Insert the method's description here.
     *
     *@return    Description of the Returned Value
     */
    public double lL() {
        return lL(phi, theta, beta);
    }


    /**
     *  Insert the method's description here.
     *
     *@param  phi    double
     *@param  theta  Description of Parameter
     *@param  beta   Description of Parameter
     *@return        double
     */
    public double lL(double phi, double[] theta, double[] beta) {
        return 0;
    }


    /**
     *  Insert the method's description here.
     *
     *@return    double
     */
    public double lQL() {
        return lQL(phi, theta, beta);
    }


    /**
     *  Insert the method's description here.
     *
     *@param  phi    double
     *@param  theta  double[]
     *@return        double
     */
    public double lQL(double phi, double[] theta) {
        return lQL(phi, theta, beta);
    }


    /**
     *  Insert the method's description here.
     *
     *@param  phi    double
     *@param  theta  double[]
     *@param  beta   double[]
     *@return        double
     */
    public double lQL(double phi, double[] theta, double[] beta) {
        T = new OrdinalSpatialQL(Y, phi, theta, Loc, beta, Z, NeighborDefinition);
        return T.lL();
    }


    /**
     *  Insert the method's description here.
     */
    public void mA() {
        mA(Type);
    }


    /**
     *  Insert the method's description here.
     *
     *@param  type  Description of Parameter
     */
    public void mA(String type) {
        type = type.toLowerCase();
        if (type.equals("mean")) {
            mMnA();
        }
        else if (type.equals("median")) {
            mMdA();
        }
        else {
            mA(Type);
        }
    }


    /**
     *  Insert the method's description here.
     */
    public void mCLETrueX() {
        if (TTrueX == null) {
            TTrueX = new OrdinalSpatialCLTrueX(TrueX, phi, Loc, beta, Z, NeighborDefinition);
        }
        else {
            TTrueX.updateX(this.TrueX);
        }
        TTrueX.mLE();
        this.phi = TTrueX.phi;
        if (K > 0) {
            this.beta = (double[]) TTrueX.beta.clone();
        }
        itotal = TTrueX.itotal;
    }


    /**
     *  Description of the Method
     *
     *@param  Phi  Description of Parameter
     */
    public void mCLETrueX(double Phi) {
        if (TTrueX == null) {
            TTrueX = new OrdinalSpatialCLTrueX(TrueX, phi, Loc, beta, Z, NeighborDefinition);
        }
        else {
            TTrueX.updateX(this.TrueX);
        }
        TTrueX.mLE(Phi);
        this.phi = TTrueX.phi;
        if (K > 0) {
            this.beta = (double[]) TTrueX.beta.clone();
        }
        itotal = TTrueX.itotal;
    }


    /**
     *  Description of the Method
     *
     *@param  Beta  Description of Parameter
     */
    public void mCLETrueX(double[] Beta) {
        if (TTrueX == null) {
            TTrueX = new OrdinalSpatialCLTrueX(TrueX, phi, Loc, beta, Z, NeighborDefinition);
        }
        else {
            TTrueX.updateX(this.TrueX);
        }
        TTrueX.mLE(null, Beta);
        this.phi = TTrueX.phi;
        if (K > 0) {
            this.beta = (double[]) TTrueX.beta.clone();
        }
        itotal = TTrueX.itotal;
    }


    /**
     *  Description of the Method
     */
    public void mCLETrueXPLE() {
        if (TTrueXPLE == null) {
            TTrueXPLE = new OrdinalSpatialCLTrueXPLE(TrueX, phi, Loc, beta, Z, NeighborDefinition);
        }
        else {
            TTrueXPLE.TrueX = this.TrueX;
        }
        TTrueXPLE.mLE();
        this.phi = TTrueXPLE.phi;
        if (K > 0) {
            this.beta = (double[]) TTrueXPLE.beta.clone();
        }
        itotal = TTrueXPLE.itotal;
    }


    /**
     *  Description of the Method
     *
     *@param  Phi  Description of Parameter
     */
    public void mCLETrueXPLE(double Phi) {
        if (TTrueXPLE == null) {
            TTrueXPLE = new OrdinalSpatialCLTrueXPLE(TrueX, phi, Loc, beta, Z, NeighborDefinition);
        }
        else {
            TTrueXPLE.TrueX = this.TrueX;
        }
        TTrueXPLE.mLE(Phi);
        this.phi = TTrueXPLE.phi;
        if (K > 0) {
            this.beta = (double[]) TTrueXPLE.beta.clone();
        }
        itotal = TTrueXPLE.itotal;
    }


    /**
     *  Description of the Method
     *
     *@param  Beta  Description of Parameter
     */
    public void mCLETrueXPLE(double[] Beta) {
        if (TTrueXPLE == null) {
            TTrueXPLE = new OrdinalSpatialCLTrueXPLE(TrueX, phi, Loc, beta, Z, NeighborDefinition);
        }
        else {
            TTrueXPLE.TrueX = this.TrueX;
        }
        TTrueXPLE.mLE(null, Beta);
        this.phi = TTrueXPLE.phi;
        if (K > 0) {
            this.beta = (double[]) TTrueXPLE.beta.clone();
        }
        itotal = TTrueXPLE.itotal;
    }


    /**
     *  put your documentation comment here
     */
    public void mExA() {
        T = new OrdinalSpatialExA(Y, phi, theta, Loc, beta, Z, NeighborDefinition);
        T.mLE();
        this.phi = T.phi;
        if (K > 0) {
            this.beta = (double[]) T.beta.clone();
        }
        if (L > 0) {
            this.theta = (double[]) T.theta.clone();
        }
        itotal = T.itotal;
    }


    /**
     *  Insert the method's description here.
     */
    public void mQLE() {
        T = new OrdinalSpatialQL(Y, phi, theta, Loc, beta, Z, NeighborDefinition);
        T.mLE();
        this.phi = T.phi;
        if (K > 0) {
            this.beta = (double[]) T.beta.clone();
        }
        if (L > 0) {
            this.theta = (double[]) T.theta.clone();
        }
        itotal = T.itotal;
    }


    /**
     *  Insert the method's description here.
     */
    public void mQLEM() {
        T = new OrdinalSpatialQLEM(Y, phi, theta, Loc, beta, Z, NeighborDefinition);
        T.mLE();
        this.phi = T.phi;
        if (K > 0) {
            this.beta = (double[]) T.beta.clone();
        }
        if (L > 0) {
            this.theta = (double[]) T.theta.clone();
        }
        itotal = T.itotal;
    }


    /**
     *  put your documentation comment here
     *
     *@param  phi    Description of Parameter
     *@param  theta  Description of Parameter
     *@param  beta   Description of Parameter
     *@return        Description of the Returned Value
     */
//	public void mSingleAveregeQL() {
//		T = new OrdinalSpatialSingleAveregeQL(Y, phi, theta, Loc, beta, Z,
//				NeighborDefinition);
//		T.mLE();
//		this.phi = T.phi;
//		if (K > 0) {
//			this.beta = (double[]) T.beta.clone();
//		}
//		if (L > 0) {
//			this.theta = (double[]) T.theta.clone();
//		}
//	}


    /**
     *  put your documentation comment here
     *
     *  put your documentation comment here put your documentation comment here
     *  put your documentation comment here put your documentation comment here
     *  put your documentation comment here put your documentation comment here
     *
     *@param  phi    Description of Parameter
     *@param  theta  Description of Parameter
     *@param  beta   Description of Parameter
     *@return        Description of the Returned Value
     *@return        Description of the Returned Value
     *@return        Description of the Returned Value
     *@return        Description of the Returned Value
     *@return        Description of the Returned Value
     *@return        Description of the Returned Value
     *@return
     */
    public double lA(double phi, double[] theta, double[] beta) {
        return lA(Type, phi, theta, beta);
    }


    /**
     *  Insert the method's description here.
     *
     *@param  phi    double
     *@param  theta  double[]
     *@param  beta   double[]
     *@param  type   Description of Parameter
     *@return        double
     */
    public double lA(String type, double phi, double[] theta, double[] beta) {
        type = type.toLowerCase();
        if (type.equals("mean")) {
            return lMnA(phi, theta, beta);
        }
        if (type.equals("median")) {
            return lMdA(phi, theta, beta);
        }
        return lA(Type, phi, theta, beta);
    }


    /**
     *  Insert the method's description here.
     */
    public void mQLENearestAverege() {
        if (TNearestAvg == null) {
            //			int[][] neighborhood = new int[2][2];
            //			TNearestAvg = new OrdinalSpatialQL[(MyMath.factorial(NeighborDefinition[0].length) / (2 * MyMath.factorial(
            //					NeighborDefinition[0].length - 2)))];
            //			int ij = 0;
            //			for (int i = 0; i < NeighborDefinition[0].length; i++) {
            //				for (int j = i + 1; j < NeighborDefinition[0].length; j++) {
            //					neighborhood[0][0] = NeighborDefinition[0][i];
            //					neighborhood[1][0] = NeighborDefinition[1][i];
            //					neighborhood[0][1] = NeighborDefinition[0][j];
            //					neighborhood[1][1] = NeighborDefinition[1][j];
            //					TNearestAvg[ij++] = new OrdinalSpatialQL(Y, xLoc, yLoc,
            //							Z, neighborhood);
            //				}
            //			}
            TNearestAvg = new OrdinalSpatialQL[2];
            int[][] neighborhood1 = {
                    {
                    -1, 1
                    }, {
                    0, 0
                    }
                    };
            TNearestAvg[0] = new OrdinalSpatialQL(Y, xLoc, yLoc, Z, neighborhood1);
            int[][] neighborhood2 = {
                    {
                    0, 0
                    }, {
                    -1, 1
                    }
                    };
            TNearestAvg[1] = new OrdinalSpatialQL(Y, xLoc, yLoc, Z, neighborhood2);
        }
        else {
            for (int i = 0; i < TNearestAvg.length; i++) {
                TNearestAvg[i].Y = Y;
            }
        }
        for (int i = 0; i < TNearestAvg.length; i++) {
            TNearestAvg[i].phi = phi;
            if (L > 0) {
                TNearestAvg[i].theta = (double[]) theta.clone();
            }
            if (K > 0) {
                TNearestAvg[i].beta = (double[]) beta.clone();
            }
            long t = System.currentTimeMillis();
            TNearestAvg[i].mLE();
            System.out.print(TNearestAvg[i].itotal + " " + (System.currentTimeMillis()
                     - t));
        }
        this.phi = TNearestAvg[0].phi;
        this.itotal = TNearestAvg[0].itotal;
        if (K > 0) {
            this.beta = (double[]) TNearestAvg[0].beta.clone();
        }
        if (L > 0) {
            this.theta = (double[]) TNearestAvg[0].theta.clone();
        }
        for (int i = 1; i < TNearestAvg.length; i++) {
            this.phi += TNearestAvg[i].phi;
            if (K > 0) {
                this.beta = La.plus(this.beta, TNearestAvg[i].beta);
            }
            if (L > 0) {
                this.theta = La.plus(this.theta, TNearestAvg[i].theta);
            }
        }
        this.phi *= NeighborDefinition[0].length / (double) (2 * TNearestAvg.length);
        if (K > 0) {
            this.beta = La.times(this.beta, 1 / (double) TNearestAvg.length);
        }
        if (L > 0) {
            this.theta = La.times(this.theta, 1 / (double) TNearestAvg.length);
        }
        //theta=(double[]) T.theta.clone();
    }


    /**
     *  Insert the method's description here.
     */
    public void dataAgumentation() {
        if (TDataAgumentation == null) {
            TDataAgumentation = new OrdinalSpatialDataAgumentation(Y, phi,
                    theta, Loc, beta, Z, NeighborDefinition);
        }
        else {
            TDataAgumentation.Y = Y;
        }
        TDataAgumentation.dataAugmentation();
        this.phi = TDataAgumentation.phi;
        //if (K > 0)
        //this.beta = (double[]) T.beta.clone();
        this.theta = (double[]) TDataAgumentation.theta;
    }


    /**
     *  Description of the Method
     *
     *@param  inteructive  Description of Parameter
     */
    public void dataAgumentation(boolean inteructive) {
        if (TDataAgumentation == null) {
            TDataAgumentation = new OrdinalSpatialDataAgumentation(Y, phi,
                    theta, Loc, beta, Z, NeighborDefinition);
        }
        else {
            TDataAgumentation.Y = Y;
        }
        TDataAgumentation.dataAugmentation(inteructive);
        this.phi = TDataAgumentation.phi;
        //if (K > 0)
        //this.beta = (double[]) T.beta.clone();
        this.theta = (double[]) TDataAgumentation.theta;
    }


    /**
     *  put your documentation comment here
     */
//	public void mQLESingleAveregee() {
//		if (TpartialAvgSingle == null) {
//			TpartialAvgSingle = new OrdinalSpatialQL[NeighborDefinition[0].length];
//			int[][] neighborgood = new int[2][1];
//			for (int i = 0; i < NeighborDefinition[0].length; i++) {
//				neighborgood[0][0] = NeighborDefinition[0][i];
//				neighborgood[1][0] = NeighborDefinition[1][i];
//				TpartialAvgSingle[i] = new OrdinalSpatialQL(Y, xLoc, yLoc,
//						Z, neighborgood);
//			}
//		}
//		else {
//			for (int i = 0; i < TpartialAvgSingle.length; i++) {
//				TpartialAvgSingle[i].updateY(Y);
//			}
//		}
//		for (int i = 0; i < TpartialAvgSingle.length; i++) {
//			TpartialAvgSingle[i].phi = phi;
//			TpartialAvgSingle[i].theta = (double[]) theta.clone();
//			TpartialAvgSingle[i].mLE();
//		}
//		//for (int i = 0; i < TpartialAvgSingle.length; i++)
//		// System.out.print(TpartialAvgSingle[i].phi + " ");
//		this.phi = TpartialAvgSingle[0].phi;
//		if (K > 0) {
//			this.beta = (double[]) TpartialAvgSingle[0].beta.clone();
//		}
//		if (L > 0) {
//			this.theta = (double[]) TpartialAvgSingle[0].theta.clone();
//		}
//		for (int i = 1; i < TpartialAvgSingle.length; i++) {
//			this.phi += TpartialAvgSingle[i].phi;
//			if (K > 0) {
//				this.beta = La.plus(this.beta, TpartialAvgSingle[i].beta);
//			}
//			if (L > 0) {
//				this.theta = La.plus(this.theta, TpartialAvgSingle[i].theta);
//			}
//		}
//		this.phi /= (double) (TpartialAvgSingle.length);
//		if (K > 0) {
//			this.beta = La.times(this.beta, 1 / (double) TpartialAvgSingle.length);
//		}
//		if (L > 0) {
//			this.theta = La.times(this.theta, 1 / (double) TpartialAvgSingle.length);
//		}
//		//theta=(double[]) T.theta.clone();
//	}


    /**
     *  put your documentation comment here
     *
     *  put your documentation comment here put your documentation comment here
     *  put your documentation comment here put your documentation comment here
     *  put your documentation comment here Insert the method's description
     *  here.
     *
     */
    public void mMdA() {
        if (TMdA == null) {
            TMdA = new OrdinalSpatialMdA(Y, phi, theta, Loc, beta, Z, NeighborDefinition);
        }
        else {
            TMdA.Y = Y;
        }
        TMdA.phi = phi;
        TMdA.beta = beta;
        TMdA.theta = theta;
        TMdA.itotal = 0;
        TMdA.mLE();
        this.phi = TMdA.phi;
        itotal = TMdA.itotal;
        //theta=(double[]) T.theta.clone();
    }


    /**
     *  Insert the method's description here.
     */
    public void mMnA() {
        if (TMnA == null) {
            TMnA = new OrdinalSpatialMnA(Y, phi, theta, Loc, beta, Z, NeighborDefinition);
        }
        else {
            TMnA.Y = Y;
        }
        TMnA.phi = phi;
        TMnA.beta = beta;
        TMnA.theta = theta;
        TMnA.itotal = 0;
        TMnA.mLE();
        this.phi = TMnA.phi;
        itotal = TMnA.itotal;
    }


    /**
     *  Description of the Method
     *
     *@param  Phi  Description of Parameter
     */
    public void mMnA(double Phi) {
        mMnA(Phi, null, null);
    }


    /**
     *  Description of the Method
     *
     *@param  Phi    Description of Parameter
     *@param  Theta  Description of Parameter
     */
    public void mMnA(double Phi, double[] Theta) {
        mMnA(Phi, Theta, null);
    }


    /**
     *  Description of the Method
     *
     *@param  Phi    Description of Parameter
     *@param  Theta  Description of Parameter
     *@param  Beta   Description of Parameter
     */
    public void mMnA(double Phi, double[] Theta, double[] Beta) {
        if (TMnA == null) {
            TMnA = new OrdinalSpatialMnA(Y, phi, theta, Loc, beta, Z, NeighborDefinition);
        }
        else {
            TMnA.Y = Y;
        }
        //TMnA.phi = phi;
        //TMnA.beta = beta;
        //TMnA.theta = theta;
        TMnA.itotal = 0;
        TMnA.mLE(Phi, Theta, Beta);
        this.phi = TMnA.phi;
        this.theta = (double[]) TMnA.theta;
        this.beta = (double[]) TMnA.beta;
        itotal = TMnA.itotal;
    }


    /**
     *  Description of the Method
     *
     *@param  Theta  Description of Parameter
     *@param  Beta   Description of Parameter
     */
    public void mMnA(double[] Theta, double[] Beta) {
        if (TMnA == null) {
            TMnA = new OrdinalSpatialMnA(Y, phi, theta, Loc, beta, Z, NeighborDefinition);
        }
        else {
            TMnA.Y = Y;
        }
        //TMnA.phi = phi;
        //TMnA.beta = beta;
        //TMnA.theta = theta;
        TMnA.itotal = 0;
        TMnA.mLE(Theta, Beta);
        this.theta = (double[]) TMnA.theta;
        this.beta = (double[]) TMnA.beta;
        this.phi = TMnA.phi;
        itotal = TMnA.itotal;
    }


    /**
     *  Description of the Method
     *
     *@param  Phi  Description of Parameter
     */
    public void mMdA(double Phi) {
        mMdA(Phi, null, null);
    }


    /**
     *  Description of the Method
     *
     *@param  Phi    Description of Parameter
     *@param  Theta  Description of Parameter
     */
    public void mMdA(double Phi, double[] Theta) {
        mMdA(Phi, Theta, null);
    }


    /**
     *  Description of the Method
     *
     *@param  Phi    Description of Parameter
     *@param  Theta  Description of Parameter
     *@param  Beta   Description of Parameter
     */
    public void mMdA(double Phi, double[] Theta, double[] Beta) {
        if (TMdA == null) {
            TMdA = new OrdinalSpatialMdA(Y, phi, theta, Loc, beta, Z, NeighborDefinition);
        }
        else {
            TMdA.Y = Y;
        }
        //TMnA.phi = phi;
        //TMnA.beta = beta;
        //TMnA.theta = theta;
        TMdA.itotal = 0;
        TMdA.mLE(Phi, Theta, Beta);
        this.phi = TMdA.phi;
        this.theta = (double[]) TMdA.theta;
        this.beta = (double[]) TMdA.beta;
        itotal = TMdA.itotal;
    }


    /**
     *  Description of the Method
     *
     *@param  Theta  Description of Parameter
     *@param  Beta   Description of Parameter
     */
    public void mMdA(double[] Theta, double[] Beta) {
        if (TMdA == null) {
            TMdA = new OrdinalSpatialMdA(Y, phi, theta, Loc, beta, Z, NeighborDefinition);
        }
        else {
            TMdA.Y = Y;
        }
        //TMnA.phi = phi;
        //TMnA.beta = beta;
        //TMnA.theta = theta;
        TMdA.itotal = 0;
        TMdA.mLE(Theta, Beta);
        this.theta = (double[]) TMdA.theta;
        this.beta = (double[]) TMdA.beta;
        itotal = TMdA.itotal;
    }


    /**
     *  Insert the method's description here. Creation date: (6/13/00 2:22:00
     *  PM)
     */
    private void generateX() {
        //double t = 0;
        if (phi == 0) {
            for (int i = 0; i < N; i++) {
                TrueX[i] = rnd.gaussian();
            }
        }
        else
                if (TimeSeries ||
                (NeighborDefinition[0].length == 1 &&
                ((NeighborDefinition[0][0]
                 == -1 && NeighborDefinition[1][0] == 0) ||
                (NeighborDefinition[0][0]
                 == 0 && NeighborDefinition[1][0] == -1)))) {
            TrueX = Rnd.nor(N, 0, Math.sqrt(1 - TruePhi * TruePhi));
            for (int i = 1; i < N; i++) {
                TrueX[i] += TruePhi * TrueX[i - 1];
            }
        }
        else {
            TrueX = Rnd.mvnor_sqrtSigma(sqrtSigma);
        }
    }


    /**
     *  Insert the method's description here.
     *
     *@return    double
     */
    private double lMdA() {
        return lMdA(phi, theta, beta);
    }


    /**
     *  Insert the method's description here.
     *
     *@param  phi    double
     *@param  theta  double[]
     *@return        double
     */
    private double lMdA(double phi, double[] theta) {
        return lMdA(phi, theta, beta);
    }


    /**
     *  Insert the method's description here.
     *
     *@param  phi    double
     *@param  theta  double[]
     *@param  beta   double[]
     *@return        double
     */
    private double lMdA(double phi, double[] theta, double[] beta) {
        T = new OrdinalSpatialMdA(Y, phi, theta, Loc, beta, Z, NeighborDefinition);
        return T.lL();
    }


    /**
     *  Insert the method's description here.
     *
     *@return    double
     */
    private double lMnA() {
        return lMdA(phi, theta, beta);
    }


    /**
     *  Insert the method's description here.
     *
     *@param  phi    double
     *@param  theta  double[]
     *@return        double
     */
    private double lMnA(double phi, double[] theta) {
        return lMnA(phi, theta, beta);
    }


    /**
     *  Insert the method's description here.
     *
     *@param  phi    double
     *@param  theta  double[]
     *@param  beta   double[]
     *@return        double
     */
    private double lMnA(double phi, double[] theta, double[] beta) {
        T = new OrdinalSpatialMnA(Y, phi, theta, Loc, beta, Z, NeighborDefinition);
        return T.lL();
    }


    /**
     *  put your documentation comment here
     *
     *@param  Neighborgood
     *@return
     */
    private static int[][] neighborgood1to2(int[] Neighborgood) {
        int[][] a = new int[2][Neighborgood.length];
        a[0] = (int[]) Neighborgood.clone();
        return a;
    }
}

