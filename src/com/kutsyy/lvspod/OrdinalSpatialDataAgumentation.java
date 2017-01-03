/*
 * put your module comment here
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 */

package com.kutsyy.lvspod;

import com.kutsyy.util.*;
import java.awt.*;
import javax.swing.*;

/**
 *  Title: com.kutsyy.lvspod Description: This library is set of rutins for
 *  Laten Variable model for Spatialy dependent Ordinal Data Copyright:
 *  Copyright (c) 2000 Company: The University of Michigan
 *
 *@author     Vadim Kutsyy
 *@created    January 9, 2001
 *@version    1.0
 */
public class OrdinalSpatialDataAgumentation extends OrdinalSpatialAbstract {
	/**
	 *  Description of the Field
	 */
	public double alpha = 2;
	/**
	 *  Description of the Field
	 */
	public double beta = 2;
	double[][] thetaSim = null;
	double[] phiSim = null;
	private int m = 1;
	private double[][] xSim = null;
	private double[] w = null;
	private OrdinalSpatialDataAgumentationPlot Plot;
	private boolean inteructive = false;


	/**
	 *  put your documentation comment here
	 */
	public OrdinalSpatialDataAgumentation() {
		super();
	}


	/**
	 *  OrdinalSpatialEM constructor comment.
	 *
	 *@param  y                   int[]
	 *@param  phi                 double
	 *@param  theta               double[]
	 *@param  loc                 com.kutsyy.lvspod.SpatialPoint[]
	 *@param  beta                cern.colt.matrix.DoubleMatrix1D
	 *@param  z                   cern.colt.matrix.DoubleMatrix2D
	 *@param  neighborDefinition  int[][]
	 */
	protected OrdinalSpatialDataAgumentation(int[] y, double phi, double[] theta,
			com.kutsyy.lvspod.SpatialPoint[] loc, double[] beta, double[][] z,
			int[][] neighborDefinition) {
		super(y, phi, theta, loc, beta, z, neighborDefinition);
	}


	/**
	 *  OrdinalSpatialPLEM constructor comment.
	 *
	 *@param  y          int[]
	 *@param  xLocation  int[]
	 *@param  yLocation  int[]
	 *@param  z          double[][]
	 */
	protected OrdinalSpatialDataAgumentation(int[] y, int[] xLocation, int[] yLocation,
			double[][] z) {
		super(y, xLocation, yLocation, z);
	}


	/**
	 *  put your documentation comment here
	 *
	 *@return
	 */
	public double lL() {

		throw new IllegalStateException("not implemented yet");
	}


	/**
	 *  Insert the method's description here. Creation date: (11/15/2000 2:23:37
	 *  AM)
	 */
	public final void dataAugmentation() {
		dataAugmentation(0, true);
	}


	/**
	 *  Description of the Method
	 *
	 *@param  Inteructive  Description of Parameter
	 */
	public final void dataAugmentation(boolean Inteructive) {
		dataAugmentation(0, Inteructive);
	}


	/**
	 *  Insert the method's description here. Creation date: (11/15/2000 2:17:51
	 *  AM)
	 *
	 *@param  m            int
	 *@param  Inteructive  Description of Parameter
	 */
	public final void dataAugmentation(int m, boolean Inteructive) {
		inteructive = Inteructive;
		if (m < 1) {
			m = N * N;
		}
		this.m = m;
		if (xSim == null || xSim.length != m || xSim[0] == null || xSim[0].length != N) {
			xSim = new double[m][N];
		}
		if (L == 0) {
			thetaSim = null;
		}
		else if (thetaSim == null || thetaSim.length != m || thetaSim[0] == null || thetaSim[0].length != L) {
			thetaSim = new double[m][L];
		}
		if (phiSim == null || phiSim.length != m) {
			phiSim = new double[m];
		}
		if (w == null || w.length != m) {
			w = new double[m];
		}
		java.util.Arrays.fill(phiSim, phi);
		int i;
		int j;
		for (i = 0; i < m; i++) {
			thetaSim[i] = (double[]) theta.clone();
		}
		double[] phiSimKeep;
		if (inteructive) {
			try {
				UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
			}
			catch (Exception e) {
				e.printStackTrace();
			}
			if (Plot == null) {
				Plot = new OrdinalSpatialDataAgumentationPlot();
				Plot.validate();
				Plot.setVisible(true);
			}
			Plot.initialize(phiSim, m);
		}
		loop :
		for (int k = 0; k < m; k++) {
			phiSimKeep = (double[]) phiSim.clone();
			//for (i = 0; i < m; i++) {
			generateX();
			//}
			for (i = 0; i < m; i++) {
				j = unifDiscrete(m - 1);
				phiSim[i] = simPhi(xSim[j], phiSimKeep[j]);
				simTheta(xSim[j], thetaSim[i]);
				if (inteructive) {
					Plot.jProgressBarPhi.setValue(i);
				}
			}
			if (inteructive) {
				Plot.update(phiSim);
				if (Plot.stopDataAugmentation) {
					break loop;
				}
			}
			//System.out.println(k+" "+Stat.mean(phiSim));
		}
		//while (Math.abs(Stat.mean(phiSim) - Stat.mean(phiSimKeep)) > epsilon);
		phi = Stat.mean(phiSim);
		thetaSim = La.t(thetaSim);
		for (i = 0; i < L; i++) {
			theta[i] = Stat.mean(thetaSim[i]);
		}
	}


	/**
	 *  Description of the Method
	 *
	 *@param  m  Description of Parameter
	 *@return    Description of the Returned Value
	 */
	public int unifDiscrete(int m) {
		return (int) Math.floor(rnd.raw() * (m + 1));
	}


	/**
	 *  Description of the Method
	 *
	 *@param  low  Description of Parameter
	 *@return      Description of the Returned Value
	 */
	public double ntrc(double low) {
		if (Double.isNaN(low)) {
			throw new IllegalArgumentException("low is NaN");
		}
//		double cut = .45;
//		if (low > cut) {
//			double z = -Math.log(rnd.raw()) / low;
//			while (rnd.raw() > Math.exp(-.5 * z * z)) {
//				z = -Math.log(rnd.raw()) / low;
//			}
//			if (z < 0) {
//				throw new IllegalStateException("z<0");
//			}
//			return z + low;
//		}
//		double x = rnd.gaussian();
//		while (x < low) {
//			x = rnd.gaussian();
//		}
//		return x;
		double LOW = Cdf.nor(low);
		return CdfInv.nor(rnd.raw() * (1 - LOW) + LOW);
	}



	/**
	 *  Insert the method's description here. Creation date: (2/10/00 2:31:07 PM)
	 *
	 *@param  low   double
	 *@param  high  double
	 *@return       double
	 */
	public double ntrc(double low, double high) {
		if (Double.isNaN(low) || Double.isNaN(high)) {
			throw new IllegalArgumentException("low or high is NaN");
		}
		if (low >= high) {
			throw new IllegalArgumentException(" low limit has be smaller than higher");
		}
		if (Double.NEGATIVE_INFINITY == low) {
			return -ntrc(-high);
		}
		if (Double.POSITIVE_INFINITY == high) {
			return ntrc(low);
		}
//		if (high - low > 0.2 &&
//				(!((high > 1.5 && low > 1.5) || (high < 1.5 && low < 1.5)))) {
//			double aa = 0;
//			if (high > 0 && low < 0) {
//				aa = Pdf.nor(0);
//			}
//			else {
//				aa = Pdf.nor(low);
//				aa = Math.max(aa, Pdf.nor(high));
//			}
//			do {
//				double y = rnd.raw()*(high-low)+low;
//				double u = rnd.raw();
//				if (u * aa <= Pdf.nor(y)) {
//					return y;
//				}
//			} while (true);
//		}
		double LOW = Cdf.nor(low);
		return CdfInv.nor(rnd.raw() * (Cdf.nor(high) - LOW) + LOW);
	}



	/**
	 *  Insert the method's description here. Creation date: (2/10/00 2:32:46 PM)
	 *
	 *@param  low   double
	 *@param  high  double
	 *@param  mu    double
	 *@param  sd    double
	 *@return       double
	 */
	public double ntrc(double low, double high, double mu, double sd) {
		low -= mu;
		high -= mu;
		low /= sd;
		high /= sd;
		return (ntrc(low, high) + mu) * sd;
	}



	/**
	 *  Insert the method's description here. Creation date: (2/10/00 2:36:04 PM)
	 *
	 *@param  n     int
	 *@param  low   double[]
	 *@param  high  double[]
	 *@return       double[]
	 */
	public double[] ntrc(int n, double low, double high) {
		double[] t = new double[n];
		for (int i = 0; i < n; i++) {
			t[i] = ntrc(low, high);
		}
		return t;
	}



	/**
	 *  Insert the method's description here. Creation date: (2/10/00 2:35:15 PM)
	 *
	 *@param  n     int
	 *@param  low   double[]
	 *@param  high  double[]
	 *@param  mu    double[]
	 *@param  sd    double[][]
	 *@return       double[]
	 */
	public double[] ntrc(
			int n,
			double low,
			double high,
			double mu,
			double sd) {
		double[] t = new double[n];
		low = -mu;
		high -= mu;
		low /= sd;
		high /= sd;
		for (int i = 0; i < n; i++) {
			t[i] = (ntrc(low, high) + mu) * sd;
		}
		return t;
	}


	/**
	 *  Insert the method's description here. Creation date: (11/15/2000 2:19:52
	 *  AM)
	 */
	protected void generateX() {
		if (thetaU == null) {
			thetaU = new double[N];
		}
		if (thetaL == null) {
			thetaL = new double[N];
		}
		if (x == null || x.length != N) {
			x = new double[N];
		}
		long tt = System.currentTimeMillis();
		for (int jj = 0; jj < m; jj++) {

			if (inteructive) {
				Plot.jProgressBarX.setValue(jj);
			}
			if (Y != null && thetaSim[jj] != null &&
					(thetaSim[jj].length == 1 || ((thetaSim[jj].length > 1)
					 && thetaSim[jj][1] > thetaSim[jj][0]))) {
				for (int i = 0; i < N; i++) {
					if (Y[i] == 0) {
						thetaU[i] = (thetaSim[jj][0]);
						thetaL[i] = Double.NEGATIVE_INFINITY;
					}
					else if (Y[i] == L) {
						thetaU[i] = Double.POSITIVE_INFINITY;
						thetaL[i] = thetaSim[jj][L - 1];
					}
					else {
						thetaL[i] = thetaSim[jj][Y[i] - 1];
						thetaU[i] = thetaSim[jj][Y[i]];
					}
				}
			}
			else {
				java.util.Arrays.fill(thetaU, Double.POSITIVE_INFINITY);
				java.util.Arrays.fill(thetaL, Double.NEGATIVE_INFINITY);
			}
			int mm = (int) Math.sqrt(N) + 1;
			int i;
			int j;
			for (int ii = 0; ii < mm; ii++) {
				for (i = 0; i < N; i++) {
					double t = 0;
					for (j = 0; j < Loc[i].neighbors.length; j++) {
						t += xSim[jj][Loc[i].neighbors[j]];
					}
					t *= phiSim[jj] / Loc[i].neighbors.length;
					xSim[jj][i] = ntrc(thetaL[i], thetaU[i], t, Math.sqrt(1 - phiSim[jj] * phiSim[jj] / (double) Loc[i].neighbors.length));
				}
			}
//			for (int i = 0; i < N; i++) {
//				x[i] = Exp.ntrc(thetaL[i], thetaU[i]);
//			}
//			for (int i = 0; i < N; i++) {
//				if (Loc[i].neighbors.length > 0) {
//					double t = 0;
//					for (int j = 0; j < Loc[i].neighbors.length; j++) {
//						t += x[Loc[i].neighbors[j]];
//					}
//					xSim[jj][i] = Rnd.ntrc(thetaL[i], thetaU[i],
//                                        t * phiSim[jj]/Loc[i].neighbors.length,
//                                        Math.sqrt(1 - phiSim[jj] * phiSim[jj] / Loc[i].neighbors.length));
//				}
//				else {
//					xSim[jj][i] = Rnd.ntrc(thetaL[i], thetaU[i]);
//				}
//			}
//                        w[i]=
		}
//		System.out.println(tt - System.currentTimeMillis());
	}


	/**
	 *  Insert the method's description here. Creation date: (11/20/2000 8:37:59
	 *  PM)
	 *
	 *@param  X    Description of Parameter
	 *@param  Phi  Description of Parameter
	 *@return      Description of the Returned Value
	 */
	private double simPhi(double[] X, double Phi) {
		double newPhi;
		if (x == null) {
			x = new double[N];
		}
		else {
			java.util.Arrays.fill(x, 0);
		}
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < Loc[j].neighbors.length; k++) {
				x[j] += X[Loc[j].neighbors[k]];
			}
		}
		double s1 = 0;
		double s2 = 0;
		double s3 = 0;
		for (int j = 0; j < N; j++) {
			s1 += X[j] * X[j];
			s2 += X[j] * x[j];
			s3 += x[j] * x[j];
		}
		int n = NeighborDefinition[0].length;
		double lu1 = Math.log(rnd.raw() * Math.exp(-1.0 / (2 - 2 * Phi * Phi / n) * (
				Phi * Phi * s1 - 2 * Phi * s2 + s3)));
		double lu2 = Math.log(rnd.uniform(1, Math.pow(1 - Phi * Phi, -N / 2.0)));
		double lu3 = Math.log(rnd.raw() * Math.pow(1 - Phi, beta - 1));
		double lu4 = Math.log(rnd.raw() * Math.pow(1 + Phi, alpha - 1));
		double t11 = -1;
		double t12 = 1;
		if (lu1 != Double.NEGATIVE_INFINITY) {
			t11 = Math.sqrt(n * n * s2 * s2 - 2 * n * n * s1 * lu1 - n * n * s1 * s3 + 4 * lu1 * lu1 * n
					 + 2 * lu1 * n * s3);
			t12 = (NeighborDefinition[0].length * s2 + t11) / (NeighborDefinition[0].length * s1
					 - 2 * lu1);
			t11 = (NeighborDefinition[0].length * s2 - t11) / (NeighborDefinition[0].length * s1
					 - 2 * lu1);
		}
		double t2 = Math.sqrt(1 - Math.exp(-2.0 * lu2 / N));
		double t3 = 1 - Math.exp(lu3 / (beta - 1));
		double t4 = Math.exp(lu4 / (alpha - 1)) - 1;
		if (t11 < t4) {
			t11 = t4;
		}
		if (t12 > t3) {
			t12 = t3;
		}
		if (-t2 > t11) {
			if (t2 < t12) {
				newPhi = rnd.uniform(0, -t2 - t11 + t12 - t2);
				if (newPhi < -t2 - t11) {
					newPhi += t11;
				}
				else {
					newPhi += t2 + t2 + t11;
				}
			}
			else {
				newPhi = rnd.uniform(t11, -t2);
			}
		}
		else {
			newPhi = rnd.uniform(t2, t12);
		}
		if (Double.isNaN(newPhi) || Math.abs(newPhi) > 1) {
			throw new IllegalStateException("theta is NAN");
		}
		return newPhi;
	}


	/**
	 *  Insert the method's description here. Creation date: (11/23/2000 10:44:05
	 *  PM)
	 *
	 *@param  X      double[]
	 *@param  Theta  double[]
	 */
	private void simTheta(double[] X, double[] Theta) {
		if (thetal == null) {
			thetal = new double[L];
			thetau = new double[L];
		}
		java.util.Arrays.fill(thetal, Double.NEGATIVE_INFINITY);
		java.util.Arrays.fill(thetau, Double.POSITIVE_INFINITY);
		for (int i = 0; i < N; i++) {
			if (Y[i] == 0) {
				thetal[0] = Math.max(thetal[0], X[i]);
			}
			else if (Y[i] == L) {
				thetau[L - 1] = Math.min(thetau[L - 1], X[i]);
			}
			else {
				thetal[Y[i]] = Math.max(thetal[Y[i]], X[i]);
				thetau[Y[i] - 1] = Math.min(thetau[Y[i] - 1], X[i]);
			}
		}
		for (int i = 0; i < L; i++) {
			Theta[i] = Exp.ntrc(thetal[i], thetau[i]);
			if (Double.isNaN(Theta[i]) || (i > 0 && Theta[i] <= Theta[i - 1])) {
				throw new IllegalStateException("theta is NAN");
			}
		}
	}


}

