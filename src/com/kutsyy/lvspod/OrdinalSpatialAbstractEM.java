/*
 * put your module comment here
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 */

package com.kutsyy.lvspod;

import com.kutsyy.util.*;

/**
 *  Insert the type's description here. Creation date: (10/21/2000 2:26:35 PM)
 *
 *@author     vadim
 *@created    December 2, 2000
 *@author:
 */
public abstract class OrdinalSpatialAbstractEM extends OrdinalSpatialAbstract {

	/**
	 *  Insert the method's description here.
	 *
	 *@param  y  int[]
	 */
	public OrdinalSpatialAbstractEM(int[] y) {
		super(y);
	}


	/**
	 *  Insert the method's description here.
	 *
	 *@param  y  int[]
	 *@param  z  double[][]
	 */
	public OrdinalSpatialAbstractEM(int[] y, double[][] z) {
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
	public OrdinalSpatialAbstractEM(int[] y, int[] xLocation, int[] yLocation,
			double[][] z, int[][] Neighborgood) {
		super(y, xLocation, yLocation, z, Neighborgood);
	}


	/**
	 *  OrdinalSpatialEM constructor comment.
	 */
	protected OrdinalSpatialAbstractEM() {
		super();
	}


	/**
	 *  OrdinalSpatialEM constructor comment.
	 *
	 *@param  y          int[]
	 *@param  xLocation  int[]
	 *@param  yLocation  int[]
	 *@param  z          double[][]
	 */
	protected OrdinalSpatialAbstractEM(int[] y, int[] xLocation, int[] yLocation,
			double[][] z) {
		super(y, xLocation, yLocation, z);
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
	protected OrdinalSpatialAbstractEM(int[] y, double phi, double[] theta,
			com.kutsyy.lvspod.SpatialPoint[] loc, double[] beta, double[][] z,
			int[][] neighborDefinition) {
		super(y, phi, theta, loc, beta, z, neighborDefinition);
	}


	/**
	 *  Insert the method's description here. Creation date: (10/21/2000 2:54:12
	 *  PM)
	 *
	 *@return    Description of the Returned Value
	 */
	public double lL() {
		throw new IllegalStateException("log lik estimator is not implemented");
	}


	/**
	 *  Insert the method's description here. Creation date: (10/21/2000 2:30:40
	 *  PM)
	 */
	public void mLE() {
		double diff = Double.POSITIVE_INFINITY;
		while (diff > EPS) {
			double[] thetaOld = (double[]) theta.clone();
			double phiOld = phi;
			double[] betaOld = null;
			if (beta != null) {
				betaOld = (double[]) beta.clone();
			}
			e();
			m();
			diff = Math.abs(phi + phiOld);
			for (int i = 0; i < theta.length; i++) {
				diff += Math.abs(thetaOld[i] - theta[i]);
			}
			for (int i = 0; i < K; i++) {
				diff += Math.abs(betaOld[i] - beta[i]);
			}
		}
	}


	/**
	 *  Insert the method's description here. Creation date: (10/21/2000 2:28:40
	 *  PM)
	 */
	protected abstract void e();


	/**
	 *  Insert the method's description here.
	 *
	 *@param  phi    double
	 *@param  theta  double[]
	 *@param  beta   double[]
	 *@return        double
	 */
	protected double lL(double phi, double[] theta, cern.colt.matrix.DoubleMatrix1D beta) {
		throw new IllegalStateException("log lik estimator is not implemented");
	}


	/**
	 *  Insert the method's description here. Creation date: (10/21/2000 2:28:54
	 *  PM)
	 */
	protected void m() {
		phiGivenX();
		thetaGivenX();
	}


	/**
	 *  Insert the method's description here. Creation date: (4/2/00 5:17:05 PM)
	 */
	protected void phiGivenX() {
		int n = x.length;
		double Sxy = 0;
		double Sxx = 0;
		for (int i = 1; i < n; i++) {
			Sxy += x[i] * x[i - 1];
			Sxx += x[i - 1] * x[i - 1];
		}
		phi = Math.min(Math.max(Sxy / Sxx, -PHI_LIMIT), PHI_LIMIT);
	}


	/**
	 *  Insert the method's description here. Creation date: (4/2/00 5:16:52 PM)
	 */
	protected void thetaGivenX() {
		int l = theta.length;
		int n = x.length;
		double[] thetaL = new double[l];
		double[] thetaU = new double[l];
		for (int i = 0; i < l; i++) {
			thetaL[i] = Double.NEGATIVE_INFINITY;
		}
		for (int i = 0; i < l; i++) {
			thetaU[i] = Double.POSITIVE_INFINITY;
		}
		for (int i = 0; i < n; i++) {
			if (Y[i] == 0) {
				thetaL[0] = Math.max(thetaL[0], x[i]);
			}
			else if (Y[i] == l) {
				thetaU[l - 1] = Math.min(thetaU[l - 1], x[i]);
			}
			else {
				thetaL[Y[i]] = Math.max(thetaL[Y[i]], x[i]);
				thetaU[Y[i] - 1] = Math.min(thetaU[Y[i] - 1], x[i]);
			}
		}
		for (int i = 0; i < l; i++) {
			if (thetaL[i] == Double.NEGATIVE_INFINITY) {
				theta[i] = thetaU[i] - Math.abs(thetaU[i]) * 0.95;
			}
			else if (thetaU[i] == Double.POSITIVE_INFINITY) {
				theta[i] = thetaL[i] + Math.abs(thetaL[i]) * 0.95;
			}
			else {
				theta[i] = (thetaL[i] + thetaU[i]) / 2;
			}
		}
	}
}


