/*
 * put your module comment here
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 */

package com.kutsyy.lvspod;

import com.kutsyy.util.*;

/**
 *  Insert the type's description here. Creation date: (10/10/2000 4:50:42 PM)
 *
 *@author     vadim
 *@created    December 2, 2000
 *@author:
 */
abstract class OrdinalSpatialAbstractSubstitution extends OrdinalSpatialAbstract {

	/**
	 *  Insert the method's description here.
	 *
	 *@param  y  int[]
	 */
	public OrdinalSpatialAbstractSubstitution(int[] y) {
		super(y);
		x = new double[N];
	}


	/**
	 *  Insert the method's description here.
	 *
	 *@param  y  int[]
	 *@param  z  double[][]
	 */
	public OrdinalSpatialAbstractSubstitution(int[] y, double[][] z) {
		super(y, z);
		x = new double[N];
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
	public OrdinalSpatialAbstractSubstitution(int[] y, int[] xLocation, int[] yLocation,
			double[][] z, int[][] Neighborgood) {
		super(y, xLocation, yLocation, z, Neighborgood);
		x = new double[N];
	}


	/**
	 *  OrdinalSpatialAbstractSubstitution constructor comment.
	 */
	protected OrdinalSpatialAbstractSubstitution() {
		super();
	}


	/**
	 *  OrdinalSpatialAbstractSubstitution constructor comment.
	 *
	 *@param  y          int[]
	 *@param  xLocation  int[]
	 *@param  yLocation  int[]
	 *@param  z          double[][]
	 */
	protected OrdinalSpatialAbstractSubstitution(int[] y, int[] xLocation,
			int[] yLocation, double[][] z) {
		super(y, xLocation, yLocation, z);
		x = new double[N];
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
	protected OrdinalSpatialAbstractSubstitution(int[] y, double phi, double[] theta,
			com.kutsyy.lvspod.SpatialPoint[] loc, double[] beta, double[][] z,
			int[][] neighborDefinition) {
		super(y, phi, theta, loc, beta, z, neighborDefinition);
		x = new double[N];
	}


	/**
	 *  Insert the method's description here. Creation date: (10/10/2000 4:52:09
	 *  PM)
	 */
	protected abstract void generateX();


	/**
	 *  Insert the method's description here.
	 */
	protected void initializeL() {
		int L = theta.length;
		thetaL = new double[N];
		thetaU = new double[N];
		for (int i = 0; i < N; i++) {
			if (Y[i] == 0) {
				thetaU[i] = (theta[0]);
				thetaL[i] = Double.NEGATIVE_INFINITY;
			}
			else if (Y[i] == L) {
				thetaU[i] = Double.POSITIVE_INFINITY;
				thetaL[i] = theta[L - 1];
			}
			else {
				thetaL[i] = theta[Y[i] - 1];
				thetaU[i] = theta[Y[i]];
			}
		}
		generateX();
//		double l = 0;
//		for (int i = 0; i < N; i++) {
//			if (Loc[i].neighbors.length > 0) {
//				double s = La.sum(La.getArray(x, Loc[i].neighbors)) * phi / Loc[i].neighbors.length;
//				double s1 = Math.sqrt((Loc[i].neighbors.length / (double) NeighborDefinition[0].length) * (
//						1 - phi * phi / (double) Loc[i].neighbors.length));
//				l += Math.log(Cdf.nor(thetaU[i], s, s1) - Cdf.nor(thetaL[i],
//						s, s1));
//			}
//		}
//		return l;
	}


	/**
	 *  Description of the Method
	 *
	 *@param  i  Description of Parameter
	 *@return    Description of the Returned Value
	 */
	protected double ll(int i) {
		if (Loc[i].neighbors.length > 0) {
			double s = 0;
			double s1 = 1;
			if (Loc[i].neighbors.length > 0) {
				for (int j = 0; j < Loc[i].neighbors.length; j++) {
					s += x[Loc[i].neighbors[j]];
				}
				s *= phi / (double) Loc[i].neighbors.length;
				s1 = Math.sqrt((
						1 - phi * phi / (double) Loc[i].neighbors.length));
			}
			double l = Math.log(Cdf.nor(thetaU[i], s, s1) - Cdf.nor(thetaL[i],
					s, s1));
			return l;
		}
		return 0;
	}

}

