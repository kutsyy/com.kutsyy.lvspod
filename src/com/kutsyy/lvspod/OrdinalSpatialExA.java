/*
 * put your module comment here
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 */

package com.kutsyy.lvspod;

import com.kutsyy.util.*;

/**
 *  Insert the type's description here. Creation date: (10/10/2000 6:08:50 PM)
 *
 *@author     vadim
 *@created    December 2, 2000
 *@author:
 */
public final class OrdinalSpatialExA extends OrdinalSpatialAbstract {

	/**
	 *  OrdinalSpatialPLEM constructor comment.
	 */
	protected OrdinalSpatialExA() {
		super();
	}


	/**
	 *  OrdinalSpatialPLEM constructor comment.
	 *
	 *@param  y          int[]
	 *@param  xLocation  int[]
	 *@param  yLocation  int[]
	 *@param  z          double[][]
	 */
	protected OrdinalSpatialExA(int[] y, int[] xLocation, int[] yLocation,
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
	protected OrdinalSpatialExA(int[] y, double phi, double[] theta, com.kutsyy.lvspod.SpatialPoint[] loc,
			double[] beta, double[][] z, int[][] neighborDefinition) {
		super(y, phi, theta, loc, beta, z, neighborDefinition);
	}


	/**
	 *  Insert the method's description here.
	 *
	 *@return    double
	 */
	public double lL() {
		int L = theta.length;
		double[] thetaL = new double[N];
		double[] thetaU = new double[N];
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
		if (K > 0) {
			x = La.times(Z, beta);
			//x.assign(Z.zMult(beta, null));
			for (int i = 0; i < N; i++) {
				thetaL[i] += x[i];
				thetaU[i] += x[i];
			}
		}
		//generateX();
		double l = 0;
		for (int i = 0; i < N; i++) {
			if (Loc[i].neighbors.length > 0) {
				System.out.println(i);
				double[][] s = new double[Loc[i].neighbors.length + 1][Loc[i].neighbors.length
						 + 1];
				double t1 = phi / Loc[i].neighbors.length;
				double t2 = t1 * t1;
				for (int j = 1; j <= Loc[i].neighbors.length; j++) {
					s[j][0] = s[0][j] = t1;
					for (int k = j + 1; k <= Loc[i].neighbors.length; k++) {
						s[j][k] = s[k][j] = t2;
					}
				}
				for (int j = 0; j <= Loc[i].neighbors.length; j++) {
					s[j][j] = 1;
				}
				double[] thetal = new double[Loc[i].neighbors.length];
				double[] thetau = new double[Loc[i].neighbors.length];
				for (int j = 0; j < Loc[i].neighbors.length; j++) {
					thetal[j] = thetaL[Loc[i].neighbors[j]];
					thetau[j] = thetaU[Loc[i].neighbors[j]];
				}
				double t = La.sum(Exp.mvntrc(thetal, thetau, La.solve(La.removeI(La.solve(s),
						0))));
				double st = Math.sqrt((Loc[i].neighbors.length / (double) NeighborDefinition[0].length) * (
						1 - phi * phi / (double) Loc[i].neighbors.length));
				l += Math.log(Cdf.nor(thetaU[i], t, st) - Cdf.nor(thetaL[i],
						t, st));
			}
		}
		return l;
	}
}


