/*
 * put your module comment here
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 */

package com.kutsyy.lvspod;

import com.kutsyy.util.*;

/**
 *  Insert the type's description here. Creation date: (10/10/2000 6:02:53 PM)
 *
 *@author     vadim
 *@created    December 1, 2000
 *@author:
 */
public final class OrdinalSpatialMnA extends OrdinalSpatialAbstractSubstitution {

	/**
	 *  OrdinalSpatialMnA constructor comment.
	 */
	protected OrdinalSpatialMnA() {
		super();
	}


	/**
	 *  OrdinalSpatialMnA constructor comment.
	 *
	 *@param  y          int[]
	 *@param  xLocation  int[]
	 *@param  yLocation  int[]
	 *@param  z          double[][]
	 */
	protected OrdinalSpatialMnA(int[] y, int[] xLocation, int[] yLocation,
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
	protected OrdinalSpatialMnA(int[] y, double phi, double[] theta, com.kutsyy.lvspod.SpatialPoint[] loc,
			double[] beta, double[][] z, int[][] neighborDefinition) {
		super(y, phi, theta, loc, beta, z, neighborDefinition);
	}


	/**
	 *  Insert the method's description here. Creation date: (10/10/2000 6:02:53
	 *  PM)
	 */
	protected void generateX() {
		if (K > 0) {
			x = La.times(Z, beta);
			for (int i = 0; i < N; i++) {
				thetaL[i] += x[i];
				thetaU[i] += x[i];
			}
		}
		for (int i = 0; i < N; i++) {
			x[i] = Exp.ntrc(thetaL[i], thetaU[i]);
		}
	}

}


