/*
 * put your module comment here
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 */

package com.kutsyy.lvspod;

import com.kutsyy.util.*;

/**
 *  Insert the type's description here. Created by <A
 *  href="http://www.kutsyy.com">Vadim Kutsyy</A> <BR>
 *
 *
 *@author     <A href="http://www.kutsyy.com">Vadim Kutsyy</A>
 *@created    December 1, 2000
 */
public final class OrdinalSpatialQL extends OrdinalSpatialAbstract
		 implements Function {
	private double sd1;
	private double[][] cov2 = new double[2][2];
	private double[][] cov3 = new double[3][3];
	private double[][] m = new double[3][3];
	private boolean makeThetaLU = true;
	//private Mvndstpack mvndstpack = new Mvndstpack();
	private LMvndstTread[] lMvndstTread;
        //private boolean[][] lock;
        //private double[][] value;


	//private Thread[] lThread ;

	/**
	 *  Insert the method's description here.
	 *
	 *@param  y  int[]
	 */
	public OrdinalSpatialQL(int[] y) {
		super(y);
		initialize();
	}


	/**
	 *  Insert the method's description here.
	 *
	 *@param  y  int[]
	 *@param  z  double[][]
	 */
	public OrdinalSpatialQL(int[] y, double[][] z) {
		super(y, z);
		initialize();
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
	public OrdinalSpatialQL(int[] y, int[] xLocation, int[] yLocation, double[][] z,
			int[][] Neighborgood) {
		super(y, xLocation, yLocation, z, Neighborgood);
		initialize();
	}


	/**
	 *  put your documentation comment here
	 */
	protected OrdinalSpatialQL() {
		super();
		initialize();
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
	protected OrdinalSpatialQL(int[] y, double phi, double[] theta, com.kutsyy.lvspod.SpatialPoint[] loc,
			double[] beta, double[][] z, int[][] neighborDefinition) {
		super(y, phi, theta, loc, beta, z, neighborDefinition);
		initialize();
	}


	/**
	 *  Insert the method's description here.
	 *
	 *@return                                              double
	 *@exception  java.lang.UnsupportedOperationException  Description of Exception
	 */
	public synchronized double lL() throws java.lang.UnsupportedOperationException {
		if (makeThetaLU) {
			makeThetaLU();
		}
		double l = 0;
		switch (NeighborDefinition[0].length) {
			default:
				cov3[0][0] = NeighborDefinition[0].length / (2.0 - 2.0 * phi * phi);
				cov3[0][1] = cov3[1][0] = cov3[0][2] = cov3[2][0] = phi * cov3[0][0];
				cov3[1][2] = cov3[2][1] = phi * cov3[0][1];
				cov3[1][1] = cov3[2][2] = (2 - phi * phi) * cov3[0][0];
			case 1:
				cov2[1][1] = cov2[0][0] = NeighborDefinition[0].length / (1 -
						phi * phi);
				cov2[0][1] = cov2[1][0] = cov2[0][0] * phi;
			case 0:
				sd1 = Math.sqrt((double) NeighborDefinition[0].length);
		}
		for (int k = NeighborDefinition[0].length; k >= 0; k--) {
			double[] lKeep = new double[(int) Math.pow(L + 1, k + 1)];
			if (k > 1) {
				int[] lLoc = new int[N];
				int[][] lNeigh = new int[lKeep.length][k + 1];
				boolean[] lExist = new boolean[lKeep.length];
				java.util.Arrays.fill(lLoc, -1);
				java.util.Arrays.fill(lExist, false);
				for (int i = 0; i < N; i++) {
					if (Loc[i].neighbors.length == k) {
						lLoc[i] = Y[i] * (int) Math.pow(L + 1, k);
						for (int i1 = 0; i1 < k; i1++) {
							lLoc[i] += Y[Loc[i].neighbors[i1]] * Math.pow(L + 1, k -
									i1 - 1);
						}
						if (!lExist[lLoc[i]]) {
							lExist[lLoc[i]] = true;
							lNeigh[lLoc[i]][0] = i;
							for (int j = 0; j < k; j++) {
								lNeigh[lLoc[i]][j + 1] = Loc[i].neighbors[j];
							}
						}
					}
				}
				for (int i = 0; i < lKeep.length; i++) {
					if (lExist[i]) {
						double[][] w = new double[k + 1][k + 1];
						double[][] m = new double[k + 1][k + 1];
						for (int j0 = 0; j0 <= k; j0++) {
							w[j0][j0] = 1;
						}
						m[0][0] = k;
						for (int j0 = 1; j0 <= k; j0++) {
							w[0][j0] = -phi / (double) k;
							w[j0][0] = -phi;
							m[j0][j0] = k;
						}
						lMvndstTread[i] = new LMvndstTread(
								La.getArray(thetaL, lNeigh[i]),
								La.getArray(thetaU, lNeigh[i]),
								La.times(La.solve(w), m),  EPS);
					}
				}
				for (int i = 0; i < lKeep.length; i++) {
					if (lExist[i]) {
						//lThread[i] = new Thread(lMvndstpack[i]);
						lMvndstTread[i].start();
					}
				}
				for (int i = 0; i < lKeep.length; i++) {
					if (lExist[i]) {
						while (lMvndstTread[i].lock) {
							try {
								wait(5);
							}
							catch (InterruptedException e) {
								/*
								 * do nothing
								 */
							}
						}
						lMvndstTread[i].yield();
						//System.out.println(i+" "+lMvndstpack[i].value);
					}
				}
				for (int i = 0; i < N; i++) {
					if (lLoc[i] > -1) {
						l += Math.log(lMvndstTread[lLoc[i]].mvndstpack.value);
                                                //System.out.println(l+" "+value[lLoc[i]][0]+" "+lMvndstTread[lLoc[i]].mvndstpack.value);
					}
				}
			}
			else {
				int thisl;
				boolean[] lNotDefined = new boolean[lKeep.length];

				java.util.Arrays.fill(lNotDefined, true);
				for (int i = 0; i < N; i++) {
					if (Loc[i].neighbors.length == k) {
						thisl = Y[i] * (int) Math.pow(L + 1, k);
						for (int i1 = 0; i1 < k; i1++) {
							thisl += Y[Loc[i].neighbors[i1]] * Math.pow(L + 1, k -
									i1 - 1);
						}
						if (lNotDefined[thisl]) {
							lNotDefined[thisl] = false;
							switch (k) {
//								case 2:
//									lKeep[thisl] = mvnor(La.getArray(thetaL,
//											Loc[i].fullneighbors), La.getArray(thetaU,
//											Loc[i].fullneighbors), cov3);
//									break;
								case 1:
									lKeep[thisl] = Cdf.bvnor(La.getArray(thetaL,
											Loc[i].fullneighbors), La.getArray(thetaU,
											Loc[i].fullneighbors), cov2);
									break;
								case 0:
									lKeep[thisl] = Cdf.nor(thetaU[i], 0, Math.sqrt(
											1 - phi * phi)) - Cdf.nor(thetaL[i],
											0, sd1);
									break;
								default:
                                                                throw new IllegalStateException("no suppose to be here");
//									double[][] w = new double[k + 1][k + 1];
//									double[][] m = new double[k + 1][k + 1];
//									for (int j0 = 0; j0 <= k; j0++) {
//										w[j0][j0] = 1;
//									}
//									m[0][0] = k;
//									for (int j0 = 1; j0 <= k; j0++) {
//										w[0][j0] = -phi / (double) k;
//										w[j0][0] = -phi;
//										m[j0][j0] = k;
//									}
//									lKeep[thisl] = mvnor(La.getArray(thetaL,
//											Loc[i].fullneighbors), La.getArray(thetaU,
//											Loc[i].fullneighbors), La.times(La.solve(w), m));
//
//									break;
							}
						}

						l += lKeep[thisl];
					}
				}
			}
		}
		//System.out.println(phi + " " + theta[0] + " " + theta[1] + " " + l);
		return l;
	}


	/**
	 *  Description of the Method
	 *
	 *@param  g  Description of Parameter
	 */
	public void grad(double[] g) {
		makeThetaLU();
		makeThetaLU = false;
		if (g.length > L) {
			throw new java.lang.UnsupportedOperationException("not implemented yet");
		}
		for (int i = 0; i < g.length; i++) {
			g[i] = 0;
		}
		switch (NeighborDefinition[0].length) {
			default:
				throw new java.lang.UnsupportedOperationException("not implemented yet");
			case 2:
				cov3[0][0] = NeighborDefinition[0].length / (2.0 - 2.0 * phi * phi);
				cov3[0][1] = cov3[1][0] = cov3[0][2] = cov3[2][0] = phi * cov3[0][0];
				cov3[1][2] = cov3[2][1] = phi * cov3[0][1];
				cov3[1][1] = cov3[2][2] = (2 - phi * phi) * cov3[0][0];
			case 1:
				cov2[1][1] = cov2[0][0] = NeighborDefinition[0].length / (1 -
						phi * phi);
				cov2[0][1] = cov2[1][0] = cov2[0][0] * phi;
			case 0:
				sd1 = Math.sqrt((double) NeighborDefinition[0].length);
		}
		for (int i = 0; i < N; i++) {
			switch (Loc[i].neighbors.length) {
				case 0:
					break;
				case 1:
					break;
				case 2:
					g[0] += Cdf.mvnor(La.getArray(thetaL, Loc[i].fullneighbors),
							La.getArray(thetaU, Loc[i].fullneighbors), cov3) * 2 * phi / (
							1 - phi * phi);
					break;
				default:
					throw new java.lang.UnsupportedOperationException("not implemented yet");
			}
		}
		makeThetaLU = true;
	}


	/**
	 *  Description of the Method
	 */
	protected void initialize() {
		int i = (int) Math.pow(L + 1, NeighborDefinition[0].length + 1);
		if (this.lMvndstTread == null || this.lMvndstTread.length < i) {
			lMvndstTread = new LMvndstTread[i];
                        //lock=new boolean[i][1];
                        //value=new double[i][1];
		}
	}


	/**
	 *  Description of the Method
	 *
	 *@param  low   Description of Parameter
	 *@param  high  Description of Parameter
	 *@param  sig   Description of Parameter
	 *@return       Description of the Returned Value
	 */
//	private final double mvnor(double[] low, double[] high, double[][] sig) {
//
//		if (low.length == 1) {
//			return Cdf.nor(high[0] / Math.sqrt(sig[0][0])) - Cdf.nor(low[0] / Math.sqrt(sig[0][0]));
//		}
//		if (low.length == 2) {
//			return Cdf.bvnor(low, high, sig);
//		}
//
//		for (int i = 0; i < low.length; i++) {
//			if (low[i] == Double.NEGATIVE_INFINITY && high[i] == Double.POSITIVE_INFINITY) {
//				return mvnor(La.removeI(low, i), La.removeI(high, i), La.removeI(sig,
//						i));
//			}
//		}
//		int max = 100 * low.length;
//		mvndstpack.Mvndstpack(low, high, null, sig, max, EPS, EPS);
//		while (mvndstpack.error > Math.min(EPS, Math.abs(mvndstpack.value) * EPS)
//				 && mvndstpack.maxpts < Integer.MAX_VALUE) {
//			mvndstpack.maxpts = mvndstpack.maxpts < Integer.MAX_VALUE / 10.0 ? mvndstpack.maxpts * 10 : Integer.MAX_VALUE;
//			//System.out.print(T.maxpts + " " + T.error + " " + T.value + " ");
//			mvndstpack.run();
//		}
//		return mvndstpack.value;
//	}


	/**
	 *  Description of the Class
	 *
	 *@author     Vadim Kutsyy
	 *@created    January 8, 2001
	 */
	private class LMvndstTread extends Thread {
		double[] low;
		double[] high;
		double[][] cor;
		double EPS;
		boolean lock;
		//double[] value;
		private Mvndstpack mvndstpack;


		/**
		 *  Constructor for the lMvndstTread object
		 *
		 *@param  Low    Description of Parameter
		 *@param  High   Description of Parameter
		 *@param  Cor    Description of Parameter
		 *@param  Eps    Description of Parameter
		 *@param  Lock   Description of Parameter
		 *@param  Value  Description of Parameter
		 */
		public LMvndstTread(double[] Low, double[] High, double[][] Cor, double Eps) {
			low = (double[]) Low.clone();
			high = (double[]) High.clone();
			cor = (double[][]) Cor.clone();
			EPS = Eps;
			lock = true;
			//value = Value;
		}


		/**
		 *  Main processing method for the lMvndstTread object
		 */
		public void run() {
			lock = true;
			this.mvndstpack = new Mvndstpack(low, high, cor, 100000, EPS, EPS);
			//value[0] = mvndstpack.value;
			lock = false;
		}
	}
}

