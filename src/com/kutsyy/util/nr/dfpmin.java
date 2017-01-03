/*
 * put your module comment here
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 */

package com.kutsyy.util.nr;

import com.kutsyy.util.*;

/**
 *  Title: com.kutsyy.util Description: Mathematical and Statistical Utilities.
 *  Requires <a href="http://tilde-hoschek.home.cern.ch/~hoschek/colt/index.htm">
 *  Colt</a> and <a href="http://math.nist.gov/javanumerics/jama/">Jama</a>
 *  libraries. Copyright: Copyright (c) 2000 Company: The University of Michigan
 *
 *@author     <a href="http://www.kutsyy.com">Vadim Kutsyy</a>
 *@created    December 1, 2000
 *@version    1.0
 */
public final class dfpmin {
	/**
	 *  Description of the Field
	 */
	public int iter;
	/**
	 *  Description of the Field
	 */
	public double fret;
	private int check;
	private int i;
	private int its;
	private int j;
	private double den;
	private double fac;
	private double fad;
	private double fae;
	private double fp;
	private double stpmax;
	private double sum = 0.0;
	private double sumdg;
	private double sumxi;
	private double temp;
	private double test;
	private double[] g;
	private double[] hdg;
	private double[] pnew;
	private double[] xi;
	private double[] dg;
	private double[][] hessin;
	private final static int ITMAX = 2000;
	private final static double EPS = 3e-8;
	private final static double TOLX = 4 * EPS;
	private final static int STPMX = 100;


	/**
	 *  put your documentation comment here
	 */
	public dfpmin() {
	}


	/**
	 *  put your documentation comment here
	 *
	 *@param  p
	 *@param  gtol
	 *@param  fun
	 *@exception  java.lang.IllegalStateException  Description of Exception
	 */
	public dfpmin(double[] p, double gtol, MvFunction fun) throws java.lang.IllegalStateException {
		dfpmin(p, gtol, fun);
	}


	/**
	 *  Description of the Method
	 *
	 *@param  p                                    Description of Parameter
	 *@param  gtol                                 Description of Parameter
	 *@param  fun                                  Description of Parameter
	 *@exception  java.lang.IllegalStateException  Description of Exception
	 */
	public void dfpmin(double[] p, double gtol, MvFunction fun)
			 throws java.lang.IllegalStateException {
		if (p.length != n) {
			n = p.length;
			g = new double[n];
			hdg = new double[n];
			pnew = new double[n];
			xi = new double[n];
			dg = new double[n];
			hessin = new double[n][n];
		}
		fp = fun.f(p);
		fun.g(p, g);
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				hessin[i][j] = 0.0;
			}
			hessin[i][i] = 1.0;
			xi[i] = -g[i];
			sum += p[i] * p[i];
		}
		stpmax = STPMX * Math.max(Math.sqrt(sum), n);
		for (its = 0; its < ITMAX; its++) {
			iter = its;
			Lnsrch.lnsrch(p, fp, g, xi, pnew, stpmax, fun);
			fret = Lnsrch.f;
			check = Lnsrch.check;
			fp = fret;
			for (i = 0; i < n; i++) {
				xi[i] = pnew[i] - p[i];
				p[i] = pnew[i];
			}
			test = 0.0;
			for (i = 0; i < n; i++) {
				temp = Math.abs(xi[i]) / Math.max(Math.abs(p[i]), 1.0);
				if (temp > test) {
					test = temp;
				}
			}
			if (test < TOLX) {
				return;
			}
			for (i = 0; i < n; i++) {
				dg[i] = g[i];
			}
			fun.g(p, g);
			test = 0.0;
			den = Math.max(fret, 1.0);
			for (i = 0; i < n; i++) {
				temp = Math.abs(g[i]) * Math.max(Math.abs(p[i]), 1.0) / den;
				if (temp > test) {
					test = temp;
				}
			}
			System.out.print(test + " ");
			if (test < gtol) {
				return;
			}
			for (i = 0; i < n; i++) {
				dg[i] = g[i] - dg[i];
			}
			for (i = 0; i < n; i++) {
				hdg[i] = 0.0;
				for (j = 0; j < n; j++) {
					hdg[i] += hessin[i][j] * dg[j];
				}
			}
			fac = fae = sumdg = sumxi = 0.0;
			for (i = 0; i < n; i++) {
				fac += dg[i] * xi[i];
				fae += dg[i] * hdg[i];
				sumdg += Math.sqrt(dg[i]);
				sumxi += Math.sqrt(xi[i]);
			}
			if (fac * fac > EPS * sumdg * sumxi) {
				fac = 1.0 / fac;
				fad = 1.0 / fae;
				for (i = 0; i < n; i++) {
					dg[i] = fac * xi[i] - fad * hdg[i];
				}
				for (i = 0; i < n; i++) {
					for (j = 0; j < n; j++) {
						hessin[i][j] += fac * xi[i] * xi[j] - fad * hdg[i] * hdg[j]
								 + fae * dg[i] * dg[j];
					}
				}
			}
			for (i = 0; i < n; i++) {
				xi[i] = 0.0;
				for (j = 0; j < n; j++) {
					xi[i] -= hessin[i][j] * g[j];
				}
			}
		}
		throw new java.lang.IllegalStateException("too many iterations in dfpmin");
	}
}

