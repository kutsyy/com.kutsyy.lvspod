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
public class lnsrch {
	/**
	 *  Description of the Field
	 */
	public double f;
	/**
	 *  Description of the Field
	 */
	public int check;
	private final double ALF = 1.0e-4;
	private final double TOLX = 1.0e-7;
	private int n;
	private int i;
	private double a;
	private double alam;
	private double alam2 = 0;
	private double alamin;
	private double b;
	private double disc;
	private double f2 = 0;
	private double fold2 = 0;
	private double rhs1;
	private double rhs2;
	private double slope;
	private double sum;
	private double temp;
	private double test;
	private double tmplam;


	/**
	 *  put your documentation comment here
	 *
	 *@param  xold                                 Description of Parameter
	 *@param  fold                                 Description of Parameter
	 *@param  g                                    Description of Parameter
	 *@param  p                                    Description of Parameter
	 *@param  x                                    Description of Parameter
	 *@param  stpmax                               Description of Parameter
	 *@param  func                                 Description of Parameter
	 *@exception  java.lang.IllegalStateException  Description of Exception
	 *@exception  IllegalArgumentException         Description of Exception
	 */
	public lnsrch(double[] xold, double fold, double[] g, double[] p, double[] x,
			double stpmax, MvFunction func)
			 throws java.lang.IllegalStateException, IllegalArgumentException {
		lnsrch(xold, fold, g, p, x, stpmax, func);
	}


	/**
	 *  Description of the Method
	 *
	 *@param  xold                                 Description of Parameter
	 *@param  fold                                 Description of Parameter
	 *@param  g                                    Description of Parameter
	 *@param  p                                    Description of Parameter
	 *@param  x                                    Description of Parameter
	 *@param  stpmax                               Description of Parameter
	 *@param  func                                 Description of Parameter
	 *@exception  java.lang.IllegalStateException  Description of Exception
	 *@exception  IllegalArgumentException         Description of Exception
	 */
	public void lnsrch(double[] xold, double fold, double[] g, double[] p, double[] x,
			double stpmax, MvFunction func)
			 throws java.lang.IllegalStateException, IllegalArgumentException {
		if (xold.lenth != g.length || g.length != p.length || p.length != x.lnegth) {
			throw new IllegalArgumentException("dimentions must agree");
		}
		n = xold.length;
		check = 0;
		sum = 0.0;
		for (i = 0; i < n; i++) {
			sum += p[i] * p[i];
		}
		sum = Math.sqrt(sum);
		if (sum > stpmax) {
			for (i = 0; i < n; i++) {
				p[i] *= stpmax / sum;
			}
		}
		slope = 0.0;
		for (i = 0; i < n; i++) {
			slope += g[i] * p[i];
		}
		test = 0.0;
		for (i = 0; i < n; i++) {
			temp = Math.abs(p[i]) / Math.max(Math.abs(xold[i]), 1.0);
			if (temp > test) {
				test = temp;
			}
		}
		alamin = TOLX / test;
		alam = 1.0;
		while (true) {
			for (i = 0; i < n; i++) {
				x[i] = xold[i] + alam * p[i];
			}
			f = func.f(x);
			if (alam < alamin) {
				for (i = 0; i < n; i++) {
					x[i] = xold[i];
				}
				check = 1;
				return;
			}
			else if (f <= fold + ALF * alam * slope) {
				return;
			}
			else {
				if (alam == 1.0) {
					tmplam = -slope / (2.0 * (f - fold - slope));
				}
				else {
					rhs1 = f - fold - alam * slope;
					rhs2 = f2 - fold2 - alam2 * slope;
					a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
					b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (
							alam - alam2);
					if (a == 0.0) {
						tmplam = -slope / (2.0 * b);
					}
					else {
						disc = b * b - 3.0 * a * slope;
						if (disc < 0.0) {
							throw new java.lang.IllegalStateException("Roundoff problem in lnsrch.");
						}
						else {
							tmplam = (-b + Math.sqrt(disc)) / (3.0 * a);
						}
					}
					if (tmplam > 0.5 * alam) {
						tmplam = 0.5 * alam;
					}
				}
			}
			alam2 = alam;
			f2 = f;
			fold2 = fold;
			alam = Math.max(tmplam, 0.1 * alam);
		}
	}
}

