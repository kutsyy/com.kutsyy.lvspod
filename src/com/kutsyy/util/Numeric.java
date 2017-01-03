/*
 * put your module comment here
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 */

package com.kutsyy.util;

/**
 *  Numeric class include some maximization and simmular numeric functions. <BR>
 *  Created by <A href="http://www.kutsyy.com">Vadim Kutsyy</A> <BR>
 *
 *
 *@author     <A href="http://www.kutsyy.com">Vadim Kutsyy</A>
 *@created    December 1, 2000
 */
public final class Numeric {
	//private static int count=0;
	private static boolean min = true;
	private static Function local_f;
	private static MvFunction local_mvF;
	private final static double tol = 1e-10;
	private final static int ITMAX = 10000;
	private final static double CGOLD = 0.3819660;
	private final static double ZEPS = 1e-20;


	/**
	 *  Function Minimizaton based on Brent method
	 *
	 *@param  ax   double - lower limit for x
	 *@param  bx   double - initial gess for x
	 *@param  cx   double -upper limit for x
	 *@param  _f   kutsyy.util.Function - inteface with function f(double) to
	 *      minimize
	 *@param  tol  double - tolerance
	 *@return      double - value of x with minimize function f
	 */
	public static double brent(double ax, double bx, double cx, Function _f,
			double tol) {
		return brent(ax, bx, cx, _f, tol, "min");
	}


	/**
	 *  Function Minimizaton based on Brent method
	 *
	 *@param  ax       double - lower limit for x
	 *@param  bx       double - initial gess for x
	 *@param  cx       double -upper limit for x
	 *@param  _f       kutsyy.util.Function - inteface with function f(double) to
	 *      minimize
	 *@param  tol      double - tolerance
	 *@param  min_max  java.lang.String - "max" or "min" -to maximize or minimize
	 *@return          double - value of x with minimize function f
	 */
	public static double brent(double ax, double bx, double cx, Function _f,
			double tol, String min_max) {
		local_f = _f;
		min = !min_max.substring(0, 3).toLowerCase().equals("max");
		return brent_(ax, bx, cx, tol);
	}


	/**
	 *  Function Minmizaton based on golden search
	 *
	 *@param  ax   double - lower limit for x
	 *@param  bx   double - initial gess for x
	 *@param  cx   double -upper limit for x
	 *@param  _f   kutsyy.util.Function - inteface with function f(double) to
	 *      minimize
	 *@param  tol  double - tolerance
	 *@return      double - value of x with minimize function f
	 */
	public static double golden(double ax, double bx, double cx, Function _f,
			double tol) {
		return golden(ax, bx, cx, _f, tol, "min");
	}


	/**
	 *  Function Minmizaton based on golden search
	 *
	 *@param  ax       double - lower limit for x
	 *@param  bx       double - initial gess for x
	 *@param  cx       double -upper limit for x
	 *@param  _f       kutsyy.util.Function - inteface with function f(double) to
	 *      minimize
	 *@param  tol      double - tolerance
	 *@param  min_max  java.lang.String - "max" or "min" -to maximize or minimize
	 *@return          double - value of x with minimize function f
	 */
	public static double golden(double ax, double bx, double cx, Function _f,
			double tol, String min_max) {
		local_f = _f;
		min = min_max.substring(0, 3).toLowerCase().equals("min");
		return golden_(ax, bx, cx, tol);
	}


	/**
	 *  Insert the method's description here.
	 *
	 *@param  ax       double
	 *@param  bx       double
	 *@param  cx       double
	 *@param  _f       com.kutsyy.util.Function
	 *@param  tol      double
	 *@param  min_max  java.lang.String
	 *@return          double
	 */
	public static double minimize(double ax, double bx, double cx, Function _f,
			double tol, String min_max) {
		local_f = _f;
		min = min_max.substring(0, 3).toLowerCase().equals("min");
		double[] x = {
				ax, bx, cx
				};
		double[] fx = {
				f(ax), f(bx), f(cx)
				};
		double fb = fx[1];
		while (Math.abs(x[0] - x[2]) > tol) {
			if (La.min(fx) > fb) {
				throw new IllegalStateException("something is wrong");
			}
			double xx = x[1] + CGOLD * ((x[1] - x[0] > x[2] - x[1] ? x[0] : x[2])
					 - x[1]);
			double ff = f(xx);
			int m = -1;
			double t = ff;
			for (int i = 0; i < 2; i++) {
				if (fx[i] < t) {
					m = i;
					t = fx[i];
				}
			}
			if (m == -1) {
				if (xx < x[1]) {
					x[2] = x[1];
					fx[2] = fx[1];
					x[1] = xx;
					fx[1] = ff;
				}
				else {
					x[0] = x[1];
					fx[0] = fx[1];
					x[1] = xx;
					fx[1] = ff;
				}
			}
			else if (m == 0) {
				if (xx < x[1]) {
					x[2] = x[1];
					fx[2] = fx[1];
					x[1] = xx;
					fx[1] = ff;
				}
				else {
					x[2] = xx;
					fx[2] = ff;
				}
			}
			else if (m == 2) {
				if (xx > x[1]) {
					x[0] = x[1];
					fx[0] = fx[1];
					x[1] = xx;
					fx[1] = ff;
				}
				else {
					x[0] = xx;
					fx[0] = ff;
				}
			}
			else if (xx < x[1]) {
				x[0] = xx;
				fx[0] = ff;
			}
			else {
				x[2] = xx;
				fx[2] = ff;
			}
		}
		return x[1];
	}


	/**
	 *  Multivariate minimizaton based on discrete grid. Next point is found by
	 *  fitting 2 degree multivariate polinomial, and maximizing it.
	 *
	 *@param  x        double[] - imput initial gues, output, value that
	 *      minimize/maximize function
	 *@param  a        double[] - lower limit for x;
	 *@param  b        double[] - upper limit for x;
	 *@param  mvF      kutsyy.util.MvFunction -inteface that contain f(double[]) to
	 *      be minimaze/maximize
	 *@param  tol      double[] - tolerance for each dimention
	 *@param  min_max  java.lang.String - "min" or "max" for minimizaton, or
	 *      maximization.
	 *@return          double - value of function at the minimum/maximum;
	 */
	public static double[] mvGrid(double[] x, double[] a, double[] b, MvFunction mvF,
			double[] tol, String min_max) {
		local_mvF = mvF;
		min = min_max.substring(0, 3).toLowerCase().equals("min");
		int n = x.length;
		if (a.length != n || b.length != n || tol.length != n) {
			throw new IllegalArgumentException("Matrix dimentions must agree");
		}
		int N = (int) Math.pow(3, n);
		double[][] xi = new double[N][n];
		double[] fx = new double[N];
		for (int k = 0; k < n; k++) {
			int I = (int) Math.pow(3, n - k - 1);
			int J = (int) Math.pow(3, k);
			for (int j = 0; j < J; j++) {
				for (int i = 0; i < I; i++) {
					xi[i + j * 3 * I][k] = a[k];
					xi[i + I + j * 3 * I][k] = x[k];
					xi[i + I + I + j * 3 * I][k] = b[k];
				}
			}
		}
		for (int i = 0; i < N; i++) {
			fx[i] = mvF(xi[i]);
		}
		boolean cont = false;
		main :
		do {
			//System.out.print(Numeric.min(fx)+" ");
			//LaIo.println(x);
			double[][] X = new double[N][2 * n + n * (n - 1) / 2 + 1];
			for (int i = 0; i < N; i++) {
				int m = 0;
				for (int j = 0; j < n; j++) {
					X[i][j] = xi[i][j];
					X[i][j + n] = xi[i][j] * xi[i][j];
					for (int k = j + 1; k < n; k++) {
						X[i][2 * n + m + k - j - 1] = xi[i][j] * xi[i][k];
					}
					m += n - j - 1;
					X[i][X[i].length - 1] = 1;
				}
			}
			double[] xOld = (double[]) x.clone();
			double[] beta = new double[N];
			try {
				beta = La.solveLS(X, fx);
			}
			catch (IllegalArgumentException e1) {
				int i = 0;
				double t = Double.POSITIVE_INFINITY;
				for (int j = 0; j < fx.length; j++) {
					if (fx[i] < t) {
						t = fx[i];
						i = j;
					}
				}
				x = xi[i];
				break main;
			}
			;
			X = null;
			X = new double[n][n];
			double[] c = new double[n];
			int m = 0;
			for (int i = 0; i < n; i++) {
				c[i] = -beta[i];
				X[i][i] = 2 * beta[i + n];
			}
			for (int i = 0; i < n; i++) {
				for (int k = i + 1; k < n; k++) {
					X[k][i] += beta[2 * n + m + k - i - 1];
					X[i][k] += X[k][i];
				}
				m += n - i - 1;
			}
			try {
				x = La.solve(X, c);
			}
			catch (RuntimeException e1) {
				int i = 0;
				double t = Double.POSITIVE_INFINITY;
				for (int j = 0; j < fx.length; j++) {
					if (fx[i] < t) {
						t = fx[i];
						i = j;
					}
				}
				x = xi[i];
				break main;
			}
			;
			X = null;
			double ff = Double.POSITIVE_INFINITY;
			loop :
			 {
				for (int i = 0; i < n; i++) {
					if (x[i] < a[i] || x[i] > b[i]) {
						break loop;
					}
				}
				ff = mvF(x);
			}
			double[] xx = (double[]) x.clone();
			double fff = ff;
			int[] M = new int[n];
			for (int i = 0; i < fx.length; i++) {
				if (ff > fx[i]) {
					ff = fx[i];
					for (int j = 0; j < n; j++) {
						M[j] = 2;
						if (xi[i][j] == a[j]) {
							M[j] = 1;
						}
						else if (xi[i][j] == b[j]) {
							M[j] = 3;
						}
					}
				}
			}
			for (int i = 0; i < n; i++) {
				switch (M[i]) {
					case 0:
						if (x[i] != xOld[i]) {
							if (x[i] <= a[i]) {
								b[i] = xOld[i];
								x[i] = (a[i] + b[i]) / 2.0;
							}
							else if (x[i] >= b[i]) {
								a[i] = xOld[i];
								x[i] = (a[i] + b[i]) / 2.0;
							}
							else if (x[i] < xOld[i]) {
								b[i] = xOld[i];
							}
							else {
								a[i] = xOld[i];
							}
						}
						break;
					case 1:
						if (x[i] > a[i]) {
							b[i] = Math.min(x[i], xOld[i]);
						}
						else {
							b[i] = xOld[i];
						}
						if (x[i] >= b[i] || x[i] <= a[i]) {
							x[i] = (a[i] + b[i]) / 2.0;
						}
						break;
					case 3:
						if (x[i] < b[i]) {
							a[i] = Math.max(x[i], xOld[i]);
						}
						else {
							a[i] = xOld[i];
						}
						if (x[i] >= b[i] || x[i] <= a[i]) {
							x[i] = (a[i] + b[i]) / 2.0;
						}
						break;
					case 2:
						if (x[i] < xOld[i]) {
							if (a[i] >= x[i]) {
								a[i] += (xOld[i] - a[i]) * 0.1;
							}
							else {
								a[i] = x[i];
							}
						}
						else if (b[i] <= x[i]) {
							b[i] -= (b[i] - xOld[i]) * 0.1;
						}
						else {
							b[i] = x[i];
						}
						x[i] = xOld[i];
						break;
				}
				if (a[i] >= x[i] || x[i] >= b[i]) {
					throw new IllegalArgumentException("x is not between a nd b");
				}
			}
			for (int i = 0; i < n; i++) {
				if (!(x[i] > a[i] && x[i] < b[i])) {
					x[i] = (a[i] + b[i]) / 2;
				}
			}
			cont = false;
			for (int i = 0; i < n; i++) {
				if (Math.max(Math.abs(x[i] - xOld[i]), Math.abs(a[i] - b[i])) > tol[i]) {
					cont = true;
					break;
				}
			}
			if (!cont) {
				break main;
			}
			double[][] xiOld = new double[xi.length][];
			for (int i = 0; i < xi.length; i++) {
				xiOld[i] = (double[]) xi[i].clone();
			}
			for (int k = 0; k < n; k++) {
				int I = (int) Math.pow(3, n - k - 1);
				int J = (int) Math.pow(3, k);
				for (int j = 0; j < J; j++) {
					for (int i = 0; i < I; i++) {
						xi[i + j * 3 * I][k] = a[k];
						xi[i + I + j * 3 * I][k] = x[k];
						xi[i + I + I + j * 3 * I][k] = b[k];
					}
				}
			}
			double[] fxOld = (double[]) fx.clone();
			if (La.min(fx) < mvF(x) - 1) {
				throw new IllegalStateException("x is not a min");
			}
			big :
			for (int k = 0; k < fx.length; k++) {
				small :
				for (int j = 0; j < fx.length; j++) {
					for (int i = 0; i < n; i++) {
						if (xi[k][i] != xiOld[j][i]) {
							continue small;
						}
					}
					fx[k] = fxOld[j];
					continue big;
				}
				small :
				 {
					for (int i = 0; i < n; i++) {
						if (xi[k][i] != xx[i]) {
							break small;
						}
					}
					fx[k] = fff;
					continue big;
				}
				fx[k] = mvF(xi[k]);
			}
			fxOld = null;
			xiOld = null;
		} while (true);
		return x;
	}


	/**
	 *  Multivariate minimizaton based on discrete grid. Next point is found by
	 *  fitting 2 degree multivariate polinomial, and maximizing it.
	 *
	 *@param  x        double[] - imput initial gues, output, value that
	 *      minimize/maximize function
	 *@param  a        double[] - lower limit for x;
	 *@param  b        double[] - upper limit for x;
	 *@param  mvF      kutsyy.util.MvFunction -inteface that contain f(double[]) to
	 *      be minimaze/maximize
	 *@param  tol      double - tolerance
	 *@param  min_max  java.lang.String - "min" or "max" for minimizaton, or
	 *      maximization.
	 *@return          double - value of function at the minimum/maximum;
	 */
	public static double[] mvGrid(double[] x, double[] a, double[] b, MvFunction mvF,
			double tol, String min_max) {
		double[] Tol = new double[a.length];
		for (int i = 0; i < a.length; i++) {
			Tol[i] = tol;
		}
		return mvGrid(x, a, b, mvF, Tol, min_max);
	}




	/**
	 *  Insert the method's description here. Creation date: (2/13/00 6:09:40 PM)
	 *
	 *@param  ax   double
	 *@param  bx   double
	 *@param  cx   double
	 *@param  tol  double
	 *@return      double
	 */
	private static double brent_(double ax, double bx, double cx, double tol) {
		double[] xx = {
				ax, bx, cx
				};
		//double[] ff = {f(ax), f(bx), f(cx)};
		double[] ff = new double[3];
		mnbrac_(xx, ff);
		ax = xx[0];
		bx = xx[1];
		cx = xx[2];
		int iter;
		double a;
		double b;
		double d = 0;
		double etemp;
		double fu;
		double fv = 0;
		double fw;
		double fx;
		double p;
		double q;
		double r;
		double tol1;
		double tol2;
		double
				u;
		double v;
		double w;
		double x;
		double xm;
		double e = 0;
		a = ax;
		b = cx;
		x = w = v = bx;
		fw = fx = fx = ff[1];
		for (iter = 1; iter <= ITMAX; iter++) {
			xm = .5 * (a + b);
			tol2 = 2 * (tol1 = tol * Math.abs(x) + ZEPS);
			if (Math.abs(x - xm) <= (tol2 - .5 * (b - a)) && iter > 1) {
				return x;
			}
			if (Math.abs(e) > tol1) {
				r = (x - w) * (fx - fv);
				q = (x - v) * (fx - fw);
				p = (x - v) * q - (x - w) * r;
				q = 2.0 * (q - r);
				if (q > 0) {
					p = -p;
				}
				q = Math.abs(q);
				etemp = e;
				e = d;
				if (Math.abs(p) >= Math.abs(.5 * q * etemp) || p <= q * (a - x) ||
						p >= q * (b - x)) {
					d = CGOLD * (e = (x >= xm ? a - x : b - x));
				}
				else {
					d = p / q;
					u = x + d;
					if (u - a < tol2 || b - u < tol2) {
						d = sign(tol1, xm - x);
					}
				}
			}
			else {
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			}
			u = (Math.abs(d) >= tol1 ? x + d : x + sign(tol1, d));
			fu = f(u);
			if (fu <= fx) {
				if (u >= x) {
					a = x;
				}
				else {
					b = x;
				}
				v = w;
				w = x;
				x = u;
				fv = fw;
				fw = fx;
				fx = fu;
			}
			else {
				if (u < x) {
					a = u;
				}
				else {
					b = x;
				}
				if (fu <= fw || w == x) {
					v = w;
					w = u;
					fv = fw;
					fw = fu;
				}
				else if (fu <= fv || v == x || v == w) {
					v = u;
					fv = fu;
				}
			}
		}
		throw new IllegalStateException("Too many itterations");
	}


	/**
	 *  Insert the method's description here. Creation date: (3/16/00 5:00:58 PM)
	 *
	 *@param  n  int
	 *@param  F  Description of Parameter
	 *@return    Description of the Returned Value
	 */
	private static double expQgauss(Function F, int n) {
		double[] x = new double[n];
		double[] w = new double[n];
		gauher(x, w);
		double[] y = new double[n];
		for (int i = 0; i < n; i++) {
			y[i] = F.f(x[i]);
		}
		double z1 = 0;
		double z2 = 0;
		for (int i = 0; i < n; i++) {
			z1 += w[i] * y[i];
			z2 += w[i] * x[i] * y[i];
		}
		return z2 / z1;
	}


	/**
	 *  Insert the method's description here. Creation date: (2/13/00 7:08:35 PM)
	 *
	 *@param  x  double
	 *@return    double
	 */
	private static double f(double x) {
		return min ? local_f.f(x) : -local_f.f(x);
	}


	/**
	 *  Insert the method's description here. Creation date: (3/16/00 4:35:14 PM)
	 *
	 *@param  x  double[]
	 *@param  w  double[]
	 */
	private static void gauher(double[] x, double[] w) {
		if (x.length != w.length) {
			throw new IllegalArgumentException("length of the vertors do no agree");
		}
		final double EPS = 3e-14;
		//final double PIM4 = 0.7511255444649425;
		final double PIM4 = 1 / Math.sqrt(Math.sqrt(Math.PI));
		final int MAXIT = 100;
		int n = x.length;
		int its;
		int m = (int) (n + 1) / 2;
		double p1;
		double p2;
		double p3;
		double pp = 0;
		double z = 0;
		double z1;
		for (int i = 0; i < m; i++) {
			switch (i) {
				case 0:
					z = Math.sqrt((double) (2 * n + 1)) - 1.85575 * Math.pow((double) (
							2 * n + 1), -0.16667);
					break;
				case 1:
					z -= 1.14 * Math.pow((double) n, 0.426) / z;
					break;
				case 2:
					z = 1.86 * z - 0.86 * x[0];
					break;
				case 3:
					z = 1.91 * z - 0.91 * x[1];
					break;
				default:
					z = 2.0 * z - x[i - 2];
					break;
			}
			do {
				p1 = PIM4;
				p2 = 0;
				for (int j = 1; j <= n; j++) {
					p3 = p2;
					p2 = p1;
					p1 = z * Math.sqrt(2.0 / j) * p2 - Math.sqrt(((double) (j - 1)) / j) * p3;
				}
				pp = Math.sqrt((double) 2 * n) * p2;
				z1 = z;
				z = z1 - p1 / pp;
			} while (Math.abs(z - z1) >= EPS);
			x[i] = z;
			x[n - i - 1] = -z;
			w[i] = 2.0 / (pp * pp);
			w[n - i - 1] = w[i];
		}
	}


	/**
	 *  Insert the method's description here. Creation date: (2/13/00 6:09:40 PM)
	 *
	 *@param  ax   double
	 *@param  bx   double
	 *@param  cx   double
	 *@param  tol  double
	 *@return      double
	 */
	private static double golden_(double ax, double bx, double cx, double tol) {
		double r = 0.61803399;
		double c = 1 - r;
		double f1;
		double f2;
		double x0;
		double x1;
		double x2;
		double x3;
		double[] x = {
				ax, bx, cx
				};
		double[] fx = new double[3];
		mnbrac_(x, fx);
		ax = x[0];
		bx = x[1];
		cx = x[2];
		x0 = ax;
		x3 = cx;
		if (Math.abs(cx - bx) > Math.abs(bx - ax)) {
			x1 = bx;
			f1 = fx[1];
			x2 = bx + c * (cx - bx);
			f2 = f(x2);
		}
		else {
			x2 = bx;
			f2 = fx[2];
			x1 = bx - c * (bx - ax);
			f1 = f(x1);
		}
		while (Math.abs(x3 - x0) > tol * (Math.abs(x1) + Math.abs(x2))) {
			if (f2 < f1) {
				x0 = x1;
				x1 = x2;
				x2 = r * x1 + c * x3;
				f1 = f2;
				f2 = f(x2);
			}
			else {
				x3 = x2;
				x2 = x1;
				x1 = r * x2 + c * x0;
				f2 = f1;
				f1 = f(x1);
			}
		}
		if (f1 < f2) {
			return x1;
		}
		else {
			return x2;
		}
	}


	/**
	 *  Insert the method's description here. Creation date: (3/21/00 10:50:39 PM)
	 *
	 *@param  x  double[]
	 *@return    double
	 */
	private static double min(double[] x) {
		double t = x[0];
		for (int i = 1; i < x.length; i++) {
			t = t < x[i] ? t : x[i];
		}
		return t;
	}


	/**
	 *  Insert the method's description here. Note Minimazation Creation date:
	 *  (2/13/00 5:20:57 PM)
	 *
	 *@param  x   double[]
	 *@param  fx  double[]
	 */
	private static void mnbrac_(double[] x, double[] fx) {
		//System.out.println(count++);
		if (Double.isNaN(x[0]) || Double.isNaN(x[1]) || Double.isNaN(x[2])) {
			throw new IllegalArgumentException("x is NAN");
		}
		if (Double.isInfinite(x[0]) || Double.isInfinite(x[1]) || Double.isInfinite(x[2])) {
			throw new IllegalArgumentException("x is Infinite");
		}
		if (x[0] >= x[1] || x[1] >= x[2]) {
			throw new IllegalArgumentException("x is not ordered");
		}
		double low = x[0];
		double high = x[2];
		x[0] = x[1] + (x[0] - x[1]) * CGOLD;
		x[2] = (x[2] - x[1]) * CGOLD + x[1];
		fx[0] = f(x[0]);
		fx[1] = f(x[1]);
		fx[2] = f(x[2]);
		int count = 0;
		while (!((fx[1] < fx[0] && fx[1] < fx[2]) || Math.abs(x[0] - x[2]) < tol * 2.0 * Math.abs(x[1]))
				 && count < ITMAX) {
			count++;
			double xx = (-fx[0] * x[1] + fx[0] * x[2] - fx[1] * x[2] + x[1] * fx[2]
					 - x[0] * fx[2] + x[0] * fx[1]) / (x[2] * x[2] * fx[0] - fx[1] * x[2] * x[2]
					 - x[1] * x[1] * fx[0] + fx[1] * x[0] * x[0] - fx[2] * x[0] * x[0] +
					fx[2] * x[1] * x[1]) / 2.0;
			//double xx = x[1] - ((x[1] - x[0]) * (x[1] - x[0]) * (fx[1] - fx[2]) - (x[1] - x[2]) * (x[1] - x[2]) * (fx[1] - x[0])) / ((x[1] - x[0]) * (fx[1] - fx[2]) - (x[1] - x[2]) * (fx[1] - fx[0]));
			if (xx >= high) {
				xx = (high - x[2]) * CGOLD + x[2];
			}
			if (xx <= low) {
				xx = x[0] + (low - x[0]) * CGOLD;
			}
			if (Double.isNaN(xx) || Double.isInfinite(xx) || xx == x[0]) {
				xx = (Rnd.unif() > .5) ? (x[1] + x[2]) / 2.0 : (x[0] + x[1]) / 2.0;
			}
			double ff = f(xx);
			if (ff < fx[1]) {
				if (xx > x[1]) {
					x[0] = x[1];
					fx[0] = fx[1];
					if (xx > x[2]) {
						x[1] = x[2];
						fx[1] = fx[2];
						x[2] = xx;
						fx[2] = ff;
					}
					else {
						x[1] = xx;
						fx[1] = xx;
					}
				}
				else {
					x[2] = x[1];
					fx[2] = fx[1];
					if (xx < x[0]) {
						x[1] = x[0];
						fx[1] = fx[0];
						x[0] = xx;
						fx[0] = ff;
					}
					else {
						x[1] = xx;
						fx[1] = xx;
					}
				}
			}
			else {
				if (xx > x[1]) {
					x[2] = xx;
					fx[2] = ff;
				}
				else {
					x[0] = xx;
					fx[0] = ff;
				}
			}
		}
	}


	/**
	 *  Insert the method's description here. Creation date: (3/19/00 8:08:16 PM)
	 *
	 *@param  x  double[]
	 *@return    double
	 */
	private static double mvF(double[] x) {
		try {
			return min ? local_mvF.f(x) : -local_mvF.f(x);
		}
		catch (Exception e) {
			e.printStackTrace();
			return 0;
		}
	}


	/**
	 *  Insert the method's description here. Creation date: (3/16/00 5:00:58 PM)
	 *
	 *@param  n  int
	 *@param  F  Description of Parameter
	 *@return    Description of the Returned Value
	 */
	private static double qgauss(Function F, int n) {
		//local_f=_f;
		double[] x = new double[n];
		double[] w = new double[n];
		gauher(x, w);
		double z = 0;
		for (int i = 0; i < n; i++) {
			z += w[i] * F.f(x[i]);
		}
		return z;
	}


	/**
	 *  Insert the method's description here. Creation date: (2/14/00 1:48:26 PM)
	 *
	 *@param  a  double
	 *@param  b  double
	 *@return    double
	 */
	private static double sign(double a, double b) {
		return a == 0 ? 0 : a > 0 ? b : -b;
	}
}


