package com.kutsyy.util;

//import Jama.*;
import java.util.Random;
import cern.jet.random.engine.*;
//import edu.cornell.lassp.houle.RngPack.*;
//import cern.jet.random.engine.MersenneTwister;
/**
 *  Collection of methods in this class, generate random variables. Some methods
 *  call other, but most are stand alone. <BR>
 *  <h2>Documentation to be added.</h2> <br>
 *  Requires <A href="http://nicewww.cern.ch/%7Ehoschek/colt/index.htm">Colt
 *  library</a> <br>
 *  The random engine is of mose 1e1000 <BR>
 *  Created by <A href="http://www.kutsyy.com">Vadim Kutsyy</A> <BR>
 *
 *
 *@author     <A href="http://www.kutsyy.com">Vadim Kutsyy</A>
 *@created    January 9, 2001
 */

public final class Rnd {
	//private static RandomElement rnd = new MersenneTwister(20);
	//private static edu.cornell.lassp.houle.RngPack.RandomElement rnd = new cern.jet.random.engine.MersenneTwister(new java.util.Date());
	//private static RandomElement rnd = new Ranmar();
	//private static Random rnd = new Random();

	private final static long[] zrng = {0, 1973272912, 281629770, 20006270, 1280689831, 2096730329, 1933576050, 913566091, 246780520, 1363774876, 604901985, 1511192140, 1259851944, 824064364, 150493284, 242708531, 75253171, 1964472944, 1202299975, 233217322, 1911216000, 726370533, 403498145, 993232223, 1103205531, 762430696, 1922803170, 1385516923, 76271663, 413682397, 726466604, 336157085, 1432650381, 1120463904, 595778810, 877722890, 1046574445, 68911991, 2088367019, 748545416, 622401386, 2122378830, 640690903, 1774806513, 2132545692, 2079249579, 78130110, 852776735, 1187867272, 1351423507, 1645973084, 1997049139, 922510944, 2045512870, 898585771, 243649545, 1004818771, 773686062, 403188473, 372279877, 1901633463, 498067494, 2087759558, 493157915, 597104727, 1530940798, 1814496276, 536444882, 1663153658, 855503735, 67784357, 1432404475, 619691088, 119025595, 880802310, 176192644, 1116780070, 277854671, 1366580350, 1142483975, 2026948561, 1053920743, 786262391, 1792203830, 1494667770, 1923011392, 1433700034, 1244184613, 1147297105, 539712780, 1545929719, 190641742, 1645390429, 264907697, 620389253, 1502074852, 927711160, 364849192, 2049576050, 638580085, 547070247};
	private final static long MODLUS = 2147483647;
	private final static long MULT1 = 24112;
	private final static long MULT2 = 26143;
	private static int remaining = 0;
	private static double remainder;
	private static MersenneTwister rnd = new MersenneTwister(new java.util.Date());



	/**
	 *  Be(alpha1, alpha2)
	 *
	 *@param  alpha1  alpha1
	 *@param  alpha2  alpha2
	 *@return         Be(alpha1, alpha2)
	 *@author:        <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double beta(double alpha1, double alpha2) {
		if (alpha1 == 1 && alpha2 == 1) {
			return unif();
		}
		double y1 = gamma(alpha1, 1);
		return y1 / (y1 + gamma(alpha2, 1));
	}



	/**
	 *  Be(alpha1, alpha2)
	 *
	 *@param  n       number of variable to generate
	 *@param  alpha1  alpha1
	 *@param  alpha2  alpha2
	 *@return         vector of Be(alpha1, alpha2)
	 *@author:        <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double[] beta(int n, double alpha1, double alpha2) {
		if (alpha1 == 1 && alpha2 == 1) {
			return unif(n);
		}
		double x[] = new double[n];
		for (int i = 0; i < n; i++) {
			x[i] = beta(alpha1, alpha2);
		}
		return x;
	}



	/**
	 *  Binomial(p)
	 *
	 *@param  p   p
	 *@return     Binomial(p)
	 *@author:    <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static boolean binomial(double p) {
		return unif() < p;
	}



	/**
	 *  Binomial(p)
	 *
	 *@param  n   number of variables to generate
	 *@param  p   p
	 *@return     vector of Binomial(p)
	 *@author:    <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static boolean[] binomial(int n, double p) {
		boolean x[] = new boolean[n];
		for (int i = 0; i < n; i++) {
			x[i] = unif() < p;
		}
		return x;
	}



	/**
	 *  Binomial(p)
	 *
	 *@param  p   p
	 *@return     Binomial(p)
	 *@author:    <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static int binomial_int(double p) {
		return unif() < p ? 0 : 1;
	}



	/**
	 *  Binomial(p)
	 *
	 *@param  n   number of variables to generate
	 *@param  p   p
	 *@return     vector of Binomial(p)
	 *@author:    <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static int[] binomial_int(int n, int p) {
		int x[] = new int[n];
		for (int i = 0; i < n; i++) {
			x[i] = unif() < p ? 0 : 1;
		}
		return x;
	}



	/**
	 *  Exponential(lambda)
	 *
	 *@param  lambda  lambda
	 *@return         Exponential(lambda)
	 *@author:        <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double exp(double lambda) {
		return -lambda * Math.log(unif());
	}



	/**
	 *  Exponential(lambda)
	 *
	 *@param  n       number of variables to generate
	 *@param  lambda  lambda
	 *@return         vector of Exponential(lambda)
	 *@author:        <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double[] exp(int n, double lambda) {
		double x[] = new double[n];
		for (int i = 0; i < n; i++) {
			x[i] = -lambda * Math.log(unif());
		}
		return x;
	}



	/**
	 *  Gamma(alpha,beta)
	 *
	 *@param  alpha  alpha
	 *@param  beta   beta
	 *@return        Gamma(alpha,beta)
	 *@author:       <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double gamma(double alpha, double beta) {
		if (alpha < 1) {
			double b = (Math.E + alpha) / Math.E;
			double y;
			double p;
			do {
				do {
					double u1 = unif();
					p = b * u1;
					if (p > 1) {
						break;
					}
					y = Math.pow(p, 1 / alpha);
					double u2 = unif();
					if (u2 <= Math.exp(-y)) {
						return beta * y;
					}
				} while (true);
				y = -Math.log((beta - p) / alpha);
				double u2 = unif();
				if (u2 <= Math.pow(y, alpha - 1)) {
					return beta * y;
				}
			} while (true);
		}
		else {
			double a = 1 / Math.sqrt(2 * alpha - 1);
			double b = alpha - Math.log(4);
			double q = alpha + 1 / a;
			double theta = 4.5;
			double d = 1 + Math.log(theta);
			do {
				double u1 = unif();
				double u2 = unif();
				double v = a * Math.log(u1 / (1 - u2));
				double y = alpha * Math.exp(v);
				double z = u1 * u1 * u2;
				double w = b + q * v + y;
				if (w + d - theta * z >= 0 && w >= Math.log(z)) {
					return beta * y;
				}
			} while (true);
		}
	}



	/**
	 *  Gamma(alpha,beta)
	 *
	 *@param  n      number of variables to generate
	 *@param  alpha  alpha
	 *@param  beta   beta
	 *@return        vector Gamma(alpha,beta)
	 *@author:       <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double[] gamma(int n, double alpha, double beta) {
		double x[] = new double[n];
		for (int i = 0; i < n; i++) {
			x[i] = gamma(alpha, beta);
		}
		return x;
	}



	/**
	 *  Geometric(p)
	 *
	 *@param  p   p
	 *@return     Geometric(p)
	 *@author:    <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static int geometric(double p) {
		return (int) Math.floor(Math.log(unif()) / Math.log(1 - p));
	}



	/**
	 *  Geometric(p)
	 *
	 *@param  n   number of variables to generate
	 *@param  p   p
	 *@return     vector of Geometric(p)
	 *@author:    <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static int[] geometric(int n, double p) {
		double a = 1 / Math.log(1 - p);
		int x[] = new int[n];
		for (int i = 0; i < n; i++) {
			x[i] = (int) Math.floor(Math.log(unif()) * a);
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 11:45:20 AM)
	 *
	 *@param  mu     double
	 *@param  sigma  double
	 *@return        double
	 *@author:       <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double lognor(double mu, double sigma) {
		return Math.exp((nor() + mu) * sigma);
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 11:46:29 AM)
	 *
	 *@param  n      number of variables to generate
	 *@param  mu     double
	 *@param  sigma  double
	 *@return        double[]
	 *@author:       <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double[] lognor(int n, double mu, double sigma) {
		double x[] = new double[n];
		for (int i = 0; i < n; i++) {
			x[i] = Math.exp((nor() + mu) * sigma);
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 1:09:14 PM)
	 *
	 *@param  sigma  Description of Parameter
	 *@return        double[]
	 *@author:       <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double[] mvnor(double[][] sigma) {
		double sqrtSigma[][] = La.choldc(sigma);
		;
		return mvnor_sqrtSigma(sqrtSigma);
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 1:32:45 PM)
	 *
	 *@param  mu     double
	 *@param  sigma  double[][]
	 *@return        double[]
	 *@author:       <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double[] mvnor(double[] mu, double[][] sigma) {
		double x[] = mvnor(sigma);
		for (int i = 0; i < sigma.length; i++) {
			x[i] += mu[i];
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 1:32:45 PM)
	 *
	 *@param  mu     double
	 *@param  sigma  double[][]
	 *@return        double[]
	 *@author:       <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double[] mvnor(double mu, double[][] sigma) {
		double x[] = mvnor(sigma);
		for (int i = 0; i < sigma.length; i++) {
			x[i] += mu;
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 1:09:14 PM)
	 *
	 *@param  n      number of variables to generate
	 *@param  sigma  Description of Parameter
	 *@return        double[]
	 *@author:       <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double[][] mvnor(int n, double[][] sigma) {
		double x[][] = new double[sigma.length][n];
		double sqrtSigma[][] = La.choldc(sigma);
		for (int i = 0; i < n; i++) {
			double y[] = mvnor_sqrtSigma(sqrtSigma);
			for (int j = 0; j < y.length; j++) {
				x[j][i] = y[j];
			}
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 1:09:14 PM)
	 *
	 *@param  n      number of variables to generate
	 *@param  mu     Description of Parameter
	 *@param  sigma  Description of Parameter
	 *@return        double[]
	 *@author:       <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double[][] mvnor(int n, double mu[], double[][] sigma) {
		double x[][] = new double[sigma.length][n];
		double sqrtSigma[][] = La.choldc(sigma);
		for (int i = 0; i < n; i++) {
			double y[] = mvnor_sqrtSigma(sqrtSigma);
			for (int j = 0; j < y.length; j++) {
				x[j][i] = y[j] + mu[j];
			}
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 1:09:14 PM)
	 *
	 *@param  n      number of variables to generate
	 *@param  mu     Description of Parameter
	 *@param  sigma  Description of Parameter
	 *@return        double[]
	 *@author:       <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double[][] mvnor(int n, double mu, double[][] sigma) {
		double x[][] = new double[sigma.length][n];
		double sqrtSigma[][] = La.choldc(sigma);
		for (int i = 0; i < n; i++) {
			double y[] = mvnor_sqrtSigma(sqrtSigma);
			for (int j = 0; j < y.length; j++) {
				x[j][i] = y[j] + mu;
			}
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 1:07:50 PM)
	 *
	 *@param  sqrtSigma  double[][]
	 *@return            double[]
	 *@author:           <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double[] mvnor_sqrtSigma(double[][] sqrtSigma) {
		int n = sqrtSigma.length;
		double y[] = nor(n);
		double x[] = new double[n];
		for (int i = 0; i < n; i++) {
			x[i] = 0;
			for (int j = 0; j <= i; j++) {
				x[i] += sqrtSigma[i][j] * y[j];
			}
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (2/10/00 3:03:10 PM)
	 *
	 *@param  low    double[]
	 *@param  high   double[]
	 *@param  sigma  Description of Parameter
	 *@return        double[]
	 */
	public static double[] mvntrc(double[] low, double[] high, double[][] sigma) {
		return mvntrcChol(low, high, La.choldc(sigma));
	}



	/**
	 *  Insert the method's description here. Creation date: (2/10/00 3:05:54 PM)
	 *
	 *@param  n      number of variables to generate
	 *@param  low    double[]
	 *@param  high   double[]
	 *@param  sigma  double[][]
	 *@return        double[][]
	 */
	public static double[][] mvntrc(
			int n,
			double[] low,
			double[] high,
			double[][] sigma) {
		double sqrtSigma[][] = La.choldc(sigma);
		int m = sigma.length;
		double[][] t = new double[n][];
		main :
		for (int j = 0; j < n; j++) {
			local :
			do {
				t[j] = (double[]) mvnor_sqrtSigma(sqrtSigma).clone();
				for (int i = 0; i < n; i++) {
					if (t[j][i] < low[i] || t[j][i] > high[i]) {
						continue local;
					}
				}
				continue main;
			} while (true);
		}
		return La.t(t);
	}



	/**
	 *  Insert the method's description here. Creation date: (2/10/00 3:03:10 PM)
	 *
	 *@param  low        double[]
	 *@param  high       double[]
	 *@param  sqrtSigma  Description of Parameter
	 *@return            double[]
	 */
	public static double[] mvntrcChol(
			double[] low,
			double[] high,
			double[][] sqrtSigma) {
		int n = sqrtSigma.length;
		main :
		do {
			double[] t = mvnor_sqrtSigma(sqrtSigma);
			for (int i = 0; i < n; i++) {
				if (t[i] < low[i] || t[i] > high[i]) {
					continue main;
				}
			}
			return t;
		} while (true);
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 12:42:53 PM)
	 *
	 *@param  s   int
	 *@param  p   double
	 *@return     int
	 *@author:    <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static int negbinomial(int s, double p) {
		int x = 0;
		for (int i = 0; i < s; i++) {
			x += geometric(p);
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 12:46:09 PM)
	 *
	 *@param  n   number of variables to generate
	 *@param  s   int
	 *@param  p   double
	 *@return     int[]
	 *@author:    <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static int[] negbinomial(int n, int s, double p) {
		int x[] = new int[n];
		for (int i = 0; i < n; i++) {
			x[i] = negbinomial(s, p);
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (1/24/00 9:00:19 AM)
	 *
	 *@return    double
	 */
	public static double nor() {
		return rnd.gaussian();
		//return gauss();
		//return rnd.nextGaussian();
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 11:01:06 AM)
	 *
	 *@param  mu  double
	 *@param  sd  double
	 *@return     double
	 *@author:    <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double nor(double mu, double sd) {
		return (nor() + mu) * sd;
	}



	/**
	 *  Insert the method's description here. Creation date: (1/25/00 4:25:28 PM)
	 *
	 *@param  n   number of variables to generate
	 *@return     double[]
	 *@author:    <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double[] nor(int n) {
		double x[] = new double[n];
		for (int i = 0; i < n; i++) {
			x[i] = nor();
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 11:02:44 AM)
	 *
	 *@param  n   int
	 *@param  mu  double
	 *@param  sd  double
	 *@return     double[]
	 *@author:    <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double[] nor(int n, double mu, double sd) {
		double x[] = new double[n];
		for (int i = 0; i < n; i++) {
			x[i] = (nor() + mu) * sd;
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (3/6/00 11:42:39 AM)
	 *
	 *@param  n  number of variables to generate
	 *@return    double[]
	 */
	public static double[] norOrdered(int n) {
		double x[] = nor(n);
		java.util.Arrays.sort(x);
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (2/10/00 4:30:30 PM)
	 *
	 *@param  low  double
	 *@return      double
	 */
	public static double ntrc(double low) {
		if (Double.isNaN(low)) {
			throw new IllegalArgumentException("low is NaN");
		}
		double cut = .45;
		if (low > cut) {
			double z = -Math.log(unif()) / low;
			while (unif() > Math.exp(-.5 * z * z)) {
				z = -Math.log(unif()) / low;
			}
			if (z < 0) {
				throw new IllegalStateException("z<0");
			}
			return z + low;
		}
		double x = nor();
		while (x < low) {
			x = nor();
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (2/10/00 2:31:07 PM)
	 *
	 *@param  low   double
	 *@param  high  double
	 *@return       double
	 */
	public static double ntrc(double low, double high) {
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
		if (high - low > 0.2 &&
				(!((high > 1.5 && low > 1.5) || (high < 1.5 && low < 1.5)))) {
			double aa = 0;
			if (high > 0 && low < 0) {
				aa = Pdf.nor(0);
			}
			else {
				aa = Pdf.nor(low);
				aa = Math.max(aa, Pdf.nor(high));
			}
			do {
				double y = unif(low, high);
				double u = unif();
				if (u * aa <= Pdf.nor(y)) {
					return y;
				}
			} while (true);
		}
		return CdfInv.nor(unif(Cdf.nor(low), Cdf.nor(high)));
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
	public static double ntrc(double low, double high, double mu, double sd) {
		low -= mu;
		high -= mu;
		low /= sd;
		high /= sd;
		return (ntrc(low, high) + mu) * sd;
	}



	/**
	 *  Insert the method's description here. Creation date: (2/10/00 2:36:04 PM)
	 *
	 *@param  n     number of variables to generate
	 *@param  low   double[]
	 *@param  high  double[]
	 *@return       double[]
	 */
	public static double[] ntrc(int n, double low, double high) {
		double[] t = new double[n];
		for (int i = 0; i < n; i++) {
			t[i] = ntrc(low, high);
		}
		return t;
	}



	/**
	 *  Insert the method's description here. Creation date: (2/10/00 2:35:15 PM)
	 *
	 *@param  n     number of variables to generate
	 *@param  low   double[]
	 *@param  high  double[]
	 *@param  mu    double[]
	 *@param  sd    double[][]
	 *@return       double[]
	 */
	public static double[] ntrc(
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
	 *  Insert the method's description here. Creation date: (1/26/00 12:47:13 PM)
	 *
	 *@param  lambda  double
	 *@return         int
	 *@author:        <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static int poisson(double lambda) {
		double a = Math.exp(-lambda);
		double b = 1;
		int i = 0;
		do {
			b *= unif();
			if (b < a) {
				return i;
			}
			i++;
		} while (true);
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 12:50:34 PM)
	 *
	 *@param  n       number of variables to generate
	 *@param  lambda  double
	 *@return         int[]
	 *@author:        <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static int[] poisson(int n, double lambda) {
		double a = Math.exp(-lambda);
		int x[] = new int[n];
		for (int i = 0; i < n; i++) {
			double b = 1;
			int j = 0;
			l1 :
			do {
				b *= unif();
				if (b < a) {
					x[i] = j;
					break l1;
				}
				j++;
			} while (true);
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (2/10/00 2:15:05 PM)
	 *
	 *@param  c  double
	 *@return    double
	 */
	public static double triang(double c) {
		double u = unif();
		return (u <= c) ? Math.sqrt(c * u) : 1.0 - Math.sqrt((1.0 - c) * (1.0 - u));
	}



	/**
	 *  Insert the method's description here. Creation date: (2/10/00 2:17:32 PM)
	 *
	 *@param  a  double
	 *@param  b  double
	 *@param  c  double
	 *@return    double
	 */
	public static double triang(double a, double b, double c) {
		return (triang(c) - a) / (b - a);
	}



	/**
	 *  Insert the method's description here. Creation date: (2/10/00 2:18:21 PM)
	 *
	 *@param  n  number of variables to generate
	 *@param  c  double
	 *@return    double[]
	 */
	public static double[] triang(int n, double c) {
		double[] t = new double[n];
		for (int i = 0; i < n; i++) {
			t[i] = triang(c);
		}
		return t;
	}



	/**
	 *  Insert the method's description here. Creation date: (2/10/00 2:17:52 PM)
	 *
	 *@param  n  number of variables to generate
	 *@param  a  double
	 *@param  b  double
	 *@param  c  double
	 *@return    double[]
	 */
	public static double[] triang(int n, double a, double b, double c) {
		double[] t = new double[n];
		for (int i = 0; i < n; i++) {
			t[i] = triang(a, b, c);
		}
		return t;
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 11:03:40 AM)
	 *
	 *@return     double
	 *@author:    <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double unif() {
		return rnd.raw();
		//return rnd.nextDouble();
		//return rand(1);
	}



	/**
	 *  Insert the method's description here. Creation date: (3/6/00 4:04:19 PM)
	 *
	 *@param  low   double[]
	 *@param  high  double[]
	 *@return       double
	 */
	public final static double unif(double[] low, double[] high) {
		int n = low.length;
		double[] cut = new double[n];
		cut[0] = Math.max(0, high[0] - low[0]);
		for (int i = 1; i < n; i++) {
			cut[i] = cut[i - 1] + Math.max(0, high[i] - low[i]);
		}
		double x = unif(0, cut[n - 1]);
		if (x < cut[0]) {
			return x + low[0];
		}
		for (int i = 0; i < n - 1; i++) {
			if (x < cut[i + 1]) {
				return (x - cut[i]) + low[i + 1];
			}
		}
		return Double.NaN;
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 11:06:43 AM)
	 *
	 *@param  a   double
	 *@param  b   double
	 *@return     double
	 *@author:    <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double unif(double a, double b) {
		if (a > b) {
			double tmp = a;
			a = b;
			b = a;
		}
		if (a == b) {
			return a;
		}
		return (unif()) * (b - a) + a;
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 11:08:33 AM)
	 *
	 *@param  n   number of variables to generate
	 *@return     double[]
	 *@author:    <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double[] unif(int n) {
		double x[] = new double[n];
		for (int i = 0; i < n; i++) {
			x[i] = unif();
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (1/26/00 11:09:26 AM)
	 *
	 *@param  n   number of variables to generate
	 *@param  a   double
	 *@param  b   double
	 *@return     double[]
	 *@author:    <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double[] unif(int n, double a, double b) {
		if (a > b) {
			double tmp = a;
			a = b;
			b = a;
		}
		double x[] = new double[n];
		for (int i = 0; i < n; i++) {
			x[i] = (unif()) * (b - a) + a;
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (3/6/00 12:57:13 PM)
	 *
	 *@param  m  int
	 *@return    int
	 */
	public static int unifDiscrete(int m) {
		return (int) Math.floor(unif() * (m + 1));
	}



	/**
	 *  generate single weibull(eta,betta) random variable; Creation date: (1/26/00
	 *  1:39:51 PM)
	 *
	 *@param  eta   double - scale parameter
	 *@param  beta  double - shape parameter
	 *@return       double - random variables
	 *@author:      <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double weibull(double eta, double beta) {
		return eta * Math.pow(-Math.log(1 - unif()), 1 / beta);
	}



	/**
	 *  generate array(n) if weibull(eta,betta) random variables; Creation date:
	 *  (1/26/00 1:39:51 PM)
	 *
	 *@param  n     number of variables to generate
	 *@param  eta   double - scale parameter
	 *@param  beta  double - shape parameter
	 *@return       double[n] - random variables
	 *@author:      <Vadum Kutsyy, kutsyy@hotmail.com>
	 */
	public static double[] weibull(int n, double eta, double beta) {
		double x[] = new double[n];
		for (int i = 0; i < n; i++) {
			x[i] = eta * Math.pow(-Math.log(1 - unif()), 1 / beta);
		}
		return x;
	}



	/**
	 *  Insert the method's description here. Creation date: (3/3/00 9:49:44 AM)
	 *
	 *@return    double
	 */
	private static double gauss() {
		double u1;
		double u2;
		if (remaining == 0) {
			u1 = rand(1);
			u2 = rand(1);
			remainder = Math.cos(2.0 * Math.PI * u2) * Math.sqrt(-2.0 * Math.log(u1));
			remaining = 1;
			return Math.sin(2.0 * Math.PI * u2) * Math.sqrt(-2.0 * Math.log(u1));
		}
		else {
			remaining = 0;
			return remainder;
		}
	}



	/**
	 *  Insert the method's description here. Creation date: (3/3/00 9:43:11 AM)
	 *
	 *@param  stream  int
	 *@return         double
	 */
	private static double rand(int stream) {
		long zi;
		long lowprd;
		long hi31;
		zi = zrng[stream];
		lowprd = (zi & 65535) * MULT1;
		hi31 = (zi >> 16) * MULT1 + (lowprd >> 16);
		zi = ((lowprd & 65535) - MODLUS) + ((hi31 & 32767) << 16) + (hi31 >> 15);
		if (zi < 0) {
			zi += MODLUS;
		}
		lowprd = (zi & 65535) * MULT2;
		hi31 = (zi >> 16) * MULT2 + (lowprd >> 16);
		zi = ((lowprd & 65535) - MODLUS) + ((hi31 & 32767) << 16) + (hi31 >> 15);
		if (zi < 0) {
			zi += MODLUS;
		}
		zrng[stream] = zi;
		return ((zi >> 7 | 1) + 1) / 16777216.0;
	}
}
