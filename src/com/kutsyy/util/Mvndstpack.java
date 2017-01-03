
package com.kutsyy.util;
import cern.jet.random.engine.*;

/**
 *  This class implements multivariate normal integration. <BR>
 *  Method is based on MCMC. class use its own random number generator, <BR>
 *  and does not depends an any other classes <BR>
 *  so no sinchronization is required for multithreading. <BR>
 *  <BR>
 *  Description from original FORTRAN code: <BR>
 *  A subroutine for computing multivariate normal probabilities. <BR>
 *  This subroutine uses an algorithm given in the paper <BR>
 *  "Numerical Computation of Multivariate Normal Probabilities", in <BR>
 *  J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by <BR>
 *  Alan Genz <BR>
 *  Department of Mathematics <BR>
 *  Washington State University <BR>
 *  Pullman, WA 99164-3113 <BR>
 *  Email : AlanGenz@wsu.edu <BR>
 *  Parameters N INTEGER, the number of variables. <BR>
 *  LOWER REAL, array of lower integration limits. <BR>
 *  UPPER REAL, array of upper integration limits. <BR>
 *  INFIN INTEGER, array of integration limits flags: <BR>
 *  if INFIN(I) < 0, Ith limits are (-infinity, infinity); <BR>
 *  if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)]; <BR>
 *  if INFIN(I) = 1, Ith limits are [LOWER(I), infinity); <BR>
 *  if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)]. <BR>
 *  CORREL REAL, array of correlation coefficients; the correlation <BR>
 *  coefficient in row I column J of the correlation matrix <BR>
 *  should be stored in CORREL( J + ((I-2)*(I-1))/2 ), for J < I. <BR>
 *  THe correlation matrix must be positive semidefinite. <BR>
 *  MAXPTS INTEGER, maximum number of function values allowed. This <BR>
 *  parameter can be used to limit the time. A sensible <BR>
 *  strategy is to start with MAXPTS = 1000*N, and then <BR>
 *  increase MAXPTS if ERROR is too large. <BR>
 *  ABSEPS REAL absolute error tolerance. <BR>
 *  RELEPS REAL relative error tolerance. <BR>
 *  ERROR REAL estimated absolute error, with 99% confidence level. <BR>
 *  VALUE REAL estimated value for the integral <BR>
 *  INFORM INTEGER, termination status parameter: <BR>
 *  if INFORM = 0, normal completion with ERROR < EPS; <BR>
 *  if INFORM = 1, completion with ERROR >EPS and MAXPTS <BR>
 *  function vaules used; increase MAXPTS to <BR>
 *  decrease ERROR; <BR>
 *  if INFORM = 2, N >100 or N < 1. <BR>
 *
 *
 *@author     kutsyy
 *@created    January 14, 2001
 */
public final class Mvndstpack {
    /**
     *  Description of the Field
     */
    public boolean running = false;
    /**
     *  Output : error of the estimator
     */
    public double error;
    /**
     *  Output: code of the computation (see description)
     */
    public int inform;
    /**
     *  Output: value of the integral
     */
    public double value = 0;
    /**
     *  Input: relative error
     */
    public double releps;
    /**
     *  Input: number of itteration
     */
    public int maxpts;
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private double[] lower;
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private double[] sd;
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private double[] upper;
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private int[] infin;
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private double[] correl;
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private double abseps;
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private int nl = 100;
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private double[] cov = new double[nl * (nl + 1) / 2];
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private double[] a = new double[nl];
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private double[] b = new double[nl];
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private int[] infi = new int[nl];
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private int n;
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private final int[][] c = {
            {
            12, 9, 9, 13, 12, 12, 12, 12, 12, 12, 12, 12, 3, 3, 3, 12, 7, 7,
            12
            }, {
            13, 11, 17, 10, 15, 15, 15, 15, 15, 15, 22, 15, 15, 6, 6, 6, 15,
            15, 9
            }, {
            27, 28, 10, 11, 11, 20, 11, 11, 28, 13, 13, 28, 13, 13, 13, 14,
            14, 14, 14
            }, {
            35, 27, 27, 36, 22, 29, 29, 20, 45, 5, 5, 5, 21, 21, 21, 21, 21,
            21, 21
            }, {
            64, 66, 28, 28, 44, 44, 55, 67, 10, 10, 10, 10, 10, 10, 38, 38,
            10, 10, 10
            }, {
            111, 42, 54, 118, 20, 31, 31, 72, 17, 94, 14, 14, 11, 14, 14, 14,
            94, 10, 10
            }, {
            163, 154, 83, 43, 82, 92, 150, 59, 76, 76, 47, 11, 11, 100, 131,
            116, 116, 116, 116
            }, {
            246, 189, 242, 102, 250, 250, 102, 250, 280, 118, 196, 118, 191,
            215, 121, 121, 49, 49, 49
            }, {
            347, 402, 322, 418, 215, 220, 339, 339, 339, 337, 218, 315, 315,
            315, 315, 167, 167, 167, 167
            }, {
            505, 220, 601, 644, 612, 160, 206, 206, 206, 422, 134, 518, 134,
            134, 518, 652, 382, 206, 158
            }, {
            794, 325, 960, 528, 247, 247, 338, 366, 847, 753, 753, 236, 334,
            334, 461, 711, 652, 381, 381
            }, {
            1189, 888, 259, 1082, 725, 811, 636, 965, 497, 497, 1490, 1490,
            392, 1291, 508, 508, 1291, 1291, 508
            }, {
            1763, 1018, 1500, 432, 1332, 2203, 126, 2240, 1719, 1284, 878,
            1983, 266, 266, 266, 266, 747, 747, 127
            }, {
            2872, 3233, 1534, 2941, 2910, 393, 1796, 919, 446, 919, 919, 1117,
            103, 103, 103, 103, 103, 103, 103
            }, {
            4309, 3758, 4034, 1963, 730, 642, 1502, 2246, 3834, 1511, 1102,
            1102, 1522, 1522, 3427, 3427, 3928, 915, 915
            }, {
            6610, 6977, 1686, 3819, 2314, 5647, 3953, 3614, 5115, 423, 423,
            5408, 7426, 423, 423, 487, 6227, 2660, 6227
            }, {
            9861, 3647, 4073, 2535, 3430, 9865, 2830, 9328, 4320, 5913, 10365,
            8272, 3706, 6186, 7806, 7806, 7806, 8610, 2563
            }, {
            10327, 7582, 7124, 8214, 9600, 10271, 10193, 10800, 9086, 2365,
            4409, 13812, 5661, 9344, 9344, 10362, 9344, 9344, 8585
            }, {
            19540, 19926, 11582, 11113, 24585, 8726, 17218, 419, 4918, 4918,
            4918, 15701, 17710, 4037, 4037, 15808, 11401, 19398, 25950
            }, {
            34566, 9579, 12654, 26856, 37873, 38806, 29501, 17271, 3663, 10763,
            18955, 1298, 26560, 17132, 17132, 4753, 4753, 8713, 18624
            }, {
            31929, 49367, 10982, 3527, 27066, 13226, 56010, 18911, 40574, 20767,
            20767, 9686, 47603, 47603, 11736, 11736, 41601, 12888,
            32948
            }, {
            40701, 69087, 77576, 64590, 39397, 33179, 10858, 38935, 43129,
            35468, 35468, 2196, 61518, 61518, 27945, 70975, 70975,
            86478, 86478
            }, {
            103650, 125480, 59978, 46875, 77172, 83021, 126904, 14541, 56299,
            43636, 11655, 52680, 88549, 29804, 101894, 113675, 48040,
            113675, 34987
            }, {
            165843, 90647, 59925, 189541, 67647, 74795, 68365, 167485, 143918,
            74912, 167289, 75517, 8148, 172106, 126159, 35867, 35867,
            35867, 121694
            }, {
            130365, 236711, 110235, 125699, 56483, 93735, 234469, 60549, 1291,
            93937, 245291, 196061, 258647, 162489, 176631, 204895,
            73353, 172319, 28881
            }
            };
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private int hisum;
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private int np;
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private double olds = 0;
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private final int[] p = {
            31, 47, 73, 113, 173, 263, 397, 593, 907, 1361, 2053, 3079, 4621, 6947,
            10427, 15641, 23473, 35221, 52837, 79259, 118891, 178349, 267523,
            401287, 601942
            };
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private final double[] prime = {
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
            67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131,
            137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
            197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263,
            269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337,
            347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409
            };
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private double[] psqt;
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private int sampls;
    //new int[plim][klim - 1];
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private double varest;

    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private int[] infis = new int[1];
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private double[] e = new double[1];
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private double[] d = new double[1];
    /**
     *  Random Number Generator (faster than Java's build in)
     */
    private MersenneTwister rnd = new MersenneTwister(new java.util.Date());
    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private int[] N;


    /**
     *  Constructor
     *
     *@param  Lower   Lower limit of integration
     *@param  Upper   Upper limit of integration
     *@param  Cor     Covariance matrix
     *@param  Maxpts  maximum number of itteration
     *@param  Abseps  absolute error
     *@param  Releps  relative arror
     */
    public Mvndstpack(double[] Lower, double[] Upper, double[][] Cor, int Maxpts,
            double Abseps, double Releps) {
        Mvndstpack(Lower, Upper, null, Cor, Maxpts, Abseps, Releps);
    }


    /**
     *  Constructor
     *
     *@param  Lower   Lower limit of integration
     *@param  Upper   Upper limit of integration
     *@param  Maxpts  maximum number of itteration
     *@param  Abseps  absolute error
     *@param  Releps  relative arror
     *@param  Correl  Description of Parameter
     */
    public Mvndstpack(double[] Lower, double[] Upper, double[] Correl, int Maxpts,
            double Abseps, double Releps) {
        Mvndstpack(Lower, Upper, null, Correl, Maxpts, Abseps, Releps);
    }


    /**
     *  Constructor
     *
     *@param  Lower   Lower limit of integration
     *@param  Upper   Upper limit of integration
     *@param  Cor     Covariance matrix
     *@param  Infin   as above
     *@param  Maxpts  maximum number of itteration
     *@param  Abseps  absolute error
     *@param  Releps  relative arror
     */
    public Mvndstpack(double[] Lower, double[] Upper, int[] Infin, double[][] Cor,
            int Maxpts, double Abseps, double Releps) {
        Mvndstpack(Lower, Upper, Infin, Cor, Maxpts, Abseps, Releps);
    }



    /**
     *  Empty constructor
     */
    public Mvndstpack() {
    }


    /**
     *  Constructor
     *
     *@param  Lower   Lower limit of integration
     *@param  Upper   Upper limit of integration
     *@param  Correl  Covariance matrix
     *@param  Infin   as above
     *@param  Maxpts  maximum number of itteration
     *@param  Abseps  absolute error
     *@param  Releps  relative arror
     */
    public Mvndstpack(double[] Lower, double[] Upper, int[] Infin, double[] Correl,
            int Maxpts, double Abseps, double Releps) {
        Mvndstpack(Lower, Upper, Infin, Correl,
                Maxpts, Abseps, Releps);
    }


    /**
     *  main method
     *
     *@param  Lower   Lower limit of integration
     *@param  Upper   Upper limit of integration
     *@param  Infin   as above
     *@param  Maxpts  maximum number of itteration
     *@param  Abseps  absolute error
     *@param  Releps  relative arror
     *@param  Cov     Description of Parameter
     */
    public void Mvndstpack(double[] Lower, double[] Upper, int[] Infin, double[][] Cov,
            int Maxpts, double Abseps, double Releps) {
        if (sd == null || sd.length != Lower.length) {
            sd = new double[Lower.length];
            lower = new double[Lower.length];
            upper = new double[Lower.length];
        }
        for (int i = 0; i < Lower.length; i++) {
            sd[i] = Cov[i][i] == 1 ? 1 : Math.sqrt(Cov[i][i]);
        }
        if (correl == null || correl.length != Cov.length * (Cov.length - 1) / 2) {
            correl = new double[Cov.length * (Cov.length - 1) / 2];
        }
        int ij = 0;
        for (int i = 0; i < Cov.length; i++) {
            for (int j = 0; j < i; j++) {
                correl[ij++] = Cov[i][j] / (sd[i] * sd[j]);
            }
            upper[i] = Upper[i] / sd[i];
            lower[i] = Lower[i] / sd[i];

        }
        Mvndstpack(lower, upper, Infin, correl, Maxpts, Abseps, Releps);
    }


    /**
     *  main method
     *
     *@param  Lower   Lower limit of integration
     *@param  Upper   Upper limit of integration
     *@param  Correl  Covariance matrix
     *@param  Infin   as above
     *@param  Maxpts  maximum number of itteration
     *@param  Abseps  absolute error
     *@param  Releps  relative arror
     */
    public void Mvndstpack(double[] Lower, double[] Upper, int[] Infin, double[] Correl,
            int Maxpts, double Abseps, double Releps) {
        n = Lower.length;
        this.lower = Lower;
        this.upper = Upper;
        this.correl = Correl;
        this.maxpts = Maxpts;
        this.abseps = Abseps;
        this.releps = Releps;
        nl = n;
        if (a == null || a.length != nl) {
            cov = new double[nl * (nl + 1) / 2];
            a = new double[nl];
            b = new double[nl];
            infi = new int[nl];
        }
        if (Infin != null) {
            this.infin = Infin;
        }
        else {
            if (infin == null || this.infin.length != n) {
                this.infin = new int[n];
            }
            for (int i = 0; i < n; i++) {
                if (lower[i] == Double.NEGATIVE_INFINITY) {
                    if (upper[i] == Double.POSITIVE_INFINITY) {
                        infin[i] = -1;
                    }
                    else {
                        infin[i] = 0;
                    }
                }
                else if (upper[i] == Double.POSITIVE_INFINITY) {
                    infin[i] = 1;
                }
                else {
                    infin[i] = 2;
                }
            }
        }
        if (upper.length != n || infin.length != n || correl.length != n * (n
                 - 1) / 2) {
            throw new IllegalArgumentException("one of the arguments is illigal");
        }
        this.correl = Correl;
        this.maxpts = Maxpts;
        this.abseps = Abseps;
        this.releps = Releps;
        compute();
    }


    /**
     *  To actualy re compute.
     */
    public void run() {
        compute();
    }


    /**
     *  Same as run
     *
     *@see    run()
     */

    private void compute() {
        running = true;
        if (n > 100 || n < 1) {
            inform = 1;
            value = 0;
            error = 1;
        }
        else {
            inform = (int) mvndnt(infis, d, e);
            if (n - infis[0] == 0) {
                value = 1;
                error = 1;
            }
            else if (n - infis[0] == 1) {
                value = e[0] - d[0];
                error = 2e-16;
            }
            else {
                dkbvrc(n - infis[0] - 1, 0, maxpts);
            }
        }
        running = false;
    }


    /**
     *  Discription from FORTRAN code: <br>
     *  Automatic Multidimensional Integration Subroutine AUTHOR: Alan Genz
     *  Department of Mathematics Washington State University Pulman, WA
     *  99164-3113 Email: AlanGenz@wsu.edu Last Change: 5/15/98 KRBVRC computes
     *  an approximation to the integral 1 1 1 I I ... I F(X)
     *  dx(NDIM)...dx(2)dx(1) 0 0 0 DKBVRC uses randomized Korobov rules for the
     *  first 20 variables. The primary references are "Randomization of Number
     *  Theoretic Methods for Multiple Integration" R. Cranley and T.N.L.
     *  Patterson, SIAM J Numer Anal, 13, pp. 904-14, and "Optimal Parameters
     *  for Multidimensional Integration", P. Keast, SIAM J Numer Anal, 10,
     *  pp.831-838. If there are more than 20 variables, the remaining variables
     *  are integrated using Richtmeyer rules. A reference is "Methods of
     *  Numerical Integration", P.J. Davis and P. Rabinowitz, Academic Press,
     *  1984, pp. 482-483. Parameters Input parameters NDIM Number of variables,
     *  must exceed 1, but not exceed 40 MINVLS Integer minimum number of
     *  function evaluations allowed. MINVLS must not exceed MAXVLS. If MINVLS <
     *  0 then the routine assumes a previous call has been made with the same
     *  integrand and continues that calculation. MAXVLS Integer maximum number
     *  of function evaluations allowed. FUNCTN EXTERNALly declared user defined
     *  function to be integrated. It must have parameters (NDIM,Z), where Z is
     *  a real array of dimension NDIM. ABSEPS Required absolute accuracy.
     *  RELEPS Required relative accuracy. Output parameters MINVLS Actual
     *  number of function evaluations used. ABSERR Estimated absolute accuracy
     *  of FINEST. FINEST Estimated value of integral. INFORM INFORM = 0 for
     *  normal exit, when ABSERR <= MAX(ABSEPS, RELEPS*ABS(FINEST)) and INTVLS
     *  <= MAXCLS. INFORM = 1 If MAXVLS was too small to obtain the required
     *  accuracy. In this case a value FINEST is returned with estimated
     *  absolute accuracy ABSERR. Creation date: (10/11/2000 2:43:21 AM)
     *
     *@param  ndim    int
     *@param  minvls  int[]
     *@param  maxvls  int
     */
    private void dkbvrc(int ndim, int minvls, int maxvls) {
        final int plim = 25;
        final int nlim = 100;
        final int klim = 20;
        final int minsmp = 8;
        int klimi;
        int intvls;
        int i = 0;
        double difint;
        double finval;
        double varsqr;
        double varprd;
        double finest = value;
        double[] x = new double[2 * nlim];
        double[] vk = new double[klim];
        inform = 1;
        intvls = 0;
        klimi = klim;
        int label = 0;
        main :
        do {
            switch (label) {
                case 0:
                    if (minvls >= 0) {
                        finest = 0;
                        varest = 0;
                        sampls = minsmp;
                        for (i = 0; i < plim; i++) {
                            np = i;
                            if (minvls < 2 * sampls * p[i]) {
                                label = 10;
                                continue main;
                            }
                        }
                        sampls = Math.max(minsmp, minvls / (2 * p[np]));
                    }
                case 10:
                    vk[0] = 1.0 / p[np];
                    for (i = 1; i < Math.min(ndim, klim); i++) {
                        vk[i] = MyMath.mod(c[np][Math.min(ndim - 1, klim -
                                1) - 1] * vk[i - 1], 1.0);
                    }
                    finval = 0;
                    varsqr = 0;
                    for (i = 0; i < sampls; i++) {
                        dksmrc(ndim, klimi, p[np], vk, x);
                        difint = (value - finval) / (i + 1.0);
                        finval += difint;
                        varsqr = (i - 1.0) * varsqr / (i + 1.0) + difint * difint;
                    }
                    intvls += 2 * sampls * p[np];
                    varprd = varest * varsqr;
                    finest += (finval - finest) / (1 + varprd);
                    if (varsqr > 0) {
                        varest = (1 + varprd) / varsqr;
                    }
                    error = 3 * Math.sqrt(varsqr / (1 + varprd));
                    if (error > Math.max(abseps, Math.abs(finest) * releps)) {
                        if (np < plim) {
                            np++;
                        }
                        else {
                            sampls = Math.min(3 * sampls / 2, (maxvls - intvls) / (2 * p[np]));
                            sampls = Math.max(minsmp, sampls);
                        }
                        if (intvls + 2 * sampls * p[np] <= maxvls) {
                            label = 10;
                            continue main;
                        }
                    }
                    else {
                        inform = 0;
                    }
            }
            minvls = intvls;
            break;
        } while (true);
        value = finest;
    }


    /**
     *  Used internaly, see original FORTRAN code for details
     */
    private void mvndst() {
        int n = lower.length;
        int[] infis = new int[1];
        double[] e = new double[1];
        double[] d = new double[1];
        if (n > 100 || n < 1) {
            inform = 1;
            value = 0;
            error = 1;
        }
        else {
            inform = (int) mvndnt(infis, d, e);
            if (n - infis[0] == 0) {
                value = 1;
                error = 1;
            }
            else if (n - infis[0] == 1) {
                value = e[0] - d[0];
                error = 2e-16;
            }
            else {
                dkbvrc(n - infis[0] - 1, 0, maxpts);
            }
        }
    }


//	/**
//	 *  Insert the method's description here. Creation date: (11/23/2000 5:20:34
//	 *  PM)
//	 *
//	 *@param  Lower   double[]
//	 *@param  Upper   double[]
//	 *@param  Infin   int[]
//	 *@param  Cor     double[][]
//	 *@param  Maxpts  int
//	 *@param  Abseps  double
//	 *@param  Releps  double
//	 */
//	public void recompute(double[] Lower, double[] Upper, int[] Infin, double[][] Cor,
//			int Maxpts, double Abseps, double Releps) {
//		recompute(Lower, Upper, Infin, cor2to1(Cor), Maxpts, Abseps, Releps);
//	}
//
//
//	/**
//	 *  Insert the method's description here. Creation date: (11/23/2000 5:20:34
//	 *  PM)
//	 *
//	 *@param  Lower   double[]
//	 *@param  Upper   double[]
//	 *@param  Infin   int[]
//	 *@param  Correl  double[]
//	 *@param  Maxpts  int
//	 *@param  Abseps  double
//	 *@param  Releps  double
//	 */
//	public void recompute(double[] Lower, double[] Upper, int[] Infin, double[] Correl,
//			int Maxpts, double Abseps, double Releps) {
//		if (n != Lower.length || upper.length != n || infin.length != n ||
//				correl.length != n * (n + 1) / 2) {
//			throw new IllegalArgumentException("one of the arguments is illigal");
//		}
//		this.lower = Lower;
//		this.upper = Upper;
//		this.correl = Correl;
//		this.maxpts = Maxpts;
//		this.abseps = Abseps;
//		this.releps = Releps;
//		nl = n;
//		for (int i = 0; i < cov.length; i++) {
//			cov[i] = 0;
//		}
//		for (int i = 0; i < nl; i++) {
//			a[i] = 0;
//			b[i] = 0;
//			infi[i] = 0;
//		}
//		//cov = new double[nl * (nl + 1) / 2];
//		//a = new double[nl];
//		//b = new double[nl];
//		//infi = new int[nl];
//		if (Infin != null) {
//			this.infin = Infin;
//		}
//		else {
//			this.infin = new int[n];
//			for (int i = 0; i < n; i++) {
//				if (lower[i] == Double.NEGATIVE_INFINITY) {
//					if (upper[i] == Double.POSITIVE_INFINITY) {
//						infin[i] = -1;
//					}
//					else {
//						infin[i] = 0;
//					}
//				}
//				else if (upper[i] == Double.POSITIVE_INFINITY) {
//					infin[i] = 1;
//				}
//				else {
//					infin[i] = 2;
//				}
//			}
//		}
//		compute();
//	}


    /**
     *  Used internaly, see original FORTRAN code for details
     *
     *@param  y      double[]
     *@param  infis  int[]
     */
    private void covsrt(double[] y, int[] infis) {
        int il;
        int jmin;
        double sumsq;
        double aj;
        double bj;
        double sum;
        double[] d = new double[1];
        double[] e = new double[1];
        double cvdiag;
        double amin = 0;
        double bmin = 0;
        double dmin;
        double emin;
        double yl;
        double yu;
        final double sqtwpi = 2.506628274631001;
        final double eps = 1e-10;
        int ij = -1;
        int ii = -1;
        infis[0] = 0;
        for (int i = 0; i < n; i++) {
            a[i] = 0;
            b[i] = 0;
            infi[i] = infin[i];
            if (infi[i] < 0) {
                infis[0]++;
            }
            else {
                if (infi[i] != 0) {
                    a[i] = lower[i];
                }
                if (infi[i] != 1) {
                    b[i] = upper[i];
                }
            }
            for (int j = 0; j <= i - 1; j++) {
                cov[++ij] = correl[++ii];
            }
            cov[++ij] = 1;
        }
        // First move any doubly infinite limits to innermost positions.
        if (infis[0] < n) {
            main10 :
            for (int i = n - 1; i >= n - infis[0]; i--) {
                if (infi[i] >= 0) {
                    for (int j = 0; j < i - 1; j++) {
                        if (infi[j] < 0) {
                            rcswp(j, i);
                            break main10;
                        }
                    }
                }
            }
        }
        //Sort remaining limits and determine Cholesky factor.
        ii = -1;
        for (int i = 0; i < n - infis[0]; i++) {
            //Determine the integration limits for variable with minimum
            //  expected probability and interchange that variable with Ith.
            dmin = 0;
            emin = 1;
            jmin = i;
            cvdiag = 0;
            ij = ii;
            for (int j = i; j < n - infis[0]; j++) {
                if (cov[ij + j + 1] > eps) {
                    sumsq = Math.sqrt(cov[ij + j + 1]);
                    sum = 0;
                    for (int k = 0; k < i; k++) {
                        sum += cov[ij + k + 1] * y[k];
                    }
                    aj = (a[j] - sum) / sumsq;
                    bj = (b[j] - sum) / sumsq;
                    mvnlms(aj, bj, infi[j], d, e);
                    if (emin + d[0] >= e[0] + dmin) {
                        jmin = j;
                        amin = aj;
                        bmin = bj;
                        dmin = d[0];
                        emin = e[0];
                        cvdiag = sumsq;
                    }
                }
                ij += j + 1;
            }
            if (jmin > i) {
                rcswp(i, jmin);
            }
            cov[ii + i + 1] = cvdiag;
            //Compute Ith column of Cholesky factor.
            //  Compute expected value for Ith integration variable and
            //  scale Ith covariance matrix row and limits.
            if (cvdiag > 0) {
                il = ii + i + 1;
                for (int l = i + 1; l < n - infis[0]; l++) {
                    cov[il + i + 1] /= cvdiag;
                    ij = ii + i + 1;
                    for (int j = i + 1; j <= l; j++) {
                        cov[il + j + 1] -= cov[il + i + 1] * cov[ij + i + 1];
                        ij += j + 1;
                    }
                    il += l + 1;
                }
                if (emin > dmin + eps) {
                    yl = 0;
                    yu = 0;
                    if (infi[i] != 0) {
                        yl = -Math.exp(-amin * amin / 2) / sqtwpi;
                    }
                    if (infi[i] != 1) {
                        yu = -Math.exp(-bmin * bmin / 2) / sqtwpi;
                    }
                    y[i] = (yu - yl) / (emin - dmin);
                }
                else {
                    switch (infi[i]) {
                        case 0:
                            y[i] = bmin;
                            break;
                        case 1:
                            y[i] = amin;
                            break;
                        case 2:
                            y[i] = (amin + bmin) / 2;
                            break;
                    }
                }
                for (int j = 0; j <= i; j++) {
                    cov[++ii] /= cvdiag;
                }
                a[i] /= cvdiag;
                b[i] /= cvdiag;
            }
            else {
                il = ii + i;
                for (int l = i + 1; l < n - infis[0]; l++) {
                    cov[il + i + 1] = 0;
                    il += l + 1;
                }
                //If the covariance matrix diagonal entry is zero,
                //permute limits and/or rows, if necessary.
                main20 :
                for (int j = i - 1; j >= 0; j--) {
                    if (Math.abs(cov[ii + j + 1]) > eps) {
                        a[i] /= cov[ii + j + 1];
                        b[i] /= cov[ii + j + 1];
                        if (cov[ii + j + 1] < 0) {
                            double tmp = a[i];
                            a[i] = b[i];
                            b[i] = tmp;
                            if (infi[i] != 2) {
                                infi[1] = 1 - infi[1];
                            }
                        }
                        for (int l = 0; l <= j; l++) {
                            cov[ii + l + 1] /= cov[ii + j + 1];
                        }
                        for (int l = j + 1; l < i; l++) {
                            if (cov[(l + 1) * l / 2 + j] > 0) {
                                ij = ii;
                                for (int k = i - 1; k >= l; k--) {
                                    for (int m = 0; m <= k; m++) {
                                        double tmp = cov[ij - k + m];
                                        cov[ij - k + m + 1] = cov[ij + m +
                                                1];
                                        cov[ij + m + 1] = tmp;
                                    }
                                    dkswap(a, k, k + 1);
                                    dkswap(b, k, k + 1);
                                    dkswap(infi, k, k + 1);
                                    ij -= k + 1;
                                }
                                break main20;
                            }
                        }
                        break main20;
                    }
                    cov[ii + j + 1] = 0;
                }
                ii += i + 1;
                y[i] = 0;
            }
        }
        double tmp1 = 0;
    }


    /**
     *  This subroutine generates a new quasi-random Richtmeyer vector. A
     *  reference is "Methods of Numerical Integration", P.J. Davis and P.
     *  Rabinowitz, Academic Press, 1984, pp. 482-483. INPUTS: S - the number of
     *  dimensions; KRRCHT is initialized for each new S or S < 1. OUTPUTS:
     *  QUASI - a new quasi-random S-vector Creation date: (10/11/2000 11:20:54
     *  AM)
     *
     *@param  s      double
     *@param  quasi  double[]
     */
    private void dkrcht(int s, double[] quasi) {
        dkrcht(s, quasi, 0);
    }


    /**
     *  This subroutine generates a new quasi-random Richtmeyer vector. A
     *  reference is "Methods of Numerical Integration", P.J. Davis and P.
     *  Rabinowitz, Academic Press, 1984, pp. 482-483. INPUTS: S - the number of
     *  dimensions; KRRCHT is initialized for each new S or S < 1. OUTPUTS:
     *  QUASI - a new quasi-random S-vector Creation date: (10/11/2000 11:20:54
     *  AM)
     *
     *@param  s      double
     *@param  quasi  double[]
     *@param  start  Description of Parameter
     */
    private void dkrcht(int s, double[] quasi, int start) {
        final int mxdim = 80;
        final int mxhsum = 48;
        final int b = 2;
        if (psqt == null) {
            psqt = new double[mxdim];
        }
        if (N == null) {
            N = new int[mxhsum + 1];
        }
        int i;
        double rn;
        if (s != olds || s < 1) {
            olds = s;
            N[0] = 0;
            hisum = 0;
            for (i = 0; i < s; i++) {
                rn = prime[i];
                psqt[i] = Math.sqrt(rn);
            }
        }
        main :
        do {
            for (i = 0; i <= hisum; i++) {
                N[i]++;
                if (N[i] < b) {
                    break main;
                }
                N[i] = 0;
            }
            hisum++;
            if (hisum > 48) {
                hisum = 0;
            }
            N[hisum] = 1;
            break;
        } while (true);
        //10
        rn = 0;
        for (i = hisum; i >= 0; i--) {
            rn = N[i] + b * rn;
        }
        for (i = 0; i < s; i++) {
            quasi[i + start] = MyMath.mod(rn * psqt[i], 1);
        }
    }


    /**
     *  Used internaly, see original FORTRAN code for details
     *
     *@param  ndim   int
     *@param  klim   int
     *@param  prime  int
     *@param  vk     double[]
     *@param  x      double[]
     */
    private void dksmrc(int ndim, int klim, int prime, double[] vk, double[] x) {
        int k;
        int j;
        int jp;
        int nk;
        double xt;
        double sumkpo = value;
        nk = Math.min(ndim, klim);
        for (j = 0; j < nk; j++) {
            jp = j + (int) rnd.raw() * (nk - j);
            xt = vk[j];
            vk[j] = vk[jp];
            vk[jp] = xt;
        }
        for (j = 0; j < ndim; j++) {
            x[ndim + j] = rnd.raw();
        }
        for (k = 1; k <= prime; k++) {
            for (j = 0; j < nk; j++) {
                x[j] = MyMath.mod(k * vk[j], 1);
            }
            if (ndim > klim) {
                dkrcht(ndim - klim, x, klim);
            }
            for (j = 0; j < ndim; j++) {
                xt = x[j] + x[ndim + j + 1];
                if (xt > 1) {
                    xt--;
                }
                x[j] = Math.abs(2 * xt - 1);
            }
            sumkpo += (functn(ndim, x) - sumkpo) / (2 * k - 1);
            for (j = 0; j < ndim; j++) {
                x[j] = 1 - x[j];
            }
            sumkpo += (functn(ndim, x) - sumkpo) / (2 * k);
        }
        value = sumkpo;
    }


    /**
     *  Swap i and j element od x
     *
     *@param  x  double[]
     *@param  i  int
     *@param  j  int
     */
    private void dkswap(double[] x, int i, int j) {
        double tmp = x[i];
        x[i] = x[j];
        x[j] = tmp;
    }


    /**
     *  Swap i and j element od x
     *
     *@param  x  int[]
     *@param  i  int
     *@param  j  int
     */
    private void dkswap(int[] x, int i, int j) {
        int tmp = x[i];
        x[i] = x[j];
        x[j] = tmp;
    }


    /**
     *  same as mvndfn(int, double[])
     *
     *@param  ndim  int
     *@param  x     double[]
     *@return       double
     *@see          mvndfn(int, double[])
     */
    private double functn(int ndim, double[] x) {
        return mvndfn(ndim, x);
    }


    /**
     *  Insert the method's description here. Creation date: (10/10/2000
     *  10:50:14 PM)
     *
     *@param  n  int
     *@param  w  double[]
     *@return    double
     */
    private double mvndfn(int n, double[] w) {
        //return mvndfn_mvndnt(0,w,  null, null, null);
        double Mvndfn = 0;
        double[] y = new double[nl];
        int ij;
        int jij;
        int ik;
        int infa;
        int infb;
        double sum;
        double ai = 0;
        double bi = 0;
        double[] di = new double[1];
        double[] ei = new double[1];
        Mvndfn = 1;
        infa = 0;
        infb = 0;
        ik = 0;
        ij = -1;
        for (int i = 0; i <= n; i++) {
            sum = 0;
            for (int j = 0; j <= i - 1; j++) {
                ij++;
                if (j < ik) {
                    sum += cov[ij] * y[j];
                }
            }
            if (infi[i] != 0) {
                if (infa == 1) {
                    ai = Math.max(ai, a[i] - sum);
                }
                else {
                    ai = a[i] - sum;
                    infa = 1;
                }
            }
            if (infi[i] != 1) {
                if (infb == 1) {
                    bi = Math.min(bi, b[i] - sum);
                }
                else {
                    bi = b[i] - sum;
                    infb = 1;
                }
            }
            ij++;
            if (i == n || (ij + ik + 2 < cov.length && cov[ij + ik + 2] > 0)) {
                mvnlms(ai, bi, 2 * infa + infb - 1, di, ei);
                if (di[0] >= ei[0]) {
                    return 0;
                }
                else {
                    Mvndfn *= ei[0] - di[0];
                    if (i <= n) {
                        y[ik] = CdfInv.nor(di[0] + w[ik] * (ei[0] - di[0]));
                    }
                    ik++;
                    infa = 0;
                    infb = 0;
                }
            }
        }
        return Mvndfn;
    }


    /**
     *  Used internaly, see original FORTRAN code for details
     *
     *@param  infis  int[]
     *@param  d      double[]
     *@param  e      double[]
     *@return        double[]
     */
    private double mvndnt(int[] infis, double[] d, double[] e) {
        //return mvndfn_mvndnt(1,null,infis,d,e);
        double Mvndnt = 0;
        double[] y = new double[nl];
        covsrt(y, infis);
        if (n - infis[0] == 1) {
            mvnlms(a[0], b[0], infi[0], d, e);
        }
        else if (n - infis[0] == 2) {
            if (Math.abs(cov[2]) > 0) {
                d[0] = Math.sqrt(1 + cov[1] * cov[1]);
                if (infi[1] != 0) {
                    a[1] /= d[0];
                }
                if (infi[1] != 1) {
                    b[1] /= d[0];
                }
                double[] aa = {
                        a[0], a[1]
                        };
                double[] bb = {
                        b[0], b[1]
                        };
                int[] iinfi = {
                        infi[0], infi[1]
                        };
                e[0] = Cdf.bvnor(aa, bb, iinfi, cov[1] / d[0]);
                d[0] = 0;
            }
            else {
                if (infi[0] != 0) {
                    if (infi[1] != 0) {
                        a[0] = Math.max(a[0], a[1]);
                    }
                    else if (infi[1] != 0) {
                        a[0] = a[1];
                    }
                }
                if (infi[1] != 1) {
                    if (infi[1] != 1) {
                        b[0] = Math.min(b[0], b[1]);
                    }
                    else if (infi[1] != 1) {
                        b[0] = b[1];
                    }
                }
                if (infi[0] != infi[1]) {
                    infi[0] = 2;
                    mvnlms(a[0], b[0], infi[0], d, e);
                }
                infis[0]++;
            }
        }
        return Mvndnt;
    }


    /**
     *  Used internaly, see original FORTRAN code for details
     *
     *@param  a      double
     *@param  b      double
     *@param  infin  int
     *@param  lower  double[]
     *@param  upper  double[]
     */
    private void mvnlms(double a, double b, int infin, double[] lower, double[] upper) {
        lower[0] = 0;
        upper[0] = 1;
        if (infin >= 0) {
            if (infin != 0) {
                lower[0] = Cdf.nor(a);
            }
            if (infin != 1) {
                upper[0] = Cdf.nor(b);
            }
        }
    }


    /**
     *  Swaps rows and columns P and Q in situ, with P <= Q. Creation date:
     *  (10/12/2000 12:36:03 AM)
     *
     *@param  p  int
     *@param  q  int
     */
    private void rcswp(int p, int q) {
        dkswap(a, p, q);
        dkswap(b, p, q);
        dkswap(infi, p, q);
        int jj = p * (p + 1) / 2 - 1;
        int ii = q * (q + 1) / 2 - 1;
        dkswap(cov, jj + p + 1, ii + q + 1);
        for (int j = 0; j < p; j++) {
            dkswap(cov, jj + j + 1, ii + j + 1);
        }
        jj += p + 1;
        for (int i = p + 1; i < q; i++) {
            dkswap(cov, jj + p + 1, ii + i + 1);
            jj += i + 1;
        }
        ii += q + 1;
        for (int i = q + 1; i < n; i++) {
            dkswap(cov, ii + p + 1, ii + q + 1);
            ii += i + 1;
        }
    }


    /**
     *  Insert the method's description here. Creation date: (10/10/2000 6:47:53
     *  PM)
     *
     *@param  cor2  Description of Parameter
     */
//	public static void main(String[] args) {
//		double abseps;
//		double releps;
//		double val;
//		double err;
//		int n;
//		int nn;
//		int i;
//		int j;
//		int k;
//		int ij;
//		int ift;
//		n = 5;
//		nn = (n - 1) * n / 2;
//		int maxpts = 5000 * n * n * n;
//		abseps = 0.00005;
//		releps = 0;
//		double[] correl = {
//				-0.707107, 0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5
//				};
//		double[] low = {
//				.0, .0, 1.7817, 1.4755, 1.5949
//				};
//		double[] up = {
//				.0, 1.5198, 1.7817, 1.4755, 1.5949
//				};
//		int[] infin = {
//				1, 2, 1, 1, 0
//				};
//		System.out.println("               Test of MVNDST");
//		System.out.println("            Requested Accuracy " + Math.max(abseps,
//				releps));
//		System.out.println("           Number of Dimensions is " + n);
//		System.out.println("     Maximum # of Function Values is " + maxpts);
//		for (k = 1; k <= 3; k++) {
//			System.out.println();
//			System.out.println("I     Limits");
//			System.out.println("    'Lower  Upper  Lower Left of Correlation Matrix");
//			ij = 0;
//			for (i = 0; i < n; i++) {
//				if (infin[i] < 0) {
//					System.out.print(i + " -infin  infin ");
//					for (j = 0; j < i - 1; j++) {
//						System.out.print(correl[ij + j] + " ");
//					}
//					System.out.println("1.0");
//				}
//				else if (infin[i] == 0) {
//					System.out.print(i + " -infin " + up[i] + " ");
//					for (j = 0; j < i - 1; j++) {
//						System.out.print(correl[ij + j] + " ");
//					}
//					System.out.println("1.0");
//				}
//				else if (infin[i] == 1) {
//					System.out.print(i + " " + low[i] + " infin ");
//					for (j = 0; j < i - 1; j++) {
//						System.out.print(correl[ij + j] + " ");
//					}
//					System.out.println("1.0");
//				}
//				else {
//					System.out.print(i + " " + low[i] + " " + up[i] + " ");
//					for (j = 0; j < i - 1; j++) {
//						System.out.print(correl[ij + j] + " ");
//					}
//					System.out.println("1.0");
//				}
//				ij += i;
//			}
//			Mvndstpack T = new Mvndstpack(low, up, infin, correl, maxpts, abseps,
//					releps);
//			err = T.error;
//			val = T.value;
//			ift = T.inform;
//			//mvndst(low,up,infin,correl,maxpts,abseps,releps,err,val,ift);
//			System.out.println(" Results for:  MVNDST");
//			System.out.println("      Value      :   " + val + " " + ift);
//			System.out.println("  Error Estimate :   (" + err + ")");
//			infin[0]--;
//		}
//	}


    /**
     *  Insert the method's description here. Creation date: (10/10/2000 6:47:53
     *  PM)
     *
     *  Transform 2 dimentional matric in to 1 dimentional
     *
     *@param  cor2  Description of Parameter
     */
    private void cor2to1(double[][] cor2) {
        if (correl == null || correl.length != cor2.length * (cor2.length - 1) / 2) {
            correl = new double[cor2.length * (cor2.length - 1) / 2];
        }
        int ij = 0;
        for (int i = 0; i < cor2.length; i++) {
            for (int j = 0; j < i; j++) {
                correl[ij++] = cor2[i][j];
            }
        }
        //return cor1;
    }
}

