/*
 * put your module comment here
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 */

package com.kutsyy.util;

import java.util.*;
import java.awt.event.*;

//import VisualNumerics.math.Statistics.*;
/**
 *  Cdf class include functions for computing Cumulative density functions or
 *  integrals of PDF (in case of multivariate) <BR>
 *  Created by <A href="http://www.kutsyy.com">Vadim Kutsyy</A> <BR>
 *
 *
 *@author     <A href="http://www.kutsyy.com">Vadim Kutsyy</A>
 *@created    January 13, 2001
 */
public final class Cdf {
    /**
     *  default minimum tolerance for mvnor
     */
    private final static double mvnor_epsmin = 1e-3;
    /**
     *  default tolerance for mvnor
     */
    private final static double mvnor_eps = 1e-3;
    /**
     *  default tolerance for ranmvn
     */
    private final static double ranmvn_eps = 1e-3;
    /**
     *  default relevant tolerance for ranmvn
     */
    private final static double ranmvn_releps = 1e-3;
    /**
     *  variable used by mvnor
     */
    private static double mvnor_ept;
    /**
     *  variable used by mvnor
     */
    private static double mvnor_eplos;
    /**
     *  variable used by mvnor
     */
    private static double mvnor_z = 100;
    /**
     *  variable used internaly
     */
    private static double cup;
    /**
     *  variable used by mvnor
     */
    private static double mvnor_cup;
    /**
     *  constant used by mvnor
     */
    private final static double mvnor_z_cup = 0.6744897501960817;
    /**
     *  default tolerance used by bvnor
     */
    private final static double bvnor_eps = 1e-10;
    /**
     *  Internal random number generator
     */
    private static Random rnd = new Random();
    /**
     *  variable used by ranmvn
     */
    private static double ranmvn_error = 0;
    /**
     *  variable used by ranmvn
     */
    private static double ranmvn_value = 0;
    /**
     *  variable used by ranmvn
     */
    private static double ranmvn_varest = 0;
    /**
     *  variable used by dmv
     */
    private final static double DMV_CUT = 10;


    /**
     *  A function for computing bivariate normal probabilities. <BR>
     *
     *
     *@param  x      upper limit
     *@param  sigma  covariate matrix
     *@return        double result
     *@see           #bvnor(double,double,double)
     */
    public final static double bvnor(double[] x, double[][] sigma) {
        if (!(x.length == 2 && sigma.length == 2 && sigma[0].length == 2)) {
            throw new IllegalArgumentException("X and mu have to be double[2], and sigma has to be double[2][2]");
        }
        if (sigma[0][1] != sigma[1][0]) {
            throw new IllegalArgumentException("sigma must be symetric");
        }
        return bvnor(x[0] / Math.sqrt(sigma[0][0]), x[1] / Math.sqrt(sigma[1][1]),
                sigma[0][1] / Math.sqrt(sigma[0][0] * sigma[1][1]));
    }


    /**
     *  A function for computing bivariate normal probabilities. <BR>
     *
     *
     *@param  lower  lower limits of integrations
     *@param  upper  upper limits of integration
     *@param  sigma  covariance matrix
     *@return        probability
     *@see           #bvnor(double,double,double)
     */
    public final static double bvnor(double[] lower, double[] upper, double[][] sigma) {
        if (lower.length == 1 || upper.length == 1 || sigma.length == 1 ||
                sigma[0].length == 1) {
            return nor(upper[0] / Math.sqrt(sigma[0][0])) - nor(lower[0] / Math.sqrt(sigma[0][0]));
        }
        if (lower.length != 2 || upper.length != 2 || sigma.length != 2 ||
                sigma[0].length != 2) {
            throw new IllegalArgumentException("all matrix are not length of 2");
        }
        double[] c = (double[]) lower.clone();
        double[] d = (double[]) upper.clone();
        double[] sd = {
                Math.sqrt(sigma[0][0]), Math.sqrt(sigma[1][1])
                };
        for (int i = 0; i < 2; i++) {
            c[i] /= sd[i];
            d[i] /= sd[i];
        }
        return bvnor(c, d, sigma[1][0] / (sd[0] * sd[1]));
    }


    /**
     *  A function for computing bivariate normal probabilities. <BR>
     *
     *
     *@param  lower  lower limits of integrations
     *@param  upper  upper limits of integration
     *@param  infin  integration limits, provided in order to work with
     *      Mvndstpack
     *@param  cor    correlation coeffitient
     *@return        probability
     *@return        result
     *@see           #bvnor(double,double,double)
     */
    public final static double bvnor(double[] lower, double[] upper, int[] infin,
            double cor) {
        for (int i = 0; i < 2; i++) {
            if (infin[i] < 1) {
                lower[i] = Double.NEGATIVE_INFINITY;
            }
            if (infin[i] < 0 || infin[i] == 1) {
                upper[i] = Double.POSITIVE_INFINITY;
            }
        }
        double[][] cov = {
                {
                1, cor
                }, {
                cor, 1
                }
                };
        return bvnor(lower, upper, cov);
    }


    /**
     *  A function for computing bivariate normal probabilities. <BR>
     *
     *
     *@param  lower  lower limit of integration
     *@param  upper  upper limit of integration
     *@param  rho    correlation coeffitient
     *@return        probability
     *@author:       <A href="http://www.kutsyy.com>Vadum Kutsyy <\A>
     *@see           #bvnor(double,double,double)
     */
    public final static double bvnor(double[] lower, double[] upper, double rho) {
        if (lower.length == 1 && upper.length == 1) {
            return nor(upper[0]) - nor(lower[0]);
        }
        if (lower.length != 2 || upper.length != 2) {
            throw new IllegalArgumentException("all matrix are not length of 2");
        }
        if (lower[0] == upper[0] || lower[1] == upper[1]) {
            return 0;
        }
        if (lower[0] > upper[0] || lower[1] > upper[1]) {
            throw new IllegalArgumentException("lower limit bigger than upper");
        }
        int[] inf = {
                3, 3
                };
        for (int i = 0; i < 2; i++) {
            if (lower[i] == Double.NEGATIVE_INFINITY) {
                if (upper[i] != Double.POSITIVE_INFINITY) {
                    inf[i] = 0;
                }
            }
            else if (upper[i] == Double.POSITIVE_INFINITY) {
                inf[i] = 1;
            }
            else {
                inf[i] = 2;
            }
        }
        switch (inf[0]) {
            case 0:
                switch (inf[1]) {
                    case 0:
                        return bvnor(upper[0], upper[1], rho);
                    case 1:
                        return bvnor(upper[0], -lower[1], -rho);
                    case 2:
                        return bvnor(upper[0], upper[1], rho) - bvnor(upper[0],
                                lower[1], rho);
                    default:
                        return nor(upper[1]);
                }
            case 1:
                switch (inf[1]) {
                    case 0:
                        return bvnor(-lower[0], upper[1], -rho);
                    case 1:
                        return bvnor(-lower[0], -lower[1], rho);
                    case 2:
                        return bvnor(-lower[0], -lower[1], rho) - bvnor(-lower[0],
                                -upper[1], rho);
                    default:
                        return 1 - nor(lower[1]);
                }
            case 2:
                switch (inf[1]) {
                    case 0:
                        return bvnor(upper[0], upper[1], rho) - bvnor(lower[0],
                                upper[1], rho);
                    case 1:
                        return bvnor(-lower[0], -lower[1], rho) - bvnor(-upper[0],
                                -lower[1], rho);
                    case 2:
                    {
                        double t = bvnor(upper[0], upper[1], rho) + bvnor(lower[0],
                                lower[1], rho) - bvnor(upper[0], lower[1],
                                rho) - bvnor(lower[0], upper[1], rho);
                        if (t < 0 || t > 1) {
                            t = bvnor(-upper[0], -upper[1], rho) + bvnor(-lower[0],
                                    -lower[1], rho) - bvnor(-upper[0],
                                    -lower[1], rho) - bvnor(-lower[0],
                                    -upper[1], rho);
                        }
                        if (t < 0 || t > 1) {
                            throw new IllegalStateException("result of probability is outside the reagen");
                        }
                        return t;
                    }
                    default:
                        return nor(upper[1]) - nor(lower[1]);
                }
            default:
                switch (inf[1]) {
                    case 0:
                        return nor(upper[0]);
                    case 1:
                        return 1 - nor(lower[0]);
                    case 2:
                        return nor(upper[0]) - nor(lower[0]);
                    default:
                        return 1;
                }
        }
    }


    /**
     *  A function for computing bivariate normal probabilities. <BR>
     *  Based on algorithms by <BR>
     *  Yihong Ge, Department of Computer Science and Electrical Engineering,
     *  Washington State University, Pullman, WA 99164-2752 <BR>
     *  and <BR>
     *  Alan Genz <a href="http://www.sci.wsu.edu/math/faculty/genz/homepage">
     *  http://www.sci.wsu.edu/math/faculty/genz/homepage</a> <BR>
     *  Department of Mathematics, Washington State University, Pullman, WA
     *  99164-3113, Email : alangenz@wsu.edu <BR>
     *  <BR>
     *
     *
     *@param  sh  integration limit
     *@param  sk  integration limit
     *@param  r   correlation coefficient
     *@return     result
     */
    public final static double bvnor(double sh, double sk, double r) {
        sh = -sh;
        sk = -sk;
        if (sh == Double.POSITIVE_INFINITY) {
            return nor(sk);
        }
        if (sh == Double.NEGATIVE_INFINITY || sk == Double.NEGATIVE_INFINITY) {
            return 0;
        }
        if (sk == Double.POSITIVE_INFINITY) {
            return nor(sh);
        }
        double[][] w = {
                {
                0.1713244923791705e+00, 0.4717533638651177e-01, 0.1761400713915212e-01
                }, {
                0.3607615730481384e+00, 0.1069393259953183e+00, 0.4060142980038694e-01
                }, {
                0.4679139345726904e+00, 0.1600783285433464e+00, 0.6267204833410906e-01
                }, {
                0, 0.2031674267230659e+00, 0.8327674157670475e-01
                }, {
                0, 0.2334925365383547e+00, 0.1019301198172404e+00
                }, {
                0, 0.2491470458134029e+00, 0.1181945319615184e+00
                }, {
                0, 0, 0.1316886384491766e+00
                }, {
                0, 0, 0.1420961093183821e+00
                }, {
                0, 0, 0.1491729864726037e+00
                }, {
                0, 0, 0.1527533871307259e+00
                }
                };
        double[][] x = {
                {
                -0.9324695142031522e+00, -0.9815606342467191e+00, -0.9931285991850949e+00
                }, {
                -0.6612093864662647e+00, -0.9041172563704750e+00, -0.9639719272779138e+00
                }, {
                -0.2386191860831970e+00, -0.7699026741943050e+00, -0.9122344282513259e+00
                }, {
                0, -0.5873179542866171e+00, -0.8391169718222188e+00
                }, {
                0, -0.3678314989981802e+00, -0.7463319064601508e+00
                }, {
                0, -0.1252334085114692e+00, -0.6360536807265150e+00
                }, {
                0, 0, -0.5108670019508271e+00
                }, {
                0, 0, -0.3737060887154196e+00
                }, {
                0, 0, -0.2277858511416451e+00
                }, {
                0, 0, -0.7652652113349733e-01
                }
                };
        double bvn;
        double as;
        double a;
        double b;
        double c;
        double d;
        double rs;
        double xs;
        double mvnphi;
        double sn;
        double asr;
        double h;
        double k;
        double bs;
        double hs;
        double hk;
        int lg = 10;
        int ng = 2;
        if (Math.abs(r) < 0.3) {
            ng = 0;
            lg = 3;
        }
        else if (Math.abs(r) < 0.75) {
            ng = 1;
            lg = 6;
        }
        ;
        h = sh;
        k = sk;
        hk = h * k;
        bvn = 0;
        if (Math.abs(r) < 0.925) {
            hs = (h * h + k * k) / 2;
            asr = Math.asin(r);
            for (int i = 0; i < lg; i++) {
                sn = Math.sin(asr * (x[i][ng] + 1) / 2.0);
                bvn += w[i][ng] * Math.exp((sn * hk - hs) / (1 - sn * sn));
                sn = Math.sin(asr * (-x[i][ng] + 1) / 2);
                bvn += w[i][ng] * Math.exp((sn * hk - hs) / (1 - sn * sn));
            }
            bvn = bvn * asr / (4 * Math.PI) + nor(-h) * nor(-k);
        }
        else {
            if (r < 0) {
                k = -k;
                hk = -hk;
            }
            if (Math.abs(r) < 1) {
                as = (1 - r) * (1 + r);
                a = Math.sqrt(as);
                bs = (h - k) * (h - k);
                c = (4 - hk) / 8;
                d = (12 - hk) / 16;
                bvn = a * Math.exp(-(bs / as + hk) / 2) * (1 - c * (bs - as) * (1 - d * bs / 5) / 3
                         + c * d * as * as / 5);
                if (hk > -160) {
                    b = Math.sqrt(bs);
                    bvn -= Math.exp(-hk / 2) * Math.sqrt(2 * Math.PI) * nor(-b / a) * b * (
                            1 - c * bs * (1 - d * bs / 5) / 3);
                }
                a /= 2;
                for (int i = 0; i < lg; i++) {
                    xs = a * (x[i][ng] + 1) * a * (x[i][ng] + 1);
                    rs = Math.sqrt(1 - xs);
                    bvn += a * w[i][ng] * (Math.exp(-bs / (2 * xs) - hk / (1 + rs)) / rs
                             - Math.exp(-(bs / xs + hk) / 2) * (1 + c * xs * (1 + d * xs)));
                    xs = as * (-x[i][ng] + 1) * (-x[i][ng] + 1) / 4;
                    rs = Math.sqrt(1 - xs);
                    bvn += a * w[i][ng] * Math.exp(-(bs / xs + hk) / 2) * (Math.exp(-hk * (
                            1 - rs) / (2 * (1 + rs))) / rs - (1 + c * xs * (1 + d * xs)));
                }
                bvn = -bvn / (Math.PI * 2);
            }
            if (r > 0) {
                bvn += nor(-Math.max(h, k));
            }
            if (r < 0) {
                bvn = -bvn + Math.max(0, nor(-h) - nor(-k));
            }
        }
        return bvn;
    }


    /**
     *  Compute Cdf of Miltivariate normal <BR>
     *  Based on algorithms by <BR>
     *  Zvi Drezner <BR>
     *  Calofornia State University, Fullerton <BR>
     *  "Computation of the Multivariate Normal Integral" <BR>
     *  ACM Transaction on Mathematical Software, Vol 18, No 4, December 1992,
     *  Pages 470-480
     *
     *@param  high  upper limits of integration
     *@param  cov   covariate matrix
     *@return       probabilty
     */
    public final static double dmv(double[] high, double[][] cov) {
        if (high.length == 1) {
            return Cdf.nor(high[0] / Math.sqrt(cov[0][0]));
        }
        for (int i = 0; i < high.length; i++) {
            if (high[i] > 0) {
                if (high[i] == Double.POSITIVE_INFINITY) {
                    return dmv(La.removeI(high, i), La.removeI(cov, i));
                }
                else {
                    return dmv(La.removeI(high, i), La.removeI(cov, i)) -
                            dmv(La.negativeI(high, i), La.negativeI(cov, i));
                }
            }
        }
        return dmv_tree(high, cov);
    }


    /**
     *  Compute Multivariate normal integral
     *
     *@param  low   lower limit of integration
     *@param  high  upper limit of integration
     *@param  cov   covariance matrix
     *@return       probability
     *@see          #dmv(double[],double[][])
     */
    public final static double dmv(double[] low, double[] high, double[][] cov) {
        for (int i = 0; i < low.length; i++) {
            if (low[i] != Double.NEGATIVE_INFINITY) {
                return dmv(La.changeI(low, i, Double.NEGATIVE_INFINITY), high,
                        cov) - dmv(La.changeI(low, i, Double.NEGATIVE_INFINITY),
                        La.changeI(high, i, low[i]), cov);
            }
        }
        return dmv(high, cov);
    }


    /**
     *  Compute multivariate normal distribution function and the probability
     *  that a multivariate normal vector falls in a rectangle in n-space <BR>
     *  routine called dmv(double[], double[], double[][]) if the result is less
     *  than (1e-10/low.length), the routine will call ranmvn(double[],
     *  double[], double[]) <BR>
     *  Note: the cut point above seems to work fine, in order to get at most 1%
     *  of relative error
     *
     *@param  low   Low limit of integration
     *@param  high  Upper limits of integration
     *@param  sig   Covariate matrix
     *@return       Probability
     *@see          #dmv(double[], double[][])
     *@see          #ranmvn(double[], double[], double[][])
     */
    public final static double mvnor(double[] low, double[] high, double[][] sig) {
        if (low.length != high.length || high.length != sig.length || sig.length
                 != sig[0].length) {
            throw new IllegalArgumentException("Matrix dimentions must agree");
        }
        if (low.length == 1) {
            return nor(high[0] / Math.sqrt(sig[0][0])) - nor(low[0] / Math.sqrt(sig[0][0]));
        }
        if (low.length == 2) {
            return bvnor(low, high, sig);
        }
        //if (low.length < DMV_CUT) {
        //double t = dmv(low, high, sig);
        //if (t > 1e-10 / low.length)
        //return t;
        //}
        //return ranmvn_main(low, high, sig);
        for (int i = 0; i < low.length; i++) {
            if (low[i] == Double.NEGATIVE_INFINITY && high[i] == Double.POSITIVE_INFINITY) {
                return mvnor(La.removeI(low, i), La.removeI(high, i), La.removeI(sig,
                        i));
            }
        }
        int max = 100 * low.length;
        long t = System.currentTimeMillis();
        Mvndstpack T = new Mvndstpack(low, high, sig, max, ranmvn_error, ranmvn_releps);
        while (T.error > Math.min(ranmvn_eps, Math.abs(T.value) * ranmvn_releps)
                 && T.maxpts < Integer.MAX_VALUE) {
            T.maxpts = T.maxpts < Integer.MAX_VALUE / 10.0 ? T.maxpts * 10 : Integer.MAX_VALUE;
            //System.out.print(T.maxpts + " " + T.error + " " + T.value + " ");
            T.run();
        }
        return T.value;
    }


    /**
     *  Compute Cdf of standart normal
     *
     *@param  x  upper limit
     *@return    probability
     */
    public final static double nor(double x) {
        //return VisualNumerics.math.Statistics.normalCdf(x);
        if (Double.isNaN(x)) {
            return 0;
        }
        if (Double.isInfinite(x)) {
            return x > 0 ? 1 : 0;
        }
        //return (1 + Stat.erf(x / Math.sqrt(2))) / 2;
        double zabs = Math.abs(x);
        if (zabs > 37) {
            return x > 0 ? 1 : 0;
        }
        double expntl = Math.exp(-(zabs * zabs) / 2);
        double p = 0;
        if (zabs < 7.071067811865475) {
            p = expntl * ((((((zabs * .03526249659989109 + .7003830644436881) * zabs
                     + 6.37396220353165) * zabs + 33.912866078383) * zabs + 112.0792914978709) * zabs
                     + 221.2135961699311) * zabs + 220.2068679123761) / (((((((zabs * .08838834764831844
                     + 1.755667163182642) * zabs + 16.06417757920695) * zabs + 86.78073220294608) * zabs
                     + 296.5642487796737) * zabs + 637.3336333788311) * zabs + 793.8265125199484) * zabs
                     + 440.4137358247522);
        }
        else {
            p = expntl / (zabs + 1 / (zabs + 2 / (zabs + 3 / (zabs + 4 / (zabs + .65))))) / 2.506628274631001;
        }
        return x > 0 ? 1 - p : p;
    }


    /**
     *  Compute Cdf of normal
     *
     *@param  x   upper limit
     *@param  mu  mean
     *@param  sd  stand diviation
     *@return     probability
     *@see        #nor(double)
     */
    public final static double nor(double x, double mu, double sd) {
        return nor((x - mu) / sd);
    }


    /**
     *  Compute Cdf of Miltivariate normal <BR>
     *  Based on algorithms by <BR>
     *  Alan Genz <a href="http://www.sci.wsu.edu/math/faculty/genz/homepage">
     *  http://www.sci.wsu.edu/math/faculty/genz/homepage</a> <BR>
     *  Department of Mathematics, Washington State University, Pullman, WA
     *  99164-3113, Email : alangenz@wsu.edu <BR>
     *  "Numerical Computation of Multivariate Normal Probabilities"" <BR>
     *
     *
     *@param  low   lower limit of integration
     *@param  high  upper limits of integration
     *@param  sig   covariate matrix
     *@return       probabilty
     */
    public final static double ranmvn(double[] low, double[] high, double[][] sig) {
        return ranmvn_main(low, high, sig);
    }


    /**
     *  Internal function used by #dmv(double[], double[][])
     *
     *@param  h  double[]
     *@param  r  double[][]
     *@return    double
     *@see       #dmv(double[], double[][])
     */
    private final static double dmv_solve(double[] h, double[][] r) {
        //double[] coef={2.484061520284430E-001,3.923310666523990E-001,2.114181930760570E-001,3.324666035134390E-002,8.248533445156280E-004};
        //double[] x={1.002421519682160E-001,4.828139660462010E-001,1.060949821525720,1.779729418520260,2.669760356087660};
        double[] x = {
                2.16869474675590e-2, 1.12684220347775e-1, 2.70492671421899e-1,
                4.86902370381935e-1, 7.53043683072978e-1, 1.06093100362236,
                1.40425495820363, 1.77864637941183, 2.18170813144494, 2.61306084533352,
                3.07461811380851, 3.57140815113714, 4.11373608977209, 4.72351306243148,
                5.46048893578335
                };
        double[] coef = {
                5.54433663102343e-2, 1.24027738987730e-1, 1.75290943892075e-1,
                1.91488340747342e-1, 1.63473797144070e-1, 1.05937637278492e-1,
                5.00270211534535e-2, 1.64429690052673e-2, 3.57320421428311e-3,
                4.82896509305201e-4, 3.74908650266318e-5, 1.49368411589636e-6,
                2.55270496934465e-8, 1.34217679136316e-10, 9.56227446736465e-14
                };
        int m = h.length;
        double det = La.choldet(r);
        double[][] r1 = La.cholsl(r);
        double[] h1 = new double[m];
        double[] h2 = new double[m];
        double[] y = new double[m];
        int[] list = new int[m];
        double prod = 1;
        for (int i = 0; i < m; i++) {
            h1[i] = Math.sqrt(2 / r1[i][i]);
            h2[i] = h[i] / h1[i];
            list[i] = 0;
            prod *= Math.PI * r1[i][i];
        }
        prod = 1 / Math.sqrt(det * prod);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                r1[i][j] *= h1[i] * h1[j] * 0.5;
            }
        }
        double solve = 0;
        L60 :
        do {
            double sum = 0;
            for (int i = 0; i < m; i++) {
                int l = list[i];
                sum += x[l] * x[l];
                y[i] = x[l] - h2[i];
            }
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    sum -= y[i] * y[j] * r1[i][j];
                }
            }
            sum = Math.exp(sum);
            for (int i = 0; i < m; i++) {
                sum *= coef[list[i]];
            }
            solve += sum;
            for (int i = 0; i < m; i++) {
                list[i]++;
                if (list[i] < x.length) {
                    continue L60;
                }
                list[i] = 0;
            }
            break L60;
        } while (true);
        solve *= prod;
        return solve;
    }


    /**
     *  Internal function used by #dmv(double[], double[][])
     *
     *@param  h  double[]
     *@param  r  double[][]
     *@return    double
     *@see       #dmv(double[], double[][])
     */
    private static double dmv_tree(double[] h, double[][] r) {
        if (h.length < 3) {
            if (h.length == 1) {
                return Cdf.nor(h[0]);
            }
            return Cdf.bvnor(h[0] / Math.sqrt(r[0][0]), h[1] / Math.sqrt(r[1][1]),
                    r[0][1] / Math.sqrt(r[0][0] * r[1][1]));
        }
        for (int i = 0; i < h.length; i++) {
            if (h[i] > 0) {
                return dmv_tree(La.removeI(h, i), La.removeI(r, i)) - dmv_tree(La.negativeI(h,
                        i), La.negativeI(r, i));
            }
        }
        return dmv_solve(h, r);
        ;
    }


    /**
     *  main function for #ranmvn(double[], double[], double[][])
     *
     *@param  sig   double[]
     *@param  low   Description of Parameter
     *@param  high  Description of Parameter
     *@return       double
     *@author:      <Vadum Kutsyy, kutsyy@hotmail.com>
     *@see          #ranmvn(double[], double[], double[][])
     */
    private final static double ranmvn_main(double[] low, double[] high, double[][] sig) {
        if (low.length != high.length || low.length != sig.length || low.length
                 != sig[0].length) {
            throw new IllegalArgumentException("Matrix dimentions must agree");
        }
        if (low.length == 1) {
            return nor(high[0], 0, sig[0][0]) - nor(low[0], 0, sig[0][0]);
        }
        if (low.length == 2) {
            return bvnor(low, high, sig);
        }
        for (int i = 0; i < low.length; i++) {
            if (low[i] >= high[i]) {
                return 0;
            }
        }
        double[] a = (double[]) low.clone();
        double[] b = (double[]) high.clone();
        double[][] v = new double[sig.length][];
        for (int i = 0; i < sig.length; i++) {
            v[i] = (double[]) sig[i].clone();
        }
        //for (int i = 0; i < low.length; i++)
        //if (v[i][i] != 1) {
        //double c = Math.sqrt(v[i][i]);
        //a[i] /= c;
        //b[i] /= c;
        //for (int j = 0; j < a.length; j++) {
        //v[i][j] /= c;
        //v[j][i] /= c;
        //}
        //}
        //for (int i = 0; i < low.length; i++)
        //for (int j = i + 1; j < low.length; j++)
        //if (b[i] - a[i] > b[j] - a[j]) {
        //double aa = a[i];
        //double bb = b[i];
        //a[i] = a[j];
        //b[i] = b[j];
        //a[j] = aa;
        //b[j] = bb;
        //for (int k = 0; k < low.length; k++) {
        //aa = v[i][k];
        //v[i][k] = v[j][k];
        //v[j][k] = aa;
        //}
        //for (int k = 0; k < low.length; k++) {
        //aa = v[k][i];
        //v[k][i] = v[k][j];
        //v[k][j] = aa;
        //}
        //}
        //ranmvn_c =
        v = La.choldc(v);
        return ranmvn_solve(a, b, v, ranmvn_eps);
    }


    /**
     *  Internal function for #ranmvn(double[], double[], double[][])
     *
     *@param  e0     double
     *@param  d0     double
     *@param  lower  double[]
     *@param  upper  double[]
     *@param  cov    double[][]
     *@param  w      Description of Parameter
     *@return        double
     */
    private final static double ranmvn_mvnfnc(double[] w, double e0, double d0,
            double[] lower, double[] upper, double[][] cov) {
        double d = d0;
        double e = e0;
        double prod = e0 - d0;
        int n = w.length;
        double[] y = new double[n];
        for (int i = 1; i <= n; i++) {
            y[i - 1] = CdfInv.nor(d + w[i - 1] * (e - d));
            double sum = 0;
            for (int j = 0; j < i; j++) {
                sum += cov[i][j] * y[j];
            }
            d = (lower[i] == Double.NEGATIVE_INFINITY) ? 0 : Cdf.nor((lower[i]
                     - sum) / cov[i][i]);
            e = (upper[i] == Double.POSITIVE_INFINITY) ? 1 : Cdf.nor((upper[i]
                     - sum) / cov[i][i]);
            prod *= (e - d);
        }
        if (Double.isNaN(prod)) {
            throw new IllegalStateException("prod is NAN");
        }
        y = null;
        return prod;
    }


    /**
     *  Internal function for #ranmvn(double[], double[], double[][])
     *
     *@param  ndim    int
     *@param  maxpts  int
     *@param  ir      int
     *@param  e0      double
     *@param  d0      double
     *@param  lower   double[]
     *@param  upper   double[]
     *@param  cov     double[][]
     *@see            #ranmvn(double[], double[], double[][])
     */
    private final static void ranmvn_rcrude(int ndim, int maxpts, int ir,
            double e0, double d0, double[] lower, double[] upper, double[][] cov) {
        double absest = ranmvn_error;
        double fun;
        double uni;
        double varprd;
        double findif;
        double finval;
        double[] x = new double[ndim];
        if (ir <= 0) {
            ranmvn_varest = 0;
            ranmvn_value = 0;
        }
        double varest = ranmvn_varest;
        double finest = ranmvn_value;
        finval = 0;
        double varsqr = 0;
        int npts = (int) maxpts / 2;
        for (int m = 1; m <= npts; m++) {
            x = Rnd.unif(ndim);
            fun = ranmvn_mvnfnc(x, e0, d0, lower, upper, cov);
            for (int k = 0; k < ndim; k++) {
                x[k] = 1 - x[k];
            }
            fun += ranmvn_mvnfnc(x, e0, d0, lower, upper, cov);
            fun /= 2.0;
            findif = (fun - finval) / (double) m;
            varsqr = (m - 2) * varsqr / m + findif * findif;
            finval += findif;
        }
        varprd = varest * varsqr;
        ranmvn_value = finest + (finval - finest) / (1 + varprd);
        if (varsqr > 0) {
            varest = (1 + varprd) / varsqr;
        }
        ranmvn_varest = varest;
        ranmvn_error = 3 * Math.sqrt(varsqr / (1 + varprd));
    }


    /**
     *  #ranmvn(double[], double[], double[][])
     *
     *@param  lower
     *@param  upper
     *@param  c
     *@param  releps
     *@return         result
     *@see            #ranmvn(double[], double[], double[][])
     */
    private final static double ranmvn_solve(double[] lower, double[] upper,
            double[][] c, double releps) {
        int n = lower.length;
        double d0 = (lower[0] == Double.NEGATIVE_INFINITY) ? 0 : Cdf.nor(lower[0] / c[0][0]);
        double e0 = (upper[0] == Double.POSITIVE_INFINITY) ? 1 : Cdf.nor(upper[0] / c[0][0]);
        int mpt = 250 + 10 * n;
        int ivls = mpt;
        ranmvn_rcrude(n - 1, mpt, 0, e0, d0, lower, upper, c);
        double error = ranmvn_error;
        double eps = Math.min(ranmvn_eps, releps * Math.abs(ranmvn_value));
        int maxItr = 100000 * n;
        while (ranmvn_error > eps && ivls < maxItr) {
            ranmvn_rcrude(n - 1, 10 * n, 1, e0, d0, lower, upper, c);
            eps = Math.min(ranmvn_eps, releps * Math.abs(ranmvn_value));
            ivls += mpt;
        }
        //System.out.println(ivls+" "+ranmvn_value);
        if (Double.isNaN(ranmvn_value)) {
            throw new IllegalMonitorStateException("NAN");
        }
        lower = null;
        upper = null;
        c = null;
        return ranmvn_value;
    }
}

