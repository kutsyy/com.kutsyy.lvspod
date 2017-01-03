/*
 * Exp class contains functions for expected value calculations
 */

package com.kutsyy.util;

/**
 *  Class Esp contains functions for computing expected value for some
 *  distributions. <BR>
 *  Created by <A href="http://www.kutsyy.com">Vadim Kutsyy</A> <BR>
 *
 *
 *@author     <A href="http://www.kutsyy.com">Vadim Kutsyy</A>
 *@created    December 2, 2000
 */
public final class Exp {
    /**
     *  constant used internaly by mvntrc
     */
    private final static double mvntrc_z = 100;
    /**
     *  variable used internaly by mvntrc
     */
    private static double mvntrc_cup;
    /**
     *  constant used internaly by mvntrc
     */
    private final static double mvntrc_eps = 1e-20;
    /**
     *  variable used internaly by mvntrc
     */

    private static double mvntrc_ept;


    /**
     *  Constructor for the mvntrc object
     *
     *@param  low                           lower limit
     *@param  high                          upper limit
     *@param  u                             mean of normal distribution
     *@param  v                             variance of normal distribution
     *@param  eps                           tolerance
     *@param  tmean                         mean of truncated normal (output)
     *@param  tvar                          covariance matrix of truncated
     *      normal (output)
     *@exception  IllegalArgumentException  Description of Exception
     *@see                                  #mvntrcMoments(double[], double[],
     *      double[][], double, double[], double[][])
     */
    public void mvntrcMoments(double[] low, double[] high, double[] u, double[][] v,
            double eps, double[] tmean, double[][] tvar) throws IllegalArgumentException {
        if (low.length != high.length || high.length != u.length || u.length
                 != v.length || v.length != v[0].length) {
            throw new IllegalArgumentException("dimentions must agree");
        }
        for (int i = 0; i < low.length; i++) {
            if (low[i] != Double.NEGATIVE_INFINITY) {
                double[] tmean1 = new double[low.length];
                double[][] tvar1 = new double[low.length][low.length];
                mvntrcMoments(La.changeI(low, i, Double.NEGATIVE_INFINITY),
                        high, u, v, eps, tmean, tvar);
                mvntrcMoments(La.changeI(low, i, Double.NEGATIVE_INFINITY),
                        La.changeI(high, i, low[i]), u, v, eps, tmean1, tvar1);
                for (int i1 = 0; i1 < low.length; i1++) {
                    tmean[i1] -= tmean1[i1];
                    for (int j = 0; j < low.length; j++) {
                        tvar[i][j] -= tvar1[i][j];
                    }
                }
                return;
            }
        }
        mvntrcMoments(high, u, v, eps, tmean, tvar);
    }


    /**
     *  Description of the Method
     *
     *@param  bounds                        lower limit
     *@param  u                             mean of original normal distribution
     *@param  v                             covariance matrix of original normal
     *      distribution
     *@param  eps                           tolerance
     *@param  tmean                         mean of truncated normal (output)
     *@param  tvar                          covariance matrix of truncated
     *      normal (output)
     *@exception  IllegalArgumentException  Description of Exception
     */
    public void mvntrcMoments(double[] bounds, double[] u, double[][] v, double eps,
            double[] tmean, double[][] tvar) throws IllegalArgumentException {
        if (bounds.length != u.length || u.length != v.length || v.length !=
                v[0].length) {
            throw new IllegalArgumentException("dimentions must agree");
        }
        int n = bounds.length;
        if (tmean == null || tmean.length != n) {
            tmean = new double[n];
        }
        if (tvar == null || tvar.length != n || tvar[0].length != n) {
            tvar = new double[n][n];
        }
        final double twopi = 6.283185308;
        final double rtwopi = 2.506628275;
        final double alphas = 1.e-3;
        int[] ind = new int[n];
        double[][] r = new double[n][n];
        double[][] rinv = new double[n][n];
        double[][] rq = new double[n - 1][n - 1];
        double[][] rqr = new double[n - 2][n - 2];
        double[][] phib = new double[n][n];
        double[][] phin2 = new double[n][n];
        double[][] ephin2 = new double[n][n];
        double[][] ww = new double[n][n];
        int[] iord = new int[n];
        double[] sd = new double[n];
        double[] anew = new double[n];
        double[] a = new double[n];
        double[] aq = new double[n - 1];
        double[] aqr = new double[n - 2];
        double[] phiu = new double[n];
        double[] phin1 = new double[n];
        double[] ephin1 = new double[n];
        double w;
        double alpha;
        double sum;
        double aa;
        double error;
        int label = 0;
        int l;
        int k;
        Mvndstpack mulnor;
        //Calculate normalised truncation points and correlation
        //matrix R
        for (int i = 0; i < n; i++) {
            iord[i] = i;
            if (v[i][i] < 0) {
                throw new IllegalArgumentException("variance can not be negative");
            }
            sd[i] = Math.sqrt(v[i][i]);
        }
        // 1
        for (int i = 0; i < n; i++) {
            a[i] = (bounds[i] - u[i]) / sd[i];
            for (int j = 0; j < n; j++) {
                r[i][j] = v[i][j] / (sd[i] * sd[j]);
            }
        }
        //Order variates by decreasing size of range of integration
        //to increase the speed of MULNOR at highest level
        //(Smallest standardised boundary = largest range, etc)
        if (n > 2) {
            for (int i = 0; i < n - 1; i++) {
                for (int j = i + 1; j < n; j++) {
                    if (a[i] > a[j]) {
                        aa = a[i];
                        a[i] = a[j];
                        a[j] = aa;
                        k = iord[i];
                        iord[i] = iord[j];
                        iord[j] = k;
                        for (k = 0; k < n; k++) {
                            aa = r[i][k];
                            r[i][k] = r[j][k];
                            r[j][k] = aa;
                        }
                        for (k = 0; k < n; k++) {
                            aa = r[k][i];
                            r[k][i] = r[k][j];
                            r[k][j] = aa;
                        }
                    }
                }
            }
        }
        //Evaluate n dimensional integral alpha = Phi(a:R) and error bound
        //alpha(1)=alpha+error   alpha(2)=alpha   alpha(3)=alpha-error
        mulnor = new Mvndstpack(a, a, ind, r, Integer.MAX_VALUE, eps, eps);
        alpha = mulnor.value;
        error = mulnor.error;
        //Evaluate univariate and bivariate normal densities
        //Initialise multinormals of dimension n-1 and n-2
        for (int i = 0; i < n; i++) {
            phiu[i] = Math.exp(-a[i] * a[i] / 2.0) / rtwopi;
            phin1[i] = 1;
            ephin1[i] = 0;
            for (int j = i + 1; j < n; j++) {
                phib[i][j] = (a[i] * a[i] - 2.0 * r[i][j] * a[i] * a[j] + a[j] * a[j]) / +(2.0 * (
                        1 - r[i][j] * r[i][j]));
                phib[i][j] = Math.exp(-phib[i][j]) / (twopi * Math.sqrt(1 - r[i][j] * r[i][j]));
                phib[j][i] = phib[i][j];
                phin2[i][j] = 1;
                phin2[j][i] = 1;
                ephin2[i][j] = 0;
                ephin2[j][i] = 0;
            }
        }
        //Calculation of n-1 dimensional integrals Phi(Aq:Rq)
        if (n > 1) {
            ww = new double[n - 1][n - 1];
            ind = new int[n - 1];
            rinv = La.cholsl(r);
            for (int iq = 0; iq < n; iq++) {
                int m1 = 0;
                //Determine bounds of integration,vector Aq
                for (int i = 0; i < n; i++) {
                    if (i != iq) {
                        m1++;
                        aq[m1] = (a[i] - r[i][iq] * a[iq]) / Math.sqrt(1 - r[i][iq] * r[i][iq]);
                        int m2 = 0;
                        for (int j = 0; j < n; j++) {
                            if (j != iq) {
                                m2++;
                                ww[m1][m2] = rinv[i][j];
                            }
                        }
                    }
                }
                //Determine n-1 x n-1 correlation matrix Rq
                ww = La.cholsl(La.removeI(ww, n - 1));
                for (int i = 0; i < m1; i++) {
                    for (int j = 0; j < m1; j++) {
                        rq[i][j] = ww[i][j] / Math.sqrt(ww[i][i] * ww[j][j]);
                    }
                }
                //Evaluate n-1 dimensional integral  Phi(Aq:Rq) and error bound
                //call mulnor(w,aq,w,eps,n-1,ind,phin1(iq),ephin1(iq),i)
                mulnor = new Mvndstpack(aq, aq, ind, rq, Integer.MAX_VALUE,
                        eps, eps);
                phin1[iq] = mulnor.value;
                ephin1[iq] = mulnor.error;
            }
        }
        //Calculation of n-2 dimensional integrals Phi(Aqr:Rqr)
        if (n > 2) {
            ww = new double[n - 2][n - 1];
            ind = new int[n - 2];
            for (int iq = 0; iq < n; iq++) {
                for (int ir = iq + 1; ir < n; ir++) {
                    //Determine n-2 x n-2 correlation matrix Rqr
                    int m1 = 0;
                    for (int i = 0; i < n; i++) {
                        if (i != iq && i != ir) {
                            m1++;
                            int m2 = 0;
                            for (int j = 0; j < n; j++) {
                                if (j != iq && j != ir) {
                                    m2++;
                                    ww[m1][m2] = rinv[i][j];
                                }
                            }
                        }
                    }
                    ww = La.solve(ww);
                    for (int i = 0; i < m1; i++) {
                        for (int j = 0; j < m1; j++) {
                            rqr[i][j] = ww[i][j] / Math.sqrt(ww[i][i] * ww[j][j]);
                        }
                    }
                    //Determine bounds of integration,vector Aqr
                    int m2 = 0;
                    for (int i = 0; i < n; i++) {
                        if (i != iq && i != ir) {
                            m2++;
                            double bsqr = (r[iq][i] - r[iq][ir] * r[ir][i]) / (
                                    1 - r[iq][ir] * r[iq][ir]);
                            double bsrq = (r[ir][i] - r[iq][ir] * r[iq][i]) / (
                                    1 - r[iq][ir] * r[iq][ir]);
                            double rsrq = (1 - r[i][iq] * r[i][iq]) * (1 - r[ir][iq] * r[ir][iq]);
                            rsrq = (r[i][ir] - r[i][iq] * r[iq][ir]) / Math.sqrt(rsrq);
                            aqr[m2] = a[i] - bsqr * a[iq] - bsrq * a[ir];
                            aqr[m2] = aqr[m2] / Math.sqrt((1 - r[i][iq] * r[i][iq]) * (
                                    1 - rsrq * rsrq));
                        }
                    }
                    //Evaluate n-2 dimensional integral Phi(Aqr:Rqr) and error bound
                    //          call mulnor(w,aqr,w,eps,n-2,ind,phin2(iq,ir),ephin2(iq,ir),i)
                    mulnor = new Mvndstpack(aqr, aqr, ind, rq, Integer.MAX_VALUE,
                            eps, eps);
                    phin2[ir][iq] = phin2[iq][ir] = mulnor.value;
                    ephin2[ir][iq] = ephin2[iq][ir] = mulnor.error;
                }
            }
        }
        //Calculation of E(Xi) ,with upper and lower bounds
        //Tallis(1961) , equation (3).
        for (int i = 0; i < n; i++) {
            sum = 0;
            for (int j = 0; j < n; j++) {
                aa = r[i][j] * phiu[j];
                w = aa * phin1[j];
                sum += w;
            }
            tmean[iord[i]] = sum / alpha;
        }
        // Calculation of E(Xi,Xj), with upper and lower bounds
        //Tallis(1961), equation (4).
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                sum = 0;
                for (k = 0; k < n; k++) {
                    aa = r[k][i] * r[k][j] * a[k] * phiu[k];
                    w = aa * phin1[k];
                    sum += w;
                    for (l = 0; l < n; l++) {
                        if (l != k) {
                            aa = (r[l][j] - r[k][l] * r[k][j]) * r[k][i] * phib[k][l];
                            w = aa * phin2[k][l];
                            sum += w;
                        }
                    }
                }
                tvar[iord[j]][iord[i]] = tvar[iord[i]][iord[j]] = r[i][j] +
                        sum / alpha;
            }
        }
        //Calculate truncated covariance matrix, with upper and
        //lower bounds. Original scale restored.
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                tvar[j][i] = tvar[i][j] = tvar[i][j] * sd[i] * sd[j];
            }
        }
        //Original location and scale restored to truncated mean
        for (int i = 0; i < n; i++) {
            tmean[i] = u[i] + sd[i] * tmean[i];
        }
        // Normal termination of subroutine
    }


    /**
     *  calulate mean and covariance matrix of truncated normal
     *
     *@param  bounds                        lower limit
     *@param  u                             mean of original normal distribution
     *@param  v                             covariance matrix of original normal
     *      distribution
     *@param  eps                           tolerance
     *@param  tmean                         mean of truncated normal (output)
     *@param  tvar                          covariance matrix of truncated
     *      normal (output)
     *@exception  IllegalArgumentException  thows if bounds, u of v are not the
     *      same lenth/dimentions
     *@see                                  #mvntrcMoments(double[], double[],
     *      double[][], double, double[], double[][])
     */
    public void mvntrc(double[] bounds, double[] u, double[][] v, double eps,
            double[] tmean, double[][] tvar) throws IllegalArgumentException {
        if (bounds.length != u.length || u.length != v.length || v.length !=
                v[0].length) {
            throw new IllegalArgumentException("dimentions must agree");
        }
        mvntrcMoments(bounds, u, v, eps, tmean, tvar);
        for (int i = 0; i < bounds.length; i++) {
            for (int j = 0; j < bounds.length; j++) {
                tvar[j][i] = tvar[i][j] = tvar[i][j] - (tmean[i] - u[i]) * (tmean[j]
                         - u[j]);
            }
        }
    }


    /**
     *  Description of the Method
     *
     *@param  low                           lower limit
     *@param  high                          upper limit
     *@param  u                             mean of normal distribution
     *@param  v                             variance of normal distribution
     *@param  eps                           tolerance
     *@param  tmean                         mean of truncated normal (output)
     *@param  tvar                          covariance matrix of truncated
     *      normal (output)
     *@exception  IllegalArgumentException  Description of Exception
     *@see                                  #mvntrcMoments(double[], double[],
     *      double[][], double, double[], double[][])
     */
    public void mvntrc(double[] low, double[] high, double[] u, double[][] v,
            double eps, double[] tmean, double[][] tvar) throws IllegalArgumentException {
        if (low.length != high.length || high.length != u.length || u.length
                 != v.length || v.length != v[0].length) {
            throw new IllegalArgumentException("dimentions must agree");
        }
        mvntrcMoments(low, high, u, v, eps, tmean, tvar);
        for (int i = 0; i < low.length; i++) {
            for (int j = 0; j < low.length; j++) {
                tvar[j][i] = tvar[i][j] = tvar[i][j] - (tmean[i] - u[i]) * (tmean[j]
                         - u[j]);
            }
        }
    }


    /**
     *  transfor one dimentional array in to two dimentional, used internaly
     *
     *@param  r  one dimentional array
     *@return    two dimentional array
     */
    private double[][] cov1to2(double[] r) {
        int n = (1 + (int) Math.sqrt(1 + 8 * r.length)) / 2;
        double[][] cov = new double[n][n];
        int l = 0;
        for (int i = 0; i < n; i++) {
            cov[i][i] = 1;
            for (int j = 0; j < n; j++) {
                cov[i][j] = cov[j][i] = r[l++];
            }
        }
        return cov;
    }


    /**
     *  Evaluate mean of trantated standart normal Creation date: (2/3/00
     *  2:37:20 PM)
     *
     *@param  low   lower trancation
     *@param  high  upper trancation
     *@return       expected value
     */
    public static double ntrc(double low, double high) {
        return ntrc(low, high, 1);
    }


    /**
     *  Evaluate mean of trantated normal
     *
     *@param  low   lower trancation
     *@param  high  upper trancation
     *@param  sd    standart diviation
     *@return       expected value
     */
    public static double ntrc(double low, double high, double sd) {
        return sd * (Pdf.nor(low / sd) - Pdf.nor(high / sd)) / (Cdf.nor(high / sd) -
                Cdf.nor(low / sd));
    }


    //private static double Infinity = Math.exp(1e200);
    /**
     *  Evaluate mean and covariance matrix of multivariate tuncated normal
     *
     *@param  a  lower truncation limit
     *@param  r  covariance matrix or original normal distribution, before
     *      truncation
     *@return    mean of tuncated normal (ie E(X_i))
     *@see       #mvntrcMoments(double[], double[], double[][], double,
     *      double[], double[][])
     */
    public static double[] mvntrc(double[] a, double[][] r) {
        return mvntrcMoment(a, r);
    }


    /**
     *  Evaluate mean and covariance matrix of multivariate tuncated normal
     *
     *@param  a  lower truncation limit
     *@param  b  upper truncation limit
     *@param  r  covariance matrix or original normal distribution, before
     *      truncation
     *@return    mean of tuncated normal (ie E(X_i))
     *@see       #mvntrcMoments(double[], double[], double[][], double,
     *      double[], double[][])
     */
    public static double[] mvntrc(double[] a, double[] b, double[][] r) {
        return mvntrcMoment(a, b, r);
    }


    /**
     *  Evaluate first and second moments of multivariate tuncated normal Based
     *  on algorithms: <BR>
     *  Leppard, P. ; Tallis, G. M. ( <BR>
     *    "Evaluation of the mean and covariance of the truncated multinormal
     *  distribution" <BR>
     *  Applied Statistics, 1989, V 38, 543-553
     *
     *@param  a  lower truncation limit
     *@param  r  covariance matrix or original normal distribution, before
     *      truncation
     *@return    mean of tuncated normal (ie E(X_i))
     *@see       #mvntrcMoments(double[], double[], double[][], double,
     *      double[], double[][])
     */
    public static double[] mvntrcMoment(double[] a, double[][] r) {
        double[] b = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            b[i] = Double.POSITIVE_INFINITY;
        }
        return La.times(mvntrcMomentPrivate(a, r), 1 / Cdf.mvnor(a, b, r));
    }


    /**
     *  Evaluate first and second moments of multivariate tuncated normal
     *
     *@param  a  lower truncation limit
     *@param  b  upper truncation limit
     *@param  r  covariance matrix or original normal distribution, before
     *      truncation
     *@return    mean of tuncated normal (ie E(X_i))
     *@see       #mvntrcMoments(double[], double[], double[][], double,
     *      double[], double[][])
     */
    public static double[] mvntrcMoment(double[] a, double[] b, double[][] r) {
        return La.times(mvntrcMomentPrivate(a, b, r), 1 / Cdf.mvnor(a, b, r));
    }


    /**
     *  transfor two dimentional array in to one dimentional, used internaly
     *
     *@param  cor2  2 dimentional array
     *@return       1 dimentional array
     */
    private final static double[] cor2to1(double[][] cor2) {
        double[] cor1 = new double[cor2.length * (cor2.length + 1) / 2];
        for (int i = 0; i < cor2.length; i++) {
            for (int j = 0; j < i; j++) {
                cor1[j + (i - 1) * i / 2] = cor2[i][j];
            }
        }
        return cor1;
    }


    /**
     *  used internaly
     *
     *@param  low  Description of Parameter
     *@param  R    Description of Parameter
     *@return      double[]
     */
    private static double[] mvntrcMomentPrivate(double[] low, double[][] R) {
        int n = low.length;
        double[] a = new double[n];
        double[][] r = new double[n][n];
        double[] sd = new double[n];
        for (int i = 0; i < n; i++) {
            sd[i] = Math.sqrt(R[i][i]);
        }
        for (int i = 0; i < n; i++) {
            a[i] = low[i] / sd[i];
            for (int j = 0; j < n; j++) {
                r[i][j] = R[i][j] / (sd[i] * sd[j]);
            }
        }
        double[] tmean = new double[n];
        if (n == 1) {
            if (a[0] == Double.NEGATIVE_INFINITY) {
                tmean[0] = 0;
            }
            else if (a[0] == Double.POSITIVE_INFINITY) {
                tmean[0] = Double.POSITIVE_INFINITY;
            }
            else {
                tmean[0] = Pdf.nor(a[0]) / (1 - Cdf.nor(a[0])) * sd[0];
            }
        }
        else if (n == 2) {
            double rho = r[1][0];
            double lambda = 1 / Math.sqrt(1 - rho * rho);
            //double p = Pdf.bvnor(a[0], a[1], rho);
            double Za = Pdf.nor(a[0]);
            double Zb = Pdf.nor(a[1]);
            double Qa = Cdf.nor(-lambda * (a[1] - rho * a[0]));
            double Qb = Cdf.nor(-lambda * (a[0] - rho * a[1]));
            tmean[0] = (Za * Qa + rho * Zb * Qb) * sd[0];
            tmean[1] = (Zb * Qb + rho * Za * Qa) * sd[1];
        }
        else {
            for (int i = 0; i < n; i++) {
                tmean[i] = 0;
            }
            double[][] rinv = La.cholsl(r);
            double[] bq = new double[n - 1];
            for (int s = 0; s < n - 1; s++) {
                bq[s] = Double.POSITIVE_INFINITY;
            }
            for (int q = 0; q < n; q++) {
                if (a[q] != Double.NEGATIVE_INFINITY) {
                    double[][] rq = La.cholsl(La.removeI(rinv, q));
                    for (int ii = 0; ii < n - 1; ii++) {
                        double aa = 1 / Math.sqrt(rq[ii][ii]);
                        for (int jj = 0; jj < n - 1; jj++) {
                            rq[ii][jj] *= aa;
                            rq[jj][ii] *= aa;
                        }
                    }
                    double[] aq = new double[n - 1];
                    for (int s = 0; s < n; s++) {
                        if (s != q) {
                            aq[s < q ? s : s - 1] = (a[s] - r[s][q] * a[q]) / Math.sqrt(
                                    1 - r[s][q] * r[s][q]);
                        }
                    }
                    double aa = Pdf.nor(a[q]) * Cdf.mvnor(aq, bq, rq);
                    for (int i = 0; i < n; i++) {
                        tmean[i] += r[i][q] * aa;
                    }
                }
            }
            for (int i = 0; i < n; i++) {
                tmean[i] *= sd[i];
            }
        }
        return tmean;
    }


    /**
     *  used internaly
     *
     *@param  a  double[]
     *@param  b  double[]
     *@param  r  double[][]
     *@return    double[]
     */
    private static double[] mvntrcMomentPrivate(double[] a, double[] b, double[][] r) {
        int i = -1;
        for (int k = 0; k < a.length; k++) {
            if (b[k] != Double.POSITIVE_INFINITY) {
                i = k;
                break;
            }
        }
        if (i == -1) {
            return mvntrcMomentPrivate(a, r);
        }
        return La.minus(mvntrcMomentPrivate(a, La.changeI(b, i, Double.POSITIVE_INFINITY),
                r), mvntrcMomentPrivate(La.changeI(a, i, b[i]), La.changeI(b,
                i, Double.POSITIVE_INFINITY), r));
    }
}

