package com.kutsyy.util;

import java.util.Arrays;
/**
 * Some additional Math/Stat functions (I could not name this class Math)<BR>
 * Created by <A href="http://www.kutsyy.com">Vadim Kutsyy</A><BR>
 * @author <A href="http://www.kutsyy.com">Vadim Kutsyy</A>
 */
public final class Stat {



    /**
         * Compute covariance matrix of colums of x
         * @return double[][] - covariance matrix
         * @param x double[][] - x
         */
    public static double[][] cov(double[][] x) {
        double[][] X = La.t(x);
        double[] m = new double[X.length];
        for (int i = 0; i < X.length; i++)
            m[i] = mean(X[i]);
        double[][] c = new double[X.length][X.length];
        for (int i = 0; i < X.length; i++)
            for (int j = 0; j < X.length; j++)
                for (int k = 0; k < X[0].length; k++)
                    c[i][j] += (X[i][k] - m[i]) * (X[j][k] - m[j]);
        for (int i = 0; i < X.length; i++)
            for (int j = 0; j < X.length; j++)
                c[i][j] /= (X[0].length - 1);
        return c;
    }



    public static double erf(double x) {
        return calerf(x, 0);
    }



    public static double erfc(double x) {
        return calerf(x, 1);
    }



    /**
         * Mean of row of the matrix x
         * @return double[] - mean(x)
         * @param x double[][] - x
         */
    public static double[] mean(double[][] x) {
        double[][] X = La.t(x);
        double[] m = new double[X.length];
        for (int i = 0; i < X.length; i++)
            m[i] = mean(X[i]);
        return m;
    }



    /**
         * Mean of the vector
         * @return double - mean(a)
         * @param a double[] - a
         */
    public static double mean(double[] a) {
        return (sum(a) / a.length);
    }



    /**
         * Median of the vector
         * @return double - median(a)
         * @param a double[] - a
         */
    public static double median(double[] a) {
        return percentile(.5, a);
    }



    /**
         * Find percetiles of given (unsorted) vector
         * @return double[] - quantiles
         * @param p double[] - percentile at with quaniles needed
         * @param a double[] - vector
         */
    public static double[] percentile(double[] p, double[] a) {
        double[] b = (double[]) a.clone();
        Arrays.sort(b);
        double[] t = new double[p.length];
        for (int i = 0; i < p.length; i++)
            t[i] =
                (b[(int) Math.floor(a.length * p[i])] + b[(int) Math.ceil(a.length * p[i])])
                    / 2.0;
        return t;
    }



    /**
         * Find quontale of given (unsorted) percentile
         * @return double - quantile
         * @param p double - percentile at with quanile needed
         * @param a double[] - vector
         */
    public static double percentile(double p, double[] a) {
        double[] b = (double[]) a.clone();
        Arrays.sort(b);
        return (b[(int) Math.floor(a.length * p)] + b[(int) Math.ceil(a.length * p)])
            / 2.0;
    }



    /**
         * Sum of the vector
         * @return double - sum(a)
         * @param a double[] - a
         */
    public static double sum(double[] a) {
        double sum = 0;
        for (int i = 0; i < a.length; i++)
            sum += a[i];
        return sum;
    }



    /**
         * Sum of squares of the vector
         * @return double - sumsq(a)
         * @param a double[] - a
         */
    public static double sumsq(double[] a) {
        double sum = 0;
        for (int i = 0; i < a.length; i++)
            sum += a[i] * a[i];
        return sum;
    }



    /**
         * Variance of the vector
         * @return double - var(a)
         * @param a double[] - a
         */
    public static double var(double[] a) {
        double sum = sum(a);
        return ((sumsq(a) - sum * sum / a.length) / a.length);
    }



    private static double calerf(double x, int jint) {
        /*
        Intro is copied from original fortran version
        C------------------------------------------------------------------
        C
        C THIS PACKET COMPUTES THE ERROR AND COMPLEMENTARY ERROR FUNCTIONS
        C   FOR REAL ARGUMENTS  ARG.  IT CONTAINS TWO FUNCTION TYPE
        C   SUBPROGRAMS,  ERF  AND  ERFC  (OR  DERF  AND  DERFC),  AND ONE
        C   SUBROUTINE TYPE SUBPROGRAM,  CALERF.  THE CALLING STATEMENTS
        C   FOR THE PRIMARY ENTRIES ARE
        C
        C                   Y=ERF(X)     (OR   Y=DERF(X) )
        C   AND
        C                   Y=ERFC(X)    (OR   Y=DERFC(X) ).
        C
        C   THE ROUTINE  CALERF  IS INTENDED FOR INTERNAL PACKET USE ONLY,
        C   ALL COMPUTATIONS WITHIN THE PACKET BEING CONCENTRATED IN THIS
        C   ROUTINE.  THE FUNCTION SUBPROGRAMS INVOKE  CALERF  WITH THE
        C   STATEMENT
        C          CALL CALERF(ARG,RESULT,JINT)
        C   WHERE THE PARAMETER USAGE IS AS FOLLOWS
        C
        C      FUNCTION                     PARAMETERS FOR CALERF
        C       CALL              ARG                  RESULT          JINT
        C     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
        C     ERFC(ARG)     ABS(ARG) .LT. XMAX        ERFC(ARG)         1
        C
        C   THE MAIN COMPUTATION EVALUATES NEAR MINIMAX APPROXIMATIONS
        C   FROM "RATIONAL CHEBYSHEV APPROXIMATIONS FOR THE ERROR FUNCTION"
        C   BY W. J. CODY, MATH. COMP., 1969, PP. 631-638.  THIS
        C   TRANSPORTABLE PROGRAM USES RATIONAL FUNCTIONS THAT THEORETICALLY
        C       APPROXIMATE  ERF(X)  AND  ERFC(X)  TO AT LEAST 18 SIGNIFICANT
        C   DECIMAL DIGITS.  THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC
        C   SYSTEM, THE COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER
        C   SELECTION OF THE MACHINE-DEPENDENT CONSTANTS.
        C
        C  AUTHOR: W. J. CODY
        C          MATHEMATICS AND COMPUTER SCIENCE DIVISION
        C          ARGONNE NATIONAL LABORATORY
        C          ARGONNE, IL 60439
        C
        C  LATEST MODIFICATION: JANUARY 8, 1985
        C
        C------------------------------------------------------------------
        */
        double a[] =
            {
                3.16112374387056560e00,
                1.13864154151050156e02,
                3.77485237685302021e02,
                3.20937758913846947e03,
                1.85777706184603153e-1};
        double b[] =
            {
                2.36012909523441209e01,
                2.44024637934444173e02,
                1.28261652607737228e03,
                2.84423683343917062e03};
        double c[] =
            {
                .564188496988670089,
                8.88314979438837594,
                66.1191906371416295,
                298.635138197400131,
                881.95222124176909,
                1712.04761263407058,
                2051.07837782607147,
                1230.33935479799725,
                2.15311535474403846e-8};
        double d[] =
            {
                15.7449261107098347,
                117.693950891312499,
                537.181101862009858,
                1621.38957456669019,
                3290.79923573345963,
                4362.61909014324716,
                3439.36767414372164,
                1230.33935480374942};
        double p[] =
            {
                .305326634961232344,
                .360344899949804439,
                .125781726111229246,
                .0160837851487422766,
                6.58749161529837803e-4,
                .0163153871373020978};
        double q[] =
            {
                2.56852019228982242,
                1.87295284992346047,
                .527905102951428412,
                .0605183413124413191,
                .00233520497626869185};
        double sqrpi = Math.sqrt(Math.PI);
        double thresh = .46875;
        double xmax = 200;
        /*for x>xmax, erf=1*/
        double xsmall = 4.2e-30;
        /*for x<xsmal erf=0*/
        double result, ysq, xnum, xden;
        double y = Math.abs(x);
        if (y > 6) {
            result = 0;
            /*200*/
            if (y <= xmax) {
                ysq = 1 / (y * y);
                /*220*/
                xnum = p[5] * ysq;
                xden = ysq;
                for (int i = 0; i < 4; i++) {
                    xnum = (xnum + p[i]) * ysq;
                    xden = (xden + q[i]) * ysq;
                };
                /*240*/
                result = ysq * (xnum + p[4]) / (xden + q[4]);
                result = Math.exp(-y * y) / y * (sqrpi - result);
            }
        } else
            if (y > thresh) {
                ysq = y * y;
                /*100*/
                xnum = c[8] * y;
                xden = y;
                for (int i = 0; i < 7; i++) {
                    xnum = (xnum + c[i]) * y;
                    xden = (xden + d[i]) * y;
                };
                /*120*/
                result = Math.exp(-ysq) * (xnum + c[7]) / (xden + d[7]);
            } else {
                ysq = 0;
                if (y > xsmall) {
                    ysq = y * y;
                };
                xnum = a[4] * ysq;
                xden = ysq;
                for (int i = 0; i < 3; i++) {
                    xnum = (xnum + a[i]) * ysq;
                    xden = (xden + b[i]) * ysq;
                }
                /*20*/
                result = x * (xnum + a[3]) / (xden + b[3]);
                if (jint == 0) {
                    return result;
                } else {
                    return (1 - result);
                };
            }
        if (jint == 0) /*300*/ {
            result = 0.5 - result + 0.5;
            /*350*/
            if (x < 0)
                result = -result;
        } else {
            if (x < 0)
                result = 2 - result;
        }
        return result;
    }
}