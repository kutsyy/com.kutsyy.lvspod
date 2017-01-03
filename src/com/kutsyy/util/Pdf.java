package com.kutsyy.util;

/**
 * This class has a collection of PDF<BR>
 * Created by <A href="http://www.kutsyy.com">Vadim Kutsyy</A><BR>
 * @author <A href="http://www.kutsyy.com">Vadim Kutsyy</A>
 */
public final class Pdf {



    /**
         * Pdf of bivariate normal
         * @return value of pdf
         * @param x  x={fist dimension, second dimension}
         * @param cov covariance matrix
         */
    public static double bvnor(double[] x, double[][] cov) {
        double[] sd = {Math.sqrt(cov[0][0]), Math.sqrt(cov[1][1])};
        return bvnor(x[0] / sd[0], x[1] / sd[1], cov[0][1] / (sd[1] * sd[0]));
    }



    /**
         * Pdf of bivariate normal
         * @return  value of pdf
         * @param x first dimension
         * @param y second dimension
         * @param r correlation coefficient
         */
    public static double bvnor(double x, double y, double r) {
        if (Double.isInfinite(x) || Double.isInfinite(y))
            return 0;
        return Math.exp(- (x * x + y * y - 2.0 * r * x * y) / (2.0 * (1 - r * r)))
            / (2.0 * Math.PI * Math.sqrt(1.0 - r * r));
    }



    /**
         * Pdf of multuvariate normal
         * @return value of pdf
         * @param x values
         * @param covariance matrix
         */
    public static double mvnor(double[] x, double[][] sigma) {
        double[][] s = new double[x.length][x.length];
        double det = La.choldet(sigma, s);
        double t = La.times(La.times(s, x), x);
        t = Math.exp(-t / 2.0);
        t /= Math.pow(2 * Math.PI, x.length / 2.0);
        t /= Math.sqrt(det);
        return t;
    }



    /**
         * Pdf of standard normal
         * @return pdf
         * @param x  x
         */
    public static double nor(double x) {
        if (Double.isInfinite(x))
            return 0;
        return Math.exp(-x * x / 2.0) / Math.sqrt(2.0 * Math.PI);
    }



    /**
         * Pdf of normal
         * @return pdf
         * @param x  x
         * @param mu mean
         */
    public static double nor(double x, double mu) {
        return nor(x, mu, 1);
    }



    /**
         * Pdf of normal
         * @return pdf
         * @param x x
         * @param mu mean
         * @param sd standart deviation
         */
    public static double nor(double x, double mu, double sd) {
        //return nor((x-mu)/sd);
        return Math.exp(- (x - mu) * (x - mu) / (2.0 * sd * sd))
            / (Math.sqrt(2.0 * Math.PI) * sd);
    }



    /**
         * Pdf of standard truncated normal
         * @return pdf
         * @param x x
         * @param low lower limit of truncation
         */
    public static double ntrc(double x, double low) {
        return nor(x) / (1 - Cdf.nor(low));
    }



    /**
         * Pdf of truncated normal
         * @return pdf
         * @param x x
         * @param low lower limit of truncation
         * @param high upper limit of truncation
         */
    public static double ntrc(double x, double low, double high) {
        return nor(x) / (Cdf.nor(high) - Cdf.nor(low));
    }



    /**
         * Pdf of truncated normal
         * @return pdf
         * @param x x
         * @param low lower limit of truncation
         * @param high upper limit of truncation
         * @param mu mean
         * @param sd standart deviation
         */
    public static double ntrc(
        double x,
        double low,
        double high,
        double mu,
        double sd) {
        low -= mu;
        high -= mu;
        low /= sd;
        high /= sd;
        return nor(x) / (Cdf.nor(high) - Cdf.nor(low));
    }
}