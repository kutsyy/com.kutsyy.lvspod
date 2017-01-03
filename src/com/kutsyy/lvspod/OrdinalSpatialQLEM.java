/*
 * put your module comment here
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 */


package  com.kutsyy.lvspod;

import  com.kutsyy.util.*;


/**
 * Insert the type's description here.
 * Creation date: (10/21/2000 2:56:53 PM)
 * @author:
 */
public final class OrdinalSpatialQLEM extends OrdinalSpatialAbstractEM {

    /**
     * put your documentation comment here
     * @param     int[] y
     */
    public OrdinalSpatialQLEM (int[] y) {
        super(y);
    }

    /**
     * put your documentation comment here
     * @param     int[] y
     * @param     double[][] z
     */
    public OrdinalSpatialQLEM (int[] y, double[][] z) {
        super(y, z);
    }

    /**
     * put your documentation comment here
     * @param     int[] y
     * @param     double phi
     * @param     double[] theta
     * @param     com.kutsyy.lvspod.SpatialPoint[] loc
     * @param     double[] beta
     * @param     double[][] z
     * @param     int[][] neighborDefinition
     */
    protected OrdinalSpatialQLEM (int[] y, double phi, double[] theta, com.kutsyy.lvspod.SpatialPoint[] loc,
            double[] beta, double[][] z, int[][] neighborDefinition) {
        super(y, phi, theta, loc, beta, z, neighborDefinition);
    }

    /**
     * put your documentation comment here
     */
    protected OrdinalSpatialQLEM () {
        super();
    }

    /**
     * put your documentation comment here
     * @param     int[] y
     * @param     int[] xLocation
     * @param     int[] yLocation
     * @param     double[][] z
     */
    protected OrdinalSpatialQLEM (int[] y, int[] xLocation, int[] yLocation,
            double[][] z) {
        super(y, xLocation, yLocation, z);
    }

    /**
     * put your documentation comment here
     */
    protected void e () {
        generateX();
    }

    /**
     * Insert the method's description here.
     * Creation date: (10/10/2000 6:08:50 PM)
     */
    private void generateX () {
        double[] thetaL = new double[N];
        double[] thetaU = new double[N];
        for (int i = 0; i < N; i++)
            if (Y[i] == 0) {
                thetaU[i] = (theta[0]);
                thetaL[i] = Double.NEGATIVE_INFINITY;
            }
            else if (Y[i] == L) {
                thetaU[i] = Double.POSITIVE_INFINITY;
                thetaL[i] = theta[L - 1];
            }
            else {
                thetaL[i] = theta[Y[i] - 1];
                thetaU[i] = theta[Y[i]];
            }
        if (K > 0) {
            x = La.times(Z, beta);
            for (int i = 0; i < N; i++) {
                thetaL[i] += x[i];
                thetaU[i] += x[i];
            }
        }
        for (int i = 0; i < N; i++) {
            double[][] s = new double[Loc[i].neighbors.length + 1][Loc[i].neighbors.length
                    + 1];
            double t1 = phi/Loc[i].neighbors.length;
            double t2 = t1*t1;
            for (int j = 1; j <= Loc[i].neighbors.length; j++) {
                s[j][0] = s[0][j] = t1;
                for (int k = j + 1; k <= Loc[i].neighbors.length; k++)
                    s[j][k] = s[k][j] = t2;
            }
            for (int j = 0; j <= Loc[i].neighbors.length; j++)
                s[j][j] = 1;
            double[] thetal = new double[Loc[i].neighbors.length + 1];
            double[] thetau = new double[Loc[i].neighbors.length + 1];
            for (int j = 0; j < Loc[i].neighbors.length; j++) {
                thetal[j + 1] = thetaL[Loc[i].neighbors[j]];
                thetau[j + 1] = thetaU[Loc[i].neighbors[j]];
            }
            thetal[0] = thetaL[i];
            thetau[0] = thetaU[i];
            x[i] = Exp.mvntrc(thetal, thetau, s)[0];
        }
    }
}



