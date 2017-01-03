/*
 * put your module comment here
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 */


package  com.kutsyy.lvspod;


/**
 * Title:        com.kutsyy.lvspod
 * Description:  This library is set of rutins for Laten Variable model for Spatialy dependent Ordinal Data
 * Copyright:    Copyright (c) 2000
 * Company:      The University of Michigan
 * @author Vadim Kutsyy
 * @version 1.0
 */
public class OrdinalSpatialSingleAveregeQL extends OrdinalSpatialAbstract {
    private OrdinalSpatialQL[] T;
    private boolean todefine = true;

    /**
     * put your documentation comment here
     * @param     int[] y
     */
    public OrdinalSpatialSingleAveregeQL (int[] y) {
        super(y);
    }

    /**
     * put your documentation comment here
     * @param     int[] y
     * @param     double[][] z
     */
    public OrdinalSpatialSingleAveregeQL (int[] y, double[][] z) {
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
    protected OrdinalSpatialSingleAveregeQL (int[] y, double phi, double[] theta,
            com.kutsyy.lvspod.SpatialPoint[] loc, double[] beta, double[][] z,
            int[][] neighborDefinition) {
        super(y, phi, theta, loc, beta, z, neighborDefinition);
    }

    /**
     * put your documentation comment here
     */
    protected OrdinalSpatialSingleAveregeQL () {
        super();
    }

    /**
     * put your documentation comment here
     * @param     int[] y
     * @param     int[] xLocation
     * @param     int[] yLocation
     * @param     double[][] z
     */
    protected OrdinalSpatialSingleAveregeQL (int[] y, int[] xLocation, int[] yLocation,
            double[][] z) {
        super(y, xLocation, yLocation, z);
    }

    /**
     * put your documentation comment here
     * @return
     */
    public double lL () {
        double l = 0;
        for (int i = 0; i < T.length; i++) {
            T[i].phi = phi*(double)T.length;
            for (int j = 0; j < L; j++)
                T[i].theta[j] = theta[j];
            for (int j = 0; j < K; j++)
                T[i].beta[j] = beta[j];
            l += T[i].lL();
        }
 return  l;
    }

    /**
     * put your documentation comment here
     */
    protected void define () {
        T = new OrdinalSpatialQL[NeighborDefinition[0].length];
        int[][] thisNeighborDefinition = new int[2][1];
        for (int i = 0; i < NeighborDefinition[0].length; i++) {
            thisNeighborDefinition[0][0] = NeighborDefinition[0][i];
            thisNeighborDefinition[1][0] = NeighborDefinition[1][i];
            T[i] = new OrdinalSpatialQL(Y, xLoc, yLoc, Z, thisNeighborDefinition);
        }
    }

    /**
     * put your documentation comment here
     * @param y
     */
    public void updateY (int[] y) {
        super.updateY(y);
        for (int i = 0; i < NeighborDefinition[0].length; i++)
            T[i].updateY(y);
    }
}



