/*
 * put your module comment here
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 */

package com.kutsyy.lvspod;

import com.kutsyy.util.*;

/**
 *  Insert the type's description here. Created by <A
 *  href="http://www.kutsyy.com">Vadim Kutsyy</A> <BR>
 *
 *
 *@author     <A href="http://www.kutsyy.com">Vadim Kutsyy</A>
 *@created    December 27, 2000
 */
public final class OrdinalSpatialCLTrueX extends OrdinalSpatialAbstract {
    double[] xl = null;
    double xx, xxl, xlxl;
    double c = -Math.log(Math.PI * 2.0) * N / 2.0;


    /**
     *  Insert the method's description here.
     */
    protected OrdinalSpatialCLTrueX() {
    }


    /**
     *  OrdinalSpatialCMPL constructor comment.
     *
     *@param  loc                 com.kutsyy.ar.util.SpatialPoint[]
     *@param  z                   double[][]
     *@param  trueX               Description of Parameter
     *@param  Phi                 Description of Parameter
     *@param  Beta                Description of Parameter
     *@param  neighborDefinition  Description of Parameter
     */
    protected OrdinalSpatialCLTrueX(double[] trueX, double Phi, SpatialPoint[] loc,
            double[] Beta, double[][] z, int[][] neighborDefinition) {
        this.NeighborDefinition = neighborDefinition;
        TrueX = trueX;
        phi = Phi;
        Loc = loc;
        N = TrueX.length;
        if (Beta != null && Beta.length > 0) {
            this.beta = (double[]) Beta.clone();
            this.Z = z;
            K = beta.length;
        }
        createXX();
    }


    /**
     *  Insert the method's description here.
     *
     *@return    double
     */
    public double lL() {
        double s = 1;
        double ss = 1;
        double v = 1;
        //-phi*phi/(double)NeighborDefinition[0].length;
        for (int i = 0; i < N; i++) {
            s *= (1 - phi * lambda[i]);
            ss *= 1 - (Loc[i].neighbors.length > 0 ? phi * phi / (double) Loc[i].neighbors.length : 0);
        }
        s = Math.log(s / ss) / 2.0;
        s += -(xx - phi * xxl) / (2.0 * (1 - phi * phi / (double) NeighborDefinition.length)) + c;
        return s;
    }


    /**
     *  Description of the Method
     */
    public void updateX(double[] X){
    TrueX=X;
    createXX();
    }
    public void createXX() {
        createEigenvalues();
        x = (double[]) TrueX.clone();
        xl = new double[N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < Loc[i].neighbors.length; j++) {
                xl[i] += x[Loc[i].neighbors[j]];
            }
        }
        for (int i = 0; i < N; i++) {
            xl[i] /= (double) Loc[i].neighbors.length;
        }
        xx = 0;
        xxl = 0;
        for (int i = 0; i < N; i++) {
            xx += x[i] * x[i];
            xxl += x[i] * xl[i];
        }

    }
}

