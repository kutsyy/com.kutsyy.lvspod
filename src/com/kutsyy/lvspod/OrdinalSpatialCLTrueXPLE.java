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
public final class OrdinalSpatialCLTrueXPLE extends OrdinalSpatialAbstract {

    /**
     *  Insert the method's description here.
     */
    protected OrdinalSpatialCLTrueXPLE() {
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
    protected OrdinalSpatialCLTrueXPLE(double[] trueX, double Phi, SpatialPoint[] loc,
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
    }


    /**
     *  Insert the method's description here.
     *
     *@return    double
     */
    public double lL() {
        double[] xz = new double[TrueX.length];
        if (K > 0) {
            xz = La.times(Z, beta);
        }
        for (int i = 0; i < N; i++) {
            xz[i] += TrueX[i];
        }
        double s = 0;
        double sd = 1 / Math.sqrt(1 - phi * phi);
        if (phi != 0) {
            for (int i = 0; i < N; i++) {
                if (Loc[i].neighbors.length > 0) {
                    double t = 0;
                    for (int j = 0; j < Loc[i].neighbors.length; j++) {
                        t += xz[Loc[i].neighbors[j]];
                    }
                    t *= phi / (double) Loc[i].neighbors.length;
                    s += Math.log(Pdf.nor(xz[i], t, Math.sqrt(1 - phi * phi / (double) Loc[i].neighbors.length)));
                }
                else {
                    s += Math.log(Pdf.nor(xz[i], 0, sd));
                }
            }
        }
        else {
            for (int i = 0; i < N; i++) {
                s += Math.log(Pdf.nor(xz[i]));
            }
        }
        return s;
    }
}

