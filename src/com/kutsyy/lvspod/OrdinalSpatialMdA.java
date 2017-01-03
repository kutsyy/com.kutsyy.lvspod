/*
 * put your module comment here
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 */


package  com.kutsyy.lvspod;

import  com.kutsyy.util.*;


/**
 * Insert the type's description here.
 * Created by <A href="http://www.kutsyy.com">Vadim Kutsyy</A><BR>
 * @author <A href="http://www.kutsyy.com">Vadim Kutsyy</A>
 */
public final class OrdinalSpatialMdA extends OrdinalSpatialAbstractSubstitution {

    /**
     * Insert the method's description here.
     */
    public OrdinalSpatialMdA () {
    }

    /**
     * put your documentation comment here
     */
    protected final void generateX () {
        if (K > 0) {
            x = La.times(Z, beta);
            for (int i = 0; i < N; i++) {
                thetaL[i] += x[i];
                thetaU[i] += x[i];
            }
        }
        for (int i = 0; i < N; i++)
            x[i] = CdfInv.nor((Cdf.nor(thetaL[i]) + Cdf.nor(thetaU[i]))/2.0);
    }

    /**
     * Insert the method's description here.
     * @param y int[]
     */
    public OrdinalSpatialMdA (int[] y) {
        super(y);
    }

    /**
     * Insert the method's description here.
     * @param y int[]
     * @param z double[][]
     */
    public OrdinalSpatialMdA (int[] y, double[][] z) {
        super(y, z);
    }

    /**
     * Insert the method's description here.
     * @param Y int[]
     * @param XLocation int[]
     * @param YLocation int[]
     */
    public OrdinalSpatialMdA (int[] y, int[] xLocation, int[] yLocation, double[][] z,
            int[][] Neighborgood) {
        super(y, xLocation, yLocation, z, Neighborgood);
    }

    /**
     * OrdinalSpatialEM constructor comment.
     * @param y int[]
     * @param phi double
     * @param theta double[]
     * @param loc com.kutsyy.lvspod.SpatialPoint[]
     * @param beta cern.colt.matrix.DoubleMatrix1D
     * @param z cern.colt.matrix.DoubleMatrix2D
     * @param neighborDefinition int[][]
     */
    protected OrdinalSpatialMdA (int[] y, double phi, double[] theta, com.kutsyy.lvspod.SpatialPoint[] loc,
            double[] beta, double[][] z, int[][] neighborDefinition) {
        super(y, phi, theta, loc, beta, z, neighborDefinition);
    }
}



