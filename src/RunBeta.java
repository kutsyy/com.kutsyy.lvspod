/*
 * put your module comment here
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 *
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 *
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 */

import com.kutsyy.util.*;
import com.kutsyy.lvspod.*;
import java.io.*;

/**
 *  Title: Description: Copyright: Copyright (c) 2000 Company:
 *
 *@author
 *@created    December 2, 2000
 *@version    1.0
 */
public class RunBeta {

    /**
     *  Starts the application.
     *
     *@param  args  an array of command-line arguments
     */
    public static void main(java.lang.String[] args) {
        String homedir = "/var/tmp/kutsyy";
        int K = 4;
        int N = 1000;
        int n = 20;
        int[][] nghb = {
                {
                -1, 1, 0, 0
                }, {
                0, 0, -1, 1
                }
                };
        double[] Phi = {0, 0.2, 0.4, 0.6};
        //double[] Phi = {0.1,0.3,0.5,0.7};
        // Double.parseDouble(args[0]);
        double[][] beta = {{0}, {.1}, {.2}, {.3}, {.4}, {.5}, {.6}, {.7}, {.8}, {.9}, {1.25}, {1.5}, {1.75}, {2}};
        double[] betaNull = {0};
        double[] theta = new double[K];
        for (int i = 0; i < K; i++) {
            theta[i] = CdfInv.nor((i + 1) / (double) (K + 1));
        }
        double[][][] z = new double[n][n][1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                z[i][j][0] = (i > (n / (double) 2)) ? 1 : 0;
            }
        }
        OrdinalSpatial T;
        OrdinalSpatial Th;
        OrdinalSpatial Tv;
        double th;
        double tv;
        long t;
        for (int i1 = 0; i1 < Phi.length; i1++) {
            double phi = Phi[i1];
            try {
                //PrintWriter summary = new PrintWriter(new FileWriter(new File(homedir, "sim." + K + "." + n + ".half.beta.summary.txt")));
                PrintWriter out = new PrintWriter(new FileWriter(new File(homedir, "sim." + K + "." + n + "." + phi + ".half.beta.txt")));
                out.print("size i phi beta ");
                for (int i = 0; i < K; i++) {
                    out.print("Theta" + i + " ");
                }
                //out.print("PhiMLE llMLE llMLE0 timeMLE PhiMMnA llMMnA llMMnAPhi0 llMMnABeta0 llMMnAPhi0Beta0");
                //out.println("timeMMnA PhiMdnA llMMdA llMMdAPhi0 llMMdABeta0 llMMdAPhi0Beta0 timeMMdA");
                for (int l = 0; l < beta.length; l++) {
                    T = new OrdinalSpatial(n, n, theta, phi, beta[l], z, nghb);
                    double meanMLETrueX = 0;
                    double meanMnA = 0;
                    double meanMdA = 0;
                    double varMLETrueX = 0;
                    double varMnA = 0;
                    double varMdA = 0;
                    for (int j = 0; j < N; j++) {
                        out.print(n + " " + j + " " + phi + " " + beta[l][0] + " ");
                        for (int jj = 0; jj < K; jj++) {
                            out.print(theta[jj] + " ");
                        }
                        T.phi = 0;
                        T.theta = (double[]) theta.clone();
                        t = System.currentTimeMillis();
                        T.mCLETrueX();
                        out.print(T.phi + " " + T.beta[0] + T.lCLTrueX() + " ");
                        T.mCLETrueX(0);
                        out.print(T.phi + " " + T.beta[0] + T.lCLTrueX() + " ");
                        T.mCLETrueX(betaNull);
                        out.print(T.phi + " " + T.beta[0] + T.lCLTrueX() + " 0 0 "
                                 + T.lCLTrueX(0, betaNull) + " ");
                        T.phi = 0;
                        T.theta = (double[]) theta.clone();
                        t = System.currentTimeMillis();
                        T.mMnA();
                        double a = T.lA("mean");
                        out.print(T.phi + " " + T.beta[0] + " " + a + " ");
                        T.mMnA(0);
                        double b = T.lA("mean");
                        out.print(T.phi + " " + T.beta[0] + " " + b + " " + (a - b) + " ");
                        T.mMnA(null, betaNull);
                        b = T.lA("mean");
                        out.print(T.phi + " " + T.beta[0] + " " + b + " " + (a - b) + " ");
                        T.mMnA(0, null, betaNull);
                        b = T.lA("mean");
                        out.print(T.phi + " " + T.beta[0] + " " + b + " " + (a - b) + " ");
                        T.generateY();
                        out.println();
                    }
                    out.close();
                }
            }
            catch (IOException e) {
            }
        }
    }
}

