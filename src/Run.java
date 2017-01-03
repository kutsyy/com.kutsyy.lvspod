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
public class Run {

    /**
     *  Starts the application.
     *
     *@param  args  an array of command-line arguments
     */
    public static void main(java.lang.String[] args) {
        String homedir = ".";
        int K = 2;
        int N = 10000;
        int n = 20;
        int[][] nghb = {
                {
                -1, 1, 0, 0
                }, {
                0, 0, -1, 1
                }
                };
        double[] phi ={0, .1, .2, .3, .4, .5, .6, .7, -.1, -.2, -.3, -.4, -.5, -.6, -.7, 1 / 3.0, 2 / 3.0};
        double[] theta = new double[K];
        for (int i = 0; i < K; i++) {
            theta[i] = CdfInv.nor((i + 1) / (double) (K + 1));
        }
        OrdinalSpatial T;
        OrdinalSpatial Th;
        OrdinalSpatial Tv;
        double th;
        double tv;
        long t;
//		try {
//			PrintWriter summary = new PrintWriter(new FileWriter(new File(homedir, "sim." + K + "." + n + ".summary.txt")));
//			PrintWriter out = new PrintWriter(new FileWriter(new File(homedir, "sim." + K + "." + n + ".txt")));
//			out.print("size i phi ");
//			for (int i = 0; i < K; i++) {
//				out.print("Theta" + i + " ");
//			}
//			out.println("PhiMLE llMLE llMLE0 timeMLE PhiMMnA llMMnA llMMnA0 timeMMnA PhiMdnA llMMdA llMMdA0 timeMMdA");
        for (int k = 0; k < phi.length; k++) {
            T = new OrdinalSpatial(n, n, theta, phi[k], nghb);
            double meanMLETrueX = 0;
            double meanMnA = 0;
            double meanMdA = 0;
            double varMLETrueX = 0;
            double varMnA = 0;
            double varMdA = 0;
            for (int j = 0; j < N; j++) {
                System.out.print(n + " " + j + " " + phi[k] + " ");
                for (int jj = 0; jj < K; jj++) {
                    System.out.print(theta[jj] + " ");
                }
                // for(double jj=-.9;jj<1;jj+=.1) System.out.println(jj+" "+T.lCLTrueX(jj));
                T.phi = 0;
                T.theta = (double[]) theta.clone();
                t = System.currentTimeMillis();
                T.mCLETrueX();
                System.out.print(T.phi + " " +  (System.currentTimeMillis()
                         - t) + " ");
                T.phi = 0;
                T.theta = (double[]) theta.clone();
                t = System.currentTimeMillis();
                T.mCLETrueXPLE();
                System.out.print(T.phi + " " + (System.currentTimeMillis()
                         - t) + " ");
//					T.phi = 0;
//					T.theta = (double[]) theta.clone();
//					t = System.currentTimeMillis();
//					T.mMnA();
//					meanMnA += T.phi;
//					varMnA += T.phi * T.phi;
//					out.print(T.phi + " " + T.lA("mean") + " ");
//					T.mMnA(0);
//					out.print(T.lA("mean") + " " + (System.currentTimeMillis()
//							 - t) + " ");
//					T.phi = 0;
//					T.theta = (double[]) theta.clone();
//					t = System.currentTimeMillis();
//					T.mMdA();
//					meanMdA += T.phi;
//					varMdA += T.phi * T.phi;
//					out.print(T.phi + " " + T.lA("median") + " ");
//					T.mMdA(0);
//					out.print(T.lA("median") + " " + (System.currentTimeMillis()
//							 - t) + " ");
                T.generateY();
                System.out.println();
            }
//				meanMLETrueX /= (double) N;
//				meanMnA /= (double) N;
//				meanMdA /= (double) N;
//				varMLETrueX /= (double) N;
//				varMnA /= (double) N;
//				varMdA /= (double) N;
//				varMLETrueX -= meanMLETrueX * meanMLETrueX;
//				varMnA -= meanMnA * meanMnA;
//				varMdA -= meanMdA * meanMdA;
//				for (int j = 0; j < K; j++) {
//					summary.print(theta[j] + " ");
//				}
//				summary.println(phi[k] + " " + meanMLETrueX + " " + varMLETrueX + " " + meanMnA + " " + varMnA + " " + meanMdA + " " + varMdA);
        }
//			summary.close();
//			out.close();
//		}
//		catch (IOException e) {
//		}

    }
}

