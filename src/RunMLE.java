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
public class RunMLE {

	/**
	 *  Starts the application.
	 *
	 *@param  args  an array of command-line arguments
	 */
	public static void main(java.lang.String[] args) {
		String homedir = "";
		//"/var/tmp/kutsyy";
		int K = 2;
		int N = 20000;
		int n = 20;
		int[][] nghb = {
				{
				-1, 1, 0, 0
				}, {
				0, 0, -1, 1
				}
				};
		double[] phi =
				{0};
		//{0, .1, .2, .3, .4, .5, .6, .7, -.1, -.2, -.3, -.4, -.5, -.6, -.7, 1 / 3.0, 2 / 3.0};
		double[] theta = new double[K];
		for (int i = 0; i < K; i++) {
			theta[i] = CdfInv.nor((i + 1) / (double) (K + 1));
		}
		long t;
		OrdinalSpatial T;
		try {
			PrintWriter out = new PrintWriter(new FileWriter(new File(homedir, "sim." + K + "." + n + ".0.QLE.txt")));
			out.print("size i phi ");
			out.println("PhiMLE llMLE llMLE0 timeMLE");
			for (int k = 0; k < phi.length; k++) {
				T = new OrdinalSpatial(n, n, theta, phi[k], nghb);
				for (int j = 0; j < N; j++) {
					out.print(n + " " + j + " " + phi[k] + " ");
					for (int i = 0; i < T.TrueX.length; i++) {
						T.TrueX[i] = T.rnd.gaussian();
					}
					T.phi = 0;
					T.theta = (double[]) theta.clone();
					t = System.currentTimeMillis();
					T.mCLETrueX();
					out.print(T.phi + " " + T.lCLTrueX() + " " + T.lCLTrueX(0) +
							" " + (System.currentTimeMillis()
							 - t) + " ");
					T.phi = 0;
					T.theta = (double[]) theta.clone();
					t = System.currentTimeMillis();
					T.mMnA();
					out.print(T.phi + " " + T.lA("mean") + " ");
					T.mMnA(0);
					out.print(T.lA("mean") + " " + (System.currentTimeMillis() /
							-t) + " ");
					out.println();
				}
			}
			out.close();
		}
		catch (IOException e) {
		}

	}
}

