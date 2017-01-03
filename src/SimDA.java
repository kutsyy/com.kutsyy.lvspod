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
public class SimDA {

	/**
	 *  Description of the Class
	 *
	 *@author     Vadim Kutsyy
	 *@created    January 9, 2001
	 */
	public SimDA() {
	}


	/**
	 *  Description of the Method
	 */
	public synchronized void main() {
		String homedir = ".";
		int thread = 1;
		int K = 2;
		int N = 100;
		int n = 20;
		int[][] nghb = {{-1, 1, 0, 0}, {0, 0, -1, 1}};
		double[] phi = {1 / 3.0};
		//{0, .1, .2, .3, .4, .5, .6, .7, -.1, -.2, -.3, -.4, -.5, -.6, -.7, 1 / 3.0, 2 / 3.0};
		double[] theta = new double[K];
		for (int i = 0; i < K; i++) {
			theta[i] = CdfInv.nor((i + 1) / (double) (K + 1));
		}
		tThread[] TThread = new tThread[thread];
		//OrdinalSpatial T;
		OrdinalSpatialPlot plotField = new OrdinalSpatialPlot(4);
		try {
			PrintWriter out = new PrintWriter(new FileWriter(
					new File(homedir, "sim." + K + "." + n + ".txt")
					));
			System.out.print("size i phi ");
			for (int i = 0; i < K; i++) {
				System.out.print("Theta" + i + " ");
			}
			System.out.println("PhiMLE PhiDA timeDA");
			for (int k = 0; k < phi.length; k++) {
				int count = 0;
				//T = new OrdinalSpatial(n, n, theta, phi[k], nghb);
				for (int i = 0; i < thread; i++) {
//					T.generateY();
					TThread[i] = new tThread();
					TThread[i].T = new OrdinalSpatial(n, n, theta, phi[k], nghb);
					TThread[i].lock = true;
                                        plotField.add(TThread[i].T);
					//TThread[i] = new tThread(T, 0.0, T.theta);
					TThread[i].start();
				}
				while (count < N) {
					try {
						wait(10);
					}
					catch (InterruptedException e) {
						// do nothing
					}
					for (int i = 0; i < thread; i++) {
						if (!TThread[i].lock) {
							System.out.print(n + " " + phi[k] + " ");
							for (int jj = 0; jj < K; jj++) {
								System.out.print(theta[jj] + " ");
							}
							// for(double jj=-.9;jj<1;jj+=.1) System.out.println(jj+" "+T.lCLTrueX(jj));
//                        T.phi = 0;
//                        T.theta = (double[]) theta.clone();
//                        long t = System.currentTimeMillis();
//                        T.mCLETrueX();
//                        System.out.print(T.phi + " ");
//                        T.phi = 0;
//                        T.theta = (double[]) theta.clone();
//                        t = System.currentTimeMillis();
//                        T.dataAgumentation();
//                        System.out.print(T.phi + " " + (System.currentTimeMillis() - t));
//                        T.generateY();
							System.out.println(TThread[i].PhimCLETrueX + " " + TThread[i].PhiMain + " " + TThread[i].t);
							count++;
							TThread[i].T.generateY();
							TThread[i].T.phi = 0;
							TThread[i].T.theta = (double[]) theta.clone();
							// TThread[i] = new tThread(T, 0.0, T.theta);
							TThread[i].lock = true;
                                                        plotField.add(TThread[i].T);
							TThread[i].start();
							try {
								wait(10);
							}
							catch (InterruptedException e) {
								// do nothing
							}
						}

					}
				}
			}
//			summary.close();
			out.close();
		}
		catch (IOException e) {
		}

	}



	/**
	 *  Description of the Class
	 *
	 *@author     kutsyy
	 *@created    January 9, 2001
	 */
	private class tThread extends Thread {
		/**
		 *  Description of the Field
		 */
		public OrdinalSpatial T;
		boolean lock = false;
		double PhimCLETrueX;
		double PhiMain;
		long t;


		/**
		 *  Description of the Method
		 *
		 */
		public tThread() {
		}


		/**
		 *  Constructor for the tThread object
		 *
		 *@param  t  Description of Parameter
		 */
		public tThread(OrdinalSpatial t) {
			T = new OrdinalSpatial(t);
			if (t.TrueX != null) {
				T.TrueX = (double[]) t.TrueX.clone();
			}
		}


		/**
		 *  Description of the Method
		 *
		 *@param  t      Description of Parameter
		 *@param  Phi    Description of Parameter
		 *@param  Theta  Description of Parameter
		 */
		public tThread(OrdinalSpatial t, double Phi, double[] Theta) {
			T = new OrdinalSpatial(t);
			if (t.TrueX != null) {
				T.TrueX = (double[]) t.TrueX.clone();
			}
			T.phi = Phi;
			T.theta = (double[]) Theta.clone();
			lock = true;
		}


		/**
		 *  Main processing method for the tTread object
		 */
		public void run() {
			lock = true;
			t = System.currentTimeMillis();
			PhimCLETrueX = 0;
			if (T.TrueX != null) {
				T.mCLETrueX();
				PhimCLETrueX = T.phi;
			}
			T.dataAgumentation(true);
			PhiMain = T.phi;
			t = System.currentTimeMillis() - t;
			lock = false;
		}
	}
}

