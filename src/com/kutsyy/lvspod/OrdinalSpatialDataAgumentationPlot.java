/*
 * put your module comment here
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 */

package com.kutsyy.lvspod;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import com.kutsyy.util.*;

//import com.borland.jbcl.layout.*;
/**
 *  Title: com.kutsyy.lvspod Description: This library is set of rutins for
 *  Laten Variable model for Spatialy dependent Ordinal Data Copyright:
 *  Copyright (c) 2000 Company: The University of Michigan
 *
 *@author     Vadim Kutsyy
 *@created    January 8, 2001
 *@version    1.0
 */
/**
 *  Title: com.kutsyy.lvspod Description: This library is set of rutins for
 *  Laten Variable model for Spatialy dependent Ordinal Data Copyright:
 *  Copyright (c) 2000 Company: The University of Michigan
 *
 *  Title: com.kutsyy.lvspod Description: This library is set of rutins for
 *  Laten Variable model for Spatialy dependent Ordinal Data Copyright:
 *  Copyright (c) 2000 Company: The University of Michigan
 *
 *@author     Vadim Kutsyy
 *@author     Vadim Kutsyy
 *@created    January 8, 2001
 *@version    1.0
 */
public class OrdinalSpatialDataAgumentationPlot extends JFrame {
	/**
	 *  Description of the Field
	 */
	public boolean stopDataAugmentation = false;
	JPanel contentPane;
	JMenuBar jMenuBar = new JMenuBar();
	JMenu jMenuFile = new JMenu();
	JMenuItem jMenuFileExit = new JMenuItem();
	JMenu jMenuHelp = new JMenu();
	JMenuItem jMenuHelpAbout = new JMenuItem();
	JLabel statusBar = new JLabel();
	BorderLayout borderLayoutMain = new BorderLayout();
	JPanel jPanelPlot = new JPanel();
	GridBagLayout gridBagLayout1 = new GridBagLayout();
	JButton jButtonStop = new JButton();
	JPanel jPanelStopBar = new JPanel();
	JProgressBar jProgressBarX = new JProgressBar();
	JProgressBar jProgressBarPhi = new JProgressBar();
	JPanel jPanelBar = new JPanel();
	CardLayout cardLayout1 = new CardLayout();
	JSplitPane jSplitPaneButton = new JSplitPane();
	JButton jButtonQuit = new JButton();
	JMenuItem jMenuFileStop = new JMenuItem();
	double q1Phi, q2Phi, q3Phi, meanPhi, xScale, yScale;
	double[] Phi;
	int current = 0, total;
	JPanel jPanelText = new JPanel();
	JTextArea jTextAreaMean = new JTextArea();
	JTextArea jTextAreaQ1 = new JTextArea();
	JTextArea jTextAreaQ2 = new JTextArea();
	JTextArea jTextAreaQ3 = new JTextArea();


	/**
	 *  Construct the frame
	 */
	public OrdinalSpatialDataAgumentationPlot() {
		enableEvents(AWTEvent.WINDOW_EVENT_MASK);
		try {
			jbInit();
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		validate();
		setVisible(true);
	}


	/**
	 *  put your documentation comment here
	 *
	 *@param  phi
	 *@param  m
	 */
	public void initialize(double[] phi, int m) {
		try {
			meanPhi = Stat.mean(phi);
			xScale = jPanelPlot.getWidth() / (double) m;
			yScale = jPanelPlot.getHeight() / 2.0;
			total = m;
			current = 0;
			statusBar.setText("0/" + total);
			jPanelPlot.repaint();
			jPanelPlot.setBackground(Color.yellow);
			jProgressBarPhi.setMaximum(m);
			jProgressBarPhi.setValue(0);
			jProgressBarX.setMaximum(m);
			jProgressBarX.setValue(0);
			stopDataAugmentation = false;
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}


	/**
	 *  put your documentation comment here
	 *
	 *@param  phi
	 */
	public void update(double[] phi) {
		//jPanelPlot.getGraphics().setColor(Color.black);
		//jPanelPlot.getGraphics().drawLine((int)((current++)*xScale),(int)((1-meanPhi)*yScale),(int)(current*xScale),(int)((1-(meanPhi=Stat.mean(phi)))*yScale));
		meanPhi = Stat.mean(phi);
		statusBar.setText(current + "/" + total);
		jTextAreaMean.setText("Mean=" + Math.round(meanPhi * 1000) / 1000.0);
		jProgressBarPhi.setValue(0);
		jProgressBarX.setValue(0);
		Phi = (double[]) phi;
		java.util.Arrays.sort(Phi);
		//jPanelPlot.getGraphics().setColor(Color.blue);
		if (current == total) {
			current = 0;
			jPanelPlot.repaint();
			jPanelPlot.setBackground(Color.yellow);
		}
		jPanelPlot.getGraphics().drawLine((int) ((current) * xScale), (int) ((1
				 - q2Phi) * yScale), (int) ((current + 1) * xScale), (int) ((1 - (q2Phi = Phi[(int) (2 * Phi.length / 4.0)])) * yScale));
		//jPanelPlot.getGraphics().setColor(Color.red);
		jPanelPlot.getGraphics().drawLine((int) ((current) * xScale), (int) ((1
				 - q1Phi) * yScale), (int) ((current + 1) * xScale), (int) ((1 - (q1Phi = Phi[(int) (Phi.length / 4.0)])) * yScale));
		jPanelPlot.getGraphics().drawLine((int) ((current) * xScale), (int) ((1
				 - q3Phi) * yScale), (int) ((current + 1) * xScale), (int) ((1 - (q3Phi = Phi[(int) (3 * Phi.length / 4.0)])) * yScale));
		jTextAreaQ2.setText("  Q2=" + Math.round(q2Phi * 1000) / 1000.0);
		jTextAreaQ1.setText("  Q1=" + Math.round(q1Phi * 1000) / 1000.0);
		jTextAreaQ3.setText("  Q3=" + Math.round(q3Phi * 1000) / 1000.0);
		current++;
	}


	/**
	 *  File | Exit action performed
	 *
	 *@param  e  Description of Parameter
	 */
	public void jMenuFileExit_actionPerformed(ActionEvent e) {
		System.exit(0);
	}


	/**
	 *  Help | About action performed
	 *
	 *@param  e  Description of Parameter
	 */
	public void jMenuHelpAbout_actionPerformed(ActionEvent e) {
		OrdinalSpatialDataAgumentationPlot_AboutBox dlg = new OrdinalSpatialDataAgumentationPlot_AboutBox(this);
		Dimension dlgSize = dlg.getPreferredSize();
		Dimension frmSize = getSize();
		Point loc = getLocation();
		dlg.setLocation((frmSize.width - dlgSize.width) / 2 + loc.x, (frmSize.height
				 - dlgSize.height) / 2 + loc.y);
		dlg.setModal(true);
		dlg.show();
	}


	/**
	 *  Overridden so we can exit when window is closed
	 *
	 *@param  e  Description of Parameter
	 */
	protected void processWindowEvent(WindowEvent e) {
		super.processWindowEvent(e);
		if (e.getID() == WindowEvent.WINDOW_CLOSING) {
			jMenuFileExit_actionPerformed(null);
		}
	}


	/**
	 *  Component initialization
	 *
	 *@exception  Exception  Description of Exception
	 */
	private void jbInit() throws Exception {
		//setIconImage(Toolkit.getDefaultToolkit().createImage(OrdinalSpatialDataAgumentationPlot.class.getResource("[Your Icon]")));
		contentPane = (JPanel) this.getContentPane();
		contentPane.setLayout(borderLayoutMain);
		this.setSize(new Dimension(500, 250));
		this.setTitle("Data Agumentation Plot");
		statusBar.setBackground(Color.lightGray);
		statusBar.setOpaque(true);
		statusBar.setToolTipText("Status bar display current information");
		statusBar.setText(" Status");
		jMenuFile.setText("File");
		jMenuFileExit.setText("Exit");
		jMenuFileExit.addActionListener(
			new ActionListener() {

				/**
				 *  put your documentation comment here
				 *
				 *@param  e
				 */
				public void actionPerformed(ActionEvent e) {
					System.exit(0);
				}
			});
		jMenuHelp.setText("Help");
		jMenuHelpAbout.setText("About");
		jMenuHelpAbout.addActionListener(
			new ActionListener() {

				/**
				 *  put your documentation comment here
				 *
				 *@param  e
				 */
				public void actionPerformed(ActionEvent e) {
					jMenuHelpAbout_actionPerformed(e);
				}
			});
		contentPane.setBackground(Color.white);
		//contentPane.setMinimumSize(new Dimension(700, 500));
		//    jPanelPlot.setBackground(Color.yellow);
		//    jPanelPlot.setBorder(BorderFactory.createLineBorder(Color.black));
		//    jPanelPlot.setSize(600,200);
		//    //jPanelPlot.setLocation(25,125);
		//jMeanPhi.setText("Mean(Phi)=0");
		jButtonStop.removeAll();
		jButtonStop.setBackground(Color.red);
		jButtonStop.setFont(new java.awt.Font("Dialog", 1, 14));
		jButtonStop.setBorder(BorderFactory.createRaisedBevelBorder());
		jButtonStop.setToolTipText("To stop Data Augmentation Click Here");
		jButtonStop.setHorizontalTextPosition(SwingConstants.CENTER);
		jButtonStop.setText("Stop");
		jButtonStop.setVerticalAlignment(SwingConstants.BOTTOM);
		jButtonStop.addMouseListener(
			new java.awt.event.MouseAdapter() {

				/**
				 *  put your documentation comment here
				 *
				 *@param  e
				 */
				public void mouseClicked(MouseEvent e) {
					stopDataAugmentation = true;
				}
			});
		jButtonQuit.addMouseListener(
			new java.awt.event.MouseAdapter() {

				/**
				 *  put your documentation comment here
				 *
				 *@param  e
				 */
				public void mouseClicked(MouseEvent e) {
					System.exit(0);
				}
			});
		jPanelStopBar.setBackground(Color.white);
		jPanelStopBar.setDoubleBuffered(false);
		jPanelStopBar.setMinimumSize(new Dimension(10, 30));
		jPanelStopBar.setLayout(new BorderLayout());
		jProgressBarX.setBackground(Color.cyan);
		jProgressBarX.setBorder(BorderFactory.createLineBorder(Color.black));
		jProgressBarX.setMinimumSize(new Dimension(10, 5));
		jProgressBarX.setToolTipText("");
		jProgressBarX.setMaximum(10000);
		jProgressBarPhi.setToolTipText("");
		jProgressBarPhi.setMaximum(10000);
		jProgressBarPhi.setMinimumSize(new Dimension(10, 5));
		jProgressBarPhi.setBorder(BorderFactory.createLineBorder(Color.black));
		jProgressBarPhi.setBackground(Color.cyan);
		jPanelBar.setBackground(Color.white);
		jPanelBar.setLayout(new GridLayout(2, 1));
		jSplitPaneButton.setOrientation(JSplitPane.VERTICAL_SPLIT);
		jSplitPaneButton.setBackground(Color.white);
		jSplitPaneButton.setForeground(Color.white);
		jSplitPaneButton.setBorder(null);
		jSplitPaneButton.setDividerSize(0);
		jSplitPaneButton.setLeftComponent(null);
		jSplitPaneButton.setTopComponent(jButtonStop);
		jButtonQuit.setBackground(Color.red);
		jButtonQuit.setFont(new java.awt.Font("Dialog", 1, 16));
		jButtonQuit.setBorder(BorderFactory.createRaisedBevelBorder());
		jButtonQuit.setText("Quit");
		jMenuFileStop.addActionListener(
			new ActionListener() {

				/**
				 *  put your documentation comment here
				 *
				 *@param  e
				 */
				public void actionPerformed(ActionEvent e) {
					stopDataAugmentation = true;
				}
			});
		jMenuFileStop.setActionCommand("Stop");
		jMenuFileStop.setText("Stop");
		//    jQ1Phi.setText("Q1(phi)=0");
		//    jQ3Phi.setText("Q3(phi)=");
		//    jQ2Phi.setText("Q2(phi)");
		//    jQ2Phi.setToolTipText("");
		jPanelText.setBackground(Color.white);
		jTextAreaMean.setText("Mean=0.000");
		jTextAreaQ1.setText("Q1=0.000");
		jTextAreaQ2.setText("Q2=0.000");
		jTextAreaQ3.setText("Q3=0.000");
		jPanelPlot.setBackground(Color.yellow);
		jMenuFile.add(jMenuFileStop);
		jMenuFile.add(jMenuFileExit);
		jMenuHelp.add(jMenuHelpAbout);
		jMenuBar.add(jMenuFile);
		jMenuBar.add(jMenuHelp);
		this.setJMenuBar(jMenuBar);
		contentPane.add(statusBar, BorderLayout.SOUTH);
		contentPane.add(jPanelPlot, BorderLayout.CENTER);
		contentPane.add(jPanelStopBar, BorderLayout.NORTH);
		jPanelStopBar.add(jPanelBar, BorderLayout.WEST);
		jPanelBar.add(jProgressBarX, BorderLayout.WEST);
		jPanelBar.add(jProgressBarPhi, BorderLayout.CENTER);
		jPanelStopBar.add(jSplitPaneButton, BorderLayout.EAST);
		jSplitPaneButton.add(jButtonStop, JSplitPane.TOP);
		jSplitPaneButton.add(jButtonQuit, JSplitPane.BOTTOM);
		jPanelText.setLayout(new GridLayout(2, 2));
		jPanelStopBar.add(jPanelText, BorderLayout.CENTER);
		jPanelText.add(jTextAreaMean, null);
		jPanelText.add(jTextAreaQ1, null);
		jPanelText.add(jTextAreaQ2, null);
		jPanelText.add(jTextAreaQ3, null);
	}
	//  protected void stop(){
	//  stopDataAugmentation=true;
	//   }
	//  protected void quit(){
	//  if(stopDataAugmentation) System.exit(0);
	//  }
	//  void jMenuFileStop_actionPerformed(ActionEvent e) {
	//  stopDataAugmentation=true;
	//
	//  }
}


