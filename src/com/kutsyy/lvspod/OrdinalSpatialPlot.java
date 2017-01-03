/*
 * put your module comment here
 * formatted with JxBeauty (c) johann.langhofer@nextra.at
 */

package com.kutsyy.lvspod;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import com.kutsyy.util.*;

/**
 *  put your documentation comment here
 *
 *@author     Vadim Kutsyy
 *@created    January 8, 2001
 */
public class OrdinalSpatialPlot extends JFrame {
	JPanel contentPane;
	JMenuBar jMenuBar1 = new JMenuBar();
	JMenu jMenuFile = new JMenu();
	JMenuItem jMenuFileExit = new JMenuItem();
	JMenu jMenuHelp = new JMenu();
	JMenuItem jMenuHelpAbout = new JMenuItem();
	GridLayout gridLayout1 = new GridLayout();
	jPanelSpatialDataPlot[] jPanelPlot;
	// = new JPanel();
	//JPanel jPanelPlot = new JPanel();
	int n = 10;
	int[] centerX, centerY, size;
	int activePanel = 0;
	Graphics[] graphics;


	/**
	 *  Construct the frame
	 */
	public OrdinalSpatialPlot() {
		this(1);
	}


	/**
	 *  put your documentation comment here
	 *
	 *@param  N    Description of Parameter
	 */
	public OrdinalSpatialPlot(int N) {
		n = (int) Math.ceil(Math.sqrt(N));
		n *= n;
		try {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		enableEvents(AWTEvent.WINDOW_EVENT_MASK);
		try {
			jbInit();
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		validate();
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		Dimension frameSize = getSize();
		setLocation((screenSize.width - frameSize.width), (screenSize.height
				 - frameSize.height) - 50);
		setVisible(true);

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
		OrdinalSpatialPlot_AboutBox dlg = new OrdinalSpatialPlot_AboutBox(this);
		Dimension dlgSize = dlg.getPreferredSize();
		Dimension frmSize = getSize();
		Point loc = getLocation();
		dlg.setLocation((frmSize.width - dlgSize.width) / 2 + loc.x, (frmSize.height
				 - dlgSize.height) / 2 + loc.y);
		dlg.setModal(true);
		dlg.show();
	}


	/**
	 *  put your documentation comment here
	 *
	 *@param  T          Description of Parameter
	 */
	public void add(OrdinalSpatial T) {
		add(T.Y, T.xLoc, T.yLoc);
	}


	/**
	 *  Description of the Method
	 *
	 *@param  Y          Description of Parameter
	 *@param  xLocation  Description of Parameter
	 *@param  yLocation  Description of Parameter
	 */
	public void add(int[] Y, int[] xLocation, int[] yLocation) {
		if (activePanel == jPanelPlot.length) {
			activePanel = 0;
			for (int i = 0; i < jPanelPlot.length; i++) {
				jPanelPlot[i].data = false;
				jPanelPlot[i].repaint();
			}
		}
		jPanelPlot[activePanel].add(Y, xLocation, yLocation);
		activePanel++;
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
		//setIconImage(Toolkit.getDefaultToolkit().createImage(OrdinalSpatialPlot.class.getResource("[Your Icon]")));
		contentPane = (JPanel) this.getContentPane();
		contentPane.setLayout(new GridLayout((int) Math.sqrt(n), (int) Math.sqrt(n)));
		//contentPane.setLayout(new FlowLayout());
		this.setSize(new Dimension(400, 300));
		this.setTitle("Plot of Ordinal Spatial Data");
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
					jMenuFileExit_actionPerformed(e);
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
		contentPane.setEnabled(true);
		contentPane.setBorder(BorderFactory.createLineBorder(Color.black));
		jMenuFile.add(jMenuFileExit);
		jMenuHelp.add(jMenuHelpAbout);
		jMenuBar1.add(jMenuFile);
		jMenuBar1.add(jMenuHelp);
		jPanelPlot = new jPanelSpatialDataPlot[n];
		graphics = new Graphics[n];
		for (int i = 0; i < n; i++) {
			jPanelPlot[i] = new jPanelSpatialDataPlot();
		}
		for (int i = 0; i < n; i++) {
			jPanelPlot[i].setBackground(Color.white);
		}
		for (int i = 0; i < n; i++) {
			contentPane.add(jPanelPlot[i], i);
		}
		//jPanelPlot.setBackground(Color.white);
		//contentPane.add(jPanelPlot,i);
		//contentPane.add(jPanelPlot);
		this.setJMenuBar(jMenuBar1);
	}
}


