package com.kutsyy.lvspod;

import javax.swing.JPanel;
import java.awt.Color;
import java.awt.*;
import java.util.*;
import java.awt.event.*;
import com.kutsyy.util.*;

/**
 *  panel of spatial data Copyright: Copyright (c) 2000
 *
 *@author     Vadim Kutsyy
 *@created    January 14, 2001
 *@version    1.0
 */
public class jPanelSpatialDataPlot extends JPanel {

    /**
     *  observations
     */
    private int[] Y;
    /**
     *  X location
     */
    private int[] xLocation;
    /**
     *  y Location
     */
    private int[] yLocation;

    /**
     *  if plot shoudl be created or not
     */
    public boolean data = false;


    /**
     *  creates panel
     */
    public jPanelSpatialDataPlot() {
        try {
            jbInit();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }


    /**
     *  make plot
     *
     *@param  g graphics
     */
    public void paint(Graphics g) {
        super.paint(g);
        if (data) {
            int xLim = MyMath.max(xLocation);
            int yLim = MyMath.max(yLocation);
            int xSize = (int) (this.getSize().width / (xLim + 2));
            int ySize = (int) (this.getSize().height / (yLim + 2));
            int xloc;
            int yloc;
            g.setColor(Color.black);
            for (int i = 0; i < Y.length; i++) {
                g.drawRect((xLocation[i] + 1) * xSize, (yLocation[i] + 1) * ySize,
                        xSize, ySize);
                switch (Y[i]) {
                    case 0:
                        g.setColor(Color.white);
                        break;
                    case 1:
                        g.setColor(Color.black);
                        break;
                    case 2:
                        g.setColor(Color.yellow);
                        break;
                }
                g.fillRect((xLocation[i] + 1) * xSize, (yLocation[i] + 1) * ySize,
                        xSize, ySize);
                g.setColor(Color.black);
                g.drawRect((xLocation[i] + 1) * xSize, (yLocation[i] + 1) * ySize,
                        xSize, ySize);
            }
        }
    }


    /**
     *  make plot of the data
     *
     *@param  Y observations
     *@param  xLocation location (x axes)
     *@param  yLocation  location (y axes)
     */
    public void add(int[] Y, int[] xLocation, int yLocation[]) {
        data = true;
        this.Y = Y;
        this.xLocation = xLocation;
        this.yLocation = yLocation;
        repaint();
    }


    /**
     *  called by constructor
     *
     */
    private void jbInit() throws Exception {
        this.setBackground(Color.white);
        this.setEnabled(true);
    }
}

