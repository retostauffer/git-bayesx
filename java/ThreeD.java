//This program is to create and run the applet for surfplot program

import java.applet.Applet;
import java.awt.*;
import java.io.*;
import java.net.URL;

public class ThreeD extends Applet implements Runnable
{
    Model3D md;
    boolean painted = true;
    double xfac;
    int prevx, prevy;
    double xtheta, ytheta;
    Matrix3D amat = new Matrix3D(), tmat = new Matrix3D();
    String mdname = null;
    String message = null;

    public void init()
    {
	mdname = getParameter("model");
	amat.yrot(20);
	amat.xrot(110);
	if (mdname == null)
	    mdname = "model.txt";
	resize(getSize().width <= 20 ? 400 : getSize().width,
	       getSize().height <= 20 ? 400 : getSize().height);
    }

    public void run()
    {
	InputStream is = null;
	try
	{
	    Thread.currentThread().setPriority(Thread.MIN_PRIORITY);
	    is = new URL(getDocumentBase(), mdname).openStream();
	    Model3D m = new Model3D (is);
	    md = m;
	    m.findBB();
	    m.compress();
	    double xw = m.xmax - m.xmin;
	    double yw = m.ymax - m.ymin;
	    double zw = m.zmax - m.zmin;
	    if (yw > xw)
		xw = yw;
	    if (zw > xw)
		xw = zw;
	    double f1 = getSize().width / xw;
	    double f2 = getSize().height / xw;
	    xfac = 0.7 * (f1 < f2 ? f1 : f2);
	}
	catch(Exception e)
	{
	    md = null;
	    message = e.toString();
	}
	try
	{
	    if (is != null)
		is.close();
	}
	catch(Exception e)
	{
	}
	repaint();
    }

    public void start()
    {
	 if (md == null && message == null)
	    new Thread(this).start();
    }

    public void stop()
    {
    }

    // For rotating the object by clicking on the applet
    public boolean mouseDown(Event e, int x, int y)
    {
	 prevx = x;
	 prevy = y;
	 return true;
    }

    public boolean mouseDrag(Event e, int x, int y)
    {
	 tmat.unit();
	 double xtheta = (prevy - y) * 360 / getSize().width;
	 double ytheta = (x - prevx) * 360 / getSize().height;
	 tmat.xrot(xtheta);
	 tmat.yrot(ytheta);
	 amat.mult(tmat);
	 if (painted)
	 {
	    painted = false;
	    repaint();
	 }
	 prevx = x;
	 prevy = y;
	 return true;
    }

    // Creating an object from Model3D class and painting it on the applet
    public void paint(Graphics g2)
    {
     Graphics2D g = (Graphics2D) g2;
	 if (md != null)
	 {
	    md.mat.unit();
	    md.mat.translate(-(md.xmin + md.xmax) / 2,
			     -(md.ymin + md.ymax) / 2,
			     -(md.zmin + md.zmax) / 2);
	    md.mat.mult(amat);
	    md.mat.scale(xfac, -xfac, 16 * xfac / getSize().width);
	    md.mat.translate(getSize().width / 2, getSize().height / 2, 8);
	    md.transformed = false;
	    md.paint(g);
	    setPainted();

	 }
	 else if (message != null)
	 {
	    g.drawString("Error in model:", 3, 20);
	    g.drawString(message, 10, 40);
	 }
	}

    private synchronized void setPainted()
    {
	 painted = true;
	 notifyAll();
    }
}
