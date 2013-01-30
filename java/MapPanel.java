import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.util.*;
import java.io.*;


//import gov.sandia.postscript.PSGrBase;
//import gov.sandia.postscript.PSGr1;
//import gov.sandia.postscript.PSGr2;

/**
 *MapPanel.java
 *
 *Created on 12.04.2001
 *
 *Last modified on 17.04.2001
 */


public class MapPanel extends JPanel
{

private BayesX b;

public static final double MAXDOUBLE = 1.7976931348623157E308;

public static final double WHITE_X = 95.047;
public static final double WHITE_Y = 100.000;
public static final double WHITE_Z = 108.883;
public static final double WHITE_u = 0.1978398;
public static final double WHITE_v = 0.4683363;

/* Gamma-Korrektur Faktor fuer hcl Farben */

public static final double GAMMA = 2.2;

/* PI/180.0 */

public static final double DEG2RAD = 0.0174532925199433;

/* Default-Werte fuer hcl Farben */

private double h1 = 130;
private double h2 =  25;
private double  c = 100;
private double l1 =  90;
private double l2 =  70;			

protected short page;

private double xstep;
private double ystep;
private double xstart;
private double ystart;

private int PAGEHEIGHT = 842;
private int PAGEWIDTH = 596;

private int xoffset = 120;
private int yoffset = 120;
private int width = 356;
private int height = 210;

private int scale = 1;
private int pointsize = 3;
private int pspointsize = 20;
private int fontsize = 12;

private double titlescale = 1.5;

private double minX;
private double maxX;
private double minY;
private double maxY;

private PSGr2 PSGr;

public MapPanel(BayesX b)
	{
	super();
	this.b=b;
	}

public void setplotparam(int xoff, int yoff, int w, int h, int ps, int fs, double ts)
        {
        xoffset = xoff;
        yoffset = yoff;
        width = w;
        height = h;
	    pspointsize = ps;
	    fontsize = fs;
	    titlescale = ts;
        }

public void paintComponent(Graphics g1)                 // fuer Bildschirmanzeige
	{

        Graphics2D g2 = (Graphics2D)g1;

	g2.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC));
        g2.setPaint(Color.white);

        Rectangle R = this.getBounds();

        R.setBounds(0,0,(int)R.getWidth(),(int)R.getHeight());
        g2.fill(R);

        g2.setColor(Color.black);

        switch(b.function){
            case 1:
                showmap(g2);    // Umrisse zeichnen (map.describe)
                break;
            case 2:
                drawmap(g1);    // Funktion drawmap
                break;
            case 3:
                plotnonp(g1);
                break;
            case 4:
                plotsample(g1);
                break;
            case 5:
                plotautocor(g1);
                break;
	   case 6:
	        plotsurf(g1);    // For plotting a surface
	        break;
            }
	}

private void showmap(Graphics2D g)
        {

        double width = getWidth();
	double height = getHeight();
	double[] d2 = new double[4];
        double[] centroid = new double[2];

        b.getboundaries(d2);
        maxX = d2[2];
	maxY = d2[3];
	minX = d2[0];
	minY = d2[1];

        setSize((int)width,(int)height);
        height = (int)(width*(maxY-minY)/(maxX-minX));

        for(int i=0; i<b.getnrregions(); i++)
		{
        	for(int j=0; j<b.getnrpoly(i); j++)
	        	{
		        for(int k=0; k<b.getnrlines(i,j); k++)
			        {
       				b.getline(d2,i,j,k);
       				g.draw(new Line2D.Double((d2[0]-minX)*0.9*width/(maxX-minX)+0.05*width,-(d2[2]-minY)*0.9*height/(maxY-minY)+0.95*height,
       					                 (d2[1]-minX)*0.9*width/(maxX-minX)+0.05*width,-(d2[3]-minY)*0.9*height/(maxY-minY)+0.95*height));
			        }
		        }
                }

        if(b.drawnames)
                {
                g.setFont(new Font("TimesRoman", Font.BOLD, b.fontsize));
                for(int i=0; i<b.getnrregions(); i++)
                        {
                        b.getcentroid(centroid,i);
                        g.drawString(b.getregionname(i),
                                    (int)((centroid[0]-minX)*0.9*width/(maxX-minX)+0.05*width),
                                    (int)(-(centroid[1]-minY)*0.9*height/(maxY-minY)+0.95*height)
                                    );
                        }
                g.setFont(new Font("TimesRoman", Font.PLAIN, b.fontsize));
                }
        }

private void drawmap(Graphics g)
        {

        double width = getWidth();
	double height = getHeight();
	double[] d2 = new double[4];

        b.getboundaries(d2);
	maxX = d2[2];
	maxY = d2[3];
	minX = d2[0];
	minY = d2[1];

        int offset = 0;
        if(!b.title.equals(""))
           offset = offset + 10 + (int)(b.titlescale*b.fontsize);

        helpdrawmap(g,height,width,offset,true);

        }

public void SaveMap(PrintWriter out)              // Zum Speichern als PostScript
	{

	double[] d2 = new double[4];
        double[] centroid = new double[2];

        b.getboundaries(d2);
	maxX = d2[2];
	maxY = d2[3];
	minX = d2[0];
	minY = d2[1];

	double width = PAGEWIDTH;
	double height = width*(maxY-minY)/(maxX-minX);

        if(height>PAGEHEIGHT)
                {
                width = width*PAGEHEIGHT/height;
                height = height*PAGEHEIGHT/height;
                }

        out.println("%!PS-Adobe-3.0");
        out.print("%%BoundingBox:");
        out.print((int)Math.round(0.05*width));
        out.print(" ");
        out.print(PAGEHEIGHT-(int)Math.round(0.95*height+2));
        out.print(" ");
        out.print((int)Math.round(0.95*width));
        out.print(" ");
        out.println(PAGEHEIGHT-((int)Math.round(0.05*height)));
        out.println("0 setlinewidth");
        out.println("%%Pages:1");
        out.println("%%Page:1 1");

        Graphics g = new PSGr2(out);

        for(int i=0; i<b.getnrregions(); i++)
		{
        	for(int j=0; j<b.getnrpoly(i); j++)
	        	{
		        for(int k=0; k<b.getnrlines(i,j); k++)
			        {
       				b.getline(d2,i,j,k);
       				g.drawLine((int)Math.round((d2[0]-minX)*0.9*width/(maxX-minX)+0.05*width),(int)Math.round(-(d2[2]-minY)*0.9*height/(maxY-minY)+0.95*height),
       					   (int)Math.round((d2[1]-minX)*0.9*width/(maxX-minX)+0.05*width),(int)Math.round(-(d2[3]-minY)*0.9*height/(maxY-minY)+0.95*height));
			        }
		        }
                }


        if(b.drawnames)
                {
                g.setFont(new Font("TimesRoman", Font.BOLD, b.fontsize));
                for(int i=0; i<b.getnrregions(); i++)
                        {
                        b.getcentroid(centroid,i);
                        g.drawString(b.getregionname(i),
                                     (int)((centroid[0]-minX)*0.9*width/(maxX-minX)+0.05*width),
                                     (int)(-(centroid[1]-minY)*0.9*height/(maxY-minY)+0.95*height)
                                    );
                        }
                g.setFont(new Font("TimesRoman", Font.PLAIN, b.fontsize));
                }

        out.println("showpage");
        out.println("grestore");
        out.println("gsave");
	out.close();

	}

public void Savedrawmap(PrintWriter out)              // Zum Speichern als PostScript
	{

	double[] d2 = new double[4];

        b.getboundaries(d2);
	maxX = d2[2];
	maxY = d2[3];
	minX = d2[0];
	minY = d2[1];

	double width = PAGEWIDTH;
	double height = width*(maxY-minY)/(maxX-minX);

        int offset = 0;
        if(!b.title.equals(""))
            offset = offset + 10 + (int)(b.titlescale*b.fontsize);

        int hlp = offset;

        if(b.legend)
                hlp = hlp + 30 + 13*b.fontsize/10;

        if(height+hlp>PAGEHEIGHT)
                {
                width = width*(PAGEHEIGHT-hlp)/(height+hlp);
                height = height*(PAGEHEIGHT-hlp)/(height+hlp);
                }

        out.println("%!PS-Adobe-3.0");
        out.print("%%BoundingBox:");
        out.print((int)Math.round(0.05*width));
        out.print(" ");
        out.print(PAGEHEIGHT-(int)Math.round(0.95*height+hlp+2));
        out.print(" ");
        out.print((int)Math.round(0.95*width));
        out.print(" ");
        out.println(PAGEHEIGHT-((int)Math.round(0.05*height)));
        out.println("0 setlinewidth");
        out.println("%%Pages:1");
        out.println("%%Page:1 1");

        Graphics g = new PSGr2(out);

        helpdrawmap(g,height,width,offset,false);

        out.println("showpage");
        out.println("grestore");
        out.println("gsave");
	out.close();

	}

private void helpdrawmap(Graphics g, double height, double width, int offset, boolean centering)
        {

        boolean NA = false;

        b.nrNA = 0;

	double[] d2 = new double[4];
        double[] help = new double[4];
        double[] centroid = new double[2];

        Polygon p;

        double x;
        double y;

        setSize((int)width,(int)height);
        height = (int)(width*(maxY-minY)/(maxX-minX));

        for(int i=0; i<b.getnrregions(); i++)
		{
                if(!b.isin(i))
                        {
	        	for(int j=0; j<b.getnrpoly(i); j++)
        			{
// Polygon einlesen
                                p = new Polygon();
                                b.getline(help,i,j,b.getnrlines(i,j)-1);
                                b.getline(d2,i,j,0);
                                if( (d2[0]==help[0] && d2[2]==help[2]) || (d2[0]==help[1] && d2[2]==help[3]))
                                        {
                                        p.addPoint( (int)Math.round((d2[0]-minX)*0.9*width/(maxX-minX)+0.05*width),
                                                    (int)Math.round(-(d2[2]-minY)*0.9*height/(maxY-minY)+0.95*height + offset)
                                                  );
                                        p.addPoint( (int)Math.round((d2[1]-minX)*0.9*width/(maxX-minX)+0.05*width),
                                                    (int)Math.round(-(d2[3]-minY)*0.9*height/(maxY-minY)+0.95*height + offset)
                                                  );
                                        x = d2[1];
                                        y = d2[3];
                                        }
                                else
                                        {
                                        p.addPoint( (int)Math.round((d2[1]-minX)*0.9*width/(maxX-minX)+0.05*width),
                                                    (int)Math.round(-(d2[3]-minY)*0.9*height/(maxY-minY)+0.95*height + offset)
                                                  );
                                        p.addPoint( (int)Math.round((d2[0]-minX)*0.9*width/(maxX-minX)+0.05*width),
                                                    (int)Math.round(-(d2[2]-minY)*0.9*height/(maxY-minY)+0.95*height + offset)
                                                  );
                                        x = d2[0];
                                        y = d2[2];
                                        }
                                for(int k=1; k<b.getnrlines(i,j); k++)
        				{
                                        b.getline(d2,i,j,k);
                                        if( d2[0]==x && d2[2]==y )
                                                {
                                                x = d2[1];
                                                y = d2[3];
                                                p.addPoint( (int)Math.round((d2[1]-minX)*0.9*width/(maxX-minX)+0.05*width),
                                                            (int)Math.round(-(d2[3]-minY)*0.9*height/(maxY-minY)+0.95*height + offset)
                                                          );
                                                }
                                        else
                                                {
                                                x = d2[0];
                                                y = d2[2];
                                                p.addPoint( (int)Math.round((d2[0]-minX)*0.9*width/(maxX-minX)+0.05*width),
                                                            (int)Math.round(-(d2[2]-minY)*0.9*height/(maxY-minY)+0.95*height + offset)
                                                          );
                                                }
                                        }
// ENDE: Polygon einlesen
                                if( b.shades > 0 )
                                	{
	                                NA = ComputeColor(g,i);
        	                        g.fillPolygon(p.xpoints,p.ypoints,b.getnrlines(i,j)+1);
                	                g.setColor(Color.black);
                        	        }
                                g.drawPolygon(p.xpoints,p.ypoints,b.getnrlines(i,j)+1);

                                if(NA)
                                        {
                                        drawNA(g,p);
                                        b.nrNA = b.nrNA + 1;
                                        }

                                }   // ENDE: for(int i=0; i<b.getnrpoly(); i++)
                        }   // ENDE: if(b.isnotin())
        	}   // ENDE: for(int i=0; i<b.getnrregions(); i++)

        for(int i=0; i<b.getnrregions(); i++)
		{
                if(b.isin(i))
                        {
	        	for(int j=0; j<b.getnrpoly(i); j++)
        			{
// Polygon einlesen
                                p = new Polygon();
                                b.getline(help,i,j,b.getnrlines(i,j)-1);
                                b.getline(d2,i,j,0);
                                if( (d2[0]==help[0] && d2[2]==help[2]) || (d2[0]==help[1] && d2[2]==help[3]))
                                        {
                                        p.addPoint( (int)Math.round((d2[0]-minX)*0.9*width/(maxX-minX)+0.05*width),
                                                    (int)Math.round(-(d2[2]-minY)*0.9*height/(maxY-minY)+0.95*height + offset)
                                                  );
                                        p.addPoint( (int)Math.round((d2[1]-minX)*0.9*width/(maxX-minX)+0.05*width),
                                                    (int)Math.round(-(d2[3]-minY)*0.9*height/(maxY-minY)+0.95*height + offset)
                                                  );
                                        x = d2[1];
                                        y = d2[3];
                                        }
                                else
                                        {
                                        p.addPoint( (int)Math.round((d2[1]-minX)*0.9*width/(maxX-minX)+0.05*width),
                                                    (int)Math.round(-(d2[3]-minY)*0.9*height/(maxY-minY)+0.95*height + offset)
                                                  );
                                        p.addPoint( (int)Math.round((d2[0]-minX)*0.9*width/(maxX-minX)+0.05*width),
                                                    (int)Math.round(-(d2[2]-minY)*0.9*height/(maxY-minY)+0.95*height + offset)
                                                  );
                                        x = d2[0];
                                        y = d2[2];
                                        }
                                for(int k=1; k<b.getnrlines(i,j); k++)
        				{
                                        b.getline(d2,i,j,k);
                                        if( d2[0]==x && d2[2]==y )
                                                {
                                                x = d2[1];
                                                y = d2[3];
                                                p.addPoint( (int)Math.round((d2[1]-minX)*0.9*width/(maxX-minX)+0.05*width),
                                                            (int)Math.round(-(d2[3]-minY)*0.9*height/(maxY-minY)+0.95*height + offset)
                                                          );
                                                }
                                        else
                                                {
                                                x = d2[0];
                                                y = d2[2];
                                                p.addPoint( (int)Math.round((d2[0]-minX)*0.9*width/(maxX-minX)+0.05*width),
                                                            (int)Math.round(-(d2[2]-minY)*0.9*height/(maxY-minY)+0.95*height + offset)
                                                          );
                                                }
                                        }
// ENDE: Polygon einlesen
                                if(b.shades > 0)
                                	{
	                                NA = ComputeColor(g,i);
        	                        g.fillPolygon(p.xpoints,p.ypoints,b.getnrlines(i,j)+1);
                	                g.setColor(Color.black);
                        	        }
                                g.drawPolygon(p.xpoints,p.ypoints,b.getnrlines(i,j)+1);

                                if(NA)
                                        {
                                        drawNA(g,p);
                                        b.nrNA = b.nrNA + 1;
                                        }

                                }   // ENDE: for(int i=0; i<b.getnrpoly(); i++)
                        }   // ENDE: if(b.isnotin())
        	}   // ENDE: for(int i=0; i<b.getnrregions(); i++)

                if(b.drawnames)
                        {
                        g.setFont(new Font("TimesRoman", Font.BOLD, b.fontsize));
                        for(int i=0; i<b.getnrregions(); i++)
                                {
                                b.getcentroid(centroid,i);
                                g.drawString(b.getregionname(i),
                                             (int)((centroid[0]-minX)*0.9*width/(maxX-minX)+0.05*width),
                                             (int)(-(centroid[1]-minY)*0.9*height/(maxY-minY)+0.95*height+offset)
                                             );
                                }
                        g.setFont(new Font("TimesRoman", Font.PLAIN, b.fontsize));
                        }

                if(b.legend)
                        drawlegend(g,height,width,centering);
                if(!b.title.equals(""))
                        drawtitle(g,height,width,centering);
/*
                if(nrNA > 0 && nrNA < b.getnrregions())
                    {
                    String str = "NOTE: " + String.valueOf(nrNA) + " missing value(s) plotted\n";
                    b.Out(str,false,false,(short)11,0,0,0);
                    }
                else if (nrNA >= b.getnrregions())
                    {
                    String str = "WARNING: only missing values plotted - map probably doesn't match data file\n";
                    b.Out(str,true,true,(short)11,255,0,0);
                    }
*/
        }

private boolean ComputeColor(Graphics g, int i)
        {

        int r;

        double min = b.lowerlimit;
        double max = b.upperlimit;

        double value = 0.0;
        double mitte = min+(max-min)/2;

        int d=0;
        while(d<b.getDRows() && b.getname(i)!=b.getDoubleValue(d,1))
                d++;

        if(d<b.getDRows())
        {
        value = b.getDoubleValue(d,0);

        if(b.color == true)
                {
		if(b.hcl)
			{
			double rval;
			int[] RGB = new int[3];

                        if(b.shades == 1)
                            rval = 0.0;
                        else
                            rval = (double)(Math.floor( (b.shades-1)*(value-min)/(max-min) ))/(b.shades-1) * 2.0 - 1.0 ;

                        if(rval > 1.0)
                                rval = 1.0;
                        if(rval < -1.0)
                                rval = -1.0;

			if(b.swap)
				rval = -rval;

			if(rval > 0.0)
				RGB = hcl2rgb(h1,c*Math.abs(rval),l1+(l2-l1)*Math.abs(rval));
			else
				RGB = hcl2rgb(h2,c*Math.abs(rval),l1+(l2-l1)*Math.abs(rval));

                        g.setColor(new Color(RGB[0],RGB[1],RGB[2]));

			}
		else
			{
        	        if(value < mitte)
                	        {
	
        	                if(b.shades == 1)
                	            r = 255;
                        	else
	                            r = (int)(b.shades*(value-min)/(max-min)) * 510/(b.shades-1);

        	                if(r>255)
                	                r = 255;
                        	if(r<0)
	                                	r = 0;
	                        if(b.swap)
        	                        g.setColor(new Color(r,255,0));
                	        else
                        	        g.setColor(new Color(255,r,0));
	                        }
        	        else
                	        {
	
        	                if(b.shades == 1)
                	            r = 255;
                        	else
	                            r = ( (b.shades-1) - (int)(b.shades*(value-min)/(max-min)) ) * 510/(b.shades-1);

        	                if(r>255)
                	                r = 255;
                        	if(r<0)
                                	r = 0;
	                        if(b.swap)
        	                        g.setColor(new Color(255,r,0));
                	        else
                        	        g.setColor(new Color(r,255,0));
	                        }
        	        }
		}
        else
                {

                if(b.shades == 1)
                    r = 127;
                else
                    r = (int)(b.shades*(value-min)/(max-min)) * 255/(b.shades-1);

                if(r>255)
                        r = 255;
                if(r<0)
                        r = 0;
                if(b.swap)
                        g.setColor(new Color(255-r,255-r,255-r));
                else
                        g.setColor(new Color(r,r,r));
                }
        return false;
        }
        else
        {
        g.setColor(new Color(255,255,255));
        return true;
        }

        }

private void drawlegend(Graphics g,double height,double width,boolean centering)
        {

        int r;
        int end;
	double rval;
        double step;

        int offset = 0;
        if(!b.title.equals(""))
            offset = offset + 10 + (int)(b.titlescale*b.fontsize);

        double xoffset = 0.55*width;
        double yoffset = 0.95*height + offset + 10;
        int legendheight = 20;

        Polygon p;

        if(b.color)
                {
		if(b.hcl)
			{
			int[] RGB = new int[3];

        	        step = width*0.3/b.shades;

                	for(int i=0;i<b.shades;i++)
	                        {
        	                p = new Polygon();
                	        p.addPoint((int)Math.round(xoffset + i*step), (int)Math.round(yoffset));
                        	p.addPoint((int)Math.round(xoffset + i*step), (int)Math.round(yoffset+legendheight));
	                        p.addPoint((int)Math.round(xoffset + (i+1)*step), (int)Math.round(yoffset+legendheight));
        	                p.addPoint((int)Math.round(xoffset + (i+1)*step), (int)Math.round(yoffset));
                	        p.addPoint((int)Math.round(xoffset + i*step), (int)Math.round(yoffset));

                	        if(b.shades == 1)
                        	    rval = 0.0;
	                        else
        	                    rval = (double)(i)/(b.shades-1) * 2.0 - 1.0;

                        	if(rval > 1.0)
                	                rval = 1.0;
        	                if(rval < -1.0)
	                                rval = -1.0;

				if(b.swap)
					rval = -rval;

				if(rval > 0.0)
					RGB = hcl2rgb(h1,c*Math.abs(rval),l1+(l2-l1)*Math.abs(rval));
				else
					RGB = hcl2rgb(h2,c*Math.abs(rval),l1+(l2-l1)*Math.abs(rval));

	                        g.setColor(new Color(RGB[0],RGB[1],RGB[2]));
	                        g.fillPolygon(p.xpoints,p.ypoints,5);
				}
			}
		else
			{
	                end = (b.shades+1)/2;
        	        step = width*0.3/b.shades;

                	for(int i=0;i<end;i++)
	                        {
        	                p = new Polygon();
                	        p.addPoint((int)Math.round(xoffset + i*step), (int)Math.round(yoffset));
                        	p.addPoint((int)Math.round(xoffset + i*step), (int)Math.round(yoffset+legendheight));
	                        p.addPoint((int)Math.round(xoffset + (i+1)*step), (int)Math.round(yoffset+legendheight));
        	                p.addPoint((int)Math.round(xoffset + (i+1)*step), (int)Math.round(yoffset));
                	        p.addPoint((int)Math.round(xoffset + i*step), (int)Math.round(yoffset));

	                        if(b.shades == 1)
       		                    r = 255;
               		        else
                       		    r = i*510/(b.shades-1);

	                        if(b.swap)
       		                        g.setColor(new Color(r,255,0));
               		        else
                       		        g.setColor(new Color(255,r,0));

	                        g.fillPolygon(p.xpoints,p.ypoints,5);
        	                }
                	for(int i=0;i<end;i++)
                        	{
	                        p = new Polygon();
        	                p.addPoint((int)Math.round(xoffset + 0.3*width - i*step), (int)Math.round(yoffset));
                	        p.addPoint((int)Math.round(xoffset + 0.3*width - i*step), (int)Math.round(yoffset+legendheight));
                        	p.addPoint((int)Math.round(xoffset + 0.3*width - (i+1)*step), (int)Math.round(yoffset+legendheight));
	                        p.addPoint((int)Math.round(xoffset + 0.3*width - (i+1)*step), (int)Math.round(yoffset));
        	                p.addPoint((int)Math.round(xoffset + 0.3*width - i*step), (int)Math.round(yoffset));

                	        if(b.shades == 1)
                       		    r = 255;
	               	        else
       		                    r = i*510/(b.shades-1);
	
        	                if(b.swap)
                	       	        g.setColor(new Color(255,r,0));
              	        	else
       	                        	g.setColor(new Color(r,255,0));

	                        g.fillPolygon(p.xpoints,p.ypoints,5);
				}
                        }
                }
        else
                {
                end = b.shades;
                step = width*0.3/end;
                for(int i=0;i<end;i++)
                        {
                        p = new Polygon();
                        p.addPoint((int)Math.round(xoffset + i*step), (int)Math.round(yoffset));
                        p.addPoint((int)Math.round(xoffset + i*step), (int)Math.round(yoffset+legendheight));
                        p.addPoint((int)Math.round(xoffset + (i+1)*step), (int)Math.round(yoffset+legendheight));
                        p.addPoint((int)Math.round(xoffset + (i+1)*step), (int)Math.round(yoffset));
                        p.addPoint((int)Math.round(xoffset + i*step), (int)Math.round(yoffset));

                        if(b.shades == 1)
                            r = 127;
                        else
                            r = i*255/(b.shades-1);

                        if(b.swap)
                                g.setColor(new Color(255-r,255-r,255-r));
                        else
                                g.setColor(new Color(r,r,r));

                        g.fillPolygon(p.xpoints,p.ypoints,5);
                        }
                }

        p = new Polygon();
        p.addPoint((int)Math.round(xoffset), (int)Math.round(yoffset));
        p.addPoint((int)Math.round(xoffset), (int)Math.round(yoffset+legendheight));
        p.addPoint((int)Math.round(xoffset+0.3*width), (int)Math.round(yoffset+legendheight));
        p.addPoint((int)Math.round(xoffset+0.3*width), (int)Math.round(yoffset));
        p.addPoint((int)Math.round(xoffset), (int)Math.round(yoffset));
        g.setColor(Color.black);
        g.drawPolygon(p.xpoints,p.ypoints,5);

        String str;
	int center = 0;
        double max = b.upperlimit;
        double min = b.lowerlimit;

        g.setFont(new Font("TimesRoman", Font.PLAIN, b.fontsize));

        str = String.valueOf(min);
	if(centering)
		center = (int)(str.length()*b.fontsize*0.25);
        g.drawString(str,(int)Math.round(xoffset)-center,(int)Math.round(yoffset+legendheight+1.3*b.fontsize));
        g.drawLine((int)Math.round(xoffset),(int)Math.round(yoffset+legendheight),(int)Math.round(xoffset),(int)Math.round(yoffset+legendheight+3));
        if(min < 0 && 0 < max)
            if(b.lowerlimit + (b.upperlimit-b.lowerlimit)/3 < 0 && 0 < b.upperlimit - (b.upperlimit-b.lowerlimit)/3)
                {
                str = String.valueOf(0);
		if(centering)
			center = (int)(str.length()*b.fontsize*0.25);
                g.drawString(str,(int)Math.round(xoffset-min/(max-min)*0.3*width)-center,(int)Math.round(yoffset+legendheight+1.3*b.fontsize));
                g.drawLine((int)Math.round(xoffset-min/(max-min)*0.3*width),(int)Math.round(yoffset+legendheight),
                           (int)Math.round(xoffset-min/(max-min)*0.3*width),(int)Math.round(yoffset+legendheight+3));
                }
        str = String.valueOf(max);
	if(centering)
		center = (int)(str.length()*b.fontsize*0.25);
        g.drawString(str,(int)Math.round(xoffset+0.3*width)-center,(int)Math.round(yoffset+legendheight+1.3*b.fontsize));
        g.drawLine((int)Math.round(xoffset+0.3*width),(int)Math.round(yoffset+legendheight),(int)Math.round(xoffset+0.3*width),(int)Math.round(yoffset+legendheight+3));

        }

private void drawNA(Graphics g, Polygon p)
        {

        int x,y;
        int dist = 5;
        Rectangle rect = p.getBounds();

        int start;

        start = rect.y;
        while(start < rect.y+rect.height)
            {
            x = rect.x;
            y = start;
            while(x <= rect.x+rect.width && y <= rect.y+rect.height)
                {
                if(p.contains(x,y,1,1))
                        g.drawLine(x,y,x+1,y+1);
                y++;
                x++;
                }
            start = start + dist;
            }

        start = rect.x;
        while(start < rect.x+rect.width)
            {
            y = rect.y;
            x = start;
            while(x <= rect.x+rect.width && y <= rect.y+rect.height)
                {
                if(p.contains(x,y,1,1))
                        g.drawLine(x,y,x+1,y+1);
                y++;
                x++;
                }
            start = start + dist;
            }

        }

private void drawtitle(Graphics g,double height,double width,boolean centering)
        {
	int center = 0;
	if(centering)
  		center = (int)(b.title.length()*b.titlescale*b.fontsize*0.5);

        g.setFont(new Font("TimesRoman", Font.BOLD, (int)(b.titlescale*b.fontsize)));
        g.drawString(b.title,(int)(width/2-center),(int)(0.05*height+b.titlescale*b.fontsize));
        g.setFont(new Font("TimesRoman", Font.PLAIN, b.fontsize));
        }


public void plotnonp(Graphics g)
        {
        setplotparam((PAGEWIDTH-b.width)/2,120,b.width,b.height,b.pointsize,b.fontsize,b.titlescale);
        plotframe(g,0,true);
        for(int i=1;i<b.getDCols();i++)
            plot(g,i);
        }

public void Saveplotnonp(PrintWriter out)
        {

        setplotparam((PAGEWIDTH-b.width)/2,120,b.width,b.height,b.pointsize,b.fontsize,b.titlescale);

	int fontwidth = (int)(0.6*fontsize);
	int fontheight = (int)(1.5*fontsize);

        out.println("%!PS-Adobe-3.0");
        out.print("%%BoundingBox:");

        if(b.ylab.equals(""))
	    out.print(xoffset-6*fontwidth-6);
	else
	    out.print(xoffset-6*fontwidth-fontheight);
	out.print(" ");

        if(b.xlab.equals(""))
            out.print(PAGEHEIGHT-yoffset-height-fontheight-2);
        else
            out.print(PAGEHEIGHT-yoffset-height-2*fontheight-2);
        out.print(" ");

        out.print(xoffset+width+2);
        out.print(" ");

	if(b.title.equals(""))
            out.println(PAGEHEIGHT-yoffset+2);
        else
            out.println(PAGEHEIGHT-yoffset+1.5*titlescale*fontheight);

        PSGr2 g = new PSGr2(out);
	PSGr = g;

        scale = 10;
        pointsize = pspointsize;

        out.println("%%Pages:(atend)");
        out.println("%%Page:1 1");
        out.println("0.1 0.1 scale");
	out.print(b.linewidth);
        out.println(" setlinewidth");

        plotframe(g,0,false);

        for(int i=1;i<b.getDCols();i++)
            plot(g,i);

        out.println("showpage");
        out.close();

        }

private void plotframe(Graphics g, int col, boolean centering)
        {

        boolean date = false;

        String str;

        int center = 0;

        int nrticks_x = 5;
        int nrticks_y = 5;

	int fontwidth = (int)(0.6*fontsize);
	int fontheight = (int)(1.5*fontsize);

// Wertebereich
        if(b.function == 3)
                {
                minX = b.getMin(0);
                maxX = b.getMax(0);
                minY = b.getMin(1);
                maxY = b.getMax(1);
                for(int i=2;i<b.getDCols();i++)
                    if(b.getMin(i)<minY)
                        minY = b.getMin(i);
                for(int i=2;i<b.getDCols();i++)
                    if(b.getMax(i)>maxY)
                        maxY = b.getMax(i);

                if(b.xmax > -MAXDOUBLE)
                    maxX = b.xmax;
                if(b.xmin < MAXDOUBLE)
                    minX = b.xmin;

                if(b.ymax > -MAXDOUBLE)
                    maxY = b.ymax;
                if(b.ymin < MAXDOUBLE)
                    minY = b.ymin;

                if(b.year != 0 && b.month != 0)
                    date = true;

                }
        else if(b.function == 4)
                {
                minX = b.getMin(0);
                maxX = b.getMax(0);
                minY = b.getMin(col);
                maxY = b.getMax(col);
                }
        else if(b.function == 5)
                {
                minX = b.getMin(0);
                maxX = b.getMax(0);
                minY = -0.1;
                maxY = 1.0;
                }

        if(maxY==minY)
                {
                maxY = maxY*1.1;
                minY = minY*0.9;
                }

	xstart = minX;
	ystart = minY;

	if(b.function == 3)
 		{
		if(b.xstart != MAXDOUBLE)
			xstart = b.xstart;
		if(b.ystart != MAXDOUBLE)
			ystart = b.ystart;
		}

        xstep = (maxX-xstart)/(nrticks_x-1);
        ystep = (maxY-ystart)/(nrticks_y-1);

        if(b.function == 3)
            {
            if(b.xstep > 0.0)
                {
                xstep = b.xstep;
                nrticks_x = (int)((maxX-xstart)/xstep)+1;
                }
            if(b.ystep > 0.0)
                {
                ystep = b.ystep;
                nrticks_y = (int)((maxY-ystart)/ystep)+1;
                }
            }
        else
            {
            if((int)(maxX%4)==0)
                {
                nrticks_x = 5;
                xstep = maxX/4;
                }
            else if((int)(maxX%3)==0)
                {
                nrticks_x = 4;
                xstep = maxX/3;
                }
            else if((int)(maxX%5)==0)
                {
                nrticks_x = 6;
                xstep = maxX/5;
                }
            else if((int)(maxX%2)==0)
                {
                nrticks_x = 3;
                xstep = maxX/2;
                }
            else
                {
                nrticks_x = 5;
                xstep = (maxX - maxX%4)/4;
                }
            }

        g.setColor(Color.black);
        g.drawRect(scale*xoffset, PAGEHEIGHT-scale*(PAGEHEIGHT-yoffset), scale*width, scale*height);
	if(centering)
	        g.setFont(new Font("Monospaced", Font.BOLD, scale*fontsize));
	else
	        g.setFont(new Font("TimesRoman", Font.BOLD, scale*fontsize));

// Beschriftung und Ticks x-Achse
        if(b.function == 3)
                {
		if(date)
                    {

                    int start = (int)(b.getDoubleValue(0,0));
                    int stop = (int)(b.getDoubleValue(b.getDRows()-1,0));
                    int year = b.year;
                    int step = 12;
                    if(b.xstep != 0.0)
                        step = (int)(b.xstep);
                    int i = start;
                    int j = b.month;

                    if(stop-start < 24)
                        {
                        String[] names = {};
                        if(xstep == 12)
                            names = new String[]{"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
                        else if(xstep == 4)
                            names = new String[]{"Jan","Apr","Jul","Oct"};
                        else if(xstep == 2)
                            names = new String[]{"Jan","Jul"};

                        while(i <= stop)
                            {
                            if(stop-start < 12 || (j-1)%2 == 0)
                                {
                                str = String.valueOf(names[j-1]);
                                str = formatLabel(str,centering);
		                if(centering)
	        		        center = (str.length()+1)*fontwidth/2;
                                g.drawString(str,translateX(i)-scale*center,PAGEHEIGHT-scale*(PAGEHEIGHT-(yoffset+height+17*fontheight/20)));
                                }
                            if((j-1)%step == 0)
                                {
                                str = String.valueOf(year);
                                str = formatLabel(str,centering);
		                if(centering)
	        		        center = (str.length()+1)*fontwidth/2;
                                g.drawString(str,translateX(i)-scale*center,PAGEHEIGHT-scale*(PAGEHEIGHT-(yoffset+height+29*fontheight/20)));
                                }
                            if(j%step == 0)
                                year++;
                            if(j == step)
                                j = 0;

                            g.drawLine(translateX(i), PAGEHEIGHT-scale*(PAGEHEIGHT-yoffset-height),
				       translateX(i), PAGEHEIGHT-scale*(PAGEHEIGHT-yoffset-height-5));

                            i++;
                            j++;
                            }
                        }   // end: <24
                    else
                        {
                        int k = (int)((stop-start)/6);
                        while(i <= stop)
                            {
                            if((j-1)%k == 0)
                                {
                                str = String.valueOf(year);
                                str = formatLabel(str,centering);
		                if(centering)
	        		        center = (str.length()+1)*fontwidth/2;
                                g.drawString(str,translateX(i)-scale*center,PAGEHEIGHT-scale*(PAGEHEIGHT-(yoffset+height+fontheight)));
                                }
                            if((j-1)%step == 0)
                                g.drawLine(translateX(i), PAGEHEIGHT-scale*(PAGEHEIGHT-yoffset-height),
					   translateX(i), PAGEHEIGHT-scale*(PAGEHEIGHT-yoffset-height-5));
                            if(i%step == 0)
                                year++;

                            i++;
                            j++;
                            }
                        }
                    }   // end: date
                else
                    {
                    for(int i=0;i<nrticks_x;i++)
                            {
                            double value = (double)(Math.round((xstart+i*xstep)*1.0E10))/1.0E10;
                            str = String.valueOf(value);
                            str = formatLabel(str,centering);
			    double pos = xstart+i*xstep;
			    if(xstart <= maxX && pos >= minX && pos <= maxX*1.0000000001)
 				    {
			            if(centering)
		             	    	center = (str.length()+1)*fontwidth/2;
	                            g.drawString(str,translateX(pos)-scale*center, PAGEHEIGHT-scale*(PAGEHEIGHT-(yoffset+height+fontheight)));
        	                    g.drawLine(translateX(pos), PAGEHEIGHT-scale*(PAGEHEIGHT-yoffset-height),
				               translateX(pos), PAGEHEIGHT-scale*(PAGEHEIGHT-yoffset-height-5));
				    }
                            }
                    }
                }       // end: function 3
        else
                {
                str = String.valueOf(minX);
                str = str.substring(0,str.indexOf('.'));
		if(centering)
	                center = str.length()*fontwidth/2;
                g.drawString(str,translateX(minX)-scale*center,PAGEHEIGHT-scale*(PAGEHEIGHT-(yoffset+height+fontheight)));
                for(int i=1;i<nrticks_x;i++)
                        {
                        str = String.valueOf(i*xstep);
                        str = str.substring(0,str.indexOf('.'));
			if(centering)
		                center = str.length()*fontwidth/2;
                        g.drawString(str,translateX(i*xstep)-scale*center, PAGEHEIGHT-scale*(PAGEHEIGHT-(yoffset+height+fontheight)));
                        g.drawLine(translateX(i*xstep), PAGEHEIGHT-scale*(PAGEHEIGHT-yoffset-height),
				   translateX(i*xstep), PAGEHEIGHT-scale*(PAGEHEIGHT-yoffset-height-5));
                        }
                g.drawLine(translateX(minX), PAGEHEIGHT-scale*(PAGEHEIGHT-yoffset-height),
			   translateX(minX), PAGEHEIGHT-scale*(PAGEHEIGHT-yoffset-height-5));
                }

// Ticks auf y-Achse
        if(b.function == 5)
                {
                for(int i=0;i<6;i++)
                        g.drawLine(scale*xoffset, translateY(i*0.2), scale*(xoffset-5), translateY(i*0.2));
                }
        else
                {
                for(int i=0;i<nrticks_y;i++)
			{
			double pos = ystart+i*ystep;
			if(ystart <= maxY && pos >= minY && pos <= maxY*1.0000000001)
	                        g.drawLine(scale*xoffset, translateY(pos), scale*(xoffset-5), translateY(pos));
			}
                }

// Beschriftung y-Achse
        if(b.function == 5)
                {
                for(int i=0;i<6;i++)
                        {
                        str = String.valueOf(i*0.2);
                        str = str.substring(0,3);
			if(centering)
	                        center = (str.length()+1)*fontwidth+5 - 3*fontwidth;
                        g.drawString(str,scale*(xoffset-center-3*fontwidth),translateY(i*0.2)+scale*fontsize/3);
                        }
                }
        else
                {
                for(int i=0;i<nrticks_y;i++)
                        {
                        double value = (double)(Math.round((ystart+i*ystep)*1.0E10))/1.0E10;
                        str = String.valueOf(value);
                        str = formatLabel(str,centering);
			if(centering)
	                        center = (str.length()+1)*fontwidth+5 - 3*fontwidth;
			double pos = ystart+i*ystep;
			if(ystart <= maxY && pos >= minY && pos <= maxY*1.0000000001)
	                        g.drawString(str,scale*(xoffset-center-3*fontwidth),translateY(pos)+scale*fontsize/3);
                        }
                }

// title
        if(b.function == 3)
                {
                if(!b.title.equals(""))
                        {
			if(centering)
	                        g.setFont(new Font("Monospaced", Font.BOLD, (int)(scale*fontsize*titlescale)));
			else
	                        g.setFont(new Font("TimesRoman", Font.BOLD, (int)(scale*fontsize*titlescale)));
                        str = b.title;
                        if(str.length()>32)
                          str = str.substring(0,32);
			if(centering)
	                        center = (int)((str.length()*titlescale*fontwidth)/2);
                        g.drawString(str,scale*((int)(xoffset+width/2)-center),PAGEHEIGHT-scale*(PAGEHEIGHT-(yoffset-(int)(titlescale*fontheight))));
			if(centering)
	                        g.setFont(new Font("Monospaced", Font.BOLD, scale*fontsize));
			else
	                        g.setFont(new Font("TimesRoman", Font.BOLD, scale*fontsize));
                        }
                if(!b.xlab.equals(""))
                        {
                        str = b.xlab;
			if(centering)
	                        center = (int)((str.length()*fontwidth)/2);
                        g.drawString(str,scale*((int)(xoffset+width/2)-center),PAGEHEIGHT-scale*(PAGEHEIGHT-(yoffset+height+2*fontheight)));
                        }
                if(!b.ylab.equals(""))
                        {
                        str = b.ylab;
			if(centering)
				{
	                        if(str.length()>14)
					center = 7*fontwidth;
				else
		                        center = (int)((str.length()*fontwidth)/2);
				}
	                if(centering)
				g.drawString(str,scale*(xoffset-center),PAGEHEIGHT-scale*(PAGEHEIGHT-(yoffset-3*fontheight/4)));
			else
				{
				PSGr.emitThis("gsave");
				PSGr.translate(scale*(xoffset-48),scale*(PAGEHEIGHT-(yoffset+height/2)));
				PSGr.emitThis("90 rotate");
				PSGr.drawString(str,0,PAGEHEIGHT);
				PSGr.emitThis("grestore");
				}
                        }
                }
        else if(b.function == 4)
                {
                str = "iteration";
		if(centering)
	                center = (int)((str.length()*fontwidth)/2);
                g.drawString(str,scale*((int)(xoffset+width/2)-center),PAGEHEIGHT-scale*(PAGEHEIGHT-(yoffset+height+2*fontheight)));

                str = "par "+col;
		if(centering)
			{
	                if(str.length()>14)
				center = 5*fontwidth;
			else
		                center = (int)((str.length()*fontwidth)/2);
			}
                g.drawString(str,scale*(xoffset-center),PAGEHEIGHT-scale*(PAGEHEIGHT-(yoffset-fontheight)));
                }
        else if(b.function == 5)
                {
                str = "lag";
		if(centering)
	                center = (int)((str.length()*fontwidth)/2);
                g.drawString(str,scale*((int)(xoffset+width/2)-center),PAGEHEIGHT-scale*(PAGEHEIGHT-(yoffset+height+2*fontheight)));

                str = b.getVarname(col);
		if(centering)
			{
	                if(str.length()>10)
				center = 5*fontwidth;
			else
		                center = (int)((str.length()*fontwidth)/2);
			}
                g.drawString(str,scale*(xoffset-center),PAGEHEIGHT-scale*(PAGEHEIGHT-(yoffset-fontheight)));
                }

        }

private int translateX(double x)
        {
        return (int)Math.round( scale*(0.05*width + (x-minX)*0.9*width/(maxX-minX)) ) + scale*xoffset;
        }

private int translateY(double y)
        {
        return (int)Math.round( PAGEHEIGHT - scale*(PAGEHEIGHT-yoffset-height + 0.05*height + (y-minY)*0.9*height/(maxY-minY)) );
        }

private void plotsample(Graphics g)
        {

        width = 165;
        height = 130;

        int start = page*b.plotsperpage+1;
        int count = 0;

        yoffset = 50-(height+100);
        while(start < b.getDCols() && count < b.plotsperpage)
                {
                if(start%2!=0)
                        {
                        xoffset = 100;
                        yoffset += height+100;
                        }
                else
                        {
                        xoffset = 100+width+66;
                        }
                plotframe(g,start,true);
                plot(g,start);
                start++;
                count++;
                }

        }

public void Saveplotsample(PrintWriter out)
	{

        width = 165;
        height = 130;
        xoffset = 100;

        out.println("%!PS-Adobe-3.0");
        out.print("%%BoundingBox:");
        out.print(xoffset-6*11*12/18-6);
        out.print(" ");
        out.print(PAGEHEIGHT-100-2*100-3*130-42);
        out.print(" ");
        out.print(xoffset+2*width+66+5);
        out.print(" ");
        out.println(PAGEHEIGHT-100+30);
        out.println("%%Pages:(atend)");

        PSGr2 g = new PSGr2(out);

        scale = 10;
        pointsize = pspointsize;

        int start = 1;
        int nrplots = b.getDCols();
        int nrpages = (b.getDCols()-2)/b.plotsperpage+1;
        int count;

        for(int i=1;i<=nrpages;i++)
                {
                out.println("%%Page:"+i+" "+i);
                out.println("0.1 0.1 scale");
                out.println("3 setlinewidth");

                yoffset = 100-(height+100);
                count = 0;
                while(start < nrplots && count < 6)
                        {
                        if(start%2!=0)
                                {
                                xoffset = 100;
                                yoffset += height+100;
                                }
                        else
                                {
                                xoffset = 100+width+66;
                                }

                        plotframe(g,start,false);
                        plot(g,start);

                        start = start+1;
                        count = count+1;
                        }
                out.println("showpage");
                }

        out.close();

        }

public void Saveplotautocor(PrintWriter out)
        {

        width = 165;
        height = 130;
        xoffset = 100;

        out.println("%!PS-Adobe-3.0");
        out.print("%%BoundingBox:");
        out.print(xoffset-4*11*12/18-6);
        out.print(" ");
        out.print(PAGEHEIGHT-100-2*100-3*130-42);
        out.print(" ");
        out.print(100+2*width+66+5);
        out.print(" ");
        out.println(PAGEHEIGHT-100+30);
        out.println("%%Pages:(atend)");

        PSGr2 g = new PSGr2(out);

        scale = 10;
        pointsize = pspointsize;

        int start = 1;
        int nrplots = b.getDCols();
        int nrpages = (b.getDCols()-2)/b.plotsperpage+1;
        int count;

        for(int i=1;i<=nrpages;i++)
                {
                out.println("%%Page:"+i+" "+i);
                out.println("0.1 0.1 scale");
                out.println("3 setlinewidth");

                yoffset = 100-(height+100);
                count = 0;
                while(start < nrplots && count < b.plotsperpage)
                        {

                        if(b.plotsperpage == 6)
                        {
                        if(count%2==0)
                                {
                                xoffset = 100;
                                yoffset += height+100;
                                }
                        else
                                {
                                xoffset = 100+width+66;
                                }
                        }
                        else if(b.plotsperpage == 3)
                        {
                        width = 396;
                        xoffset = 100;
                        yoffset += height+100;
                        }

                        plotframe(g,start,false);
                        plot(g,start);
                        drawLine(g);

                        start = start+1;
                        count = count+1;
                        }
                out.println("showpage");
                }

        out.close();

        }


private void plotautocor(Graphics g)
        {
        width = 165;
        height = 130;

        int start = page*b.plotsperpage+1;
        int count = 0;

        yoffset = 50-(height+100);
        while(start < b.getDCols() && count < b.plotsperpage)
                {

                if(b.plotsperpage == 6)
                {
                if(count%2==0)
                        {
                        xoffset = 100;
                        yoffset += height+100;
                        }
                else
                        {
                        xoffset = 100+width+66;
                        }
                }
                else if(b.plotsperpage == 3)
                {
                width = 396;
                xoffset = 100;
                yoffset += height+100;
                }

                plotframe(g,start,true);
                plot(g,start);
                drawLine(g);
                start++;
                count++;
                }
        }



// Function to plot surface from the given 3D-points data
private void plotsurf(Graphics g)
	{
	setplotparam((PAGEWIDTH-500)/2,120,500,500,b.pointsize,b.fontsize,b.titlescale);
	double x[][] = new double [b.getDRows()][3];
	for(int i=0;i<b.getDRows();i++)
		{
		x[i][0] = b.getDoubleValue(i,0);
		x[i][1] = b.getDoubleValue(i,1);
		x[i][2] = b.getDoubleValue(i,2);
		}

	char color = 'G';
	if(b.linecolor.length()>0)
          color = b.linecolor.charAt(0);

    Plot3D plot1 = new Plot3D(x,b.gridsize,color,b.title,b.xlab,b.ylab,b.zlab,b.xstart,b.xstep,b.ystart,b.ystep,b.zstart,b.zstep,b.xrot,b.yrot,b.zrot);
	plot1.Plot(g);
	return;
	}



public void Saveplotsurf(PrintWriter out)
	{
	 out.println("%!PS-Adobe-3.0");
     out.println("%%Pages:1");
     out.println("%%Page:1 1");

	 PSGr2 g = new PSGr2(out);
	 plotsurf(g);

	 out.println("showpage");
     out.close();
    }

private void plot(Graphics g,int col)
        {

	if(b.linecolor!=null && b.linecolor.length()>=col)
		{
		if(b.linecolor.charAt(col-1)=='B')
			{
			g.setColor(Color.black);
			}
		else if(b.linecolor.charAt(col-1)=='G')
			{
			g.setColor(Color.gray);
			}
		else if(b.linecolor.charAt(col-1)=='r')
			{
			g.setColor(Color.red);
			}
		else if(b.linecolor.charAt(col-1)=='g')
			{
			g.setColor(Color.green);
			}
		else if(b.linecolor.charAt(col-1)=='b')
			{
			g.setColor(Color.blue);
			}
		else if(b.linecolor.charAt(col-1)=='c')
			{
			g.setColor(Color.cyan);
			}
		else if(b.linecolor.charAt(col-1)=='m')
			{
			g.setColor(Color.magenta);
			}
		else if(b.linecolor.charAt(col-1)=='o')
			{
			g.setColor(Color.orange);
			}
		else if(b.linecolor.charAt(col-1)=='y')
			{
			g.setColor(Color.yellow);
			}
		}

        if(b.function==5)
                {
                if(b.connect.equals("points") || (b.connect.length()>=col&&b.connect.charAt(col-1)=='p'))
                        {
                        for(int i=0;i<b.getDRows();i++)
                                if(b.getDoubleValue(i,col)>=-0.16)
                                        g.fillOval(translateX(b.getDoubleValue(i,0)),translateY(b.getDoubleValue(i,col)),pointsize,pointsize);
                        }
                else
                        {
                        int x[] = new int[b.getDRows()];
                        int y[] = new int[b.getDRows()];

                        for(int i=0;i<b.getDRows();i++)
                                {
                                x[i] = translateX(b.getDoubleValue(i,0));
                                if(b.getDoubleValue(i,col)<-0.16)
                                    y[i] = translateY(-0.16);
                                else
                                    y[i] = translateY(b.getDoubleValue(i,col));
                                }
                        g.drawPolyline(x,y,b.getDRows());
                        }
                }
        else
                {
                if( (b.connect.length()>=col&&b.connect.charAt(col-1)=='p') ||
		    (b.connect.length()>=col&&b.connect.charAt(col-1)=='5') )
                        {
                        double x;
                        double y;
                        for(int i=0;i<b.getDRows();i++)
				{
                                x = b.getDoubleValue(i,0);
                                y = b.getDoubleValue(i,col);
                                if(minX<=x&x<=maxX&minY<=y&y<=maxY)
					g.fillOval(translateX(x),translateY(y),pointsize,pointsize);
				}
                        }
		else if( (b.connect.length()>=col&&b.connect.charAt(col-1)=='-') ||
                         (b.connect.length()>=col&&b.connect.charAt(col-1)=='4') )
			{
			dashedLine(g,30,0.5,col);
			}
		else if( (b.connect.length()>=col&&b.connect.charAt(col-1)=='_') ||
                         (b.connect.length()>=col&&b.connect.charAt(col-1)=='3') )
			{
			dashedLine(g,25,0.7,col);
			}
		else if( (b.connect.length()>=col&&b.connect.charAt(col-1)=='d') ||
                         (b.connect.length()>=col&&b.connect.charAt(col-1)=='2') )
			{
			dashedLine(g,12,0.7,col);
			}
		else
                        {
                        double x;
                        double y;
                        double x2;
                        double y2;
                        for(int i=1;i<b.getDRows();i++)
                                {
                                x = b.getDoubleValue(i-1,0);
                                y = b.getDoubleValue(i-1,col);
                                x2 = b.getDoubleValue(i,0);
                                y2 = b.getDoubleValue(i,col);
                                if(minX<=x&x<=maxX&minX<=x2&x2<=maxX&minY<=y&y<=maxY&minY<=y2&y2<=maxY)
                                    g.drawLine(translateX(x),translateY(y),translateX(x2),translateY(y2));
                                }
                        }
                }

	g.setColor(Color.black);

        }

private void dashedLine(Graphics g, int intervals, double frac, int col)
        {
	double xi1;
	double xi2;
	double yi1;
	double yi2;
        double x1;
        double y1;
        double x2;
        double y2;
	double xstart;
	double xend1;
	double xend2;
	double xi1help;
	double yi1help;
	double xi2help;
	double yi2help;
	double help;

	int nrobs = b.getDRows();	
	int start;
	int stop;

	int i;
	double dist = (maxX-minX)/intervals;

        for(int j=0;j<intervals;j++)
	        {

		xstart = minX + j*dist;
        	xend1 = minX + (j+frac)*dist;
		xend2 = minX + (j+1.0)*dist;

		start = -1;
		stop = -2;

		xi1 = b.getDoubleValue(0,0);
		xi2 = b.getDoubleValue(nrobs-1,0);
		yi1 = b.getDoubleValue(0,col);
		yi2 = b.getDoubleValue(nrobs-1,col);

		if(xstart<xi1)
			xstart=xi1;
		if(xend1>xi2)
			xend1=xi2;

		xi1help = minX-1;	
		yi1help = minX-1;
		xi2help = maxX+1;	
		yi2help = maxX+1;

		i=0;
		help = b.getDoubleValue(i,0);

		while(i+1<nrobs && help<xstart)
			{
			i++;
			help = b.getDoubleValue(i,0);
			}
		if(xstart<=help && help<=xend1)
			{
   			xi1help = help;
			yi1help = b.getDoubleValue(i,col);
			start = i;
			}

		while(i+1<nrobs && help<xend1)
			{
			i++;
			help = b.getDoubleValue(i,0);
			}
		xi2 = help;
		yi2 = b.getDoubleValue(i,col);


		i=nrobs-1;
		help = b.getDoubleValue(i,0);

		while(i>0 && help>xend1)
			{
			i--;
			help = b.getDoubleValue(i,0);
			}
		if(xstart<=help && help<=xend1)
			{
   			xi2help = help;
			yi2help = b.getDoubleValue(i,col);
			stop = i;
			}

		while(i>0 && help>xstart)
			{
			i--;
			help = b.getDoubleValue(i,0);
			}
		xi1 = help;
		yi1 = b.getDoubleValue(i,col);


		if(xi1help<minX && xi2help>maxX) // keine Beobachtung in [xstart,xend1]
			{
			x1 = xstart;			
			y1 = yi1 + (xstart-xi1) * (yi2-yi1)/(xi2-xi1);
			x2 = xend1;
			y2 = yi1 + (xend1-xi1) * (yi2-yi1)/(xi2-xi1);
		
			if(x1!=x2 && minX<=x1&x1<=maxX&minX<=x2&x2<=maxX&minY<=y1&y1<=maxY&minY<=y2&y2<=maxY)
       		                g.drawLine(translateX(x1),translateY(y1),translateX(x2),translateY(y2));
			else if(x1==x2)
       		                g.fillOval(translateX(x1),translateY(y1),3,3);

			}
		else // mindestens eine Beobachtung in [xstart,xend1]
			{
			x1 = xstart;			
			y1 = yi1 + (xstart-xi1) * (yi1help-yi1)/(xi1help-xi1);
			x2 = xi1help; 
			y2 = yi1help;
		
			if(x1!=x2 && minX<=x1&x1<=maxX&minX<=x2&x2<=maxX&minY<=y1&y1<=maxY&minY<=y2&y2<=maxY)
       		                g.drawLine(translateX(x1),translateY(y1),translateX(x2),translateY(y2));
			else if(x1==x2)
       		                g.fillOval(translateX(x1),translateY(y1),3,3);

			x1 = xi2help;
			y1 = yi2help;
			x2 = xend1;
			y2 = yi1 + (xend1-xi1) * (yi2-yi2help)/(xi2-xi2help);
		
			if(x1!=x2 && minX<=x1&x1<=maxX&minX<=x2&x2<=maxX&minY<=y1&y1<=maxY&minY<=y2&y2<=maxY)
       		                g.drawLine(translateX(x1),translateY(y1),translateX(x2),translateY(y2));
			else if(x1==x2)
       		                g.fillOval(translateX(x1),translateY(y1),3,3);

			}	


		for(int k=start;k<stop;k++)
			{
			x1 = b.getDoubleValue(k,0); 
			y1 = b.getDoubleValue(k,col);
			x2 = b.getDoubleValue(k+1,0);
			y2 = b.getDoubleValue(k+1,col);

			if(x1!=x2 && minX<=x1&x1<=maxX&minX<=x2&x2<=maxX&minY<=y1&y1<=maxY&minY<=y2&y2<=maxY)
       		                g.drawLine(translateX(x1),translateY(y1),translateX(x2),translateY(y2));
			else if(x1==x2)
       		                g.fillOval(translateX(x1),translateY(y1),3,3);

			}

		}

        }

private void drawLine(Graphics g)
        {
        double start;
        double end;
        for(int i=0;i<100;i++)
                {
                if(i%3==0)
                    {
                    start = minX + i*(maxX-minX)/101;
                    end = start + (maxX-minX)/101;
                    g.drawLine(translateX(start),translateY(0.1),translateX(end),translateY(0.1));
                    }
                }
        }

private String formatLabel(String str, boolean centering)
        {

        if(str.indexOf('E')>-1)
            {
            if(str.indexOf('.')>-1)
                str = str.substring(0,str.indexOf('.')+2) + str.substring(str.indexOf('E'));
            }
        else if(str.indexOf('e')>-1)
            {
            if(str.indexOf('.')>-1)
                str = str.substring(0,str.indexOf('.')+2) + str.substring(str.indexOf('e'));
            }
        else if(str.indexOf('.')>-1)
                {

                int prec;

                if(str.charAt(0)!='-')          // um Position von '.' richtig festzustellen
                    str = ' '+str;

                if(str.indexOf('.') < 3)        // weniger als 2 Ziffern vor dem '.'
                    {
                    prec = 2;
                    while(Math.pow(10.0,(prec-2)) * ystep < 0.1)
                        prec++;
                    double d = round(Double.parseDouble(str),prec);
                    str = String.valueOf(d);
                    }
                else if(str.indexOf('.') == 3)  // genau 2 Ziffern vor dem '.'
                    {
                    prec = 1;
                    double d = round(Double.parseDouble(str),prec);
                    str = String.valueOf(d);
                    }
                else                            // mehr als zwei Ziffern vor dem '.'
                    {
                    double d = round(Double.parseDouble(str),0);
                    str = String.valueOf(d);
                    }

                // '.0' am Ende entfernen
                if(str.charAt(str.length()-2)=='.'&str.charAt(str.length()-1)=='0')
                    str = str.substring(0,str.indexOf('.'));

                }
        else
                {
                // keine Formatierung fuer ganze Zahlen
                }

	if(centering)
	        if(str.charAt(0)!='-')
        	    str = ' '+str;          // erstes Zeichen '-' oder ' '

        return str;
        }


private double round(double d,int n)
        {
        double x = Math.pow(10.0,n);
        d = Math.round(d*x)/x;
        return d;
        }


static double gtrans(double u)
	{
	if (u > 0.00304)
	        return 1.055 * Math.pow(u, (1 / GAMMA)) - 0.055;
	else
        	return 12.92 * u;
	}

static int[] hcl2rgb(double h, double c, double l)
	{
	double L, U, V;
	double u, v;
	double X, Y, Z;
	double R, G, B;

    /* Step 1 : Convert to CIE-LUV */

	h = DEG2RAD * h;
	L = l;
    	U = c * Math.cos(h);
	V = c * Math.sin(h);

    /* Step 2 : Convert to CIE-XYZ */

	if (L <= 0 && U == 0 && V == 0) 
		{
	        X = 0; Y = 0; Z = 0;
	    	}
	else 
		{
	        Y = WHITE_Y * ((L > 7.999592) ? Math.pow((L + 16)/116, 3) : L / 903.3);
        	u = U / (13 * L) + WHITE_u;
	        v = V / (13 * L) + WHITE_v;
        	X =  9.0 * Y * u / (4 * v);
	        Z =  - X / 3 - 5 * Y + 3 * Y / v;
		}

    /* Step 4 : CIE-XYZ to sRGB */

	R = 255.0 * gtrans(( 3.240479 * X - 1.537150 * Y - 0.498535 * Z) / WHITE_Y);
	G = 255.0 * gtrans((-0.969256 * X + 1.875992 * Y + 0.041556 * Z) / WHITE_Y);
	B = 255.0 * gtrans(( 0.055648 * X - 0.204043 * Y + 1.057311 * Z) / WHITE_Y);

	if(R>255.0)
 		R = 255.0;
	if(R<0.0)
		R = 0.0;

	if(G>255.0)
 		G = 255.0;
	if(G<0.0)
		G = 0.0;

	if(B>255.0)
 		B = 255.0;
	if(B<0.0)
		B = 0.0;

	int[] RGB = new int[3];
        RGB[0] = (int)R;
        RGB[1] = (int)G;
        RGB[2] = (int)B;
	return RGB;
	}



}        // END: class MapPanel
