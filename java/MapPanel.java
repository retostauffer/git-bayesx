import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.util.*;
import java.io.*;

import gov.sandia.postscript.PSGr1;
import gov.sandia.postscript.PSGr2;

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

protected short page;

private double xstep;
private double ystep;
private double xstart;
private double ystart;

private int xoffset = 120;
private int yoffset = 120;
private int width = 356;
private int height = 210;

private int scale = 1;
private int pointsize = 3;
private int pspointsize = 20;
private int fontsize = 12;

private double minX;
private double maxX;
private double minY;
private double maxY;

public MapPanel(BayesX b)
	{
	super();
	this.b=b;
	}

public void setplotparam(int xoff, int yoff, int w, int h, int ps, int fs)
        {
        xoffset = xoff;
        yoffset = yoff;
        width = w;
        height = h;
	pspointsize = ps;
	fontsize = fs;
        }
        
public void paintComponent(Graphics g1)                 // f�r Bildschirmanzeige        
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
        
//        resize((int)width,(int)height);        
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
            offset = offset + 10 + 3*b.fontsize/2;
        
        helpdrawmap(g,height,width,offset);        
        
        }

public void SaveMap(PrintWriter out)              // Zum Speichern als PostScript   
	{
        
	double[] d2 = new double[4];

        b.getboundaries(d2);
	maxX = d2[2];
	maxY = d2[3];
	minX = d2[0];
	minY = d2[1];

	double width = 596;
	double height = width*(maxY-minY)/(maxX-minX);

        if(height>842)
                {
                width = width*842/height;
                height = height*842/height;             
                }        

        out.println("%!PS-Adobe-3.0");      
        out.print("%%BoundingBox:");
        out.print((int)Math.round(0.05*width));
        out.print(" ");
        out.print(842-(int)Math.round(0.95*height+2));
        out.print(" ");
        out.print((int)Math.round(0.95*width));
        out.print(" ");
        out.println(842-((int)Math.round(0.05*height)));

        Graphics g = new PSGr2(out);                    
        
        out.println("0 setlinewidth");
        out.println("%%Pages:1");
        out.println("%%Page:1 1");                 
            
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

	double width = 596;
	double height = width*(maxY-minY)/(maxX-minX);
        
        int offset = 0;
        if(!b.title.equals(""))
            offset = offset + 10 + 3*b.fontsize/2;
        
        int hlp = offset;
        
        if(b.legend)
                hlp = hlp + 30 + 13*b.fontsize/10;

        if(height+hlp>842)
                {
                width = width*(842-hlp)/(height+hlp);
                height = height*(842-hlp)/(height+hlp);             
                }        
        
        out.println("%!PS-Adobe-3.0");      
        out.print("%%BoundingBox:");
        out.print((int)Math.round(0.05*width));
        out.print(" ");
        out.print(842-(int)Math.round(0.95*height+hlp+2));
        out.print(" ");
        out.print((int)Math.round(0.95*width));
        out.print(" ");
        out.println(842-((int)Math.round(0.05*height)));

        Graphics g = new PSGr2(out);                    
        
        out.println("0 setlinewidth");
        out.println("%%Pages:1");
        out.println("%%Page:1 1");                 
            
        helpdrawmap(g,height,width,offset);
        
        out.println("showpage");    
        out.println("grestore");
        out.println("gsave");        
	out.close();
       
	}         
        
private void helpdrawmap(Graphics g, double height, double width, int offset)
        {

        boolean NA = false;   
//        int nrNA = 0;
        b.nrNA = 0;
        
	double[] d2 = new double[4];
        double[] help = new double[4];
        double[] centroid = new double[2];

        Polygon p;

        double x;
        double y;

//        resize((int)width,(int)height);        
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
                        drawlegend(g,height,width);
                if(!b.title.equals(""))
                        drawtitle(g,height,width);
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
      
private void drawlegend(Graphics g,double height,double width)
        {
        
        int r;
        int end;
        double step;
        
        int offset = 0;
        if(!b.title.equals(""))
            offset = offset + 10 + 3*b.fontsize/2;
        
        double xoffset = 0.55*width;
        double yoffset = 0.95*height + offset + 10;        
        int legendheight = 20;
        
        Polygon p;    

        if(b.color)
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
        double max = b.upperlimit;
        double min = b.lowerlimit;
        
        g.setFont(new Font("TimesRoman", Font.PLAIN, b.fontsize));                                    
 
        str = String.valueOf(min);
        g.drawString(str,(int)Math.round(xoffset)-str.length()*9*b.fontsize/36,(int)Math.round(yoffset+legendheight+1.3*b.fontsize));
        g.drawLine((int)Math.round(xoffset),(int)Math.round(yoffset+legendheight),(int)Math.round(xoffset),(int)Math.round(yoffset+legendheight+3));
        if(min < 0 && 0 < max)                
            if(b.lowerlimit + (b.upperlimit-b.lowerlimit)/3 < 0 && 0 < b.upperlimit - (b.upperlimit-b.lowerlimit)/3)        
                {
                str = String.valueOf(0);
                g.drawString(str,(int)Math.round(xoffset-min/(max-min)*0.3*width)-str.length()*9*b.fontsize/36,(int)Math.round(yoffset+legendheight+1.3*b.fontsize));
                g.drawLine((int)Math.round(xoffset-min/(max-min)*0.3*width),(int)Math.round(yoffset+legendheight),
                           (int)Math.round(xoffset-min/(max-min)*0.3*width),(int)Math.round(yoffset+legendheight+3));
                }
        str = String.valueOf(max);
        g.drawString(str,(int)Math.round(xoffset+0.3*width)-str.length()*9*b.fontsize/36,(int)Math.round(yoffset+legendheight+1.3*b.fontsize));
        g.drawLine((int)Math.round(xoffset+0.3*width),(int)Math.round(yoffset+legendheight),(int)Math.round(xoffset+0.3*width),(int)Math.round(yoffset+legendheight+3));        
        
        }

private void drawNA(Graphics g, Polygon p)
        {

        int x,y;            
        int dist = 5;        
        Rectangle rect = p.getBounds();
/*
        x = rect.x;
        while(x < rect.x+rect.width)
                {
                y = rect.y;    
                while(y < rect.y+rect.height) 
                       {
                       if(p.contains(x,y,dist,dist))
                           g.drawLine(x,y,x+dist,y+dist);
                       y = y+dist;    
                       }
                x = x+dist;
                }                
*/        

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
        
private void drawtitle(Graphics g,double height,double width)
        {
        g.setFont(new Font("TimesRoman", Font.BOLD, 2*b.fontsize));            
        g.drawString(b.title,(int)(width/2-b.title.length()*9*b.fontsize/18),(int)(0.05*height)+3*b.fontsize/2);
        g.setFont(new Font("TimesRoman", Font.PLAIN, b.fontsize));            
        }

        
public void plotnonp(Graphics g)        
        {
        setplotparam((596-b.width)/2,120,b.width,b.height,b.pointsize,b.fontsize);
        plotframe(g,0);      
        for(int i=1;i<b.getDCols();i++)
            plot(g,i);
        }

public void Saveplotnonp(PrintWriter out)
        {
            
        setplotparam((596-b.width)/2,120,b.width,b.height,b.pointsize,b.fontsize);        
        
	int fontwidth = 11*fontsize/18;
	int fontheight = 3*fontsize/2;

        out.println("%!PS-Adobe-3.0");      
        out.print("%%BoundingBox:");
        out.print(xoffset-6*fontwidth-6);
        out.print(" ");
        if(b.xlab.equals(""))
            out.print(842-yoffset-height-fontheight-2);
        else
            out.print(842-yoffset-height-2*fontheight-2);
        out.print(" ");
        out.print(xoffset+width+2);
        out.print(" ");
        if(b.ylab.equals("") && b.title.equals(""))        
            out.println(842-yoffset+2);
        else if(!b.title.equals(""))         
            out.println(842-yoffset+9*fontheight/4);
        else    
            out.println(842-yoffset+3*fontheight/2);            
            
        PSGr2 g = new PSGr2(out);        

        scale = 10;        
        pointsize = pspointsize;   

        out.println("%%Pages:(atend)");
        out.println("%%Page:1 1");
        out.println("0.1 0.1 scale");
	out.println(b.linewidth);
        out.println(" setlinewidth");                
        
        plotframe(g,0);
        
        for(int i=1;i<b.getDCols();i++)
            plot(g,i);
       
        out.println("showpage");  
        out.close();
        
        }

private void plotframe(Graphics g, int col)
        {

        boolean date = false;
        
        String str;        
        
        int center = 0;
 
        int nrticks_x = 5;
        int nrticks_y = 5;

	int fontwidth = 11*fontsize/18;
	int fontheight = 3*fontsize/2;

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
                
//                if(nrticks_x*xstep < (maxX-xstart)*19/18)
//                    nrticks_x++;
                }
            if(b.ystep > 0.0)
                {
                ystep = b.ystep;
                nrticks_y = (int)((maxY-ystart)/ystep)+1;                                      

//                if(nrticks_y*ystep < (maxY-ystart)*19/18)
//                    nrticks_y++;
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
        g.drawRect(scale*xoffset, 842-scale*(842-yoffset), scale*width, scale*height);        
        g.setFont(new Font("Monospaced", Font.BOLD, scale*fontsize));                                                                 
        
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
                                str = formatLabel(str);                                                                
                                center = (str.length()+1)*fontwidth/2;     
                                g.drawString(str,translateX(i)-scale*center,842-scale*(842-(yoffset+height+17*fontheight/20)));                                                                                
                                }
                            if((j-1)%step == 0)
                                {
                                str = String.valueOf(year);
                                str = formatLabel(str);                                
                                center = (str.length()+1)*fontwidth/2;     
                                g.drawString(str,translateX(i)-scale*center,842-scale*(842-(yoffset+height+29*fontheight/20)));                                                                                
                                }
                            if(j%step == 0)
                                year++;
                            if(j == step)
                                j = 0;

                            g.drawLine(translateX(i), 842-scale*(842-yoffset-height), 
				       translateX(i), 842-scale*(842-yoffset-height-5));                                                                                    
                            
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
                                str = formatLabel(str);
                                center = (str.length()+1)*fontwidth/2;     
                                g.drawString(str,translateX(i)-scale*center,842-scale*(842-(yoffset+height+fontheight)));                                                                                
                                }
                            if((j-1)%step == 0)
                                g.drawLine(translateX(i), 842-scale*(842-yoffset-height), 
					   translateX(i), 842-scale*(842-yoffset-height-5));                                                                                    
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
                            str = formatLabel(str);
                            center = (str.length()+1)*fontwidth/2;
			    double pos = xstart+i*xstep;	
			    if(xstart <= maxX && pos >= minX && pos <= maxX*1.0000000001)
 				    {	
	                            g.drawString(str,translateX(pos)-scale*center,
				  	 842-scale*(842-(yoffset+height+fontheight)));
        	                    g.drawLine(translateX(pos), 842-scale*(842-yoffset-height), 
				       translateX(pos), 842-scale*(842-yoffset-height-5));                                
				    }	
                            }
                    }
                }       // end: function 3              
        else        
                {
                str = String.valueOf(minX);
                str = str.substring(0,str.indexOf('.'));                
                center = str.length()*fontwidth/2; 
                g.drawString(str,translateX(minX)-scale*center,842-scale*(842-(yoffset+height+fontheight)));
                for(int i=1;i<nrticks_x;i++)                    
                        {
                        str = String.valueOf(i*xstep);
                        str = str.substring(0,str.indexOf('.'));
                        center = str.length()*fontwidth/2; 
                        g.drawString(str,translateX(i*xstep)-scale*center,
					842-scale*(842-(yoffset+height+fontheight)));
                        g.drawLine(translateX(i*xstep), 842-scale*(842-yoffset-height),
				   translateX(i*xstep), 842-scale*(842-yoffset-height-5));    
                        }
                g.drawLine(translateX(minX), 842-scale*(842-yoffset-height),
			   translateX(minX), 842-scale*(842-yoffset-height-5));                    
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
                        center = (str.length()+1)*fontwidth+5;
                        g.drawString(str,scale*(xoffset-center),translateY(i*0.2)+scale*fontsize/3);
                        }
                }              
        else            
                {
                for(int i=0;i<nrticks_y;i++)                    
                        {
                        double value = (double)(Math.round((ystart+i*ystep)*1.0E10))/1.0E10;   
                        str = String.valueOf(value);                            
                        str = formatLabel(str);
                        center = (str.length()+1)*fontwidth+5;
			double pos = ystart+i*ystep;
			if(ystart <= maxY && pos >= minY && pos <= maxY*1.0000000001)
	                        g.drawString(str,scale*(xoffset-center),translateY(pos)+scale*fontsize/3); 
                        }
                }               
                
// title            
        if(b.function == 3)
                {
                if(!b.title.equals(""))
                        {
                        g.setFont(new Font("Monospaced", Font.BOLD, (int)(scale*fontsize*1.5)));                                        
                        str = b.title;
                        if(str.length()>32)
                          str = str.substring(0,32);
                        center = (int)((str.length()*1.5*fontwidth)/2);
                        g.drawString(str,scale*((int)(xoffset+width/2)-center),842-scale*(842-(yoffset-3*fontheight/2)));
                        g.setFont(new Font("Monospaced", Font.BOLD, scale*fontsize));                                                                                                           
                        }
                if(!b.xlab.equals(""))
                        {
                        str = b.xlab;
                        center = (int)((str.length()*fontwidth)/2);
                        g.drawString(str,scale*((int)(xoffset+width/2)-center),
					 842-scale*(842-(yoffset+height+2*fontheight)));
                        }
                if(!b.ylab.equals(""))
                        {
                        str = b.ylab;
                        center = (int)((str.length()*fontwidth)/2);
                        if(str.length()>14)
                            g.drawString(str,scale*(xoffset-7*fontwidth),842-scale*(842-(yoffset-3*fontheight/4)));
                        else
                            g.drawString(str,scale*(xoffset-center),842-scale*(842-(yoffset-3*fontheight/4)));
                        }
                }               
        else if(b.function == 4)
                {
                str = "iteration";
                center = (int)((str.length()*fontwidth)/2);
                g.drawString(str,scale*((int)(xoffset+width/2)-center),842-scale*(842-(yoffset+height+2*fontheight)));

                str = "par "+col;
                center = (int)((str.length()*fontwidth)/2);
                if(str.length()>14)
                    g.drawString(str,scale*(xoffset-5*fontwidth),842-scale*(842-(yoffset-fontheight)));
                else
                    g.drawString(str,scale*(xoffset-center),842-scale*(842-(yoffset-fontheight)));
                }
        else if(b.function == 5)
                {
                str = "lag";
                center = (int)((str.length()*fontwidth)/2);
                g.drawString(str,scale*((int)(xoffset+width/2)-center),842-scale*(842-(yoffset+height+2*fontheight)));

                str = b.getVarname(col);
                center = (int)((str.length()*fontwidth)/2);
                if(str.length()>10)
                    g.drawString(str,scale*(xoffset-5*fontwidth),842-scale*(842-(yoffset-fontheight)));                    
                else
                    g.drawString(str,scale*(xoffset-center),842-scale*(842-(yoffset-fontheight)));
                }
                
        }                

private int translateX(double x)
        { 
        return (int)Math.round( scale*(0.05*width + (x-minX)*0.9*width/(maxX-minX)) ) + scale*xoffset;        
        }
 
private int translateY(double y)
        { 
        return (int)Math.round( 842 - scale*(842-yoffset-height + 0.05*height + (y-minY)*0.9*height/(maxY-minY)) );        
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
                plotframe(g,start);
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
        out.print(842-100-2*100-3*130-42);
        out.print(" ");
        out.print(xoffset+2*width+66+5);
        out.print(" ");
        out.println(842-100+30);
            
        PSGr2 g = new PSGr2(out);        

        out.println("%%Pages:(atend)");

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

                        plotframe(g,start);
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
        out.print(842-100-2*100-3*130-42);
        out.print(" ");
        out.print(100+2*width+66+5);
        out.print(" ");
        out.println(842-100+30);
            
        PSGr2 g = new PSGr2(out);        

        out.println("%%Pages:(atend)");
                
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
                        
                        plotframe(g,start);
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
                
                plotframe(g,start);                
                plot(g,start);     
                drawLine(g);
                start++;
                count++;
                }
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
                if(b.connect.equals("points") || (b.connect.length()>=col&&b.connect.charAt(col-1)=='p'))
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

private String formatLabel(String str)
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
                // keine Formatierung f�r ganze Zahlen
                }

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
        
    
       
}        // END: class MapPanel