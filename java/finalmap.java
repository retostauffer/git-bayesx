import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.util.*;
import java.io.*;

class finalmap 
{
    double p[][];
    double vert[];
    double xmin, xmax, ymin, ymax, zmin, zmax;
    int lin[][];
    int tvert[];
    int nvert, maxvert;
    int con[];
    int ncon, maxcon,cflag;
    boolean transformed;
    String Title;
    char c;
    Matrix3D mat;
    
    finalmap()
    {
	 mat = new Matrix3D ();
	 mat.xrot(20);
	 mat.yrot(30);
    }

    /** Constructor for the program*/ 
    finalmap (double[][] p1, int[][] lin1, char c1, String t1) 
    {
	     p = p1;
	   lin = lin1;
	     c = c1;
     Title = t1;
	 
	 // Parsing the point array
	 for (int i=0; i<p.length; i++)
	 {
	  addVert(p[i][0], p[i][1], p[i][2]);
	 }
 
     // Parsing the lines array
 	 for (int i=0; i<p.length; i++)
	 {
	  add(lin[i][0], lin[i][1]);
	 }	
     
//     transform();
     compress();
    }
   
   
   /** Adding a vertex to this model */
    int addVert(double x, double y, double z)
    {
	 int i = nvert;
	 if (i >= maxvert)
	    if (vert == null)
	    {
		 maxvert = 100;
		 vert = new double[maxvert * 3];
	    }
	    else
	    {
		 maxvert *= 2;
		 double nv[] = new double[maxvert * 3];
		 System.arraycopy(vert, 0, nv, 0, vert.length);
		 vert = nv;
	    }
	i *= 3;
	vert[i] = x;
	vert[i + 1] = y;
	vert[i + 2] = z;
	return nvert++;
    }


    /** Add a line from vertex p1 to vertex p2 */
    void add(int p1, int p2)
    {
	 int i = ncon;
	 if (p1 >= nvert || p2 >= nvert)
	     return;
	 if (i >= maxcon)
	    if (con == null)
	    {
		 maxcon = 100;
		 con = new int[maxcon];
	    }
	    else
	    {
		 maxcon *= 2;
		 int nv[] = new int[maxcon];
		 System.arraycopy(con, 0, nv, 0, con.length);
		 con = nv;
	    }
   	 if (p1 > p2)
   	 {
	    int t = p1;
	    p1 = p2;
	    p2 = t;
	 }
	 con[i] = (p1 << 16) | p2;
	 ncon = i + 1;
    }


    /** Transform all the points in this model */
   /** void transform()
    {
	 if (transformed || nvert <= 0)
	    return;
	 if (tvert == null || tvert.length < nvert * 3)
	    tvert = new int[nvert*3];
	 mat.transform(vert, tvert, nvert);
	 transformed = true;
	 }*/


    private void sort(int lo0, int hi0)
    {
	 int a[] = con;
	  int lo = lo0;
	  int hi = hi0;
	 if (lo >= hi)
	     return;
	 int mid = a[(lo + hi) / 2];
	 while (lo < hi) 
	 {
	    while (lo < hi && a[lo] < mid) 
	    {
		lo++;
	    }
	    while (lo < hi && a[hi] >= mid) 
	    {
		hi--;
	    }
	    if (lo < hi) 
	    {
		int T = a[lo];
		a[lo] = a[hi];
		a[hi] = T;
	    }
	}
	
	if (hi < lo) 
	{
	    int T = hi;
	    hi = lo;
	    lo = T;
	}
	sort(lo0, lo);
	sort(lo == lo0 ? lo + 1 : lo, hi0);
    }


    /** Eliminating duplicate lines */
    void compress() 
    {
	int limit = ncon;
	int c[] = con;
	sort(0, ncon - 1);
	int d = 0;
	int pp1 = -1;
	for (int i = 0; i < limit; i++) 
	{
	    int p1 = c[i];
	    if (pp1 != p1) 
	    {
		c[d] = p1;
		d++;
	    }
	    pp1 = p1;
	}
	ncon = d;
    }

    static Color gr[];

    /** Paint this model to a graphics context.It uses the matrix associated
	*   with this model to map from model space to screen space*/

    void paintcomp(Graphics g) 
    {
    Graphics2D g2 = (Graphics2D) g;    
	if (vert == null || nvert <= 0)
	    return;
//	transform();
	if (gr == null) 
	{
	    gr = new Color[16];
	    for (int i = 0; i < 16; i++) 
	    {
		int grey = (int) (Math.pow(3*(1-(i/16)),5));
		
	    //Fixing the color of the 3D plot
	    if (c == 'r') gr[i] = new Color(grey, 0, 0);
		if (c == 'g') gr[i] = new Color(0, grey, 0);
		else gr[i] = new Color(0, 0, grey);
	    }
	}
	
	int lg = 0;
	int lim = ncon;
	int c[] = con;
	int v[] = tvert;
	if (lim <= 0 || nvert <= 0)
	    return;
	for (int i = 0; i < lim-6; i++) 
	{
	    int T = c[i];
	    int p1 = ((T >> 16) & 0xFFFF) * 3;
	    int p2 = (T & 0xFFFF) * 3;
	    int grey = v[p1 + 2] + v[p2 + 2];
	    if (grey < 0)
		grey = 0;
	    if (grey > 15)
		grey = 15;
	    if (grey != lg) 
	    {
		lg = grey;
		g2.setColor(gr[grey]);
	    }
	    g2.drawLine(v[p1], v[p1 + 1],
		       v[p2], v[p2 + 1]);
	}
	
	for (int i = lim-12; i < lim-9; i++) 
	{
	    int T = c[i];
	    int p1 = ((T >> 16) & 0xFFFF) * 3;
	    int p2 = (T & 0xFFFF) * 3;
	    g2.setColor(Color.BLACK);
	    g2.drawLine(v[p1], v[p1 + 1],
		       v[p2], v[p2 + 1]);
    }
	
	 int temp = 0;
	int temp1 = 0; 
	int temp2 = 0;
	for (int i = lim-9; i < lim-6; i++) 
	{
	    int T1 = c[i];
	    int p1 = ((T1 >> 16) & 0xFFFF) * 3;
	    int p2 = (T1 & 0xFFFF) * 3;
	    g2.setColor(Color.BLACK);
	    g2.drawLine(v[p1], v[p1 + 1],
		       v[p2], v[p2 + 1]);
	    temp = temp+1;	       
		if (temp == 1)g2.drawString("X", v[p2]+15,v[p2+1]+15);
		else if (temp == 2)g2.drawString("Y", v[p2]+15,v[p2+1]+15);
		else g2.drawString("Z", v[p2]+15,v[p2+1]+15);
		g2.drawString("0", v[p1]+5,v[p1+1]+5);
	}
	
	for (int i = lim-6; i < lim; i++) 
	{
	    int T = c[i];
	    int p1 = ((T >> 16) & 0xFFFF) * 3;
	    int p2 = (T & 0xFFFF) * 3;
	    g2.setColor(Color.BLACK);
	    g2.drawLine(v[p1], v[p1 + 1],
		       v[p2], v[p2 + 1]);
    }
    
    Font font = new Font("Franklin Gothic Medium", Font.BOLD, 25);
    g2.setFont(font);
    g2.drawString(Title, 3, 40);
	}

    /** Find the bounding box of this model 
    void findBB() 
    {
	if (nvert <= 0)
	    return;
	 double v[] = vert;
	double xmin = v[0], xmax = xmin;
	double ymin = v[1], ymax = ymin;
	double zmin = v[2], zmax = zmin;
	for (int i = nvert * 3; (i -= 3) > 0;)
	{
	    double x = v[i];
	    if (x < xmin)
		xmin = x;
	    if (x > xmax)
		xmax = x;
	    double y = v[i + 1];
	    if (y < ymin)
		ymin = y;
	    if (y > ymax)
		ymax = y;
	    double z = v[i + 2];
	    if (z < zmin)
		zmin = z;
	    if (z > zmax)
		zmax = z;
	}
	this.xmax = xmax;
	this.xmin = xmin;
	this.ymax = ymax;
	this.ymin = ymin;
	this.zmax = zmax;
	this.zmin = zmin;
    }*/
     
}

