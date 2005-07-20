// This is the program where "model.txt" file is parsed
// and it is painted into a graphics content

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.io.*;

// The representation of a 3D model
class Model3D
{

	double xmin, xmax, ymin, ymax, zmin, zmax;
    double vert[];
    int tvert[],con[];
    int nvert, maxvert, ncon, maxcon, k1, k2, k3, gridsize;
    boolean transformed;
    char col;
    String Title, xlab, ylab, zlab;
    Matrix3D mat;


    Model3D ()
    {
	mat = new Matrix3D ();
	mat.xrot(20);
	mat.yrot(30);
    }

    // Create a 3D model by parsing an input stream
    Model3D (InputStream is) throws IOException, FileFormatException
    {
	this();
	          Reader r = new BufferedReader(new InputStreamReader(is));
    StreamTokenizer st = new StreamTokenizer(r);
	st.eolIsSignificant(true);

	// Fixing the color of the 3D plot
	st.nextToken();
	st.nextToken();
	if('r' == st.sval.charAt(0)) col = 'r';
	else if('g' == st.sval.charAt(0)) col = 'g';
	else if('b' == st.sval.charAt(0)) col = 'b';
	else if('y' == st.sval.charAt(0)) col = 'y';
	else if('c' == st.sval.charAt(0)) col = 'c';
	else if('m' == st.sval.charAt(0)) col = 'm';
	else if('o' == st.sval.charAt(0)) col = 'o';
	else if('G' == st.sval.charAt(0)) col = 'G';
	else col = 'B';

	st.nextToken();

	// Writing the title
	st.nextToken();
	st.nextToken();
	Title = " ";
	while(st.sval!= null)
	{
	Title = Title +" "+ st.sval;
	st.nextToken();
	}

    // Assigning the Axes Labels
    st.nextToken();
	st.nextToken();
	xlab = " ";
	while(st.sval!= null)
	{
	 xlab = xlab +" "+ st.sval;
	 st.nextToken();
	}

	st.nextToken();
	st.nextToken();
	ylab = " ";
	while(st.sval!= null)
	{
	 ylab = ylab +" "+ st.sval;
	 st.nextToken();
	}

	st.nextToken();
	st.nextToken();
	zlab = " ";
	while(st.sval!= null)
	{
	 zlab = zlab +" "+ st.sval;
	 st.nextToken();
	}


    // Assigning the gridsize and number of tick marks
    st.nextToken();
    st.nextToken();
          k1 = (int)st.nval;
    st.nextToken();
          k2 = (int)st.nval;
    st.nextToken();
          k3 = (int)st.nval;
    st.nextToken();
    gridsize = (int)st.nval;

    st.nextToken();


	// Parsing the remaining part of the input file
scan:
	while (true)
	{
	    switch (st.nextToken())
	    {
	    default:
		break scan;

	    case StreamTokenizer.TT_EOL:
		break;

	    case StreamTokenizer.TT_WORD:
	    if ("v".equals(st.sval))
		{
		    double x = 0, y = 0, z = 0;
		    if (st.nextToken() == StreamTokenizer.TT_NUMBER)
		    {
			 x = st.nval;
			 if (st.nextToken() == StreamTokenizer.TT_NUMBER)
			 {
			    y = st.nval;
			    if (st.nextToken() == StreamTokenizer.TT_NUMBER)
				z = st.nval;
			 }
		    }

		    addVert((double) x, (double) y, (double) z);
		    while (st.ttype != StreamTokenizer.TT_EOL &&
			    st.ttype != StreamTokenizer.TT_EOF)
			st.nextToken();

		}

		else if ("l".equals(st.sval))
		{
		    int start = -1;
		    int prev = -1;
		    int n = -1;
		    while (true)
			if (st.nextToken() == StreamTokenizer.TT_NUMBER)
			{
			    n = (int) st.nval;
			    if (prev >= 0)
				add(prev - 1, n - 1);
			    if (start < 0)
				start = n;
			    prev = n;
			}
			else if (st.ttype == '/')
			    st.nextToken();
			else
			    break;
		    if (start >= 0)
			add(start - 1, prev - 1);
		    if (st.ttype != StreamTokenizer.TT_EOL)
			break scan;
		}

		else
		{
		    while (st.nextToken() != StreamTokenizer.TT_EOL
			          && st.ttype != StreamTokenizer.TT_EOF);
		}
	    }
	}

	is.close();
	if (st.ttype != StreamTokenizer.TT_EOF)
	    throw new FileFormatException(st.toString());
    }

    // Adding a vertex to this model
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

    // Adding a line from vertex p1 to vertex p2
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

    // Transform all the points in this model
    void transform()
    {
	 if (transformed || nvert <= 0)
	    return;
	 if (tvert == null || tvert.length < nvert * 3)
	    tvert = new int[nvert*3];
	 mat.transform(vert, tvert, nvert);
	 transformed = true;
	 }

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

    // Eliminating duplicate lines
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
	   with this model to map from model space to screen space*/

    void paint(Graphics g)
    {
    if (vert == null || nvert <= 0)
			    return;

	transform();

	if (gr == null)
	{
	    gr = new Color[16];
	    for (int i = 0; i < 16; i++)
	    {
		int grey = (int) (Math.pow(3*(1-(i/16)),5));
		if (col == 'B') gr[i] = new Color(grey, grey, grey);
		if (col == 'G') gr[i] = new Color( 3*grey/4, 3*grey/4, 3*grey/4);
		if (col == 'r') gr[i] = new Color(grey, 0, 0);
		if (col == 'g') gr[i] = new Color(0, grey, 0);
		if (col == 'b') gr[i] = new Color(0, 0, grey);
		if (col == 'y') gr[i] = new Color(grey, grey, 0);
		if (col == 'c') gr[i] = new Color(0, grey, grey);
	    if (col == 'm') gr[i] = new Color(grey, 0, grey);
    	if (col == 'o') gr[i] = new Color(grey, grey/3, 0);
        }
	}

	int lg = 0;
	int lim = ncon;
	int c[] = con;
	int v[] = tvert;

	if (lim <= 0 || nvert <= 0)
	return;

    // Drawing the main plot
	for (int i = 0; i < lim-(12+k1+k2+k3); i++)
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
		g.setColor(gr[grey]);
		}
	    g.drawLine(v[p1], v[p1 + 1], v[p2], v[p2 + 1]);
	}

    // Drawing the axes and writing the axes-labels
	g.setColor(Color.BLACK);
	int temp = 0;
	for (int i = lim-(12+k1+k2+k3); i < lim-(9+k1+k2+k3); i++)
	{
	    int T1 = c[i];
	    int p1 = ((T1 >> 16) & 0xFFFF) * 3;
	    int p2 = (T1 & 0xFFFF) * 3;
	    g.drawLine(v[p1], v[p1 + 1], v[p2], v[p2 + 1]);
	    temp = temp+1;
		if (temp == 1)
		{
		 if (xlab.length() == 0)	g.drawString("X" , v[p2]+15,v[p2+1]+10);
		 else 	g.drawString(xlab , v[p2]+15,v[p2+1]+10);
		}
		else if (temp == 2)
		{
		 if (ylab.length() == 0)	g.drawString("Y" , v[p2]+15,v[p2+1]+10);
		 else 	g.drawString(ylab , v[p2]+15,v[p2+1]+10);
		}
		else
		{
		 if (zlab.length() == 0)	g.drawString("Z" , v[p2]+15,v[p2+1]+10);
		 else 	g.drawString(zlab , v[p2]+15,v[p2+1]+10);
		}
		g.drawString("O", v[p1]+5,v[p1+1]+5);
	}

	// Drawing the arrows at the end of the axes
	for (int i = lim-(9+k1+k2+k3); i < lim-(k1+k2+k3+3); i++)
	{
	    int T = c[i];
	    int p1 = ((T >> 16) & 0xFFFF) * 3;
	    int p2 = (T & 0xFFFF) * 3;
	    g.setColor(Color.BLACK);
	    g.drawLine(v[p1], v[p1 + 1], v[p2], v[p2 + 1]);
   }

   //Drawing the tick marks
	    for (int i = lim-(k1+k2+k3+3); i < lim; i++)
    	{
		int T = c[i];
	    int p1 = ((T >> 16) & 0xFFFF) * 3;
	    int p2 = (T & 0xFFFF) * 3;
	    g.setColor(Color.BLACK);
		g.drawLine(v[p1], v[p1 + 1], v[p2], v[p2 + 1]);
        }

	   Font font = new Font("Franklin Gothic Medium", Font.BOLD, 25);
	   g.setFont(font);
	   g.drawString(Title, 3, 40);
	}

    // Find the bounding box of this model
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
    }
}

