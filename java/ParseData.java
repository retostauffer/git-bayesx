/** Similar to the Model3D.java program.The difference is that here we directly need to
  * enter point and line array to form the model*/

import java.awt.*;
import java.awt.geom.*;
import java.util.*;
import java.lang.Comparable.*;

class ParseData
{
    double p[][], p_org[][];
    double vert[];
    double xmin, xmax, ymin, ymax, zmin, zmax;
    int lin[][];
    int tvert[],con[];
    int nvert, maxvert,ncon, maxcon, gridsize, cflag, k1, k2, k3;
    boolean transformed;
    String Title, xlab, ylab, zlab;
    char col;
    Matrix3D mat;


    ParseData (double[][] p_1, int[][] lin1, int g1, char col1, String t1, String x1, String y1, String z1,double xst,double xs,double yst,double ys,double zst,double zs)
    {

        	 mat = new Matrix3D();
	         mat.xrot(20);
             mat.yrot(20);
           p_org = p_1;
               p = new double[p_org.length][3];
        	 lin = lin1;
        gridsize = g1;
	         col = col1;
           Title = t1;
            xlab = x1;
            ylab = y1;
            zlab = z1;

       double[] xcord = new double[p_org.length];
	   double[] ycord = new double[p_org.length];
	   double[] zcord = new double[p_org.length];

	   for(int i=0; i<p_org.length; i++)
	   {
	    xcord[i] = p_org[i][0];
	    ycord[i] = p_org[i][1];
	    zcord[i] = p_org[i][2];
	   }

	   double[]temp1 = new double[p_org.length];
	         Sort s1 = new Sort(xcord);
	           temp1 = s1.mergeSort();

	   double[]temp2 = new double[p_org.length];
	         Sort s2 = new Sort(ycord);
	           temp2 = s2.mergeSort();

	   double[]temp3 = new double[p_org.length];
	         Sort s3 = new Sort(zcord);
	           temp3 = s3.mergeSort();

       if(xs==0||ys==0||zs==0)
	   {
	   	 k1 = 4;
	     k2 = 4;
	     k3 = 4;
	   }

       else
       {
	     k1 = (int)((p_org[(gridsize+1)*(gridsize+1)-1][0] - xst) / xs);
         k2 = (int)((p_org[(gridsize+1)*(gridsize+1)-1][1] - yst) / ys);
         k3 = (int)((temp3[p_org.length-1]/1.4 - zst) / zs);
	   }

      double[] scale = new double[3];
	 double[] scalet = new double[3];
	        scale[0] = temp1[p_org.length-1]- temp1[0];
	        scale[1] = temp2[p_org.length-1]- temp2[0];
	        scale[2] = temp3[p_org.length-1]- temp3[0];

	  for(int i=0; i<3; i++)
	  {
	   scalet[i] = scale[i];
	  }
	         Sort s4 = new Sort(scale);
	   double[]temp4 = new double[3];
	           temp4 = s4.mergeSort();
	  double scale_f = temp4[2];

	  for(int i=0; i<3; i++)
	  {
	   if((scale_f/scalet[i]-(int)(scale_f/scalet[i])) >= 0.5) scalet[i] = (int)(scale_f/scalet[i])+1;
	   else scalet[i] = (int)(scale_f/scalet[i]);
      }

      for(int i=0; i<p.length; i++)
	  {
	    for(int j=0; j<3; j++)
	    {
	      p[i][j] = scalet[j]*p_org[i][j];
	    }
      }


   // Parsing the point array
	  for (int i=0; i<p.length; i++)
	  {
	   addVert(p[i][0], p[i][1], p[i][2]);
	  }

   // Parsing the lines array
 	  for (int i=0; i<lin.length; i++)
	  {
	   add(lin[i][0]-1, lin[i][1]-1);
	  }

     }


  // Adding a vertex to the model
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


    static Color gr[];

 // Paint this model to a graphics context.It uses the matrix associated
 //	with this model to map from model space to screen space
    void paintc(Graphics g)
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
			g.drawLine(v[p1], v[p1 + 1],v[p2], v[p2 + 1]);
		}

        // Drawing the axes and writing the axes-labels
	    int temp = 0;
		for (int i = lim-(12+k1+k2+k3); i < lim-(9+k1+k2+k3); i++)
		{
		    int T1 = c[i];
		    int p1 = ((T1 >> 16) & 0xFFFF) * 3;
		    int p2 = (T1 & 0xFFFF) * 3;
		    g.setColor(Color.BLACK);
		    g.drawLine(v[p1], v[p1 + 1], v[p2], v[p2 + 1]);
		    temp = temp+1;
			if (temp == 1)
			{
			 if (xlab.length() == 0)	g.drawString("" , v[p2]+15,v[p2+1]+10);
			 else 	g.drawString(xlab , v[p2]+25,v[p2+1]+10);
			}
			else if (temp == 2)
			{
			 if (ylab.length() == 0)	g.drawString("" , v[p2]+15,v[p2+1]+10);
			 else 	g.drawString(ylab , v[p2]+15,v[p2+1]+15);
			}
			else
			{
			 if (zlab.length() == 0)	g.drawString("" , v[p2]+15,v[p2+1]+10);
			 else 	g.drawString(zlab , v[p2]+15,v[p2+1]+15);
			}
		}

        // Drawing arrows at the end of each axis,denoting positive direction
		for (int i = lim-(9+k1+k2+k3); i < lim-(k1+k2+k3+3); i++)
		{
		    int T = c[i];
		    int p1 = ((T >> 16) & 0xFFFF) * 3;
		    int p2 = (T & 0xFFFF) * 3;
		    g.setColor(Color.BLACK);
		    g.drawLine(v[p1], v[p1 + 1], v[p2], v[p2 + 1]);
        }

        // Drawing tick marks and writing the labels
        int flag1 = 0;
        for (int i = lim-(k1+k2+k3+3); i < lim; i++)
		{
			int T = c[i];
		    int p1 = ((T >> 16) & 0xFFFF) * 3;
		    int p2 = (T & 0xFFFF) * 3;
		    g.setColor(Color.BLACK);
			g.drawLine(v[p1], v[p1 + 1], v[p2], v[p2 + 1]);

			if (0<= flag1 && flag1<k1+1)
			{
				 int flag2 = (Double.toString(p_org[((gridsize+1)*(gridsize+1)+10)+2*flag1][0])).indexOf('.');
				 char[] tempc = new char[ Math.min( flag2+3, (Double.toString(p_org[((gridsize+1)*(gridsize+1)+10)+2*flag1][0])).length()) ];
				 for(int j=0; j<tempc.length; j++)
				 {
				  tempc[j] = (Double.toString(p_org[((gridsize+1)*(gridsize+1)+10)+2*flag1][0])).charAt(j);
				 }

				 String stemp = new String(tempc);
				 g.drawString(stemp, v[p2]+15, v[p2+1]+10);
			}

			if (k1+1<= flag1 && flag1<k1+k2+2)
			{
				 int flag2 = (Double.toString(p_org[((gridsize+1)*(gridsize+1)+10)+2*flag1][1])).indexOf('.');
				 char[] tempc = new char[ Math.min( flag2+3, (Double.toString(p_org[((gridsize+1)*(gridsize+1)+10)+2*flag1][1])).length()) ];
				 for(int j=0; j<tempc.length; j++)
				 {
				  tempc[j] = (Double.toString(p_org[((gridsize+1)*(gridsize+1)+10)+2*flag1][1])).charAt(j);
				 }

                 String stemp = new String(tempc);
				 g.drawString(stemp, v[p2], v[p2+1]+10);
			}

			if (k1+k2+2<= flag1 && flag1<k1+k2+k3+3)
			{
				 int flag2 = (Double.toString(p_org[((gridsize+1)*(gridsize+1)+10)+2*flag1][2])).indexOf('.');
				 char[] tempc = new char[ Math.min( flag2+3, (Double.toString(p_org[((gridsize+1)*(gridsize+1)+10)+2*flag1][2])).length()) ];
				 for(int j=0; j<tempc.length; j++)
				 {
				  tempc[j] = (Double.toString(p_org[((gridsize+1)*(gridsize+1)+10)+2*flag1][2])).charAt(j);
				 }

				 String stemp = new String(tempc);
				 g.drawString(stemp, v[p2]+10, v[p2+1]+10);
			}

			flag1 = flag1+1;
        }

        Font font = new Font("TimesRoman", Font.BOLD, 25);
	    g.setFont(font);
	    g.drawString(Title,250-(Title.length()/2),40);
	}



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
   	 this.xmin = 1.2*xmin;
   	 this.ymax = ymax;
   	 this.ymin = 1.2*ymin;
   	 this.zmax = zmax;
   	 this.zmin = 1.2*zmin;
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


}