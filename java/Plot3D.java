/** This program constructs a 3D surface from the scattered point data entered as an array.*/

import java.awt.*;
import java.awt.geom.*;
import java.io.*;
import java.awt.*;

public class Plot3D
{
  double[][] p;
  int[][] linf;
  String mdname, xlab, ylab, zlab, title;
  char c;
  int gridsize,k1,k2,k3;
  double xstart,xstep,ystart,ystep,zstart,zstep,angx,angy,angz;
  InputStream is;


 /** The constructor here takes the 3D data set, gridsize,the color of the plot, title, axes labels,
  *  starting points and length of tick marks and the angle of rotation about each of the three axes as arguments*/
    public Plot3D(double[][] p_1, int g1, char c1, String title1, String xlab1, String ylab1, String zlab1,double x1,double xs,double y1,double ys,double z1,double zs,double angx1, double angy1, double angz1)
    {
     Dense d1 = new Dense(p_1,g1,x1,xs,y1,ys,z1,zs);
     gridsize = g1;
            p = d1.densegen();
            c = c1;
        title = title1;
         xlab = xlab1;
         ylab = ylab1;
         zlab = zlab1;
       xstart = x1;
        xstep = xs;
       ystart = y1;
        ystep = ys;
       zstart = z1;
        zstep = zs;
         angx = angx1;
         angy = angy1;
         angz = angz1;
    }

    void Plot(Graphics g)
     {
      int l1 = (gridsize*gridsize*3)+(2*gridsize)+9;
      int l2 = (gridsize+1)*(gridsize+1);
      int l3 = (gridsize*gridsize*3)+(2*gridsize);

      int[][] lin = new int[l1][2];
            int f = 0;

     for (int i=0; i<gridsize; i++)
     {
       for (int j=0; j<gridsize; j++)
       {
       lin[f][0] = (gridsize+1)*i+j+1;
       lin[f][1] = lin[f][0]+1;
               f = f+1;
       lin[f][0] = (gridsize+1)*i+j+1;
       lin[f][1] = lin[f][0]+gridsize+1;
               f = f+1;
       lin[f][0] = (gridsize+1)*i+j+1;
       lin[f][1] = lin[f][0]+gridsize+2;
               f = f+1;
       }
     }

     for (int j=0; j<gridsize; j++)
     {
      lin[f][0] = (gridsize+1)*j+(gridsize+1);
      lin[f][1] = (gridsize+1)*j+2*(gridsize+1);
              f = f+1;
     }

     for (int j=0; j<gridsize; j++)
     {
      lin[f][0] = (gridsize+1)*gridsize+j+1;
      lin[f][1] = lin[f][0]+1;
              f = f+1;
     }


   // Adding lines to draw axes and arrows at the ends of each axis denoting the positive side
     lin[l3][0] = l2+1;
     lin[l3][1] = l2+2;
     lin[l3+1][0] = l2+1;
     lin[l3+1][1] = l2+3;
     lin[l3+2][0] = l2+1;
     lin[l3+2][1] = l2+4;
     lin[l3+3][0] = l2+2;
     lin[l3+3][1] = l2+5;
     lin[l3+4][0] = l2+2;
     lin[l3+4][1] = l2+6;
     lin[l3+5][0] = l2+3;
     lin[l3+5][1] = l2+7;
     lin[l3+6][0] = l2+3;
     lin[l3+6][1] = l2+8;
     lin[l3+7][0] = l2+4;
     lin[l3+7][1] = l2+9;
     lin[l3+8][0] = l2+4;
     lin[l3+8][1] = l2+10;

   // Adding lines to draw tick marks on each axis
     double[] xcord = new double[p.length];
	 double[] ycord = new double[p.length];
	 double[] zcord = new double[p.length];

	 for(int i=0; i<p.length; i++)
	 {
	   xcord[i] = p[i][0];
	   ycord[i] = p[i][1];
	   zcord[i] = p[i][2];
	 }

     double[]temp1 = new double[p.length];
	       Sort s1 = new Sort(xcord);
	         temp1 = s1.mergeSort();

	 double[]temp2 = new double[p.length];
	       Sort s2 = new Sort(ycord);
	         temp2 = s2.mergeSort();

	 double[]temp3 = new double[p.length];
	       Sort s3 = new Sort(zcord);
	         temp3 = s3.mergeSort();

     if(xstep == 0|| ystep == 0||zstep == 0)
     {
	  linf = new int[l1+15][2];
	  for(int i= 0; i<l1; i++)
	  {
	   linf[i] = lin[i];
	  }

      for(int i= l1; i<l1+15; i++)
	  {
	   linf[i][0] = l2+10+2*(i-l1+1)-1;
	   linf[i][1] = linf[i][0]+1;
	  }
     }

     else
     {
	   k1 = (int)((p[l2-1][0]-xstart)/xstep);
	   k2 = (int)((p[l2-1][1]-ystart)/ystep);
	   k3 = (int)(((temp3[p.length-1]/1.4)-zstart)/zstep);

	  linf = new int[l1+k1+k2+k3+3][2];
	  for(int i= 0; i<l1; i++)
	  {
	   linf[i] = lin[i];
	  }

	  for(int i= l1; i<linf.length; i++)
	  {
	   linf[i][0] = l2+10+2*(i-l1+1)-1;
	   linf[i][1] = linf[i][0]+1;
	  }
	 }

    // Implementing the main program
     Surfplot s = new Surfplot(p,linf,gridsize,c,title,xlab,ylab,zlab,xstart,xstep,ystart,ystep,zstart,zstep,angx,angy,angz);
	 s.run();
	 s.paint(g);


    // Adding few additional things for the applet part of the program
    double[] scale = new double[3];
   double[] scalet = new double[3];
          scale[0] = Math.max(Math.abs(temp1[0]), Math.abs(temp1[p.length-1]));
          scale[1] = Math.max(Math.abs(temp2[0]), Math.abs(temp2[p.length-1]));
          scale[2] = Math.max(Math.abs(temp3[0]), Math.abs(temp3[p.length-1]));

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

    System.out.println ("The scaling factors corresponding to each axis in the plot are as following:");
    System.out.println ("X-Axis: " + " " +scalet[0]);
    System.out.println ("Y-Axis: " + " " +scalet[1]);
    System.out.println ("Z-Axis: " + " " +scalet[2]);

    for(int i=0; i<p.length; i++)
    {
      for(int j=0; j<3; j++)
      {
       p[i][j] = scalet[j]*p[i][j];
      }
    }

    if (linf.length == l1+15)
    {
	 k1 = 4;
	 k2 = 4;
	 k3 = 4;
	}

    File_Writer fil1= new File_Writer(p,linf,c,title,xlab,ylab,zlab,k1,k2,k3,gridsize);
    try
	{
	 fil1.writing();
	}
	catch (IOException e)
	{
    }
   }
}
