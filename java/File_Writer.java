/** Given a point and line set this program writes a file in .obj format.
  * This is the standard file format for many other plotting softwares*/

import java.lang.System;
import java.io.*;

class File_Writer
{
 String t,xlab,ylab,zlab;
 int k1,k2,k3,g;
 double poi [][];
 int lin [][];
 char s;

/**First of the two 2D arrays in the argument here is the 3D point data set
 * whereas the other gives information about which points are connected by a
 * line*/
  File_Writer(double[][] poi_1, int[][] lin_1,char s1, String t1, String x1, String y1, String z1, int k11, int k21, int k31,int g1)
  {
   poi = poi_1;
   lin = lin_1;
     s = s1;
     t = t1;
  xlab = x1;
  ylab = y1;
  zlab = z1;
    k1 = k11;
    k2 = k21;
    k3 = k31;
     g = g1;
  }

 /**This method gives a text file in the standard .obj format as output*/
  public void writing() throws IOException
  {
   try
   {
      BufferedWriter out = new BufferedWriter( new FileWriter("model.txt") );

      /**For writing vertex into the output file*/
       out.write('#');
       out.write(' ');
       out.write(s);
       out.newLine();
       out.write('#');
       out.write(' ');
       out.write(t);
       out.newLine();
       out.write('#');
	   out.write(' ');
	   out.write(xlab);
	   out.newLine();
       out.write('#');
	   out.write(' ');
	   out.write(ylab);
	   out.newLine();
       out.write('#');
	   out.write(' ');
	   out.write(zlab);
	   out.newLine();
       out.write('#');
       out.write(' ');
       out.write(Integer.toString(k1));
       out.write(' ');
       out.write(Integer.toString(k2));
       out.write(' ');
       out.write(Integer.toString(k3));
       out.write(' ');
       out.write(Integer.toString(g));
       out.newLine();

        for(int i=0;i<poi.length;++i)
       {
        out.write('v');
        for(int j=0;j<3;j++)
        {
         out.write(' ');
         if(Math.abs(poi[i][j]) < 0.001) poi[i][j] = 0;
         out.write(Double.toString(poi[i][j]));
        }
        out.newLine();
       }

       /**For writing information about lines in the output files*/
       for(int i=0;i<lin.length;++i)
       {
        out.write('l');
        for(int j=0;j<2;j++)
        {
         out.write(' ');
         out.write(Integer.toString(lin[i][j]));
        }
        out.newLine();
       }

       out.close();
   }
   catch (IOException e)
   {
   }

  }

}
