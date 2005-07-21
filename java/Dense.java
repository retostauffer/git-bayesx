/** This program creates a grid of size gridsize x gridsize and evaluates
  * the function value at new points using Shepherd interpolation technique.
  */

class Dense
{
  double[][] p, finalm;
  double xstart, xstep, ystart, ystep, zstart, zstep;
  int gridsize;

    Dense(double[][] p_1, int g1, double x1, double xs, double y1, double ys, double z1, double zs)
    {
           p = p_1;
    gridsize = g1;
      xstart = x1;
       xstep = xs;
      ystart = y1;
       ystep = ys;
      zstart = z1;
       zstep = zs;
    }

    double[][] densegen()
    {
     double[] xcord = new double[p.length];
     double[] ycord = new double[p.length];
     double[] zcord = new double[p.length];
     for(int i=0; i<p.length; i++)
     {
      xcord[i] = p[i][0];
      ycord[i] = p[i][1];
      zcord[i] = p[i][2];
     }

     double[]temp1 = new double[xcord.length];
           Sort s1 = new Sort(xcord);
             temp1 = s1.mergeSort();

     double[]temp2 = new double[ycord.length];
           Sort s2 = new Sort(ycord);
             temp2 = s2.mergeSort();

     double[]temp3 = new double[zcord.length];
           Sort s3 = new Sort(zcord);
             temp3 = s3.mergeSort();



    double[] scale = new double[3];
   double[] scalet = new double[3];
          scale[0] = temp1[p.length-1]- temp1[0];
		  scale[1] = temp2[p.length-1]- temp2[0];
	      scale[2] = temp3[p.length-1]- temp3[0];

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

           double h1 =(temp1[p.length-1] - temp1[0] )/gridsize;
           double h2 =(temp2[p.length-1] - temp2[0] )/gridsize;

      double[] meshx = new double[gridsize+1];
      double[] meshy = new double[gridsize+1];


     // Generating the gird on the X-Y plane
     for(int i=0; i<gridsize+1; i++)
      {
       meshx [i] = temp1[0] + i*h1;
      }

      for(int i=0; i<gridsize+1; i++)
      {
       meshy [i] = temp2[0] + i*h2;
      }

      int l1 = (gridsize+1)*(gridsize+1);
      double final1[][] = new double [l1][2];
      for(int i=0; i<gridsize+1; i++)
      {
       for(int j=0; j<gridsize+1; j++)
       {
        final1[i*(gridsize+1)+j][0] = meshx[i];
        final1[i*(gridsize+1)+j][1] = meshy[j];
       }
      }

      int l2 = l1+10;
      double[][] final2 = new double[l2][3];

      for(int i=0; i<l1; i++)
      {
       Interpolate i1 = new Interpolate(p);                //Finding the value at the particular point by interpolation
         final2[i][0] = final1[i][0];
         final2[i][1] = final1[i][1];
         final2[i][2] = i1.zval(final1[i]);
      }


      // Adding points for drawing axes on the plot
         double xax,yax,zax;
         if(temp1[0] >0) xax = 0.9*temp1[0];
         else xax = 1.1*temp1[0];
         if(temp2[0] >0) yax = 0.9*temp2[0];
         else yax = 1.1*temp2[0];
         if(temp3[0] >0) zax = 0.9*temp3[0];
         else zax = 1.1*temp3[0];

         final2[l1][0] = xax;
		 final2[l1][1] = yax;
		 final2[l1][2] = zax;

         final2[l1+1][0] = 1.2*temp1[p.length-1];
         final2[l1+1][1] = yax;
         final2[l1+1][2] = zax;

         final2[l1+2][0] = xax;
		 final2[l1+2][1] = 1.2*temp2[p.length-1];
         final2[l1+2][2] = zax;

         final2[l1+3][0] = xax;
		 final2[l1+3][1] = yax;
		 final2[l1+3][2] = 1.2*temp3[p.length-1];


     //  Adding points for drawing arrows at the end of the axes
               double h3 = (temp1[p.length-1] - temp1[0])*0.0075;
               double h4 = (temp2[p.length-1] - temp2[0])*0.0075;
               double h5 = (temp3[p.length-1] - temp3[0])*0.0075;
         final2[l1+4][0] = 1.2*temp1[p.length-1]-(8*h3);
         final2[l1+4][1] = yax+2*h4;
         final2[l1+4][2] = zax;
         final2[l1+5][0] = 1.2*temp1[p.length-1]-(8*h3);
         final2[l1+5][1] = yax-2*h4;
         final2[l1+5][2] = zax;

         final2[l1+6][0] = xax+2*h3;
         final2[l1+6][1] = 1.2*temp2[p.length-1]-(8*h4);
         final2[l1+6][2] = zax;
         final2[l1+7][0] = xax-2*h3;
         final2[l1+7][1] = 1.2*temp2[p.length-1]-(8*h4);
         final2[l1+7][2] = zax;

         final2[l1+8][0] = xax+2*h3;
         final2[l1+8][1] = yax;
         final2[l1+8][2] = 1.2*temp3[p.length-1]-(8*h5);
         final2[l1+9][0] = xax-2*h3;
         final2[l1+9][1] = yax;
         final2[l1+9][2] = 1.2*temp3[p.length-1]-(8*h5);


     //  Adding points for drawing the tick marks
         if(xstep == 0|| ystep == 0 || zstep == 0)
         {

	        finalm = new double[l2+30][3];
            for(int i =0; i<l2; i++)
            {
		     finalm[i] = final2[i];
	        }

            int f1 = 0;
            double h6 = (temp1[p.length-1] - temp1[0])/4;
		    for(int i=l2; i<l2+9; i=i+2)
		    {
		      finalm[i][0] = temp1[0] + f1*h6;
		      finalm[i][1] = yax+2*h4;
		      finalm[i][2] = zax;
		      finalm[i+1][0] = finalm[i][0];
		      finalm[i+1][1] = yax-2*h4;
		      finalm[i+1][2] = zax;
		      f1++;
		    }

            int f2 = 0;
            double h7 = (temp2[p.length-1] - temp2[0])/4;
	        for(int i=l2+10; i<l2+19; i=i+2)
		    {
		      finalm[i][0] = xax+2*h3;
		      finalm[i][1] = temp2[0] + f2*h7;
		      finalm[i][2] = zax;
		      finalm[i+1][0] = xax-2*h3;
		      finalm[i+1][1] = finalm[i][1];
		      finalm[i+1][2] = zax;
		      f2++;
		    }

	        int f3 = 0;
	        double h8 = (temp3[p.length-1] - temp3[0])/4;
		    for(int i=l2+20; i<l2+29; i=i+2)
		    {
		      finalm[i][0] = xax+2*h3;
		      finalm[i][1] = yax;
		      finalm[i][2] = temp3[0] + f3*h8;
		      finalm[i+1][0] = xax-2*h3;
		      finalm[i+1][1] = yax;
		      finalm[i+1][2] = finalm[i][2];
		      f3++;
           }
	   }

	   else
	   {
		   int k1 = (int)((temp1[p.length-1] -  xstart)/xstep);
		   int k2 = (int)((temp2[p.length-1] -  ystart)/ystep);
		   int k3 = (int)((temp3[p.length-1] -  zstart)/zstep);

		   finalm = new double[l2 + 2*(k1+k2+k3+3)][3];

           for(int i=0; i<l2; i++)
		   {
		    finalm[i] = final2[i];
	       }

	       int f1 = 0;
	       for(int i=l2; i<l2+2*k1+1; i=i+2)
		   {
		      finalm[i][0] = xstart + f1*xstep;
		      finalm[i][1] = yax+2*h4;
		      finalm[i][2] = zax;
		    finalm[i+1][0] = finalm[i][0];
		    finalm[i+1][1] = yax-2*h4;
		    finalm[i+1][2] = zax;
		    f1++;
	       }

           int f2 = 0;
           for(int i=l2+2*(k1+1); i<l2+2*(k1+k2)+3; i=i+2)
		   {
		      finalm[i][0] = xax+2*h3;
		      finalm[i][1] = ystart + f2*ystep;
		      finalm[i][2] = zax;
		    finalm[i+1][0] = xax-2*h3;
		    finalm[i+1][1] = finalm[i][1];
		    finalm[i+1][2] = zax;
		    f2++;
		   }

           int f3 = 0;
	       for(int i=l2+2*(k1+k2+2); i<l2+2*(k1+k2+k3)+5; i=i+2)
		   {
		      finalm[i][0] = xax+2*h3;
			  finalm[i][1] = yax;
			  finalm[i][2] = zstart + f3*zstep;
		    finalm[i+1][0] = xax-2*h3;
		    finalm[i+1][1] = yax;
		    finalm[i+1][2] = finalm[i][2];
		    f3++;
           }


	   }


  // Returning the final point array
    return finalm;
   }

}
