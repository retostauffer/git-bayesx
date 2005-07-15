/** Given a set of 3D scattered data points
*   this Program is used for interpolation using
*   "Shepard's Method".
*   It takes the data set as a 2-d array
*   and the point where we need to evaluate the
*   function value,as input and returns function
*   value as output.
*/

 import java.math.*;
 import java.util.*;

 public class Interpolate
 {
	double[][] a;

/** The constructor for the class takes the available data
  * as the argument*/
	public Interpolate(double[][] d)
	{
	 a = d;
	}

/** This part of the program implements the Shepherd's algorithm
 *  @param The coordinates of the point where the function value is to be evaluated
 *  @return The function value at point poi
 */
    double zval(double[] poi)
	{
	  /**Evaluating distance of the interpolation point from each of the scatter points.*/
	  double[] h = new double [a.length];
	  for( int i=0; i<a.length; i++ )
	  {
	      double flag = 0;
	      for( int j=0; j<2; j++ )
	      {
	       flag = flag+Math.pow((a[i][j]-poi[j]),2);
	      }
	       h[i] = Math.sqrt(flag);
	      if(h[i]==0)
	      {
	       return a[i][2];
	      }
	  }

	  /** Finding the 5 nearest points to the point of interest.*/
	     int temp1[] = new int[a.length];
	     Sort_ind s1 = new Sort_ind(h);
	           temp1 = s1.mergeSort();

	  /** Finding the distance from the interpolation point to the most distant scatter point in the 5-point neighbourhood.*/
	         Sort s2 = new Sort(h);
	  double temp2[] = new double [a.length];
	           temp2 = s2.mergeSort();
	        double r = temp2[4];

	  /** Evaluating the weights for each of the 5neighbouring points.*/
       double flag_1 = 0;
       for ( int i=0; i<5; i++ )
       {
        flag_1 = flag_1 + Math.pow(((r-temp2[i])/(r*temp2[i])),2);
       }

	  double[] w = new double[5];
	  if(flag_1 <=0.0001)
	  {
	    for ( int i=0; i<5; i++ )
         {
          w[i] = 0.20;
         }
      }
      else
      {
	    for ( int i=0; i<5; i++ )
         {
          w[i] = Math.pow(((r-temp2[i])/(r*temp2[i])),2)/flag_1;
         }
      }

      /** Evaluating the function value at the interpolation point.*/
	  double val = 0;
	  for( int i=0; i<5; i++)
       {
        val = val+w[i]*a[temp1[i]-1][2];
       }

      return val;
   }
}
