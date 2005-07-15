/** This program implements MergeSort algorithm and returns the
 *  original indices of the entries of the sorted array.
 */

class Sort_ind
{
  double[] numbers;

 /** Enter the array to be sorted here in the argument.*/
   public Sort_ind(double numb[])
    {
     numbers = numb;
    }


    int[] mergeSort()
    {
     int[] flag = new int[numbers.length];
     for(int i=0; i < numbers.length; i++)
     {
      flag[i] = i+1;
     }

     m_sort(numbers, flag, 0, numbers.length - 1);
     return flag;
    }

  /** The main algorithm is implemented in this part of the program*/
   void m_sort(double a[], int b[], int low, int high)
   {
      if(low == high)
           return;
      int length = high-low+1;
       int pivot = (low+high) / 2;
      m_sort(a, b, low, pivot);
      m_sort(a, b, pivot+1, high);

      double working[] = new double[length];
      int working_1[]  = new int [length];

      for(int i = 0; i < length; i++)
      {
       working[i]  = a[low+i];
       working_1[i]= b[low+i];
      }

      int m1 = 0;
      int m2 = pivot-low+1;

      for(int i = 0; i < length; i++)
      {
        if(m2 <= high-low)
        {
          if(m1 <= pivot-low)
           {
             if(working[m1] > working[m2])
             {
              a[i+low] = working[m2++];
              b[i+low] = working_1[m2-1];
             }
             else
             {
              a[i+low] = working[m1++];
              b[i+low] = working_1[m1-1];
             }
           }
          else
           {
            a[i+low] = working[m2++];
            b[i+low] = working_1[m2-1];
           }
        }
        else
        {
         a[i+low] = working[m1++];
         b[i+low] = working_1[m1-1];
        }
     }
   }
}

