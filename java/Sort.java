/** This program implements MergeSort algorithm and returns the
  * sorted array.*/

public class Sort
{
  double[] numbers;

 /** Enter the array to be sorted here in the argument.*/
   public Sort(double numb[])
    {
     numbers = numb;
    }


    double[] mergeSort()
    {
     m_sort(numbers, 0, numbers.length - 1);
     return numbers;
    }

  /** The main algorithm is implemented in this part of the program*/
   void m_sort(double a[], int low, int high)
   {
      if(low == high)
           return;
      int length = high-low+1;
       int pivot = (low+high) / 2;
      m_sort(a, low, pivot);
      m_sort(a, pivot+1, high);

      double working[] = new double[length];
      for(int i = 0; i < length; i++)
      {
       working[i]  = a[low+i];
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
             }
             else
             {
              a[i+low] = working[m1++];
             }
           }
          else
           {
            a[i+low] = working[m2++];
           }
        }
        else
        {
         a[i+low] = working[m1++];
        }
     }
   }
}


