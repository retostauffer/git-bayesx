// DATE: 20.01.98


#include "statmat.h"
#include <fstream.h>

//------------------------------------------------------------------------------
//----------- CLASS statmatrix: implementation of member functions -------------
//------------------------------------------------------------------------------

template<class T>
statmatrix<T>::statmatrix(const SparseMatrix & m)
                            : Matrix<T>(m.get_rows(),m.get_cols())
  {
  register unsigned i,j;
  double * work = getV();
  for(i=0;i<rows();i++)
    for(j=0;j<cols();j++,work++)
      *work = m(i,j);
  }


template<class T>
statmatrix<T>::statmatrix(const vector<T> & v)
                            : Matrix<T>(v.size(),1)
  {
  register unsigned i;
  double * work = getV();
  for(i=0;i<rows();i++,work++)
    *work = v[i];
  }


template<class T>
const statmatrix<T> & statmatrix<T>::operator=(const SparseMatrix & m)
  {
  statmatrix<T> res(m.get_rows(),m.get_cols());
  double * work = res.getV();
  register unsigned i,j;
  for(i=0;i<res.rows();i++)
    for(j=0;j<res.cols();j++,work++)
      *work = T(m(i,j));
  return res;
  }


// OVERLOADED ASSIGNMENT OPERATOR
	template<class T>
   statmatrix<T>  statmatrix<T>::operator= (const adja& a)
	{
		datamatrix res(a.rows(),a.cols());

		double * work = res.getV();
		register unsigned i,j;
		for(i=0;i<res.rows();i++)
			for(j=0;j<res.cols();j++,work++)
				 *work = a(i,j);
		return res;
	}


template<class T>
void statmatrix<T>::assign(const statmatrix & A)
  {
  assert(rows()==A.rows());
  unsigned size = rows( ) * cols( );

  T *workA = A.getV( );
  T *workR = getV( );
  register unsigned i;
  for ( i = 0;i < size;i++, workA++,workR++ )
    *workR = *workA;
  }


template<class T>
void statmatrix<T>::plus(const statmatrix & A,const statmatrix & B)
  {

  unsigned size = rows( ) * cols( );

  T *workA = A.getV( );
  T *workB = B.getV( );
  T *workR = getV( );
  register unsigned i;

  for ( i = 0;i < size;i++, workA++, workB++, workR++ )
        *workR = *workA + *workB;
  }

template<class T>
void statmatrix<T>::plus(const statmatrix & A)
  {

  unsigned size = rows( ) * cols( );

  T *workA = A.getV( );
  T *workR = getV( );
  register unsigned i;

  for ( i = 0;i < size;i++, workA++, workR++ )
        *workR += *workA;
  }


template<class T>
void statmatrix<T>::minus(const statmatrix & A,const statmatrix & B)
  {
  unsigned size = rows( ) * cols( );

  T *workA = A.getV( );
  T *workB = B.getV( );
  T *workR = getV( );
  register unsigned i;

  for ( i = 0;i < size;i++, workA++, workB++, workR++ )
    *workR = *workA - *workB;
  }


template<class T>
void statmatrix<T>::minus(const statmatrix & A,const statmatrix & B,
                         const unsigned & colA,const unsigned & colB)
  {

  unsigned size = rows();
  unsigned sizeA = A.cols();
  unsigned sizeB = B.cols();
  register unsigned i;
  T * workA = A.getV()+colA;
  T * workB = B.getV()+colB;
  T* workR = getV();
  for (i=0;i<size;i++,workA+=sizeA,workB+=sizeB,workR++)
    *workR = *workA- *workB;

  }


template<class T>
void statmatrix<T>::mult(const statmatrix & A,const statmatrix & B)
  {

  assert(A.cols() == B.rows());
  assert(rows() == A.rows());
  assert(cols() == B.cols());

  T * workA;
  T * workB;
  T * workR = getV();
  unsigned n = cols();
  unsigned size = rows()*n;
  register unsigned i, k;

  for (i=0;i <size;i++,workR++)
    {

    *workR=T(0);
    workA = A.getV( ) +  (i / n) * A.cols();
    workB = B.getV( ) +  i % n;

    for (k = 0; k < A.cols(); ++k, ++workA, workB += n )
      if (!(*workA == T(0) || *workB == T(0)))
        *workR += *workA * *workB;

    }
  }



template<class T>
void statmatrix<T>::addmult(const statmatrix & A,const statmatrix & B)
  {

  assert(A.cols() == B.rows());
  assert(rows() == A.rows());
  assert(cols() == B.cols());

  T * workA;
  T * workB;
  T * workR = getV();
  unsigned n = cols();
  unsigned size = rows()*n;
  register unsigned i, k;

  for (i=0;i <size;i++,workR++)
    {

    workA = A.getV( ) +  (i / n) * A.cols();
    workB = B.getV( ) +  i % n;

    for (k = 0; k < A.cols(); ++k, ++workA, workB += n )
      if (!(*workA == T(0) || *workB == T(0)))
        *workR += *workA * *workB;

    }
  }


template<class T>
void statmatrix<T>::addmultsym(const statmatrix & A,const statmatrix & B)
  {

  assert(A.cols() == B.rows());
  assert(rows() == A.rows());
  assert(cols() == B.cols());

  T * workA;
  T * workB;
  T * workR = getV();
  unsigned n = cols();
  unsigned m = A.cols();

  register unsigned i,j,k;

  for (i=0;i<rows();i++)
    for(j=0;j<n;j++,workR++)
      {

      workA = A.getV( ) +  i*m;
      workB = B.getV( ) +  j;

      for(k=0;k<=i;k++,workB+=n,workA++)
        *workR += * workA * *workB;
//        *workR += A(i,k)*B(j,k);

      workA = A.getV() + (i+1)*m + i;

      for(k=i+1;k<rows();k++,workB+=n,workA+=m)
        *workR += *workA * *workB;
//        *workR += A(k,i)*B(j,k);

      }

  }

template<class T>
statmatrix<T> statmatrix<T>::inverse(void)
  {
  assert(rows()==cols());
  if (rows() == 1)
    {
    assert( *getV() != T(0) );
    return statmatrix<T>(1,1,T(1)/(*getV()));
    }
  else if (rows()==2)
    {
    T det = get(0,0)*get(1,1)-get(0,1)*get(1,0);
    assert(det !=  T(0));
    statmatrix<T> result(2,2);
    T* work = result.getV();
    *work =  get(1,1)/det;                 // result(0,0)
    work++;
    *work =  -get(0,1)/det;                // result(0,1)
    work++;
    *work =  -get(1,0)/det;                // result(1,0)
    work++;
    *work =  get(0,0)/det;                 // result(0,0)
    return result;
    }
  else
    return Matrix<T>::inverse();
  }

template<class T>
void statmatrix<T>::multdiagback(const statmatrix & d)
  {
  T* dpoint;
  T* thispoint=getV();
  unsigned i, j;
  for(i=0; i<rows(); i++)
    {
    for(j=0, dpoint=d.getV(); j<cols(); j++, dpoint++, thispoint++)
      {
      *thispoint *= *dpoint;
      }
    }
  }

template<class T>
void statmatrix<T>::multdiagfront(const statmatrix & d)
  {
  T* dpoint=d.getV();
  T* thispoint=getV();
  unsigned i, j;
  for(i=0; i<rows(); i++, dpoint++)
    {
    for(j=0; j<cols(); j++, thispoint++)
      {
      *thispoint *= *dpoint;
      }
    }
  }

template<class T>
void statmatrix<T>::multdiagfront(const statmatrix & A, const statmatrix & d)
  {
  assert(A.rows()==rows());
  assert(A.cols()==cols());
  assert(A.rows()==d.rows());
  T* dpoint=d.getV();
  T* apoint=A.getV();
  T* thispoint=getV();
  unsigned i, j;
  for(i=0; i<rows(); i++, dpoint++)
    {
    for(j=0; j<cols(); j++, thispoint++, apoint++)
      {
      *thispoint = *apoint * *dpoint;
      }
    }
  }

template<class T>
void statmatrix<T>::addtodiag(const statmatrix & d, unsigned first,
                              unsigned last)
  {
  assert(rows()==cols());
  T* dpoint=d.getV();
  T* thispoint=getV()+first*cols()+first;
  unsigned i;
  for(i=first; i<last; i++, dpoint++, thispoint+=cols()+1)
    {
    *thispoint += *dpoint;
    }
  }

template<class T>
void statmatrix<T>::subfromdiag(const statmatrix & d, unsigned first,
                              unsigned last)
  {
  assert(rows()==cols());
  T* dpoint=d.getV();
  T* thispoint=getV()+first*cols()+first;
  unsigned i;
  for(i=first; i<last; i++, dpoint++, thispoint+=cols()+1)
    {
    *thispoint -= *dpoint;
    }
  }

template<class T>
void statmatrix<T>::elemmult(const statmatrix<T> & A)
  {
  assert(A.cols()==cols());
  assert(A.rows()==rows());
  T* Apoint = A.getV();
  T* thispoint = getV();
  unsigned i;
  for(i=0; i<cols()*rows(); i++, thispoint++, Apoint++)
    {
    *thispoint *= *Apoint;
    }
  }

template<class T>
void statmatrix<T>::elemquot(const statmatrix<T> & A)
  {
  assert(A.cols()==cols());
  assert(A.rows()==rows());
  T* Apoint = A.getV();
  T* thispoint = getV();
  unsigned i;
  for(i=0; i<cols()*rows(); i++, thispoint++, Apoint++)
    {
    *thispoint /= *Apoint;
    }
  }

template<class T>
void statmatrix<T>::weightedsscp(const statmatrix<T> & X,
                                 const statmatrix<T> & w)
  {
  unsigned xcols=X.cols();
  unsigned n=X.rows();

  assert(cols()==xcols);
  assert(rows()==xcols);
  assert(w.rows()==n);

  T* xpointi;
  T* xpointj;
  T* wpoint;

  register unsigned i,j,k;
  double sum;

  for(i=0; i<xcols; i++)
    {
    for(j=i; j<xcols; j++)
      {
      sum=0;
      xpointi=X.getV()+i;
      xpointj=X.getV()+j;
      wpoint=w.getV();

//      for(k=0; k<n; k++)
      for(k=0; k<n; k++, wpoint++, xpointi+=xcols, xpointj+=xcols)
        {
//        sum += X(k,i)*X(k,j)*w(k,0);
        if(!(*xpointi==T(0)||*xpointj==T(0)))
          {
          sum += *xpointi * *xpointj * *wpoint;
          }
        }
      put(i,j,sum);
      if(i!=j)
        {
        put(j,i,sum);
        }
      }
    }
  }

template<class T>
void statmatrix<T>::weightedsscp2(const statmatrix<T> & X, const statmatrix<T> & Z,
                                  const statmatrix<T> & w)
  {
  unsigned xcols=X.cols();
  unsigned zcols=Z.cols();
  unsigned n=Z.rows();

  assert(cols()==xcols+zcols);
  assert(rows()==xcols+zcols);
  assert(w.rows()==n);
  assert(X.rows()==n);

  T* xpointi;
  T* xpointj;
  T* zpointi;
  T* zpointj;
  T* wpoint;

  register unsigned i,i1,j,j1,k;
  double sum;

// compute X'WX, X'WZ and Z'WX
  for(i=0; i<xcols; i++)
    {
    for(j=i; j<xcols; j++)
      {
      sum=0;
      xpointi=X.getV()+i;
      xpointj=X.getV()+j;
      wpoint=w.getV();

//      for(k=0; k<n; k++)
      for(k=0; k<n; k++, wpoint++, xpointi+=xcols, xpointj+=xcols)
        {
//        sum += X(k,i)*X(k,j)*w(k,0);
        if(!(*xpointi==T(0)||*xpointj==T(0)))
          {
          sum += *xpointi * *xpointj * *wpoint;
          }
        }
      put(i,j,sum);
      if(i!=j)
        {
        put(j,i,sum);
        }
      }
    for(j=0, j1=xcols; j<zcols; j++, j1++)
      {
      sum=0;
      xpointi=X.getV()+i;
      zpointj=Z.getV()+j;
      wpoint=w.getV();

//      for(k=0; k<n; k++)
      for(k=0; k<n; k++, wpoint++, xpointi+=xcols, zpointj+=zcols)
        {
//        sum += X(k,i)*Z(k,j)*w(k,0);
        if(!(*xpointi==T(0)||*zpointj==T(0)))
          {
          sum += *xpointi * *zpointj * *wpoint;
          }
        }
      put(i,j1,sum);
      put(j1,i,sum);
      }
    }

// compute Z'WZ
  for(i=0, i1=xcols; i<zcols; i++, i1++)
    {
    for(j=i, j1=i+xcols; j<zcols; j++, j1++)
      {
      sum=0;
      zpointi=Z.getV()+i;
      zpointj=Z.getV()+j;
      wpoint=w.getV();

//      for(k=0; k<n; k++)
      for(k=0; k<n; k++, wpoint++, zpointi+=zcols, zpointj+=zcols)
        {
//        sum += Z(k,i)*Z(k,j)*w(k,0);
        if(!(*zpointi==T(0)||*zpointj==T(0)))
          {
          sum += *zpointi * *zpointj * *wpoint;
          }
        }
      put(i1,j1,sum);
      if(i!=j)
        {
        put(j1,i1,sum);
        }
      }
    }
  }

template<class T>
void statmatrix<T>::weightedsscp_resp(const statmatrix & X, const statmatrix & y,
                          const statmatrix & w)
  {
  unsigned xcols=X.cols();
  unsigned n=X.rows();

  assert(rows()==xcols);
  assert(w.rows()==n);
  assert(y.rows()==n);

  register unsigned i,k;
  double sum;

  T* wpoint=w.getV();
  T* ypoint=y.getV();

  statmatrix<T>wy(n,1);
  T* wypoint=wy.getV();

  for(i=0; i<n; i++, ++wpoint, ++ypoint, ++wypoint)
    {
    *wypoint = *ypoint * *wpoint;
    }

  T* xpointi;
  T* thispoint=getV();

  for(i=0; i<xcols; i++, ++thispoint)
    {
    sum=0;
    wypoint=wy.getV();
    xpointi=X.getV()+i;
    for(k=0; k<n; k++, xpointi+=xcols, ++wypoint)
      {
      if(*xpointi!=T(0))
        {
        sum += *xpointi * *wypoint;
        }
      }
    *thispoint = sum;
    }
  }

template<class T>
void statmatrix<T>::weightedsscp_resp2(const statmatrix<T> & X,
                                       const statmatrix<T> & Z,
                                       const statmatrix<T> & y,
                                       const statmatrix<T> & w)
  {
  unsigned xcols=X.cols();
  unsigned zcols=Z.cols();
  unsigned n=Z.rows();

  assert(rows()==xcols+zcols);
  assert(w.rows()==n);
  assert(y.rows()==n);
  assert(X.rows()==n);

  register unsigned i,k;
  double sum;

  T* wpoint=w.getV();
  T* ypoint=y.getV();

  statmatrix<T>wy(n,1);
  T* wypoint=wy.getV();

  for(i=0; i<n; i++, ++wpoint, ++ypoint, ++wypoint)
    {
    *wypoint = *ypoint * *wpoint;
    }

  T* xpointi;
  T* zpointi;
  T* thispoint=getV();

// compute X'Wy
  for(i=0; i<xcols; i++, ++thispoint)
    {
    sum=0;
    wypoint=wy.getV();
    xpointi=X.getV()+i;
    for(k=0; k<n; k++, xpointi+=xcols, ++wypoint)
      {
      if(*xpointi!=T(0))
        {
        sum += *xpointi * *wypoint;
        }
      }
    *thispoint = sum;
    }

// compute Z'Wy
  for(i=0; i<zcols; i++, ++thispoint)
    {
    sum=0;
    zpointi=Z.getV()+i;
    wypoint=wy.getV();

    for(k=0; k<n; k++, wypoint++, zpointi+=zcols)
      {
      if(*zpointi!=T(0))
        {
        sum += *zpointi * *wypoint;
        }
      }
    *thispoint = sum;
    }
  }

template<class T>
void statmatrix<T>::sort(int start,int ende,int col)
  {
  int i = start;
  int j = ende;
  T x = get((start+ende)/2,col);
  statmatrix<T> hilfe;
  do
	 {
	 while (get(i,col) < x)
		i++;
	 while (x < get(j,col))
		j--;
	 if (i <= j)
		{
		hilfe = getRow(i);
		putRow(i,getRow(j));
		putRow(j,hilfe);
		i++;
		j--;
		}
	 }
  while ( i <= j );
	 if (start < j)
		sort(start,j,col);
	 if (i < ende)
		sort(i,ende,col);
  }


template<class T>
void statmatrix<T>::sortcol(int start,int ende,int col)
  {
  int i = start;
  int j = ende;
  T x = get((start+ende)/2,col);
  T hilfe;
  do
	 {
	 while (get(i,col) < x)
		i++;
	 while (x < get(j,col))
		j--;
	 if (i <= j)
		{
		hilfe = get(i,col);
		put(i,col,get(j,col));
		put(j,col,hilfe);
		i++;
		j--;
		}
	 }
  while ( i <= j );
	 if (start < j)
		sortcol(start,j,col);
	 if (i < ende)
		sortcol(i,ende,col);
  }


template<class T>
void statmatrix<T>::indexinit(void)
  {
  unsigned i;
  unsigned j;
  for (i=0;i<cols();i++)
	 for (j=0;j<rows();j++)
		put(j,i,j);
  }


template<class T>
void statmatrix<T>::indexsort(statmatrix<int> & index,int start,int ende,
										int col,int indexcol) const
  {
  int i = start;
  int j = ende;
  T x = get(index((start+ende)/2,indexcol),col);
  int hilfe;
  do
	 {
	 while (get(index(i,indexcol),col) < x)
		i++;
	 while (x < get(index(j,indexcol),col))
		j--;
	 if (i <= j)
		{
		hilfe = index(i,indexcol);
		index(i,indexcol) = index(j,indexcol);
		index(j,indexcol) = hilfe;
		i++;
		j--;
		}
	 }
  while ( i <= j );
	 if (start < j)
		indexsort(index,start,j,col,indexcol);
	 if (i < ende)
		indexsort(index,i,ende,col,indexcol);
  }


template<class T>
statmatrix<T> statmatrix<T>::sum (void) const
  {
  statmatrix<T> s(cols(),1,0);
  unsigned col;
  for (col=0;col<cols();col++)
    s(col,0) = sum(col);
  return s;
  }


template<class T>
T statmatrix<T>::sum (const unsigned & col) const
  {

  assert(col < cols());

  T sum = 0;
  register unsigned i;
  T* work = getV()+col;
  for (i=0;i<rows();i++,work+=cols())
    sum += *work;
  return sum;
  }


template<class T>
T statmatrix<T>::sum2 (const unsigned & col) const
  {

  assert(col < cols());

  T sum = 0;
  register unsigned i;
  T* work = getV()+col;
  for (i=0;i<rows();i++,work+=cols())
    sum += *work * *work;
  return sum;
  }


template<class T>
T  statmatrix<T>::sum2(const unsigned & col,const statmatrix<T> & weight) const
  {

  assert(col < cols());

  T sum = 0;
  T* work = getV()+col;
  T* workweight = weight.getV();
  register unsigned i;
  for (i=0;i<rows();i++,work+=cols(),workweight++)
    {
    sum += *workweight* (*work) * (*work);
    }

  return sum;

  }

template<class T>
statmatrix<T> statmatrix<T>::sum2()
  {
  statmatrix<T>res(cols(),1,0);
  unsigned i;
  for(i=0; i<cols(); i++)
    {
    res(i,0)=sum2(i);
    }
  return res;
  }

template<class T>
T statmatrix<T>::mean(const unsigned & col,
                      const statmatrix<T> & weight) const
  {
  assert(col < cols());

  T sum = 0;
  T sumweight = 0;
  register unsigned i;
  T* work = getV()+col;
  T* workweight = weight.getV();
  for (i=0;i<rows();i++,work+=cols(),workweight++)
    {
    sumweight+= *workweight;
    sum += *workweight * *work;
    }
  return sum/sumweight;
  }


template<class T>
T statmatrix<T>::min (const unsigned & c) const
  {
  T* work = getV()+c;
  T minv = *work;
  work+=cols();
  unsigned i;
  for (i=1;i<rows();i++,work+=cols())
    {
    if ((*work) < minv)
      minv = *work;
    }
  return minv;
  }


template<class T>
T statmatrix<T>::max (const unsigned & c) const
  {
  T* work = getV()+c;
  T maxv = *work;
  work+=cols();
  unsigned i;
  for (i=1;i<rows();i++,work+=cols())
    {
    if ((*work) > maxv)
      maxv = *work;
    }
  return maxv;
  }


template<class T>
T statmatrix<T>::sumcomplete(void) const
  {
  register unsigned i;
  unsigned size = rows()*cols();
  T* work = getV();
  T sum = 0;
  for (i=0;i<size;i++,work++)
	 sum += *work;
  return sum;
  }

template<class T>
T statmatrix<T>::norm(unsigned col)
  {
  T norm=0;
  norm = sqrt(sum2(col));
  return norm;
  }

template<class T>
statmatrix<T> statmatrix<T>::norm()
  {
  statmatrix<T>res(cols(),1,0);
  unsigned i;
  for(i=0; i<cols(); i++)
    {
    res(i,0)=norm(i);
    }
  return res;
  }

template<class T>
statmatrix<T> statmatrix<T>::mean() const
  {
  statmatrix<T> m(cols(),1);
  for (unsigned col=0;col<cols();col++)
	 m(col,0) = mean(col);
  return m;
  }

template<class T>
T statmatrix<T>::var(const unsigned & col) const
  {
  T m = mean(col);
  return T(1)/T(rows())*sum2(col)-m*m;
  }

template<class T>
T statmatrix<T>::var(const unsigned & col,
                     const statmatrix<double> & weight) const
  {
  T m = mean(col,weight);
  T ws = weight.sum(0);
  T s2 = sum2(col,weight);
  return T(1)/ws*s2-m*m;
  }

template<class T>
T statmatrix<T>::quantile(const T & percent,const unsigned & col) const
  {

  T k = rows()*(percent/100.0);           // (alpha * Anzahl der Beobachtungen)
  unsigned kganz = unsigned(k);

  statmatrix<int> index(rows(),1);
  index.indexinit();
  indexsort(index,0,rows()-1,col,0);

  if (k == kganz)                             // Falls k ganzzahlig ist
	 return (get(index(kganz-1,0),col) + get(index(kganz,0),col))/2.0;
  else                                       // Falls k nicht ganzzahlig ist
	 return get(index(kganz,0),col);

  }


template<class T>
statmatrix<T> statmatrix<T>::quantile(T percent)
  {
  statmatrix<T> quant(cols(),1);
  for (int col=0;col<cols();col++)
	 quant(col,0) = quantile(percent,col);
  return quant;
  }


template<class T>
T statmatrix<T>::autocorr (const unsigned & lag,const unsigned & col) const
  {

  T sum = 0;                               // Summe der Werte
  T sum2 = 0;                              // Quadratsumme der Werte
  T sum_lag = 0;                           // Summe der verzögerten Werte
  T sum_lag2 = 0;                          // Quadratsumme der verzögerten Werte
  T sum_wertlag = 0;                       // Summe Wert * verzögerter Wert
  T mean,mean_lag;                         // Mittelwert, verzögerter Mittelwert
  T anz = rows()-lag;                      // Anzahl Beobachtungen

  for (unsigned k=lag;k<rows();k++)
	 {
	 sum = sum + get(k,col);
	 sum2 = sum2 + get(k,col)*get(k,col);
	 sum_lag = sum_lag + get(k-lag,col);
	 sum_lag2 = sum_lag2 + get(k-lag,col)*get(k-lag,col);
	 sum_wertlag = sum_wertlag + get(k,col)*get(k-lag,col);
	 }

  mean = (1.0/anz)*sum;
  mean_lag = (1.0/anz)*sum_lag;

  return	(sum_wertlag - anz*mean*mean_lag)/
			 sqrt( (sum2-anz*mean*mean)*(sum_lag2-anz*mean_lag*mean_lag) );

  }

//------------------------------------------------------------------------------

template<class T>
statmatrix<T> statmatrix<T>::autocorr (const unsigned & lag) const
  {

  statmatrix<T> corr(lag,cols());

  for (unsigned i=1;i<=lag;i++)
	 for (unsigned j=0;j<cols();j++)
		corr(i-1,j) = autocorr(i,j);

  return corr;
  }

//------------------------------------------------------------------------------


template<class T>
statmatrix<T> statmatrix<T>::autocorr(const unsigned & beginlag,
                                      const unsigned & endlag,
                                      const unsigned & col) const
  {



  unsigned rowstot = endlag-beginlag+1;
  unsigned i;
  statmatrix corr(rowstot,1);


  T sum = 0;                               // Summe der Werte
  T sum2 = 0;                              // Quadratsumme der Werte
  T sum_lag = 0;                           // Summe der verzögerten Werte
  T sum_lag2 = 0;                          // Quadratsumme der verzögerten Werte
  T sum_wertlag = 0;                       // Summe Wert * verzögerter Wert
  T mean,mean_lag;                         // Mittelwert, verzögerter Mittelwert
  T anz = rows() - beginlag;                   // Anzahl Beobachtungen

  for (unsigned k=beginlag;k<rows();k++)
    {
	 sum = sum + get(k,col);
	 sum2 = sum2 + get(k,col)*get(k,col);
	 sum_lag = sum_lag + get(k-beginlag,col);
	 sum_lag2 = sum_lag2 + get(k-beginlag,col)*get(k-beginlag,col);
	 sum_wertlag = sum_wertlag + get(k,col)*get(k-beginlag,col);
	 }

  mean = (1.0/anz)*sum;
  mean_lag = (1.0/anz)*sum_lag;

  if ((sum2-anz*mean*mean)*(sum_lag2-anz*mean_lag*mean_lag) <= 0)
    corr(0,0) = 2;
  else
    corr(0,0) = (sum_wertlag - anz*mean*mean_lag)/
			     sqrt( (sum2-anz*mean*mean)*(sum_lag2-anz*mean_lag*mean_lag) );


  for(i=beginlag+1;i<=endlag;i++)
    {

    sum -= get(i-1,col);
    sum2 -= get(i-1,col)*get(i-1,col);
    sum_lag -= get(rows()-1-(i-1),col);
    sum_lag2 -=  get(rows()-1-(i-1),col)*get(rows()-1-(i-1),col);

    sum_wertlag = 0;
    for (unsigned k=i;k<rows();k++)
      sum_wertlag += get(k,col)*get(k-i,col);

    anz--;
    mean = (1.0/anz)*sum;
    mean_lag = (1.0/anz)*sum_lag;


    if ((sum2-anz*mean*mean)*(sum_lag2-anz*mean_lag*mean_lag) <= 0)
      corr(i-beginlag,0) = 2;
    else
      corr(i-beginlag,0) =  (sum_wertlag - anz*mean*mean_lag)/
			     sqrt( (sum2-anz*mean*mean)*(sum_lag2-anz*mean_lag*mean_lag) );

//    corr(i-beginlag,0) = autocorr(i,0);

    }


  return corr;

  }


template<class T>
statmatrix<T> statmatrix<T>::cov()
  {
  statmatrix<T> one(rows(),1,1);

  return (1.0/(rows()-1))*( (*this).transposed()*(*this) -
			(1.0/rows())*(*this).transposed()*one*one.transposed()*(*this) );
  }


template<class T>
statmatrix<T> statmatrix<T>::corr()
  {
  int i,j;
  statmatrix<T> c = cov();
  statmatrix<T> co(cols(),cols());
  for (i=0;i<c.rows();i++)
	 for(j=0;j<c.cols();j++)
		co(i,j) = c(i,j)/sqrt(c(i,i)*c(j,j));
  return co;
  }


template<class T>
T statmatrix<T>::compute_quadform(const statmatrix<T> & x,const unsigned & c)
  {
  unsigned i,j;
  T res=0;
  T * xi=x.getV()+c;
  T * xj;
  T * workm=getV();
  unsigned d = x.cols();
  for (i=0;i<rows();i++,xi+=d)
    {
    workm+=i;
//    res+= x(i,0)*x(i,0)*get(i,i);
    res+= *xi * *xi * *workm;
    xj = xi+d;
    workm++;
    for(j=i+1;j<cols();j++,xj+=d,workm++)
      {
//      res+=2*x(i,0)*x(j,0)*get(i,j);
      res+=2* *xi * *xj * *workm;

      }

    }

  return res;

  }


template<class T>
statmatrix<T> statmatrix<T>::strike (unsigned int k)
{
	unsigned int i,j;
	unsigned rows_new = rows()-1;

	statmatrix<T> matrix_new (rows_new,rows_new);

	if(k==0)
	{
		for(i=0; i<rows_new; i++)
			for(j=0; j<rows_new; j++)
				matrix_new(i,j)=get(i+1,j+1);
	}
	else if(k==rows_new+1)
	{
		for(i=0; i<rows_new; i++)
			for(j=0; j<rows_new; j++)
				matrix_new(i,j)=get(i,j);
	}

	else
	{
		for(i=0; i<rows_new; i++)
		{
			for(j=0; j<rows_new; j++)
			{
				if(i<k && j<k)
					 matrix_new(i,j)=get(i,j);
				else if(i<k && j>k-1)
					matrix_new(i,j)=get(i,j+1);
				else if(i>k-1 && j<k)
					matrix_new(i,j)=get(i+1,j);
				else if (i>k-1 && j>k-1)
					matrix_new(i,j)=get(i+1,j+1); 
			}
		}
	}

	return matrix_new;				
}

template<class T>
statmatrix<T> statmatrix<T>::get_cov_iX (int i, int j)
{
	assert(rows()==cols());
	int k,l;
	datamatrix res (1,rows()-2);

	l=0;
	for(k=0; k<rows(); k++)
	{
		if(k==i)
			l--;
		else if(k==j)
			l--;
		else
			res(0,l) = get(i,k);
		l++;
	}
	
	


/*	if(i<j)
	{
		for(k=0; k<rows()-2; k++)
		{
			if(k<i)
				res(0,k) = get(i,k);

			else if(k>i-1 && k<j && i!=j-1)
				res(0,k) = get(i,k+1);

			else if( k>i-1  && k<j && i==j-1)
				res(0,k) = get(i,k+2);

			else if(k>j-1)
				res(0,k) = get(i,k+2);

			else 
				res(0,k) = get(i,k+2);
		}
	}
	else if (j<i)
	{
		for(k=0; k<rows()-2; k++)
		{
			if(k<j)
				res(0,k) = get(i,k);

			else if(k>j-1 && k<i && j!=i-1)
				res(0,k) = get(i,k+1);

			else if(k>j-1 && k<i && j==i-1)
				res(0,k) = get(i,k+2);

			else if(k>i-1)
				res(0,k) = get(i,k+2);

			else 
				res(0,k) = get(i,k+2);
		}
	} */
	

	return res;
}



template<class T>
statmatrix<T> statmatrix<T>::partial_var(void) 
{
	unsigned i,j,m,n;
	unsigned nvar = cols();
	unsigned nobs = rows(); 

	double numerator, denominator;

	datamatrix cov_all (nvar,nvar);
	datamatrix par_var (nvar,nvar,-999);

	cov_all.assign(cov());

	datamatrix var_iX (nvar-1,nvar-1);
	datamatrix var_jX (nvar-1,nvar-1);
	datamatrix cov_iX (1,nvar-2);
	datamatrix cov_jX (1,nvar-2);
	datamatrix var_X (nvar-2, nvar-2);

	datamatrix help1 (nvar-2,1);
	datamatrix help2 (1, 1);

	datamatrix test1 (1,nvar-2);
	datamatrix test2 (1,1);


	double var_i_X;
	double var_j_X;

	for(i=0; i<nvar; i++)
	{
		for(j=0; j<nvar; j++)
		{
			if(i<j)
			{
				cov_iX.assign(cov_all.get_cov_iX(i,j)); 
				cov_jX.assign(cov_all.get_cov_iX(j,i)); 

			//	cout<<cov_iX<<endl;
			//	cout<<cov_jX<<endl;

				datamatrix help (nvar-1, nvar-1);

				var_X.assign((cov_all.strike(i)).strike(j-1));

			//	cout<<var_X<<endl;

				test1.mult(cov_iX,var_X.inverse());
				test2.mult(test1,cov_iX.transposed());

				var_i_X = cov_all(i,i) - test2(0,0);

				test1.mult(cov_jX,var_X.inverse());
				test2.mult(test1,cov_jX.transposed());

				var_j_X = cov_all(j,j) - test2(0,0);

				help1.mult(var_X.inverse(), cov_jX.transposed());
				help2.mult(cov_iX,help1);

				numerator = cov_all(i,j) - help2(0,0);

				denominator = sqrt(var_i_X *var_j_X );

				par_var(i,j)=numerator/denominator;
				par_var(j,i)=numerator/denominator;
			}

			else if (i==j)
				par_var(i,j)=1;
		}
	}

	return par_var;
}

/*
template<class T>
statmatrix<T> statmatrix<T>::diag_one(void)
{
	assert(rows()==cols());
	unsigned i,j;
	statmatrix<T> matrix_new (rows(),cols(),0);

	for(i=0; i<rows();i++)
		for(j=0; j<rows();j++)
			matrix_new(i,j) = get(i,j)/ get(i,i);

	return matrix_new;
}


*/

statmatrix<double> multdiagback(datamatrix X, const datamatrix & d)
  {
  double* dpoint;
  double* Xpoint=X.getV();
  unsigned i, j;
  for(i=0; i<X.rows(); i++)
    {
    for(j=0, dpoint=d.getV(); j<X.cols(); j++, dpoint++, Xpoint++)
      {
      *Xpoint *= *dpoint;
      }
    }
  return X;
  }

statmatrix<double> multdiagfront(datamatrix X, const datamatrix & d)
  {
  double* dpoint=d.getV();
  double* Xpoint=X.getV();
  unsigned i, j;
  for(i=0; i<X.rows(); i++, dpoint++)
    {
    for(j=0; j<X.cols(); j++, Xpoint++)
      {
      *Xpoint *= *dpoint;
      }
    }
  return X;
  }

statmatrix<double> diffmat(const int k, const int d)
  {
  assert(k>=1);
  assert(k<=2);
  int i;
  statmatrix<double>res(d-k,d,0);
  for(i=0; i<d-k; i++)
    {
    if(k==1)
      {
      res(i,i)=-1;
      res(i,i+1)=1;
      }
    else if(k==2)
      {
      res(i,i)=1;
      res(i,i+1)=-2;
      res(i,i+2)=1;
      }
    }
  return res;
  }


statmatrix<double> diffmat_k(const int k, const int d)
  {
  assert(k>=0);
  assert(k<d);
  int i,j,m;
  statmatrix<double>res(d-k,d,0);
  for(i=0;i<res.rows();i++)
    res(i,i+k) = 1.0;
  for(i=1;i<=k;i++)
    {
    for(j=0;j<res.rows();j++)
      for(m=k-i;m<k;m++)
        res(j,j+m) = res(j,j+m)-res(j,j+m+1);
    }
  return res;
  }


statmatrix<double> weighteddiffmat(const int k, const vector<double> & weight)
  {
  assert(k>=1);
  assert(k<=2);
  unsigned i;
  unsigned d = weight.size();
  statmatrix<double>res(d-k,d,0);
  if(k==1)
    {
    for(i=0; i<d-1; i++)
      {
      res(i,i)=1/sqrt(weight[i+1]);
      res(i,i+1)=-res(i,i);
      }
    }
  else
    {
    for(i=0; i<d-2; i++)
      {
      res(i,i+2)=1/sqrt(weight[i+2]*(1+weight[i+2]/weight[i+1]));
      res(i,i+1)=-(1+weight[i+2]/weight[i+1])*res(i,i+2);
      res(i,i)=weight[i+2]/weight[i+1]*res(i,i+2);
      }
    }
  return res;
  }

statmatrix<double> seasonalfactor(const unsigned & per, const unsigned & s)
  {
  unsigned i,j;
  statmatrix<double> res(s,s-per+1,0);
  for(i=0; i<s-per+1; i++)
    {
    for(j=0; j<per; j++)
      {
      res(i+j,i)=1;
      }
    }
  return res;
  }

statmatrix<double> seasonalX(const unsigned & per, const unsigned & s)
  {
  unsigned i,k;
  statmatrix<double> res(s,per-1,0);
  for(i=0; i<s; i++)
    {
    k=(i+1)%per;
    if(k>0)
      {
      res(i,k-1)=1;
      }
    else
      {
      for(k=0; k<per-1; k++)
        {
        res(i,k)=-1;
        }
      }
    }
  return res;
  }

void rotate(statmatrix<double> & a, const double & s, const double & tau,
            const int & i, const int & j, const int & k, const int & l)
  {
  double g,h;

  g=a(i,j);
  h=a(k,l);
  a(i,j)=g-s*(h+g*tau);
  a(k,l)=h+s*(g-h*tau);
  }

void tridiag(statmatrix<double> & a, statmatrix<double> & d,
             statmatrix<double> & e)
  {
  assert(a.rows()==a.cols());
  assert(a.rows()==d.rows());
  assert(d.cols()==1);
  assert(e.rows()==a.rows());
  assert(e.cols()==1);

  int l, k, j, i;
  double scale, hh, h, g, f;

  int n = d.rows();
  for(i=n-1; i>0; i--)
    {
    l=i-1;
    h=scale=0;
    if(l>0)
      {
      for(k=0; k<l+1; k++)
        {
        scale += fabs(a(i,k));
        }
      if(scale==0)
        {
        e(i,0) = a(i,l);
        }
      else
        {
        for(k=0; k<l+1; k++)
          {
          a(i,k) /= scale;
          h += a(i,k)*a(i,k);
          }
        f = a(i,l);
        g = (f >= 0 ? -sqrt(h) : sqrt(h));
        e(i,0) = scale*g;
        h -= f*g;
        a(i,l) = f-g;
        f=0;
        for(j=0; j<l+1; j++)
          {
          a(j,i) = a(i,j)/h;
          g=0;
          for(k=0; k<j+1; k++)
            {
            g += a(j,k)*a(i,k);
            }
          for(k=j+1; k<l+1; k++)
            {
            g += a(k,j)*a(i,k);
            }
          e(j,0) = g/h;
          f += e(j,0)*a(i,j);
          }
        hh = f/(h+h);
        for(j=0; j<l+1; j++)
          {
          f=a(i,j);
          e(j,0)=g=e(j,0)-hh*f;
          for(k=0; k<j+1; k++)
            {
            a(j,k) -= (f*e(k,0) + g*a(i,k));
            }
          }
        }
      }
    else
      {
      e(i,0) = a(i,l);
      }
    d(i,0)=h;
    }
  d(0,0)=0;
  e(0,0)=0;
  for(i=0; i<n; i++)
    {
    l=i;
    if(d(i,0) != 0)
      {
      for(j=0; j<l; j++)
        {
        g=0;
        for(k=0; k<l; k++)
          {
          g += a(i,k)*a(k,j);
          }
        for(k=0; k<l; k++)
          {
          a(k,j) -= g*a(k,i);
          }
        }
      }
    d(i,0) = a(i,i);
    a(i,i) = 1;
    for(j=0; j<l; j++)
      {
      a(j,i) = a(i,j) = 0;
      }
    }
  }

bool eigentridiag(statmatrix<double> & d, statmatrix<double> & e,
                  statmatrix<double> & z)
  {
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;

  int n=d.rows();
  for(i=1; i<n; i++)
    {
    e(i-1,0)=e(i,0);
    }
  e(n-1,0)=0;

  for(l=0; l<n; l++)
    {
    iter=0;
    do
      {
      for(m=l; m<n-1; m++)
        {
        dd=fabs(d(m,0))+fabs(d(m+1,0));
        if(fabs(e(m,0))+dd == dd)
          {
          break;
          }
        }
      if(m != l)
        {
        if(iter++ == 30)
          {
          return false;
          }
        g=(d(l+1,0)-d(l,0))/(2*e(l,0));
        r=pythag(g,1);
        g=d(m,0)-d(l,0)+e(l,0)/(g+SIGN(r,g));
        s=c=1;
        p=0;
        for(i=m-1; i>=l; i--)
          {
          f=s*e(i,0);
          b=c*e(i,0);
          e(i+1,0)=(r=pythag(f,g));
          if(r==0)
            {
            d(i+1,0) -= p;
            e(m,0)=0;
            break;
            }
          s=f/r;
          c=g/r;
          g=d(i+1,0)-p;
          r=(d(i,0)-g)*s+2*c*b;
          d(i+1,0)=g+(p=s*r);
          g=c*r-b;
          for(k=0; k<n; k++)
            {
            f=z(k,i+1);
            z(k,i+1)=s*z(k,i)+c*f;
            z(k,i)=c*z(k,i)-s*f;
            }
          }
        if(r==0 && i>=1)
          {
          continue;
          }
        d(l,0) -= p;
        e(l,0) = g;
        e(m,0) = 0;
        }
      }
    while (m!=l);
    }
  return true;
  }

double pythag(const double & a, const double & b)
  {
  double absa, absb;

  absa=fabs(a);
  absb=fabs(b);
  if(absa > absb)
    {
    return absa*sqrt(1+sqr(absb/absa));
    }
  else
    {
    return (absb == 0 ? 0 : absb*sqrt(1+sqr(absa/absb)));
    }
  }
double sqr(const double & a)
  {
  return a*a;
  }
double SIGN(const double & a, const double & b)
  {
  return b>=0 ? (a>=0 ? a : -a) : (a>=0 ? -a : a);
  }


int eigen(statmatrix<double> & a, statmatrix<double> & values,
           statmatrix<double> & vectors)
  {
  assert(a.cols()==vectors.cols());
  assert(a.rows()==vectors.rows());
  assert(a.rows()==values.rows());
  assert(values.cols()==1);
  assert(a.cols()==a.rows());

  int i,j,ip,iq;
  double tresh,theta,tau,t,sm,s,h,g,c;

  int n=a.cols();

  vectors=statmatrix<double>::diag(n,1);
  statmatrix<double> b = a.diag();
  values = a.diag();
  statmatrix<double>z(n,1,0);

  for(i=1; i<=50; i++)
    {
    sm=0;
    for(ip=0; ip<n-1; ip++)
      {
      for(iq=ip+1; iq<n; iq++)
        {
        sm += fabs(a(ip,iq));
        }
      }

    if(sm==0)
      {
      return i;
      }

    if(i<4)
      {
      tresh=0.2*sm/(n*n);
      }
    else
      {
      tresh=0;
      }

    for(ip=0; ip<n-1; ip++)
      {
      for(iq=ip+1; iq<n; iq++)
        {
        g=100*fabs(a(ip,iq));
        if(i>4 && (fabs(values(ip,0))+g)==fabs(values(ip,0)) &&
                  (fabs(values(iq,0))+g)==fabs(values(iq,0)))
          {
          a(ip,iq)=0;
          }
        else if(fabs(a(ip,iq))>tresh)
          {
          h=values(iq,0)-values(ip,0);
          if((fabs(h)+g)==fabs(h))
            {
            t=a(ip,iq)/h;
            }
          else
            {
            theta=0.5*h/(a(ip,iq));
            t=1/(fabs(theta)+sqrt(1+theta*theta));
            if(theta<0)
              {
              t=-t;
              }
            }
          c=1/sqrt(1+t*t);
          s=t*c;
          tau=s/(1+c);
          h=t*a(ip,iq);
          z(ip,0) -= h;
          z(iq,0) += h;
          values(ip,0) -= h;
          values(iq,0) += h;
          a(ip,iq)=0;
          for(j=0; j<ip; j++)
            {
            rotate(a,s,tau,j,ip,j,iq);
            }
          for(j=ip+1; j<iq; j++)
            {
            rotate(a,s,tau,ip,j,j,iq);
            }
          for(j=iq+1; j<n; j++)
            {
            rotate(a,s,tau,ip,j,iq,j);
            }
          for(j=0; j<n; j++)
            {
            rotate(vectors,s,tau,j,ip,j,iq);
            }
          }
        }
      }
    for(ip=0; ip<n; ip++)
      {
      b(ip,0) += z(ip,0);
      values(ip,0) = b(ip,0);
      z(ip,0) = 0;
      }
    }
  return i;
  }

bool eigen2(statmatrix<double> & a, statmatrix<double> & d)
  {
  statmatrix<double> e(d.rows(),1,0);
  tridiag(a,d,e);
  return eigentridiag(d,e,a);
  }


void eigensort(statmatrix<double> & values, statmatrix<double> & vectors)
  {
  int i,j,k;
  double p;

  int n=values.rows();
  for(i=0; i<n-1; i++)
    {
    p=values(k=i,0);
    for(j=i; j<n; j++)
      {
      if(values(j,0)>=p)
        {
        p=values(k=j,0);
        }
      }
    if(k!=i)
      {
      values(k,0)=values(i,0);
      values(i,0)=p;
      for(j=0; j<n; j++)
        {
        p=vectors(j,i);
        vectors(j,i)=vectors(j,k);
        vectors(j,k)=p;
        }
      }
    }
  }

statmatrix<double> kronecker(const statmatrix<double> & A, const statmatrix<double> & B)
  {
  statmatrix<double> res(A.rows()*B.rows(),A.cols()*B.cols(),0);
  unsigned i,j;
  unsigned brows=B.rows();
  unsigned bcols=B.cols();
  for(i=0; i<A.rows(); i++)
    {
    for(j=0; j<A.cols(); j++)
      {
      res.putBlock(A(i,j)*B,i*brows,j*bcols,(i+1)*brows,(j+1)*bcols);
      }
    }
  return res;
  }

void compare(const datamatrix & refdata, const datamatrix & neudata, double limit,
             unsigned col, const ST::string & colname, vector<ST::string> & out)
  {

  double diff;
  datamatrix help = datamatrix(refdata.rows(),1,0.0);

  help.minus(neudata.getCol(col),refdata.getCol(col));
  diff = help.norm(0)/refdata.norm(col);

// Ausgabe

  if(diff > limit)
    out.push_back("WARNUNG:");

  out.push_back("  '" + colname + "': " + ST::doubletostring(diff,4) );

  if(diff > limit)
    out.push_back("  Toleranzgrenze: " + ST::doubletostring(limit) );

  }


void compare_nonp(const ST::string & ref, const ST::string & neu, double limit,
                  vector<ST::string> & out)
  {

/* pmean   pqu2p5   pqu10   pmed   pqu90   pqu97p5 */
// pmode   ci95lower   ci80lower   std   ci80upper   ci95upper
// pstd
// variance
// sigma2
// linpred   mu   saturated_deviance   leverage
// eta

  datamatrix refdata;
  datamatrix neudata;

//  Dateien lesen

  ST::string header;

  ifstream in(ref.strtochar());
  if(!in.fail())
    {
    ST::getline(in,50000,header);
    refdata.prettyScan(in);
    }

  ifstream in2(neu.strtochar());
  if(!in2.fail())
    {
    ST::getline(in2,50000,header);
    neudata.prettyScan(in2);
    }

  list<ST::string> colnames = header.strtokenlist(" \t",false);

// diff berechnen

  out.push_back("Relative quadratische Abweichung zur Referenz in der Datei");
  out.push_back("  '" + neu + "':");
  out.push_back("\n");

  unsigned i = 0;
  list<ST::string>::iterator it = colnames.begin();
  while(it != colnames.end())
    {
    if(*it == "pmean")
      compare(refdata,neudata,limit,i,*it,out);
    i++;
    it++;
    }

//  compare(refdata,neudata,limit,2,"pmean",out);
  compare(refdata,neudata,limit,5,"pmed",out);

  compare(refdata,neudata,2.0*limit,4,"pqu10",out);
  compare(refdata,neudata,2.0*limit,6,"pqu90",out);

  compare(refdata,neudata,2.5*limit,3,"pqu2p5",out);
  compare(refdata,neudata,2.5*limit,7,"pqu97p5",out);

  out.push_back("\n");

  }

