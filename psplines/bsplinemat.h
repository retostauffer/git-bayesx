
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE
#else
#define __EXPORT_TYPE __import
#endif

#ifndef bsplinematH
#define bsplinematH

#include<deque.h>
#include "statmat.h"
#include "mcmc_nonpbasis.h"

namespace MCMC
{


//---------------------------------------------------------------------------
//----------------------- class: bsplinemat -------------------------------
//---------------------------------------------------------------------------


class __EXPORT_TYPE bsplinemat
  {

  protected:

  datamatrix B;
  datamatrix BS;

  unsigned nrknots;
  unsigned degree;
  unsigned nrdiffobs;
  unsigned nrpar;

  knotpos knpos;

  vector<int> freq;
  vector<int> freqoutput;
  vector<int> index2;
  vector<int> begcol;

  deque<int> firstnonzero;
  deque<int> lastnonzero;
  deque<double> knot;

  statmatrix<int> index;

  datamatrix Bcolmean;


  void make_index(const datamatrix & md);

  void make_Bspline(const bool & deriv, const datamatrix & md, const bool & minnull = false);

  datamatrix bspline(const double & x);
  datamatrix bspline_derivative(const double & x);


  public:

  // DEFAULT CONSTRUCTOR

  bsplinemat(void)
    {
    }

  // CONSTRUCTOR

  bsplinemat(const datamatrix & data, const unsigned & nrk, const unsigned & degr,
                const MCMC::knotpos & kp, const bool & minnull = false, const deque<double> & k = deque<double>());

  // CONSTRUCTOR 2 (for derivatives)

  bsplinemat(const bool & deriv, const datamatrix & data, const unsigned & nrk, const unsigned & degr,
                const MCMC::knotpos & kp, const bool & minnull = false, const deque<double> & k = deque<double>());


  // COPY CONSTRUCTOR

  bsplinemat(const bsplinemat & bmat);

  // OVERLOADED ASSIGNMENT OPERATOR

  const bsplinemat & operator=(const bsplinemat & bmat);

  void mult(datamatrix & res, const datamatrix & beta);

  void mult_index(datamatrix & res, const datamatrix & beta);

  // DESTRUCTOR

  ~bsplinemat(){}

  };


}   // END: namespace MCMC

//---------------------------------------------------------------------------
#endif
