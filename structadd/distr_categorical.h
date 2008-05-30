

#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (DISTRcategorical_INCLUDED)
#define DISTRcategorical_INCLUDED

#include"statmat.h"
#include"random.h"
#include"GENERAL_OPTIONS.h"
#include"fc.h"
#include"distr.h"


namespace MCMC
{


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_binomial ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_binomial : public DISTR
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

   DISTR_binomial(void) : DISTR()
     {
     }

   // CONSTRUCTOR

   DISTR_binomial(GENERAL_OPTIONS * o, const datamatrix & r,
                  const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_binomial(const DISTR_binomial & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_binomial & operator=(const DISTR_binomial & nd);

   // DESTRUCTOR

   ~DISTR_binomial() {}

  double loglikelihood(double * response, double * linpred,
                       double * weight) const;

  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);

  void outoptions(void);

  };

} // end: namespace MCMC


#endif
