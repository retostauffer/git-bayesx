
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (MCMC_RIDGE_INCLUDED)

#define MCMC_RIDGE_INCLUDED

#include "fullcond.h"
#include "mcmc_const.h"

namespace MCMC
{
class __EXPORT_TYPE FULLCOND_ridge : public FULLCOND
  {

  protected:

  vector<double> variances;  // vector of variances for the ridge penalty
  datamatrix X1;                   // (X'WX)^-0.5
  datamatrix X2;                   // (X'WX)^-1X'W

  datamatrix linold;
  datamatrix mu1;

  DISTRIBUTION * likep;

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FULLCOND_ridge(void) : FULLCOND()
    {
    }

  // CONSTRUCTOR

  FULLCOND_ridge(MCMCoptions * o, DISTRIBUTION * dp, const datamatrix & d,
                 const ST::string & t, const ST::string & fs,
                 const ST::string & fr, const vector<double> & vars,
                 const unsigned & c);

  // COPY CONSTRUCTOR

  FULLCOND_ridge(const FULLCOND_ridge & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_ridge & operator=(const FULLCOND_ridge & m);

  // DESTRUCTOR

  ~FULLCOND_ridge()
    {
    }

//-------------------------- UPDATE and related methods-------------------------

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  // FUNCTION: outresults
  // TASK: - write results to output window and files

  void outresults(void);

  };

} // end: namespace MCMC

#endif

