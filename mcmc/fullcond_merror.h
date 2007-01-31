
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FULLCOND_MERROR_INCLUDED)

#define FULLCOND_MERROR_INCLUDED

#include "fullcond.h"
#include "mcmc_nonpbasis.h"
#include "fullcond_nonp_gaussian.h"
#include "spline_basis.h"

namespace MCMC
{
class __EXPORT_TYPE fullcond_merror : public FULLCOND
  {

  protected:

  FULLCOND_nonp_gaussian * designp;
  DISTRIBUTION * likep;

// BEGIN: merror
  spline_basis * splinep;
  bool varcoeff;

  double maxx;
  double minx;

  datamatrix meandata;      // mean of the observed covariate values
  datamatrix old;           // sampled values from the previous iteration

  datamatrix currentspline; // stores current f(x) (i.e. f(x^{old})
  datamatrix diffspline;    // stores f(x^{prop})-f(x^{old})
  datamatrix logfcold;      // full conditional evaluated for old values
  datamatrix logfcnew;      // full conditional evaluated for propsed values

  int merror;               // counts number of replicated measurements for each covariate value

  statmatrix<int> index;

  FULLCOND fc_merrorvar;    // full conditional for the variance of the measurement error
  FULLCOND fc_ximu;         // full conditional for the expectation of the true covariate values
  FULLCOND fc_xivar;        // full conditional for the variance of the true covariate values

  unsigned generrcount;      // counts the number of generated values out of range
  unsigned generrtrial;

  ST::string pathresults;
// END: merror

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  fullcond_merror(void) : FULLCOND()
    {
    }

  // CONSTRUCTOR : Susi (Measurement error in the interaction variable of a VCM)
  // o    : pointer to MCMCoptions object
  // t    : title of the full conditional (for example "fixed effects")
  //        (i.e. number of categories of the response variable)
  // fp   : file path for storing sampled parameters

  fullcond_merror(MCMCoptions * o, FULLCOND_nonp_gaussian * p, DISTRIBUTION * dp,
           const datamatrix & d, const ST::string & t, const ST::string & fp);

// BEGIN: merror
  // CONSTRUCTOR : Thomas (Measurement error in a nonparametric effect)
  fullcond_merror(MCMCoptions * o, spline_basis * p, DISTRIBUTION * dp,
           const datamatrix & d, const ST::string & t, const ST::string & fp,
           const ST::string & pres, const double & lk, const double & uk);
// END: merror

  // COPY CONSTRUCTOR

  fullcond_merror(const fullcond_merror & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const fullcond_merror & operator=(const fullcond_merror & m);

  // DESTRUCTOR

  ~fullcond_merror()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  void posteriormode_set_beta_mode(void);

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void outoptions(void);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(void);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void);

  vector<ST::string> & get_results_latex(void);

  };

} // end: namespace MCMC

#endif

