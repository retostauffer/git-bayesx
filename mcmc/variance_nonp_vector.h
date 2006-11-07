
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (VARIANCENONP_VECTOR_INCLUDED)
#define VARIANCENONP_VECTOR_INCLUDED

#include"mcmc_const.h"

namespace MCMC
{


class __EXPORT_TYPE FULLCOND_variance_nonp_vector : public FULLCOND
  {

  protected:

  FULLCOND_const * Cp;

  bool update_sigma2;

  DISTRIBUTION * distrp;

  vector<double> a_invgamma;
  vector<double> b_invgamma;

  ST::string pathresults;

  unsigned column;

  vector<double> tau;
  vector<double> lambda;

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_variance_nonp_vector(void) : FULLCOND()
    {
    }


  // CONSTRUCTOR1

  FULLCOND_variance_nonp_vector(MCMCoptions * o, FULLCOND_const * p,
                         DISTRIBUTION * d, const vector<double> & a,
                         const vector<double> & b, const ST::string & ti,
                         const ST::string & fp, const ST::string & fr,
                         const unsigned & c);

  // COPY CONSTRUCTOR

  FULLCOND_variance_nonp_vector(const FULLCOND_variance_nonp_vector & t);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_variance_nonp_vector & operator=(const FULLCOND_variance_nonp_vector & t);

  void update(void);

  void outresults(void);

  void outoptions(void);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    FULLCOND::reset();
    setbeta(beta.rows(),1,0.1);
    }

  // DESTRUCTOR

  ~FULLCOND_variance_nonp_vector() {}

  }; // end: class FULLCOND_variance_nonp_vector



} // end: namespace MCMC

#endif


