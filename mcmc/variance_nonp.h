#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (VARIANCENONP_INCLUDED)
#define VARIANCENONP_INCLUDED

#include<mcmc_nonpbasis.h>
#include<randomeffect.h>
#include<mcmc_nonp.h>

namespace MCMC
{


class __EXPORT_TYPE FULLCOND_variance_nonp : public FULLCOND
  {

  protected:

  bool constlambda;
  bool uniformprior;

  FULLCOND_nonp_basis * Kp;

  FULLCOND_random * REp;
  bool randomeffect;

  FULLCOND_nonp * Fnp;
  bool fullcondnonp;

  bool update_sigma2;

  DISTRIBUTION * distrp;

  double a_invgamma;
  double b_invgamma;
  unsigned rankK;

  ST::string pathresults;

  bool average;

  unsigned column;

  public:


  // DEFAULT CONSTRUCTOR

  FULLCOND_variance_nonp(void) : FULLCOND()
    {
    }


  // CONSTRUCTOR1

  FULLCOND_variance_nonp(MCMCoptions * o,FULLCOND_nonp_basis * p,
                         DISTRIBUTION * d,
                         const double & a,const double & b,
                         const ST::string & ti, const ST::string & fp,
                         const ST::string & fr,const bool & av,
                         const unsigned & c);

  // CONSTRUCTOR2

  FULLCOND_variance_nonp(MCMCoptions * o,FULLCOND_random * p,
                         DISTRIBUTION * d,
                         const double & a,const double & b,
                         const ST::string & ti, const ST::string & fp,
                         const ST::string & fr,const bool & av,
                         const unsigned & c);

  // CONSTRUCTOR3

  FULLCOND_variance_nonp(MCMCoptions * o,FULLCOND_nonp * p,
                         DISTRIBUTION * d,
                         const double & a,const double & b,
                         const ST::string & ti, const ST::string & fp,
                         const ST::string & fr,const bool & av,
                         const unsigned & c);


  // COPY CONSTRUCTOR

  FULLCOND_variance_nonp(const FULLCOND_variance_nonp & t);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_variance_nonp & operator=(const FULLCOND_variance_nonp & t);

  void update(void);

  bool posteriormode(void);

  void outresults(void);

  void outoptions(void);

  void set_update_sigma2(void)
    {
    update_sigma2 = false;
    average = false;
    }

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    FULLCOND::reset();
    setbeta(1,1,0.1);
    }

  // FUNCTION: predict (virtual)
  // TASK: predicts the mean for a new observation Xnew

  void predict(const datamatrix & newX, datamatrix & linpred)
    {
    }

  // FUNCTION: set_constlambda
  // TASK: variance/smoothingparameter is held constant

  void set_constlambda(void)
    {
    constlambda=true;
    }

  void set_uniformprior(void)
    {
    uniformprior=true;
    }

  // DESTRUCTOR

  ~FULLCOND_variance_nonp() {}

  }; // end: class FULLCOND_variance_nonp



} // end: namespace MCMC

#endif


