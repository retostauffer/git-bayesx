
#if !defined (TVARIANCE_INCLUDED)
#define TVARIANCE_INCLUDED

#include"../export_type.h"
#include"fullcond_nonp_gaussian.h"

namespace MCMC
{



class __EXPORT_TYPE FULLCOND_tvariance : public FULLCOND
  {


  protected:

  FULLCOND_nonp_basis * Kp;

  datamatrix u;

  unsigned nu;

  unsigned start;

  public:


  // DEFAULT CONSTRUCTOR

  FULLCOND_tvariance(void) : FULLCOND()
    {
    }


  // CONSTRUCTOR

  FULLCOND_tvariance(MCMCoptions * o,FULLCOND_nonp_basis * p,unsigned & v,
                     const ST::string & ti, const ST::string & fp,
                     const ST::string & pres);


  // COPY CONSTRUCTOR

  FULLCOND_tvariance(const FULLCOND_tvariance & t)
    : FULLCOND(FULLCOND(t))
    {
    Kp = t.Kp;
    u = t.u;
    nu = t.nu;
    start=t.start;
    }

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_tvariance & operator=(const FULLCOND_tvariance & t)
    {
    if (this == &t)
      return *this;
    FULLCOND::operator=(FULLCOND(t));
    Kp = t.Kp;
    u = t.u;
    nu = t.nu;
    start=t.start;
    return *this;
    }


  void update(void);


  bool posteriormode(void)
    {
    return true;
    }

  void outresults(void);

  void outoptions(void);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    FULLCOND::reset();
    setbeta(nrpar,1,1);
    }

  // FUNCTION: predict (virtual)
  // TASK: predicts the mean for a new observation Xnew

  void predict(const datamatrix & newX, datamatrix & linpred)
    {

    }

  // DESTRUCTOR

  ~FULLCOND_tvariance() {}

  }; // end: class FULLCOND_tvariance






} // end: namespace MCMC

#endif





