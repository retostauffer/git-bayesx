//---------------------------------------------------------------------------
#ifndef IWLS_psplineH
#define IWLS_psplineH

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif

#include<mcmc.h>
#include<fullcond.h>
#include "time.h"
#include <deque>
#include <sparsemat.h>
#include<mcmc_nonpbasis.h>
#include<spline_basis.h>
#include<fullcond_nonp_gaussian.h>

namespace MCMC
{


//------------------------------------------------------------------------------
//---------------------------- class: IWLS_pspline -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE IWLS_pspline : public spline_basis
  {

  protected:

  updatetype utype;

  double a_invgamma;
  double b_invgamma;
  double kappa;
  double kappaprop;
  double kappamode;
  double kappamean;

  bool diagtransform;

  unsigned updateW;

  void create_iwls(void);

  void update_IWLS(void);

  void update_IWLS_mode(void);

  void update_IWLS_hyperblock(void);

  void update_IWLS_hyperblock_mode(void);

  void update_isotonic(void);

  void update_diagtransform(void);


  public:

  // DEFAULT CONSTRUCTOR

  IWLS_pspline(void) : spline_basis()
    {
    }

  // CONSTRUCTOR 1

  IWLS_pspline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & d,const bool & mode,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l, const fieldtype & ft, const ST::string & monotone,
                    const unsigned & upW, const bool & updatetau, const double & fstart,
                    const double & a, const double & b, const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const bool & diag, const unsigned & c);

  // CONSTRUCTOR 2 (f�r variierende Koeffizienten)

  IWLS_pspline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & effmod,const datamatrix & intact,const bool & mode,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l, const fieldtype & ft, const ST::string & monotone,
                    const unsigned & upW, const bool & updatetau, const double & fstart,
                    const double & a, const double & b, const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const bool & diag, const unsigned & c);

  // CONSTRUCTOR 3 (f�r Cox)

  IWLS_pspline(MCMCoptions * o, DISTRIBUTION * dp,FULLCOND_const * fcc,
                const fieldtype & ft,const ST::string & ti,
                const unsigned & nrk, const unsigned & degr, const MCMC::knotpos & kp,
                const int & gs, const ST::string & fp,
                const ST::string & pres, const bool & deriv, const unsigned & c);

  // COPY CONSTRUCTOR

  IWLS_pspline(const IWLS_pspline & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const IWLS_pspline & operator=(const IWLS_pspline & fc);

  void update(void);

  void outresults(void);

  void outoptions(void);

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    f = 100.0;
    oldacceptance = 0;
    oldnrtrials = 0;
    FULLCOND_nonp_basis::reset();
    }

  // FUNCTION: predict
  // TASK: predicts the mean for a new observation Xnew

  void predict(const datamatrix & newX, datamatrix & linpred);

  // FUNCTION: compute_quadform
  // TASK: returns beta(.,v)' K beta(.,v) where K is the penalty matrix

  double compute_quadform(void);


  // DESTRUCTOR

  ~IWLS_pspline() {}

  };


}   // end: namespace MCMC


//---------------------------------------------------------------------------
#endif
