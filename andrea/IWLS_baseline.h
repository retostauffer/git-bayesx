//---------------------------------------------------------------------------
#ifndef IWLS_baselineH
#define IWLS_baselineH

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
//---------------------------- class: IWLS_baseline -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE IWLS_baseline : public spline_basis
  {

  protected:

     bool begin0;
   datamatrix int_knots;
   datamatrix int_D;
   MCMC::bsplinemat testmat;
   vector<MCMC::bsplinemat> gaussmat;
   vector<IWLS_baseline*> baselinep;
   datamatrix zi;
   unsigned gauss_n;
   datamatrix coeff;
   datamatrix z_vc;
   datamatrix zi_ges;
   datamatrix beg_i;
   statmatrix<int> zi_index;
   statmatrix<int> ges_index;
   datamatrix spline_ges;
   datamatrix spline_zi;
   datamatrix gaussspline;
   datamatrix int_ti_help;
   bool vc_dummy1;

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

  IWLS_baseline(void) : spline_basis()
    {
    }

  // CONSTRUCTOR 1

  IWLS_baseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & d,const bool & mode,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l, const fieldtype & ft, const ST::string & monotone,
                    const unsigned & upW, const bool & updatetau, const double & fstart,
                    const double & a, const double & b, const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const bool & diag, const unsigned & c, const datamatrix & anfang);

  // CONSTRUCTOR 2 (für variierende Koeffizienten)

  IWLS_baseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & effmod,const datamatrix & intact,const bool & mode,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l, const fieldtype & ft, const ST::string & monotone,
                    const unsigned & upW, const bool & updatetau, const double & fstart,
                    const double & a, const double & b, const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const bool & diag, const unsigned & c, const datamatrix & anfang);

  // COPY CONSTRUCTOR

  IWLS_baseline(const IWLS_baseline & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const IWLS_baseline & operator=(const IWLS_baseline & fc);

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

  void set_fcconst(FULLCOND_const * fcc)
    {
    fcconst = fcc;
    }

  // DESTRUCTOR

  ~IWLS_baseline() {}

  };


}   // end: namespace MCMC


//---------------------------------------------------------------------------
#endif
