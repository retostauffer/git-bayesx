#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FULLCOND_mult_INCLUDED)

#define FULLCOND_mult_INCLUDED

#include "randomeffect.h"
#include "mcmc_nonpbasis.h"
#include "statmat_penalty.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_mult ---------------------------------
//------------------------------------------------------------------------------


enum multtype {re_rw,mrf_rw,re_pspline,mrf_pspline}; 

class __EXPORT_TYPE FULLCOND_mult : public FULLCOND
  {

  protected:

  FULLCOND_nonp_basis * basis1p;
  FULLCOND_nonp_basis * basis2p;
  FULLCOND_random * reffectp;

  multtype ttype;

  bool first;


  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_mult(void) : FULLCOND()
    {
    }

  // CONSTRUCTOR 1  (for i.i.d. RE * nonlinear)
  // o    : pointer to MCMCoptions object
  // dp   : pointer to distribution object
  // fp   : file where sampled parameters are stored
  // pres : file where results are stored

  FULLCOND_mult(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_random * rp,
                         FULLCOND_nonp_basis * ba,
                         bool fi,
                         const ST::string & ti,
                         const ST::string & fp, const ST::string & pres,
                         const unsigned & c);

  // CONSTRUCTOR 2  (for i.i.d. mrf * nonlinear)
  // o    : pointer to MCMCoptions object
  // dp   : pointer to distribution object
  // fp   : file where sampled parameters are stored
  // pres : file where results are stored

  FULLCOND_mult(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_nonp_basis * ba1,
                         FULLCOND_nonp_basis * ba2,
                         bool fi,
                         const ST::string & ti,
                         const ST::string & fp, const ST::string & pres,
                         const unsigned & c);


  // COPY CONSTRUCTOR

  FULLCOND_mult(const FULLCOND_mult & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_mult & operator=(const FULLCOND_mult & fc);

  void update(void);

//  void update_linpred(const bool & add)
//    {
//    }


  bool posteriormode(void);

//  bool posteriormode_converged(const unsigned & itnr);

  void outresults(void);

  void get_effectmatrix(datamatrix & e,vector<ST::string> & enames,
                        unsigned be, unsigned en,effecttype t);


  unsigned get_nreffects(effecttype t);


  void outoptions(void);

  ST::string getinfo(void);

  void init_name(const ST::string & na);

  void init_names(const vector<ST::string> & na);

  void init_priorassumptions(const ST::string & na);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void);

  // DESTRUCTOR

  ~FULLCOND_mult() {}


  };


} // end: namespace MCMC

#endif

 