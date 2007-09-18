

#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (MCMCadaptiv_INCLUDED)

#define MCMCadaptiv_INCLUDED

#include"mcmc.h"
#include"fullcond.h"
#include "mcmc_nonp.h"
#include "mcmc_nonpbasis.h"
#include "variance_nonp.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- class: FULLCOND_adaptiv --------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_adaptiv : public FULLCOND
  {

  protected:
  fieldtype type;                    // order of variance random walk (RW1 or RW2)
  long double sigma2;                // variance parameter sigma_j^2 for
                                     // variance random walk
  double a_invgamma;                 // Hyperparameter a_j' of sigma_j^2
  double b_invgamma;                 // Hyperparameter b_j' of sigma_j^2
  double sigma2sum;
  bool unifb;

  unsigned minblocksize;             // Minimum blocksize for blockmove
  unsigned maxblocksize;             // Maximum blocksize for blockmove

  unsigned start;

  FULLCOND_nonp_basis * Gp;          // pointer

  ST::string pathresults;            // path for variance-parameter-results

  datamatrix h;                      // matrix of h_t's
  PenaltyMatrix  Pm;                 // Penalty matrix object

   // FUNCTION: compute_denquot
   // TASK: computes ln{P(beta_t|beta_[s<t],h*_t)/P(beta_t|beta_[s<t],h_t)}

  double compute_denquot(unsigned i,double hp);

//  double compute_hprop(unsigned & i);

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_adaptiv(void) : FULLCOND()
    {
    }

  FULLCOND_adaptiv(MCMCoptions * o,FULLCOND_nonp_basis * p,
                   const fieldtype & ft,
                   const ST::string & ti, const double & a, const double & b,
                   const bool & uniformb,
                   const double & startvalue,const unsigned & minb,
                   const unsigned & maxb, const ST::string & fp,
                   const ST::string & pres);

  // COPY CONSTRUCTOR

  FULLCOND_adaptiv(const FULLCOND_adaptiv & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_adaptiv & operator=(const FULLCOND_adaptiv & fc);

  void update(void);

  void outresults(void);

  void outoptions(void);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    FULLCOND::reset();
    sigma2 = 1;
    sigma2sum = 0;
    }

  // FUNCTION: predict (virtual)
  // TASK: predicts the mean for a new observation Xnew

  void predict(const datamatrix & newX, datamatrix & linpred)
    {
    }


  // DESTRUCTOR

  ~FULLCOND_adaptiv() {}

  };


}   // end: namespace MCMC

#endif



