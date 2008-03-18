
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (DISTR_INCLUDED)
#define DISTR_INCLUDED

#include"statmat.h"
#include"random.h"
#include"GENERAL_OPTIONS.h"
#include"fc.h"


namespace MCMC
{

using randnumbers::rand_invgamma;
using randnumbers::rand_normal;
using randnumbers::uniform;
using randnumbers::trunc_normal;
using randnumbers::trunc_normal2;
using randnumbers::truncnormal;
using randnumbers::kssample;
using randnumbers::rand_gamma;



class __EXPORT_TYPE DISTR
  {

  protected:


  GENERAL_OPTIONS * optionsp;         // pointer to general MCMC options object


  public:

  ST::string family;              // name of the distribution

  unsigned nrobs;                 // Number of observations

  datamatrix response;            // Response
  ST::string responsename;        // Name of the response


  datamatrix weight;              // Weightvariable for weighted regression
  ST::string weightname;          // Name of the weightvariable

  datamatrix workingweight;       // Working weight
  bool changingweight;


  datamatrix linearpred;          // Linear predictor
  datamatrix linearpredprop;      // Proposed linear predictor
  datamatrix * linpred_current;   // Pointer that contains adress of current
                                  // predictor
  datamatrix * linpred_proposed;  // Pointer that contains adress of proposed
                                  // predictor


  // FUNCTION: swap_linearpred
  // TASK: swaps the pointers linpred_current and linpred_proposed, that is
  //       the proposed new predictor will be the current

  void swap_linearpred(void);


  double trmult;                   // multiplicative constant with which
                                   // the response has been transformed


//------------------------------------------------------------------------------
//--------------------------- CONSTRUCTORS -------------------------------------
//------------------------------------------------------------------------------

  // DEFAULT CONSTRUCTOR

  DISTR(void)
    {
    }

  // CONSTRUCTOR1
  // TASK: initializes data
  //       response = r
  //       weight = w
  //       nrobs = r.rows()

  DISTR(GENERAL_OPTIONS * o,const datamatrix & r,
               const datamatrix & w=datamatrix());


  // COPY CONSTRUCTOR

  DISTR(const DISTR & d);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR & operator=(const DISTR & d);

  // DESTRUCTOR

  ~DISTR() {}

  //----------------------------------------------------------------------------
  //------------------------------ WRITING OPTIONS -----------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: outoptions
  // TASK: writing options

  virtual void outoptions(void);


  //----------------------------------------------------------------------------
  //----------------------- COMPUTING THE LOGLIKELIHOOD ------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: loglikelihood
  // TASK: computes the loglikelihood for a single observation

  virtual double loglikelihood(double * res,double * lin,double * weight) const
    {
    return 0;
    }

  // FUNCTION: loglikelihood
  // TASK: computes the complete loglikelihood for all observations

  double loglikelihood(const bool & current=true) const;

  // FUNCTION: loglikelihood
  // TASK: computes the loglikelihood for observations between 'beg' and 'end'

  double loglikelihood(const unsigned & beg,const unsigned & end,
                       const statmatrix<int> & index,
                       const bool & current=true) ;

  //----------------------------------------------------------------------------
  //----------------------- ACCESS TO SCALE PARAMETER --------------------------
  //----------------------------------------------------------------------------

  virtual double get_scale(void);

  //----------------------------------------------------------------------------
  //----------------------- POSTERIORMODE FUNCTIONS ----------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  virtual bool posteriormode(void)
    {
    return true;
    }

  //----------------------------------------------------------------------------
  //--------------------------- UPDATE FUNCTIONS -------------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: update
  // TASK: base function for inherited classes,
  //       should update the scale parameter
  //       the base function updates the estimated mean and variance
  //       of the scale parameter only

  virtual void update(void);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: outresults
  // TASK: writes estimation results for the scale parameter
  //       estimated mean and variance

  virtual void outresults(ST::string pathresults="");

  // FUNCTION: reset
  // TASK: resets linpred (all values to 0)

  void reset(void);


  }; // end: class DISTR



//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_gaussian ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gaussian : public DISTR
  {

  protected:

  double a_invgamma;                    // hyperparameter a (for the inverse
                                        // gamma distribution of the scale
                                        // parameter, i.e. sigma^2
  double b_invgamma;                    // hyperparameter b

  FC FCsigma2;
  double sigma2;

  // FUNCTION: standardise
  // TASK: standardises the response and the offset
  //       sets scalesave.transform = trmult*trmult (!!!)

  void standardise(void);


  public:


   // DEFAULT CONSTRUCTOR

   DISTR_gaussian(void) : DISTR()
     {
     }

   // CONSTRUCTOR1
   // a_invgamma = a
   // b_invgamma = b

   DISTR_gaussian(const double & a,const double & b,GENERAL_OPTIONS * o,
                  const datamatrix & r,const ST::string & ps,
                         const datamatrix & w=datamatrix());




   // COPY CONSTRUCTOR

   DISTR_gaussian(const DISTR_gaussian & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_gaussian & operator=(const DISTR_gaussian & nd);

   // DESTRUCTOR

   ~DISTR_gaussian() {}

   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * res,
                       double * lin,
                       double * w) const;


  virtual double get_scale(void);

  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  bool posteriormode(void);

  void outresults(ST::string pathresults="");


  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_gaussian_re -------------------------
//------------------------------------------------------------------------------

/*
class __EXPORT_TYPE DISTRIBUTION_gaussian_re : public DISTRIBUTION
  {

  protected:

  double a_invgamma;                   // hyperparameter a (for the inverse
                                        // gamma distribution of the scale
                                        // parameter, i.e. sigma^2
  double b_invgamma;                   // hyperparameter b

  bool constscale;
  bool uniformprior;

  DISTRIBUTION * distrp;

  public:


   // DEFAULT CONSTRUCTOR

   DISTRIBUTION_gaussian_re(void) : DISTRIBUTION()
     {
     family = "gaussian";
     a_invgamma = 1;
     b_invgamma = 0.005;
     constscale = false;
     uniformprior = false;
     }

   // CONSTRUCTOR1
   // a_invgamma = a
   // b_invgamma = b

   DISTRIBUTION_gaussian_re(const double & a,
                         const double & b,
                         MCMCoptions * o,
                         const datamatrix & r,const ST::string & p,
                         const ST::string & ps,
                         const datamatrix & w=datamatrix());

   // CONSTRUCTOR2
   // a_invgamma = a
   // b_invgamma = b

   DISTRIBUTION_gaussian_re(const datamatrix & offset, const double & a,
                         const double & b,
                         MCMCoptions * o,
                         const datamatrix & r,
                         const ST::string & p,const ST::string & ps,
                         const datamatrix & w=datamatrix());



   // COPY CONSTRUCTOR

   DISTRIBUTION_gaussian_re(const DISTRIBUTION_gaussian_re & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_gaussian_re & operator=(const DISTRIBUTION_gaussian_re & nd);

   // DESTRUCTOR

   ~DISTRIBUTION_gaussian_re() {}

   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * res,
                       double * lin,
                       double * w,
                       const int & i) const;


  // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred'

  void compute_mu(const double * linpred,double * mu) const;

  void compute_mu_notransform(const double * linpred,double * mu) const;

  // FUNCTION: compute_deviance
  // TASK: computes the retransformed individual deviance
  //       scale and response is assumed to be NOT RETRANSFORMED
  //       but will be retransformed when computing the residual
  //       mu is assumed to be already restransformed

  void compute_deviance(const double * response,const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           const datamatrix & scale,const int & i) const;

  double compute_weight(double * linpred, double * weight,
                        const int & i, const unsigned & col=0) const;

  void compute_iwls(void)
    {
    tildey.assign(response);
    DISTRIBUTION::compute_weight(weightiwls,0);
    }

  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  void update_predict(void)
    {
    DISTRIBUTION::update_predict();
    }

  bool posteriormode(void);

  bool posteriormode_converged_fc(const datamatrix & beta,
                                  const datamatrix & beta_mode,
                                  const unsigned & itnr)
    {
    return true;
    }

  void outresults(void)
    {
    DISTRIBUTION::outresults();
    }




  void get_residuals(datamatrix & r);

  void tr_nonlinear(vector<double *> b,vector<double *> br,
                    vector<FULLCOND*> & fcp,unsigned & nr,
                    unsigned & it,ST::string & trtype);

  void set_constscale(double s);

  void undo_constscale(void);

  void set_uniformprior(void);

  void set_distrpointer(DISTRIBUTION * dp);


  };
*/


} // end: namespace MCMC


#endif
