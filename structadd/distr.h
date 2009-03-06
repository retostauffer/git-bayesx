

#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (DISTR_INCLUDED)
#define DISTR_INCLUDED

#include"statmat.h"
#include"Random.h"
#include"GENERAL_OPTIONS.h"
#include"FC.h"


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

/*
1. workingweights ändern sich, weights ungleich 1

2. workingweights ändern sich, weights gleich eins

3. workingweights ändern sich nicht und sind konstant

4. workingweights ändern sich nicht und sind eins
*/

enum weighttype{wweightschange_weightsneqone,wweightschange_weightsone,
wweightsnochange_constant,wweightsnochange_one};

class __EXPORT_TYPE DISTR
  {

  protected:

  // FUNCTION: check_workingweights_one
  // TASK: checks if all workingweights are one (returns true if this is the
  //       case)

  bool check_weightsone(void);


  GENERAL_OPTIONS * optionsp;         // pointer to general MCMC options object


  public:

  bool optionbool1;
  ST::string option1;


  double sigma2;

  bool updateIWLS;
  ST::string family;              // name of the distribution

  unsigned nrobs;                 // Number of observations

  datamatrix response;                // Response
  datamatrix response_untransformed;  // untransformed response, i.e.
                                      // response in the original measurement
                                      // scale
  datamatrix workingresponse;         // Working response, tilde y
  ST::string responsename;            // Name of the response


  datamatrix weight;              // Weightvariable for weighted regression
  ST::string weightname;          // Name of the weightvariable

  datamatrix workingweight;       // Working weight (workingweight = weight
                                  // in the constructor)

  weighttype wtype;               // weight type: default is
                                  // wweightschange_weightsneqone, i.e.
                                  // workingweights change and weights are
                                  // not equal to one
  bool weightsone;                // true if weights are one for all
                                  // observations


  datamatrix linearpred1;          // Linear predictor
  datamatrix linearpred2;          // Proposed linear predictor
  int linpred_current;

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

  virtual double loglikelihood_weightsone(double * res,double * lin) const
    {
    return 0;
    }

  // FUNCTION: loglikelihood
  // TASK: computes the complete loglikelihood for all observations

  double loglikelihood(const bool & current=true) const;

  // FUNCTION: loglikelihood
  // TASK: computes the loglikelihood for observations between begin and end
  //       response, weights, predicor stored in responsep,workingweightp,
  //       linpredp

  double loglikelihood(int & begin,
                       int & end, statmatrix<double *> & responsep,
                       statmatrix<double *> & workingweightp,
                       statmatrix<double *> & linpredp) const;


  //----------------------------------------------------------------------------
  //------------------------------- COMPUTE mu ---------------------------------
  //----------------------------------------------------------------------------

  virtual void compute_mu(const double * linpred,double * mu,
                          bool notransform=false);

  virtual void compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           double * scale) const;


  //----------------------------------------------------------------------------
  //----------------------------- IWLS Algorithm -------------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: compute_IWLS (for one observation)
  // TASK: computes the iwls weights (will be stored in workingweight),
  //       tildey=predicor+(y-mu)g'(mu) (stored in workingresponse) and
  //       the loglikelihood (will be returned)

  //       type: wweightschange_weightsneqone

  virtual double compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse,const bool & like)
    {
    return 0;
    }

  // FUNCTION: compute_IWLS (for one observation)
  // TASK: computes the iwls weights (will be stored in workingweight),
  //       tildey=predicor+(y-mu)g'(mu) (stored in workingresponse) and
  //       the loglikelihood stored in like (only if compute_like = true)
  //       assumes that weighs=1 (for all observations)

  virtual void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
    {
    }


  // FUNCTION: compute_IWLS (for one observation)
  // TASK: computes tildey=predicor+(y-mu)g'(mu) (stored in workingresponse) and
  //       the loglikelihood stored in like (only if compute_like = true)
  //       assumes that workingweighs=constant (for all observations), i.e.
  //       they are not recomputed in the function

  //       wweightsnochange_constant

  virtual void compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
    {
    }


  // FUNCTION: compute_IWLS (for one observation)
  // TASK: computes tildey=predicor+(y-mu)g'(mu) (stored in workingresponse) and
  //       the loglikelihood stored in like (only if compute_like = true)
  //       assumes that workingweighs=1 (for all observations), must be set
  //       to one in advance

  virtual void compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
    {
    }


  // FUNCTION: compute_IWLS (for the whole dataset
  // TASK:

  double compute_iwls(const bool & current,const bool & like);


  // FUNCTION: compute_IWLS (
  // TASK: computes the iwls weights (will be stored in workingweight),
  //       tildey=(y-mu)g'(mu) (stored in workingresponse) and
  //       the loglikelihood (will be returned) for the begin - end observation
  //       in the pointer vectors

  double compute_iwls_loglikelihood(int & begin,
                                 int & end, statmatrix<double *> & responsep,
                                 statmatrix<double *> & workingresponsep,
                                 statmatrix<double *> & weightp,
                                 statmatrix<double *> & workingweightp,
                                 statmatrix<double *> & linpredp);


  double compute_iwls_loglikelihood_sumworkingweight(
         int & begin,int & end, statmatrix<double *> & responsep,
         statmatrix<double *> & workingresponsep,statmatrix<double *> & weightp,
         statmatrix<double *> & workingweightp, statmatrix<double *> & linpredp,
         datamatrix & intvar2,double & sumworkingweight);


  //----------------------------------------------------------------------------
  //----------------------- ACCESS TO SCALE PARAMETER --------------------------
  //----------------------------------------------------------------------------

  virtual double get_scale(bool tranform=false);
  virtual double get_scalemean(void);

  //----------------------------------------------------------------------------
  //----------------------- POSTERIORMODE FUNCTIONS ----------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  virtual bool posteriormode(void);


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

  // FUNCTION: standardise
  // TASK: standardises the response and the offset
  //       sets scalesave.transform = trmult*trmult (!!!)

  virtual void standardise(void);


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

  void compute_mu(const double * linpred,double * mu, bool notransform);

  void compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           double * scale) const;

  double loglikelihood(double * res,
                       double * lin,
                       double * w) const;

  double loglikelihood_weightsone(double * res,double * lin) const;

  double compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse, const bool & like);

  void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like);

  void compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  bool posteriormode(void);

  void outresults(ST::string pathresults="");

  double get_scalemean(void);


  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_gaussian_exp ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gaussian_exp : public DISTR_gaussian
  {

  protected:

  // FUNCTION: standardise
  // TASK: standardises the response and the offset
  //       sets scalesave.transform = trmult*trmult (!!!)

  void standardise(void);


  public:

   // DEFAULT CONSTRUCTOR

   DISTR_gaussian_exp(void) : DISTR_gaussian()
     {
     }

   // CONSTRUCTOR1
   // a_invgamma = a
   // b_invgamma = b

   DISTR_gaussian_exp(const double & a,const double & b,GENERAL_OPTIONS * o,
                  const datamatrix & r,const ST::string & ps,
                         const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_gaussian_exp(const DISTR_gaussian_exp & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_gaussian_exp & operator=(const DISTR_gaussian_exp & nd);

   // DESTRUCTOR

   ~DISTR_gaussian_exp() {}

  void compute_mu(const double * linpred,double * mu, bool notransform);


  double loglikelihood(double * res,
                       double * lin,
                       double * w) const;

  double compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse, const bool & like);

  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  bool posteriormode(void);

  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_gaussian_mult -----------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gaussian_mult : public DISTR_gaussian_exp
  {

  protected:

  // FUNCTION: standardise
  // TASK: standardises the response and the offset
  //       sets scalesave.transform = trmult*trmult (!!!)

  void standardise(void);


  public:

  void set_mult(bool & m);

   // DEFAULT CONSTRUCTOR

   DISTR_gaussian_mult(void) : DISTR_gaussian_exp()
     {
     }

   // CONSTRUCTOR1
   // a_invgamma = a
   // b_invgamma = b

   DISTR_gaussian_mult(const double & a,const double & b,GENERAL_OPTIONS * o,
                  const datamatrix & r,const ST::string & ps,
                         const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_gaussian_mult(const DISTR_gaussian_mult & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_gaussian_mult & operator=(const DISTR_gaussian_mult & nd);

   // DESTRUCTOR

   ~DISTR_gaussian_mult() {}

  void compute_mu(const double * linpred,double * mu, bool notransform);


  double loglikelihood(double * res,
                       double * lin,
                       double * w) const;

  double compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse, const bool & like);

  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  bool posteriormode(void);

  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_gaussian_re -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gaussian_re : public DISTR_gaussian
  {

  protected:

  public:

   // DEFAULT CONSTRUCTOR

   DISTR_gaussian_re(void) : DISTR_gaussian()
     {
     }

   // CONSTRUCTOR1

   DISTR_gaussian_re(GENERAL_OPTIONS * o,
                  const datamatrix & r,
                  const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_gaussian_re(const DISTR_gaussian_re & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_gaussian_re & operator=(const DISTR_gaussian_re & nd);

   // DESTRUCTOR

   ~DISTR_gaussian_re() {}

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  bool posteriormode(void);

  void outresults(ST::string pathresults="");

  void outoptions(void);

  };


} // end: namespace MCMC


#endif
