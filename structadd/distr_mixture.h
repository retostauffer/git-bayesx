

#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (DISTRmixture_INCLUDED)
#define DISTRmixture_INCLUDED

#include"statmat.h"
#include"Random.h"
#include"GENERAL_OPTIONS.h"
#include"FC.h"
#include"distr.h"


namespace MCMC
{

//------------------------------------------------------------------------------
//-------------------- CLASS: DISTR_gaussianmixture ----------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE DISTR_gaussianmixture : public DISTR_gaussian
  {

  protected:

  void define_knots(double left, double dist);

  unsigned nrknots;

  datamatrix alpha;          // transformed weights according to logistic
                             // transformation
  FC alphasample;

  datamatrix alpha_prob;     // probabilities

  datamatrix m;              // knots
  double s2;                  // variance of mixture components

  statmatrix<int> rho;         // index of mixture component




  public:

   // DEFAULT CONSTRUCTOR

   DISTR_gaussianmixture(void) : DISTR_gaussian()
     {
     }

   // CONSTRUCTOR1
   // a_invgamma = a
   // b_invgamma = b

   DISTR_gaussianmixture(const double & a,const double & b,GENERAL_OPTIONS * o,
                  const datamatrix & r,const ST::string & ps,
                  const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_gaussianmixture(const DISTR_gaussianmixture & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_gaussianmixture & operator=(const DISTR_gaussianmixture & nd);

   // DESTRUCTOR

   ~DISTR_gaussianmixture() {}

/*
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

  */

  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the scale parameter, alpha's

  void update(void);


  bool posteriormode(void);

  /*
  void outresults(ST::string pathresults="");

  double get_scalemean(void);

  void sample_responses(unsigned i,datamatrix & sr);

  void outresults_predictive_check(datamatrix & D,datamatrix & sr);
  */
  };






} // end: namespace MCMC


#endif
