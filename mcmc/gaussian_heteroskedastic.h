
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (GAUSSIANHETEROSKEDASTIC_INCLUDED)

#define GAUSSIANHETEROSKEDASTIC_INCLUDED

#include"distribution.h"


namespace MCMC
{

class __EXPORT_TYPE DISTRIBUTION_gaussianh : public DISTRIBUTION
  {

  protected:

  unsigned nrcat;               // number of categories of the response

  // FUNCTION: standardise
  // TASK: standardises the response and the offset
  //       sets scalesave.transform = trmult*trmult (!!!)

  void standardise(void);



  public:

  // DEFAULT CONSTRUCTOR

  DISTRIBUTION_gaussianh(void) : DISTRIBUTION() {}

   // CONSTRUCTOR

   DISTRIBUTION_gaussianh(const double & a,const datamatrix & b,
                   MCMCoptions * o, const datamatrix & r,
                   const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTRIBUTION_gaussianh(const DISTRIBUTION_gaussianh & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_gaussianh &
   operator=(const DISTRIBUTION_gaussianh & nd);

   // DESTRUCTOR

   ~DISTRIBUTION_gaussianh() {}


   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * response,double * linpred,
                       double * weight,const int & i) const;

  // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred' and stores
  //       the result in 'mu'

  void compute_mu(const double * linpred,double * mu) const;

  void compute_mu_notransform(const double * linpred,double * mu) const;


  // FUNCTION: compute_deviance

  void compute_deviance(const double * response, const double * weight,
                        const double * mu,double * deviance,
                        double * deviancesat,
                        const datamatrix & scale,const int & i) const;


  double compute_weight(double * linpred, double * weight,
                        const int & i, const unsigned & col=0) const;


 double compute_IWLS(double * response,double * linpred, double * weight,
                      const int & i,double * weightiwls,double * tildey,
                      bool weightyes, const unsigned & col=0);

 void compute_IWLS_weight_tildey(double * response,double * linpred,
                              double * weight,const int & i,
                              double * weightiwls,double * tildey,
                              const unsigned & col=0);

double compute_gmu(double * linpred,const unsigned & col=0) const;


  // FUNCTION: outoptions
  // TASK: writing options of the distribution

  void outoptions(void);


  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  void update_predict(void);

  void outresults(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode for the scale parameter and the
  //       intercept

  bool posteriormode(void);

/*
  bool posteriormode_converged_fc(const datamatrix & beta,
                                  const datamatrix & beta_mode,
                                  const unsigned & itnr);
*/

  void compute_iwls(void);


  };


 }




 // end: namespace MCMC

#endif

