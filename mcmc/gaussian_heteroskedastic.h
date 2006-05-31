
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


  public:

  // DEFAULT CONSTRUCTOR

  DISTRIBUTION_gaussianh(void) : DISTRIBUTION() {}

   // CONSTRUCTOR

   DISTRIBUTION_gaussianh(const double & a,const datamatrix & b,
                   MCMCoptions * o, const datamatrix & r,
                   const ST::string & fp,const ST::string & fs,
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
                       double * weight,const int & i) const
      {

      /*
      datamatrix V;
      datamatrix likeli;
      V=Sigma.inverse();
      likeli = 0.5*nrobs * V.det() * log(V) - 0.5 * diff.transposed()*V*diff;

       */
    // include program code to compute the loglikelihood for the distribution

    return 0;
    }

  // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred' and stores
  //       the result in 'mu'

  void compute_mu(const double * linpred,double * mu) const;

  // FUNCTION: compute_deviance

  void compute_deviance(const double * response, const double * weight,
                        const double * mu,double * deviance,
                        double * deviancesat,
                        const datamatrix & scale,const int & i) const;


  // FUNCTION: outoptions
  // TASK: writing options of the distribution

  void outoptions(void)
    {
    DISTRIBUTION::outoptions();

    optionsp->out("\n");

    }


  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  void outresults(void)
    {
    DISTRIBUTION::outresults();
    }

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode for the scale parameter and the
  //       intercept

  bool posteriormode(void);

  bool posteriormode_converged_fc(const datamatrix & beta,
                                  const datamatrix & beta_mode,
                                  const unsigned & itnr)
    {
    return true;
    }

  void compute_iwls(void)
    {
    tildey.assign(response);
    }


  };


 }




 // end: namespace MCMC

#endif

