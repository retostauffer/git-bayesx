
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (MULTGAUSSIAN_INCLUDED)

#define MULTGAUSSIAN_INCLUDED

#include"distribution.h"


namespace MCMC
{

class __EXPORT_TYPE DISTRIBUTION_multgaussian : public DISTRIBUTION
  {

   protected:

   double A;                    // Hyperparameter for Inverse Wishart
   datamatrix B;                // Hyperparameter for Inverse Wishart
   datamatrix sumB;
   datamatrix diff;

   datamatrix SIGMA_mr;         // Sigma_-r
   datamatrix SIGMA_rmr;        // Sigma_r|-r (row vector)
   datamatrix sigma_rmr;        // column vector of sigma^2_r | -r

   datamatrix offset;           // stores the offsets in a nrobs x nrcat matrix
                                // i -th row corresponds to i-th observation
                                // r - th column stores o_r|-r 's

   unsigned nrcat;              // number of categories of the response


   // FUNCTION: compute_IWproduct
   // TASK: computes 0.5*SUM_{i=1}^{n} (y_i-eta_i)(y_i-eta_i)' and stores
   //       the result in sumB
   //       y_i = (y_i1,...y_ik)'

   void compute_IWproduct(void);

   // FUNCTION: compute_SIGMA_mr
   // TASK: computes the submatrix of Sigma where the elements corresponding to
   //       the r-th component of Sigma are removed
   //       stores the INVERSE of SIGMA_mr

   void compute_SIGMA_mr(unsigned r);

   // FUNCTION: compute_SIGMA_rmr
   // TASK: computes SIGMA_r|-r as a row vector

   void compute_SIGMA_rmr(unsigned r);

  // FUNCTION: compute_sigmarmr
  // TASK: computes for r=1,...,k: sigma2_r|-r =
  //       sigma2_r - SIGMA_r|-r * (SIGMA_-r)^-1 * (SIGMA_r|-r)'
  //       stores the results in the column vector sigma_rmr

  void compute_sigmarmr(void);

  // FUNCTION: compute_offset
  // TASK computes the offset and stores it in the matrix 'offset'

  void compute_offset(void);

  // FUNCTION: standardise
  // TASK: standarizes the response

  void standardise(void);

  public:

  // DEFAULT CONSTRUCTOR

  DISTRIBUTION_multgaussian(void) : DISTRIBUTION() {}

   // CONSTRUCTOR

   DISTRIBUTION_multgaussian(const double & a,const datamatrix & b,
                   MCMCoptions * o, const datamatrix & r,
                   const ST::string & fp,const ST::string & fs,
                   const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTRIBUTION_multgaussian(const DISTRIBUTION_multgaussian & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_multgaussian &
   operator=(const DISTRIBUTION_multgaussian & nd);

   // DESTRUCTOR

   ~DISTRIBUTION_multgaussian() {}


  void compute_respminuslinpred(datamatrix & res, const unsigned & co);

  const double & get_scale(const unsigned & r=0,const unsigned & c=0) const
    {
    return sigma_rmr(r,0);
    }

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

