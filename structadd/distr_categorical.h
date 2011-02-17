
#if !defined (DISTRcategorical_INCLUDED)
#define DISTRcategorical_INCLUDED

#include"../export_type.h"
#include"statmat.h"
#include"Random.h"
#include"GENERAL_OPTIONS.h"
#include"FC.h"
#include"distr.h"


namespace MCMC
{


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_binomial ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_binomial : public DISTR
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

   DISTR_binomial(void) : DISTR()
     {
     }

   // CONSTRUCTOR

   DISTR_binomial(GENERAL_OPTIONS * o, const datamatrix & r,
                  const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_binomial(const DISTR_binomial & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_binomial & operator=(const DISTR_binomial & nd);

   // DESTRUCTOR

   ~DISTR_binomial() {}

  void compute_mu(const double * linpred,double * mu);

  void compute_deviance(const double * response, const double * weight,
                        const double * mu,double * deviance,
                        double * deviancesat, double * scale) const;

  double loglikelihood(double * response, double * linpred,
                       double * weight) const;

  double loglikelihood_weightsone(double * response, double * linpred) const;

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

  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_binomialprobit ----------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_binomialprobit : public DISTR
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

   DISTR_binomialprobit(void) : DISTR()
     {
     }

   // CONSTRUCTOR

   DISTR_binomialprobit(GENERAL_OPTIONS * o, const datamatrix & r,
                  const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_binomialprobit(const DISTR_binomialprobit & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_binomialprobit & operator=(const DISTR_binomialprobit & nd);

   // DESTRUCTOR

   ~DISTR_binomialprobit() {}

  void compute_mu(const double * linpred,double * mu);

  void compute_deviance(const double * response, const double * weight,
                        const double * mu,double * deviance,
                        double * deviancesat, double * scale) const;

  double loglikelihood(double * response, double * linpred,
                       double * weight) const;

  double loglikelihood_weightsone(double * response, double * linpred) const;

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

  void update(void);

  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_poisson -----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_poisson : public DISTR
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

   DISTR_poisson(void) : DISTR()
     {
     }

   // CONSTRUCTOR

   DISTR_poisson(GENERAL_OPTIONS * o, const datamatrix & r,
                  const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_poisson(const DISTR_poisson & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_poisson & operator=(const DISTR_poisson & nd);

   // DESTRUCTOR

   ~DISTR_poisson() {}

  void compute_mu(const double * linpred,double * mu);

  void compute_deviance(const double * response, const double * weight,
                        const double * mu,double * deviance,
                        double * deviancesat, double * scale) const;

  double loglikelihood(double * response, double * linpred,
                       double * weight) const;

  double loglikelihood_weightsone(double * response, double * linpred) const;

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

  };



} // end: namespace MCMC


#endif
