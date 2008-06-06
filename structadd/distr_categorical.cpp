

#include "distr_categorical.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_binomial --------------------------
//------------------------------------------------------------------------------


DISTR_binomial::DISTR_binomial(GENERAL_OPTIONS * o, const datamatrix & r,
                               const datamatrix & w)
  : DISTR(o,r,w)

  {

  family = "Binomial";
  changingweight = true;
  updateIWLS = true;
  }


const DISTR_binomial & DISTR_binomial::operator=(
                                      const DISTR_binomial & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  return *this;
  }


DISTR_binomial::DISTR_binomial(const DISTR_binomial & nd)
   : DISTR(DISTR(nd))
  {
  }


void DISTR_binomial::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: logistic distribution function\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTR_binomial::loglikelihood(double * response, double * linpred,
                                     double * weight) const
  {

  if (*linpred >= 10)
    return *weight *(*response * *linpred - *linpred);
  else
    return *weight *(*response * *linpred - log(1+exp(*linpred)));
  }


void DISTR_binomial::compute_mu(const double * linpred,double * mu,
                                bool notransform)
  {
  double el = exp(*linpred);
  *mu = el/(1+el);
  }


void DISTR_binomial::compute_deviance(const double * response,
                   const double * weight,const double * mu,double * deviance,
                   double * deviancesat, double * scale) const
  {

  if (*response==0)
    {
    *deviance = -2* *weight * log(1-*mu);
    *deviancesat = *deviance;
    }
  else if (*response == 1)
    {
    *deviance = -2* *weight*log(*mu);
    *deviancesat = *deviance;
    }
  else
    {
    *deviance = -2* *weight*( *response*log(*mu)+(1-*response)*log(1-*mu) );
    *deviancesat = *deviance +
    2* *weight*( *response*log(*response)+(1-*response)*log(1-*response) );
    }

  }


double DISTR_binomial::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {

  double el = exp(*linpred);
  double mu = el/(1+el);
  if(mu > 0.999)
    mu = 0.999;
  if(mu < 0.001)
    mu = 0.001;
  double v = mu*(1-mu);

  *workingweight = *weight * v;

  *workingresponse = *linpred + (*response - mu)/v;

  if (like)
    {
    if (*linpred >= 10)
      return *weight *(*response * *linpred - *linpred);
    else
      return *weight *(*response * *linpred - log(1+el));
    }
  else
    {
    return 0;
    }

  }


} // end: namespace MCMC



