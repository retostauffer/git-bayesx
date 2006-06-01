
#include "first.h"

#include "gaussian_heteroskedastic.h"

namespace MCMC
{



DISTRIBUTION_gaussianh::DISTRIBUTION_gaussianh(const double & a,
                   const datamatrix & b, MCMCoptions * o, const datamatrix & r,
                   const ST::string & fp,const ST::string & fs,
                   const datamatrix & w)
  : DISTRIBUTION(o,r,w,fp,fs)
  {

  }


DISTRIBUTION_gaussianh::DISTRIBUTION_gaussianh(
const DISTRIBUTION_gaussianh & nd)
: DISTRIBUTION(DISTRIBUTION(nd))
  {
  }


const DISTRIBUTION_gaussianh & DISTRIBUTION_gaussianh::operator=(
const DISTRIBUTION_gaussianh & nd)
  {

  if (this==&nd)
    return *this;
  DISTRIBUTION::operator=(DISTRIBUTION(nd));

  return *this;
  }


void DISTRIBUTION_gaussianh::compute_mu(const double * linpred,double * mu)
                                           const
  {

  }


void DISTRIBUTION_gaussianh::compute_deviance(const double * response,
                             const double * weight,const double * mu,
                             double * deviance,double * deviancesat,
                             const datamatrix & scale,const int & i) const
  {

  }




bool DISTRIBUTION_gaussianh::posteriormode(void)
  {

  return true;

  }


void DISTRIBUTION_gaussianh::update(void)
  {

//  DISTRIBUTION::update();

  }






} // end: namespace MCMC

