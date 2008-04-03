
#include "FC_hrandom_variance.h"


//------------------------------------------------------------------------------
//--------- CLASS: FC_hrandom_variance implementation of member functions ------
//------------------------------------------------------------------------------


namespace MCMC
{


FC_hrandom_variance::FC_hrandom_variance(void)
  {

  }


FC_hrandom_variance::FC_hrandom_variance(GENERAL_OPTIONS * o,DISTR * lp,
                  DISTR_gaussian_re * lpRE,
                 const ST::string & t,const ST::string & fp,
                 DESIGN * Dp,FC_nonp * FCn,double la)
     : FC_nonp_variance(o,lp,t,fp,Dp,FCn,la)
  {
  likepRE = lpRE;
  }


FC_hrandom_variance::FC_hrandom_variance(const FC_hrandom_variance & m)
  : FC_nonp_variance(FC_nonp_variance(m))
  {
  likepRE = m.likepRE;
  }


const FC_hrandom_variance & FC_hrandom_variance::operator=(const FC_hrandom_variance & m)
  {

  if (this==&m)
	 return *this;
  FC_nonp_variance::operator=(FC_nonp_variance(m));
  likepRE = m.likepRE;
  return *this;
  }


double FC_hrandom_variance::compute_quadform(void)
  {

  unsigned n;
  double sum = 0;
  double * workbeta = FCnonpp->beta.getV();
  register unsigned i;

  n = FCnonpp->beta.rows();

  datamatrix * linpredRE = likepRE->linpred_current;
  double * linpredREp = (*linpredRE).getV();


  for(i=0;i<n;i++,workbeta++,linpredREp++)
    {
    sum += pow(*workbeta-(*linpredREp),2);
    }

  return sum;

  }


void FC_hrandom_variance::update(void)
  {

  beta(0,0) = rand_invgamma(a_invgamma+0.5*designp->rankK,
                                  b_invgamma+0.5*compute_quadform());

  beta(0,1) = likep->get_scale()/beta(0,0);

  FCnonpp->tau2 = beta(0,0);
  FCnonpp->lambda = beta(0,1);

  transform(0,0) = pow(likep->trmult,2);
  acceptance++;
  FC::update();

  likepRE->sigma2 = beta(0,0);
  }


bool FC_hrandom_variance::posteriormode(void)
  {

  bool conv = FC::posteriormode();
  likepRE->sigma2 = beta(0,0);
  return conv;
  }


} // end: namespace MCMC



