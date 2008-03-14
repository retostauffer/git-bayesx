
#include "FC_nonp.h"


//------------------------------------------------------------------------------
//----------------- CLASS: FC implementation of member functions ---------------
//------------------------------------------------------------------------------


namespace MCMC
{


FC_nonp::FC_nonp(void)
  {
  }


FC_nonp::FC_nonp(GENERAL_OPTIONS * o,DISTR * lp,
                 const ST::string & t,const ST::string & fp,
                 DESIGN * Dp)
     : FC(o,t,Dp->Zout.rows(),1,fp)
  {
  likep = lp;
  designp = Dp;
  partres = datamatrix(likep->nrobs,1,0);
  param = datamatrix(designp->nrpar,1,0);
  }


FC_nonp::FC_nonp(const FC_nonp & m)
  {
  likep = m.likep;
  partres = m.partres;
  designp = m.designp;
  param = m.param;
  }


const FC_nonp & FC_nonp::operator=(const FC_nonp & m)
  {

  if (this==&m)
	 return *this;
  likep = m.likep;
  partres = m.partres;  
  designp = m.designp;
  param = m.param;
  return *this;

  }


void FC_nonp::update(void)
  {

  FC::update();
  }


bool FC_nonp::posteriormode(void)
  {

  double lambda = 1;
  bool lambdaconst = false;
  datamatrix * linp = likep->linpred_current;


  designp->update_linpred(beta,false);
  partres.minus(likep->response,*linp);


  if ((likep->changingweight) || (changingdesign))
    designp->compute_XtransposedWX_XtransposedWres(partres);
  else
    designp->compute_XtransposedWres(partres);

  if ((likep->changingweight) || (changingdesign) || (!lambdaconst))
    designp->compute_precision(lambda);


  designp->precision.solve(designp->XWres,param);

//  if(center)
//    centerparam();

  designp->compute_f(param,beta);

  designp->update_linpred(beta,true);


  return FC::posteriormode();
  }



void FC_nonp::outresults(void)
  {
  FC::outresults();
  }


void FC_nonp::centerparam(void)
  {

  unsigned i;
  double sum=0;
  double * workparam = param.getV();
  unsigned nrparam = param.rows();

  for (i=0;i<nrparam;i++,workparam++)
    {
    sum+= *workparam;
    }

  workparam = param.getV();

  sum /= double(nrparam);

  for (i=0;i<nrparam;i++,workparam++)
    *workparam-= sum;

  }

  
void FC_nonp::reset(void)
  {

  }


} // end: namespace MCMC



