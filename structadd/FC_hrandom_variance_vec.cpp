
#include "FC_hrandom_variance_vec.h"


//------------------------------------------------------------------------------
//--------- CLASS: FC_hrandom_variance implementation of member functions ------
//------------------------------------------------------------------------------


namespace MCMC
{


void FC_hrandom_variance_vec::read_options(vector<ST::string> & op,
vector<ST::string> & vn)
  {

  int f;

  /*
  1       degree
  2       numberknots
  3       difforder
  4       lambda
  5       a
  6       b
  7       center
  8       map
  9       lambda_re
  10      a_re
  11      b_re
  12      internal_mult
  */


  if (op[12] == "true")
    {
    mult =true;
    f = op[9].strtodouble(lambdastart);
    f = op[10].strtodouble(a_invgamma);
    f = op[11].strtodouble(b_invgamma);
    }
  else
    {
    mult = false;
    }


  }


FC_hrandom_variance_vec::FC_hrandom_variance_vec(void)
  {

  }


FC_hrandom_variance_vec::FC_hrandom_variance_vec(MASTER_OBJ * mp,
                 GENERAL_OPTIONS * o,DISTR * lp,
                  DISTR * lpRE,
                 const ST::string & t,const ST::string & fp,
                 DESIGN * Dp,FC_nonp * FCn,vector<ST::string> & op,
                 vector<ST::string> & vn)
     : FC_nonp_variance_vec(mp,o,lp,t,fp,Dp,FCn,op,vn)
  {
  read_options(op,vn);
  likepRE = lpRE;

  b_invgamma = masterp->level1_likep->trmult*b_invgamma_orig;

//  hyperLambda=rand_gamma(a_invgamma,b_invgamma);
  hyperLambda=0.1;
  }


FC_hrandom_variance_vec::FC_hrandom_variance_vec(const FC_hrandom_variance_vec & m)
  : FC_nonp_variance_vec(FC_nonp_variance_vec(m))
  {
  likepRE = m.likepRE;
  mult = m.mult;
  hyperLambda=m.hyperLambda;
  }


const FC_hrandom_variance_vec & FC_hrandom_variance_vec::operator=(
const FC_hrandom_variance_vec & m)
  {

  if (this==&m)
	 return *this;
  FC_nonp_variance_vec::operator=(FC_nonp_variance_vec(m));
  likepRE = m.likepRE;
  mult = m.mult;
  hyperLambda=m.hyperLambda;
  return *this;
  }



void FC_hrandom_variance_vec::update(void)
  {

  b_invgamma = masterp->level1_likep->trmult*b_invgamma_orig;

  register unsigned i;
  double * workbeta = beta.getV();
  double * workbetafcn = FCnonpp->beta.getV();

  double hyperLambda2 = hyperLambda*hyperLambda;
  double sumtau2 = 0;

  double * linpredREp;
  if (likepRE->linpred_current==1)
    linpredREp = likepRE->linearpred1.getV();
  else
    linpredREp = likepRE->linearpred2.getV();

  double * ww = likepRE->workingweight.getV();

  for (i=0;i<beta.rows();i++,workbeta++,workbetafcn++,ww++)
    {
    *workbeta = rand_inv_gaussian(fabs(hyperLambda)/
    fabs((*workbetafcn - (*linpredREp))),hyperLambda2);
    *ww = 1/(*workbeta);
    sumtau2 += (*workbeta);
    }

  hyperLambda = sqrt(rand_gamma(a_invgamma+beta.rows(),b_invgamma+0.5*sumtau2));

  FCnonpp->tau2 = 1;
  designp->compute_penalty2(beta);

  likepRE->sigma2 = 1;



  acceptance++;
  FC::update();

  }


bool FC_hrandom_variance_vec::posteriormode(void)
  {
  likepRE->wtype=wweightschange_weightsone;
  return  FC_nonp_variance_vec::posteriormode();
  }





} // end: namespace MCMC


