
#include "FC_linear.h"


//------------------------------------------------------------------------------
//----------------- CLASS: FC implementation of member functions ---------------
//------------------------------------------------------------------------------


namespace MCMC
{


void FC_linear::read_options(vector<ST::string> & op,vector<ST::string> & vn)
  {

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
  13      samplemult
  14      constraints
  */

  }


FC_linear::FC_linear(void)
  {
  }


FC_linear::FC_linear(GENERAL_OPTIONS * o,DISTR * lp,datamatrix & d,
                 const ST::string & t,const ST::string & fp,
                 vector<ST::string> & op, vector<ST::string> & vn)
     : FC(o,t,1,1,fp)
  {
  read_options(op,vn);
  likep = lp;
  design = d;
  }


FC_linear::FC_linear(const FC_linear & m)
  : FC(FC(m))
  {
  IWLS = m.IWLS;
  likep = m.likep;
  design = m.design;
  XWX = m.XWX;
  }


const FC_linear & FC_linear::operator=(const FC_linear & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  IWLS = m.IWLS;
  likep = m.likep;
  design = m.design;
  XWX = m.XWX;  
  return *this;
  }


void FC_linear::update_IWLS(void)
  {
  FC::update();
  }

void FC_linear::update(void)
  {
  if (IWLS)
    update_IWLS();
  else
    update_gaussian();
  }

void FC_linear::update_gaussian(void)
  {
  FC::update();
  }

/*
  likep->fisher(XWX,data,column);
  // BEGIN: shrinkage
  if(shrinkage)
    {
    unsigned j;
    for(j=0; j<nrconst; j++)
      {
      XWX(j,j) += 1/variances(j,0);
      }
    }
  // END: shrinkage
  XWX.assign(XWX.cinverse());

  unsigned i;
  double * worklinold=linold.getV();
  for(i=0;i<linold.rows();i++,worklinold++)
    *worklinold += interceptadd;
  interceptadd=0;


  likep->substr_linearpred_m(linold,column);

  likep->compute_weightiwls_workingresiduals(column);

  beta = XWX*data.transposed()*likep->get_workingresiduals();

  linold.mult(data,beta);                   // updates linold
  likep->add_linearpred_m(linold,column);   // updates linpred

  return FULLCOND_const::posteriormode();
*/

/*
void FC_linear::compute_XWX(void)
  {

  unsigned nrconst = beta.rows();
  nrobs = Xt.cols();
  double * XWXp = XWX.getV();
  double * Xt_ip;
  double * Xt_jp;
  double * workingweightp;
  for (i=0;i<nrconst;i++)
    for (j=i;j<nrconst;j++,XWXp++)
      {
      Xt_ip = Xt.getV()+i*nrobs;
      Xt_jp = Xt.getV()+j*nrobs;
      workingweightp->likep->workingweight.getV();

      for (k=0;k<nrobs;k++,Xt_ip++,Xt_jp++,workingweightp++)
        *XWXp = (*workingweightp) * (*Xt_ip)*(*Xt_jp);

      }

  }
 */ 

bool FC_linear::posteriormode(void)
  {
/*
  unsigned i;
  double * worklinold=linold.getV();        // linold = data * beta




  likep->fisher(X1,data,column);            // recomputes X1 = (data' W data)^{-1}


  X1.assign((X1.cinverse()));               // continued
  likep->substr_linearpred_m(linold,column);  // substracts linold from linpred
  likep->compute_weightiwls_workingresiduals(column); // computes W(y-linpred)
  beta = X1*data.transposed()*likep->get_workingresiduals();
  linold.mult(data,beta);                   // updates linold
  likep->add_linearpred_m(linold,column);   // updates linpred
  return FULLCOND_const::posteriormode();


  return FC::posteriormode();
  */
  }


void FC_linear::outoptions(void)
  {
//  optionsp->out("  OPTIONS FOR TERM: " + title + "\n",true);
//  optionsp->out("\n");
  }

void FC_linear::outresults(const ST::string & pathresults)
  {

  }


void FC_linear::reset(void)
  {

  }


} // end: namespace MCMC



