
#include "mcmc_const.h"
#include "mcmc_const_stepwise.h"

namespace MCMC
{

ST::string FULLCOND_const::get_effect(void)
  {
  ST::string h="";
  unsigned i;
  if (fctype==MCMC::factor)
    {
    for(i=0;i<datanames.size();i++)
      h = h + " + " + datanames[i];
    }
  else
    {
    h = datanames[0];
    for(i=1;i<datanames.size();i++)
      h = h + " + " + datanames[i];
    }
  return h;
  }

double FULLCOND_const::compute_df(void)
  {
  if (lambda==-1)
    return data.cols();
  else
    return 0;
  }


void FULLCOND_const::update_stepwise(double la)
  {
  lambda=la;
  }


void FULLCOND_const::include_effect(vector<ST::string> & names, datamatrix & newx)
  {
  unsigned i,j;

  nrconst+=names.size();
  datamatrix dataold = data;

  data = datamatrix(data.rows(),nrconst);

  double * workold = dataold.getV();
  double * workdata = data.getV();
  double * worknew = newx.getV();

  for(i=0;i<dataold.rows();i++)
    {

    for (j=0;j<dataold.cols();j++,workold++,workdata++)
      {
      *workdata = *workold;
      }

    for (j=0;j<newx.cols();j++,workdata++,worknew++)
      *workdata = *worknew;

    }

  for (j=0;j<names.size();j++)
    {
    datanames.push_back(names[j]);
    }

  datamatrix betao = beta;

  setbeta(nrconst,1,0);

  double * workbeta = beta.getV();
  double * workbetao = betao.getV();
  double * workbetameanold = betameanold.getV();
  for(i=0;i<betao.rows();i++,workbetao++,workbeta++,workbetameanold++)
    {
    *workbeta=*workbetao;
    *workbetameanold = *workbeta;
    }
  }


void FULLCOND_const::reset_effect(unsigned & pos)
  {

  unsigned i,j;

  nrconst--;
  datamatrix dataold = data;

  data = datamatrix(data.rows(),nrconst);

  double * workold = dataold.getV();
  double * workdata = data.getV();

  for(i=0;i<dataold.rows();i++)
    {

    for (j=0;j<dataold.cols();j++,workold++)
      {
      if (j!=pos)
        {
        *workdata = *workold;
        workdata++;
        }
      }

    }

  vector<ST::string> dn = datanames;

  datanames.erase(datanames.begin(),datanames.end());
//  datanames.reserve(dn.size()-1);
  for(i=0;i<dn.size();i++)
    {
    if (i!=pos)
      datanames.push_back(dn[i]);
    }


  datamatrix betao = beta;

  setbeta(nrconst,1,0);

  double * workbeta = beta.getV();
  double * workbetao = betao.getV();
  double * workbetameanold = betameanold.getV();
  for(i=0;i<betao.rows();i++,workbetao++)
    {
    if (i!=pos)
      {
      *workbeta=*workbetao;
      *workbetameanold = *workbeta;
      workbeta++;
      workbetameanold++;
      }
    }


  likep->substr_linearpred_m(linold,column);
  linold.mult(data,beta);
  likep->add_linearpred_m(linold,column);
  }


void FULLCOND_const::make_design(datamatrix & d)
  {
  vector<unsigned> zaehlen;

  statmatrix<int> index(d.rows(),1);
  index.indexinit();
  d.indexsort(index,0,d.rows()-1,0,0);

  unsigned i = 0;
  unsigned j;
  while(i<index.rows())
     {
     double anz=0;
     int* p = index.getV() + i;
     int* q = index.getV() + i;
     for(j=i;j<index.rows();j++,p++)
        {
        if (d.get(*p,0) == d.get(*q,0))
                  anz = anz+1;
               }
               zaehlen.push_back(anz);
               diff_categories.push_back(d.get(*q,0));
               i = i + anz;
          }

  bool gefunden = false;
  for(i=0;i<diff_categories.size();i++)
     {
     if(diff_categories[i]==reference)
        gefunden = true;
     }
  if(gefunden==false)
     {
     optionsp->out("WARNING: The value for the reference category doesn't exist!\n");
     optionsp->out("Category " + ST::doubletostring(diff_categories[0]) + " used instead\n");
     reference = diff_categories[0];
     }

  data = datamatrix(d.rows(),diff_categories.size()-1);
  unsigned spalte = 0;
  unsigned zeile = 0;
  unsigned ref;
  for(i=0;i<diff_categories.size();i++)
     {
     if(diff_categories[i]!=reference)
       {
       int* q = index.getV();
       for(j=0;j<zeile;j++,q++)
          data(*q,spalte) = 0;
       q = index.getV() + zeile;
       for(j=zeile;j<zeile+zaehlen[i];j++,q++)
          data(*q,spalte) = 1;
       q = index.getV() + zeile + zaehlen[i];
       for(j=zeile+zaehlen[i];j<d.rows();j++,q++)
          data(*q,spalte) = 0;
       spalte = spalte + 1;
       }

     else if(diff_categories[i]==reference && coding=="effect")
       ref = i;

     zeile = zeile + zaehlen[i];
     }

  if(coding=="effect")
     {
     zeile = 0;
     for(i=0;i<ref;i++)
        zeile = zeile + zaehlen[i];
     int* q = index.getV() + zeile;
     for(i=zeile;i<zeile+zaehlen[ref];i++,q++)
        {
        for(j=0;j<data.cols();j++)
           data(*q,j) = -1;
        }
     }
  }


void FULLCOND_const::compute_lambdavec(vector<double> & lvec,unsigned & number)
  {

  assert(fctype == MCMC::factor);
  lvec.push_back(-1);
  get_forced();
  if(forced_into==false)
     lvec.push_back(0);

  }

// -----------------------------------------------------------------------------
// ------------------- STEPWISE-FACTOR -----------------------------------------
//------------------------------------------------------------------------------

FULLCOND_const::FULLCOND_const(MCMCoptions * o,DISTRIBUTION * dp,
                 const datamatrix & d,const ST::string & code, int
                 & ref,
                 const ST::string & t,
                 const ST::string & fs,const ST::string & fr,
                 const unsigned & c)
  {

  lambda=-1;

  fctype = MCMC::factor;

  reference = ref;

  coding = code;

  make_design(d);

  interceptadd=0;

  likep = dp;

  sumold = 0;

  datamatrix w = likep->get_weight();

  column = c;

  nrconst = data.cols();

  interceptyes = false;

  pathresult = fr;
  pathcurrent = fr;

  results_type="fixed";

//  negbin=false;


  }

FULLCOND_const::FULLCOND_const(MCMCoptions * o,DISTRIBUTION * dp,
                                const datamatrix & d, const ST::string & t,
                                const int & constant, const ST::string & fs,
                                const ST::string & fr,
                                const unsigned & c)
          : FULLCOND(o,d,t,d.cols(),1,fs)
  {

  lambda=-1;

  interceptadd=0;

  fctype = MCMC::fixed;

  likep = dp;

  sumold = 0;

  datamatrix w = likep->get_weight();

  column = c;

  nrconst = data.cols();

  linold = datamatrix(likep->get_nrobs(),1,0);
  linnew = linold;
  linoldp = &linold;
  linnewp = &linnew;

  if (constant>=0)
    {
    identifiable = false;
    interceptpos = constant;
    interceptyes = true;
    }
  else
    {
    interceptyes = false;
    interceptpos = constant;
    }

  pathresult = fr;
  pathcurrent = fr;

  results_type="fixed";

//  negbin=false;

  }


FULLCOND_const::FULLCOND_const(const FULLCOND_const & m) : FULLCOND(FULLCOND(m))
  {
  lambda=m.lambda;
  reference = m.reference;
  coding = m.coding;
  diff_categories = m.diff_categories;
  interceptadd = m.interceptadd;
//  negbin=m.negbin;
  likep = m.likep;
  nrconst = m.nrconst;
  linold = m.linold;
  linnew = m.linnew;
  sumold = m.sumold;
  pathresult = m.pathresult;
  table_results = m.table_results;
  interceptpos = m.interceptpos;
  interceptyes = m.interceptyes;
  }


const FULLCOND_const & FULLCOND_const::operator=(const FULLCOND_const & m)
  {
  if (this==&m)
	 return *this;
  FULLCOND::operator=(FULLCOND(m));
  lambda=m.lambda;
  reference = m.reference;
  diff_categories = m.diff_categories;
  interceptadd = m.interceptadd;
//  negbin=m.negbin;
  likep = m.likep;
  nrconst = m.nrconst;
  linold = m.linold;
  linnew = m.linnew;
  sumold = m.sumold;
  pathresult = m.pathresult;
  table_results = m.table_results;
  interceptpos = m.interceptpos;
  interceptyes = m.interceptyes;
  return *this;
  }

//------------------------------------------------------------------------------
//------------------ CLASS: FULLCOND_const_stepwise_gaussian_special -----------
//------------------------------------------------------------------------------

double FULLCOND_const_gaussian_special::compute_df(void)
  {
  if (lambda==0)
    return 0;
  else
    return 1.0;
  }


FULLCOND_const_gaussian_special::FULLCOND_const_gaussian_special(
                                  MCMCoptions * o,DISTRIBUTION * dp,
                                  const datamatrix & d,const ST::string & t,
                                  const ST::string & fs,const ST::string & fr,
                                  const unsigned & c)
  : FULLCOND_const(o,dp,d,t,false,fs,fr,c)
  {

  fctype = MCMC::nonlinearf;
  lambda = -1;
  transform = likep->get_trmult(c);
  datatransformed=data;
  mu = datamatrix(data.rows(),1);
  }


FULLCOND_const_gaussian_special::FULLCOND_const_gaussian_special(
const FULLCOND_const_gaussian_special & m)  : FULLCOND_const(FULLCOND_const(m))
  {
  datatransformed = m.datatransformed;
  mu = m.mu;
  }


const FULLCOND_const_gaussian_special &
FULLCOND_const_gaussian_special::operator=(
const FULLCOND_const_gaussian_special & m)
  {
  if (this==&m)
	 return *this;
  FULLCOND_const::operator=(FULLCOND_const(m));
  datatransformed = m.datatransformed;
  mu = m.mu;
  return *this;
  }


void FULLCOND_const_gaussian_special::update_stepwise(double la)
  {
  lambda=la;
  unsigned i=0;
  reset_effect(i);
  compute_datatransformed(lambda);
  }


bool FULLCOND_const_gaussian_special::posteriormode(void)
  {
  /*
  double xwx =0;
  unsigned i;
  double * workdatatransformed = datatransformed.getV();
  likep->set_weightp(0);
  double * workweight = likep->get_weightp();
  double * worklinold=linold.getV();
  for(i=0;i<data.rows();i++,workdatatransformed++,workweight++,worklinold++)
    {
    *worklinold = beta(0,0) * *workdatatransformed;
    xwx += *workdatatransformed * *workdatatransformed * *workweight;

    }

  likep->substr_linearpred_m(linold,column);

  mu.minus(likep->get_response(),likep->get_linearpred());

  double sumy = 0;
  likep->set_weightp(0);
  workweight = likep->get_weightp();
  double * workmu = mu.getV();
  workdatatransformed = datatransformed.getV();
  for(i=0;i<data.rows();i++,workdatatransformed++,workweight++,workmu++)
    {
    sumy += *workweight * *workdatatransformed * *workmu;
    }

  beta(0,0) = sumy/xwx;

  double * worklinnew=linnew.getV();
  workdatatransformed = datatransformed.getV();
  for(i=0;i<data.rows();i++,worklinnew++,workdatatransformed++)
    *worklinnew = beta(0,0) * *workdatatransformed;

  likep->add_linearpred_m(linnew,column);
  */
  return FULLCOND_const::posteriormode();
  }


bool FULLCOND_const_gaussian_special::posteriormode_converged(
const unsigned & itnr)
  {
  return likep->posteriormode_converged_fc(beta,beta_mode,itnr);
  }



void FULLCOND_const_gaussian_special::outresults(void)
  {
  FULLCOND_const::outresults();
  }


void FULLCOND_const_gaussian_special::compute_datatransformed(double lambda)
  {
  double * workdata=data.getV();
  double * workdatatransformed = datatransformed.getV();
  unsigned i;
  for (i=0;i<data.rows();i++,workdata++,workdatatransformed++)
    {
    if (lambda == -1)
      *workdatatransformed = *workdata;
    else if (lambda == 1)
      *workdatatransformed = log(*workdata);
    else if ( lambda == 2)
      *workdatatransformed = 1.0/(*workdata+1);
    }

  }


void FULLCOND_const_gaussian_special::compute_lambdavec(
  vector<double> & lvec,unsigned & number)
  {

  lvec.push_back(2);       // 1/x+1
  lvec.push_back(1);       // ln
  lvec.push_back(-1);
  get_forced();
  if(forced_into==false)
     lvec.push_back(0);

  }


ST::string  FULLCOND_const_gaussian_special::get_effect(void)
  {
  ST::string h;
  if (lambda==-1)
    h = datanames[0];
  else if (lambda==1)
    h = "log(" + datanames[0] + ")";
  else if (lambda == 2)
    h = "1/(" + datanames[0] + "+1)";

  return h;
  }


void FULLCOND_const_gaussian_special::reset_effect(unsigned & pos)
  {
  double * worklinnew=linnew.getV();
  double * workdatatransformed = datatransformed.getV();
  unsigned i;
  for(i=0;i<data.rows();i++,worklinnew++,workdatatransformed++)
    *worklinnew = -beta(0,0) * *workdatatransformed;

  likep->add_linearpred_m(linnew,column);

  beta(0,0) = 0;
  }

//------------------------------------------------------------------------------
//------------------ CLASS: FULLCOND_const_gaussian ----------------------------
//------------------------------------------------------------------------------


void FULLCOND_const_gaussian::include_effect(vector<ST::string> & names,
datamatrix & newx)
  {
  FULLCOND_const::include_effect(names,newx);

  X1 = datamatrix(nrconst,nrconst,0);
  X2 = datamatrix(nrconst,likep->get_nrobs());
  }


void FULLCOND_const_gaussian::reset_effect(unsigned & pos)
  {
  FULLCOND_const::reset_effect(pos);

  X1 = datamatrix(nrconst,nrconst,0);
  X2 = datamatrix(nrconst,likep->get_nrobs());
  }

//------------------------------------------------------------------------------
//--- CLASS FULLCOND_const_nongaussian: implementation of member functions -----
//------------------------------------------------------------------------------

void FULLCOND_const_nongaussian::include_effect(vector<ST::string> & names,
datamatrix & newx)
  {
  FULLCOND_const::include_effect(names,newx);

  XWX = datamatrix(nrconst,nrconst);
  }

  
void FULLCOND_const_nongaussian::reset_effect(unsigned & pos)
  {
  FULLCOND_const::reset_effect(pos);

  XWX = datamatrix(nrconst,nrconst);
  }

//------------------------------------------------------------------------------
//--- CLASS FULLCOND_const_nbinomial: implementation of member functions -------
//------------------------------------------------------------------------------


} // end: namespace MCMC







