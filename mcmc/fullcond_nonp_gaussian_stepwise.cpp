#include "fullcond_nonp_gaussian.h"
#include "fullcond_nonp_gaussian_stepwise.h"

namespace MCMC
{

const datamatrix & FULLCOND_nonp_gaussian::get_data_forfixedeffects(void)
  {

  if ( (data_forfixed.rows() < index.rows()) &&
       (!varcoeff) &&
       ( (type==RW1) || (type==RW2) )
     )
    {
    data_forfixed=datamatrix(index.rows(),1);
    unsigned i,j;
    int * workindex = index.getV();
    double h;
    for(i=0;i<posbeg.size();i++)
      {
      h = effectvdouble[i];
      if (posbeg[i] != -1)
        for(j=posbeg[i];j<=posend[i];j++,workindex++)
          {
          data_forfixed(*workindex,0) = h;
          }
      }

    }
  else if ( (data_forfixed.rows() < data.rows()) && (varcoeff==true) )
    {
    data_forfixed=datamatrix(data.rows(),1);
    unsigned i;
    int * workindex = index.getV();
    double * workdata = data.getV();
    for (i=0;i<data.rows();i++,workindex++,workdata++)
      {
      data_forfixed(*workindex,0) = *workdata;
      }
    }

  return data_forfixed;

  }


double FULLCOND_nonp_gaussian::compute_df(void)
  {

  if ( (lambda_prec != lambda) || (likep->iwlsweights_constant() == false) )
    {

    if (likep->iwlsweights_constant() == false)
      {
      if (varcoeff)
        compute_XWX_varcoeff_env(likep->get_weightiwls(),column);
      else
        compute_XWX_env(likep->get_weightiwls(),column);
      }

    precenv.addtodiag(XXenv,Kenv,1.0,lambda);
    lambda_prec = lambda;
    }


  if (type==MCMC::mrf)
    invprec = envmatdouble(precenv.getXenv(),0,precenv.getDim());
  else
    invprec = envmatdouble(0,nrpar,Kenv.getBandwidth());

  precenv.inverse_envelope(invprec);

  if (varcoeff)
    return invprec.traceOfProduct(XXenv);  
  else
    return invprec.traceOfProduct(XXenv)-1;
  }


ST::string  FULLCOND_nonp_gaussian::get_effect(void)
  {
  ST::string h;

  ST::string t;
  if (type==MCMC::RW1)
    t = "rw1";
  else if (type==MCMC::RW2)
    t = "rw2";
  else if (type==MCMC::seasonal)
    t = "seasonal";
  else if (type==MCMC::mrf)
    t = "spatial";

  if(varcoeff)
    h = datanames[1] + "*" + datanames[0];
  else
    h = datanames[0];

  h = h + "(" + t + ",df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";

  return h;
  }


void FULLCOND_nonp_gaussian::reset_effect(unsigned & pos)
  {

  update_linpred(false);

  unsigned i;
  double * work;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;

  }

  
void FULLCOND_nonp_gaussian::hierarchie_rw1(vector<double> & untervector)
  {

  unsigned number = untervector.size()-1;

  update_stepwise(untervector[0]);
  double df_max = compute_df();

  update_stepwise(untervector[number]);
  double df_min = compute_df();

  if(df_max > 1 && df_min < 1)
     {
     bool geordnet = false;
     unsigned stelle_oben = number;
     unsigned stelle_unten = 0;
     while(geordnet==false)
        {
        unsigned stelle = stelle_oben + stelle_unten;
        update_stepwise(untervector[stelle/2]);
        double df_mitteunten = compute_df();
        update_stepwise(untervector[stelle/2 + 1]);
        double df_mitteoben = compute_df();

        if(df_mitteunten > 1 && df_mitteoben > 1)
          stelle_unten = stelle/2;
        else if(df_mitteunten < 1 && df_mitteoben < 1)
          stelle_oben = stelle/2 + 1;
        else
          {
          geordnet = true;
          vector<double> hilf;
          unsigned i;
          stelle_unten = stelle/2;
          stelle_oben = stelle/2 + 1;
          for(i=0;i<=stelle_unten;i++)
             hilf.push_back(untervector[i]);
          hilf.push_back(-1);
          for(i=stelle_oben;i<untervector.size();i++)
            hilf.push_back(untervector[i]);
          untervector = hilf;
          }
        }
     }
  else if(df_min >= 1)
     {
     untervector.push_back(-1);
     }
  else
     {
     vector<double> hilf;
     hilf.push_back(-1);
     unsigned i;
     for(i=0;i<untervector.size();i++)
        hilf.push_back(untervector[i]);
     untervector = hilf;
     }
  }


void FULLCOND_nonp_gaussian::compute_lambdavec(
vector<double> & lvec,unsigned & number)
  {
  if (get_df_equidist()==true)
     FULLCOND::compute_lambdavec_equi(lvec,number);
  else
     FULLCOND::compute_lambdavec(lvec,number);

  if ( (type==RW1) && (!varcoeff) )
    {
    hierarchie_rw1(lvec);
    }
  else if ( (type==RW1) && (varcoeff) )
    {
    lvec.push_back(-1);
    }
  else if (type==RW2)
    {
    lvec.push_back(-1);
    }
  else if ( (type==mrf) && (varcoeff) )
    {
    lvec.push_back(-1);
    }


  get_forced();
  if(forced_into==false)
     lvec.push_back(0);
  }

 
} // end: namespace MCMC
 