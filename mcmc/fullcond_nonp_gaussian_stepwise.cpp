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
    unsigned i;
    int j;
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


void FULLCOND_nonp_gaussian::reset_effect(const unsigned & pos)
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
vector<double> & lvec, int & number)
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


// BEGIN: MODEL-AVERAGING ------------------------------------------------------

void FULLCOND_nonp_gaussian::save_betas(vector<double> & modell, unsigned & anzahl)
  {
  vector<double> beta_neu;
  unsigned i;
  if(anzahl == -1)
    {
    double * workbeta = beta.getV();
    for(i=0;i<beta.rows();i++,workbeta++)
      beta_neu.push_back(*workbeta);
    }
  else if(anzahl >= 1)
    {
    double fix = fcconst->get_betafix(anzahl);
    beta_neu.push_back(fix);
    }
  // else
  // Vektor "beta_neu" bleibt leer!

  beta_average.push_back(beta_neu);
  }


void FULLCOND_nonp_gaussian::average_posteriormode(vector<double> & crit_weights)
  {
  unsigned i,k;
  int j;
  vector<double> beta_spline;
  for(i=0;i<nrpar;i++)
    beta_spline.push_back(0);
  double beta_fix = 0;

  for(i=0;i<crit_weights.size();i++)
    {
    if(beta_average[i].size()>1)
      {
      for(k=0;k<beta_average[i].size();k++)
        beta_spline[k] += beta_average[i][k] * crit_weights[i];
      }
    else if(beta_average[i].size()==1)
      beta_fix += beta_average[i][0] * crit_weights[i];
    }

  update_linpred(false);
  setbeta(beta_spline.size(),1,0);
  double * workbeta = beta.getV();
  for(i=0;i<beta_spline.size();i++,workbeta++)
    *workbeta = beta_spline[i];
    
  datamatrix pmean_spline = datamatrix(likep->get_nrobs(),1,0);
  workbeta = beta.getV();
  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();
  if(varcoeff)
    {
    int * workindex = index.getV();
    double * workdata = data.getV();
    for(i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        {
        for(j=*itbeg;j<=*itend;j++,workindex++,workdata++)
          effect_sort(pmean_spline,*workbeta*(*workdata),unsigned(*workindex));
        }
      }
    }
  else
    {
    for(i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        effect_sort(pmean_spline,*workbeta,*itbeg,*itend,index);
      }
    }
  datamatrix pmean_fix = datamatrix(likep->get_nrobs(),1,0);
  datamatrix beta_fixx = datamatrix(1,1,beta_fix);
  if(beta_fix != 0)
    pmean_fix.mult(data_forfixed,beta_fixx);     // berechnet den Anteil der fixen Effekte
  pmean_spline.plus(pmean_spline,pmean_fix);     

  // für Ausgabe: Vektor "pmean_spline" muß für Ausgabe sortiert werden!
  workbeta = beta.getV();
  double * workeff = pmean_spline.getV();
  itbeg = posbeg.begin();
  itend = posend.begin();
  int * workindex = index.getV(); 
  for(i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
    {
    if(*itbeg != -1)
      {
      if(!varcoeff)
        *workbeta = pmean_spline(*workindex,0);
      else
        *workbeta = pmean_spline(*workindex,0) / data_forfixed(*workindex,0);
      for(j=*itbeg;j<=*itend;j++)
        workindex++;
      }
    }

  if(!varcoeff)      // Zentrieren
    {
    double intercept = centerbeta();
    fcconst->set_intercept_for_center(intercept);
    datamatrix inter = datamatrix(likep->get_nrobs(),1,-intercept);
    pmean_spline.plus(pmean_spline,inter);         // zentrierte durchschnittliche Fkt, Einträge in Reihenfolge x_1,...,x_n
    }

  likep->add_linearpred_m(pmean_spline,column);     // addiert durchschnittl. Fkt. zum Gesamtprädiktor

  workbeta = beta.getV();
  workeff = betamean.getV();
  for(i=0;i<nrpar;i++,workbeta++,workeff++)
    *workeff = *workbeta * transform;
  }


void FULLCOND_nonp_gaussian::effect_sort(datamatrix & effect, const double & m,
            const unsigned & beg, const unsigned & end,const statmatrix<int> & index)
  {
  unsigned register i;
  int * workindex = index.getV() + beg;
  for (i=beg;i<=end;i++,workindex++)
    effect(*workindex,0)+=m;
  }

void FULLCOND_nonp_gaussian::effect_sort(datamatrix & effect, const double & m, unsigned & row)
  {
  double * workl = effect.getV() + row;
  *workl += m;
  }


 
} // end: namespace MCMC
