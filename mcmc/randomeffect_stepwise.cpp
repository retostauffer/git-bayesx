
#include "randomeffect.h"
#include "randomeffect_stepwise.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//----------------- class FULLCOND_random --------------------------------------
//------------------------------------------------------------------------------

void FULLCOND_random::compute_lambdavec(vector<double> & lvec, int & number)
  {
  if (get_df_equidist()==true)
     FULLCOND::compute_lambdavec_equi(lvec,number);
  else
     FULLCOND::compute_lambdavec(lvec,number);
  if (randomslope)
    lvec.push_back(-1);
  get_forced();
  if(forced_into==false)
     lvec.push_back(0);
  }


const datamatrix & FULLCOND_random::get_data_forfixedeffects(void)
  {
  // useful for randomslopes only
  if ( (data_forfixed.rows() < data.rows()) && (randomslope==true) )
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
  

double FULLCOND_random::compute_df(void)
  {
  unsigned i;
  double df=0;
  double * workXX=XX.getV();
  unsigned n;

  if (randomslope && includefixed)
    {
    n = nrpar-1;
    df=1;
    }
  else
    n = nrpar;


  if ((lambdaold1==lambda) && (likep->iwlsweights_constant() == true) )
    {
    df = df_lambdaold1;
    }
  else if ((lambdaold2==lambda) && (likep->iwlsweights_constant() == true) )
    {
    df = df_lambdaold2;
    }
  else
    {
    for(i=0;i<n;i++,workXX++)
      {
      df += (*workXX)/(*workXX+lambda);
      }

    df_lambdaold2 = df_lambdaold1;
    lambdaold2 = lambdaold1;
    df_lambdaold1 = df;
    lambdaold1 = lambda;

    }

  return df;
  }


void FULLCOND_random::update_stepwise(double la)
  {
  lambda = la;
  }


ST::string FULLCOND_random::get_effect(void)
  {

  ST::string h;

  if(randomslope)
    h = datanames[1] + "*" + datanames[0];
  else
    h = datanames[0];

  h = h + "(random,df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";

  return h;

  }


void FULLCOND_random::reset_effect(const unsigned & pos)
  {

  update_linpred(false);

  unsigned i;
  double * work;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;

  }


// BEGIN: MODEL-AVERAGING ------------------------------------------------------

void FULLCOND_random::save_betas(vector<double> & modell, unsigned & anzahl)
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


void FULLCOND_random::average_posteriormode(vector<double> & crit_weights)
  {
  unsigned n = nrpar;
  if (includefixed)
    n = nrpar-1;

  unsigned i;
  unsigned j;
  vector<double> beta_spline;
  for(j=0;j<nrpar;j++)
    beta_spline.push_back(0);
  double beta_fix = 0;

  for(i=0;i<crit_weights.size();i++)
    {
    if(beta_average[i].size()>1)
      {
      for(j=0;j<beta_average[i].size();j++)
        beta_spline[j] += beta_average[i][j] * crit_weights[i];
      }
    else if(beta_average[i].size()==1)
      beta_fix += beta_average[i][0] * crit_weights[i];
    }

  update_linpred(false);
  setbeta(nrpar,1,0);
  double * workbeta = beta.getV();
  for(i=0;i<beta_spline.size();i++,workbeta++)
    *workbeta = beta_spline[i];

  datamatrix pmean_spline = datamatrix(likep->get_nrobs(),1,0);
  workbeta = beta.getV();
  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();
  if(randomslope)
    {
    int * workindex = index.getV();
    double * workdata = data.getV();
    for(i=0;i<n;i++,workbeta++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        {
        for(j=*itbeg;j<=*itend;j++,workindex++,workdata++)
          effect_sort(pmean_spline,(*workbeta + beta_spline[nrpar-1])*(*workdata),unsigned(*workindex));
        }
      }
    }
  else
    {
    for(i=0;i<n;i++,workbeta++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        effect_sort(pmean_spline,*workbeta,*itbeg,*itend,index);
      }
    }

  datamatrix pmean_fix = datamatrix(likep->get_nrobs(),1,0);
  datamatrix beta_fixx = datamatrix(1,1,beta_fix);
  if(beta_fix != 0)
    pmean_fix.mult(data_forfixed,beta_fixx);     // berechnet den Anteil der fixen Effekte
  pmean_spline.plus(pmean_spline,pmean_fix);     // durchschnittliche Funktion
  if(includefixed)
    {
    beta_spline[nrpar-1] += beta_fix;
    workbeta = beta.getV() + nrpar-1;
    *workbeta = beta_spline[nrpar-1];
    }

  // für Ausgabe: Vektor "pmean_spline" muß für Ausgabe sortiert werden!
  workbeta = beta.getV();
  double * workeff = pmean_spline.getV();
  itbeg = posbeg.begin();
  itend = posend.begin();
  int * workindex = index.getV(); 
  for(i=0;i<n;i++,workbeta++,++itbeg,++itend)
    {
    if(*itbeg != -1)
      {
      if(!randomslope)
        *workbeta = pmean_spline(*workindex,0);
      else
        *workbeta = (pmean_spline(*workindex,0) / data_forfixed(*workindex,0)) - beta_spline[nrpar-1];
      for(j=*itbeg;j<=*itend;j++)
        workindex++;
      }
    }

  likep->add_linearpred_m(pmean_spline,column);      // addiert die durchschnittl. Funktion zum Gesamtprädiktor

  workbeta = beta.getV();
  workeff = betamean.getV();
  for(i=0;i<nrpar;i++,workbeta++,workeff++)
    *workeff = *workbeta * transform;
  }


void FULLCOND_random::effect_sort(datamatrix & effect, const double & m,
            const unsigned & beg, const unsigned & end,const statmatrix<int> & index)
  {
  unsigned register i;
  int * workindex = index.getV() + beg;
  for (i=beg;i<=end;i++,workindex++)
    effect(*workindex,0)+=m;
  }

void FULLCOND_random::effect_sort(datamatrix & effect, const double & m, unsigned & row)
  {
  double * workl = effect.getV() + row;
  *workl += m;
  }


/*
//------------------------------------------------------------------------------
//------------------- class FULLCOND_random_gaussian ---------------------------
//------------------------------------------------------------------------------

double FULLCOND_random_gaussian::compute_df(void)
  {
  unsigned i;
  double df=0;
  double * workXX=XX.getV();
  unsigned n;

  if (randomslope && includefixed)
    {
    n = nrpar-1;
    df=1;
    }
  else
    n = nrpar;


  if ((lambdaold1==lambda) && (likep->iwlsweights_constant() == true) )
    {
    df = df_lambdaold1;
    }
  else if ((lambdaold2==lambda) && (likep->iwlsweights_constant() == true) )
    {
    df = df_lambdaold2;
    }
  else
    {
    for(i=0;i<n;i++,workXX++)
      {
      df += (*workXX)/(*workXX+lambda);
      }

    df_lambdaold2 = df_lambdaold1;
    lambdaold2 = lambdaold1;
    df_lambdaold1 = df;
    lambdaold1 = lambda;

    }

  return df;
  }


void FULLCOND_random_gaussian::update_stepwise(double la)
  {
  lambda = la;
  }


ST::string FULLCOND_random_gaussian::get_effect(void)
  {

  ST::string h;

  if(randomslope)
    h = datanames[1] + "*" + datanames[0];
  else
    h = datanames[0];

  h = h + "(random,df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";

  return h;

  }


void FULLCOND_random_gaussian::reset_effect(const unsigned & pos)
  {

  update_linpred(false);

  unsigned i;
  double * work;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;

  }

//------------------------------------------------------------------------------
//----------------- class FULLCOND_random_nongaussian --------------------------
//------------------------------------------------------------------------------

double FULLCOND_random_nongaussian::compute_df(void)
  {
  unsigned i;
  double df=0;
  double * workXX=XX.getV();
  unsigned n;

  if (randomslope && includefixed)
    {
    n = nrpar-1;
    df=1;
    }
  else
    n = nrpar;


  if ((lambdaold1==lambda) && (likep->iwlsweights_constant() == true) )
    {
    df = df_lambdaold1;
    }
  else if ((lambdaold2==lambda) && (likep->iwlsweights_constant() == true) )
    {
    df = df_lambdaold2;
    }
  else
    {
    for(i=0;i<n;i++,workXX++)
      {
      df += (*workXX)/(*workXX+lambda);
      }

    df_lambdaold2 = df_lambdaold1;
    lambdaold2 = lambdaold1;
    df_lambdaold1 = df;
    lambdaold1 = lambda;

    }

  return df;
  }


void FULLCOND_random_nongaussian::update_stepwise(double la)
  {
  lambda = la;
  }


ST::string FULLCOND_random_nongaussian::get_effect(void)
  {

  ST::string h;

  if(randomslope)
    h = datanames[1] + "*" + datanames[0];
  else
    h = datanames[0];

  h = h + "(random,df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";

  return h;

  }


void FULLCOND_random_nongaussian::reset_effect(const unsigned & pos)
  {

  update_linpred(false);

  unsigned i;
  double * work;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;

  }

*/
} // end: namespace MCMC



