
#include "randomeffect_stepwise.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//----------------- class FULLCOND_random_stepwise -----------------------------
//------------------------------------------------------------------------------

FULLCOND_random_stepwise::FULLCOND_random_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                              FULLCOND_const * fcc,
                              const datamatrix & d, const ST::string & t,
                              const ST::string & fp,const ST::string & pr,
                              const double & la, const unsigned & c)
                            : FULLCOND_random(o,dp,fcc,d,t,fp,pr,la,c)
  {
  identifiable = false;
  grenzfall = 0;
  intercept = 0.0;
  }

// randomslope
FULLCOND_random_stepwise::FULLCOND_random_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                  FULLCOND_const * fcc,
                  const datamatrix & intvar,const datamatrix & effmod,
                  const ST::string & t,
                  const ST::string & fp,const ST::string & pr,
                  const ST::string & prf,
                  const double & la,
                  const bool & inclfixed,const unsigned & c)
                  : FULLCOND_random(o,dp,fcc,intvar,effmod,t,fp,pr,prf,la,false,c)
  {
  grenzfall = 1;

  if(inclfixed == false)
    {
    identifiable = false;
    grenzfall -= 1;
    }

  includefixed = false;
  intercept = 0.0;
  }


FULLCOND_random_stepwise::FULLCOND_random_stepwise(const FULLCOND_random_stepwise & fc)
                            : FULLCOND_random(FULLCOND_random(fc))
  {
  intercept = fc.intercept;
  beta_average = fc.beta_average;
  interactions_pointer = fc.interactions_pointer;
  data_varcoeff_fix = fc.data_varcoeff_fix;
  effmodi = fc.effmodi;
  fbasisp = fc.fbasisp;
  }


const FULLCOND_random_stepwise & FULLCOND_random_stepwise::
         operator=(const FULLCOND_random_stepwise & fc)
  {
  if (this==&fc)
    return *this;

  FULLCOND_random::operator=(FULLCOND_random(fc));

  intercept = fc.intercept;
  beta_average = fc.beta_average;
  interactions_pointer = fc.interactions_pointer;
  data_varcoeff_fix = fc.data_varcoeff_fix;
  effmodi = fc.effmodi;
  fbasisp = fc.fbasisp;

  return *this;
  }


bool FULLCOND_random_stepwise::posteriormode(void)
  {
  unsigned n = nrpar;

  update_linpred(false);

  compute_XWX(likep->get_weightiwls(),column);

  likep->compute_weightiwls_workingresiduals(column);

  unsigned i,j;
  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();
  int * workindex2 = index2.getV();
  itbeg = posbeg.begin();
  itend = posend.begin();
  double * workmuy = muy.getV();
  likep->set_workingresp();

  if (!randomslope)
    {
    for(i=0;i<nrpar;i++,workmuy++,++itbeg,++itend)
      {
      *workmuy = 0;
      for(j=*itbeg;j<=*itend;j++,workindex2++)
        {
        *workmuy+= likep->get_workingres(*workindex2);
        }
      }
    }
  else
    {
    double * datap = data.getV();
    for(i=0;i<n;i++,workmuy++,++itbeg,++itend)
      {
      *workmuy = 0;
      for(j=*itbeg;j<=*itend;j++,workindex2++,datap++)
        {
        *workmuy+= likep->get_workingres(*workindex2)* (*datap);
        }
      }
    }


  itbeg = posbeg.begin();
  itend = posend.begin();
  workmuy = muy.getV();
  double * workbeta = beta.getV();
  double * workXX = XX.getV();

  for(i=0;i<n;i++,workmuy++,++itbeg,++itend,workbeta++,workXX++)
    {
    *workbeta = (*workmuy)/(*workXX+lambda);
    }


  if (randomslope && center)
    {
    double * workbeta = beta.getV();
    double sum=0;
    for (i=0;i<n;i++,workbeta++)
      {
      sum += *workbeta;
      }

    intercept = sum/double(n);

    workbeta = beta.getV();
    for (i=0;i<n;i++,workbeta++)
      *workbeta -= intercept;

    update_fix_effect(intercept);
    }

  update_linpred(true);

  transform = likep->get_trmult(column);

  return FULLCOND::posteriormode();

  }


// BEGIN: For Varying Coefficients ---------------------------------------------

void FULLCOND_random_stepwise::update_fix_effect(double & intercept)
  {
  bool raus = false;
  unsigned j = 1;
  ST::string name_richtig = datanames[1];
  while(j<fcconst->get_datanames().size() && raus==false)
     {
     if(fcconst->get_datanames()[j] == datanames[1])
        {
        raus = true;
        }
     if(fcconst->get_datanames()[j] == (datanames[1]+"_1"))
        {
        raus = true;
        name_richtig = datanames[1] + "_1";
        }
     j = j + 1;
     }
  if(raus == true)
    {
    fcconst->update_fix_effect(j-1,intercept,data_forfixed);
    }
  else
    {
    vector<ST::string> names;
    names.push_back(name_richtig);
    fcconst->include_effect(names,data_forfixed);
    fcconst->update_fix_effect(j,intercept,data_forfixed);
    }
  }


void FULLCOND_random_stepwise::set_pointer_to_interaction(FULLCOND * inter)
  {
  interactions_pointer.push_back(inter);
  }


void FULLCOND_random_stepwise::hierarchical(ST::string & possible)
  {
  unsigned i;
  bool spline = false;
  bool spline1, fix1;
  if(!randomslope)
    {
    for(i=0;i<interactions_pointer.size();i++)
      {
      interactions_pointer[i]->get_inthemodel(spline1,fix1);
      if(spline1 == true)
        spline = true;
      }

    if(spline == true)
      possible = "vfix";
    else
      possible = "alles";
    }
  else
    {
    possible = "alles";
    }
  }


void FULLCOND_random_stepwise::const_varcoeff(void)
  {
  if(randomslope)
    fcconst->posteriormode_const_varcoeff(data_forfixed);
  }

// END: For Varying Coefficients -----------------------------------------------


void FULLCOND_random_stepwise::compute_lambdavec(vector<double> & lvec, int & number)
  {
  if (get_df_equidist()==true)
     FULLCOND::compute_lambdavec_equi(lvec,number);
  else
     FULLCOND::compute_lambdavec(lvec,number);
  if (randomslope && identifiable)
    lvec.push_back(-1);
  get_forced();
  if(forced_into==false)
     lvec.push_back(0);
  }


const datamatrix & FULLCOND_random_stepwise::get_data_forfixedeffects(void)
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
  

double FULLCOND_random_stepwise::compute_df(void)
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

double c = 0;
double w = 0;
for(i=0;i<n;i++,workXX++)
  {
  c += (*workXX * *workXX)/(*workXX+lambda);
  w += *workXX;
  }
c = 1/(w - c);

workXX = XX.getV();
for(i=0;i<n;i++,workXX++)
  {
  df += (*workXX * (lambda + *workXX * (-c * (*workXX + 2*lambda) + 1)))/((*workXX+lambda) * (*workXX+lambda));
  }
df += w*c - 1;

    /*for(i=0;i<n;i++,workXX++)
      {
      df += (*workXX)/(*workXX+lambda);
      } */

    df_lambdaold2 = df_lambdaold1;
    lambdaold2 = lambdaold1;
    df_lambdaold1 = df;
    lambdaold1 = lambda;

    }

  if(identifiable)
    return df;
  else
    //return df-1;
    return df;
  }


void FULLCOND_random_stepwise::update_stepwise(double la)
  {
  lambda = la;
  }


ST::string FULLCOND_random_stepwise::get_effect(void)
  {

  ST::string h;

  if(randomslope)
    h = datanames[1] + "*" + datanames[0];
  else
    h = datanames[0];

  h = h + "(random,df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";

  return h;

  }


void FULLCOND_random_stepwise::reset_effect(const unsigned & pos)
  {

  update_linpred(false);

  unsigned i;
  double * work;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;

  if(randomslope && center)
    update_fix_effect(-intercept);
  intercept = 0.0;
  }


// BEGIN: MODEL-AVERAGING ------------------------------------------------------

void FULLCOND_random_stepwise::save_betas(vector<double> & modell, unsigned & anzahl)
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


void FULLCOND_random_stepwise::average_posteriormode(vector<double> & crit_weights)
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


void FULLCOND_random_stepwise::effect_sort(datamatrix & effect, const double & m,
            const unsigned & beg, const unsigned & end,const statmatrix<int> & index)
  {
  unsigned register i;
  int * workindex = index.getV() + beg;
  for (i=beg;i<=end;i++,workindex++)
    effect(*workindex,0)+=m;
  }

void FULLCOND_random_stepwise::effect_sort(datamatrix & effect, const double & m, unsigned & row)
  {
  double * workl = effect.getV() + row;
  *workl += m;
  }


void FULLCOND_random_stepwise::init_spatialtotal(FULLCOND_nonp_basis * sp,
                                        const ST::string & pnt,
                                        const ST::string & prt)
  {

  fbasisp = sp;
  vector<ST::string> ev = sp->get_effectvalues();

  FULLCOND_random::init_spatialtotal(ev,pnt,prt);

  }


ST::string FULLCOND_random_stepwise::get_befehl(void)
  {
  ST::string h;

  if(randomslope)
    h = datanames[1] + "*" + datanames[0];
  else
    h = datanames[0];

  h = h + "(random,lambda=" + ST::doubletostring(lambda,6) + ")";
  return h;
  }


} // end: namespace MCMC



