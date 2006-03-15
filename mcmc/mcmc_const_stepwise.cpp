
#include "mcmc_const.h"
#include "mcmc_const_stepwise.h"


namespace MCMC
{

  // CONSTRUCTOR_1 linear effects

FULLCOND_const_stepwise::FULLCOND_const_stepwise(
                      MCMCoptions * o,DISTRIBUTION * dp,const datamatrix & d,
                          const ST::string & t, const int & constant,
                          const ST::string & fs,const ST::string & fr,
                          const unsigned & c)
  : FULLCOND_const(o,dp,d,t,constant,fs,fr,c)
  {

  transform = likep->get_trmult(c);

  changingweight = likep->get_changingweight();

  mu1 = datamatrix(likep->get_nrobs(),1);

  X1 = datamatrix(nrconst,nrconst,0);

  help = datamatrix(nrconst,likep->get_nrobs(),0);

  compute_matrices();

  beta_average.erase(beta_average.begin(),beta_average.end());
  betas_aktuell.erase(betas_aktuell.begin(),betas_aktuell.end());

  intercept_for_center = 0;

  if (X1.rows() < nrconst)
    errors.push_back("ERROR: design matrix for fixed effects is rank deficient\n");


  }

  // COPY CONSTRUCTOR

FULLCOND_const_stepwise::FULLCOND_const_stepwise(
const FULLCOND_const_stepwise & m) : FULLCOND_const(FULLCOND_const(m))
  {

  fcconst = m.fcconst;
  diff_categories = m.diff_categories;
  reference = m.reference;
  coding=m.coding;

  X1 = m.X1;
  mu1 = m.mu1;
  help = m.help;
  changingweight = m.changingweight;
  beta_average = m.beta_average;
  betas_aktuell = m.betas_aktuell;
  intercept_for_center = m.intercept_for_center;
  interactions_pointer = m.interactions_pointer;

  }

  // OVERLOADED ASSIGNMENT OPERATOR

const FULLCOND_const_stepwise & FULLCOND_const_stepwise::
                             operator=(const FULLCOND_const_stepwise & m)
  {
  if (this==&m)
	 return *this;
  FULLCOND_const::operator=(FULLCOND_const(m));

  fcconst = m.fcconst;
  diff_categories = m.diff_categories;
  reference = m.reference;
  coding=m.coding;

  X1 = m.X1;
  mu1 = m.mu1;
  help = m.help;
  changingweight = m.changingweight;
  beta_average = m.beta_average;
  betas_aktuell = m.betas_aktuell;
  intercept_for_center = m.intercept_for_center;
  interactions_pointer = m.interactions_pointer;

  return *this;
  }

void FULLCOND_const_stepwise::compute_matrices(void)
  {
  // computing X1

  unsigned i,j,k;
  double * workweight;
  double * workdata_i_j;
  double * workdata_i_k;

  for (j=0;j<nrconst;j++)
    for(k=0;k<=j;k++)
      {
      X1(j,k) = 0;
      workweight = likep->get_weightp();
      workdata_i_j = data.getV()+j;
      workdata_i_k = data.getV()+k;
      for (i=0;i<likep->get_nrobs();i++,workweight++,workdata_i_j+=nrconst
                         ,workdata_i_k+=nrconst)
        X1(j,k) += *workweight * (*workdata_i_j) * (*workdata_i_k);
      if (k!=j)
        X1(k,j) = X1(j,k);
      }

  X1 = X1.cinverse();
  if (X1.rows() == nrconst)
    X1 = X1.root();
  }


void FULLCOND_const_stepwise::update_intercept(double & m)
  {
  interceptadd+=m;
  beta(interceptpos,0) +=m;
  }

void FULLCOND_const_stepwise::update_interceptold(double & m)
  {
  betameanold(interceptpos,0) -= m;
  }


void FULLCOND_const_stepwise::posteriormode_intercept(double & m)   // wird bei Posteriormode(title,false)
  {
  interceptadd+=m;
  beta(interceptpos,0) +=m;
  betameanold(interceptpos,0) +=m;
  }


void FULLCOND_const_stepwise::update_linold(void)
  {
  if(fabs(interceptadd) >= pow10(-9))
    {
    unsigned i;
    double * worklinold=linold.getV();        // linold = data * beta
    for(i=0;i<linold.rows();i++,worklinold++) // add "interceptadd" to linold
      *worklinold += interceptadd;
    interceptadd = 0;

    double * workbetamean = betamean.getV();
    double * workbeta = beta.getV();
    *workbetamean = *workbeta * transform;
    }
  }


bool FULLCOND_const_stepwise::posteriormode(void)
  {
  unsigned i;
  double * worklinold=linold.getV();        // linold = data * beta
  for(i=0;i<linold.rows();i++,worklinold++) // add interceptadd to linold
    *worklinold += interceptadd;            // interceptadd contains numbers
                                            // from centering other terms
  interceptadd=0;
  likep->fisher(X1,data,column);            // recomputes X1 = (data' W data)^{-1}
  X1.assign((X1.cinverse()));               // continued
  likep->substr_linearpred_m(linold,column);  // substracts linold from linpred
  likep->compute_weightiwls_workingresiduals(column); // computes W(y-linpred)
  beta = X1*data.transposed()*likep->get_workingresiduals();
  linold.mult(data,beta);                   // updates linold
  likep->add_linearpred_m(linold,column);   // updates linpred
  return FULLCOND_const::posteriormode();
  }


void FULLCOND_const_stepwise::posteriormode_single(const vector<ST::string> & names, datamatrix newx,
                                                   const bool include)
  {
  unsigned i;
  if(interceptadd != 0)
    {
    likep->substr_linearpred_m(linold,column);  // substracts linold from linpred
    double * worklinold=linold.getV();        // linold = data * beta
    for(i=0;i<linold.rows();i++,worklinold++) // add interceptadd to linold
      *worklinold += interceptadd;            // interceptadd contains numbers
    interceptadd=0;                           // from centering other terms
    likep->add_linearpred_m(linold,column);
    }

  X2 = datamatrix(names.size()+1,names.size()+1,0);
  datamatrix beta_neu = datamatrix(names.size()+1,1,0);

  datamatrix newx2 = datamatrix(newx.rows(),newx.cols()+1,1);
  double * alt = newx.getV();
  double * neu = newx2.getV();
  unsigned j;
  for(i=0;i<newx.rows();i++)
    {
    neu++;
    for(j=0;j<newx.cols();j++,alt++,neu++)
      *neu = *alt;
    }

  likep->fisher(X2,newx2,column);            // recomputes X1 = (newx' W newx)^{-1}
  X2.assign((X2.cinverse()));               // continued
  likep->compute_weightiwls_workingresiduals(column); // computes W(y-linpred)
  beta_neu = X2*newx2.transposed()*likep->get_workingresiduals();
  likep->substr_linearpred_m(linold,column);  // substracts linold from linpred
  datamatrix linold_single = datamatrix(newx2.rows(),1,0);
  linold_single.mult(newx2,beta_neu);
  linold = linold + linold_single;            // updates linold
  likep->add_linearpred_m(linold,column);     // updates linpred

  if(include)
    include_effect(names,newx);

  //unsigned i;
  double * workbeta = beta.getV();
  double * workbeta_neu = beta_neu.getV();
  double * workbetameanold = betameanold.getV();
  *workbeta = *workbeta + *workbeta_neu;
  *workbetameanold = *workbeta;

  if(include)
    {
    workbeta = beta.getV() + beta.rows()-names.size();
    workbeta_neu = beta_neu.getV() + 1;
    workbetameanold = betameanold.getV() + beta.rows()-names.size();
    for(i=beta.rows()-names.size();i<beta.rows();i++,workbeta_neu++,workbeta++,workbetameanold++)
      {
      *workbeta=*workbeta_neu;
      *workbetameanold = *workbeta;
      }
    }
  else
    {
    bool raus = false;
    i = 1;
    while(i<datanames.size() && raus==false)
      {
      if(datanames[i] == names[0])
        raus = true;
      i = i + 1;
      }
    workbeta = beta.getV() + i-1;
    workbeta_neu = beta_neu.getV() + 1;
    workbetameanold = betameanold.getV() + i-1;
    for(i=0;i<names.size();i++,workbeta_neu++,workbeta++,workbetameanold++)
      {
      *workbeta+=*workbeta_neu;
      *workbetameanold = *workbeta;
      }
    }
  }


void FULLCOND_const_stepwise::safe_const(void)
  {
  double * workbeta = beta.getV();
  const_alt = *workbeta;
  }


void FULLCOND_const_stepwise::set_const_old(void)
  {
  double * workbeta = beta.getV();
  likep->substr_linearpred_m(linold,column);  // substracts linold from linpred
  datamatrix linold_const = datamatrix(linold.rows(),1,*workbeta);
  linold = linold - linold_const;            // updates linold (subtracts the constant)
  *workbeta = const_alt;
  linold_const = datamatrix(linold.rows(),1,*workbeta);
  linold = linold + linold_const;
  likep->add_linearpred_m(linold,column);     // updates linpred
  }


void FULLCOND_const_stepwise::posteriormode_const(void)
  {
  unsigned i;    
  if(interceptadd !=0 )
    {
    double * worklinold=linold.getV();        // linold = data * beta
    for(i=0;i<linold.rows();i++,worklinold++) // add interceptadd to linold
      *worklinold += interceptadd;            // interceptadd contains numbers
    interceptadd=0;                           // from centering other terms
    }
  
  double * workbeta = beta.getV();
  likep->substr_linearpred_m(linold,column);  // substracts linold from linpred
  datamatrix linold_const = datamatrix(linold.rows(),1,*workbeta);
  *workbeta = 0;
  linold = linold - linold_const;            // updates linold (subtracts the constant)
  likep->add_linearpred_m(linold,column);     // updates linpred

  X2 = datamatrix(1,1,1);
  datamatrix beta_const = datamatrix(1,1,0);
  datamatrix newx = datamatrix(linold.rows(),1,1);
  likep->fisher(X2,newx,column);
  X2.assign((X2.cinverse()));
  likep->compute_weightiwls_workingresiduals(column);
  beta_const = X2*newx.transposed()*likep->get_workingresiduals();
  double * workconst = beta_const.getV();
  *workbeta = *workconst;
  likep->substr_linearpred_m(linold,column);  // substracts linold from linpred
  linold_const = datamatrix(linold.rows(),1,*workbeta);
  linold = linold + linold_const;            // updates linold (subtracts the constant)
  likep->add_linearpred_m(linold,column);     // updates linpred
  }


bool FULLCOND_const_stepwise::posteriormode_converged(const unsigned & itnr)
  {
  return likep->posteriormode_converged_fc(beta,beta_mode,itnr);
  }

// For updating the fixed effect with the fixed part of a varying coefficient
void FULLCOND_const_stepwise::update_fix_effect(const unsigned & pos, double & value, datamatrix fix)
  {
  double * work = beta.getV()+pos;
  *work += value;

  double * lin = linold.getV();
  double * dat = fix.getV();
  //linold.mult(data,beta);     // nur wenn "spline" den zentrierten Effekt enthält!
  for(unsigned i=0;i<linold.rows();i++,lin++,dat++)
    {
    *lin += value * *dat;
    }
  }


void FULLCOND_const_stepwise::init_name(const ST::string & na)
  {

  if (fctype == MCMC::factor)
    {
    vector<ST::string> nam;
    unsigned i;
    if(diff_categories.size()==2)
      nam.push_back(na);
    else
      {
      for (i=0;i<diff_categories.size();i++)
        {
        if (diff_categories[i] != reference)
          nam.push_back(na+"_" + ST::doubletostring(diff_categories[i]));
        }
      }
    FULLCOND::init_names(nam);

    char charh ='_';
    ST::string stringh = "\\_";

    ST::string helpname;

    for(i=0;i<nam.size();i++)
      {
      helpname = nam[i].insert_string_char(charh,stringh);
      term_symbolic = term_symbolic + "\\gamma_{"+helpname+"}"+helpname;
      if (i+1<nam.size())
        term_symbolic = term_symbolic + " + ";
      }
    int c = column;
    if(c==0)
      {
      priorassumptions.push_back("Factor $" + na.insert_string_char(charh,stringh) + "$:");
      ST::string liste = "";
      for(i=0;i<nam.size()-1;i++)
        liste = liste + "$" + nam[i].insert_string_char(charh,stringh) + "$, ";
      liste = liste +  + "$" + nam[nam.size()-1].insert_string_char(charh,stringh) + "$";
      priorassumptions.push_back("Resulting variables: " + liste);
      priorassumptions.push_back("diffuse priors");
      priorassumptions.push_back("Coding: " + coding);
      priorassumptions.push_back("\\\\");
      }
    if(c>0)
      {
      priorassumptions.push_back(
      "Factor " + na + " (" + ST::inttostring(c+1) + ". response category):");
      priorassumptions.push_back("diffuse priors");
      priorassumptions.push_back("\\\\");
      }

    }
  else
    {
    FULLCOND_const::init_name(na);
    }
  }


void FULLCOND_const_stepwise::init_names(const vector<ST::string> & na)
  {
  FULLCOND_const::init_names(na);
  }


void FULLCOND_const_stepwise::outresults(void)
  {
  FULLCOND_const::outresults();
  }

void FULLCOND_const_stepwise::outoptions(void)
  {
  FULLCOND_const::outoptions();

  }


ST::string FULLCOND_const_stepwise::get_effect(void)
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


ST::string FULLCOND_const_stepwise::get_befehl(void)
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
    if(datanames.size()>1)
      {
      h = datanames[1];
      for(i=2;i<datanames.size();i++)
        h = h + " + " + datanames[i];
      }
    }
  return h;
  }


void FULLCOND_const_stepwise::update_stepwise(double la)
  {
  lambda=la;
  }


void FULLCOND_const_stepwise::include_effect(const vector<ST::string> & names, const datamatrix & newx)
  {
  if(fctype != factor)
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

    X1 = datamatrix(nrconst,nrconst,0);
    }
  }


void FULLCOND_const_stepwise::reset_effect(const unsigned & pos)
  {
  unsigned i;
  if(interceptadd != 0)
    {
    likep->substr_linearpred_m(linold,column);  // substracts linold from linpred
    double * worklinold=linold.getV();        // linold = data * beta
    for(i=0;i<linold.rows();i++,worklinold++) // add interceptadd to linold
      *worklinold += interceptadd;            // interceptadd contains numbers
    interceptadd=0;                           // from centering other terms
    likep->add_linearpred_m(linold,column);
    }

  if(fctype != factor)
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
    X1 = datamatrix(nrconst,nrconst,0);
    }
  }


void FULLCOND_const_stepwise::set_effect_zero(void)
  {
  /*double * workbeta = beta.getV();
  double * workbetameanold = betameanold.getV();
  unsigned i;
  for(i=0;i<beta.rows();i++,workbeta++,workbetameanold++)
    {
    //if(i==pos)
    //  {
      *workbeta = 0;
      *workbetameanold = 0;
    //  }
    }
  likep->substr_linearpred_m(linold,column);   // zieht den Anteil der fixen Effekte von Gesamtprädiktor ab
  linold.mult(data,beta);                      // berechnet den Anteil der fixen Effekte neu
  likep->add_linearpred_m(linold,column);      // addiert den Anteil der fixen Effekte zum Gesamtprädiktor */

  unsigned i;    
  if(interceptadd !=0 )
    {
    double * worklinold=linold.getV();        // linold = data * beta
    for(i=0;i<linold.rows();i++,worklinold++) // add interceptadd to linold
      *worklinold += interceptadd;            // interceptadd contains numbers
    interceptadd=0;                           // from centering other terms
    }
  
  double * workbeta = beta.getV();
  datamatrix linold_const = datamatrix(linold.rows(),1,*workbeta);
  likep->substr_linearpred_m(linold_const,column);  // substracts linold from linpred
  *workbeta = 0;
  linold = linold - linold_const;            // updates linold (subtracts the constant)

  *workbeta -= likep->get_linearpred().mean(0);

  double mean = likep->get_response().mean(0);
  *workbeta += log(mean);  // speziell für Log-Link, sonst anders!

  linold_const = datamatrix(linold.rows(),1,*workbeta);
  linold = linold + linold_const;            // updates linold (subtracts the constant)
  likep->add_linearpred_m(linold_const,column);     // updates linpred

  }


// BEGIN: MODEL-AVERAGING ------------------------------------------------------

void FULLCOND_const_stepwise::save_betas(vector<double> & modell, unsigned & anzahl)
  {
  betas_aktuell.erase(betas_aktuell.begin(),betas_aktuell.end());
  double * workbeta = beta.getV();
  unsigned i;
  if(fctype != MCMC::factor)
    {
    betas_aktuell.push_back(*workbeta);
    workbeta++;
    for(i=1;i<anzahl;i++)
      {
      if(modell[i-1]!=0)
        {
        betas_aktuell.push_back(*workbeta);
        workbeta++;
        }
      else
        betas_aktuell.push_back(0);
      }
    }
  else     
    {
    if(anzahl >= 1)
      {
      workbeta = fcconst->getbeta().getV() + anzahl;
      betas_aktuell.push_back(*workbeta);
      workbeta++;
      for(i=1;i<data.cols();i++,workbeta++)
        betas_aktuell.push_back(*workbeta);
      }
    else
      {
      for(i=0;i<data.cols();i++)
        betas_aktuell.push_back(0);
      }
    fcconst->beta_to_fix(betas_aktuell);
    }
  }

void FULLCOND_const_stepwise::save_betas2(void)
  {
  beta_average.push_back(betas_aktuell);
  }

void FULLCOND_const_stepwise::beta_to_fix(const vector<double> & betas)
  {
  unsigned i;
  for(i=0;i<betas.size();i++)
    betas_aktuell.push_back(betas[i]);
  }


void FULLCOND_const_stepwise::average_posteriormode(vector<double> & crit_weights)
  {
  if(fctype != MCMC::factor)
    {
    unsigned i;
    unsigned j;
    vector<double> beta_sum;
    for(j=0;j<beta_average[0].size();j++)
      beta_sum.push_back(0);
    for(i=0;i<crit_weights.size();i++)
      {
      for(j=0;j<beta_average[i].size();j++)
        beta_sum[j] += beta_average[i][j] * crit_weights[i];
      }
    setbeta(beta_sum.size(),1,0);
    double * workbeta = beta.getV();
    double * workbetamean = betamean.getV();
    *workbeta = beta_sum[0] + intercept_for_center;
    *workbetamean = *workbeta * transform;
    workbeta++;
    workbetamean++;
    for(i=1;i<beta_sum.size();i++,workbetamean++,workbeta++)
      {
      *workbeta = beta_sum[i];
      *workbetamean = beta_sum[i] * transform;    //(hier fehlt der Term "+addon" aus fullcond.cpp, posteriormode; "addon" wird aber nie !=0)
      }
    likep->substr_linearpred_m(linold,column);   // zieht den Anteil der fixen Effekte von Gesamtprädiktor ab
    linold.mult(data,beta);                      // berechnet den Anteil der fixen Effekte neu
    likep->add_linearpred_m(linold,column);      // addiert den Anteil der fixen Effekte zum Gesamtprädiktor
    }
  }


double FULLCOND_const_stepwise::get_betafix(unsigned & welches)
  {
  double * workbeta = beta.getV() + welches;
  return *workbeta;
  }


void FULLCOND_const_stepwise::set_intercept_for_center(double & dazu)
  {
  intercept_for_center += dazu;
  }

// END: MODEL-AVERAGING --------------------------------------------------------


// For Varying Coefficient Model -----------------------------------------------

void FULLCOND_const_stepwise::posteriormode_const_varcoeff(datamatrix newx)
  {  // Diese Funktion schätzt b0 und b1, damit b0 zum Variierenden Koeffizienten passt! (D.h., es wird nur b0 verwendet)
  unsigned i;
  if(interceptadd != 0)
    {
    likep->substr_linearpred_m(linold,column);  // substracts linold from linpred
    double * worklinold=linold.getV();        // linold = data * beta
    for(i=0;i<linold.rows();i++,worklinold++) // add interceptadd to linold
      *worklinold += interceptadd;            // interceptadd contains numbers
    interceptadd=0;                           // from centering other terms
    likep->add_linearpred_m(linold,column);
    }

  X2 = datamatrix(newx.cols()+1,newx.cols()+1,0);
  datamatrix beta_neu = datamatrix(newx.cols()+1,1,0);

  datamatrix newx2 = datamatrix(newx.rows(),newx.cols()+1,1);
  double * alt = newx.getV();
  double * neu = newx2.getV();
  unsigned j;
  for(i=0;i<newx.rows();i++)
    {
    neu++;
    for(j=0;j<newx.cols();j++,alt++,neu++)
      *neu = *alt;
    }

  likep->fisher(X2,newx2,column);            // recomputes X1 = (newx' W newx)^{-1}
  X2.assign((X2.cinverse()));               // continued
  likep->compute_weightiwls_workingresiduals(column); // computes W(y-linpred)
  beta_neu = X2*newx2.transposed()*likep->get_workingresiduals();
  likep->substr_linearpred_m(linold,column);  // substracts linold from linpred
  datamatrix linold_single = datamatrix(newx2.rows(),1,beta_neu(0,0));
  linold = linold + linold_single;            // updates linold
  likep->add_linearpred_m(linold,column);     // updates linpred

  double * workbeta = beta.getV();
  double * workbeta_neu = beta_neu.getV();
  double * workbetameanold = betameanold.getV();
  *workbeta = *workbeta + *workbeta_neu;
  *workbetameanold = *workbeta;
  }

// End: For Varying Coefficient Model ----------------------------------------------------------


void FULLCOND_const_stepwise::make_design(datamatrix & d)
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

  if(diff_categories.size()>20)
     errors.push_back("ERROR: There are too many different categories!\n");

  bool gefunden = false;
  for(i=0;i<diff_categories.size();i++)
     {
     if(diff_categories[i]==reference)
        gefunden = true;
     }
  if(gefunden==false)
     {
     optionsp->out("WARNING: The value for the reference category does not exist\n");
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


// -----------------------------------------------------------------------------
// ------------------- STEPWISE-FACTOR -----------------------------------------
//------------------------------------------------------------------------------

void FULLCOND_const_stepwise::compute_lambdavec(vector<double> & lvec, int & number)
  {

  assert(fctype == MCMC::factor);
  lvec.push_back(-1);
  get_forced();
  if(forced_into==false)
     lvec.push_back(0);

  }


// Für VCM-Modell  
void FULLCOND_const_stepwise::set_pointer_to_interaction(FULLCOND * inter)
  {
  interactions_pointer.push_back(inter);
  }


void FULLCOND_const_stepwise::hierarchical(ST::string & possible)
  {
  unsigned i;
  bool spline1, fix1;
  bool spline = false;
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


FULLCOND_const_stepwise::FULLCOND_const_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                 FULLCOND_const * fcc,
                 const datamatrix & d,const ST::string & code, int & ref,
                 const ST::string & t,
                 const ST::string & fs,const ST::string & fr,
                 const unsigned & c)
//  : FULLCOND_const(o,dp,d,t,0,fs,fr,c)          //reintun und die Zuweisung unten weglassen! -> ausprobieren!
  {

  fcconst = fcc;
  interactions_pointer.erase(interactions_pointer.begin(),interactions_pointer.end());

  optionsp = o;
  pathresult = fr;
  pathcurrent = fr;
  likep = dp;

  lambda=-1;

  fctype = MCMC::factor;

  reference = ref;

  coding = code;

  make_design(d);

  interceptadd=0;

  sumold = 0;

  datamatrix w = likep->get_weight();

  column = c;

  nrconst = data.cols();

  interceptyes = false;


  results_type="fixed";

  beta_average.erase(beta_average.begin(),beta_average.end());

  /*           braucht man hier nicht!!!
  transform = likep->get_trmult(c);

  changingweight = likep->get_changingweight();

  X1 = datamatrix(nrconst,nrconst,0);

  help = datamatrix(nrconst,likep->get_nrobs(),0);

  compute_matrices();

  if (X1.rows() < nrconst)
    errors.push_back("ERROR: design matrix for fixed effects is rank deficient\n");
  */
  }

//------------------------------------------------------------------------------
//------------------ CLASS: FULLCOND_const_gaussian_special --------------------
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
  vector<double> & lvec, int & number)
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


void FULLCOND_const_gaussian_special::reset_effect(const unsigned & pos)
  {
  double * worklinnew=linnew.getV();
  double * workdatatransformed = datatransformed.getV();
  unsigned i;
  for(i=0;i<data.rows();i++,worklinnew++,workdatatransformed++)
    *worklinnew = -beta(0,0) * *workdatatransformed;

  likep->add_linearpred_m(linnew,column);

  beta(0,0) = 0;
  }

} // end: namespace MCMC







