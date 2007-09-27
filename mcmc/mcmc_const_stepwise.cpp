
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
  changed_data = true;

  mu1 = datamatrix(likep->get_nrobs(),1);

  X1 = datamatrix(nrconst,nrconst,0);

  help = datamatrix(nrconst,likep->get_nrobs(),0);

//  compute_matrices();

  beta_average.erase(beta_average.begin(),beta_average.end());
  betas_aktuell.erase(betas_aktuell.begin(),betas_aktuell.end());

  intercept_for_center = 0;
  datanames_fixed_only.erase(datanames_fixed_only.begin(),datanames_fixed_only.end());

  if (X1.rows() < nrconst)
    errors.push_back("ERROR: design matrix for fixed effects is rank deficient\n");

  conditional = true;
  utype = "gauss";
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
  X1root = m.X1root;
  X1X = m.X1X;
  mu1 = m.mu1;
  help = m.help;
  changingweight = m.changingweight;
  changed_data = m.changed_data;
  beta_average = m.beta_average;
  betas_aktuell = m.betas_aktuell;
  intercept_for_center = m.intercept_for_center;
  interactions_pointer = m.interactions_pointer;

  data_rest = m.data_rest;
  beta_rest = m.beta_rest;
  linold_rest = m.linold_rest;
  names_rest = m.names_rest;
  datanames_fixed_only = m.datanames_fixed_only;
  fc_df = m.fc_df;
  conditional = m.conditional;

  utype = m.utype;
  proposal = m.proposal;
  weightiwls = m.weightiwls;
  diff = m.diff;
  tildey = m.tildey;
  mode = m.mode;
  linmode = m.linmode;
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
  X1root = m.X1root;
  X1X = m.X1X;
  mu1 = m.mu1;
  help = m.help;
  changingweight = m.changingweight;
  changed_data = m.changed_data;
  beta_average = m.beta_average;
  betas_aktuell = m.betas_aktuell;
  intercept_for_center = m.intercept_for_center;
  interactions_pointer = m.interactions_pointer;
  data_rest = m.data_rest;
  beta_rest = m.beta_rest;
  linold_rest = m.linold_rest;
  names_rest = m.names_rest;
  datanames_fixed_only = m.datanames_fixed_only;
  fc_df = m.fc_df;  
  conditional = m.conditional;

  utype = m.utype;
  proposal = m.proposal;
  weightiwls = m.weightiwls;
  diff = m.diff;
  tildey = m.tildey;
  mode = m.mode;
  linmode = m.linmode;

  return *this;
  }

void FULLCOND_const_stepwise::compute_matrices(void)
  {
  // nur für MCMC

  if(utype == "gauss")
    {
    unsigned i,j;
    double * workweight;
    double * workdata_i_j;
    // computing X1X
    X1X = datamatrix(nrconst,likep->get_nrobs());
    X1root = datamatrix(nrconst,nrconst);
    help = datamatrix(nrconst,likep->get_nrobs(),0);
    double * workhelp = help.getV();
    for (j=0;j<nrconst;j++)
      {
      workweight = likep->get_weightp();
      workdata_i_j = data.getV()+j;
      for(i=0;i<likep->get_nrobs();i++,workweight++,workhelp++,
                      workdata_i_j+=nrconst)
        *workhelp = *workweight * (*workdata_i_j);
      }
    if (X1.rows() == nrconst)
      {
      X1X.mult(X1,help);
      X1root.assign(X1.root());
      }
    }
  else
    {
    unsigned p,k,i;
    double * workw = weightiwls.getV();
    double * workXp;
    double * workXk;
    if(X1X.rows() != nrconst)
      {
      X1X = datamatrix(nrconst,nrconst);
      X1root = datamatrix(nrconst,nrconst);
      }
    double scale = likep->get_scale(column);
    for (p=0;p<nrconst;p++)
      for (k=p;k<nrconst;k++)
        {
        X1X(p,k)=0;
        workw = weightiwls.getV();
        workXp = data.getV()+p;
        workXk = data.getV()+k;
        for(i=0;i<weightiwls.rows();i++,workw++,workXp+=nrconst,workXk+=nrconst)
          X1X(p,k)+= *workw  *  *workXp * *workXk / scale;
        X1X(k,p) = X1X(p,k);
        }
    }
  }


void FULLCOND_const_stepwise::update_intercept(double & m)
  {
  interceptadd+=m;
  beta(interceptpos,0) +=m;
  }

/*void FULLCOND_const_stepwise::update_interceptold(double & m)
  {
  betameanold(interceptpos,0) -= m;
  } */


void FULLCOND_const_stepwise::posteriormode_intercept(double & m)   // wird bei Posteriormode(title,false)
  {
  interceptadd+=m;
  beta(interceptpos,0) +=m;
  betameanold(interceptpos,0) +=m;
  }


void FULLCOND_const_stepwise::update_linold(void)
  {
  if(fabs(interceptadd) >= std::pow(10,-9.0))
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


void FULLCOND_const_stepwise::update_linold_vc(void)
  {
  unsigned i;
  for(i=1;i<nrconst;i++)
    {
    if (effectsadd(i,0) !=0)
      {
      double v = effectsadd(i,0);
      unsigned j;
      double * worklinold=linold.getV();
      double * datap = data.getV()+i;
      for(j=0;j<linold.rows();j++,worklinold++,datap+=nrconst)
        *worklinold += v* (*datap);

      effectsadd(i,0) = 0;
      }
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
  if (calculate_xwx == true || changed_data==true)
    {
    calculate_xwx = false;
    likep->fisher(X1,data,column);            // recomputes X1 = (data' W data)^{-1}
    X1.assign((X1.cinverse()));               // continued
    changed_data = false;
    }
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
  if(fctype != MCMC::factor)
    datanames_fixed_only = datanames;
  else if(fctype == MCMC::factor)
    {
    fcconst->set_datanames_fixed_only(datanames);
    }
  }


void FULLCOND_const_stepwise::init_names(const vector<ST::string> & na)
  {
  FULLCOND_const::init_names(na);
  if(fctype != MCMC::factor)
    datanames_fixed_only = datanames;
  else if(fctype == MCMC::factor)
    {
    fcconst->set_datanames_fixed_only(datanames);
    }
  }


void FULLCOND_const_stepwise::set_datanames_fixed_only(const vector<ST::string> & na)
  {
  if(fctype != MCMC::factor)
    {
    for(unsigned i=0;i<na.size();i++)
      datanames_fixed_only.push_back(na[i]);
    }
  }


void FULLCOND_const_stepwise::outresults(void)
  {
  if(fctype != MCMC::factor)
    {
    if(betamean.rows() == datanames_fixed_only.size())
      {
      datanames = datanames_fixed_only;
      nrpar = datanames.size();
      nrconst = nrpar;
      }
    if (optionsp->get_samplesize() > 0 && flags[0] != 1)
      {
      setbeta(betamean);
      samplestream.close();
      datamatrix sample(optionsp->get_samplesize(),1);
      unsigned i;
      for(i=0;i<nrpar;i++)
        {
        readsample(sample,i);
        betavar(i,0) = sample.var(0);
        } 
      }
    FULLCOND_const::outresults();
    }
  }

void FULLCOND_const_stepwise::outoptions(void)
  {
  if(fctype != MCMC::factor)
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
    nrpar = nrconst;
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
    changed_data = true;
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
    nrpar = nrconst;
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
    changed_data = true;
    }
  }


void FULLCOND_const_stepwise::reset(void)
  {
  linold = datamatrix(linold.rows(),linold.cols(),0);
  interceptadd = 0;
  FULLCOND::reset();
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


double FULLCOND_const_stepwise::compute_df(void)
  {
  double df = 0;
  if(fctype != factor)
    {
    if(data_rest.rows()>1)
      df = df + data_rest.cols();
    df += FULLCOND_const::compute_df();
    }
  return df;
  }


// -----------------------------------------------------------------------------
// --------------------------- Mini-Backfitting --------------------------------
// -----------------------------------------------------------------------------

void FULLCOND_const_stepwise::split_data(const vector<ST::string> & names)
  {
  if(fctype != factor && data.cols() != names.size()+1)
    {
    unsigned i,j;

    nrconst = names.size()+1;  // +1 wegen Intercept
    nrpar = nrconst;
    datamatrix dataold = data;

    data = datamatrix(data.rows(),nrconst);
    datamatrix betao = beta;
    setbeta(nrconst,1,0);

    data_rest = datamatrix(data.rows(),dataold.cols()-nrconst);
    beta_rest = datamatrix(dataold.cols()-nrconst,1,0);
    linold_rest = datamatrix(data.rows(),1,0);

    double * workold = dataold.getV();
    double * workdata = data.getV();
    double * workrest = data_rest.getV();
    double * workbr = beta_rest.getV();
    double * workbeta = beta.getV();
    double * workbetao = betao.getV();
//    double * workbetameanold = betameanold.getV();


    vector<ST::string> dnames;
    vector<unsigned> reihenfolge;

    for(j=0;j<dataold.cols();j++,workbetao++)
      {
      unsigned z = 0;
      bool gefunden = false;
      while(z<names.size() && gefunden == false)
        {
        if(datanames[j] == names[z])
          gefunden = true;
        z++;
        }
      if(gefunden == true || j==0)
        {
        reihenfolge.push_back(1);
        dnames.push_back(datanames[j]);
        *workbeta = *workbetao;
        workbeta++;
        }
      else
        {
        reihenfolge.push_back(0);
        names_rest.push_back(datanames[j]);
        *workbr = *workbetao;
        workbr++;
        }
      }


    for(i=0;i<dataold.rows();i++)
      {
      for (j=0;j<dataold.cols();j++,workold++)
        {
        if(reihenfolge[j] == 1)
          {
          *workdata = *workold;
          workdata++;
          }
        else
          {
          *workrest = *workold;
          workrest++;
          }
        }
      }

    datanames = dnames;
    X1 = datamatrix(nrconst,nrconst,0);
    changed_data = true;    
    linold_rest.mult(data_rest,beta_rest);
    linold -= linold_rest;
    }
  }



void FULLCOND_const_stepwise::merge_data(void)
  {
  if(fctype != factor && names_rest.size()>0)
    {
    unsigned i,j;

    nrconst = names_rest.size() + datanames.size();
    nrpar = nrconst;
    datamatrix dataold = data;

    data = datamatrix(data.rows(),nrconst);
    datamatrix betao = beta;
    setbeta(nrconst,1,0);

    double * workold = dataold.getV();
    double * workdata = data.getV();
    double * workrest = data_rest.getV();
    double * workbr = beta_rest.getV();
    double * workbeta = beta.getV();
    double * workbetao = betao.getV();
//    double * workbetameanold = betameanold.getV();

    for(j=0;j<dataold.cols();j++,workbeta++,workbetao++)
      *workbeta = *workbetao;
    for(j=0;j<data_rest.cols();j++,workbeta++,workbr++)
      *workbeta = *workbr;

    for(i=0;i<data.rows();i++)
      {
      for (j=0;j<dataold.cols();j++,workold++,workdata++)
        *workdata = *workold;

      for (j=0;j<data_rest.cols();j++,workdata++,workrest++)
        *workdata = *workrest;
      }

    for(i=0;i<names_rest.size();i++)
      {
      datanames.push_back(names_rest[i]);
      }

    X1 = datamatrix(nrconst,nrconst,0);
    changed_data = true;    
    linold += linold_rest;
    names_rest.erase(names_rest.begin(),names_rest.end());
    data_rest = datamatrix(1,1,0);
    }
  }


// BEGIN: MODEL-AVERAGING ------------------------------------------------------

void FULLCOND_const_stepwise::save_betas(vector<double> & modell, int & anzahl)
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


/*double FULLCOND_const_stepwise::get_betafix(unsigned & welches)
  {
  double * workbeta = beta.getV() + welches;
  return *workbeta;
  }    */


void FULLCOND_const_stepwise::set_intercept_for_center(double & dazu)
  {
  intercept_for_center += dazu;
  }


void FULLCOND_const_stepwise::update_bootstrap(const bool & uncond)
  {
  interceptadd = 0;
  if(fctype != factor)
    {
    datamatrix betaold = beta;
    if(!uncond)
      {
      nrconst = datanames_fixed_only.size();
      nrpar = nrconst;
      }
    else
      datanames_fixed_only = datanames;

    if(optionsp->get_nriter()<=1)
      {
      ST::string path = samplepath.substr(0,samplepath.length()-4)+"_df.raw";
      fc_df = FULLCOND(optionsp,datamatrix(1,1),"title?",nrconst,1,path);
      fc_df.setflags(MCMC::norelchange | MCMC::nooutput);
      }
      
    if(!uncond)
      {
      if(optionsp->get_nriter()<=1)
        setbeta(nrconst,1,0);
      else
        beta = datamatrix(nrconst,1,0);
      unsigned i,j;
      for(i=0;i<datanames_fixed_only.size();i++)
        {
        bool gefunden = false;
        j = 0;
        while(j<datanames.size() && gefunden==false)
          {
          if(datanames_fixed_only[i] == datanames[j])
            gefunden = true;
          j += 1;
          }
        if(gefunden==true)
          {
          beta(i,0) = betaold(j-1,0);
          fc_df.setbetavalue(i,0,1);
          }
        else
          fc_df.setbetavalue(i,0,0);
        }
      }

    FULLCOND::update_bootstrap();
    fc_df.update_bootstrap_df();
    beta = betaold;
    nrpar = beta.rows();
    nrconst = nrpar;
    }
  }


void FULLCOND_const_stepwise::update_bootstrap_df(void)
  {
  if(fctype != factor)
    {
    conditional = false;  // wird bei MCMCbootstrap aufgerufen, nicht bei MCMCselect
    //effectsadd = datamatrix(nrpar,1,0);

    if(X1root.rows()>1)     // nötig oder nicht?
      {
      X1root = datamatrix(1,1,0);
      X1X = datamatrix(1,1,0);
      }
      
    nrconst = datanames_fixed_only.size();
    nrpar = nrconst;

    if(optionsp->get_nriter()<=1)
      {
      ST::string path = samplepath.substr(0,samplepath.length()-4)+"_df.raw";
      fc_df = FULLCOND(optionsp,datamatrix(1,1),"title?",nrconst,1,path);
      fc_df.setflags(MCMC::norelchange | MCMC::nooutput);
      }

    unsigned i,j;
    for(i=0;i<datanames_fixed_only.size();i++)
      {
      bool gefunden = false;
      j = 0;
      while(j<datanames.size() && gefunden==false)
        {
        if(datanames_fixed_only[i] == datanames[j])
          gefunden = true;
        j += 1;
        }
      if(gefunden==true)
        fc_df.setbetavalue(i,0,1);
      else
        fc_df.setbetavalue(i,0,0);
      }

    fc_df.update_bootstrap_df();
    nrpar = beta.rows();
    nrconst = nrpar;
    }
  }


void FULLCOND_const_stepwise::save_betamean(void)
  {
  interceptadd = 0;
  if(fctype != factor)
    {
    datamatrix betaold = beta;
    nrconst = datanames_fixed_only.size();
    nrpar = nrconst;
    beta = datamatrix(nrconst,1,0);

    unsigned i,j;
    for(i=0;i<datanames_fixed_only.size();i++)
      {
      bool gefunden = false;
      j = 0;
      while(j<datanames.size() && gefunden==false)
        {
        if(datanames_fixed_only[i] == datanames[j])
          gefunden = true;
        j += 1;
        }
      if(gefunden==true)
        beta(i,0) = betaold(j-1,0);
      }

    FULLCOND::save_betamean();
    beta = betaold;
    nrpar = beta.rows();
    nrconst = nrpar;
    }
  }

void FULLCOND_const_stepwise::update_bootstrap_betamean(void)
  {
  if(fctype != MCMC::factor)
    {
    FULLCOND::update_bootstrap_betamean();
    nrpar = betamean.rows();
    nrconst = nrpar;
    beta = datamatrix(nrpar,beta.cols(),1);
    }
  }

void FULLCOND_const_stepwise::outresults_df(unsigned & size)
  {
  if(fctype != MCMC::factor)
    {
    fc_df.update_bootstrap_betamean();
    nrpar = betamean.rows();
    nrconst = nrpar;
    double * workmean;

    ST::string pathdf = pathcurrent.substr(0,pathcurrent.length()-4)+"_df.res";
    ofstream outres(pathdf.strtochar());

    outres << "varname  ";
    outres << "value   ";
    outres << "frequency  ";
    outres << "pmean   " << endl;

    //samplestream.close();
    datamatrix sample(size,1);

    unsigned v,i;
    for(v=1;v<nrconst;v++)
      {
      fc_df.readsample_df(sample,v);
      workmean = fc_df.get_betameanp() + v;
      unsigned anz1 = 0;
      unsigned anz0 = 0;
      double * works = sample.getV();
      for(i=0;i<sample.rows();i++,works++)
        {
        if(*works==0)
          anz0 += 1;
        else if(*works==1)
          anz1 += 1;
        }

      outres << datanames_fixed_only[v] << "   ";
      outres << ST::inttostring(0) << "   " << ST::inttostring(anz0) << "   ";
      if(*workmean == 0)
        outres << "selected"; // ST::doubletostring(*workmean,6);
      else
        outres << "-";
      outres << endl;
      outres << datanames_fixed_only[v] << "   ";
      outres << ST::inttostring(1) << "   " << ST::inttostring(anz1) << "   ";
      if(*workmean == 1)
        outres << "selected"; // ST::doubletostring(*workmean,6);
      else
        outres << "-";
      outres << endl << endl;
      }
    }  // END: if(fctype != MCMC::factor)
  }


void FULLCOND_const_stepwise::update(void)
  {
  if(effectsadd.rows()!=nrpar)
    effectsadd = datamatrix(nrpar,1,0);

  if(utype == "gauss")
    update_gauss();
  else
    update_nongauss();
  }


void FULLCOND_const_stepwise::update_gauss(void)
  {
  if(optionsp->get_nriter()<=1 && conditional)   // nur bei MCMCselect
    effectsadd = datamatrix(nrpar,1,0);

  if(!conditional)
    {
    nrconst = datanames_fixed_only.size();
    nrpar = nrconst;
    datamatrix betaold = beta;
    setbeta(nrconst,1,0);

    unsigned i,j;
    for(i=0;i<datanames_fixed_only.size();i++)
      {
      bool gefunden = false;
      j = 0;
      while(j<datanames.size() && gefunden==false)
        {
        if(datanames_fixed_only[i] == datanames[j])
          gefunden = true;
        j += 1;
        }
      if(gefunden==true)
        {
        beta(i,0) = betaold(j-1,0);
        }
      }
    FULLCOND_const::update();
    nrpar = betaold.rows();
    nrconst = nrpar;
    setbeta(nrconst,1,0);
    beta = betaold;
    }
  else
    {
    if(optionsp->get_nriter()<=1)
      {
      datamatrix betaold = beta;
      setbeta(nrconst,1,0);
      beta = betaold;
      }
    FULLCOND_const::update();
    }
              // FEHLT: Wann neu berechnen???
  if (X1root.rows()==1 || changingweight || optionsp->get_nriter()==1)
    {
    likep->fisher(X1,data,column);            // recomputes X1 = (data' W data)^{-1}
    X1.assign((X1.cinverse()));               // continued
    compute_matrices();
    }

  unsigned i;
  double * worklinold=linold.getV();
  for(i=0;i<linold.rows();i++,worklinold++)
    *worklinold += interceptadd;
  interceptadd=0;

  for(i=1;i<nrconst;i++)
    {
    if (effectsadd(i,0) !=0)
      {
      double v = effectsadd(i,0);
      unsigned j;
      double * worklinold=linold.getV();
      double * datap = data.getV()+i;
      for(j=0;j<linold.rows();j++,worklinold++,datap+=nrconst)
        *worklinold += v* (*datap);

      effectsadd(i,0) = 0;
      }
    }

  likep->substr_linearpred_m(linold,column);

  likep->compute_respminuslinpred(mu1,column);

  beta.mult(X1X,mu1);
  beta+= sqrt(likep->get_scale(column))*X1root*rand_normvek(nrconst);

  linold.mult(data,beta);
  likep->add_linearpred_m(linold,column);
  acceptance++;
  transform = likep->get_trmult(column);
  }

  
void FULLCOND_const_stepwise::update_nongauss(void)
  {
  if(optionsp->get_nriter()==1 && conditional)   // nur bei MCMCselect
    effectsadd = datamatrix(nrpar,1,0);

  if(!conditional)
    {
    nrconst = datanames_fixed_only.size();
    nrpar = nrconst;
    datamatrix betaold = beta;
    datamatrix betamodeold = beta_mode;
    setbeta(nrconst,1,0);

    unsigned i,j;
    for(i=0;i<datanames_fixed_only.size();i++)
      {
      bool gefunden = false;
      j = 0;
      while(j<datanames.size() && gefunden==false)
        {
        if(datanames_fixed_only[i] == datanames[j])
          gefunden = true;
        j += 1;
        }
      if(gefunden==true)
        {
        beta(i,0) = betaold(j-1,0);
        }
      }
    FULLCOND_const::update();
    nrpar = betaold.rows();
    nrconst = nrpar;
    setbeta(nrconst,1,0);
    beta = betaold;
    beta_mode = betamodeold;
    }
  else
    {
    datamatrix betaold = beta;
    if(optionsp->get_nriter()<=1)
      setbeta(nrconst,1,0);
    beta = betaold;
    FULLCOND_const::update();
    }

  double qoldbeta;
  double qnewbeta;

  if (optionsp->get_nriter() == 1)
    {
    diff = linnew;
    weightiwls = datamatrix(likep->get_nrobs(),1,1);
    tildey = weightiwls;
    proposal = beta;
    help = beta;
    mu1=datamatrix(nrconst,1);
    linoldp = &linold;
    linnewp = &linnew;
    linold.mult(data,beta);
    linmode = datamatrix(data.rows(),1);
    //mode = beta;
    }
  if(proposal.rows()!=beta.rows() || help.rows()!=beta.rows())
    {
    proposal = datamatrix(beta.rows(),1);
    help = datamatrix(beta.rows(),1);
    }

  unsigned i;
  if (interceptyes && interceptadd!=0)
    {
    double * work = linoldp->getV();
    for (i=0;i<linoldp->rows();i++,work++)
      *work += interceptadd;
    interceptadd=0;
    }

  for(i=1;i<nrconst;i++)
    {
    if(effectsadd(i,0) !=0)
      {
      double v = effectsadd(i,0);
      unsigned j;
      double * worklinold=linold.getV();
      double * datap = data.getV()+i;
      for(j=0;j<linold.rows();j++,worklinold++,datap+=nrconst)
        *worklinold += v* (*datap);

      effectsadd(i,0) = 0;
      }
    }

  double invscale = 1.0/likep->get_scale(column);
  double logold = likep->loglikelihood();

  linmode.mult(data,beta_mode);  // linmode.mult(data,mode);
  diff.minus(linmode,*linoldp);
  likep->add_linearpred_m(diff,column);

  likep->compute_IWLS_weight_tildey(weightiwls,tildey,column);
  compute_matrices();
  X1root.assign((X1X.cinverse()).root());

  compute_XWtildey(&linmode,invscale);

  beta_mode = X1X.solve(mu1);  // mode = X1X.solve(mu1);

  help.minus(beta,beta_mode);  // help.minus(beta,mode);
  qoldbeta = -0.5*X1X.compute_quadform(help);

  proposal.plus(beta_mode,X1root*rand_normvek(nrconst));  // proposal.plus(mode,X1root*rand_normvek(nrconst));
  help.minus(proposal,beta_mode);  // help.minus(proposal,mode);
  qnewbeta = -0.5*X1X.compute_quadform(help);
  linnewp->mult(data,proposal);
  diff.minus(*linnewp,linmode);
  likep->add_linearpred_m(diff,column);          // (mit proposed)
  double logprop = likep->loglikelihood();     // mit proposed
  double u = log(uniform());

  if (u <= (logprop + qoldbeta - logold - qnewbeta) )
    {
    //datamatrix * mp = linoldp;
    //linoldp = linnewp;
    //linnewp = mp;
    beta.assign(proposal);
    linold.assign(linnew);
    acceptance++;
    }
  else
    {
    diff.minus(*linoldp,*linnewp);
    likep->add_linearpred_m(diff,column);
    }
  }


void FULLCOND_const_stepwise::compute_XWtildey(datamatrix * linb,double & invscale)
  {
  unsigned i,j;
  double * worktildey=tildey.getV();
  double * worklinb = linb->getV();
  double * workw = weightiwls.getV();
  double * workdata = data.getV();
  double h;
  mu1 = datamatrix(nrconst,1,0);

  for(i=0;i<tildey.rows();i++,worktildey++,worklinb++,workw++)
    {
    h = *workw * (*worktildey + *worklinb);
    for (j=0;j<nrconst;j++,workdata++)
      mu1(j,0) += *workdata * h * invscale;
    }
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

              
void FULLCOND_const_stepwise::make_design(const datamatrix & d)
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
  unsigned ref;
  for(i=0;i<diff_categories.size();i++)
     {
     if(diff_categories[i]==reference)
        {
        gefunden = true;
        ref = i;
        }
     }
  if(gefunden==false)
     {
     optionsp->out("WARNING: The value for the reference category does not exist\n");
     optionsp->out("Category " + ST::doubletostring(diff_categories[0]) + " used instead\n");
     reference = diff_categories[0];
     ref = 0;
     }

  data = datamatrix(d.rows(),diff_categories.size()-1);
  unsigned spalte = 0;
  unsigned zeile = 0;
  for(i=0;i<diff_categories.size();i++)
     {
     if(diff_categories[i]!=reference)
       {
       int* q = index.getV();
       for(j=0;j<zeile;j++,q++)
          data(*q,spalte) = 0;
       q = index.getV() + zeile;
       double wert = 1;
       if(coding=="userdef")
         //wert = double(d.rows())/double(zaehlen[i]);
         //wert = double(zaehlen[i])/double(d.rows());
         wert = double(zaehlen[ref])/double(d.rows());
       for(j=zeile;j<zeile+zaehlen[i];j++,q++)
          data(*q,spalte) = wert;
       q = index.getV() + zeile + zaehlen[i];
       for(j=zeile+zaehlen[i];j<d.rows();j++,q++)
          data(*q,spalte) = 0;
       spalte = spalte + 1;
       }

     //else if(diff_categories[i]==reference && (coding=="effect" || coding=="userdef"))
      //ref = i;

     zeile = zeile + zaehlen[i];
     }

  if(coding=="effect" || coding=="userdef")
     {
     zeile = 0;
     double wert = -1;
     vector<double> werte;
     if(coding=="userdef")
       {
       for(i=0;i<diff_categories.size();i++)
         {
         if(diff_categories[i]!=reference)
           werte.push_back(-1* double(zaehlen[i])/double(d.rows()));
         }
       }
     //if(coding=="userdef")
       //wert = wert* double(d.rows())/double(zaehlen[ref]);
       //wert = wert* double(zaehlen[ref])/double(d.rows());
     for(i=0;i<ref;i++)
        zeile = zeile + zaehlen[i];
     int* q = index.getV() + zeile;
     for(i=zeile;i<zeile+zaehlen[ref];i++,q++)
        {
        for(j=0;j<data.cols();j++)
          {
          if(coding=="userdef")
            wert = werte[j];
          data(*q,j) = wert;
          }
        }
     }
  }


// ------------ Funktionen für gesamte Hatmatrix: ------------------------------


unsigned FULLCOND_const_stepwise::get_nrvar(void)
  {
  return data.cols();
  }


void FULLCOND_const_stepwise::createreml(datamatrix & X,datamatrix & Z,
                                const unsigned & Xpos, const unsigned & Zpos)
  {
  unsigned i,j;

  double * workdata= data.getV();
  double * workX = X.getV()+Xpos;
  unsigned s = X.cols()-data.cols();

  for (i=0;i<data.rows();i++,workX+=s)
    for(j=0;j<data.cols();j++,workdata++,workX++)
      {
      *workX = *workdata;
      }

  }


// -----------------------------------------------------------------------------
// ------------------- STEPWISE-FACTOR -----------------------------------------
//------------------------------------------------------------------------------

void FULLCOND_const_stepwise::compute_lambdavec(vector<double> & lvec, int & number)
  {

  assert(fctype == MCMC::factor);
  lvec.push_back(-1);
  if(forced_into==false)
     lvec.push_back(0);

  // Startwert für lambda aus df:
  if(spfromdf==true && fctype == MCMC::factor)
    {
    if(dfstart!=0)
      lambdastart = -1;
    else if(dfstart==0)
      lambdastart = 0;
    }
  }


// Für VCM-Modell  
void FULLCOND_const_stepwise::set_pointer_to_interaction(FULLCOND * inter)
  {
  interactions_pointer.push_back(inter);
  }


void FULLCOND_const_stepwise::get_interactionspointer(vector<FULLCOND*> & inter)
  {
  inter = interactions_pointer;
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
                 FULLCOND_const_stepwise * fcc,
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







