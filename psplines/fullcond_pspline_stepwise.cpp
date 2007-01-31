
#include "first.h"

#include "fullcond_pspline_stepwise.h"


namespace MCMC
{

  // CONSTRUCTOR 1  (for additive models)

FULLCOND_pspline_stepwise::FULLCOND_pspline_stepwise(MCMCoptions * o,
                      DISTRIBUTION * dp, FULLCOND_const * fcc, const datamatrix & d,
                      const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                      const fieldtype & ft, const ST::string & monotone, const ST::string & ti,
                      const ST::string & fp, const ST::string & pres, const bool & deriv,
                      const double & l, const int & gs, const bool & diag, const unsigned & c)
  : FULLCOND_pspline_gaussian(o,dp,fcc,d,nrk,degr,kp,ft,monotone,ti,fp,pres,deriv,l,gs,diag,0.0,0.0,c)
  {

  beta_average.erase(beta_average.begin(),beta_average.end());
  lambda_nr = -1;
  //lambdas_local = datamatrix(rankK,1,1);

  if (type==RW1)
    grenzfall = 0;
  else if (type == RW2)
    grenzfall = 1;

  dimX = nrpar - rankK - 1;
  dimZ = rankK;
  }


// CONSTRUCTOR 2  (for varying coefficients term)

FULLCOND_pspline_stepwise::FULLCOND_pspline_stepwise(MCMCoptions * o, DISTRIBUTION * dp,
                      FULLCOND_const * fcc,const datamatrix & effmod, const datamatrix & intact,
                      const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                      const fieldtype & ft, const ST::string & monotone, const ST::string & ti,
                      const ST::string & fp, const ST::string & pres, const bool & deriv,
                      const double & l, const int & gs, const bool & nofixed, const unsigned & c)
  : FULLCOND_pspline_gaussian(o,dp,fcc,effmod,intact,nrk,degr,kp,ft,monotone,ti,fp,pres,false,l,gs,c)
  {

  beta_average.erase(beta_average.begin(),beta_average.end());
  lambda_nr = 0;

  data_forfixed = intact;
  effmodi = effmod;

  unsigned nrobs = index.rows();

  if (data_varcoeff_fix.rows() < nrobs)
    {
    data_varcoeff_fix = datamatrix(nrobs,2,1);
    int * workindex = index.getV();
    vector<int>::iterator freqwork = freqoutput.begin();
    for(unsigned i=0;i<nrobs;i++,workindex++,freqwork++)
      {
      data_varcoeff_fix(i,0) = intact(i,0);
      data_varcoeff_fix(i,1) = effmod(i,0)*intact(i,0);
      }
    }

  if (type==RW1)
    grenzfall = 1;
  else if (type == RW2)
    grenzfall = 2;

  //VCM_neu
  if(nofixed == true)
    identifiable = false;
  }

  // COPY CONSTRUCTOR

FULLCOND_pspline_stepwise::FULLCOND_pspline_stepwise(const FULLCOND_pspline_stepwise & fc)
  : FULLCOND_pspline_gaussian(FULLCOND_pspline_gaussian(fc))
  {
  beta_average = fc.beta_average;
  lambda_nr = fc.lambda_nr;
  lambdas_local = fc.lambdas_local;
  data_varcoeff_fix = fc.data_varcoeff_fix;
  effmodi = fc.effmodi;
  df_lambdaold = fc.df_lambdaold;
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const FULLCOND_pspline_stepwise & FULLCOND_pspline_stepwise::operator=(
                                            const FULLCOND_pspline_stepwise & fc)
  {
  if (this == &fc)
    return *this;
  FULLCOND_pspline_gaussian::operator=(FULLCOND_pspline_gaussian(fc));

  beta_average = fc.beta_average;
  lambda_nr = fc.lambda_nr;
  lambdas_local = fc.lambdas_local;
  data_varcoeff_fix = fc.data_varcoeff_fix;
  effmodi = fc.effmodi;
  df_lambdaold = fc.df_lambdaold;

  return *this;
  }


bool FULLCOND_pspline_stepwise::posteriormode(void)
  {

  unsigned i;
  transform = likep->get_trmult(column);
  fchelp.set_transform(transform);

  likep->substr_linearpred_m(spline,column,true);

  if(varcoeff && lambda == -2)
    {
    datamatrix X = datamatrix(2,2,0);
    datamatrix betas = datamatrix(2,1,0);

    likep->fisher(X,data_varcoeff_fix,column);            // recomputes X1 = (newx' W newx)^{-1}
    X.assign((X.cinverse()));               // continued
    likep->compute_weightiwls_workingresiduals(column); // computes W(y-linpred)
    betas = X*data_varcoeff_fix.transposed()*likep->get_workingresiduals();
    spline.mult(data_varcoeff_fix,betas);
    likep->substr_linearpred_m(-spline,column,true);

    if(center)
      {
      intercept = betas(0,0) + 0.5*betas(1,0)*(effmodi.max(0)+effmodi.min(0));
      int * workindex = index.getV();
      for(i=0;i<spline.rows();i++,workindex++)
        spline(*workindex,0) -= intercept*data_forfixed(*workindex,0);
      betas(0,0) -= intercept;
      update_fix_effect();
      intercept = 0.0;
      }
      
    double * fchelpbetap = fchelp.getbetapointer();
    vector<int>::iterator freqwork = freqoutput.begin();
    int * workindex = index.getV();
    for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
      {
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        {
        *fchelpbetap = betas(0,0) + effmodi(*workindex,0)*betas(1,0);
        fchelpbetap++;
        }
      }

    return fchelp.posteriormode();
    //return FULLCOND_nonp_basis::posteriormode();
    }
  else
    {
    if ( (lambda_prec != lambda) || (likep->iwlsweights_constant() == false) )
      {
      if (likep->iwlsweights_constant() == false)
        {
        compute_XWXenv(likep->get_weightiwls(),column);
        }
      prec_env.addto(XX_env,Kenv,1.0,lambda);
      lambda_prec = lambda;
      }
    likep->compute_workingresiduals(column);
    compute_XWtildey(likep->get_weightiwls(),likep->get_workingresiduals(),1.0,column);

    prec_env.solve(muy,beta);

  // monotone Regression!
    if(decreasing)
      {
      bool ok = false;
      while(!ok)
        {
        bool ok2 = true;
        for(unsigned i=1;i<nrpar;i++)
          {
          double diff = beta(i,0)-beta(i-1,0);
          if(diff > 0.0001)
            {
            ok2 = false;
            double mean = 0.5*( beta(i-1,0)+beta(i,0) );
            beta(i-1,0) = mean;
            beta(i,0) = mean;
            }
          if(diff > 0)
            {
            double help = beta(i,0);
            beta(i,0) = beta(i-1,0);
            beta(i-1,0) = help;
            }
          }
        ok = ok2;
        }
      beta.sortcol(0,nrpar-1,0);
      datamatrix bsort = beta;
      for(unsigned j=0;j<nrpar;j++)
        beta(j,0) = bsort(nrpar-1-j,0);
      }

    if(increasing)
      {
      bool ok = false;
      while(!ok)
        {
        bool ok2 = true;
        for(unsigned i=1;i<nrpar;i++)
          {
          double diff = beta(i-1,0)-beta(i,0);
          if(diff > 0.0001)
            {
            ok2 = false;
            double mean = 0.5*(beta(i-1,0)+beta(i,0));
            beta(i-1,0) = mean;
            beta(i,0) = mean;
            }
          if(diff > 0)
            {
            double help = beta(i,0);
            beta(i,0) = beta(i-1,0);
            beta(i-1,0) = help;
            }
          }
        ok = ok2;
        }
      beta.sortcol(0,nrpar-1,0);
      }
  // ENDE: monoton

/*center = false;
datamatrix eins = datamatrix(likep->get_nrobs(),1,1);
datamatrix hilf = datamatrix(nrpar,1,0);
compute_XWtildey(eins,eins,1.0,column);
// prec_env.solve(betaweight,hilf);
prec_env.solve(muy,hilf);
double sum1 = 0;
double sum2 = 0;
double * s1 = hilf.getV();
double * s2 = beta.getV();
//double * weight = betaweight.getV();
double * weight = muy.getV();
for(unsigned i=0;i<nrpar;i++,s1++,s2++,weight++)
  {
  sum1 += *weight * *s1;
  sum2 += *weight * *s2;
  }
double * work = beta.getV();
s1 = hilf.getV();
for(unsigned i=0;i<nrpar;i++,work++,s1++)
  *work -= *s1*sum2/sum1;
*/

    add_linearpred_multBS();

    bool interaction2 = false;
    if(interactions_pointer.size()>0)
      interaction2 = search_for_interaction();

    if(center)
      {
      compute_intercept();
      if(!varcoeff)
        {
        if(interaction==false || interaction2 == false)
          fcconst->posteriormode_intercept(intercept);
        else
          fcconst->update_intercept(intercept);
        }
      else
        {
        update_fix_effect();
        }
      }

    if(interaction == false || interaction2 == false)
      {
      if(center)
        {
        if(!varcoeff)
          {
          int * workindex = index.getV();
          for(i=0;i<spline.rows();i++,workindex++)
            spline(*workindex,0) -= intercept;
          }
        else
          {
          int * workindex = index.getV();
          for(i=0;i<spline.rows();i++,workindex++)
            spline(*workindex,0) -= intercept*data_forfixed(*workindex,0);
          }
        }
      double * fchelpbetap = fchelp.getbetapointer();

      if(gridsize < 0)
        {
        if(varcoeff)
          {
          multBS(splinehelp,beta);
          if(center)
            {
            int * workindex = index.getV();
            for(i=0;i<splinehelp.rows();i++,workindex++)
              splinehelp(i,0) -= intercept;
            }
          }

        vector<int>::iterator freqwork = freqoutput.begin();
        int * workindex = index.getV();
        for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
          {
          if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
            {
            if(varcoeff)
              *fchelpbetap = splinehelp(i,0);
            else
              *fchelpbetap = spline(*workindex,0);
            fchelpbetap++;
            }
          }
        }
      else
        {
        multDG(splinehelp,beta);
        for(i=0;i<gridsize;i++,fchelpbetap++)
          *fchelpbetap = splinehelp(i,0) - intercept;
        }

      intercept = 0.0;

      write_derivative();

      if(derivative)
        fcderivative.posteriormode();

      fchelp.posteriormode();
      return FULLCOND_nonp_basis::posteriormode();
      }  // end: if(interaction == false)
    else
      {
      return true;
      }
    } // END: else if(lambda != 2)
  }


void FULLCOND_pspline_stepwise::create_weight(datamatrix & w)
  {
  vector<int>::iterator freqwork = freqoutput.begin();
  int * workindex = index.getV();
  unsigned i;
  for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
    {
    if(freqwork==freqoutput.begin() || (*freqwork!=*(freqwork-1) && *freqwork == *(freqoutput.end()-1) ))
      w(*workindex,0) = 1;
    }
  }


// BEGIN: For Varying Coefficients ---------------------------------------------

void FULLCOND_pspline_stepwise::update_fix_effect(void)
  {
  bool raus = false;
  unsigned j = 1;
  ST::string name_richtig = datanames[0];
  while(j<fcconst->get_datanames().size() && raus==false)
     {
     if(fcconst->get_datanames()[j] == datanames[0])
        {
        raus = true;
        }
     if(fcconst->get_datanames()[j] == (datanames[0]+"_1"))
        {
        raus = true;
        name_richtig = datanames[0] + "_1";
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
    interactions_pointer[1]->set_inthemodel(-1);
    fcconst->update_fix_effect(j,intercept,data_forfixed);
    }
  }


void FULLCOND_pspline_stepwise::const_varcoeff(void)
  {
  if(varcoeff)
    fcconst->posteriormode_const_varcoeff(data_forfixed);
  }

// END: For Varying Coefficients -----------------------------------------------

// BEGIN: F�R INTERAKTIONEN-----------------------------------------------------

bool FULLCOND_pspline_stepwise::changeposterior(const datamatrix & main,const double & inter)
  {
  unsigned i;

  vector<int>::iterator freqwork = freq.begin();
  int * workindex = index.getV();

// spline �ndern
  for(i=0;i<spline.rows();i++,freqwork++,workindex++)
    spline(*workindex,0) += main(*freqwork,0) - intercept - inter;

  intercept = 0.0;

// fchelp �ndern
  double * fchelpbetap = fchelp.getbetapointer();
  freqwork = freq.begin();
  workindex = index.getV();
  for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
    {
    if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
      {
      *fchelpbetap = spline(*workindex,0);
      fchelpbetap++;
      }
    }

  write_derivative();
// posteriormode
  if(derivative)
    fcderivative.posteriormode();

  fchelp.posteriormode();
  return FULLCOND_nonp_basis::posteriormode();

  }


bool FULLCOND_pspline_stepwise::changeposterior2(const datamatrix & main, const double & inter)
  {

  unsigned i;

  vector<int>::iterator freqwork = freq.begin();
  int * workindex = index.getV();

  datamatrix hilf = datamatrix(spline.rows(),1,0);

// spline �ndern
  for(i=0;i<spline.rows();i++,freqwork++,workindex++)
    {
    spline(*workindex,0) += main(*freqwork,0) + inter;  
    hilf(*workindex,0) += main(*freqwork,0) + inter;
    }

  likep->add_linearpred(hilf);

// Intercept �ndern
  intercept = 0.0;

  return FULLCOND_nonp_basis::posteriormode();
  }


void FULLCOND_pspline_stepwise::set_pointer_to_interaction(FULLCOND * inter)
  {
  interactions_pointer.push_back(inter);
  }

void FULLCOND_pspline_stepwise::get_interactionspointer(vector<FULLCOND*> & inter)
  {
  inter = interactions_pointer;
  }


bool FULLCOND_pspline_stepwise::search_for_interaction(void)
  {
  unsigned i;
  bool thereis = false;
  for(i=0;i<interactions_pointer.size();i++)
    {
    bool drin, fix;
    interactions_pointer[i]->get_inthemodel(drin,fix);
    if(drin == true)
      thereis = true;
    }
/*  if(thereis == true)
    interaction = true;
  else
    interaction = false;   */
  return thereis;
  }


void FULLCOND_pspline_stepwise::wiederholen(FULLCOND * haupt, bool konst)
  {
  unsigned i;
  for(i=0;i<interactions_pointer.size();i++)
    {
    bool drin, fix;
    interactions_pointer[i]->get_inthemodel(drin,fix);
    if(drin == true)
      interactions_pointer[i]->get_zentrierung(haupt,konst);
    }
  }


void FULLCOND_pspline_stepwise::wiederholen_fix(FULLCOND * haupt, int vorzeichen, bool inter)
  {
  unsigned i;
  for(i=0;i<interactions_pointer.size();i++)
    {
    bool drin, fix;
    interactions_pointer[i]->get_inthemodel(drin,fix);
    if(drin == true)
      interactions_pointer[i]->set_zentrierung(haupt,vorzeichen,inter);
    }
  }


void FULLCOND_pspline_stepwise::hierarchical(ST::string & possible)
  {
  unsigned i;
  bool spline = false;
  bool fix = false;
  bool spline1, fix1;
  if(!varcoeff)
    {
    for(i=0;i<interactions_pointer.size();i++)
      {
      interactions_pointer[i]->get_inthemodel(spline1,fix1);
      if(spline1 == true)
        spline = true;
      if(fix1 == true)
        fix = true;
      }

    if(interaction)
      {
      if(spline == true)
        possible = "spline";
      else if(fix == true && spline == false)
        possible = "spfix";
      else
        possible = "alles";
      }
    else     // VC
      {
      if(spline == true && number == -1)
        possible = "vfix";
      else
        possible = "alles";
      }
    }
  else
    {
    possible = "valles";
    }
  }


// BEGIN: MODEL-AVERAGING ------------------------------------------------------

void FULLCOND_pspline_stepwise::save_betas(vector<double> & modell, int & anzahl)
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


void FULLCOND_pspline_stepwise::average_posteriormode(vector<double> & crit_weights)
  {
  unsigned i;
  unsigned j;
  vector<double> beta_spline;
  for(j=0;j<nrpar;j++)
    beta_spline.push_back(0);
  double beta_fix = 0;
  double alpha_fix = 0;
  for(i=0;i<crit_weights.size();i++)
    {
    if(beta_average[i].size()>1)
      {
      for(j=0;j<beta_average[i].size();j++)
        beta_spline[j] += beta_average[i][j] * crit_weights[i];
      }
    else if(beta_average[i].size()==1)
      {
      beta_fix += beta_average[i][0] * crit_weights[i];
      if(!varcoeff)          // Gerade so zentrieren, da� Integral=0, d.h. G((x_max-x_min)/2) = 0!
        alpha_fix += - beta_average[i][0] * crit_weights[i] * (data_forfixed.max(0)+data_forfixed.min(0))/2;
      }
    }
  double helpdouble = -alpha_fix;
  fcconst->set_intercept_for_center(helpdouble);

  setbeta(beta_spline.size(),1,0);
  double * workbeta = beta.getV();
  for(i=0;i<beta_spline.size();i++,workbeta++)
    *workbeta = beta_spline[i];
  datamatrix pmean_spline = datamatrix(likep->get_nrobs(),1,0);
  if(gridsize < 0)
    multBS_sort(pmean_spline,beta);
  else    // Vorsicht: hier passt es wahrscheinlich �berhaupt nicht mehr, weil nicht alle x_1,...,x_n verwendet werden!!!
    multDG(pmean_spline,beta);        // FEHLT NOCH!!!
  if(varcoeff)
    {
    double * splinep = pmean_spline.getV();
    double * workdata = data_forfixed.getV();
    for(i=0;i<likep->get_nrobs();i++,splinep++,workdata++)
      *splinep *= *workdata;
    }
  else
    {
    compute_intercept(beta);
    datamatrix inter = datamatrix(likep->get_nrobs(),1,-intercept);
    pmean_spline.plus(pmean_spline,inter);         // zentrierter durchschnittlicher Spline, Eintr�ge in Reihenfolge x_1,...,x_n
    // fcconst->set_intercept_for_center(intercept);  // beta_0 wurde vorher nicht an zentrierten Spline angepsst, deshalb hier!
    intercept = 0.0;
    }
  datamatrix pmean_fix = datamatrix(likep->get_nrobs(),1,0);
  datamatrix beta_fixx = datamatrix(1,1,beta_fix);
  datamatrix intercept_fix = datamatrix(likep->get_nrobs(),1,alpha_fix);
  if(beta_fix != 0)
    {
    pmean_fix.mult(data_forfixed,beta_fixx);     // berechnet den Anteil der fixen Effekte
    pmean_fix.plus(pmean_fix,intercept_fix);     // zentrierter linearer Effekt; Eintr�ge in Reihenfolge x_1,...,x_n
    }

  subtr_spline();
  spline.plus(pmean_spline,pmean_fix);         //zentrierte durchschnittliche Funktion
  likep->add_linearpred_m(spline,column);      // addiert den Anteil der fixen Effekte zum Gesamtpr�diktor

  // f�r Ausgabe: Vektor "spline" mu� f�r Ausgabe sortiert werden!
  double * splinep;
  double * fchelpbetap = fchelp.get_betameanp();

  if(gridsize < 0)
    {
    vector<int>::iterator freqwork = freqoutput.begin();
    int * workindex = index.getV();
    for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
      {
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        {
        if(!varcoeff)
          *fchelpbetap = spline(*workindex,0) * transform;
        else
          *fchelpbetap = spline(*workindex,0) * transform / data_forfixed(*workindex,0);
        fchelpbetap++;
        }
      }
    }
  else      //???
    {
    for(i=0;i<gridsize;i++,fchelpbetap++,splinep++)
      *fchelpbetap = *splinep;
    }
  }


void FULLCOND_pspline_stepwise::multBS_sort(datamatrix & res, const datamatrix & beta)      // soll f_dach berechnen in Reihenfolge wie Daten f�r fixen Effekt
  {

  double *workres;
  double *workbeta;
  double *workBS;

  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator workindex2 = index2.begin();

  if(varcoeff)
    workBS = B.getV();
  else
    workBS = BS.getV();

  unsigned col = degree+1;
  unsigned j,k;
  int i,stop;

  workres = res.getV();
  for(j=0;j<res.rows()*res.cols();j++,workres++)
    *workres = 0.0;

  i = 0;
  k = 0;
  workres = res.getV() + *workindex2;
  while (k<nrpar)
    {
    stop = lastnonzero[k];
    while (i <= stop)
      {
      workbeta = beta.getV() + k;
      for(j=0;j<col;j++,workBS++,workbeta++)
        *workres += *workBS * *workbeta;
      if((freqwork+1)!=freq.end() && *freqwork==*(freqwork+1))
        {
        workBS -= col;
        workbeta -= col;
        }
      i++;
      freqwork++;
      workindex2++;
      workres += *workindex2;
      }
    k++;
    }

  }

// END: MODEL-AVERAGING --------------------------------------------------------


void FULLCOND_pspline_stepwise::reset_effect(const unsigned & pos)
  {
  subtr_spline();
  unsigned i;
  double * work;
  work = spline.getV();
  for(i=0;i<spline.rows();i++,work++)
    *work = 0.0;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;

  intercept = 0.0;
  }


double FULLCOND_pspline_stepwise::compute_df(void)
  {
  double df = 0;
  //if(inthemodel == false && fixornot == true)
  //  df = 1;
  if(inthemodel == true)
    {
    // falls �bergang: leer --> VC kann fixer Effekt noch nicht enthalten sein, d.h. es mu� ein df addiert werden!
    /*if(varcoeff && !identifiable && !center)
      {
      bool raus = false;
      unsigned j = 1;
      while(j<fcconst->get_datanames().size() && raus==false)
        {
        if(fcconst->get_datanames()[j] == datanames[0]
          || fcconst->get_datanames()[j] == (datanames[0]+"_1"))
          {
          raus = true;
          }
        j = j + 1;
        }
      if(raus == false)
        {
        df += 1;
        }
      }*/

    if(varcoeff && lambda == -2)
      {
      if(identifiable)
        df = 2;
      else
        df = df + 1;
      }
    else
      {
      if (lambdaold == lambda && likep->get_iwlsweights_notchanged() == true)
        {
        df = df_lambdaold;
        }
      else if (lambdaold != lambda || likep->get_iwlsweights_notchanged() == false)
        {
        if (likep->get_iwlsweights_notchanged() == false)
          compute_XWXenv(likep->get_weightiwls(),column);
        if(lambda != lambda_prec || likep->iwlsweights_constant() == false)
          {
          prec_env.addto(XX_env,Kenv,1.0,lambda);
          lambda_prec = lambda;
          }
        invprec = envmatdouble(0,nrpar,prec_env.getBandwidth());
        prec_env.inverse_envelope(invprec);
        df = df + invprec.traceOfProduct(XX_env);
        if(!identifiable)
          df -= 1;
        df_lambdaold = df;
        lambdaold = lambda;
        }
      }
    }

  return df;
  }


void FULLCOND_pspline_stepwise::hierarchie_rw1(vector<double> & untervector, int dfo)
  {

  unsigned number = untervector.size()-1;

  update_stepwise(untervector[0]);
  double df_max = compute_df();

  update_stepwise(untervector[number]);
  double df_min = compute_df();

  if(df_max > dfo && df_min < dfo)
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

        if(df_mitteunten > dfo && df_mitteoben > dfo)
          stelle_unten = stelle/2;
        else if(df_mitteunten < dfo && df_mitteoben < dfo)
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
          if(!varcoeff)
            hilf.push_back(-1);
          else
            hilf.push_back(-2);
          for(i=stelle_oben;i<untervector.size();i++)
            hilf.push_back(untervector[i]);
          untervector = hilf;
          }
        }
     }
  else if(df_min >= dfo)
     {
     if(!varcoeff)
       untervector.push_back(-1);
     else
       untervector.push_back(-2);
     }
  else
     {
     vector<double> hilf;
     if(!varcoeff)
       hilf.push_back(-1);
     else
       hilf.push_back(-2);
     unsigned i;
     for(i=0;i<untervector.size();i++)
        hilf.push_back(untervector[i]);
     untervector = hilf;
     }
  }


void FULLCOND_pspline_stepwise::compute_lambdavec(
vector<double> & lvec, int & number)
  {
  if(number>0)
    {
    if (get_df_equidist()==true && number>1)
       FULLCOND::compute_lambdavec_equi(lvec,number);
    else
       FULLCOND::compute_lambdavec(lvec,number);
    }

  if(!varcoeff)
    {
    if(type==RW1 && number>0)
      hierarchie_rw1(lvec,1);
    else  // if(type==RW2 || (type==RW1 && number==-1))
      lvec.push_back(-1);
    }
  else
    {
    if(type==RW1 && number>0)
      {
      if(identifiable)
        {
        hierarchie_rw1(lvec,2);
        lvec.push_back(-1);
        }
      else
        hierarchie_rw1(lvec,1);
      }
    else  // if(type==RW2 || (type==RW1 && number==-1))
      {
      lvec.push_back(-2);
      if(identifiable)    //VCM_neu
        lvec.push_back(-1);
      }
    }

  get_forced();
  if(forced_into==false)
     lvec.push_back(0);
  }


const datamatrix & FULLCOND_pspline_stepwise::get_data_forfixedeffects(void)
  {

  unsigned nrobs = index.rows();

  if (data_forfixed.rows() < nrobs)
    {
    data_forfixed = datamatrix(nrobs,1);
    int * workindex = index.getV();
    vector<int>::iterator freqwork = freqoutput.begin();
    for(unsigned i=0;i<nrobs;i++,workindex++,freqwork++)
      data_forfixed(*workindex,0) = xvalues(*freqwork,0);
    }

  return data_forfixed;

  }


ST::string FULLCOND_pspline_stepwise::get_effect(void)
  {
  ST::string h;

  if(varcoeff)
    h = datanames[1] + "(psplinerw";
  else
    h = datanames[0] + "(psplinerw";
  if (type== MCMC::RW1)
    h = h + "1,df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";
  else
    h = h + "2,df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";

  return h;
  }


ST::string FULLCOND_pspline_stepwise::get_befehl(void)
  {
  ST::string h;

  if(varcoeff)
    h = datanames[1] + "(psplinerw";
  else
    h = datanames[0] + "(psplinerw";
  if (type== MCMC::RW1)
    h = h + "1";
  else
    h = h + "2";

  h = h + ", lambda=" + ST::doubletostring(lambda,6);
  if(degree!=3)
    h = h + ", degree=" + ST::inttostring(degree);
  if(nrknots!=20)
    h = h + ",nrknots=" + ST::inttostring(nrknots);
  h = h + ")";

  return h;
  }


void FULLCOND_pspline_stepwise::createreml(datamatrix & X,datamatrix & Z,
                                const unsigned & Xpos, const unsigned & Zpos)
  {
  unsigned i,j,k;

  double * workdata;
  double * workZ;
  datamatrix spline1 = datamatrix(spline.rows(),1,0);

// X berechen

  if(type == RW1)
    {
    if(varcoeff)
      {
      double * workX;
      unsigned Xcols = X.cols();

      workX = X.getV()+Xpos;

      double * workintact = data_forfixed.getV();
      for (i=0;i<spline1.rows();i++,workintact++,workX+=Xcols)
        {
        *workX = *workintact;
        }
      }
    }
  else if(type == RW2)
    {
    double * workX;
    unsigned Xcols = X.cols();

    datamatrix knoten = datamatrix(nrpar,1,0.0);
    for(i=0;i<nrpar;i++)
      knoten(i,0) = knot[i];

    multBS_index(spline1,knoten);

    if(!varcoeff)
      {
      double splinemean=spline1.mean(0);
      for(i=0; i<spline1.rows(); i++)
        {
        spline1(i,0)=spline1(i,0)-splinemean;
        }
      }

    if(X.rows()<spline1.rows())
    // category specific covariates
      {
      X_VCM = spline1;

      unsigned nrcat2 = spline1.rows()/X.rows()-1;
      unsigned nrobs=X.rows();
      for(i=0; i<X.rows(); i++)
        {
        for(j=0; j<nrcat2; j++)
          {
          X(i,Xpos+j) = spline1(j*nrobs+i,0) - spline1(nrcat2*nrobs+i,0);
          }
        }
      }
    else
    // no category specific covariates
      {
      workdata = spline1.getV();
      workX = X.getV()+Xpos;

      if(varcoeff)
        {
        double * workintact = data_forfixed.getV();
        double * workX_VCM = X_VCM.getV()+1;
        for (i=0;i<spline1.rows();i++,workdata++,workintact++,workX+=Xcols,workX_VCM+=2)
          {
          *workX = *workintact;
          *(workX+1) = *workdata**workintact;
          *workX_VCM = *workdata;
          }
        }
      else
        {
        for (i=0;i<spline1.rows();i++,workdata++,workX+=Xcols)
          {
          *workX = *workdata;
          }
        }
      }
    }

// Z berechnen

  compute_Kweights();
//  datamatrix diffmatrix = diffmat(dimX+1,nrpar);
  datamatrix diffmatrix = weighteddiffmat(type==RW1?1:2,weight);
  diffmatrix = diffmatrix.transposed()*diffmatrix.transposed().sscp().inverse();

  unsigned Zcols = Z.cols();

  if(Z.rows()<spline1.rows())
    // category specific covariates
    {
    Z_VCM = datamatrix(spline1.rows(), dimZ, 0);
    }

  for(j=0;j<dimZ;j++)
    {
    multBS_index(spline1,diffmatrix.getCol(j));

    if(Z.rows()<spline1.rows())
      // category specific covariates
      {
      unsigned nrcat2 = spline1.rows()/Z.rows()-1;
      unsigned nrobs=Z.rows();
      for(i=0; i<nrobs; i++)
        {
        for(k=0; k<nrcat2; k++)
          {
          Z(i, Zpos + k*dimZ + j) = spline1(k*nrobs+i,0) - spline1(nrcat2*nrobs+i,0);
          }
        }
      for(i=0; i<spline1.rows(); i++)
        {
        Z_VCM(i,j) = spline1(i,0);
        }
      }
    else
      // no category specific covariates
      {
      workdata = spline1.getV();
      workZ = Z.getV()+Zpos+j;

      if(varcoeff)
        {
        double * workintact = data_forfixed.getV();
        double * workZ_VCM = Z_VCM.getV()+j;
        for (i=0;i<spline1.rows();i++,workdata++,workintact++,workZ+=Zcols,workZ_VCM+=dimZ)
          {
          *workZ = *workdata**workintact;
          *workZ_VCM = *workdata;
          }
        }
      else
        {
        for (i=0;i<spline1.rows();i++,workdata++,workZ+=Zcols)
          {
          *workZ = *workdata;
          }
        }
      }
    } 

  }

  
} // end: namespace MCMC
