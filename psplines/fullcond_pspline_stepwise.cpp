
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
  : FULLCOND_pspline_gaussian(o,dp,fcc,d,nrk,degr,kp,ft,monotone,ti,fp,pres,deriv,l,gs,diag,c)
  {

  beta_average.erase(beta_average.begin(),beta_average.end());
  interactions_pointer.erase(interactions_pointer.begin(),interactions_pointer.end());
  lambda_nr = 0;
  //lambdas_local = datamatrix(rankK,1,1);

  if (type==RW1)
    grenzfall = 0;
  else if (type == RW2)
    grenzfall = 1;

  }


// CONSTRUCTOR 2  (for varying coefficients term)

FULLCOND_pspline_stepwise::FULLCOND_pspline_stepwise(MCMCoptions * o, DISTRIBUTION * dp,
                      FULLCOND_const * fcc,const datamatrix & effmod, const datamatrix & intact,
                      const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                      const fieldtype & ft, const ST::string & monotone, const ST::string & ti,
                      const ST::string & fp, const ST::string & pres, const bool & deriv,
                      const double & l, const int & gs, const unsigned & c)
  : FULLCOND_pspline_gaussian(o,dp,fcc,effmod,intact,nrk,degr,kp,ft,monotone,ti,fp,pres,false,l,gs,c)
  {

  beta_average.erase(beta_average.begin(),beta_average.end());
  interactions_pointer.erase(interactions_pointer.begin(),interactions_pointer.end());
  lambda_nr = 0;

  data_forfixed = intact;

  if (type==RW1)
    grenzfall = 1;
  else if (type == RW2)
    grenzfall = 2;

  }

  // COPY CONSTRUCTOR

FULLCOND_pspline_stepwise::FULLCOND_pspline_stepwise(const FULLCOND_pspline_stepwise & fc)
  : FULLCOND_pspline_gaussian(FULLCOND_pspline_gaussian(fc))
  {
  beta_average = fc.beta_average;
  interactions_pointer = fc.interactions_pointer;
  lambda_nr = fc.lambda_nr;
  lambdas_local = fc.lambdas_local;
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const FULLCOND_pspline_stepwise & FULLCOND_pspline_stepwise::operator=(
                                            const FULLCOND_pspline_stepwise & fc)
  {
  if (this == &fc)
    return *this;
  FULLCOND_pspline_gaussian::operator=(FULLCOND_pspline_gaussian(fc));

  beta_average = fc.beta_average;
  interactions_pointer = fc.interactions_pointer;
  lambda_nr = fc.lambda_nr;
  lambdas_local = fc.lambdas_local;

  return *this;
  }


bool FULLCOND_pspline_stepwise::posteriormode(void)
  {

  unsigned i;

  transform = likep->get_trmult(column);
  fchelp.set_transform(transform);

  likep->substr_linearpred_m(spline,column,true);

  compute_XWXenv(likep->get_weightiwls(),column);
  prec_env.addto(XX_env,Kenv,1.0,lambda);
  lambda_prec = lambda;

  likep->compute_workingresiduals(column);
  compute_XWtildey(likep->get_weightiwls(),likep->get_workingresiduals(),1.0,column);

  prec_env.solve(muy,beta);
  add_linearpred_multBS();

  if(interactions_pointer.size()>0)
      search_for_interaction();

  if(center)
    {
    compute_intercept();
    if(interaction==false)
      fcconst->posteriormode_intercept(intercept);
    else                                           
      fcconst->update_intercept(intercept);
    }

  if(interaction == false)
    {
    if(!varcoeff)
      {
      int * workindex = index.getV();
      for(i=0;i<spline.rows();i++,workindex++)
        spline(*workindex,0) -= intercept;
      }
    double * fchelpbetap = fchelp.getbetapointer();

    if(gridsize < 0)
      {
      if(varcoeff)
        multBS(splinehelp,beta);

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

  }


// BEGIN: FÜR INTERAKTIONEN-----------------------------------------------------

bool FULLCOND_pspline_stepwise::changeposterior(const datamatrix & main,const double & inter)
  {
  unsigned i;

  vector<int>::iterator freqwork = freq.begin();
  int * workindex = index.getV();

// spline ändern
  for(i=0;i<spline.rows();i++,freqwork++,workindex++)
    spline(*workindex,0) += main(*freqwork,0) - intercept - inter;

  intercept = 0.0;

// fchelp ändern
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


bool FULLCOND_pspline_stepwise::changeposterior2(const datamatrix & main,const double & inter)       
  {
  unsigned i;

  vector<int>::iterator freqwork = freq.begin();
  int * workindex = index.getV();

  datamatrix hilf = datamatrix(spline.rows(),1,0);

// spline ändern
  for(i=0;i<spline.rows();i++,freqwork++,workindex++)
    {
    spline(*workindex,0) += main(*freqwork,0) + inter;
    hilf(*workindex,0) += main(*freqwork,0) + inter;
    }

  likep->add_linearpred(hilf);

// Intercept ändern
  intercept = 0.0;     

  return FULLCOND_nonp_basis::posteriormode();
  }


void FULLCOND_pspline_stepwise::set_pointer_to_interaction(FULLCOND * inter)
  {
  interactions_pointer.push_back(inter);
  }

void FULLCOND_pspline_stepwise::search_for_interaction(void)
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
  if(thereis == true)
    interaction = true;
  else
    interaction = false;
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


// BEGIN: MODEL-AVERAGING ------------------------------------------------------

void FULLCOND_pspline_stepwise::save_betas(vector<double> & modell, unsigned & anzahl)
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
      if(!varcoeff)          // Gerade so zentrieren, daß Integral=0, d.h. G((x_max-x_min)/2) = 0!
        alpha_fix += - beta_average[i][0] * crit_weights[i] * (data_forfixed.max(0)+data_forfixed.min(0))/2;
      }
    }
  fcconst->set_intercept_for_center(-alpha_fix);

  setbeta(beta_spline.size(),1,0);
  double * workbeta = beta.getV();
  for(i=0;i<beta_spline.size();i++,workbeta++)
    *workbeta = beta_spline[i];
  datamatrix pmean_spline = datamatrix(likep->get_nrobs(),1,0);
  if(gridsize < 0)
    multBS_sort(pmean_spline,beta);
  else    // Vorsicht: hier passt es wahrscheinlich überhaupt nicht mehr, weil nicht alle x_1,...,x_n verwendet werden!!!
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
    pmean_spline.plus(pmean_spline,inter);         // zentrierter durchschnittlicher Spline, Einträge in Reihenfolge x_1,...,x_n
    // fcconst->set_intercept_for_center(intercept);  // beta_0 wurde vorher nicht an zentrierten Spline angepsst, deshalb hier!
    intercept = 0.0;
    }
  datamatrix pmean_fix = datamatrix(likep->get_nrobs(),1,0);
  datamatrix beta_fixx = datamatrix(1,1,beta_fix);
  datamatrix intercept_fix = datamatrix(likep->get_nrobs(),1,alpha_fix);
  if(beta_fix != 0)
    {
    pmean_fix.mult(data_forfixed,beta_fixx);     // berechnet den Anteil der fixen Effekte
    pmean_fix.plus(pmean_fix,intercept_fix);     // zentrierter linearer Effekt; Einträge in Reihenfolge x_1,...,x_n
    }

  subtr_spline();
  spline.plus(pmean_spline,pmean_fix);         //zentrierte durchschnittliche Funktion
  likep->add_linearpred_m(spline,column);      // addiert den Anteil der fixen Effekte zum Gesamtprädiktor

  // für Ausgabe: Vektor "spline" muß für Ausgabe sortiert werden!
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


void FULLCOND_pspline_stepwise::multBS_sort(datamatrix & res, const datamatrix & beta)      // soll f_dach berechnen in Reihenfolge wie Daten für fixen Effekt
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


void FULLCOND_pspline_stepwise::hierarchie_rw1(vector<double> & untervector)
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


void FULLCOND_pspline_stepwise::compute_lambdavec(      
vector<double> & lvec, int & number)
  {
  if(number>0)
    {
    if (get_df_equidist()==true)
       FULLCOND::compute_lambdavec_equi(lvec,number);
    else
       FULLCOND::compute_lambdavec(lvec,number);
    }

  if(type==RW1 && number>0)
    hierarchie_rw1(lvec);
  else  // if(varcoeff || type==RW2 || (type==RW1 && number==-1))
    lvec.push_back(-1);

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
    h = datanames[1] + "*" + datanames[0] + "(psplinerw";
  else
    h = datanames[0] + "(psplinerw";
  if (type== MCMC::RW1)
    h = h + "1,df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";
  else
    h = h + "2,df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";

  return h;
  }

  
} // end: namespace MCMC
