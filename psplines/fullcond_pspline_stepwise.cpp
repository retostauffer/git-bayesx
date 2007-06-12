
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

  utype = gaussian;
  data_forfixed = d;

  beta_average.erase(beta_average.begin(),beta_average.end());
  lambda_nr = -1;
  //lambdas_local = datamatrix(rankK,1,1);

  all_precenv.erase(all_precenv.begin(),all_precenv.end());
  lambdavec.erase(lambdavec.begin(),lambdavec.end());

  if(increasing || decreasing)
    {
    Menv = Krw1env(weight);
    g = datamatrix(nrpar,1,0.0);
    updateMenv();
    }
  if(monotone == "convex")
    convex = true;
  else
    convex = false;
  if(monotone == "concave")
    concave = true;
  else
    concave = false;
  if(concave || convex)
    {
    Menv = Krw2env(weight);
    g = datamatrix(nrpar,1,0.0);
    //g = weight;
    F1 = datamatrix(nrpar,1,0);
    F2 = datamatrix(nrpar,1,0);
    unsigned s;
    for (s=2;s<nrpar;s++)
      {
      double w1 = weight[s-1];
      double w2 = weight[s];
      if(type==RW1)
        {
        w1 *= 0.5;
        w2 *= 0.5;
        }
      F1(s,0) = -(1+w2/w1);
      F2(s,0) = w2/w1;
      //g(s,0) = w2*(1+w2/w1);
      }
    updateMenv();
    }

  if(increasing || decreasing)
    {
    unsigned nrobs = index.rows();
    if (data_varcoeff_fix.rows() < nrobs)
      {
      data_varcoeff_fix = datamatrix(nrobs,2,1);
      int * workindex = index.getV();
      vector<int>::iterator freqwork = freqoutput.begin();
      for(unsigned i=0;i<nrobs;i++,workindex++,freqwork++)
        {
        data_varcoeff_fix(i,0) = 1;
        data_varcoeff_fix(i,1) = d(i,0);
        }
      }
    XVX = datamatrix(2,2,0);
    }

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
                      const double & l, const int & gs, const bool & vccent, const unsigned & c)
  : FULLCOND_pspline_gaussian(o,dp,fcc,effmod,intact,nrk,degr,kp,ft,monotone,ti,fp,pres,false,l,gs,c)
  {

  utype = gaussian;  
  beta_average.erase(beta_average.begin(),beta_average.end());
  lambda_nr = 0;

  all_precenv.erase(all_precenv.begin(),all_precenv.end());
  lambdavec.erase(lambdavec.begin(),lambdavec.end());

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

  if(increasing || decreasing)
    {
    Menv = Krw1env(weight);
    g = datamatrix(nrpar,1,0.0);
    updateMenv();
    }
  if(monotone == "convex")
    convex = true;
  else
    convex = false;
  if(monotone == "concave")
    concave = true;
  else
    concave = false;
  if(concave || convex)
    {
    Menv = Krw2env(weight);
    g = datamatrix(nrpar,1,0.0);
    F1 = datamatrix(nrpar,1,0);
    F2 = datamatrix(nrpar,1,0);
    unsigned s;
    for (s=2;s<nrpar;s++)
      {
      double w1 = weight[s-1];
      double w2 = weight[s];
      if(type==RW1)
        {
        w1 *= 0.5;
        w2 *= 0.5;
        }
      F1(s,0) = -(1+w2/w1);
      F2(s,0) = w2/w1;
      }
    updateMenv();
    }

  XVX = datamatrix(2,2,0);

  if (type==RW1)
    grenzfall = 1;
  else if (type == RW2)
    grenzfall = 2;

  //VCM_neu
  if(vccent == true)
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
  XVX = fc.XVX;
  df_lambdaold = fc.df_lambdaold;
  lambdaold = fc.lambdaold;
  all_precenv = fc.all_precenv;
  lambdavec = fc.lambdavec;
  Menv = fc.Menv;
  convex = fc.convex;
  concave = fc.concave;
  lambdamono = fc.lambdamono;
  fc_df = fc.fc_df;
  utype = fc.utype;
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
  XVX = fc.XVX;
  df_lambdaold = fc.df_lambdaold;
  lambdaold = fc.lambdaold;
  all_precenv = fc.all_precenv;
  lambdavec = fc.lambdavec;
  Menv = fc.Menv;
  convex = fc.convex;
  concave = fc.concave;
  lambdamono = fc.lambdamono;
  fc_df = fc.fc_df;
  utype = fc.utype;

  return *this;
  }


bool FULLCOND_pspline_stepwise::posteriormode(void)
  {

  unsigned i;
  transform = likep->get_trmult(column);
  fchelp.set_transform(transform);

  bool interaction2 = false;
  if(interactions_pointer.size()>0)
    interaction2 = search_for_interaction();

    //  Spezialfall "lambda = -2": - linearer Effekt bei Monotoniebedingung
    //                             - Gerade als Funktion g(x) bei VC g(x)*z
  if((varcoeff || decreasing || increasing) && lambda == -2)
    {
    likep->substr_linearpred_m(spline,column,true);

    datamatrix betas = datamatrix(2,1,0);
    //if(likep->iwlsweights_constant() == false || (XVX(1,1)==0 && XVX(0,0)==0))
    if(calculate_xwx_vc == true || (XVX(1,1)==0 && XVX(0,0)==0))
      {
      calculate_xwx_vc = false;
      likep->fisher(XVX,data_varcoeff_fix,column);
      XVX.assign((XVX.cinverse()));
      }
    likep->compute_weightiwls_workingresiduals(column);
    betas = XVX*data_varcoeff_fix.transposed()*likep->get_workingresiduals();

    // wenn das Vorzeichen bei Monotoniebedingung falsch ist, wird Steigung = 0 gesetzt!
    if(!varcoeff && ((increasing && betas(1,0) < 0) || (decreasing && betas(1,0) > 0)))
      {
      betas = datamatrix(2,1,0);
      double hilf = likep->get_workingresiduals().sum(0);
      hilf = hilf / likep->get_weightiwls().sum(0);
      betas(0,0) = hilf;
      spline = datamatrix(spline.rows(),1,hilf);
      }
    else
      {
      spline.mult(data_varcoeff_fix,betas);
      }

     likep->substr_linearpred_m(-spline,column,true);

    if(center && varcoeff)
      {
      intercept = betas(0,0) + 0.5*betas(1,0)*(effmodi.max(0)+effmodi.min(0));
      int * workindex = index.getV();
      for(i=0;i<spline.rows();i++,workindex++)
        spline(*workindex,0) -= intercept*data_forfixed(*workindex,0);
      betas(0,0) -= intercept;
      update_fix_effect();
      intercept = 0.0;
      }
    else if(center && !varcoeff)
      {
      int * workindex = index.getV();
      for(i=0;i<spline.rows();i++,workindex++)
        spline(*workindex,0) -= betas(0,0);
      fcconst->update_intercept(betas(0,0));
      betas(0,0) = 0;
      }

    double * fchelpbetap = fchelp.getbetapointer();
    datamatrix help = datamatrix(beta.rows(),1,0);
    unsigned j;
    if(gridsize<0)
      {
      vector<int>::iterator freqwork = freqoutput.begin();
      int * workindex = index.getV();
      j = 0;
      for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
        {
        if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
          {
          if(varcoeff)
            *fchelpbetap = betas(0,0) + effmodi(*workindex,0)*betas(1,0);
          else
            *fchelpbetap = betas(0,0) + data_forfixed(*workindex,0)*betas(1,0);
          while(j<beta.rows())
            {
            help(j,0) = *fchelpbetap;
            j += 1;
            }
          fchelpbetap++;
          }
        }
      }
    else
      {
      j = 0;
      for(i=0;i<unsigned(gridsize);i++,fchelpbetap++)
        {
        *fchelpbetap = betas(0,0) + xvalues(i,0)*betas(1,0);   // keine Fallunterscheidung VC / nicht VC nötig!
        while(j<beta.rows())
          {
          help(j,0) = *fchelpbetap;
          j += 1;
          }
        }
      }

    beta.assign(help);
    return fchelp.posteriormode();
    //return FULLCOND_nonp_basis::posteriormode();
    }
  else // if(lambda != -2)
    {
    if(!interaction || !interaction2)  // nur, wenn nicht Teil einer Interaktion, sonst wird "posteriormode"
      {                                // in "fullcond_pspline_surf_stepwise" aufgerufen!
      likep->substr_linearpred_m(spline,column,true);

      if(!decreasing && !increasing && !convex && !concave)   // "normale" Funktion ohne Restriktion
        {
        if ( (lambda_prec != lambda) || (calculate_xwx == true))
          {
          if(calculate_xwx == true)
            {
            calculate_xwx = false;
            compute_XWXenv(likep->get_weightiwls(),column);
            }
          prec_env.addto(XX_env,Kenv,1.0,lambda);
          lambda_prec = lambda;
          }
        likep->compute_workingresiduals(column);
        compute_XWtildey(likep->get_weightiwls(),likep->get_workingresiduals(),1.0,column);

        prec_env.solve(muy,beta);
        add_linearpred_multBS();
        }

      else if(decreasing || increasing || convex || concave)
        {
        g = datamatrix(nrpar,1,0.0);
        updateMenv();
        datamatrix bold = datamatrix(nrpar,1,0);
        double nold = 0.01;
        bool first = true;
        unsigned count = 0;

        while(count < 200 && (bold==datamatrix(nrpar,1,0) || norm(beta-bold)/nold > 0.00001))
          {
          count += 1;
          if(!first)
            likep->substr_linearpred_m(spline,column,true);

          bold = beta;
          if(calculate_xwx == true && first==true)
            {
            calculate_xwx = false;
            compute_XWXenv(likep->get_weightiwls(),column);
            }
          prec_env.addto(XX_env,Kenv,1.0,lambda);
          prec_env.addto(prec_env,Menv,1.0,lambdamono);
          lambda_prec = lambda;

          likep->compute_workingresiduals(column);
          compute_XWtildey(likep->get_weightiwls(),likep->get_workingresiduals(),1.0,column);
          prec_env.solve(muy,beta);
          add_linearpred_multBS();

          // Versuch mit Monotonie:
          g = datamatrix(nrpar,1,0.0);
          if(increasing || decreasing)
            {
            for(i=1;i<beta.rows();i++)
              {
              double diff = beta(i,0)-beta(i-1,0);
              if(increasing && diff < 0.0)
                g(i,0) = 1.0;
              if(decreasing && diff > 0.0)
                g(i,0) = 1.0;
              }
            }
          else
            {
            for(i=2;i<beta.rows();i++)
              {
              double diff = beta(i,0) - 2*beta(i-1,0) + beta(i-2,0);
              if(knpos == quantiles)
                diff = (beta(i,0) - beta(i-1,0))/weight[i] - (beta(i-1,0) - beta(i-2,0))/weight[i-1];
              if(convex && diff < 0.0)
                 g(i,0) = 1/(weight[i]*weight[i]);
                 //g(i,0) = weight[i] * (1+weight[i]/weight[i-1]);
                 //g(i,0) = weight[i];
                 //g(i,0) = 1.0;
              if(concave && diff > 0.0)
                 g(i,0) = 1/(weight[i]*weight[i]);
                 //g(i,0) = weight[i] * (1+weight[i]/weight[i-1]);
                 //g(i,0) = weight[i];
                 //g(i,0) = 1.0;
              }
            }

          updateMenv();
          nold = norm(bold);
          if(nold == 0)
            nold = 0.001;
          first = false;
          }
        if(count == 200)
          optionsp->out("  NOTE: Algorithm for ensuring monotonicity restriction did not converge!");
        }

      if(center)
        {
        compute_intercept();
        if(!varcoeff)
          fcconst->posteriormode_intercept(intercept);
        else
          update_fix_effect();

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

      bool converged = FULLCOND_nonp_basis::posteriormode();
      if(converged)
        {
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
          for(i=0;i<unsigned(gridsize);i++,fchelpbetap++)
            *fchelpbetap = splinehelp(i,0) - intercept;
          }
        fchelp.posteriormode();
        }

      intercept = 0.0;

      return converged;
      }  // end: if(interaction == false)
    else if(!center && interaction && interaction2)  // Zum Bestimmen der Lambdas (es gilt: "!center")
      {                                              // sonst überflüssig!
      return interactions_pointer[interactions_pointer.size()-1]->posteriormode();
      }
    else
      {
      return true;
      }
    }  // END: else if(lambda != 2)
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
  ST::string name_richtig = datanames[1];  // (VC alt) 0
  while(j<fcconst->get_datanames().size() && raus==false)
     {
     if(fcconst->get_datanames()[j] == datanames[1])  // (VC alt) 0
        {
        raus = true;
        }
     if(fcconst->get_datanames()[j] == (datanames[1]+"_1"))   // (VC alt) 0
        {
        raus = true;
        name_richtig = datanames[1] + "_1";    // (VC alt) 0
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
    interactions_pointer[0]->set_inthemodel(-1);
    fcconst->update_fix_effect(j,intercept,data_forfixed);
    }
  }


void FULLCOND_pspline_stepwise::const_varcoeff(void)
  {
  if(varcoeff)
    fcconst->posteriormode_const_varcoeff(data_forfixed);
  }

// END: For Varying Coefficients -----------------------------------------------

// BEGIN: FÜR INTERAKTIONEN-----------------------------------------------------

bool FULLCOND_pspline_stepwise::changeposterior3(const datamatrix & betamain,const datamatrix & main,const double & inter)
   {
/*  unsigned i;
  vector<int>::iterator freqwork = freq.begin();
  int * workindex = index.getV();

// spline ändern
  for(i=0;i<spline.rows();i++,freqwork++,workindex++)
    spline(*workindex,0) += main(*freqwork,0) - inter;

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

  return fchelp.posteriormode();     */

  unsigned i;

// beta, spline ändern
  beta.assign(betamain);
  for(i=0;i<nrpar;i++)
    beta(i,0) -= inter;

  vector<int>::iterator freqwork = freq.begin();
  int * workindex = index.getV();
  for(i=0;i<spline.rows();i++,freqwork++,workindex++)
    spline(*workindex,0) += main(*freqwork,0) - inter;

  //betaold.assign(beta);

  bool converged = FULLCOND_nonp_basis::posteriormode();
  if(converged)
    {
    // fchelp ändern
    write_spline();
    fchelp.posteriormode();
    }

  return converged;
  }


bool FULLCOND_pspline_stepwise::changeposterior_varcoeff(const datamatrix & betamain,
                                        const datamatrix & main,const double & inter)
  {
  unsigned i;

// beta ändern
  beta.assign(betamain);
  for(i=0;i<nrpar;i++)
    beta(i,0) -= inter;

  vector<int>::iterator freqwork = freqoutput.begin();
  int * workindex = index.getV();

// splinehelp ändern
  for(i=0;i<splinehelp.rows();i++,freqwork++,workindex++)
    {
    splinehelp(*workindex,0) = main(*freqwork,0) - inter;
    }
  double * worksp = spline.getV();
  double * worksph = splinehelp.getV();
  double * workint = data_forfixed.getV();
  for(i=0;i<spline.rows();i++,worksp++,worksph++,workint++)
    {
    *worksp = *worksph * *workint;
    }

  bool converged = FULLCOND_nonp_basis::posteriormode();
  if(converged)
    {
    // fchelp ändern
    double * fchelpbetap = fchelp.getbetapointer();
    freqwork = freqoutput.begin();
    workindex = index.getV();
    for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
      {
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        {
        *fchelpbetap = splinehelp(*workindex,0);
        fchelpbetap++;
        }
      }
    fchelp.posteriormode();      
    }

  return converged;
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
//  unsigned i;
  bool thereis = false;
//  for(i=0;i<interactions_pointer.size();i++)
//    {
    bool drin, fix;         // 2d Interaktion muß an der letzten Stelle sein!!!
    interactions_pointer[interactions_pointer.size()-1]->get_inthemodel(drin,fix);
    if(drin == true)  
      thereis = true;
//    }
  return thereis;
  }


void FULLCOND_pspline_stepwise::hierarchical(ST::string & possible)
  {
  unsigned i;
  bool spline = false;
  bool fix = false;
  bool spline1, fix1;
  if(!varcoeff)
    {
    if(interaction)
      {       // Zeiger auf 2d-Interaktion ist an der letzten Stelle!!!
      interactions_pointer[interactions_pointer.size()-1]->get_inthemodel(spline,fix);
      if(spline == true)
        possible = "spline";
      else if(fix == true && spline == false)
        possible = "spfix";
      else
        possible = "alles";
      }
    else if(number == -1)  // für Interaktionsvariable "z" bei VC (g(x)*z + h(y)*z + b*z)
      {                    // number = -1 als Kennzeichen für "z" (theoretisch nicht eindeutig)
      for(i=0;i<interactions_pointer.size();i++)
        {
        interactions_pointer[i]->get_inthemodel(spline1,fix1);
        if(spline1 == true)
          spline = true;
        if(fix1 == true)
          fix = true;
        }

      if(spline == true)
        possible = "vfix";
      else
        possible = "alles";
      }
    }
  else    // nur für Unterscheidung VC / !VC
    {
    if(interaction)
      {       // Zeiger auf 2d-Interaktion ist an der letzten Stelle!!!
      interactions_pointer[interactions_pointer.size()-1]->get_inthemodel(spline,fix);
      if(spline == true)
        possible = "vspline";
      else if(fix == true && spline == false)
        possible = "vspfix";
      else
        possible = "valles";
      }
    else
      possible = "valles";
    }
  }


// BEGIN: MODEL-AVERAGING ------------------------------------------------------

/*void FULLCOND_pspline_stepwise::save_betas(vector<double> & modell, int & anzahl)
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
  } */


/*void FULLCOND_pspline_stepwise::average_posteriormode(vector<double> & crit_weights)
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
  double helpdouble = -alpha_fix;
  fcconst->set_intercept_for_center(helpdouble);

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
  }     */


/*void FULLCOND_pspline_stepwise::multBS_sort(datamatrix & res, const datamatrix & beta)      // soll f_dach berechnen in Reihenfolge wie Daten für fixen Effekt
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
           */
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


void FULLCOND_pspline_stepwise::reset(void)
  {
  spline = datamatrix(spline.rows(),spline.cols(),0);
  FULLCOND::reset();
  }

void FULLCOND_pspline_stepwise::update_stepwise(double la)         // neu!!!
  {
  if(smoothing == "global")
    {
    lambda=la;

    if(likep->iwlsweights_constant() == true)
      {
      bool gefunden = false;
      unsigned i = 0;
      while(i<lambdavec.size() && gefunden == false)
        {
        if(lambda == lambdavec[i])
          gefunden = true;
        i++;
        }
      if(gefunden == true)
        {
        prec_env = all_precenv[i-1];
        lambda_prec = lambda;
        }
      } 
    }
  else
    {
    lambda_prec = -1;
    lambdaold = -1;
    lambda = 1;
    lambdas_local(lambda_nr+nrpar-rankK,0) = 1/la;
// lineare Interpolation der log. Lambdas für dazwischenliegende Lambdas:
/*if( (lambda_nr+nrpar-rankK) < rankK )
  {
  double lambda0 = log10(1/lambdas_local(lambda_nr+nrpar-rankK+2,0));
  lambdas_local(lambda_nr+nrpar-rankK+1,0) = 1/(pow(10,0.5*(lambda0+log10(la))));
  }
if( (lambda_nr+nrpar-rankK) > (nrpar-rankK))
  {
  double lambda0 = log10(1/lambdas_local(lambda_nr+nrpar-rankK-2,0));
  lambdas_local(lambda_nr+nrpar-rankK-1,0) = 1/(pow(10,0.5*(lambda0+log10(la))));
  }
*/
    updateK(lambdas_local);
    }
  }


double FULLCOND_pspline_stepwise::compute_df(void)
  {
  double df = 0;
  //if(inthemodel == false && fixornot == true)
  //  df = 1;
  if(inthemodel == true)
    {
    // falls Übergang: leer --> VC kann fixer Effekt noch nicht enthalten sein, d.h. es muß ein df addiert werden!
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
    else if(!varcoeff && lambda == -2)
      {
      df = 1;
      }
    else
      {
      if (lambdaold == lambda && likep->get_iwlsweights_notchanged() == true && !increasing && !decreasing && !convex && !concave)
        {
        df = df_lambdaold;
        }
      else if (lambdaold != lambda || likep->get_iwlsweights_notchanged() == false || increasing || decreasing || convex || concave)
        {
        if (calculate_xwx == true)
          compute_XWXenv(likep->get_weightiwls(),column);
        if(lambda != lambda_prec || calculate_xwx == true || increasing || decreasing || convex || concave)
          {
          calculate_xwx = false;
          prec_env.addto(XX_env,Kenv,1.0,lambda);
          if(increasing || decreasing || convex || concave)
            prec_env.addto(prec_env,Menv,1.0,lambdamono);
          lambda_prec = lambda;
          }
        invprec = envmatdouble(0,nrpar,prec_env.getBandwidth());
        prec_env.inverse_envelope(invprec);
        df = df + invprec.traceOfProduct(XX_env);
        if(!identifiable)
          df -= 1;

if(smoothing == "local")
  df = df/rankK;

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
          if(varcoeff || increasing || decreasing)
            hilf.push_back(-2);
          else
            hilf.push_back(-1);
          for(i=stelle_oben;i<untervector.size();i++)
            hilf.push_back(untervector[i]);
          untervector = hilf;
          }
        }
     }
  else if(df_min >= dfo)
     {
     if(varcoeff || increasing || decreasing)
       untervector.push_back(-2);
     else
       untervector.push_back(-1);
     }
  else
     {
     vector<double> hilf;
     if(varcoeff || increasing || decreasing)
       hilf.push_back(-2);
     else
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
  lambdaold = -1;
  if(number>0)
    {
    if (df_equidist==true && spfromdf==true)
       FULLCOND::compute_lambdavec_equi(lvec,number);
    else
       FULLCOND::compute_lambdavec(lvec,number);
    }

  if(likep->iwlsweights_constant() == true)
    {
    lambdavec = lvec;
    compute_XWXenv(likep->get_weightiwls(),column);
    for(unsigned i=0;i<lambdavec.size();i++)
      {
      prec_env.addto(XX_env,Kenv,1.0,lambdavec[i]);
      prec_env.decomp();
      all_precenv.push_back(prec_env);
      }
    }

  if(!nofixed && !varcoeff && !decreasing && !increasing)
    {
    if(type==RW1 && number>0)
      hierarchie_rw1(lvec,1);
    else  // if(type==RW2 || (type==RW1 && number==-1))
      lvec.push_back(-1);
    }
  else if(!nofixed && !varcoeff && (decreasing || increasing))
    {
    if(type==RW1 && number>0)
      hierarchie_rw1(lvec,1);
    else  // if(type==RW2 || (type==RW1 && number==-1))
      lvec.push_back(-2);
    }
  else if(!nofixed && varcoeff)
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

  if(forced_into==false)
     lvec.push_back(0);

  // Startwert für lambda aus df:
  if(spfromdf==true)
    {
    double lambdavorg = 1000;
    if(!nofixed && !varcoeff && !increasing && !decreasing)
      {
      if(dfstart==1)
        lambdastart = -1;
      else if(dfstart==0)
        lambdastart = 0;
      else
        lambdastart = lambda_from_df(dfstart,lambdavorg);
      }
    else if(!nofixed && !varcoeff && (increasing || decreasing))
      {
      if(dfstart==1)
        lambdastart = -2;
      else if(dfstart==0)
        lambdastart = 0;
      else
        lambdastart = lambda_from_df(dfstart,lambdavorg);
      }
    else
      {
      if(!nofixed && dfstart==1 && identifiable)
        lambdastart = -1;
      else if(!nofixed && (dfstart==2 && identifiable) || (dfstart==1 && !identifiable))
        lambdastart = -2;
      else if(dfstart==0)
        lambdastart = 0;
      else
        lambdastart = lambda_from_df(dfstart,lambdavorg);
      }
    if(lambdastart==-9 || lambdastart==1000000000)    // falls dfstart nicht erreicht werden kann
      lambdastart = 0;
    }
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
    h = datanames[1] + "*" + datanames[0] + "(psplinerw";   // (VC alt) 0 - 1
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
    h = datanames[1] + "*" + datanames[0] + "(psplinerw";    // (VC alt) 0 - 1
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


void FULLCOND_pspline_stepwise::init_names(const vector<ST::string> & na)
  {

  FULLCOND::init_names(na);

  ST::string underscore = "\\_";
  ST::string helpname0 = na[0].insert_string_char('_',underscore);   // (VC alt) 1
  ST::string helpname1 = na[1].insert_string_char('_',underscore);   // (VC alt) 0
  term_symbolic = "f_{" + helpname0 + "}(" + helpname0 + ")" + " \\cdot " + helpname1;
  priorassumptions.push_back("$" + term_symbolic + "$");

  if(column > 0)
    //term_symbolic = term_symbolic + " (" + ST::inttostring(column+1) + ". response category)";
    priorassumptions.push_back("$" + term_symbolic
           + " (" + ST::inttostring(column+1) + ". \\mbox{ } response \\mbox{ } category)$:");

  if(type==MCMC::RW1)
     priorassumptions.push_back("P-spline with first order random walk penalty");
  else if(type==MCMC::RW2)
     priorassumptions.push_back("P-spline with second order random walk penalty");

  ST::string knotstr;
  if (knpos == equidistant)
    knotstr = "equidistant";
  else if (knpos == quantiles)
    knotstr = "quantiles";
  priorassumptions.push_back("Number of knots: " + ST::inttostring(nrknots));
  priorassumptions.push_back("Knot choice: " + knotstr);
  priorassumptions.push_back("Degree of Splines: " + ST::inttostring(degree));
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


void FULLCOND_pspline_stepwise::updateMenv(void)
  {
  unsigned i;  
  if(decreasing || increasing)
    {
    vector<double>::iterator workdiag = Menv.getDiagIterator();
    vector<double>::iterator workenv = Menv.getEnvIterator();
    double * workg = g.getV()+1;
    double workq = weight[1];
    if(type==RW2)
      workq *= 2;
    double wold = 0.0;
    wold=1.0/workq * *workg;
    double wnew;

    *workdiag = wold;
    *workenv = -wold;

    workdiag++;
    workenv++;
    workq++;
    workg++;
    for (i=1;i<nrpar-1;i++,workdiag++,workenv++,workg++)
      {
      workq = weight[i+1];
      if(type==RW2)
        workq *= 2;
      wnew = 0.0;
      wnew = 1.0/workq * *workg;
      *workdiag = wold+wnew;
      *workenv = -wnew;
      wold = wnew;
      }

    *workdiag = wold;
    }
  else
    {
    vector<double>::iterator workdiag=Menv.getDiagIterator();
    vector<double>::iterator workenv = Menv.getEnvIterator();

    *workdiag = (F2(2,0)*F2(2,0))*g(2,0);
    workdiag++;
    *workdiag = (F1(2,0)*F1(2,0))*g(2,0) + (F2(3,0)*F2(3,0))*g(3,0);
    workdiag++;

    *workenv = (F1(2,0)*F2(2,0))*g(2,0);                 //(2,1)
    workenv++;
    *workenv = F2(2,0)*g(2,0);                           //(3,1)
    workenv++;

    *workenv = F1(2,0)*g(2,0) + (F1(3,0)*F2(3,0))*g(3,0);     //(3,2)
    workenv++;

    *workenv = F2(3,0)*g(3,0);                           //(4,2)
    workenv++;

    for(i=2;i<nrpar-2;i++,workdiag++,workenv++)
      {
      *workdiag = (F1(i+1,0)*F1(i+1,0))*g(i+1,0) + (F2(i+2,0)*F2(i+2,0))*g(i+2,0) + 1*g(i,0);

      *workenv = F1(i+1,0)*g(i+1,0) + (F1(i+2,0)*F2(i+2,0))*g(i+2,0);      //(i+1,i)
      workenv++;

      *workenv = F2(i+2,0)*g(i+2,0);                    //(i+2,i)

      }

    *workdiag = (F1(nrpar-1,0)*F1(nrpar-1,0))*g(nrpar-1,0) + 1*g(nrpar-2,0);
    workdiag++;
    *workdiag = 1*g(nrpar-1,0);

    *workenv = F1(nrpar-1,0)*g(nrpar-1,0);
    }

/*ofstream out1("c:\\cprog\\test\\results\\Menv.txt");
Menv.print2(out1);
out1 << endl;
beta.prettyPrint(out1);
out1 << endl;
for(unsigned z=0;z<weight.size();z++)
  out1 << ST::doubletostring(weight[z]) << "  ";
out1 << endl;
Kenv.print2(out1);
out1.close(); */

    }

void FULLCOND_pspline_stepwise::update_bootstrap(const bool & uncond)
  {
  update_bootstrap_df();

  if(fixornot==true)
    {
    bool raus = false;
    unsigned j = 1;
    ST::string name_richtig;
    if(!varcoeff)
       name_richtig = datanames[0];
    else
      name_richtig = datanames[1];  // (VC alt) 0
    while(j<fcconst->get_datanames().size() && raus==false)
      {
      if(fcconst->get_datanames()[j] == name_richtig)
        raus = true;
      j = j + 1;
      }
    unsigned index_fix = j-1;
    double fix = fcconst->getbeta(index_fix,0);
    unsigned i;
    double korrektur = 0;
    if(center)
      korrektur = -0.5*fix*(data_forfixed.max(0) + data_forfixed.min(0));
    double * fchelpbetap = fchelp.getbetapointer();
    if(gridsize < 0)                              // alle verschiedene Beobachtungen
      {
      vector<int>::iterator freqwork = freqoutput.begin();
      int * workindex = index.getV();
      for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
        {
        if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
          {
          if(!varcoeff)
            *fchelpbetap = fix * data_forfixed(*workindex,0) + korrektur;
          else
            *fchelpbetap = fix;
          fchelpbetap++;
          }
        }
      }
    else //if(gridsize>0) // Gitterpunkte
      {
      for(i=0;i<unsigned(gridsize);i++,fchelpbetap++)
        {
        if(!varcoeff)
          *fchelpbetap = fix * xvalues(i,0) + korrektur;
        else
          *fchelpbetap = fix;
        }
      }
    double help = -korrektur;
    fcconst->update_intercept(help);
    fchelp.update_bootstrap();
    }
  else if(inthemodel==false && fixornot==false)
    {
    double * fchelpbetap = fchelp.getbetapointer();
    for(unsigned i=0;i<fchelp.getbeta().rows();i++,fchelpbetap++)
        *fchelpbetap = 0;
    fchelp.update_bootstrap();
    }
  else
    {
    fchelp.update_bootstrap();
    }
  // FULLCOND::update_bootstrap();     // wie bei Gerade???
  }


void FULLCOND_pspline_stepwise::update_bootstrap_df(void)
  {
  if(optionsp->get_nriter()<=1)
    {
    ST::string path = samplepath.substr(0,samplepath.length()-4)+"_df.raw";
    fc_df = FULLCOND(optionsp,datamatrix(1,1),"title?",1,1,path);
    fc_df.setflags(MCMC::norelchange | MCMC::nooutput);
    }

  if(fixornot==true)
    {
    fc_df.setbetavalue(0,0,-1.0);
    fc_df.update_bootstrap_df();
    }
  else if(inthemodel==false && fixornot==false)
    {
    fc_df.setbetavalue(0,0,0.0);
    fc_df.update_bootstrap_df();
    }
  else // if(inthemodel==true && fixornot==false)
    {
    fc_df.setbetavalue(0,0,lambda);
    fc_df.update_bootstrap_df();
    }
  }


void FULLCOND_pspline_stepwise::save_betamean(void)
  {
  if(fixornot==true)
    {
    bool raus = false;
    unsigned j = 1;
    ST::string name_richtig;
    if(!varcoeff)
      name_richtig = datanames[0];
    else
      name_richtig = datanames[1];  // (VC alt) 0
    while(j<fcconst->get_datanames().size() && raus==false)
      {
      if(fcconst->get_datanames()[j] == name_richtig)
        raus = true;
      j = j + 1;
      }
    unsigned index_fix = j-1;
    double fix = fcconst->getbeta(index_fix,0);
    unsigned i;
    double korrektur = 0;
    if(center)
      korrektur = -0.5*fix*(data_forfixed.max(0) + data_forfixed.min(0));
    double * fchelpbetap = fchelp.getbetapointer();
    if(gridsize < 0)                              // alle verschiedene Beobachtungen
      {
      vector<int>::iterator freqwork = freqoutput.begin();
      int * workindex = index.getV();
      for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
        {
        if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
          {
          if(!varcoeff)
            *fchelpbetap = fix * data_forfixed(*workindex,0) + korrektur;
          else
            *fchelpbetap = fix;
          fchelpbetap++;
          }
        }
      }
    else //if(gridsize>0) // Gitterpunkte
      {
      for(i=0;i<unsigned(gridsize);i++,fchelpbetap++)
        {
        if(!varcoeff)
          *fchelpbetap = fix * xvalues(i,0) + korrektur;
        else
          *fchelpbetap = fix;
        }
      }
    double help = -korrektur;
    fcconst->update_intercept(help);
    fchelp.save_betamean();
    }
  else if(inthemodel==false && fixornot==false)
    {
    double * fchelpbetap = fchelp.getbetapointer();
    for(unsigned i=0;i<fchelp.getbeta().rows();i++,fchelpbetap++)
        *fchelpbetap = 0;
    fchelp.save_betamean();
    }
  else
    {
    fchelp.save_betamean();
    }
  }


void FULLCOND_pspline_stepwise::update_bootstrap_betamean(void)
  {
  fchelp.update_bootstrap_betamean();
  FULLCOND::setflags(MCMC::norelchange);
  }

void FULLCOND_pspline_stepwise::outresults_df(unsigned & size)
  {
  fc_df.update_bootstrap_betamean();
  //fc_df.outresults();
  double * workmean = fc_df.get_betameanp();

  ST::string pathdf = pathcurrent.substr(0,pathcurrent.length()-4)+"_df.res";
  ofstream outres(pathdf.strtochar());

  outres << "df_value   ";
  outres << "sp_value  ";
  outres << "frequency  ";
  outres << "pmean   " << endl;

// Häufigkeitstabelle:

  // samplestream.close();
  datamatrix sample(size,1);
  fc_df.readsample_df(sample,0);
  unsigned i;

//pathdf = pathcurrent.substr(0,pathcurrent.length()-4)+"_df_sample.raw";
//ofstream out(pathdf.strtochar());
//sample.prettyPrint(out);

  vector<unsigned> number;
  vector<unsigned> number1;
  vector<unsigned> number2;
  vector<unsigned> cumnumber1;
  vector<unsigned> cumnumber;

  statmatrix<int> index(sample.rows(),1);
  index.indexinit();
  sample.indexsort(index,0,sample.rows()-1,0,0);

  i = 0;
  unsigned j,anz;
  while(i<index.rows())
     {
     anz=0;
     int* p = index.getV() + i;
     int* q = index.getV() + i;
     j=i;
     while(j<index.rows() && (sample.get(*p,0) == sample.get(*q,0)))
        {
        anz = anz+1;
        j++;
        p++;  
        }
     if(sample.get(*q,0) <= 0)
       number1.push_back(anz);
     else if(sample.get(*q,0) > 0)
       number2.push_back(anz);
     if(cumnumber1.size()>0)
       cumnumber1.push_back(cumnumber1[cumnumber1.size()-1]+anz);
     else
       cumnumber1.push_back(anz);
     i = i + anz;
     }

  int k;
  for(k=number1.size()-1;k>=0;k--)
    {
    cumnumber.push_back(cumnumber1[k]);
    number.push_back(number1[k]);
    }
  for(k=number2.size()-1;k>=0;k--)
    {
    cumnumber.push_back(cumnumber1[k+number1.size()]);
    number.push_back(number2[k]);
    }

  for(i=0;i<number.size();i++)
    {
    double help = sample.get(index(cumnumber[i]-1,0),0);
    double dfs = -1*help;
    if(help!=0 && help!=-1)
      {
      update_stepwise(help);
      set_inthemodel(help);
      dfs = compute_df();
      }
    outres << ST::doubletostring(dfs,6) << "   " << ST::doubletostring(help,6) << "   " << ST::inttostring(number[i]) << "   ";
    if(*workmean == help)
      outres << "selected"; // ST::doubletostring(*workmean,6);
    else
      outres << "-";
    outres << endl;
    }
  }


void FULLCOND_pspline_stepwise::get_samples(const ST::string & filename,const unsigned & step) const
  {
  fchelp.get_samples(filename,step);
  }


void FULLCOND_pspline_stepwise::update(void)
  {
  if(utype != gaussian && optionsp->get_nriter()==1 && interaction==true)
   betaold = beta;

  if(lambda==0)
    {
    beta = datamatrix(beta.rows(),beta.cols(),0);
    FULLCOND::update();
    double * fchelpbetap = fchelp.getbetapointer();
    for(unsigned i=0;i<fchelp.getbeta().rows();i++,fchelpbetap++)
        *fchelpbetap = 0;
    fchelp.update();
    }
  else
    {
    bool interaction_save = interaction;
    bool interaction2 = false;
    if(interactions_pointer.size()>0)
      interaction2 = search_for_interaction();

    bool update = true;
    if(interaction2 == true && interaction == true)
      {
      if(interactions_pointer[interactions_pointer.size()-1]->get_rankK2()==(nrpar-1)*(nrpar-1) &&
        interactions_pointer[interactions_pointer.size()-1]->get_lambda() != 0)
        {
        update = false;
        }
      else if(interactions_pointer[interactions_pointer.size()-1]->get_lambda() == 0)
        interaction = false;
      }
    else if(interaction2 == false && interaction == true)
      interaction = false;

    if(update==true) //(!interaction || !interaction2))  // nur, wenn nicht Teil einer Interaktion
      {
      if(utype == gaussian)
        FULLCOND_pspline_gaussian::update();
      else if(utype != gaussian)
        update_IWLS_mode();

      if(center && intercept!=0.0)
        {
        unsigned i;
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
      intercept = 0.0;
      }
    interaction = interaction_save;
    }
  }


void FULLCOND_pspline_stepwise::change_Korder(double lamb)
  {
  set_lambdaconst(1000000000);
  if(lamb==-1)
    {
    if(!varcoeff)
      {
      if(type==RW1)
        {
        Kenv = Krw2env(weight);
        rankK = nrpar-2;
        }
      }
    else
      {
      if(type==RW2)
        {
        Kenv = Krw1env(weight);
        rankK = nrpar-1;
        }
      }
    }
  else if(lamb==-2)
    {
    if(!varcoeff)
      {
      // MONOTON: fehlt noch!!!
      }
    else
      {
      if(type==RW1)
        {
        Kenv = Krw2env(weight);
        rankK = nrpar-2;
        }
      }
    }
  }

void FULLCOND_pspline_stepwise::undo_Korder(void)
  {
  if(type==RW1 && rankK==nrpar-2)
    {
    Kenv = Krw1env(weight);
    rankK = nrpar-1;
    }
  else if(type==RW2 && rankK==nrpar-1)
    {
    Kenv = Krw2env(weight);
    rankK = nrpar-2;
    }
  }


void FULLCOND_pspline_stepwise::change_varcoeff(const datamatrix & main,const double & inter)
  {
  unsigned i;

  multBS(splinehelp,beta);

  vector<int>::iterator freqwork = freqoutput.begin();
  int * workindex = index.getV();

// splinehelp ändern
  for(i=0;i<splinehelp.rows();i++,freqwork++,workindex++)
    splinehelp(*workindex,0) += main(*freqwork,0);

// Intercept ändern
  intercept += inter;

  double * worksp = spline.getV();
  double * worksph = splinehelp.getV();
  double * workint = data_forfixed.getV();
  for(i=0;i<spline.rows();i++,worksp++,worksph++,workint++)
    {
    *worksp = *worksph * *workint;
    }

// fchelp ändern
  if( (optionsp->get_nriter() > optionsp->get_burnin()) &&      ((optionsp->get_nriter()-optionsp->get_burnin()-1) % (optionsp->get_step()) == 0) )
    {

    double * fchelpbetap = fchelp.getbetapointer();

    freqwork = freqoutput.begin();
    workindex = index.getV();
    for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
      {
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        {
        *fchelpbetap = splinehelp(*workindex,0) - intercept;
        fchelpbetap++;
        }
      }

    }

  fchelp.update();
  //FULLCOND::update();
  }


void FULLCOND_pspline_stepwise::update_IWLS_mode(void)
  {

unsigned updateW = optionsp->get_iterations();
  unsigned i;

  if(lambdaconst == true)
    sigma2 = likep->get_scale(column)/lambda;

  if(optionsp->get_nriter()==1)                             // posterior mode Schätzung übernehmen
    {
    betaold = datamatrix(nrpar,1,0);
    W = datamatrix(likep->get_nrobs(),1,0);
    }
  betaold.assign(beta);

  double logold = likep->loglikelihood(true)
                - 0.5*Kenv.compute_quadformblock(betaold,0,0,nrpar-1)/sigma2;

  add_linearpred_multBS(beta_mode,betaold,true);

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
    likep->compute_IWLS_weight_tildey(W,mu,column,true);
    mu.plus(spline,mu);
    compute_XWXenv_XWtildey(W,1.0);
    }
  else
    {
    likep->tilde_y(mu,spline,column,true,W);
    compute_XWtildey(W,1.0);
    }

  prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);
  prec_env.solve(muy,betahelp);

  double * work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  prec_env.solveU(beta,betahelp);

  add_linearpred_multBS(beta,beta_mode,true);
  beta_mode.assign(betahelp);

  betahelp.minus(beta,beta_mode);
  double qold = - 0.5*prec_env.compute_quadformblock(betahelp,0,0,nrpar-1);

  betahelp.minus(betaold,beta_mode);
  double lognew = likep->loglikelihood(true)
                - 0.5*Kenv.compute_quadformblock(beta,0,0,nrpar-1)/sigma2;
  double qnew = - 0.5*prec_env.compute_quadformblock(betahelp,0,0,nrpar-1);

  double alpha = lognew + qnew - logold - qold;
  double u = log(uniform());

  if(center)
    {
    compute_intercept();
    for(i=0;i<nrpar;i++)
      beta_mode(i,0) -= intercept;
    intercept = 0.0;
    }

  if(u<=alpha)
    {
    acceptance++;
    if(center)
      {
      compute_intercept();
      for(i=0;i<nrpar;i++)
        beta(i,0) -= intercept;
        
      if (varcoeff)
        fcconst->update_fix_varcoeff(intercept,datanames[1]);
      else
        fcconst->update_intercept(intercept);

      for(i=0;i<likep->get_nrobs();i++)
        spline(i,0) -= intercept;
      intercept = 0.0;
      }
    betaold.assign(beta);
    }
  else
    {
    add_linearpred_multBS(betaold,beta,true);
    beta.assign(betaold);
    }
  write_spline();
  fchelp.update();
  FULLCOND::update();
  }


} // end: namespace MCMC














