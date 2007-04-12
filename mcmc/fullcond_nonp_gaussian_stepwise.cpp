
#include "first.h"

#include "fullcond_nonp_gaussian_stepwise.h"

namespace MCMC
{

// additive Effekte, RW1 RW2 und season

FULLCOND_nonp_gaussian_stepwise::FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,
                      DISTRIBUTION * dp,
                      const datamatrix & d,
                      FULLCOND_const * fcc,
                      const unsigned & maxint,const fieldtype & ft,
                      const ST::string & ti,
                      const ST::string & fp, const ST::string & pres,
                      const unsigned & c,const double & l,
                      const unsigned & per)
  : FULLCOND_nonp_gaussian(o,dp,d,fcc,maxint,ft,ti,fp,pres,c,l,per)
  {

  intercept = 0.0;

  beta_average.erase(beta_average.begin(),beta_average.end());

  all_precenv.erase(all_precenv.begin(),all_precenv.end());
  lambdavec.erase(lambdavec.begin(),lambdavec.end());

  if (type == RW1)
    {
    grenzfall = 0;
    dimX = 0;
    dimZ = nrpar-1;
    }
  else if (type == RW2)
    {
    grenzfall = 1;
    dimX = 1;
    dimZ = nrpar-2;
    }
  else if (type == seasonal)
    {
    grenzfall = period - 2;
    dimX = per-1;         // ?
    dimZ = nrpar-per+1;
    }

  spatialtotal = false;
  //gleichwertig = true;
  }

// varying coefficients , RW1 RW2 und season

FULLCOND_nonp_gaussian_stepwise::FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                       const datamatrix & d,
                       const datamatrix & intvar,
                       FULLCOND_const * fcc,
                       const unsigned & maxint,
                       const fieldtype & ft,const ST::string & ti,
                       const ST::string & fp, const ST::string & pres,
                       const unsigned & c,const double & l, const bool & vccent,
                       const unsigned & per)
  : FULLCOND_nonp_gaussian(o,dp,d,intvar,fcc,maxint,ft,ti,fp,pres,c,l,per)

  {

  intercept = 0.0;
  //VCM_neu
  if(vccent == true)
    identifiable = false;

  beta_average.erase(beta_average.begin(),beta_average.end());

  all_precenv.erase(all_precenv.begin(),all_precenv.end());
  lambdavec.erase(lambdavec.begin(),lambdavec.end());

  get_data_forfixedeffects();
  effmodi = d;
  unsigned nrobs = index.rows();

  if (data_varcoeff_fix.rows() < nrobs)
    {
    data_varcoeff_fix = datamatrix(nrobs,2,1);
    for(unsigned i=0;i<nrobs;i++)
      {
      data_varcoeff_fix(i,0) = intvar(i,0);
      data_varcoeff_fix(i,1) = d(i,0)*intvar(i,0);
      }
    }

  XVX = datamatrix(2,2,0);

  if (type == RW1)
    {
    grenzfall = 1;
    dimX = 1;
    dimZ = nrpar-1;
    }
  else if (type == RW2)
    {
    grenzfall = 2;
    dimX = 2;
    dimZ = nrpar-2;
    }
  else if (type == seasonal)
    {
    grenzfall = period - 1;
    dimX = per-1;
    dimZ = nrpar-per+1;
    }

  if(identifiable == false)
    {
    grenzfall -= 1;
    dimX -= 1;
    }

  spatialtotal = false;
  //gleichwertig = true;
  }

// spatial covariates

FULLCOND_nonp_gaussian_stepwise::FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,
                        DISTRIBUTION * dp,const datamatrix & d,
                        FULLCOND_const * fcc,
                        const MAP::map & m, const ST::string & mn,
                        const ST::string & ti,
                        const ST::string & fp, const ST::string & pres,
                        const unsigned & c,const double & l)
  : FULLCOND_nonp_gaussian(o,dp,d,fcc,m,mn,ti,fp,pres,c,l)

  {

  intercept = 0.0;
  beta_average.erase(beta_average.begin(),beta_average.end());

  all_precenv.erase(all_precenv.begin(),all_precenv.end());
  lambdavec.erase(lambdavec.begin(),lambdavec.end());

  grenzfall = 0;

  spatialtotal = false;
  //gleichwertig = true;

  dimX = 0;
  dimZ = rankK;
  }

// varying coefficients , spatial covariates as effect modifier

FULLCOND_nonp_gaussian_stepwise::FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,
                        DISTRIBUTION * dp,
                        FULLCOND_const * fcc,
                        const MAP::map & m,
                        const ST::string & mn,
                        const datamatrix & d,
                        const datamatrix & d2,
                        const ST::string & ti,
                        const ST::string & fp, const ST::string & pres,
                        const unsigned & c, const double & l, const bool & vccent)
  : FULLCOND_nonp_gaussian(o,dp,fcc,m,mn,d,d2,ti,fp,pres,c,l)

  {

  intercept = 0.0;
  //VCM_neu
  if(vccent == true)
    identifiable = false;

  get_data_forfixedeffects();
  beta_average.erase(beta_average.begin(),beta_average.end());

  all_precenv.erase(all_precenv.begin(),all_precenv.end());
  lambdavec.erase(lambdavec.begin(),lambdavec.end());

  grenzfall = 1;

  dimX = 1;
  dimZ = rankK;

  if(identifiable == false)
    {
    grenzfall -= 1;
    dimX -= 1;
    }

  spatialtotal = false;
  //gleichwertig = true;
  }


 void FULLCOND_nonp_gaussian_stepwise::init_spatialtotal(FULLCOND * unstructp)
  {
  spatialtotal = true;
  fcunstruct = unstructp;
  }

  // COPY CONSTRUCTOR

FULLCOND_nonp_gaussian_stepwise::FULLCOND_nonp_gaussian_stepwise(const FULLCOND_nonp_gaussian_stepwise & fc)
  : FULLCOND_nonp_gaussian(FULLCOND_nonp_gaussian(fc))
  {
  intercept = fc.intercept;
  beta_average = fc.beta_average;
  //lambda_nr = fc.lambda_nr;
  //lambdas_local = fc.lambdas_local;
  data_varcoeff_fix = fc.data_varcoeff_fix;
  effmodi = fc.effmodi;
  XVX = fc.XVX;
  fcunstruct = fc.fcunstruct;
  spatialtotal = fc.spatialtotal;
  //gleichwertig = fc.gleichwertig;
  df_lambdaold = fc.df_lambdaold;
  lambdaold = fc.lambdaold;
  df_lambdaold_unstr = fc.df_lambdaold_unstr;
  lambdaold_unstr = fc.lambdaold_unstr;
  lambdavec = fc.lambdavec;
  all_precenv = fc.all_precenv;
  fc_df = fc.fc_df;
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const FULLCOND_nonp_gaussian_stepwise & FULLCOND_nonp_gaussian_stepwise::operator=(
                                            const FULLCOND_nonp_gaussian_stepwise & fc)
  {
  if (this == &fc)
    return *this;
  FULLCOND_nonp_gaussian::operator=(FULLCOND_nonp_gaussian(fc));

  intercept = fc.intercept;
  beta_average = fc.beta_average;
  //lambda_nr = fc.lambda_nr;
  //lambdas_local = fc.lambdas_local;
  data_varcoeff_fix = fc.data_varcoeff_fix;
  effmodi = fc.effmodi;
  XVX = fc.XVX;
  fcunstruct = fc.fcunstruct;
  spatialtotal = fc.spatialtotal;
  //gleichwertig = fc.gleichwertig;
  df_lambdaold = fc.df_lambdaold;
  lambdaold = fc.lambdaold;
  df_lambdaold_unstr = fc.df_lambdaold_unstr;
  lambdaold_unstr = fc.lambdaold_unstr;
  lambdavec = fc.lambdavec;
  all_precenv = fc.all_precenv;
  fc_df = fc.fc_df;  
      
  return *this;
  }


bool FULLCOND_nonp_gaussian_stepwise::posteriormode(void)
  {
  int j;
  unsigned i;

  int * workindex;

  update_linpred(false);

  // NEU!!!
  if(varcoeff && lambda == -2)
    {
    //datamatrix X = datamatrix(2,2,0);
    datamatrix betas = datamatrix(2,1,0);

    if(calculate_xwx_vc == true || (XVX(1,1)==0 && XVX(0,0)==0))
      {
      calculate_xwx_vc = false;
      likep->fisher(XVX,data_varcoeff_fix,column);            // recomputes X1 = (newx' W newx)^{-1}
      XVX.assign((XVX.cinverse()));               // continued
      }
    likep->compute_weightiwls_workingresiduals(column); // computes W(y-linpred)
    betas = XVX*data_varcoeff_fix.transposed()*likep->get_workingresiduals();

    datamatrix spline = datamatrix(likep->get_nrobs(),1,0);
    spline.mult(data_varcoeff_fix,betas);
    double * workbeta = beta.getV();
    vector<int>::iterator itbeg = posbeg.begin();
    vector<int>::iterator itend = posend.begin();
    int * workindex = index.getV();
    for(i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        {
        *workbeta = betas(0,0) + effmodi(*workindex,0)*betas(1,0);
        for(j=*itbeg;j<=*itend;j++)
          workindex++;
        }
      }
    likep->add_linearpred_m(spline,column);     // addiert durchschnittl. Fkt. zum Gesamtprädiktor

    if(center)
      {
      intercept = centerbeta();
      update_fix_effect(intercept);
      intercept = 0.0;
      }
    }
  else
    {
    if (lambda_prec != lambda || calculate_xwx == true)
      {
      if (calculate_xwx == true)
        {
        calculate_xwx = false;
        if (varcoeff)
          compute_XWX_varcoeff_env(likep->get_weightiwls(),column);
        else
          compute_XWX_env(likep->get_weightiwls(),column);
        }
      precenv.addtodiag(XXenv,Kenv,1.0,lambda);
      lambda_prec = lambda;
      }

    likep->compute_weightiwls_workingresiduals(column);

    workindex = index.getV();
    double * workmuy = beta.getV();

    if (varcoeff)
      {
      double * workdata=data.getV();
      for(i=0;i<nrpar;i++,workmuy++)
        {
        *workmuy = 0;
        if (posbeg[i] != -1)
          for(j=posbeg[i];j<=posend[i];j++,workindex++,workdata++)
            *workmuy+=
            likep->get_workingresiduals()(*workindex,0)*(*workdata);
        }
      }
    else  // else additive
      {
      for(i=0;i<nrpar;i++,workmuy++)
        {
        *workmuy = 0;
        if (posbeg[i] != -1)
          for(j=posbeg[i];j<=posend[i];j++,workindex++)
            *workmuy+= likep->get_workingresiduals()(*workindex,0);
        }
      }

    precenv.solve(beta);
    update_linpred(true);

    if (center)
      {
      intercept = centerbeta();
      if(varcoeff == false)
        {
        fcconst->posteriormode_intercept(intercept);
        }
      else if(varcoeff == true)
        {
        update_fix_effect(intercept);
        intercept = 0.0;
        }
      }
    } // END: else if(!varcoeff || lambda != -2)

  transform = likep->get_trmult(column);

  return FULLCOND_nonp_basis::posteriormode();

  }


const datamatrix & FULLCOND_nonp_gaussian_stepwise::get_data_forfixedeffects(void)
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


void FULLCOND_nonp_gaussian_stepwise::update_stepwise(double la)
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
      precenv = all_precenv[i-1];
      lambda_prec = lambda;
      }
    }
  }


double FULLCOND_nonp_gaussian_stepwise::compute_df(void)
  {
  double df = 0;
  //if(inthemodel == false && fixornot == false)
  //  {
  //  if(spatialtotal)
  //    fcunstruct->compute_df_andererteil();
  //  }
  //if(inthemodel == false && fixornot == true)
  //  df = 1;
  if(inthemodel == true)
    {
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
      }  */

    if(varcoeff && lambda == -2)
      {
      if(identifiable)
        df = 2;
      else
        df = df + 1;
      }
    else
      {
      bool unstr_included = false;
      if(spatialtotal)
        {
        bool fix;
        fcunstruct->get_inthemodel(unstr_included,fix);
        if(unstr_included == true)
          {
          double lambda_unstr = fcunstruct->get_lambda();
          double df_str = 0;
          double df_unstr = 0;

          if (lambdaold == lambda && lambdaold_unstr == lambda_unstr && likep->get_iwlsweights_notchanged() == true)
            {
            df = df_lambdaold;
            fcunstruct->set_dfunstruct(df_lambdaold_unstr);
            }
          else
            {
            if (calculate_xwx == true)
              {
              calculate_xwx = false;
              compute_XWX_env(likep->get_weightiwls(),column);
              }

            envmatdouble Diag_neu = envmatdouble(0,nrpar);
            vector<double>::iterator d = XXenv.getDiagIterator();
            vector<double>::iterator d2 = Diag_neu.getDiagIterator();
            unsigned i;
            for(i=0;i<nrpar;i++,d++,d2++)
              *d2 = *d * lambda_unstr / (*d + lambda_unstr);
            envmatdouble precenv_neu = envmatdouble(Kenv.getXenv(),0,nrpar);
            precenv_neu.addtodiag(Diag_neu,Kenv,1.0,lambda);
            invprec = envmatdouble(precenv_neu.getXenv(),0,precenv_neu.getDim());
            precenv_neu.inverse_envelope(invprec);

            df_str += invprec.traceOfProduct(XXenv);

            d = XXenv.getDiagIterator();
            d2 = Diag_neu.getDiagIterator();
            for(i=0;i<nrpar;i++,d++,d2++)
              *d2 = *d * *d / (*d + lambda_unstr);

            df_unstr -= invprec.traceOfProduct(Diag_neu);
            df_str += df_unstr;

            d = XXenv.getDiagIterator();
            d2 = Diag_neu.getDiagIterator();
            for(i=0;i<nrpar;i++,d++,d2++)
              *d2 = *d2 * *d / (*d + lambda_unstr);

            df_unstr += invprec.traceOfProduct(Diag_neu);

            d = XXenv.getDiagIterator();
            for(i=0;i<nrpar;i++,d++,d2++)
              df_unstr += *d / (*d + lambda_unstr);

          //if(!gleichwertig)
          //  {
          //  df = df_str-1+df_unstr;
          //  }
          //else
          //  {
            fcunstruct->set_dfunstruct(df_unstr);
            df = df_str-1;
            df_lambdaold = df;
            df_lambdaold_unstr = df_unstr;
            lambdaold = lambda;
            lambdaold_unstr = lambda_unstr;
          //  }
            }
          }
        }

      if(!spatialtotal || !unstr_included)
        {

        if (lambdaold == lambda && likep->get_iwlsweights_notchanged() == true && !spatialtotal)
          {
          df = df_lambdaold;
          }
        else if (spatialtotal || lambdaold != lambda || likep->get_iwlsweights_notchanged() == false)
          {
          if (calculate_xwx == true)
            {
            if (varcoeff)
              compute_XWX_varcoeff_env(likep->get_weightiwls(),column);
            else
              compute_XWX_env(likep->get_weightiwls(),column);
            }
          if(lambda != lambda_prec || calculate_xwx == true)
            {
            calculate_xwx = false;
            precenv.addtodiag(XXenv,Kenv,1.0,lambda);
            lambda_prec = lambda;
            }

          if (type==MCMC::mrf)
            invprec = envmatdouble(precenv.getXenv(),0,precenv.getDim());
          else
            invprec = envmatdouble(0,nrpar,Kenv.getBandwidth());

          precenv.inverse_envelope(invprec);

          if(identifiable)
            df = invprec.traceOfProduct(XXenv);
          else
            df = df + invprec.traceOfProduct(XXenv)-1;

          df_lambdaold = df;
          lambdaold = lambda;
          }
        }
      }
    }

  return df;
  }


/*double FULLCOND_nonp_gaussian_stepwise::compute_df_andererteil(void)    // für spatialtotal
  {
  double df = 0;
  if(spatialtotal && !gleichwertig)
    {
    bool fix,unstr_included;
    fcunstruct->get_inthemodel(unstr_included,fix);
    if(unstr_included == true)
      df = fcunstruct->compute_df();
    }

  return df;
  } */

/*void FULLCOND_nonp_gaussian_stepwise::set_gleichwertig(const bool & gleich, bool weiter)       // für spatialtotal
  {
  if(spatialtotal)
    {
    gleichwertig = gleich;
    if(weiter == true)
      fcunstruct->set_gleichwertig(gleich,false);
    }
  } */


ST::string FULLCOND_nonp_gaussian_stepwise::get_effect(void)
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
    h = datanames[0] + "*" + datanames[1];
  else
    h = datanames[0];

  h = h + "(" + t + ",df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";

  return h;
  }


ST::string FULLCOND_nonp_gaussian_stepwise::get_befehl(void)
  {
  ST::string h;

  ST::string t;
  if (type==MCMC::RW1)
    t = "rw1";
  else if (type==MCMC::RW2)
    t = "rw2";
  else if (type==MCMC::seasonal)
    t = "season,period=" + ST::inttostring(period);
  else if (type==MCMC::mrf)
    t = "spatial,map=" + mapname;

  if(varcoeff)
    h = datanames[0] + "*" + datanames[1];
  else
    h = datanames[0];

  h = h + "(" + t + ",lambda=" + ST::doubletostring(lambda,6) + ")";

  return h;
  }


void FULLCOND_nonp_gaussian_stepwise::init_names(const vector<ST::string> & na)
    {
    FULLCOND::init_names(na);

    char hchar = '_';
    ST::string hstring =  "\\_";

    if (na.size()==1)
      {
      ST::string helpname = na[0].insert_string_char(hchar,hstring);
      if (type==MCMC::seasonal)
        term_symbolic = "f^{Season}_{" +  helpname + "}("+helpname+")";
      else
        term_symbolic = "f_{" +  helpname + "}("+helpname+")";
      }
    else
      {
      ST::string helpname1 = na[1].insert_string_char(hchar,hstring);
      ST::string helpname2 = na[0].insert_string_char(hchar,hstring);
      if (type==MCMC::seasonal)
        term_symbolic = "f^{Season}_{" +  helpname1 + "}("+helpname1+") \\cdot "
                        + helpname2;
      else
        term_symbolic = "f_{" +  helpname1 + "}("+helpname1+") \\cdot "
                        + helpname2;
      }

    if (column > 0)
      priorassumptions.push_back("$" + term_symbolic + "$" +
       " (" + ST::inttostring(column+1) + ". \\mbox{ } response \\mbox{ } category)");
    else
      priorassumptions.push_back("$" + term_symbolic + "$");
//    priorassumptions.push_back(term_symbolic);
    init_priorassumptions(na[0]);
    }


void FULLCOND_nonp_gaussian_stepwise::reset_effect(const unsigned & pos)
  {

  update_linpred(false);

  unsigned i;
  double * work;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;

  intercept = 0.0;
  }


void FULLCOND_nonp_gaussian_stepwise::hierarchie_rw1(vector<double> & untervector, int dfo)
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


void FULLCOND_nonp_gaussian_stepwise::compute_lambdavec(
vector<double> & lvec, int & number)
  {
  lambdaold = -1;
  if (df_equidist==true && spfromdf==true)
     FULLCOND::compute_lambdavec_equi(lvec,number);
  else
     FULLCOND::compute_lambdavec(lvec,number);

  if(likep->iwlsweights_constant() == true)
    {
    lambdavec = lvec;
    if (varcoeff)
      compute_XWX_varcoeff_env(likep->get_weightiwls(),column);
    else
      compute_XWX_env(likep->get_weightiwls(),column);
    for(unsigned i=0;i<lambdavec.size();i++)
      {
      precenv.addtodiag(XXenv,Kenv,1.0,lambdavec[i]);
      precenv.decomp();
      all_precenv.push_back(precenv);
      }
    }

  if (!nofixed && (type==RW1) && (!varcoeff) )
    {
    hierarchie_rw1(lvec,1);
    }
  else if (!nofixed && (type==RW1) && (varcoeff) )
    {
    if(identifiable)
      {
      hierarchie_rw1(lvec,2);
      lvec.push_back(-1);
      }
    else
      hierarchie_rw1(lvec,1);
    }
  else if (!nofixed && (type==RW2) && (!varcoeff) )
    {
    lvec.push_back(-1);
    }
  else if (!nofixed && (type==RW2) && (varcoeff) )
    {
    lvec.push_back(-2);
    if(identifiable)    //VCM_neu
      lvec.push_back(-1);
    }
  else if ( (type==mrf) && (varcoeff) )
    {
    if(!nofixed && identifiable)
      lvec.push_back(-1);
    }

  if(forced_into==false)
     lvec.push_back(0);

  // Startwert für lambda aus df:
  if(spfromdf==true)
    {
    double lambdavorg = 1000;
    if(!varcoeff)
      {
      if(!nofixed && dfstart==1 && (type==RW1 || type==RW2))
        lambdastart = -1;
      else if(dfstart==0)
        lambdastart = 0;
      else
        lambdastart = lambda_from_df(dfstart,lambdavorg);
      }
    else
      {
      if(!nofixed && dfstart==1 && identifiable)
        lambdastart = -1;
      else if(!nofixed && ((dfstart==2 && identifiable) || (dfstart==1 && !identifiable)) && (type==RW1 || type==RW2))
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


// BEGIN: For Varying Coefficients ---------------------------------------------

void FULLCOND_nonp_gaussian_stepwise::update_fix_effect(double & intercept)
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
    interactions_pointer[0]->set_inthemodel(-1);
    fcconst->update_fix_effect(j,intercept,data_forfixed);
    }
  }


void FULLCOND_nonp_gaussian_stepwise::set_pointer_to_interaction(FULLCOND * inter)
  {
  interactions_pointer.push_back(inter);
  }

void FULLCOND_nonp_gaussian_stepwise::get_interactionspointer(vector<FULLCOND*> & inter)
  {
  inter = interactions_pointer;
  }


void FULLCOND_nonp_gaussian_stepwise::hierarchical(ST::string & possible)
  {
  if(!varcoeff)
    {
    possible = "alles";
    }
  else
    {
    possible = "valles";
    }
  }


void FULLCOND_nonp_gaussian_stepwise::const_varcoeff(void)
  {
  if(varcoeff)
    fcconst->posteriormode_const_varcoeff(data_forfixed);
  }

// END: For Varying Coefficients -----------------------------------------------


void FULLCOND_nonp_gaussian_stepwise::create_weight(datamatrix & w)
  {
  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();
  int * workindex = index.getV();
  unsigned i;

  if(type != mrf)
    {
    w(*workindex,0) = 1;

    itbeg = posbeg.begin() + posbeg.size()-1;
    workindex = index.getV() + *itbeg;
    w(*workindex,0) = 1;
    }
  else
    {
    int j;
    for(i=0;i<nrpar;i++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        {
        w(*workindex,0) = 1;
        for(j=*itbeg;j<=*itend;j++)
          workindex++;
        }
      }
    }
  }


// BEGIN: MODEL-AVERAGING ------------------------------------------------------

void FULLCOND_nonp_gaussian_stepwise::update_bootstrap(const bool & uncond)
  {
  if(optionsp->get_samplesize()==1)
    {
    ST::string path = samplepath.substr(0,samplepath.length()-4)+"_df.raw";
    fc_df = FULLCOND(optionsp,datamatrix(1,1),"title?",1,1,path);
    fc_df.setflags(MCMC::norelchange | MCMC::nooutput);
    }

  datamatrix betaold = beta;

  if(fixornot==true)
    {
    bool raus = false;
    unsigned j = 1;
    ST::string name_richtig = datanames[0];
    while(j<fcconst->get_datanames().size() && raus==false)
      {
      if(fcconst->get_datanames()[j] == datanames[0])
        raus = true;
      j = j + 1;
      }
    unsigned index_fix = j-1;
    double fix = fcconst->getbeta(index_fix,0);
    unsigned i;
    double * workbeta = beta.getV();
    vector<int>::iterator itbeg = posbeg.begin();
    vector<int>::iterator itend = posend.begin();
    int * workindex = index.getV();
    double sum = 0;
    int k;
    for(i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        {
        if(!varcoeff)
          {
          *workbeta = data_forfixed(*workindex,0)*fix;
          sum += *workbeta;
          }
        else
          *workbeta = fix;
        for(k=*itbeg;k<=*itend;k++)
          workindex++;
        }
      }
    workbeta = beta.getV();
    sum /= double(nrpar);
    if(!center)
      sum = 0;
    for(i=0;i<nrpar;i++,workbeta++)
      *workbeta -= sum;
    double chelp = sum*double(nrpar);
    fcconst->update_intercept(chelp);

    FULLCOND::update_bootstrap();
    //fc_df.setbetavalue(0,0,1);
    fc_df.setbetavalue(0,0,-1);
    fc_df.update_bootstrap();
    }
  else if(inthemodel==false && fixornot==false)
    {
    beta = datamatrix(nrpar,1,0);
    FULLCOND::update_bootstrap();
    fc_df.setbetavalue(0,0,0);
    fc_df.update_bootstrap();
    }
  else
    {
    FULLCOND::update_bootstrap();
    //fc_df.setbetavalue(0,0,compute_df());
    fc_df.setbetavalue(0,0,lambda);
    fc_df.update_bootstrap();
    }
  beta = betaold;
  }


void FULLCOND_nonp_gaussian_stepwise::update_bootstrap_betamean(void)
  {
  FULLCOND::update_bootstrap_betamean();
  FULLCOND::setflags(MCMC::norelchange);

  fc_df.update_bootstrap_betamean();
  //fc_df.outresults();
  double * workmean = fc_df.get_betameanp();

  ST::string pathdf = pathcurrent.substr(0,pathcurrent.length()-4)+"_df.res";
  ofstream outres(pathdf.strtochar());

  outres << "value   ";
  outres << "frequency  ";
  outres << "pmean   " << endl;

// Häufigkeitstabelle:

  //samplestream.close();
  datamatrix sample(optionsp->get_samplesize(),1);
  fc_df.readsample(sample,0);
  unsigned i;

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
    if(help>0)
      {
      update_stepwise(help);
      set_inthemodel(help);
      dfs = compute_df();
      }
    outres << ST::doubletostring(dfs,6) << "   " << ST::inttostring(number[i]) << "   ";
    if(*workmean == help)
      outres << "selected"; // ST::doubletostring(*workmean,6);
    else
      outres << "-";
    outres << endl;
    }
  }


void FULLCOND_nonp_gaussian_stepwise::createreml(datamatrix & X,datamatrix & Z,
                                  const unsigned & Xpos, const unsigned & Zpos)
  {
/*  unsigned i,k;
  int j;
  if(!varcoeff)
    {
    if(fctype==MCMC::spatial)
      {
      if(remlspatialdesign.rows() < index.rows())
        {
        ST::string test = "hallo";
       /* datamatrix Kstat=STATMAT_PENALTY::Kmrf(m);    // wie???
        datamatrix vals = datamatrix(Kstat.rows(),1,0);

        bool eigentest=eigen2(Kstat,vals);
        if(eigentest==false)
          {    // andere Möglichkeit?!
          errors.push_back("ERROR: Unable to compute eigen decomposition for structured spatial effect.\n");
          }
        else
          {
          eigensort(vals,Kstat);

          for(i=0; i<vals.rows()-1; i++)
            {
            vals(i,0)=1/sqrt(vals(i,0));
            }
          vals(vals.rows()-1,0)=0;
          remlspatialdesign = multdiagback(Kstat,vals).getColBlock(0,Kstat.cols()-1);
          } */
/*        }
      for(i=0;i<posbeg.size();i++)
        {
        for (j=posbeg[i];j<=posend[i];j++)
          {
          if(j!=-1)
            {
            for(k=0; k<remlspatialdesign.cols(); k++)
              {
              Z(index(j,0),Zpos+k)=remlspatialdesign(i,k);
              }
            }
          }
        }
      }
    else
      FULLCOND_nonp_gaussian::createreml(X,Z,Xpos,Zpos);
    }
  else
    {
    if(fctype==MCMC::spatial)
      {
      for(i=0; i<X.rows(); i++)
        {
        X(i,Xpos)=data_forfixed(i,0);
        }
      for(i=0;i<posbeg.size();i++)
        {
        for (j=posbeg[i];j<=posend[i];j++)
          {
          if(j!=-1)
            {
            for(k=0; k<remlspatialdesign.cols(); k++)
              {
              Z(index(j,0),Zpos+k)=remlspatialdesign(i,k)*data_forfixed(index(j,0),0);
              }
            }
          }
        }
      }
    else
      FULLCOND_nonp_gaussian::createreml(X,Z,Xpos,Zpos);
    }*/
  }


} // end: namespace MCMC




