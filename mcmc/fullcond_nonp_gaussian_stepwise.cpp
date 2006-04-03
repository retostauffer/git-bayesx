
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
  interactions_pointer.erase(interactions_pointer.begin(),interactions_pointer.end());

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
  }

// varying coefficients , RW1 RW2 und season

FULLCOND_nonp_gaussian_stepwise::FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                       const datamatrix & d,
                       const datamatrix & intvar,
                       FULLCOND_const * fcc,
                       const unsigned & maxint,
                       const fieldtype & ft,const ST::string & ti,
                       const ST::string & fp, const ST::string & pres,
                       const unsigned & c,const double & l, const bool & nofixed,
                       const unsigned & per)
  : FULLCOND_nonp_gaussian(o,dp,d,intvar,fcc,maxint,ft,ti,fp,pres,c,l,per)

  {

  intercept = 0.0;
  //VCM_neu
  if(nofixed == true)
    identifiable = false;

  interactions_pointer.erase(interactions_pointer.begin(),interactions_pointer.end());
  beta_average.erase(beta_average.begin(),beta_average.end());

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
  interactions_pointer.erase(interactions_pointer.begin(),interactions_pointer.end());

  beta_average.erase(beta_average.begin(),beta_average.end());

  grenzfall = 0;

  spatialtotal = false;

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
                        const unsigned & c, const double & l, const bool & nofixed)
  : FULLCOND_nonp_gaussian(o,dp,fcc,m,mn,d,d2,ti,fp,pres,c,l)

  {

  intercept = 0.0;
  //VCM_neu
  if(nofixed == true)
    identifiable = false;

  interactions_pointer.erase(interactions_pointer.begin(),interactions_pointer.end());
  beta_average.erase(beta_average.begin(),beta_average.end());

  grenzfall = 1;

  dimX = 1;
  dimZ = rankK;

  if(identifiable == false)
    {
    grenzfall -= 1;
    dimX -= 1;
    }

  spatialtotal = false;
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
  interactions_pointer = fc.interactions_pointer;
  //lambda_nr = fc.lambda_nr;
  //lambdas_local = fc.lambdas_local;
  data_varcoeff_fix = fc.data_varcoeff_fix;
  effmodi = fc.effmodi;
  fcunstruct = fc.fcunstruct;
  spatialtotal = fc.spatialtotal;
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
  interactions_pointer = fc.interactions_pointer;
  //lambda_nr = fc.lambda_nr;
  //lambdas_local = fc.lambdas_local;
  data_varcoeff_fix = fc.data_varcoeff_fix;
  effmodi = fc.effmodi;
  fcunstruct = fc.fcunstruct;
  spatialtotal = fc.spatialtotal;

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
    datamatrix X = datamatrix(2,2,0);
    datamatrix betas = datamatrix(2,1,0);

    likep->fisher(X,data_varcoeff_fix,column);            // recomputes X1 = (newx' W newx)^{-1}
    X.assign((X.cinverse()));               // continued
    likep->compute_weightiwls_workingresiduals(column); // computes W(y-linpred)
    betas = X*data_varcoeff_fix.transposed()*likep->get_workingresiduals();

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


double FULLCOND_nonp_gaussian_stepwise::compute_df(void)
  {

  double df = 0;
  if(varcoeff && !identifiable && !center)
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
    }

  if(varcoeff && lambda == -2)
    {
    if(identifiable)
      return 2;
      //return 1;
    else
      return df + 1;
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
        double lambda_unstr = 0;
        double df_str = 0;
        double df_unstr = 0;
        fcunstruct->get_lambda(lambda_unstr);
        compute_XWX_env(likep->get_weightiwls(),column);

        envmatdouble Diag_neu = envmatdouble(0,nrpar);
        vector<double>::iterator d = XXenv.getDiagIterator();
        vector<double>::iterator d2 = Diag_neu.getDiagIterator();
        unsigned i;
        for(i=0;i<nrpar;i++,d++,d2++)
          *d2 = *d * lambda_unstr / (*d + lambda_unstr);
        precenv.addtodiag(Diag_neu,Kenv,1.0,lambda);
        invprec = envmatdouble(precenv.getXenv(),0,precenv.getDim());
        precenv.inverse_envelope(invprec);

        df_str += invprec.traceOfProduct(XXenv);

        d = XXenv.getDiagIterator();
        d2 = Diag_neu.getDiagIterator();
        for(i=0;i<nrpar;i++,d++,d2++)
          *d2 = *d * *d / (*d + lambda_unstr);

        df_str -= invprec.traceOfProduct(Diag_neu);
        df_unstr -= invprec.traceOfProduct(Diag_neu);

        d = XXenv.getDiagIterator();
        d2 = Diag_neu.getDiagIterator();
        for(i=0;i<nrpar;i++,d++,d2++)
          *d2 = *d2 * *d / (*d + lambda_unstr);

        df_unstr += invprec.traceOfProduct(Diag_neu);

        d = XXenv.getDiagIterator();
        for(i=0;i<nrpar;i++,d++,d2++)
          df_unstr += *d / (*d + lambda_unstr);

        if(identifiable)
          return df_str;
        else
          return df_str-1;
        }
      }

    if(!spatialtotal || !unstr_included)
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

      if(identifiable)
        return invprec.traceOfProduct(XXenv);
      else
        return df + invprec.traceOfProduct(XXenv)-1;
      }
    }
  }


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
    h = datanames[1];
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
    h = datanames[1];
  else
    h = datanames[0];

  h = h + "(" + t + ",lambda=" + ST::doubletostring(lambda,6) + ")";

  return h;
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
  if (get_df_equidist()==true)
     FULLCOND::compute_lambdavec_equi(lvec,number);
  else
     FULLCOND::compute_lambdavec(lvec,number);

  if ( (type==RW1) && (!varcoeff) )
    {
    hierarchie_rw1(lvec,1);
    }
  else if ( (type==RW1) && (varcoeff) )
    {
    if(identifiable)
      {
      hierarchie_rw1(lvec,2);
      lvec.push_back(-1);
      }
    else
      hierarchie_rw1(lvec,1);
    }
  else if ( (type==RW2) && (!varcoeff) )
    {
    lvec.push_back(-1);
    }
  else if ( (type==RW2) && (varcoeff) )
    {
    lvec.push_back(-2);
    if(identifiable)    //VCM_neu
      lvec.push_back(-1);
    }
  else if ( (type==mrf) && (varcoeff) )
    {
    if(identifiable)
      lvec.push_back(-1);
    }


  get_forced();
  if(forced_into==false)
     lvec.push_back(0);
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
    fcconst->update_fix_effect(j,intercept,data_forfixed);
    }
  }


void FULLCOND_nonp_gaussian_stepwise::set_pointer_to_interaction(FULLCOND * inter)
  {
  interactions_pointer.push_back(inter);
  }


void FULLCOND_nonp_gaussian_stepwise::hierarchical(ST::string & possible)
  {
  unsigned i;
  bool spline = false;
  //bool fix = false;
  bool spline1, fix1;
  if(!varcoeff)
    {
    for(i=0;i<interactions_pointer.size();i++)
      {
      interactions_pointer[i]->get_inthemodel(spline1,fix1);
      if(spline1 == true)
        spline = true;
      //if(fix1 == true)
      //  fix = true;
      }

    /*if(interaction)
      {
      if(spline == true)
        possible = "spline";
      else if(fix == true && spline == false)
        possible = "spfix";
      else
        possible = "alles";
      }*/
    //else     // VC
    //  {
      if(spline == true)
        possible = "vfix";
      else
        possible = "alles";
    //  }
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

// BEGIN: MODEL-AVERAGING ------------------------------------------------------

void FULLCOND_nonp_gaussian_stepwise::save_betas(vector<double> & modell, int & anzahl)
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


void FULLCOND_nonp_gaussian_stepwise::average_posteriormode(vector<double> & crit_weights)
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
    {
    pmean_fix.mult(data_forfixed,beta_fixx);     // berechnet den Anteil der fixen Effekte

    if(!varcoeff)      // Zentrieren
      {
      workbeta = beta.getV();
      itbeg = posbeg.begin();
      itend = posend.begin();
      int * workindex = index.getV();
      for(i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
        {
        if(*itbeg != -1)
          {
          *workbeta = pmean_fix(*workindex,0);
          for(j=*itbeg;j<=*itend;j++)
            workindex++;
          }
        }

      double intercept = centerbeta();
      fcconst->set_intercept_for_center(intercept);
      datamatrix inter = datamatrix(likep->get_nrobs(),1,-intercept);
      pmean_fix.plus(pmean_fix,inter);         // zentrierte durchschnittliche Fkt, Einträge in Reihenfolge x_1,...,x_n
      }
    }

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

  likep->add_linearpred_m(pmean_spline,column);     // addiert durchschnittl. Fkt. zum Gesamtprädiktor

  workbeta = beta.getV();
  workeff = betamean.getV();
  for(i=0;i<nrpar;i++,workbeta++,workeff++)
    *workeff = *workbeta * transform;
  }


void FULLCOND_nonp_gaussian_stepwise::effect_sort(datamatrix & effect, const double & m,
            const unsigned & beg, const unsigned & end,const statmatrix<int> & index)
  {
  unsigned register i;
  int * workindex = index.getV() + beg;
  for (i=beg;i<=end;i++,workindex++)
    effect(*workindex,0)+=m;
  }

// Vorschlag:
//void FULLCOND_nonp_gaussian_stepwise::effect_sort(datamatrix & effect, const double & m, unsigned & row)
void FULLCOND_nonp_gaussian_stepwise::effect_sort(datamatrix & effect, double m, unsigned row)
  {
  double * workl = effect.getV() + row;
  *workl += m;
  }


 
} // end: namespace MCMC


