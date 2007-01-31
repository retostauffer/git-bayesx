
#include "first.h"

#include "fullcond_merror.h"

//------------------------------------------------------------------------------
//------------ CLASS: FULLCOND_merror implementation of member functions -------
//------------------------------------------------------------------------------


namespace MCMC
{

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // CONSTRUCTOR (Susi)
  // o    : pointer to MCMCoptions object
  // t    : title of the full conditional (for example "fixed effects")
  //        (i.e. number of categories of the response variable)
  // fp   : file path for storing sampled parameters

  fullcond_merror::fullcond_merror(MCMCoptions * o, FULLCOND_nonp_gaussian * p,
           DISTRIBUTION * dp, const datamatrix & d, const ST::string & t, const ST::string & fp)
                         : FULLCOND(o,d,t,d.rows(),d.cols(),fp)
    {
    designp = p;
    likep = dp;

    varcoeff=true;
    unsigned i, j;
    meandata = datamatrix(d.rows(),1,0);
    unsigned dcols = d.cols();
    for(i=0; i<d.rows(); i++)
      {
      for(j=0; j<dcols; j++)
        {
        meandata(i,0) = d(i,j);
        }
      meandata(i,0) /= dcols;
      }
    setbeta(meandata);
    }

// BEGIN: merror
  // CONSTRUCTOR (Thomas)
  fullcond_merror::fullcond_merror(MCMCoptions * o, spline_basis * p,
           DISTRIBUTION * dp, const datamatrix & d, const ST::string & t,
           const ST::string & fp, const ST::string & pres, const double & lk,
           const double & uk)
           : FULLCOND(o,d,t,d.rows(),1,fp)
  {
  splinep = p;
  likep= dp;
  varcoeff=false;
  unsigned i, j;
  data = d;
  meandata = datamatrix(d.rows(),1,0);
  merror = d.cols();
  for(i=0; i<d.rows(); i++)
    {
    for(j=0; j<merror; j++)
      {
      meandata(i,0) += d(i,j);
      }
    meandata(i,0) /= merror;
    }
  minx = lk;
  maxx = uk;
  setbeta(meandata);

  logfcold = datamatrix(d.rows(),1,0);
  logfcnew = datamatrix(d.rows(),1,0);

  generrcount=0;
  generrtrial=0;

  pathresults = pres;
  ST::string path = pathresults.substr(0,pathresults.length()-4);

  fc_merrorvar = FULLCOND(o,datamatrix(1,1,1.0),title+"_merror_var",1,1,path+"_merror_var.res");
  fc_merrorvar.setflags(MCMC::norelchange | MCMC::nooutput);
//  double * merrorvarp = fc_merrorvar.getbetapointer()
//  *merrorvarp = 1.0;

  fc_ximu = FULLCOND(o,datamatrix(1,1,0.0),title+"_truecov_expectation",1,1,path+"_truecov_expectation.res");
  fc_ximu.setflags(MCMC::norelchange | MCMC::nooutput);
//  double * ximup = fc_ximu.getbetapointer()
//  *ximup = 0.0;

  fc_xivar = FULLCOND(o,datamatrix(1,1,1.0),title+"_truecov_var",1,1,path+"truecov_var.res");
  fc_xivar.setflags(MCMC::norelchange | MCMC::nooutput);
//  double * xivarp = fc_xivar.getbetapointer()
//  *xivarp = 1;

  index = statmatrix<int>(beta.rows(),1,0);
  for(i=0; i<beta.rows(); i++)
    {
    index(i,0)=i;
    }
  }
// END: merror

  // COPY CONSTRUCTOR

  fullcond_merror::fullcond_merror(const fullcond_merror & m)
    : FULLCOND(FULLCOND(m))
    {
    designp = m.designp;
    likep = m.likep;
// BEGIN: merror
    splinep = m.splinep;
    minx = m.minx;
    maxx = m.maxx;
    meandata = m.meandata;
    old = m.old;
    index = m.index;
    merror = m.merror;
    currentspline = m.currentspline;
    diffspline = m.diffspline;
    logfcold = m.logfcold;
    logfcnew = m.logfcnew;
    fc_merrorvar = m.fc_merrorvar;
    fc_ximu = m.fc_ximu;
    fc_xivar = m.fc_xivar;
    generrcount = m.generrcount;
    pathresults = m.pathresults;
// END: merror
    }

  // OVERLOADED ASSIGNMENT OPERATOR

  const fullcond_merror & fullcond_merror::operator=(const fullcond_merror & m)
    {
    if (this == &m)
      return *this;
    FULLCOND::operator=(FULLCOND(m));

    designp = m.designp;
    likep = m.likep;
// BEGIN: merror
    splinep = m.splinep;
    minx = m.minx;
    maxx = m.maxx;
    meandata = m.meandata;
    old = m.old;
    index = m.index;
    merror = m.merror;
    fc_merrorvar = m.fc_merrorvar;
    fc_ximu = m.fc_ximu;
    fc_xivar = m.fc_xivar;
    generrcount = m.generrcount;
    pathresults = m.pathresults;
// END: merror

    return *this;
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void fullcond_merror::update(void)
    {
    if(varcoeff)
      {
      designp->init_data_varcoeff(beta);
      }
    else
      {
      // some miscellaneous variables
      unsigned i, j;
      double u;

      // store results from previous iteration in old
      old = beta;

      // sampling step:
      // standard deviation of measurement error
      double * merrorvarp = fc_merrorvar.getbetapointer();
      double anew, bnew;
      anew = 0.001 + 0.5*beta.rows()*data.cols();
      bnew = 0.0;
      for(i=0; i<beta.rows(); i++)
        for(j=0; j<data.cols(); j++)
          bnew += (data(i,j)-beta(i,0))*(data(i,j)-beta(i,0));
      bnew = 0.001 + 0.5*bnew;
      *merrorvarp = rand_invgamma(anew,bnew);
      double mesd = sqrt(*merrorvarp);
      fc_merrorvar.update();

      // sampling step:
      // standard deviation of true covariate values
      // sampling step:
      // standard deviation of true covariate values
      double * xivarp = fc_xivar.getbetapointer();
      double * ximup = fc_ximu.getbetapointer();

      anew = 0.001 + 0.5*beta.rows();
      bnew = 0.0;
      for(i=0; i<beta.rows(); i++)
        bnew += (beta(i,0)-*ximup)*(beta(i,0)-*ximup);
      bnew = 0.001 + 0.5*bnew;
      *xivarp = rand_invgamma(anew,bnew);
      double priorsd = sqrt(*xivarp);
      fc_xivar.update();

      double muhelp = beta.sum(0)*1000/(beta.rows()*1000*1000 + *xivarp);
      double sdhelp = *xivarp *1000*1000 / (beta.rows()*1000*1000 + *xivarp);
      *ximup = muhelp + sdhelp*rand_normal();
      double priormean = *ximup;
      fc_ximu.update();

      // generate proposed values (random walk proposal according to Berry et al.)
      double * work = beta.getV();
      double * workold = old.getV();
      for(i=0;i<beta.rows();i++,work++,workold++)
        {
        *work = *workold + 2*mesd*rand_normal()/(double)merror;
        generrtrial++;
        while(*work<minx || *work>maxx)
          {
          generrcount++;
//          optionsp->out("  WARNING in "+title+":");
//          optionsp->out("          Generated true covariate value out of range!");
          *work = *workold + 2*mesd*rand_normal()/(double)merror;
          generrtrial++;
          }
        }

      // extract current f(x) from spline_basis
      currentspline = splinep->get_spline();

      // call update_merror and compute new values for f(x).
      splinep->update_merror(beta);
      diffspline = splinep->get_spline()-currentspline;

      // full conditional for old values
      for(i=0; i<beta.rows(); i++)
        {
        logfcold(i,0) = likep->loglikelihood(i,i,index) - 0.5*((old(i,0)-priormean)/priorsd)*((old(i,0)-priormean)/priorsd);
        for(j=0; j<merror; j++)
          {
          logfcold(i,0) -= 0.5*((data(i,j)-old(i,0))/mesd)*((data(i,j)-old(i,0))/mesd);
          }
        }

      // full conditional for proposed values
      likep->addtocurrent(diffspline);
      for(i=0; i<beta.rows(); i++)
        {
        logfcnew(i,0) = likep->loglikelihood(i,i,index,false) - 0.5*((beta(i,0)-priormean)/priorsd)*((beta(i,0)-priormean)/priorsd);
        for(j=0; j<merror; j++)
          {
          logfcnew(i,0) -= 0.5*((data(i,j)-beta(i,0))/mesd)*((data(i,j)-beta(i,0))/mesd);
          }
        }

      // compute acceptance probabilities; overwrite non-accepted values with old values
      for(i=0; i<beta.rows(); i++)
        {
        u = log(uniform());
        nrtrials++;
        if(u <= (logfcnew(i,0)-logfcold(i,0)))
          {
          acceptance += 1.0;
          }
        else
          {
          beta(i,0) = old(i,0);
          }
        }

      // update spline_basis with the final information
      splinep->update_merror(beta);

      // Update the linear predictor
      diffspline = splinep->get_spline()-currentspline;
      likep->addtocurrent(diffspline);
      likep->swap_linearpred();

      FULLCOND::update();
      }
    }


  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool fullcond_merror::posteriormode(void)
    {
    }

  bool fullcond_merror::posteriormode_converged(const unsigned & itnr)
    {
    }

  void fullcond_merror::posteriormode_set_beta_mode(void)
    {
    }

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void fullcond_merror::outoptions(void)
    {
    }

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void fullcond_merror::outresults(void)
    {
    FULLCOND::outresults();

    ST::string l1 = ST::doubletostring(lower1,4);
    ST::string l2 = ST::doubletostring(lower2,4);
    ST::string u1 = ST::doubletostring(upper1,4);
    ST::string u2 = ST::doubletostring(upper2,4);

    ST::string nl1 = ST::doubletostring(lower1,4);
    ST::string nl2 = ST::doubletostring(lower2,4);
    ST::string nu1 = ST::doubletostring(upper1,4);
    ST::string nu2 = ST::doubletostring(upper2,4);
    nl1 = nl1.replaceallsigns('.','p');
    nl2 = nl2.replaceallsigns('.','p');
    nu1 = nu1.replaceallsigns('.','p');
    nu2 = nu2.replaceallsigns('.','p');

    ST::string vstr;

    ofstream ou(pathresults.strtochar());

    unsigned i;
    ou << "pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
    nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
    for(i=0; i<beta.rows(); i++)
      {
      ou << betamean(i,0) << "  ";
      ou << (betavar(i,0)<0.0?0.0:sqrt(betavar(i,0))) << "  ";
      ou << betaqu_l1_lower(i,0) << "  ";
      ou << betaqu_l2_lower(i,0) << "  ";
      ou << betaqu50(i,0) << "  ";
      ou << betaqu_l2_upper(i,0) << "  ";
      ou << betaqu_l1_upper(i,0) << "  ";
      ou << betamin(i,0) << "  ";
      ou << betamax(i,0) << "  " << endl;
      }

    optionsp->out("  Results for the covariate values are stored in file\n");
    optionsp->out("  " + pathresults + "\n");

    optionsp->out("\n");

    if(generrcount>0)
      {
      optionsp->out("  WARNING: "+ST::doubletostring((double)generrcount/(double)generrtrial,5)+"% of the generated covariate values were out of the specified range!");
      optionsp->out("           Consider making the grid more wide.\n\n");
      }

// Variance of the measurement error

    fc_merrorvar.outresults();
    ST::string pathhelp = pathresults.substr(0,pathresults.length()-7)+"merrorvar_sample.raw";
    fc_merrorvar.get_samples(pathhelp);

    optionsp->out("\n");
    optionsp->out("  "+fc_merrorvar.get_title()+"\n");
    optionsp->out("\n");
    optionsp->out("\n");

    vstr = "  Mean:         ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_merrorvar.get_betamean(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  Std. dev.:    ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring((fc_merrorvar.get_betavar(0,0)<0.0?0.0:sqrt(fc_merrorvar.get_betavar(0,0))),6) + "\n");

    vstr = "  " + l1 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_merrorvar.get_beta_lower1(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + l2 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_merrorvar.get_beta_lower2(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  50% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_merrorvar.get_betaqu50(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + u1 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_merrorvar.get_beta_upper2(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + u2 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_merrorvar.get_beta_upper1(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    optionsp->out("\n");

    ST::string merrorvar_pathresults = pathresults.substr(0,pathresults.length()-7) + "merror_var.res";

    ofstream ou1(merrorvar_pathresults.strtochar());

    ou1 << "pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
    nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
    ou1 << fc_merrorvar.get_betamean(0,0) << "  ";
    ou1 << (fc_merrorvar.get_betavar(0,0)<0.0?0.0:sqrt(fc_merrorvar.get_betavar(0,0))) << "  ";
    ou1 << fc_merrorvar.get_beta_lower1(0,0) << "  ";
    ou1 << fc_merrorvar.get_beta_lower2(0,0) << "  ";
    ou1 << fc_merrorvar.get_betaqu50(0,0) << "  ";
    ou1 << fc_merrorvar.get_beta_upper2(0,0) << "  ";
    ou1 << fc_merrorvar.get_beta_upper1(0,0) << "  ";
    ou1 << fc_merrorvar.get_betamin(0,0) << "  ";
    ou1 << fc_merrorvar.get_betamax(0,0) << "  " << endl;

    optionsp->out("  Results for the variance of the measurement error are also stored in file\n");
    optionsp->out("  " + merrorvar_pathresults + "\n");

    optionsp->out("\n");

// mean of the true covariate values    

    fc_ximu.outresults();
    pathhelp = pathresults.substr(0,pathresults.length()-7)+"ximu_sample.raw";
    fc_ximu.get_samples(pathhelp);

    optionsp->out("\n");
    optionsp->out("  "+fc_ximu.get_title()+"\n");
    optionsp->out("\n");
    optionsp->out("\n");

    vstr = "  Mean:         ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_ximu.get_betamean(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  Std. dev.:    ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring((fc_ximu.get_betavar(0,0)<0.0?0.0:sqrt(fc_ximu.get_betavar(0,0))),6) + "\n");

    vstr = "  " + l1 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_ximu.get_beta_lower1(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + l2 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_ximu.get_beta_lower2(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  50% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_ximu.get_betaqu50(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + u1 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_ximu.get_beta_upper2(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + u2 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_ximu.get_beta_upper1(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    optionsp->out("\n");

    ST::string ximu_pathresults = pathresults.substr(0,pathresults.length()-7) + "truecov_expectation.res";

    ofstream ou2(ximu_pathresults.strtochar());

    ou2 << "pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
    nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
    ou2 << fc_ximu.get_betamean(0,0) << "  ";
    ou2 << (fc_ximu.get_betavar(0,0)<0.0?0.0:sqrt(fc_ximu.get_betavar(0,0))) << "  ";
    ou2 << fc_ximu.get_beta_lower1(0,0) << "  ";
    ou2 << fc_ximu.get_beta_lower2(0,0) << "  ";
    ou2 << fc_ximu.get_betaqu50(0,0) << "  ";
    ou2 << fc_ximu.get_beta_upper2(0,0) << "  ";
    ou2 << fc_ximu.get_beta_upper1(0,0) << "  ";
    ou2 << fc_ximu.get_betamin(0,0) << "  ";
    ou2 << fc_ximu.get_betamax(0,0) << "  " << endl;

    optionsp->out("  Results for the expectation of the true covariate values are also stored in file\n");
    optionsp->out("  " + ximu_pathresults + "\n");

    optionsp->out("\n");

// variance of the true covariate values    

    fc_xivar.outresults();
    pathhelp = pathresults.substr(0,pathresults.length()-7)+"xivar_sample.raw";
    fc_xivar.get_samples(pathhelp);

    optionsp->out("\n");
    optionsp->out("  "+fc_xivar.get_title()+"\n");
    optionsp->out("\n");
    optionsp->out("\n");

    vstr = "  Mean:         ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_xivar.get_betamean(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  Std. dev.:    ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring((fc_xivar.get_betavar(0,0)<0.0?0.0:sqrt(fc_xivar.get_betavar(0,0))),6) + "\n");

    vstr = "  " + l1 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_xivar.get_beta_lower1(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + l2 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_xivar.get_beta_lower2(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  50% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_xivar.get_betaqu50(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + u1 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_xivar.get_beta_upper2(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + u2 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_xivar.get_beta_upper1(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    optionsp->out("\n");

    ST::string xivar_pathresults = pathresults.substr(0,pathresults.length()-7) + "truecov_var.res";

    ofstream ou3(xivar_pathresults.strtochar());

    ou3 << "pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
    nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
    ou3 << fc_xivar.get_betamean(0,0) << "  ";
    ou3 << (fc_xivar.get_betavar(0,0)<0.0?0.0:sqrt(fc_xivar.get_betavar(0,0))) << "  ";
    ou3 << fc_xivar.get_beta_lower1(0,0) << "  ";
    ou3 << fc_xivar.get_beta_lower2(0,0) << "  ";
    ou3 << fc_xivar.get_betaqu50(0,0) << "  ";
    ou3 << fc_xivar.get_beta_upper2(0,0) << "  ";
    ou3 << fc_xivar.get_beta_upper1(0,0) << "  ";
    ou3 << fc_xivar.get_betamin(0,0) << "  ";
    ou3 << fc_xivar.get_betamax(0,0) << "  " << endl;

    optionsp->out("  Results for the variance of the true covariate values are also stored in file\n");
    optionsp->out("  " + xivar_pathresults + "\n");

    optionsp->out("\n");

    }

  // FUNCTION: reset
  // TASK: resets all parameters

  void fullcond_merror::reset(void)
    {
    }

  vector<ST::string> & fullcond_merror::get_results_latex(void)
    {
    }

} // end: namespace MCMC

