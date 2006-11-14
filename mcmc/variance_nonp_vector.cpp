#include "first.h"
#include"variance_nonp_vector.h"



//------------------------------------------------------------------------------
//------------ CLASS: FULLCOND_ridge implementation of member functions --------
//------------------------------------------------------------------------------



namespace MCMC
{

using randnumbers::rand_inv_gaussian;

//______________________________________________________________________________
//
// CONSTRUCTOR with Parameters
//
// o         : pointer to MCMCoptions object
// p         : pointer to FULLCOND_const object
// d        : pointer to DISTRIBUTION object
// a
// b
// ti        : reference to title of the full conditional
//             (for example "fixed effects")
// fp        : reference to file path for storing sampled parameters
// fr        : reference to filename for storing results
// c         : reference to responsecategory
//             (important in the case of multivariate response)
//______________________________________________________________________________

FULLCOND_variance_nonp_vector::FULLCOND_variance_nonp_vector(MCMCoptions * o,
                         FULLCOND_const * p,DISTRIBUTION * d,
                         const vector<double> & a, const vector<double> & b,
                         const ST::string & ti, const ST::string & fp,
                         const ST::string & fr, const unsigned & c)
                         : FULLCOND(o,datamatrix(1,1),ti,1,1,fp)
    {
    fctype = MCMC::variance;

//  tau = Startwerte fuer Varianzen aus Optionen holen
    
    update_sigma2 = true;
    
    column = c;
    
    pathresults = fr;
    
    Cp = p;
    
    distrp = d;

    a_invgamma = a;
    b_invgamma = b;

//    priorassumptions.push_back("Inverse gamma prior for variance component with hyperparameters a="
//    +ST::doubletostring(a,6)+ " and b=" + ST::doubletostring(b,6) );

    priorassumptions.push_back("\\\\");

    ST::string path = samplepath.substr(0,samplepath.length()-4)+"_lasso.raw";

    fc_lasso = FULLCOND(o,datamatrix(1,1),Cp->get_title()+"_lasso",1,1,path);

    fc_lasso.setflags(MCMC::norelchange | MCMC::nooutput);

    lasso = 1;                    // Startwert, aus Option holen
    a_lassogamma = 3.33;          //  Hyperparameter for lasso aus Option holen
    b_lassogamma = 4.44;          //  Hyperparameter for lasso aus Option holen
  
    setbeta(Cp->get_variances());
    }


//______________________________________________________________________________
//
// COPY CONSTRUCTOR
//______________________________________________________________________________

FULLCOND_variance_nonp_vector::FULLCOND_variance_nonp_vector(const FULLCOND_variance_nonp_vector & t)
    : FULLCOND(FULLCOND(t))
  {
  update_sigma2 = t.update_sigma2;
  column = t.column;
  pathresults = t.pathresults;
  Cp = t.Cp;
  distrp = t.distrp;
  a_invgamma = t.a_invgamma;
  b_invgamma = t.b_invgamma;
  fc_lasso = t.fc_lasso;
  lasso = t.lasso;
  a_lassogamma = t.a_lassogamma;
  b_lassogamma = t.b_lassogamma;
  }


//______________________________________________________________________________
//
// OVERLOADED ASSIGNMENT OPERATOR
//______________________________________________________________________________

const FULLCOND_variance_nonp_vector & FULLCOND_variance_nonp_vector::operator=(
                                        const FULLCOND_variance_nonp_vector & t)
  {
  if (this == &t)
    return *this;
  FULLCOND::operator=(FULLCOND(t));
  update_sigma2 = t.update_sigma2;
  column = t.column;
  pathresults = t.pathresults;
  Cp = t.Cp;
  distrp = t.distrp;
  a_invgamma = t.a_invgamma;
  b_invgamma = t.b_invgamma;
  fc_lasso = t.fc_lasso;
  lasso = t.lasso;
  a_lassogamma = t.a_lassogamma;
  b_lassogamma = t.b_lassogamma;
  return *this;
  }


//______________________________________________________________________________
//
// FUNCTION: update
// TASK: - stores sampled parameters in file 'samplepath'
//         storing order: first row, second row, ...
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector::update(void)
  {
  acceptance++;
  
  double * lassop = fc_lasso.getbetapointer();        // current value of lassoparameter

  // Gibbs-Update of varianceparameters 1/tau^2 with Inverse Normaldistribution
  //--------------------------------------------------------------------
  double * workbeta = Cp->getbetapointer();           // current value of first regressionparameter

  ofstream output("c:/bayesx/test/test.txt", ios::out|ios::app);
  ofstream outputlasso("c:/bayesx/test/lasso.txt", ios::out|ios::app);
  outputlasso << *lassop << "\n";
  output << "fall " << "index " << "randinvgaussian " << "beta " << "lasso " <<"\n";

  for(unsigned int i=0; i<nrpar; i++, workbeta++)
    {
    if (*workbeta>0 && *lassop>0)
      {
      double rand_invgaussian = rand_inv_gaussian((*lassop)/(*workbeta), (*lassop)*(*lassop));

      output << "*workbeta>0 && *lassop>0 " << i << " " << rand_invgaussian << " " << *workbeta << " " << *lassop << "\n" ;

      beta(i,0) = 1.0/rand_invgaussian;
      }
    if (*workbeta<0 && *lassop>0)
      {
      double rand_invgaussian = rand_inv_gaussian(-1.0*(*lassop)/(*workbeta), (*lassop)*(*lassop));

      output << "*workbeta<0 && *lassop>0 " << i << " " << rand_invgaussian << " " << *workbeta << " " << *lassop << "\n" ;

      beta(i,0) = 1.0/(rand_invgaussian + 1E-8);
      }
    if (*workbeta==0 || *lassop<=0)
      {
      output <<"else" << i << " " << *workbeta << " " << *lassop << "\n" ;
      beta(i,0) = 1E-6;
      }
    }
  // Ergebnisse zu den Regressionskoeffizienten zurueckschreiben
  Cp->update_variances(beta);

  // Gibbs-Update of lassoparameter with Gammadistribution
  //--------------------------------------------------------

  double sumvariances = 0;
  double * workbetamean = betamean.getV();               // mean of variances

  datamatrix variances = Cp->get_variances();          // current variances
  for(unsigned int i=0; i<nrpar; i++)
    {
    sumvariances = sumvariances + variances(i,0);      // sum of current variances
    }
  *lassop = sqrt(rand_gamma(nrpar + a_lassogamma, b_lassogamma + 0.5*sumvariances));

  outputlasso << "beta1 " << beta(0,0) << "\n";
  outputlasso << "betamean " << *workbetamean << "\n";
  outputlasso << "betamean " << *(workbetamean + 1) << "\n";
  outputlasso << "sample " << optionsp->get_samplesize() << "\n";
  outputlasso << "burnin " << optionsp->get_burnin() << "\n";
  outputlasso << "iterationen " << optionsp->get_nriter() << "\n";
  outputlasso << "lassonew " << *lassop << "\n";

  // Update lassoparameter
  fc_lasso.update();

  // Update varianceparameter
  FULLCOND::update();
  }


//______________________________________________________________________________
//
// FUNCTION: outresults
// TASK: - write results to output window and files
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector::outresults(void)
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

  optionsp->out("  Results for the variance components are stored in file\n");
  optionsp->out("  " + pathresults + "\n");

  optionsp->out("\n");

  outresults_lasso();
  }


//______________________________________________________________________________
//
// FUNCTION: outresults_lasso
// TASK: - write results to output window and files
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector::outresults_lasso(void)
  {

  fc_lasso.outresults();

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

  optionsp->out("\n");
  optionsp->out("  f_ridge_lasso \n");
  optionsp->out("\n");
  optionsp->out("\n");

  vstr = "  Mean:         ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_lasso.get_betamean(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  vstr = "  Std. dev.:    ";
  optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring((fc_lasso.get_betavar(0,0)<0.0?0.0:sqrt(fc_lasso.get_betavar(0,0))),6) + "\n");

  vstr = "  " + l1 + "% Quantile: ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_lasso.get_beta_lower1(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  vstr = "  " + l2 + "% Quantile: ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_lasso.get_beta_lower2(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  vstr = "  50% Quantile: ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_lasso.get_betaqu50(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  vstr = "  " + u1 + "% Quantile: ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_lasso.get_beta_upper2(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  vstr = "  " + u2 + "% Quantile: ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_lasso.get_beta_upper1(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  optionsp->out("\n");

  ST::string lasso_pathresults = pathresults.substr(0,pathresults.length()-7) + "lasso.res";

  ofstream ou(lasso_pathresults.strtochar());

  ou << "pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
  nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
  ou << fc_lasso.get_betamean(0,0) << "  ";
  ou << (fc_lasso.get_betavar(0,0)<0.0?0.0:sqrt(fc_lasso.get_betavar(0,0))) << "  ";
  ou << fc_lasso.get_beta_lower1(0,0) << "  ";
  ou << fc_lasso.get_beta_lower2(0,0) << "  ";
  ou << fc_lasso.get_betaqu50(0,0) << "  ";
  ou << fc_lasso.get_beta_upper2(0,0) << "  ";
  ou << fc_lasso.get_beta_upper1(0,0) << "  ";
  ou << fc_lasso.get_betamin(0,0) << "  ";
  ou << fc_lasso.get_betamax(0,0) << "  " << endl;

  optionsp->out("  Results for the lasso parameter are also stored in file\n");
  optionsp->out("  " + lasso_pathresults + "\n");

  optionsp->out("\n");
  }


//______________________________________________________________________________
//
// FUNCTION: outoptions
// TASK: - write options to output window
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector::outoptions(void)
  {
  optionsp->out("  Hyperprior a for variance parameter: " +
                   ST::doubletostring(a_invgamma[0]) + "\n" );
  optionsp->out("  Hyperprior b for variance parameter: " +
                   ST::doubletostring(b_invgamma[0]) + "\n" );
  optionsp->out("  Hyperprior a for lasso parameter: " +
                   ST::doubletostring(a_lassogamma) + "\n" );
  optionsp->out("  Hyperprior b for lasso parameter: " +
                   ST::doubletostring(b_lassogamma) + "\n" );
  optionsp->out("\n");
  }


} // end: namespace MCMC






