#include "first.h"
#include"variance_nonp_vector.h"



//------------------------------------------------------------------------------
//-- CLASS: FULLCOND_variance_nonp_vector implementation of member functions ---
//------------------------------------------------------------------------------


namespace MCMC
{

using randnumbers::rand_inv_gaussian;

//______________________________________________________________________________
//
// CONSTRUCTOR with Parameters
//
// o                 : pointer to MCMCoptions object
// p                 : pointer to FULLCOND_const object
// d                 : pointer to DISTRIBUTION object
// ti                : reference to title of the full conditional(for example "fixed effects")
// fp                : reference to file path for storing sampled parameters
// fr                : reference to filename for storing results
// shrinkage_start   : Starting value for the shrinkageparameter
// a_shrinkage_gamma : Hyperparameter for gammaprior of shrinkageparameter
// b_shrinkagegamma  : Hyperparameter for gammaprior of shrinkageparameter
// shrinkage_fix     : Should the shrinkageparameter be fixed at the value "shrinkage_start"
// is.ridge          : variable that indicates if L2- or L1- penalty for regressioncoefficients is uesd
// ct                : blocks of regression coefficients
// c                 : reference to responsecategory (important in the case of multivariate response)
//______________________________________________________________________________

FULLCOND_variance_nonp_vector::FULLCOND_variance_nonp_vector(MCMCoptions * o,
                         vector<FULLCOND_const*> & p,DISTRIBUTION * d,
                         const ST::string & ti, const ST::string & fp,
                         const ST::string & fr, const double & shrinkage_start,
                         const double & a_shrinkage_gamma, const double & b_shrinkage_gamma,
                         const bool & shrinkage_fix, const bool & isridge,
                         const vector<unsigned> & ct, const unsigned & c)
                         : FULLCOND(o,datamatrix(1,1),ti,1,1,fp)
    {

    fctype = MCMC::variance;

    update_sigma2 = true;

    column = c;

    pathresults = fr;

    Cp = p;

    distrp = d;

    cut = ct;
    is_ridge = isridge;

    priorassumptions.push_back("\\\\");

    ST::string path = pathresults.substr(0,pathresults.length()-4)+"_shrinkage.raw";

    fc_shrinkage = FULLCOND(o,datamatrix(1,1),Cp[0]->get_title()+"_shrinkage",1,1,path);
    fc_shrinkage.setflags(MCMC::norelchange | MCMC::nooutput);

    vector<ST::string> varnames = fc_shrinkage.get_datanames();

    // current value of shrinkageparameter lambda
    double * shrinkagep = fc_shrinkage.getbetapointer();

    // Startwerte setzen aus Option
    * shrinkagep = shrinkage_start;
    shrinkagefix = shrinkage_fix;
    a_shrinkagegamma = a_shrinkage_gamma;
    b_shrinkagegamma = b_shrinkage_gamma;


    //Initialisieren der Betamatrizen für die Varianzan + Übergabe der Startwerte
    datamatrix help;
    if (is_ridge == 0)                                       // L1-penalty
      {
      lassosum = 0;
      help = datamatrix(cut[cut.size()-1],1,0);
      }
    if (is_ridge == 1)                                       // L2-penalty
      {
      ridgesum = 0;
      help = datamatrix(cut[cut.size()-1],1,0);
      }

    unsigned i;
    for(i=0; i<cut.size()-1; i++)
      help.putRowBlock(cut[i],cut[i+1],Cp[i]->get_variances());
    setbeta(help);
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
  fc_shrinkage = t.fc_shrinkage;
  shrinkagefix = t.shrinkagefix;
  a_shrinkagegamma = t.a_shrinkagegamma;
  b_shrinkagegamma = t.b_shrinkagegamma;
  lassosum = t.lassosum;
  ridgesum = t.ridgesum;
  cut = t.cut;
  is_ridge = t.is_ridge;
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
  fc_shrinkage = t.fc_shrinkage;
  shrinkagefix = t.shrinkagefix;
  a_shrinkagegamma = t.a_shrinkagegamma;
  b_shrinkagegamma = t.b_shrinkagegamma;
  lassosum = t.lassosum;
  ridgesum = t.ridgesum;
  cut = t.cut;
  is_ridge = t.is_ridge;
  return *this;
  }


  // Pointer auf das Shrinkagearameter lambda Fullcond-Objekt
FULLCOND * FULLCOND_variance_nonp_vector::get_shrinkagepointer()
  {
  return &fc_shrinkage;
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
  unsigned i, j, k;

  // get current value of shrinkagearameter
  double * shrinkagep = fc_shrinkage.getbetapointer();

  // get current varianceparameters
  datamatrix variances = datamatrix(beta.rows(),1,0);
  for(i=0; i<cut.size()-1; i++)
    variances.putRowBlock(cut[i],cut[i+1],Cp[i]->get_variances());

  // getcurrent value of sqrt(scale) parameter
  double help = sqrt(distrp->get_scale(column));

  // variable for current value regressionparameters
  double * workbeta;

  // Vaeiables for summs
  double sumvariances = 0;
  double sumregcoeff = 0;
  lassosum=0;
  ridgesum=0;


  // Gibbs-Update of varianceparameters 1/tau^2 with Inverse Normaldistribution
  // if L1-penalty is used
  //---------------------------------------------------------------------------
  i=0;
  for(j=0; j<cut.size()-1; j++)
    {
    workbeta = Cp[j]->getbetapointer();               // current value of first regressionparameter
    for(k=cut[j]; k<cut[j+1]; k++, i++, workbeta++)
      {
      if (is_ridge == 0)                              // L1-penalty
        {
        if (*workbeta>0 && *shrinkagep>0)
          {
          double rand_invgaussian = rand_inv_gaussian(help*(*shrinkagep)/(*workbeta), (*shrinkagep)*(*shrinkagep));
          beta(i,0) = 1.0/rand_invgaussian;
          }
        if (*workbeta<0 && *shrinkagep>0)
          {
          double rand_invgaussian = rand_inv_gaussian(-1.0*help*(*shrinkagep)/(*workbeta), (*shrinkagep)*(*shrinkagep));
          beta(i,0) = 1.0/(rand_invgaussian);
          }
        if (*workbeta==0 || *shrinkagep<=0)
          {
          beta(i,0) = 1E-6;
          }
       }

      if (is_ridge == 1)                              // L2-penalty
        {
         beta(i,0) = 1/(2*(*shrinkagep));
        }

      sumregcoeff = sumregcoeff + (*workbeta)*(*workbeta);
      ridgesum = ridgesum + ((*workbeta)*(*workbeta))/variances(i,0);  // sum(beta^2/tau^2)
      lassosum = ridgesum;
      }
    }


  // Update varianceparameter tau^2
  datamatrix temp;
  for(i=0; i<cut.size()-1; i++)
    {
    temp = beta.getRowBlock(cut[i],cut[i+1]);
    Cp[i]->update_variances(temp);
    }

  // Ergebnisse zu ridgesum in Distributionobjekt zurückschreiben
  if (is_ridge == 0)                              // L1-penalty
    {
    distrp->update_lasso(lassosum);
    }
  if (is_ridge == 1)                              // L2-penalty
    {
    distrp->update_ridge(ridgesum);
    }



  // Gibbs-Update of Shrinkageparameter with Gammadistribution
  //----------------------------------------------------------
  if(shrinkagefix==false && is_ridge == 0)            // L1-penalty
    {
    for(i=0; i<nrpar; i++)
      {
      sumvariances = sumvariances + variances(i,0);      // sum(tau^2) of current variances
      }
    *shrinkagep = sqrt(rand_gamma(nrpar + a_shrinkagegamma, b_shrinkagegamma + 0.5*sumvariances));
    }

  if(shrinkagefix==false && is_ridge == 1)            // L2-penalty
    {
    *shrinkagep = rand_gamma(0.5*nrpar + a_shrinkagegamma, b_shrinkagegamma + sumregcoeff/(help*help));
    }

  // Update Shrinkageparameter
  fc_shrinkage.update();

  FULLCOND::update();
  }



//______________________________________________________________________________
//
// FUNCTION: outresults
// TASK: - write results for varianceparameters to output window and files
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

  outresults_shrinkage();
  bool shrinkagesamples=true;
  if(shrinkagesamples)
    {
    ST::string pathhelp = pathresults.substr(0,pathresults.length()-7)+"shrinkage_sample.raw";
    fc_shrinkage.get_samples(pathhelp);
    }

  }


//______________________________________________________________________________
//
// FUNCTION: outresults_shrinkage
// TASK: - write results for shrinkaageparameter to output window and files
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector::outresults_shrinkage(void)
  {

  fc_shrinkage.outresults();

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
  optionsp->out("  f_shrinkage_shrinkage \n");
  optionsp->out("\n");
  optionsp->out("\n");

  vstr = "  Mean:         ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_shrinkage.get_betamean(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  vstr = "  Std. dev.:    ";
  optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring((fc_shrinkage.get_betavar(0,0)<0.0?0.0:sqrt(fc_shrinkage.get_betavar(0,0))),6) + "\n");

  vstr = "  " + l1 + "% Quantile: ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_shrinkage.get_beta_lower1(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  vstr = "  " + l2 + "% Quantile: ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_shrinkage.get_beta_lower2(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  vstr = "  50% Quantile: ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_shrinkage.get_betaqu50(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  vstr = "  " + u1 + "% Quantile: ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_shrinkage.get_beta_upper2(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  vstr = "  " + u2 + "% Quantile: ";
  vstr = vstr + ST::string(' ',20-vstr.length())
              + ST::doubletostring(fc_shrinkage.get_beta_upper1(0,0),6);
  optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

  optionsp->out("\n");

  ST::string shrinkage_pathresults = pathresults.substr(0,pathresults.length()-7) + "shrinkage.res";

  ofstream ou(shrinkage_pathresults.strtochar());

  ou << "pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
  nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
  ou << fc_shrinkage.get_betamean(0,0) << "  ";
  ou << (fc_shrinkage.get_betavar(0,0)<0.0?0.0:sqrt(fc_shrinkage.get_betavar(0,0))) << "  ";
  ou << fc_shrinkage.get_beta_lower1(0,0) << "  ";
  ou << fc_shrinkage.get_beta_lower2(0,0) << "  ";
  ou << fc_shrinkage.get_betaqu50(0,0) << "  ";
  ou << fc_shrinkage.get_beta_upper2(0,0) << "  ";
  ou << fc_shrinkage.get_beta_upper1(0,0) << "  ";
  ou << fc_shrinkage.get_betamin(0,0) << "  ";
  ou << fc_shrinkage.get_betamax(0,0) << "  " << endl;

  optionsp->out("  Results for the shrinkage parameter are also stored in file\n");
  optionsp->out("  " + shrinkage_pathresults + "\n");

  optionsp->out("\n");
  }


//______________________________________________________________________________
//
// FUNCTION: outoptions
// TASK: - write options to output window
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector::outoptions(void)
  {
  optionsp->out("  Hyperparameter a for gamma-shrinkagepriori: " +
                   ST::doubletostring(a_shrinkagegamma) + "\n" );
  optionsp->out("  Hyperparameter b for gamma-shrinkagepriori: " +
                   ST::doubletostring(b_shrinkagegamma) + "\n" );
  if(shrinkagefix==true)
  {
  optionsp->out("  shrinkageparameter is fixed at value shrinkage = " +
                   ST::doubletostring(fc_shrinkage.getbeta(0,0)) + "\n" );
  }
  optionsp->out("\n");
  }


} // end: namespace MCMC





