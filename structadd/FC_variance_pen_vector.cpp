
#include"FC_variance_pen_vector.h"



//------------------------------------------------------------------------------
//-- CLASS: FC_variance_pen_vector implementation of member functions -----------
//------------------------------------------------------------------------------


namespace MCMC
{

using randnumbers::rand_inv_gaussian;


void FC_variance_pen_vector::add_variable(datamatrix & x,
                    double la, double shrink,
                    bool sfix, double a, double b)
  {

  }                    



void FC_variance_pen_vector::read_options(vector<ST::string> & op,vector<ST::string> & vn)
  {

  /*
  1       degree
  2       numberknots
  3       difforder
  4       lambda
  5       a
  6       b
  7       center
  8       map
  9       lambda_re
  10      a_re
  11      b_re
  12      internal_mult
  13      samplemult
  14      constraints
  15      round
  16      centermethod
  17      internal_mult
  18      pvalue
  19      meaneffect
  20      binning
  21      update
  */


//  is_ridge = isridge;


  }


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

FC_variance_pen_vector::FC_variance_pen_vector(MASTER_OBJ * mp,
                                          GENERAL_OPTIONS * o, FC_linear_pen * p,
                                          DISTR * d,const ST::string & ti,
                         const ST::string & fp, bool isr,
                          vector<ST::string> & op,
                         vector<ST::string> & vn)
                         : FC(o,ti,1,1,fp)
    {

    read_options(op,vn);

    is_ridge = isr;

    update_sigma2 = true;

    Cp = p;

    distrp = d;

    priorassumptions.push_back("\\\\");


    // vector<ST::string> varnames = fc_shrinkage.get_datanames();



    //Initialisieren der Betamatrizen für die Varianzan + Übergabe der Startwerte

    lassosum = 0;
    ridgesum = 0;

    }
//______________________________________________________________________________
//
// COPY CONSTRUCTOR
//______________________________________________________________________________

FC_variance_pen_vector::FC_variance_pen_vector(const FC_variance_pen_vector & t)
    : FC(FC(t))
  {
  tau = t.tau;
  lambda = t.lambda;
  update_sigma2 = t.update_sigma2;
  Cp = t.Cp;
  distrp = t.distrp;
  fc_shrinkage = t.fc_shrinkage;
  shrinkagefix = t.shrinkagefix;
  a_shrinkagegamma = t.a_shrinkagegamma;
  b_shrinkagegamma = t.b_shrinkagegamma;
  shrinkagestart = t.shrinkagestart;
  lassosum = t.lassosum;
  ridgesum = t.ridgesum;
  is_ridge = t.is_ridge;
  }


//______________________________________________________________________________
//
// OVERLOADED ASSIGNMENT OPERATOR
//______________________________________________________________________________

const FC_variance_pen_vector & FC_variance_pen_vector::operator=(
                                        const FC_variance_pen_vector & t)
  {
  if (this == &t)
    return *this;
  FC::operator=(FC(t));
  tau = t.tau;
  lambda = t.lambda;
  update_sigma2 = t.update_sigma2;
  Cp = t.Cp;
  distrp = t.distrp;
  fc_shrinkage = t.fc_shrinkage;
  shrinkagefix = t.shrinkagefix;
  a_shrinkagegamma = t.a_shrinkagegamma;
  b_shrinkagegamma = t.b_shrinkagegamma;
  shrinkagestart = t.shrinkagestart;
  lassosum = t.lassosum;
  ridgesum = t.ridgesum;
  is_ridge = t.is_ridge;
  return *this;
  }


bool FC_variance_pen_vector::posteriormode(void)
  {
  return true;
  }


//______________________________________________________________________________
//
// FUNCTION: update
// TASK: - stores sampled parameters in file 'samplepath'
//         storing order: first row, second row, ...
//______________________________________________________________________________

void FC_variance_pen_vector::update(void)
  {




  fc_shrinkage = FC(optionsp,"",shrinkagefix.size(),1,samplepath + ".shrinkage");

  /*
  double * shrinkagep fc_skrinkage.beta.getV();

  acceptance++;
  unsigned i, j, k;

  // get current value of shrinkagearameter
  double * shrinkagep = fc_shrinkage.getbetapointer();

// get current varianceparameters
//  datamatrix variances = datamatrix(beta.rows(),1,0);
//  for(i=0; i<cut.size()-1; i++)
//    variances.putRowBlock(cut[i],cut[i+1],Cp[i]->get_variances());
//  variances = beta;

  // getcurrent value of sqrt(scale) parameter
  double help = sqrt(distrp->get_scale(column));


  // variable for current value regressionparameters
  double * workbeta;

  // Vaeiables for summs
  double rand_invgaussian = 0;
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
          rand_invgaussian = rand_inv_gaussian(help*(*shrinkagep)/(*workbeta), (*shrinkagep)*(*shrinkagep));
          beta(i,0) = 1.0/rand_invgaussian;
          }
        if (*workbeta<0 && *shrinkagep>0)
          {
          rand_invgaussian = rand_inv_gaussian(-1.0*help*(*shrinkagep)/(*workbeta), (*shrinkagep)*(*shrinkagep));
          beta(i,0) = 1.0/(rand_invgaussian);
          }
        if (*workbeta==0 || *shrinkagep<=0)
          {
          beta(i,0) = 1E-6;
          }
        lassosum = lassosum + ((*workbeta)*(*workbeta))/beta(i,0);  // sum(beta^2/tau^2)
       }

      if (is_ridge == 1)                              // L2-penalty
        {
         beta(i,0) = 1/(2*(*shrinkagep));
         ridgesum = ridgesum + ((*workbeta)*(*workbeta))/beta(i,0);  // sum(beta^2/tau^2)

        }

      sumregcoeff = sumregcoeff + (*workbeta)*(*workbeta);
      }
    }


// Update varianceparameter tau^2
//  datamatrix temp;
//  for(i=0; i<cut.size()-1; i++)
//    {
//    temp = beta.getRowBlock(cut[i],cut[i+1]);
//    Cp[i]->update_variances(temp);
//    }

  double * varp;
  k=0;
  for(i=0; i<cut.size()-1; i++)
    {
    varp = Cp[i]->getvariancespointer();
    for(j=cut[i]; j<cut[i+1]; j++, varp++, k++)
      {
      *varp = beta(k,0);
      }
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
      sumvariances = sumvariances + beta(i,0);      // sum(tau^2) of current variances
      }
    *shrinkagep = sqrt(rand_gamma(nrpar + a_shrinkagegamma, b_shrinkagegamma + 0.5*sumvariances));
    }

  if(shrinkagefix==false && is_ridge == 1)            // L2-penalty
    {
    *shrinkagep = rand_gamma(0.5*nrpar + a_shrinkagegamma, b_shrinkagegamma + sumregcoeff/(help*help));
    }

  // Update Shrinkageparameter
  fc_shrinkage.update();

  FC::update();
  */
  }



//______________________________________________________________________________
//
// FUNCTION: outresults
// TASK: - write results for varianceparameters to output window and files
//______________________________________________________________________________

void FC_variance_pen_vector::outresults(ofstream & out_stata, ofstream & out_R,
                  const ST::string & pathresults)
  {
  FC::outresults(out_stata,out_R,pathresults);

  /*
  unsigned int i,j,k  = 0;
  vector<ST::string> varnames(nrpar);
  vector<ST::string> helpvarnames;

  for(j=0; j<cut.size()-1; j++)
    {
    helpvarnames = Cp[j]->get_datanames();
    i = 0;
    for(k=cut[j]; k<cut[j+1]; k++, i++)
      {
      varnames[k] = helpvarnames[i];
      }
    }

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

//  unsigned i;
  ou << "varname  pmean  pstd  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
  nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
  for(i=0; i<beta.rows(); i++)
    {
    ou << varnames[i] << "   ";
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
  */
  }

void FC_variance_pen_vector::get_samples(const ST::string & filename,ofstream & outg) const
  {

  FC::get_samples(filename,outg);


  ST::string pathhelp = filename.substr(0,filename.length()-7)+"shrinkage_sample.raw";

  fc_shrinkage.get_samples(pathhelp,outg);

  }

//______________________________________________________________________________
//
// FUNCTION: outresults_shrinkage
// TASK: - write results for shrinkaageparameter to output window and files
//______________________________________________________________________________

void FC_variance_pen_vector::outresults_shrinkage(void)
  {

  // fc_shrinkage.outresults();

  /*
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

  ou << "varname  pmean  pstd  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
  nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
  ou << "lambda" << "  ";
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
  */
  }


//______________________________________________________________________________
//
// FUNCTION: outoptions
// TASK: - write options to output window
//______________________________________________________________________________

void FC_variance_pen_vector::outoptions(void)
  {
  /*
  if(is_ridge == 0)
    {
    optionsp->out("  Hyperparameter a for lassoshrinkage-gammapriori: " +
                     ST::doubletostring(a_shrinkagegamma) + "\n" );
    optionsp->out("  Hyperparameter b for lassoshrinkage-gammapriori: " +
                     ST::doubletostring(b_shrinkagegamma) + "\n" );
    if(shrinkagefix==true)
    {
    optionsp->out("  lassoshrinkageparameter is fixed at value shrinkage_start = " +
                     ST::doubletostring(fc_shrinkage.getbeta(0,0)) + "\n" );
    }
    optionsp->out("\n");
    }

  if(is_ridge == 1)
    {
    optionsp->out("  Hyperparameter a for ridgeshrinkage-gammapriori: " +
                     ST::doubletostring(a_shrinkagegamma) + "\n" );
    optionsp->out("  Hyperparameter b for ridgeshrinkage-gammapriori: " +
                     ST::doubletostring(b_shrinkagegamma) + "\n" );
    if(shrinkagefix==true)
    {
    optionsp->out("  ridgeshrinkageparameter is fixed at value shrinkage_start = " +
                     ST::doubletostring(fc_shrinkage.getbeta(0,0)) + "\n" );
    }
    optionsp->out("\n");
    }
  */
  }


} // end: namespace MCMC





