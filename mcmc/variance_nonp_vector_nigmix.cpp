
#include"variance_nonp_vector_nigmix.h"



//------------------------------------------------------------------------------
//-- CLASS: FULLCOND_variance_nonp_vector_nigmix implementation of member functions ---
//------------------------------------------------------------------------------

namespace MCMC
{

using randnumbers::rand_inv_gaussian;
using randnumbers::bernoulli;
using randnumbers::rand_beta;

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
// ins               : Starting value for varianceparametercomponent indicator
// vv0               : Point for Indicatorfunction
// vv1               : Point for Indicatorfunction
// t2s               : Starting value for varianceparametercomponent t2
// at2               : Hyperparameter for inverse gammaprior of vainaceparametercomponent t2
// bt2               : Hyperparameter for inverse gammaprior of vainaceparametercomponent t2
// omegastart        : starting value for omega
// omf               : Should the omega be fixed at the value "omegastart"    
// ct                : blocks of regression coefficients
// c                 : reference to responsecategory (important in the case of multivariate response)
//______________________________________________________________________________

FULLCOND_variance_nonp_vector_nigmix::FULLCOND_variance_nonp_vector_nigmix(MCMCoptions * o,
                         vector<FULLCOND_const*> & p,DISTRIBUTION * d,
                         const ST::string & ti, const ST::string & fp,const ST::string & fr, 
                         const vector<double> & ins, const double & vv0, const double & vv1,
                         const vector<double> & t2s, const double & at2, const double & bt2,
                         const double & omegastart, const bool & omf,
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

    priorassumptions.push_back("\\\\");

    ST::string path = pathresults.substr(0,pathresults.length()-4)+"_shrinkage.raw";

    // Fullcondobject for shrinkageparameter omege
    fc_shrinkage = FULLCOND(o,datamatrix(1,1),Cp[0]->get_title()+"_shrinkage",1,1,path);
    fc_shrinkage.setflags(MCMC::norelchange | MCMC::nooutput);

    vector<ST::string> varnames = fc_shrinkage.get_datanames();

    // current value of shrinkageparameter omega
    double * shrinkagep = fc_shrinkage.getbetapointer();

    // Startwerte setzen aus Option
    * shrinkagep = omegastart;
    indicatorstart = ins;
    v0 = vv0;
    v1 = vv1;
    t2start = t2s;
    double test1 = t2start[0];
    double test2 = indicatorstart[0];
    a_t2 = at2;
    b_t2 = bt2;

    omegafix = omf;

    nigmixsum = 0;

    // Fullcondobject for Variancekomponents tau=indicator*t2
    fc_indicator = FULLCOND(o,datamatrix(nrpar,1),Cp[0]->get_title()+"_indicator",1,1,path);
    fc_indicator.setflags(MCMC::norelchange | MCMC::nooutput);
    fc_t2 = FULLCOND(o,datamatrix(nrpar,1),Cp[0]->get_title()+"_t2",1,1,path);
    fc_t2.setflags(MCMC::norelchange | MCMC::nooutput);

    
    // Initialisieren der Matrizen für die Varianzen
    datamatrix help;
    help = datamatrix(cut[cut.size()-1],1,0);

    double * workind = fc_indicator.getbetapointer();
    double * workt2 = fc_t2.getbetapointer();

    unsigned i;
    for(i=0; i<nrpar; i++, workind++, workt2++)
      {
      *workind = indicatorstart[i];
      *workt2 = t2start[i];
      }
    fc_indicator.update();
    fc_t2.update();

    for(i=0; i<cut.size()-1; i++)
      {
      help.putRowBlock(cut[i],cut[i+1],Cp[i]->get_variances());
      }
    setbeta(help);

    }
//______________________________________________________________________________
//
// COPY CONSTRUCTOR
//______________________________________________________________________________

FULLCOND_variance_nonp_vector_nigmix::FULLCOND_variance_nonp_vector_nigmix(const FULLCOND_variance_nonp_vector_nigmix & t)
    : FULLCOND(FULLCOND(t))
  {
  update_sigma2 = t.update_sigma2;
  column = t.column;
  pathresults = t.pathresults;
  Cp = t.Cp;
  distrp = t.distrp;
  fc_shrinkage = t.fc_shrinkage;
  indicatorstart = t.indicatorstart;
  v0 = t.v0;
  v1 = t.v1;
  t2start = t.t2start;
  a_t2 = t.a_t2;
  b_t2 = t.b_t2;
  omegafix = t.omegafix;    
  fc_indicator = t.fc_indicator;
  fc_t2 = t.fc_t2;
  nigmixsum = t.nigmixsum;
  cut = t.cut;

  }


//______________________________________________________________________________
//
// OVERLOADED ASSIGNMENT OPERATOR
//______________________________________________________________________________

const FULLCOND_variance_nonp_vector_nigmix & FULLCOND_variance_nonp_vector_nigmix::operator=(
                                        const FULLCOND_variance_nonp_vector_nigmix & t)
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
  indicatorstart = t.indicatorstart;
  v0 = t.v0;
  v1 = t.v1;
  t2start = t.t2start;
  a_t2 = t.a_t2;
  b_t2 = t.b_t2;
  omegafix = t.omegafix;    
  fc_indicator = t.fc_indicator;
  fc_t2 = t.fc_t2;
  nigmixsum = t.nigmixsum;
  cut = t.cut;
  return *this;
  }


  // Pointer auf das Shrinkagearameter lambda Fullcond-Objekt
FULLCOND * FULLCOND_variance_nonp_vector_nigmix::get_shrinkagepointer()
  {
  return &fc_shrinkage;
  }

//______________________________________________________________________________
//
// FUNCTION: update
// TASK: - stores sampled parameters in file 'samplepath'
//         storing order: first row, second row, ...
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector_nigmix::update(void)
  {
  acceptance++;
  unsigned i, j, k;

  // get current value of shrinkagearameter omega
  double * shrinkagep = fc_shrinkage.getbetapointer();

  // getcurrent value of sigma=sqrt(scale) parameter
  double help = sqrt(distrp->get_scale(column));

  // variable for current value of regressioncoefficients
  double * workbeta;

  // Variables for summs
  double probv1 = 0;
  unsigned int rand_bernoulli = 0;
  double sumvariances = 0;
  double sumregcoeff = 0;
  nigmixsum=0;

  double * workind = fc_indicator.getbetapointer();
  double * workt2 = fc_t2.getbetapointer();

  // Gibbs-Update of varianceparameters t2 with Inverse Gammadistribution and
  // indicator with Binomialdistribution
  //---------------------------------------------------------------------------
  i=0;
  for(j=0; j<cut.size()-1; j++)
    {
    workbeta = Cp[j]->getbetapointer();               // current value of first regressionparameter
    for(k=cut[j]; k<cut[j+1]; k++, i++, workbeta++, workind++, workt2++)
      {
      probv1 = 1/(1+((1-(*shrinkagep))/(*shrinkagep)*sqrt(v1/v0)*exp(-(1/v0-1/v1)*(*workbeta)*(*workbeta)/(2*help*help*(*workt2)))));
      rand_bernoulli = bernoulli(probv1);
      if(rand_bernoulli==0)
        {
        *workind = v0;
        }
      if(rand_bernoulli==1)
        {
        *workind = v1;
        }

      *workt2 = rand_invgamma(0.5*a_t2,b_t2+(*workbeta)*(*workbeta)/(2*help*help* *workind));

      beta(i,0) = *workind * *workt2;

      nigmixsum = nigmixsum + ((*workbeta)*(*workbeta))/beta(i,0);  // sum(beta^2/tau^2)
      }
    }
  fc_indicator.update();
  fc_t2.update();


  // Update varianceparameter tau^2
  //----------------------------------
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

  // Transfer of nigmixsum to the distribution objekt
  //-------------------------------------------------
    distrp->update_nigmix(nigmixsum);



  // Gibbs-Update of Shrinkageparameter omega with Betadistribution
  //----------------------------------------------------------------
  unsigned int sumindicatorv0 = 0;
  unsigned int sumindicatorv1 = 0;

  workind = fc_indicator.getbetapointer();
  if(omegafix==false)
    { for(unsigned int i=0; i<nrpar; i++, workind++)
        { if(*workind==v1)
          { sumindicatorv1 = sumindicatorv1 + 1;
          }
        }
    sumindicatorv0 = distrp->get_nrobs()-sumindicatorv1;
    *shrinkagep = rand_beta(1+sumindicatorv1,1+sumindicatorv0);
    }

  // Transfer of the updated Shrinkageparameter omega
  //--------------------------------------------------
  fc_shrinkage.update();

  FULLCOND::update();
  }



//______________________________________________________________________________
//
// FUNCTION: outresults
// TASK: - write results for varianceparameters to output window and files
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector_nigmix::outresults(void)
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

void FULLCOND_variance_nonp_vector_nigmix::outresults_shrinkage(void)
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

void FULLCOND_variance_nonp_vector_nigmix::outoptions(void)
  {
//  if(is_ridge == 0)
//    {
//    optionsp->out("  Hyperparameter a for shrinkage-gammapriori: " +
//                     ST::doubletostring(a_shrinkagegamma) + "\n" );
//    optionsp->out("  Hyperparameter b for lassoshrinkage-gammapriori: " +
//                     ST::doubletostring(b_shrinkagegamma) + "\n" );
//    if(shrinkagefix==true)
//    {
//    optionsp->out("  lassoshrinkageparameter is fixed at value shrinkage_start = " +
//                     ST::doubletostring(fc_shrinkage.getbeta(0,0)) + "\n" );
//    }
//    optionsp->out("\n");
//    }
//    
//  if(is_ridge == 1)
//    {
//    optionsp->out("  Hyperparameter a for ridgeshrinkage-gammapriori: " +
//                     ST::doubletostring(a_shrinkagegamma) + "\n" );
//    optionsp->out("  Hyperparameter b for ridgeshrinkage-gammapriori: " +
//                     ST::doubletostring(b_shrinkagegamma) + "\n" );
//    if(shrinkagefix==true)
//    {
//    optionsp->out("  ridgeshrinkageparameter is fixed at value shrinkage_start = " +
//                     ST::doubletostring(fc_shrinkage.getbeta(0,0)) + "\n" );
//    }
//    optionsp->out("\n");
//    }
  }


} // end: namespace MCMC





