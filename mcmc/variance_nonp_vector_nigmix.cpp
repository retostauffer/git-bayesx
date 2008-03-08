

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

    ST::string path1 = pathresults.substr(0,pathresults.length()-4)+"_shrinkage.raw";
    ST::string path2 = pathresults.substr(0,pathresults.length()-4)+"_indicator.raw";
    ST::string path3 = pathresults.substr(0,pathresults.length()-4)+"_t2.raw";

    // Initialize the matrix for the variances
    datamatrix help;
    help = datamatrix(cut[cut.size()-1],1,0);

    for(unsigned i=0; i<cut.size()-1; i++)
      {
      help.putRowBlock(cut[i],cut[i+1],Cp[i]->get_variances());
      }
    setbeta(help);


    // Fullcondobject for shrinkageparameter omega
    fc_shrinkage = FULLCOND(o,datamatrix(1,1),Cp[0]->get_title()+"_shrinkage",1,1,path1);
    fc_shrinkage.setflags(MCMC::norelchange | MCMC::nooutput);


    // Current value of shrinkageparameter omega
    double * shrinkagep = fc_shrinkage.getbetapointer();

    // set the starting values 
    * shrinkagep = omegastart;
    indicatorstart = ins;
    v0 = vv0;
    v1 = vv1;
    t2start = t2s;
    a_t2 = at2;
    b_t2 = bt2;

    omegafix = omf;

    nigmixsum = 0;
    
    // Fullcondobject for Variancekomponents indicator and t2
    fc_indicator = FULLCOND(o,datamatrix(nrpar,1),Cp[0]->get_title()+"_indicator",nrpar,1,path2);
    fc_indicator.setflags(MCMC::norelchange | MCMC::nooutput);
    fc_t2 = FULLCOND(o,datamatrix(nrpar,1),Cp[0]->get_title()+"_t2",nrpar,1,path3);
    fc_t2.setflags(MCMC::norelchange | MCMC::nooutput);

    
    // set the starting values of the Variancekomponents
    double * workind = fc_indicator.getbetapointer();
    double * workt2 = fc_t2.getbetapointer();
    for(unsigned int i=0; i<nrpar; i++, workind++, workt2++)
      {
      *workind = indicatorstart[i];
      *workt2 = t2start[i];
      }


  // set the variablenames of the fullcond objects
  unsigned int i, k  = 0;
  vector<ST::string> varnames_indicator(nrpar);
  vector<ST::string> varnames_t2(nrpar);
  vector<ST::string> helpvarnames;
  vector<unsigned int> helpfi(nrpar);

  for(unsigned j=0; j<cut.size()-1; j++)
    {
    helpvarnames = Cp[j]->get_datanames(); 
    i = 0;   
    for(k=cut[j]; k<cut[j+1]; k++, i++)
      { 
      varnames_t2[k] = "t2." + helpvarnames[i]; 
      varnames_indicator[k] = "I." + helpvarnames[i];
      helpfi[k] = 0; 
      }
    }
    
    fc_t2.init_names(varnames_t2);
    fc_indicator.init_names(varnames_indicator); 
    fc_shrinkage.init_name("w"); 

      

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

  nigmixsum=0;

  double * workind = fc_indicator.getbetapointer();
  double * workt2 = fc_t2.getbetapointer();

  // Gibbs-Update of varianceparameters t2 with Inverse Gammadistribution and
  // indicator with Binomialdistribution
  //---------------------------------------------------------------------------


//TEMP:BEGIN--------------------------------------------------------------------
ofstream output("c:/bayesx/test/test_nigmix.txt", ios::out|ios::app);
//TEMP:END----------------------------------------------------------------------


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

//TEMP:BEGIN--------------------------------------------------------------------
//output << "indikator " << *workind << " t2 " << *workt2 
//       << " var " << beta(i,0) << "\n" ;
//TEMP:END----------------------------------------------------------------------
      if(rand_bernoulli==0)
        {
        *workind = 0.0;
        }
      if(rand_bernoulli==1)
        {
        *workind = 1.0;
        }

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
    sumindicatorv0 = nrpar-sumindicatorv1;
    *shrinkagep = rand_beta(1+sumindicatorv1,1+sumindicatorv0);

//TEMP:BEGIN--------------------------------------------------------------------
//ofstream outputv1("c:/bayesx/test/w_nigmix.txt", ios::out|ios::app);
//for(unsigned int i=0; i<nrpar; i++)
//{outputv1 << "sumv1 " << sumindicatorv1 << " sumv0 " << sumindicatorv0 
//          << " w " << *shrinkagep << "\n" ;  
//}
//TEMP:END----------------------------------------------------------------------


    }

  // Transfer of the updated Shrinkageparameter omega
  //--------------------------------------------------
  fc_shrinkage.update();

  FULLCOND::update();
  }



//______________________________________________________________________________
//
// FUNCTION: outresults
// TASK: - write results for varianceparameters=idicator*t2 
//         to output window and files
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

  // Bildschirm-Ausgabe
  ST::string vstr;
  ST::string vmstr;
  ST::string vsdstr;
  ST::string v025str;
  ST::string v500str;
  ST::string v975str;

  optionsp->out("\n");
  optionsp->out("  f_shrinkage_variances \n");
  optionsp->out("\n");
  optionsp->out("\n");
  
  // Bildschirm-Ausgabe
  unsigned i;
//  optionsp->out("  Variable     mean           Std. Dev.      " + l1 +"% quant.    median         "+ u2 +"% quant.");
//  for(i=0; i<beta.rows(); i++)
//    {
//    vstr = "  var." + ST::inttostring(i+1) + " ";
//    vmstr = ST::doubletostring(betamean(i,0),6);
//    vsdstr = ST::doubletostring((betavar(i,0)<0.0?0.0:sqrt(betavar(i,0))),6);
//    v025str = ST::doubletostring(betaqu_l1_lower(i,0),6);
//    v500str = ST::doubletostring(betaqu50(i,0),6);
//    v975str = ST::doubletostring(betaqu_l1_upper(i,0),6);
//
//    optionsp->out(vstr + ST::string(' ',15-vstr.length())
//                  + vmstr + ST::string(' ',15-vmstr.length())
//                  + vsdstr + ST::string(' ',15-vsdstr.length())
//                  + v025str + ST::string(' ',15-v025str.length())
//                  + v500str + ST::string(' ',15-v500str.length())
//                  + v975str + ST::string(' ',15-v975str.length()) + "\n");
//    }
//
//  optionsp->out("\n");


  // Datei-Ausgabe Ergebnisse
  ofstream ou(pathresults.strtochar());
  ou << "varname  pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
  nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
  for(i=0; i<beta.rows(); i++)
    {
    ou << "var." << (i+1) << "   ";
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


  outresults_indicator();
  outresults_t2();
  outresults_shrinkage();

  // Samples-Ausgebepfad
  bool nosamples = distrp->get_nosamples();

/*  if(!nosamples)
    {
    ST::string pathhelp1 = pathresults.substr(0,pathresults.length()-7)+"indicator_sample.raw";
    ST::string pathhelp2 = pathresults.substr(0,pathresults.length()-7)+"t2_sample.raw";
    ST::string pathhelp3 = pathresults.substr(0,pathresults.length()-7)+"shrinkage_sample.raw";

    fc_indicator.get_samples(pathhelp1);
    fc_t2.get_samples(pathhelp2);
    fc_shrinkage.get_samples(pathhelp3);


    optionsp->out("  Sampled indicator parameters are stored in file\n");
    optionsp->out("  " + pathhelp1 + "\n");
    optionsp->out("\n");
    optionsp->out("  Sampled t2 parameters are stored in file\n");
    optionsp->out("  " + pathhelp2 + "\n");
    optionsp->out("\n");
    optionsp->out("  Sampled w parameters are stored in file\n");
    optionsp->out("  " + pathhelp3 + "\n");
    optionsp->out("\n");
    optionsp->out("\n");
    }*/

  }

void FULLCOND_variance_nonp_vector_nigmix::get_samples(const ST::string & filename,const unsigned & step) const
  {
  FULLCOND::get_samples(filename, step);

  ST::string pathhelp1 = pathresults.substr(0,pathresults.length()-7)+"indicator_sample.raw";
  ST::string pathhelp2 = pathresults.substr(0,pathresults.length()-7)+"t2_sample.raw";
  ST::string pathhelp3 = pathresults.substr(0,pathresults.length()-7)+"shrinkage_sample.raw";

  optionsp->out(pathhelp1 + "\n");
  optionsp->out("\n");
  fc_indicator.get_samples(pathhelp1);

  optionsp->out(pathhelp2 + "\n");
  optionsp->out("\n");
  fc_t2.get_samples(pathhelp2);

  optionsp->out(pathhelp3 + "\n");
  optionsp->out("\n");
  fc_shrinkage.get_samples(pathhelp3);
  }

//______________________________________________________________________________
//
// FUNCTION: outresults_indicator
// TASK: - write results for indicator to output window and files
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector_nigmix::outresults_indicator(void)
  {

  fc_indicator.outresults();
  vector<ST::string> varnames_indicator = fc_indicator.get_datanames();

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


  // Bildschirm-Ausgabe
  ST::string vstr;
  ST::string vmstr;
  ST::string vsdstr;
  ST::string v025str;
  ST::string v500str;
  ST::string v975str;

  optionsp->out("\n");
  optionsp->out("  f_shrinkage_indicator \n");
  optionsp->out("\n");
  optionsp->out("\n");

  unsigned i;
  optionsp->out("  Variable     frequency");
  for(i=0; i<beta.rows(); i++)
    {
    vstr = "  " + varnames_indicator[i] + " ";
    vmstr = ST::doubletostring(fc_indicator.get_betamean(i,0),6);
//    vsdstr = ST::doubletostring((fc_indicator.get_betavar(i,0)<0.0?0.0:sqrt(fc_indicator.get_betavar(i,0))),6);
//    v025str = ST::doubletostring(fc_indicator.get_beta_lower1(i,0),6);
//    v500str = ST::doubletostring(fc_indicator.get_betaqu50(i,0),6);
//    v975str = ST::doubletostring(fc_indicator.get_beta_upper1(i,0),6);

    optionsp->out(vstr + ST::string(' ',15-vstr.length())
                  + vmstr + ST::string(' ',15-vmstr.length())
//                  + vsdstr + ST::string(' ',15-vsdstr.length())
//                  + v025str + ST::string(' ',15-v025str.length())
//                  + v500str + ST::string(' ',15-v500str.length())
//                  + v975str + ST::string(' ',15-v975str.length()) 
                  + "\n");
    }

  optionsp->out("\n");


  // Datei-Ausgabe Ergebnisse
  ST::string indicator_pathresults = pathresults.substr(0,pathresults.length()-7) + "indicator.res";

  ofstream ou(indicator_pathresults.strtochar());


//  ou << "varname  pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
//  nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
  ou << "varname  frequency" << endl;  
  for(i=0; i<beta.rows(); i++)
    {
    ou << varnames_indicator[i] << "   ";
    ou << fc_indicator.get_betamean(i,0) << "  "<< endl;
//    ou << (fc_indicator.get_betavar(i,0)<0.0?0.0:sqrt(fc_indicator.get_betavar(i,0))) << "  ";
//    ou << fc_indicator.get_beta_lower1(i,0) << "  ";
//    ou << fc_indicator.get_beta_lower2(i,0) << "  ";
//    ou << fc_indicator.get_betaqu50(i,0) << "  ";
//    ou << fc_indicator.get_beta_upper2(i,0) << "  ";
//    ou << fc_indicator.get_beta_upper1(i,0) << "  ";
//    ou << fc_indicator.get_betamin(i,0) << "  ";
//    ou << fc_indicator.get_betamax(i,0) << "  " << endl;
    }
  
  optionsp->out("  Results for the indicator parameter are also stored in file\n");
  optionsp->out("  " + indicator_pathresults + "\n");

  optionsp->out("\n");
  }


//______________________________________________________________________________
//
// FUNCTION: outresults_shrinkage
// TASK: - write results for t2 to output window and files
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector_nigmix::outresults_t2(void)
  {

  fc_t2.outresults();
  vector<ST::string> varnames_t2 = fc_t2.get_datanames();
  
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


  // Bildschirm-Ausgabe
  ST::string vstr;
  ST::string vmstr;
  ST::string vsdstr;
  ST::string v025str;
  ST::string v500str;
  ST::string v975str;

  optionsp->out("\n");
  optionsp->out("  f_shrinkage_t2 \n");
  optionsp->out("\n");
  optionsp->out("\n");

  unsigned i;

//  optionsp->out("  Variable     mean           Std. Dev.      " + l1 +"% quant.    median         "+ u2 +"% quant.");
//  for(i=0; i<beta.rows(); i++)
//    {
//    vstr = varnames_t2[i] + " ";
//    vmstr = ST::doubletostring(fc_t2.get_betamean(i,0),6);
//    vsdstr = ST::doubletostring((fc_t2.get_betavar(i,0)<0.0?0.0:sqrt(fc_t2.get_betavar(i,0))),6);
//    v025str = ST::doubletostring(fc_t2.get_beta_lower1(i,0),6);
//    v500str = ST::doubletostring(fc_t2.get_betaqu50(i,0),6);
//    v975str = ST::doubletostring(fc_t2.get_beta_upper1(i,0),6);
//
//    optionsp->out(vstr + ST::string(' ',15-vstr.length())
//                  + vmstr + ST::string(' ',15-vmstr.length())
//                  + vsdstr + ST::string(' ',15-vsdstr.length())
//                  + v025str + ST::string(' ',15-v025str.length())
//                  + v500str + ST::string(' ',15-v500str.length())
//                  + v975str + ST::string(' ',15-v975str.length()) + "\n");
//    }
//
//  optionsp->out("\n");


  // Datei-Ausgabe Ergebnisse
  ST::string t2_pathresults = pathresults.substr(0,pathresults.length()-7) + "t2.res";

  ofstream ou(t2_pathresults.strtochar());


  ou << "varname  pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
  nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
  for(i=0; i<beta.rows(); i++)
    {
    ou << varnames_t2[i] << "   ";
    ou << fc_t2.get_betamean(i,0) << "  ";
    ou << (fc_t2.get_betavar(i,0)<0.0?0.0:sqrt(fc_t2.get_betavar(i,0))) << "  ";
    ou << fc_t2.get_beta_lower1(i,0) << "  ";
    ou << fc_t2.get_beta_lower2(i,0) << "  ";
    ou << fc_t2.get_betaqu50(i,0) << "  ";
    ou << fc_t2.get_beta_upper2(i,0) << "  ";
    ou << fc_t2.get_beta_upper1(i,0) << "  ";
    ou << fc_t2.get_betamin(i,0) << "  ";
    ou << fc_t2.get_betamax(i,0) << "  " << endl;
    }
  
  optionsp->out("  Results for the t2 parameter are also stored in file\n");
  optionsp->out("  " + t2_pathresults + "\n");

  optionsp->out("\n");
  }


//______________________________________________________________________________
//
// FUNCTION: outresults_shrinkage
// TASK: - write results for shrinkaageparameter omega to output window and files
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector_nigmix::outresults_shrinkage(void)
  {

  fc_shrinkage.outresults();
  vector<ST::string> varname_shrinkage = fc_shrinkage.get_datanames();

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

  // Bildschirm-Ausgabe
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


  // Datei-Ausgabe Ergebnisse
  ST::string shrinkage_pathresults = pathresults.substr(0,pathresults.length()-7) + "shrinkage.res";

  ofstream ou(shrinkage_pathresults.strtochar());

  ou << "varname  pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
  nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
  ou << varname_shrinkage[0] << "  ";
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





