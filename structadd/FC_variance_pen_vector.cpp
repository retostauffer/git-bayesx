
#include"FC_variance_pen_vector.h"



//------------------------------------------------------------------------------
//-- CLASS: FC_variance_pen_vector implementation of member functions ----------
//------------------------------------------------------------------------------


namespace MCMC
{

using randnumbers::rand_inv_gaussian;
using randnumbers::rand_gamma;


void FC_variance_pen_vector::add_variable(datamatrix & x,
                         vector<ST::string> & op,
                         vector<ST::string> & vn)
  {

  int f;
  double la;
  double shrink;
  bool sfix;
//  bool adsh;
  double a,b;
  double shw;

  f = op[4].strtodouble(la);
  f = op[29].strtodouble(shrink);
  if (op[30]=="true")
    sfix = true;
  else
    sfix = false;
//  if (op[32]=="true")
//    adsh = true;
//  else
//    adsh = false;
  f = op[5].strtodouble(a);
  f = op[6].strtodouble(b);
  f = op[31].strtodouble(shw);

  lambda.push_back(la);
  shrinkagestart.push_back(shrink);
  shrinkagefix.push_back(sfix);
//  adaptiveshrinkage.push_back(adsh);
  a_shrinkagegamma.push_back(a);
  b_shrinkagegamma.push_back(b);
  shrinkageweight.push_back(shw);
  }


//______________________________________________________________________________
//
// CONSTRUCTOR with Parameters
//
// mp                : pointer to MASTER_object
// o                 : pointer to GENERAL_options object
// p                 : pointer to FULLCOND_lin_pen object
// d                 : pointer to DISTRIBUTION object
// ti                : reference to title of the full conditional(for example "fixed effects")
// fp                : reference to file path for storing sampled parameters
// isr               : variable,indicates if L2- or L1- penalty is uesd
// fr                : reference to filename for storing results
// shrinkage_start   : Starting value for the shrinkageparameter
// a_shrinkage_gamma : Hyperparameter for gammaprior of shrinkageparameter
// b_shrinkagegamma  : Hyperparameter for gammaprior of shrinkageparameter
// shrinkage_fix     : Should the shrinkageparameter be fixed at the value "shrinkage_start"
// ct                : blocks of regression coefficients
// c                 : reference to responsecategory (important in the case of multivariate response)
//______________________________________________________________________________
//==============================================================================
// KOMMENTAR Susanne:
//-------------------
// Derzeit ist die scale-dependent Variante implementiert, d.h. die 
// Shrinkageprioris hängen in Gaussfall vom Skalenparameter sigma2 ab.
// Dies hat zur Folge, das zum update von sigma2 der Wert sum(beta^2/tau^2),
// der in lassosum und ridgesum berechnet wird, benötigt wird.
// Für den Fall das kein adaptives shrinkage gewählt wird, werden die
// ersten Einträge aus den Vektoroptionen für das Verfahren verwendet.
// Einentsprechender Hinweis sollte in das Handbuch.
//==============================================================================

FC_variance_pen_vector::FC_variance_pen_vector(MASTER_OBJ * mp,
                        GENERAL_OPTIONS * o, FC_linear_pen * p,
                        DISTR * d,const ST::string & ti,
                        const ST::string & fp, bool isr)
                        :FC(o,ti,1,1,fp)
  {

  is_ridge = isr;

  update_sigma2 = true;

  Cp = p;

  distrp = d;

  priorassumptions.push_back("\\\\");

  // vector<ST::string> vnames = FC_shrinkage.get_datanames();
  //Initialisieren der Betamatrizen für die Varianzan + Übergabe der Startwerte

  lassosum = 0;
  ridgesum = 0;
  
  is_fix = true;          // = shrinkagefix[1]
  is_adaptive = true;     // = adaptiveshrinkage[1]

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
  FC_shrinkage = t.FC_shrinkage;
  shrinkagefix = t.shrinkagefix;
  adaptiveshrinkage = t.adaptiveshrinkage;
  a_shrinkagegamma = t.a_shrinkagegamma;
  b_shrinkagegamma = t.b_shrinkagegamma;
  shrinkagestart = t.shrinkagestart;
  shrinkageweight = t.shrinkageweight;  
  lassosum = t.lassosum;
  ridgesum = t.ridgesum;
  is_ridge = t.is_ridge;
  is_fix = t.is_fix; 
  is_adaptive = t.is_adaptive;
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
  FC_shrinkage = t.FC_shrinkage;
  shrinkagefix = t.shrinkagefix;
  adaptiveshrinkage = t.adaptiveshrinkage;
  a_shrinkagegamma = t.a_shrinkagegamma;
  b_shrinkagegamma = t.b_shrinkagegamma;
  shrinkagestart = t.shrinkagestart;
  shrinkageweight = t.shrinkageweight;  
  lassosum = t.lassosum;
  ridgesum = t.ridgesum;
  is_ridge = t.is_ridge;
  is_fix = t.is_fix; 
  is_adaptive = t.is_adaptive;
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
  acceptance++;
  unsigned i;
  double helpshrinkage;
  
  
  // reset Variables for summs
  int nrpar = beta.rows();
  double rand_invgaussian = 0;
  double sumvariances = 0;
  double sumregcoeff = 0;
  lassosum = 0;
  ridgesum = 0;
  

//TEMP:BEGIN--------------------------------------------------------------------
// KOMMENTAR: Susanne 
// das kann rausgenommen werden sobald die shrinkageparameter startwerte
// an das FC_shrinkage übergeben wurden
FC_shrinkage = FC(optionsp,"",shrinkagefix.size(),1,samplepath + ".shrinkage");
unsigned k;
double * helpshrinkagep = FC_shrinkage.beta.getV();
for(k=0; k<nrpar; k++, helpshrinkagep++)
  {
  *helpshrinkagep = 1;
  }
//TEMP:END----------------------------------------------------------------------  

  // get current value of first shrinkagearameter 
  double * shrinkagep = FC_shrinkage.beta.getV();


  // get current value of first regressionparameter
  double * workbeta = Cp->beta.getV();

  // getcurrent value of sqrt(scale) parameter
  double sigma = sqrt(distrp->get_scale(column));

  // Gibbs-Update of varianceparameters tau^2 
  //-----------------------------------------
  if (is_ridge == 0 && is_adaptive == false)   // Lasso L1-Penalty
   {
   for(i=0; i<nrpar; i++, workbeta++, shrinkagep++)
    {
    if (*workbeta>0 && *shrinkagep>0)
      {
      rand_invgaussian = rand_inv_gaussian(sigma*(*shrinkagep)/(*workbeta * shrinkageweight[i]), (*shrinkagep)*(*shrinkagep)/(shrinkageweight[i]*shrinkageweight[i]));
      beta(i,0) = 1.0/rand_invgaussian;
      }
    if (*workbeta<0 && *shrinkagep>0)
      {
      rand_invgaussian = rand_inv_gaussian(-1.0*sigma*(*shrinkagep)/(*workbeta * shrinkageweight[i]), (*shrinkagep)*(*shrinkagep)/(shrinkageweight[i]*shrinkageweight[i]));
      beta(i,0) = 1.0/(rand_invgaussian);
      }
    if (*workbeta==0 || *shrinkagep<=0)
      {
      beta(i,0) = 1E-6;
      }
    lassosum = lassosum + ((*workbeta)*(*workbeta))/beta(i,0);                          // sum(beta^2/tau^2)
    sumvariances = sumvariances + beta(i,0)/(shrinkageweight[k] * shrinkageweight[k]);  // sum(tau^2/weights^2) of current variances
   }
  }
  
  if (is_ridge == 0 && is_adaptive == true)   // Lasso L1-Penalty adaptive
   {
   for(i=0; i<nrpar; i++, workbeta++, shrinkagep++)
    {
    if (*workbeta>0 && *shrinkagep>0)
      {
      rand_invgaussian = rand_inv_gaussian(sigma*(*shrinkagep)/(*workbeta), (*shrinkagep)*(*shrinkagep));
      beta(i,0) = 1.0/rand_invgaussian;
      }
    if (*workbeta<0 && *shrinkagep>0)
      {
      rand_invgaussian = rand_inv_gaussian(-1.0*sigma*(*shrinkagep)/(*workbeta), (*shrinkagep)*(*shrinkagep));
      beta(i,0) = 1.0/(rand_invgaussian);
      }
    if (*workbeta==0 || *shrinkagep<=0)
      {
      beta(i,0) = 1E-6;
      }
    lassosum = lassosum + ((*workbeta)*(*workbeta))/beta(i,0);  // sum(beta^2/tau^2)
    //sumvariances = sumvariances + beta(i,0);                    // sum(tau^2) of current variances
   }
  }

  if (is_ridge == 1 && is_adaptive == false)   // Ridge L2-penalty
    {
    for(i=0; i<nrpar; i++, workbeta++, shrinkagep++)
      {
      beta(i,0) = shrinkageweight[i]/(2*(*shrinkagep));
      ridgesum = ridgesum + ((*workbeta)*(*workbeta))/beta(i,0);  // sum(beta^2/tau^2)
      sumregcoeff = sumregcoeff + (*workbeta)*(*workbeta)/(shrinkageweight[i]);
      }
    }
  if (is_ridge == 1 && is_adaptive == true)   // Ridge L2-penalty adaptive
    {
    for(i=0; i<nrpar; i++, workbeta++, shrinkagep++)
      {
      beta(i,0) = 1/(2*(*shrinkagep));
      ridgesum = ridgesum + ((*workbeta)*(*workbeta))/beta(i,0);  // sum(beta^2/tau^2)
      //sumregcoeff = sumregcoeff + (*workbeta)*(*workbeta);
      }
    }

  // Gibbs-Update of Shrinkageparameter with Gammadistribution
  //----------------------------------------------------------
  if(is_fix==false)              
    {
    if(is_ridge == 0 && is_adaptive == false)            // Lasso L1-penalty
      {
      helpshrinkage = sqrt(rand_gamma(nrpar + a_shrinkagegamma[1], b_shrinkagegamma[1] + 0.5*sumvariances));
      for(i=0; i<nrpar; i++, shrinkagep++)
        {
        *shrinkagep = helpshrinkage;
        }
      }
    if(is_ridge == 0 && is_adaptive == true)            // Lasso L1-penalty adaptive
      {
      for(i=0; i<nrpar; i++, shrinkagep++)
        {
        *shrinkagep = sqrt(rand_gamma(1 + a_shrinkagegamma[i], b_shrinkagegamma[i] + 0.5*beta(i,0)));
        }
      }
      
      
    if(is_ridge == 1 && is_adaptive == false)            // Ridge L2-penalty
      {  
      helpshrinkage = rand_gamma(0.5*nrpar + a_shrinkagegamma[1], b_shrinkagegamma[1] + sumregcoeff/(sigma*sigma));
      for(i=0; i<nrpar; i++, shrinkagep++)
        {
        *shrinkagep = helpshrinkage;
        }
      }
    if(is_ridge == 1 && is_adaptive == true)            // Ridge L2-penalty adaptive
      {  
      for(i=0; i<nrpar; i++, shrinkagep++, workbeta++)
        {
        *shrinkagep = rand_gamma(0.5*1 + a_shrinkagegamma[i], b_shrinkagegamma[i] + ((*workbeta)*(*workbeta))/(sigma*sigma));
        }
      }
    }

//==============================================================================
// KOMMENTAR: Susanne2Stefan
//----------------------------
// Update des varianceparameter tau^2 -> muss an FC_linear_pen übergeben werden
// zum updaten der Designmatrix zur Berechnung der Erwartungswerts und der
// Varianz: betapen ~ NV(E,V)
// E = ((X'X + diag(1/tau_1^2,...,1/tau_p^2)^-1) *X'ytilde,  
// V = sigma2*(X'X + diag(1/tau_1^2,...,1/tau_p^2)^-1)
// Siehe auch \mcmc\mcmc_const.cpp : ...::compute_matrices(void)
//==============================================================================
  //double * varp = Cp->getvariancespointer();
  //for(i=0; i<nrpar; i++, varp++)
  //  {
  //  *varp = beta(i,0);
  //  }


//==============================================================================
// KOMMENTAR Susanne2Stefan :
//---------------------------
// lassosum und ridgesum braucht man zum updaten des Skalenparameters
// sigma2, je nch dem ob die prioris für die regularisierten Effekte
// von sigma2 abhängen, was sie derzeit tun, vgl. \mcmc\distribution.cpp Z.5722
// Die Werte müssen also dahin übergeben werden wo das sigma2 update stattfindet
//==============================================================================
  // Ergebnisse zu ridgesum in Distributionobjekt zurückschreiben
  //if (is_ridge == 0)                              // L1-penalty
  //  {
  //  distrp->update_lasso(lassosum);
  //  }
  //if (is_ridge == 1)                              // L2-penalty
  //  {
  //  distrp->update_ridge(ridgesum);
  //  }


  // Update Shrinkageparameter
  //FC_shrinkage.update();

  //FC::update();

  }


//______________________________________________________________________________
//
// FUNCTION: get_samples
// TASK: - write samples to files
//______________________________________________________________________________

void FC_variance_pen_vector::get_samples(const ST::string & filename,ofstream & outg) const
  {

  FC::get_samples(filename,outg);

  ST::string pathhelp = filename.substr(0,filename.length()-7)+"shrinkage_sample.raw";

  FC_shrinkage.get_samples(pathhelp,outg);

  }

//______________________________________________________________________________
//
// FUNCTION: outresults
// TASK: - write results for varianceparameters to output window and files
//______________________________________________________________________________

void FC_variance_pen_vector::outresults(ofstream & out_stata, ofstream & out_R,
                             const ST::string & pathresults)
  {
  // FC.outresults(out_stata,out_R,pathresults);
  // FC.outresults_help(out_stata,out_R,pathresults,datanames);


  optionsp->out("\n");
  optionsp->out("    Results for variances are also stored in file\n");
  optionsp->out("    " + pathresults + "\n");
  optionsp->out("\n");

  outresults_shrinkage();

  }

//______________________________________________________________________________
//
// FUNCTION: outresults_shrinkage
// TASK: - write results for shrinkageparameter to output window and files
//______________________________________________________________________________

void FC_variance_pen_vector::outresults_shrinkage(void)
  {
  
  if(is_adaptive==true)
    {
    // FC_shrinkage.outresults(out_stata,out_R,pathresults);
    // FC_shrinkage.outresults_help(out_stata,out_R,pathresults,datanames);
    }

  if(is_adaptive==false)
    {
    // KOMMENTAR Susanne:  Hie können prinzipiell die selben funktionen aufgerufen 
    // werden nur das dann nr=1 bzw nrpar=1 gesetzt werden muss, da sonst
    // nrpar mal das selbe rausgeschrieben wird
    
    // FC_shrinkage.outresults(out_stata,out_R,pathresults);
    // FC_shrinkage.outresults_help(out_stata,out_R,pathresults,datanames);

    }
    
  optionsp->out("\n");
  optionsp->out("    Results for the shrinkage parameter are also stored in file\n");
//  optionsp->out("    " + shrinkage_pathresults + "\n");
  optionsp->out("\n");

  }


//______________________________________________________________________________
//
// FUNCTION: outoptions
// TASK: - write options to output window
//______________________________________________________________________________

void FC_variance_pen_vector::outoptions(void)
  {

//TEMP:BEGIN--------------------------------------------------------------------
// KOMMENTAR Susanne:  FC_shrinkage kann dann hier wieder rausgenommen weden
   FC_shrinkage = FC(optionsp,"",shrinkagefix.size(),1,samplepath + ".shrinkage");
//TEMP:END----------------------------------------------------------------------

  int i;
  int nrpar = beta.rows();
  vector<ST::string> vnames;   // = FC_shrinkage.datanames.getV(); 

  for (i=0;i<nrpar;i++)
      vnames.push_back("_x_" + ST::inttostring(i));

  double * shrinkagep = FC_shrinkage.beta.getV();


  if(is_ridge == 0)
    {
    optionsp->out("  LINEAR EFFECTS WITH LASSO PENALTY: \n");
    }
  if(is_ridge == 1)
    {
    optionsp->out("  LINEAR EFFECTS WITH RIDGE PENALTY: \n");
    }
    
  if(is_adaptive==false)
    {
    optionsp->out("  Hyperparameter a for shrinkage: " +
                     ST::doubletostring(a_shrinkagegamma[1]) + "\n" );
    optionsp->out("  Hyperparameter b for shrinkage: " +
                     ST::doubletostring(b_shrinkagegamma[1]) + "\n" );
    if(is_fix==true)
      {
      optionsp->out("  Shrinkage is fixed at value = " +
                       ST::doubletostring(*shrinkagep) + "\n" );
      } 
    optionsp->out("\n");
    }

  if(is_adaptive==true)
    {
    for(i=0; i<nrpar; i++, shrinkagep++)
      {
      optionsp->out("  Hyperparameter a for shrinkage of " + vnames[i] + ": " +
                       ST::doubletostring(a_shrinkagegamma[i]) + "\n" );
      optionsp->out("  Hyperparameter b for shrinkage of " + vnames[i] + ": " +
                       ST::doubletostring(b_shrinkagegamma[i]) + "\n" );
      if(is_fix==true)
        {
        optionsp->out("  Shrinkage of " + vnames[i] + " is fixed at value = " +
                         ST::doubletostring(*shrinkagep) + "\n" );
        } 
      }
    optionsp->out("\n");
    
    }
  
  }


} // end: namespace MCMC





