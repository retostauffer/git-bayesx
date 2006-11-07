
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (MCMC_RIDGE_INCLUDED)

#define MCMC_RIDGE_INCLUDED

#include "fullcond.h"
#include "mcmc_const.h"


namespace MCMC
{
using randnumbers::rand_inv_gaussian;          
          
          
class __EXPORT_TYPE FULLCOND_ridge : public FULLCOND
  {

  protected:     
  
  vector<double> variances;    // vector of variances for the ridge penalty 
                               // (tau^2=:variances)
                                                                                     
  double lasso;                // lassoparameter (lambda=:lasso) 
  
  DISTRIBUTION * likep;
  
  datamatrix linold;           // linold=data*beta
  datamatrix mu1;              // mu1=response-predictorpart  
  datamatrix XX;               // XX=X'X
  datamatrix X1;               // X1=(X'WX)^-0.5
  datamatrix X2;               // X2=(X'WX)^-1X'W

//-Temorär ---------------------------------------------------------------------
  //  Hyperparameter für Lasso & Outputmatrix für Varianzschätzer, Lassoschätzer
  double hypr;
  double hyps;
  datamatrix estimVariances;
  datamatrix estimLasso;
//-Temorär Ende-----------------------------------------------------------------


  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------
  //____________________________________________________________________________
  //
  // DEFAULT CONSTRUCTOR
  //____________________________________________________________________________

  FULLCOND_ridge(void) : FULLCOND()
    {
    }
    
  //____________________________________________________________________________
  //
  // CONSTRUCTOR
  //____________________________________________________________________________

  FULLCOND_ridge(MCMCoptions * o, DISTRIBUTION * dp, const datamatrix & d,
                 const ST::string & t, const ST::string & fs,
                 const ST::string & fr, const vector<double> & vars,
                 const unsigned & c);

  //____________________________________________________________________________
  //
  // COPY CONSTRUCTOR
  //____________________________________________________________________________

  FULLCOND_ridge(const FULLCOND_ridge & m);

  //____________________________________________________________________________
  //
  // OVERLOADED ASSIGNMENT OPERATOR
  //____________________________________________________________________________

  const FULLCOND_ridge & operator=(const FULLCOND_ridge & m);

  //____________________________________________________________________________
  //
  // DESTRUCTOR
  //____________________________________________________________________________

  ~FULLCOND_ridge()
    {
    }


//-------------------------- UPDATE and related methods-------------------------
  //____________________________________________________________________________
  //
  // FUNCTION: update
  // TASK:     - stores sampled parameters in file 'samplepath'
  //           - storing order: first row, second row, ...
  //____________________________________________________________________________

  void update(void);

  //____________________________________________________________________________
  //
  // FUNCTION: outresults
  // TASK:     - write results to output window and files
  //____________________________________________________________________________

  void outresults(void);

  };


} // end: namespace MCMC

#endif
