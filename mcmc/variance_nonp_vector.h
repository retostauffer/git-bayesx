
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (VARIANCENONP_VECTOR_INCLUDED)
#define VARIANCENONP_VECTOR_INCLUDED

#include"mcmc_const.h"

namespace MCMC
{


class __EXPORT_TYPE FULLCOND_variance_nonp_vector : public FULLCOND
  {

  protected:

  vector<double> tau;           //  Wozu sind diese Vektoren (löschen veursacht Probleme in BayesX)!!!!!!!!!!!!!
  vector<double> lambda;        //  Wozu sind diese Vektoren (löschen veursacht Probleme in BayesX)!!!!!!!!!!!!!!

  bool update_sigma2;

  unsigned column;              //

  ST::string pathresults;       //  file path for storing sampled parameters
  
  vector<FULLCOND_const *> Cp;
  
  DISTRIBUTION * distrp;

  vector<double> a_invgamma;    //  Hyperparameters, not used yet
  vector<double> b_invgamma;    //  Hyperparameters, not used yet

  FULLCOND fc_lasso;
  bool lassofix;                //  Lassoparameter fix
  double a_lassogamma;          //  Hyperparameter for lasso
  double b_lassogamma;          //  Hyperparameter for lasso
  double lassomin;              //  Lower interval limit
  double lassomax;              //  Upper interval limit
  int lassogrid;                //  Interval gridpoints
  double ridgesum;               //  sum(beta^2/variances^2)

  vector<unsigned> cut;         // blocks of regression coefficients

  void outresults_lasso(void);  //  Function to write results to output window and files

  public:

  //____________________________________________________________________________
  //
  // DEFAULT CONSTRUCTOR
  //____________________________________________________________________________

  FULLCOND_variance_nonp_vector(void) : FULLCOND()
    {
    }
    

  //____________________________________________________________________________
  //
  // CONSTRUCTOR
  //____________________________________________________________________________

  FULLCOND_variance_nonp_vector(MCMCoptions * o, vector<FULLCOND_const*> & p,
                         DISTRIBUTION * d, const vector<double> & a,
                         const vector<double> & b, const ST::string & ti,
                         const ST::string & fp, const ST::string & fr,
                         const double & lasso_start, const double & a_lasso_gamma,
                         const double & b_lasso_gamma, const bool & shrinkage_fix,
                         const vector<bool> & isridge, const vector<unsigned> & ct,
                         const unsigned & c);
                         
  //____________________________________________________________________________
  //
  // COPY CONSTRUCTOR
  //____________________________________________________________________________

  FULLCOND_variance_nonp_vector(const FULLCOND_variance_nonp_vector & t);
  
  
  //____________________________________________________________________________
  //
  // OVERLOADED ASSIGNMENT OPERATOR
  //____________________________________________________________________________

  const FULLCOND_variance_nonp_vector & operator=(const FULLCOND_variance_nonp_vector & t);


  // Pointer auf das lasso-Parameter Fullcond-Objekt
  FULLCOND * get_lassopointer();

  //____________________________________________________________________________
  //
  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //____________________________________________________________________________

  void update(void);
  
  
  //____________________________________________________________________________
  //
  // FUNCTION: outresults
  // TASK: - write results to output window and files
  //____________________________________________________________________________

  void outresults(void);
  
  
  //____________________________________________________________________________
  //
  // FUNCTION: outoptions
  // TASK: - write options to output window
  //____________________________________________________________________________

  void outoptions(void);
  
  
  //____________________________________________________________________________
  //
  // FUNCTION: reset
  // TASK: resets all parameters
  //____________________________________________________________________________

  void reset(void)
    {
    FULLCOND::reset();
    setbeta(beta.rows(),1,0.1);
    }
    
    
  //____________________________________________________________________________
  //
  // DESTRUCTOR
  //____________________________________________________________________________

  ~FULLCOND_variance_nonp_vector() {}

  }; // end: class FULLCOND_variance_nonp_vector



} // end: namespace MCMC

#endif

