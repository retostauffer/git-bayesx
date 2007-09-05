
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

  vector<double> tau;               // Varianceparameters
  vector<double> lambda;            // Inverse Varianceparameter: lambda=1/tau^2

  bool update_sigma2;

  unsigned column;                  //  Category for for fullcond if multivariate response

  ST::string pathresults;           //  File path for storing sampled parameters
  
  vector<FULLCOND_const *> Cp;
  
  DISTRIBUTION * distrp;

  FULLCOND fc_shrinkage;
  bool shrinkagefix;                //  Shrinkageparameter fix
  double a_shrinkagegamma;          //  Hyperparameter for Shrinkageparameter
  double b_shrinkagegamma;          //  Hyperparameter for Shrinkageparameter

  double ridgesum;                  //  sum(beta^2/tau^2)

  vector<unsigned> cut;             //  Blocks of regression coefficients
  bool is_ridge;                    //  The Components indicates if "true" the L2-penalty
                                    //  and if "false" the L1-penalty is used
                                    
  void outresults_shrinkage(void);  //  Function to write results to output window and files

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
                         DISTRIBUTION * d,const ST::string & ti,
                         const ST::string & fp, const ST::string & fr,
                         const double & shrinkage_start, const double & a_shrinkage_gamma,
                         const double & b_shrinkage_gamma, const bool & shrinkage_fix,
                         const bool & isridge, const vector<unsigned> & ct,
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


  // Pointer auf das shrinkage-Parameter Fullcond-Objekt
  FULLCOND * get_shrinkagepointer();

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

