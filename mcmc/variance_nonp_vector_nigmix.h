
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (VARIANCENONP_VECTOR_NIG_INCLUDED)
#define VARIANCENONP_VECTOR_NIG_INCLUDED

#include"mcmc_const.h"

namespace MCMC
{


class __EXPORT_TYPE FULLCOND_variance_nonp_vector_nigmix : public FULLCOND
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
  FULLCOND fc_indicator;
  FULLCOND fc_t2;

  vector<double> indicatorstart;    //  Blocks of regression coefficients
  double v0;                        //  Hyperparameter for Shrinkageparameter
  double v1;                        //  Hyperparameter for Shrinkageparameter
  vector<double> t2start;           //  Blocks of regression coefficients
  double a_t2;                      //  Hyperparameter for Shrinkageparameter
  double b_t2;                      //  Hyperparameter for Shrinkageparameter
  bool omegafix;                    //  Mixingparameter fix

  datamatrix indicator;             // Matrix for 1. Varianceparameterkomponent: Indicators
  datamatrix t2;                    // Matrix for 1. Varianceparameterkomponent: t2
 
  double nigmixsum;                 //  sum(beta^2/tau^2) for update of scaleparameter

  vector<unsigned> cut;             //  Blocks of regression coefficients
                                  
  void outresults_shrinkage(void);  //  Function to write results of omega to output window and files
  void outresults_indicator(void);  //  Function to write results of indicator to output window and files
  void outresults_t2(void);         //  Function to write results of t2 to output window and files

  public:

  //____________________________________________________________________________
  //
  // DEFAULT CONSTRUCTOR
  //____________________________________________________________________________

  FULLCOND_variance_nonp_vector_nigmix(void) : FULLCOND()
    {
    }
    

  //____________________________________________________________________________
  //
  // CONSTRUCTOR
  //____________________________________________________________________________

  FULLCOND_variance_nonp_vector_nigmix(MCMCoptions * o, vector<FULLCOND_const*> & p,
                         DISTRIBUTION * d,const ST::string & ti,
                         const ST::string & fp, const ST::string & fr,
                         const vector<double> & ins, const double & vv0, const double & vv1,
                         const vector<double> & t2s, const double & at2, const double & bt2,
                         const double & omegastart, const bool & omf,
                         const vector<unsigned> & ct,
                         const unsigned & c);

  //____________________________________________________________________________
  //
  // COPY CONSTRUCTOR
  //____________________________________________________________________________

  FULLCOND_variance_nonp_vector_nigmix(const FULLCOND_variance_nonp_vector_nigmix & t);
  
  
  //____________________________________________________________________________
  //
  // OVERLOADED ASSIGNMENT OPERATOR
  //____________________________________________________________________________

  const FULLCOND_variance_nonp_vector_nigmix & operator=(const FULLCOND_variance_nonp_vector_nigmix & t);


  // Pointer auf das shrinkage-Parameter Fullcond-Objekt
  FULLCOND * get_shrinkagepointer();


  void get_samples(const ST::string & filename, const unsigned & step = 1) const;

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

  ~FULLCOND_variance_nonp_vector_nigmix() {}

  }; // end: class FULLCOND_variance_nonp_vector_nigmix



} // end: namespace MCMC

#endif

