
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (DISTR_INCLUDED)
#define DISTR_INCLUDED

#include"statmat.h"
#include"random.h"
#include"GENERAL_OPTIONS.h"
#include"fc.h"


namespace MCMC
{

using randnumbers::rand_invgamma;
using randnumbers::rand_normal;
using randnumbers::uniform;
using randnumbers::trunc_normal;
using randnumbers::trunc_normal2;
using randnumbers::truncnormal;
using randnumbers::kssample;
using randnumbers::rand_gamma;



class __EXPORT_TYPE DISTR
  {

  bool check_workingweights_one(void);

  protected:


  GENERAL_OPTIONS * optionsp;         // pointer to general MCMC options object


  public:

  ST::string family;              // name of the distribution

  unsigned nrobs;                 // Number of observations

  datamatrix response;            // Response
  ST::string responsename;        // Name of the response


  datamatrix weight;              // Weightvariable for weighted regression
  ST::string weightname;          // Name of the weightvariable

  datamatrix workingweight;       // Working weight
  bool changingweight;
  bool weights_one;

  datamatrix partres;             // partial residual


  datamatrix linearpred1;          // Linear predictor
  datamatrix linearpred2;          // Proposed linear predictor
  int linpred_current;

  void swap_linearpred(void);


  double trmult;                   // multiplicative constant with which
                                   // the response has been transformed


//------------------------------------------------------------------------------
//--------------------------- CONSTRUCTORS -------------------------------------
//------------------------------------------------------------------------------

  // DEFAULT CONSTRUCTOR

  DISTR(void)
    {
    }

  // CONSTRUCTOR1
  // TASK: initializes data
  //       response = r
  //       weight = w
  //       nrobs = r.rows()

  DISTR(GENERAL_OPTIONS * o,const datamatrix & r,
               const datamatrix & w=datamatrix());


  // COPY CONSTRUCTOR

  DISTR(const DISTR & d);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR & operator=(const DISTR & d);

  // DESTRUCTOR

  ~DISTR() {}

  //----------------------------------------------------------------------------
  //------------------------------ WRITING OPTIONS -----------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: outoptions
  // TASK: writing options

  virtual void outoptions(void);


  //----------------------------------------------------------------------------
  //----------------------- COMPUTING THE LOGLIKELIHOOD ------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: loglikelihood
  // TASK: computes the loglikelihood for a single observation

  virtual double loglikelihood(double * res,double * lin,double * weight) const
    {
    return 0;
    }

  // FUNCTION: loglikelihood
  // TASK: computes the complete loglikelihood for all observations

  double loglikelihood(const bool & current=true) const;

  // FUNCTION: loglikelihood
  // TASK: computes the loglikelihood for observations between 'beg' and 'end'

  double loglikelihood(const unsigned & beg,const unsigned & end,
                       const statmatrix<int> & index,
                       const bool & current=true) ;

  //----------------------------------------------------------------------------
  //----------------------- ACCESS TO SCALE PARAMETER --------------------------
  //----------------------------------------------------------------------------

  virtual double get_scale(void);

  //----------------------------------------------------------------------------
  //----------------------- POSTERIORMODE FUNCTIONS ----------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  virtual bool posteriormode(void)
    {
    return true;
    }

  //----------------------------------------------------------------------------
  //--------------------------- UPDATE FUNCTIONS -------------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: update
  // TASK: base function for inherited classes,
  //       should update the scale parameter
  //       the base function updates the estimated mean and variance
  //       of the scale parameter only

  virtual void update(void);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: outresults
  // TASK: writes estimation results for the scale parameter
  //       estimated mean and variance

  virtual void outresults(ST::string pathresults="");

  // FUNCTION: reset
  // TASK: resets linpred (all values to 0)

  void reset(void);


  }; // end: class DISTR



//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_gaussian ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gaussian : public DISTR
  {

  protected:

  double a_invgamma;                    // hyperparameter a (for the inverse
                                        // gamma distribution of the scale
                                        // parameter, i.e. sigma^2
  double b_invgamma;                    // hyperparameter b

  FC FCsigma2;

  // FUNCTION: standardise
  // TASK: standardises the response and the offset
  //       sets scalesave.transform = trmult*trmult (!!!)

  void standardise(void);


  public:

  double sigma2;


   // DEFAULT CONSTRUCTOR

   DISTR_gaussian(void) : DISTR()
     {
     }

   // CONSTRUCTOR1
   // a_invgamma = a
   // b_invgamma = b

   DISTR_gaussian(const double & a,const double & b,GENERAL_OPTIONS * o,
                  const datamatrix & r,const ST::string & ps,
                         const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_gaussian(const DISTR_gaussian & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_gaussian & operator=(const DISTR_gaussian & nd);

   // DESTRUCTOR

   ~DISTR_gaussian() {}

   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * res,
                       double * lin,
                       double * w) const;


  virtual double get_scale(void);

  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  bool posteriormode(void);

  void outresults(ST::string pathresults="");


  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_gaussian_re -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gaussian_re : public DISTR_gaussian
  {

  protected:

  DISTR * distrp;

  public:

   // DEFAULT CONSTRUCTOR

   DISTR_gaussian_re(void) : DISTR_gaussian()
     {
     }

   // CONSTRUCTOR1

   DISTR_gaussian_re(GENERAL_OPTIONS * o,DISTR * dp,
                  const datamatrix & r,const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_gaussian_re(const DISTR_gaussian_re & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_gaussian_re & operator=(const DISTR_gaussian_re & nd);

   // DESTRUCTOR

   ~DISTR_gaussian_re() {}

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  bool posteriormode(void);


  };


} // end: namespace MCMC


#endif
