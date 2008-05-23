
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FCNONPVARIANCEINCLUDED)

#define FCNONPVARIANCEINCLUDED

#include"../values.h"
#include<fstream.h>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include"FC_nonp.h"
#include"design.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_nonp_variance --------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_nonp_variance  : public FC
  {

  protected:

  FC_nonp * FCnonpp;                         // Pointer to corresponding
                                             // FC_nonp object  
  DISTR * likep;                             // Pointer to DISTR obejct
  DESIGN * designp;                          // Pointer to design object

  double a_invgamma;
  double b_invgamma;
  double lambdastart;

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_nonp_variance(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_nonp_variance(GENERAL_OPTIONS * o,DISTR * lp, const ST::string & t,
           const ST::string & fp,DESIGN * dp,FC_nonp * FCn,
           vector<ST::string> & op,vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  FC_nonp_variance(const FC_nonp_variance & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_nonp_variance & operator=(const FC_nonp_variance & m);

  // DESTRUCTOR

  ~FC_nonp_variance()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void outoptions(void);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(const ST::string & pathresults);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  virtual void transform_beta(void);

  };


} // end: namespace MCMC

#endif


