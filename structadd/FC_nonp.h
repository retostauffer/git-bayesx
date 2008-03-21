
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FCNONPINCLUDED)

#define FCNONPINCLUDED

#include"../values.h"
#include<fstream.h>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include"FC.h"
#include"design.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_nonp -----------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_nonp  : public FC
  {

  protected:

  DISTR * likep;                             // Pointer to DISTR obejct
  DESIGN * designp;                          // Pointer to design object

  datamatrix param;                          // Parameters


  datamatrix partres;                        // stores partial residuals

  datamatrix betahelp;

  void centerparam(void);

  public:

  double lambda;
  double tau2;


//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_nonp(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_nonp(GENERAL_OPTIONS * o,DISTR * lp, const ST::string & t,
           const ST::string & fp,DESIGN * dp);

  // COPY CONSTRUCTOR

  FC_nonp(const FC_nonp & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_nonp & operator=(const FC_nonp & m);

  // DESTRUCTOR

  ~FC_nonp()
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

  void outoptions(void)
    {
    }

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(const ST::string & pathresults);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void);

  };


} // end: namespace MCMC

#endif


