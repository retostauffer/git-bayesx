
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FCNONPnongaussianINCLUDED)

#define FCNONPnongaussianINCLUDED

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
//----------------------- CLASS: FC_nonp_nongaussian ---------------------------
//------------------------------------------------------------------------------



class __EXPORT_TYPE FC_nonp_nongaussian  : public FC_nonp
  {



  protected:


  public:


//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_nonp_nongaussian(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_nonp_nongaussian(GENERAL_OPTIONS * o,DISTR * lp, const ST::string & t,
           const ST::string & fp,DESIGN * dp,vector<ST::string> & op,
             vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  FC_nonp_nongaussian(const FC_nonp_nongaussian & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_nonp_nongaussian & operator=(const FC_nonp_nongaussian & m);

  // DESTRUCTOR

  ~FC_nonp_nongaussian()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  void update_IWLS(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void outoptions(void);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void);

  };


} // end: namespace MCMC

#endif


