
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FChrandomINCLUDED)

#define FChrandomINCLUDED

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
//--------------------------- CLASS: FC_hrandom --------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_hrandom  : public FC_nonp
  {

  protected:

  DISTR * likep_RE;

  FC FCrcoeff;

  void set_response(void);

  void set_rcoeff(void);

  public:



//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_hrandom(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_hrandom(GENERAL_OPTIONS * o,DISTR * lp, DISTR * lp_RE,const ST::string & t,
           const ST::string & fp, const ST::string & fp2, DESIGN * dp);

  // COPY CONSTRUCTOR

  FC_hrandom(const FC_hrandom & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_hrandom & operator=(const FC_hrandom & m);

  // DESTRUCTOR

  ~FC_hrandom()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

    // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(const ST::string & pathresults);


  };


} // end: namespace MCMC

#endif


