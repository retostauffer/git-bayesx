
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FChrandomVARIANCEINCLUDED)

#define FChrandomVARIANCEINCLUDED

#include"../values.h"
#include<fstream.h>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include"FC_nonp_variance.h"
#include"design.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_hrandom_variance -----------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_hrandom_variance  : public FC_nonp_variance
  {

  protected:

  DISTR * likepRE;

  bool mult;

  double compute_quadform(void); 

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_hrandom_variance(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_hrandom_variance(GENERAL_OPTIONS * o,DISTR * lp, DISTR * lpRE,
                      const ST::string & t, const ST::string & fp,DESIGN * dp,
                      FC_nonp * FCn,vector<ST::string> & op);

  // COPY CONSTRUCTOR

  FC_hrandom_variance(const FC_hrandom_variance & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_hrandom_variance & operator=(const FC_hrandom_variance & m);

  // DESTRUCTOR

  ~FC_hrandom_variance()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

  void transform_beta(void);

  void read_options(vector<ST::string> & op);

  };


} // end: namespace MCMC

#endif


