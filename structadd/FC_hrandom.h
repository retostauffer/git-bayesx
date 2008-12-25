
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE  __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FChrandomINCLUDED)

#define FChrandomINCLUDED

#include"../values.h"
#include<fstream>
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


enum hrandomtype {additive,mult,multexp};

class __EXPORT_TYPE FC_hrandom  : public FC_nonp
  {

  protected:

  DISTR * likep_RE;

  FC FCrcoeff;

  hrandomtype rtype;

  datamatrix response_o;
  datamatrix linpred_o;

  void set_rcoeff(void);

  void update_linpred_multexp(void);
  void update_response_multexp(void);

  bool posteriormode_multexp(void);
  bool posteriormode_additive(void);

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_hrandom(void);

  // CONSTRUCTOR

  FC_hrandom(GENERAL_OPTIONS * o,DISTR * lp, DISTR * lp_RE,const ST::string & t,
           const ST::string & fp, const ST::string & fp2, DESIGN * dp,
           vector<ST::string> & op,vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  FC_hrandom(const FC_hrandom & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_hrandom & operator=(const FC_hrandom & m);

  // DESTRUCTOR

  ~FC_hrandom()
    {
    }

  void update_linpred(int & begin, int & end, double  & value);

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  void update_IWLS(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

  void transform_beta(void);

    // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(const ST::string & pathresults);


  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  };


} // end: namespace MCMC

#endif


