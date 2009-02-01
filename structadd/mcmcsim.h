
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (MCMCsim_INCLUDED)

#define MCMCsim_INCLUDED

#include"GENERAL_OPTIONS.h"
#include"distr.h"
#include"FC.h"

namespace MCMC
{


class __EXPORT_TYPE equation
  {

  protected:

  public:

  int hlevel;
  int equationnr;
  ST::string equationtype;

  unsigned nrfc;

  ST::string header;
  ST::string paths;

  DISTR * distrp;
  ST::string pathd;

  vector<FC*> FCpointer;
  vector<ST::string> FCpaths;

  // DEFAULT CONSTRUCTOR

  equation(void);

  // CONSTRUCTOR1

  equation(int enr, int hl,ST::string t);

  // CONSTRUCTOR2

  equation(const ST::string & h, DISTR * dp, const vector<FC*> fcp,
           const ST::string & pd, const vector<ST::string> & ps);

  // COPY CONSTRUCTOR

  equation(const equation & s);

  // OVERLOADED ASSIGNMENT CONSTRUCTOR

  const equation & operator=(const equation & s);

  void add_FC(FC * FCp,const ST::string & p);

  // DESTRUCTOR

  ~equation() {}

  };


class __EXPORT_TYPE MCMCsim
  {

  protected:

  GENERAL_OPTIONS * genoptions;

  vector<equation> equations;

  unsigned maxiterations;             // for posteriormode, maximum number of
                                      // iterations, default = 1000

  public:

  // DEFAULT CONSTRUCTOR

  MCMCsim(void)
    {
    }

  // CONSTRUCTOR
  // TASK: initializes the MCMC simulation object with general MCMC options 'go'
  //       a vector of equations 'equ'

  MCMCsim(GENERAL_OPTIONS * go,vector<equation> & equ);

  // COPY CONSTRUCTOR

  MCMCsim(const MCMCsim & s);

  // OVERLOADED ASSIGNMENT CONSTRUCTOR

  const MCMCsim & operator=(const MCMCsim & s);

  // FUNCTION: simulate
  // TASK: runs a MCMC simulation
  //       returns true, if simulation error or user break occured

  bool simulate(const int & seed, const bool & computemode=true);

  bool posteriormode(const bool & presim=false);

  void out_effects(const vector<ST::string> & paths);


  // FUNCTION: autocorr
  // TASK: computes autocorrelations for all samples parameters
  //      (i.e. for all beta matrices) and stores the result in datamatrix
  //      'cmat'

  void autocorr(const unsigned & lag,datamatrix & cmat);

  // FUNCTION: autocorr
  // TASK: computes autocorrelations for all samples parameters
  //      (i.e. for all beta matrices) and stores the result in file 'path'

  void autocorr(const unsigned & lag,const ST::string & path);

  // FUNCTION: compute_nrpar
  // TASK: computes the total number of parameters

  unsigned compute_nrpar(void);

  // FUNCTION: get_samples
  // TASK: stores sampled parameters of all full conditionals in ASCII format
  //       for each full conditional one file will be created with filename
  //       'path' + title of the full conditional + "_sample.raw"

  void get_samples(
  #if defined(JAVA_OUTPUT_WINDOW)
  vector<ST::string> & newc,
  #endif
  const ST::string & path,const unsigned & step=1);

  // DESTRUCTOR

  ~MCMCsim() {}

  };



} // end: namespace MCMC

#endif
