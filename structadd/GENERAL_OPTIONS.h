
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (GENERALOPTIONS)

#define GENERALOPTIONS

#include<fstream>
#include<vector>
#include"clstring.h"

#if defined(JAVA_OUTPUT_WINDOW)
#include"adminparse_basic.h"
#endif

using std::cout;

namespace MCMC
{


//------------------------------------------------------------------------------
//--------------------------- CLASS: GENERAL_OPTIONS ---------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE GENERAL_OPTIONS
  {



  protected:


  ostream * logout;               // Pointer to filestream for writing
                                  // output

  public:

  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adminb_p;
  #endif

  unsigned iterations;            // Number of iterations
  unsigned burnin;                // Number of burnin iterations
  unsigned step;                  // Thinning parameter
                                  // (step = 100: every 100th sampled parameter
                                  //  will be stored and used for computing
                                  //  charkteristics of the posterior)

  double level1;                  // level1 for credible intervals
                                  // lower1 - upper2
  double level2;                  // level2 for credible intervals
                                  // lower2 - upper1

  double lower1;                  
  double lower2;

  double upper1;
  double upper2;


  unsigned nrbetween;
  unsigned nrout;

  unsigned nriter;                // current iteration number
  unsigned samplesize;            // current number of parameters used to
                                  // compute characteristics of the
                                  // posterior


  // DEFAULT CONSTRUCTOR
  // Defines:
  // iterations = 22000
  // burnin = 2000
  // step = 20
  // nrbetween = 1000    (may be changes via set_nrbetween)
  // nrout = 100         (may be changes via set_nrout)
  // nriter = 0
  // samplesize = 0
  // logout = &cout

  GENERAL_OPTIONS(void);

  // CONSTRUCTOR
  // Defines
  // iterations = it
  // burnin = bu
  // step = st
  // nrbetween = (iterations-burnin)/10
  // if (nrbetween < 1000) nrbetween = 1000
  // nrout = 100
  // nriter = 0
  // samplesize = 0
  // logout = lo

  #if defined(JAVA_OUTPUT_WINDOW)
  GENERAL_OPTIONS(administrator_basic * abp,
              const unsigned & it,const unsigned & bu,const unsigned & st,
              ostream * lo=&cout,const double & l1=95,const double & l2=80);
  #else
  GENERAL_OPTIONS(const unsigned & it,const unsigned & bu,const unsigned & st,
                  ostream * lo=&cout,const double & l1=95,const double & l2=80);
  #endif

  // COPY CONSTRUCTOR

  GENERAL_OPTIONS(const GENERAL_OPTIONS & o);

  // OVERLOADED ASSIGNMENT OPERATOR

  const GENERAL_OPTIONS & operator=(const GENERAL_OPTIONS & o);

  // FUNCTION: out
  // TASK: writes results to outputstream or
  //       in Output window if BORLAND_OUTPUT_WINDOW is defined

  void out(const ST::string & s,bool thick=false,bool italic = false,
           unsigned size = 12,int r=0,int g=0, int b=0);

  unsigned compute_samplesize(void);

  // FUNCTION: reset
  // TASK: resets the current number of iterations and the samplesize

  void reset(void);

  void update(void);

  void outoptions(void);

  void set_level1(double l1);

  void set_level2(double l2);

  // DESTRUCTOR

  ~GENERAL_OPTIONS() {}


  };

//------------------------------------------------------------------------------
//------------------------ End: CLASS GENERAL_OPTIONS --------------------------
//------------------------------------------------------------------------------

} // end: namespace MCMC

#endif
