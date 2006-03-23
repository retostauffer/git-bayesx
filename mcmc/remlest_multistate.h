#if !defined (remlest_multistate_INCLUDED)

#define remlest_multistate_INCLUDED

#include "statmat.h"
#include "fullcond.h"
#include "distribution.h"

#if defined(JAVA_OUTPUT_WINDOW)
#include"adminparse_basic.h"
#endif


//------------------------------------------------------------------------------
//------------------------ CLASS: remlest_multistate ---------------------------
//------------------------------------------------------------------------------

/*class remlest_multistate
  {

//------------------------------------------------------------------------------
//----------------------------- Variables --------------------------------------
//------------------------------------------------------------------------------

  private:

  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adminb_p;
  #endif

  ostream * logout;                         // Pointer to filestream for writing
                                            // output

  vector<MCMC::FULLCOND*> fullcond;

  ST::string respfamily;                         // distribution of the response

  ST::string outfile;
  int maxit;
  double lowerlim;
  double eps;
  double maxchange;

  unsigned nrtransitions;

  statmatrix<double> X;                         // fixed effects
  statmatrix<double> Z;                         // random effects

  vector<unsigned> xcut;                        // partition of X
  vector<unsigned> zcut;                        // partition of Z
  vector<unsigned> xcuttrans;            // partition of X corresponding to different transitions
  vector<unsigned> zcuttrans;            // partition of Z corresponding to different transitions

  statmatrix<double> beta;                          // regression coefficients
  statmatrix<double> theta;                         // variance parameters

  public:

//------------------------------------------------------------------------------
//----------------------------- Constructor ------------------------------------
//------------------------------------------------------------------------------


  // DEFAULT CONSTRUCTOR

  remlest_multistate(void) {}

  // CONSTRUCTOR 1: Initialize from fullcond-objects

  remlest_multistate(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adb,
  #endif
  vector<MCMC::FULLCOND*> & fc,datamatrix & re,
          const ST::string & family, const ST::string & ofile,
          const int & maxiter, const double & lowerlimit, const double & epsi,
          const double & maxch, const vector<unsigned> & nrfullconds,
          const datamatrix & weight, ostream * lo=&cout);

//------------------------------------------------------------------------------
//----------------------------- REML estimation --------------------------------
//------------------------------------------------------------------------------


  // Function: estimate
  // Task: Perform REML-estimation of multi state models
  //       returns true if an error or user break occured

  bool estimate(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight);


//------------------------------------------------------------------------------
//----------------------------- Object description -----------------------------
//------------------------------------------------------------------------------

  void outoptions();

//------------------------------------------------------------------------------
//----------------------------- Writing results --------------------------------
//------------------------------------------------------------------------------

  void make_plots(ofstream & outtex,ST::string path_batch,
                  ST::string path_splus);

  void make_model(ofstream & outtex, const ST::string & rname);

  void make_predictor(ofstream & outtex);

  void make_prior(ofstream & outtex);

  void make_options(ofstream & outtex);

  void make_fixed_table(ofstream & outtex);

  void make_graphics(const ST::string & title,
                     const ST::string & path_batch,
                     const ST::string & path_tex,
                     const ST::string & path_splus,
                     const ST::string & rname,
                     const bool & dispers);

  bool check_pause();

  // FUNCTION: out
  // TASK: writes results to outputstream or
  //       in Output window if BORLAND_OUTPUT_WINDOW is defined

  void out(const ST::string & s,bool thick=false,bool italic = false,
           unsigned size = 12,int r=0,int g=0, int b=0);

  void outerror(const ST::string & s);

  };*/

#endif

