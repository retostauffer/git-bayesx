#if !defined (remlest_multi_INCLUDED)

#define remlest_multi_INCLUDED

#include <statmat.h>
#include <fullcond.h>

#if defined(JAVA_OUTPUT_WINDOW)
#include<adminparse_basic.h>
#endif


//------------------------------------------------------------------------------
//------------------------ CLASS: remlest_multinomial --------------------------
//------------------------------------------------------------------------------

class remlest_multinomial
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

  double refcat;                                   // Referenz-Kategorie
  unsigned nrcat;                                    // Anzahl Kategorien
  unsigned nrcat2;                                   // Anzahl nicht redundanter Kategorien
  statmatrix<double> cats;                         // Mögliche Kategorien

  unsigned nrobs;
  unsigned nrobspos;

  unsigned partialnrpar;                        // Parameter pro Kategorie
  unsigned partialnrfixed;                      // fixe Effekte pro Kategorien
  unsigned partialnrrandom;
  unsigned totalnrfixed;                        // Anzahl fixer Effekte für alle Kategorien!
  unsigned totalnrpar;                          // Parameter in allen Kategorien

  unsigned partialvar;                          // Varianzparameter pro Kategorie

  statmatrix<double> X;                         // fixed effects
  statmatrix<double> Z;                         // random effects

  vector<unsigned> xcut;                        // partition of X
  vector<unsigned> zcut;                        // partition of Z

  vector<unsigned> xcutbeta;                    // partition of fixed part in beta
  vector<unsigned> zcutbeta;                    // partition of random part in beta

  statmatrix<double> beta;                          // regression coefficients
  statmatrix<double> theta;                         // variance parameters

  public:

//------------------------------------------------------------------------------
//----------------------------- Constructor ------------------------------------
//------------------------------------------------------------------------------


  // DEFAULT CONSTRUCTOR

  remlest_multinomial(void) {}

  // Initialize from fullcond-objects

  remlest_multinomial(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adb,
  #endif
  vector<MCMC::FULLCOND*> & fc,datamatrix & re,
          const ST::string & family, const ST::string & ofile,
          const int & maxiter, const double & lowerlimit, const double & epsi,
          const double & maxch, const datamatrix & categories,
          const datamatrix & weight, ostream * lo=&cout);

//------------------------------------------------------------------------------
//----------------------------- REML estimation --------------------------------
//------------------------------------------------------------------------------


  // Function: estimate
  // Task: Perform REML-estimation with nonparametric terms

  bool estimate(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight);

  // Function: estimate_glm
  // Task: Perform REML-estimation without nonparametric terms

  bool estimate_glm(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight);

//------------------------------------------------------------------------------
//------------- Weights, expectation, linear predictor, etc --------------------
//------------------------------------------------------------------------------

  // FUNCTION: compute_respind
  // TASK: Computes the indicator version of the response

  void compute_respind(const datamatrix & re, datamatrix & respind);

  // FUNCTION: compute_weights
  // TASK: Computes mu, weights and the working observations

  void compute_weights(datamatrix & mu, datamatrix & workweights,
                       datamatrix & worky, datamatrix & eta,
                       datamatrix & respind, const datamatrix & weights);

  // FUNCTION: compute_eta
  // TASK: Computes the linear predictor X*beta

  void compute_eta(datamatrix & eta);

  // FUNCTION: compute_eta
  // TASK: Computes the linear predictor X*beta+Z*b

  void compute_eta2(datamatrix & eta);

  // FUNCTION: compute_sscp
  // TASK: Computes the SSCP matrix X'WX

  void compute_sscp(datamatrix & H, datamatrix & workweight);

  // FUNCTION: compute_sscp2
  // TASK: Computes the weighted SSCP matrix for X and Z

  void compute_sscp2(datamatrix & H, datamatrix & workweight);

  // FUNCTION: compute_sscp_resp
  // TASK: Computes the SSCP matrix X'W worky

  void compute_sscp_resp(datamatrix & H1, datamatrix & workweight, datamatrix & worky);

  // FUNCTION: compute_sscp_resp2
  // TASK: Computes the SSCP matrix (X Z)'W worky

  void compute_sscp_resp2(datamatrix & H1, datamatrix & workweight, datamatrix & worky);

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
                     const ST::string & rname);

  bool check_pause();

  // FUNCTION: out
  // TASK: writes results to outputstream or
  //       in Output window if BORLAND_OUTPUT_WINDOW is defined

  void out(const ST::string & s,bool thick=false,bool italic = false,
           unsigned size = 12,int r=0,int g=0, int b=0);

  void outerror(const ST::string & s);

  };

#endif

//---------------------------------------------------------------------------
#ifndef remlest_multiH
#define remlest_multiH
//---------------------------------------------------------------------------
#endif
 