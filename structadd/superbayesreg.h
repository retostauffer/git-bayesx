
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (superBAYESREG_INCLUDED)

#define superBAYESREG_INCLUDED

#include"statobj.h"
#include"dataobj.h"
#include"MASTER_obj.h"
#include"GENERAL_OPTIONS.h"

#include"model_parameters.h"

#include"distr.h"
#include"distr_categorical.h"

#include"design.h"
#include"design_pspline.h"
#include"design_hrandom.h"
#include"design_mrf.h"
#include"design_kriging.h"

#include"FC.h"
#include"FC_predict.h"
#include"FC_predictive_check.h"
#include"FC_nonp.h"
#include"FC_linear.h"
#include"FC_hrandom.h"
#include"FC_mult.h"
#include"FC_nonp_variance.h"
#include"FC_hrandom_variance.h"

#include"mcmcsim.h"

using MCMC::MASTER_OBJ;

using randnumbers::uniform;
using randnumbers::rand_normvek;
using randnumbers::rand_normal;
using randnumbers::rand_invgamma;
using MCMC::GENERAL_OPTIONS;

using MCMC::DISTR;
using MCMC::DISTR_gaussian;
using MCMC::DISTR_quantreg;
using MCMC::DISTR_loggaussian;
using MCMC::DISTR_gaussian_re;
using MCMC::DISTR_gaussian_exp;
using MCMC::DISTR_gaussian_mult;
using MCMC::DISTR_binomial;
using MCMC::DISTR_poisson;
using MCMC::DISTR_binomialprobit;

using MCMC::DESIGN_pspline;
using MCMC::DESIGN_hrandom;
using MCMC::DESIGN_mrf;
using MCMC::DESIGN_kriging;
using MCMC::equation;

using MCMC::FC;
using MCMC::FC_predict;
using MCMC::FC_predictive_check;
using MCMC::FC_nonp;
using MCMC::FC_linear;
using MCMC::FC_mult;
using MCMC::FC_hrandom;

using MCMC::FC_nonp_variance;
using MCMC::FC_hrandom_variance;


using MCMC::MCMCsim;


class __EXPORT_TYPE superbayesreg : public statobject
  {

  private :

  ST::string pathres;
  ST::string title;
  ST::string pathnonp;

 void make_header(unsigned & modnr);

  void make_paths(ST::string & pathnonp, ST::string & pathres,
                  ST::string & title, vector<ST::string> vn,
                  ST::string  endingraw,
                  ST::string  endingres, ST::string  endingtitle);

  void extract_data(unsigned i, datamatrix & d,datamatrix & iv,unsigned dim_dm);

  bool findREdistr(ST::string & na,equation & maine,unsigned & fnr);

  bool check_errors(void);

  void clear(void);
//  void initpointers(void);

  //------------------------- PRIVATE VARIABLES --------------------------------

  // vector of pointers to current statobjects

  vector<statobject*> * statobj;

  // pointer to functions

  typedef void (* runpointer )(superbayesreg & b);

  runpointer functions[10];

  datamatrix D;

  vector<ST::string> modelvarnamesv;

  // global options

  fileoption outfile;

  optionlist globaloptions;

  // ---------------------  for method 'regress'  ------------------------------

  vector<equation> equations;          // Vector of equations
  MCMCsim simobj;                      // Simulation object;

  GENERAL_OPTIONS generaloptions;
  bool generaloptions_yes;
  bool create_generaloptions(void);

  // OPTIONS for method regress

  simpleoption modeonly;               // Computes the posterior mode only
  intoption setseed;

  // general MCMC options
  intoption iterations;                // Number of iterations
  intoption burnin;                    // Number of burnin iterations
  intoption step;                      // Thinning parameter
  doubleoption level1;                 // Nominal level 1 of credible intervals
  doubleoption level2;                 // Nominal level 2 of credible intervals

  // general MCMC options

  vector<ST::string> families;          // Response families
  stroption family;                     // specifies the response distribution
  doubleoption aresp;                   // Hyperparameter a of overal variance
                                        // (Gaussian response)
  doubleoption bresp;                   // Hyperparameter b of overal variance
                                        // (Gaussian response)

  intoption equationnr;                 // Equationnumber for multivariate
                                        // responses
  intoption hlevel;                     // hierarchy
  vector<ST::string> equationtypes;
  stroption equationtype;               // type of equation, e.g. mean, variance

  // prediction
  vector<ST::string> predictop;
  stroption predict;

  simpleoption pred_check;

  vector<ST::string> MSEop;
  stroption mse;

  doubleoption quantile;

  // linear effects

  simpleoption centerlinear;

  optionlist regressoptions;

  // end: OPTIONS for method regress

 // ------------------------------- MASTER_OBJ ---------------------------------

 MASTER_OBJ master;

//---------------------------------- DISTR  ------------------------------------

  vector<DISTR_gaussian> distr_gaussians;
  vector<DISTR_quantreg> distr_quantregs;  
  vector<DISTR_loggaussian> distr_loggaussians;
  vector<DISTR_gaussian_re> distr_gaussian_res;
  vector<DISTR_gaussian_exp> distr_gaussian_exps;
  vector<DISTR_gaussian_mult> distr_gaussian_mults;
  vector<DISTR_binomial> distr_binomials;
  vector<DISTR_poisson> distr_poissons;
  vector<DISTR_binomialprobit> distr_binomialprobits;

  bool create_distribution(void);

  bool resultsyesno;
  bool posteriormode;


  //----------------------------------------------------------------------------

  use udata;

  modelterm modreg;
  vector<basic_termtype*> termtypes;

  vector<term> terms;

  term_nonp tnonp;

  //---------------------------- for predict -----------------------------------

  vector<FC_predict> FC_predicts;

  void create_predict(void);

  //----------------------------------------------------------------------------


  //-------------------------- for predictive_checks ---------------------------

  vector<FC_predictive_check> FC_predictive_checks;

  void create_predictive_check(void);

  //----------------------------------------------------------------------------


  //---------------------------- for linear terms ------------------------------

  basic_termtype lineareffects;

  vector<FC_linear> FC_linears;

  bool create_linear(void);

  //----------------------------------------------------------------------------

  //----------------------- for nonparametric terms ----------------------------

  vector<DESIGN_pspline> design_psplines;
  vector<DESIGN_mrf> design_mrfs;
  vector<DESIGN_kriging> design_krigings;
  vector<FC_nonp> FC_nonps;
  vector<FC_nonp_variance> FC_nonp_variances;

  bool create_nonp(void);
  void create_pspline(unsigned i);
  bool create_mrf(unsigned i);
  bool create_kriging(unsigned i);

//------------------------ end for nonparametric terms -------------------------

//----------------------- hierarchical random effects --------------------------

  vector<DESIGN_hrandom>  design_hrandoms;
  vector<FC_hrandom> FC_hrandoms;
  vector<FC_hrandom_variance> FC_hrandom_variances;

  bool create_hrandom(unsigned i);

//------------------- end for hierarchical random effects ----------------------

//---------------------- multiplicative random effects -------------------------

  vector<FC_mult> FC_mults;

  bool create_random_pspline(unsigned i);

//-------------------- end: multiplicative random effects ----------------------

  friend void __EXPORT_TYPE hregressrun(superbayesreg & b);

  bool run_yes;

//------------------------ end: for method regress -----------------------------

//--------------------------- for method autocorr ------------------------------

  intoption maxlag;

  optionlist autocorroptions;

  modelStandard ma;

  useDataset ad;

  friend void __EXPORT_TYPE autocorrrun(superbayesreg & b);

//------------------------ end: for method autocorr ----------------------------

// ------------------------- for method getsample ------------------------------

  optionlist getsampleoptions;

  modelStandard mgetsample;

  useDataset usegetsample;

  friend void __EXPORT_TYPE getsamplerun(superbayesreg & b);

//----------------------- end: for method getsample ----------------------------

  void create(void);
  void create_hregress(void);
  void create_autocorr(void);
  void create_getsample(void);

  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_pointer * adminp_p;
  #endif

  public:


  //------------------------- PUBLIC FUNCTIONS ---------------------------------

  // DEFAULT CONSTRUCTOR

  superbayesreg (void)  : statobject()
    {
    type = "mcmcreg";
    resultsyesno = false;
    }

  // CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - name = n

  #if defined(JAVA_OUTPUT_WINDOW)
  superbayesreg (administrator_basic * adb, administrator_pointer * adp,
                 const ST::string & n,ofstream * lo,istream * i,
                                ST::string p,vector<statobject*> * st);
  #else
  superbayesreg (const ST::string & n,ofstream * lo,istream * i,
                 ST::string p,vector<statobject*> * st);
  #endif

  // COPY CONSTRUCTOR

  superbayesreg (const superbayesreg & b);

  // DESTRUCTOR

  ~superbayesreg()
         {
         }

  // OVERLOADED ASSIGNMENT OPERATOR

  const superbayesreg & operator=(const superbayesreg & b);


  int parse(const ST::string & c);

  void describe(const optionlist & globaloptions = optionlist());

  };

#endif

