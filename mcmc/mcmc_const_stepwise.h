
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (MCMCconststepwise_INCLUDED)

#define MCMCconststepwise_INCLUDED

#include"mcmc_const.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//----------------------- CLASS: FULLCOND_const_stepwise -----------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_const_stepwise : public FULLCOND_const
  {

  protected:

  FULLCOND_const_stepwise * fcconst;
//  vector<FULLCOND*> interactions_pointer;

  vector<ST::string> datanames_fixed_only;
  bool conditional;

  vector<double> diff_categories;
  double reference;
  ST::string coding;

  double const_alt; 
  
  bool changed_data;
  bool changingweight;

  datamatrix X1;
  datamatrix X2;

  datamatrix X1root;     // für MCMC: (X'WX)^-0.5
  datamatrix X1X;        // für MCMC: Gauss: (X'WX)^-1X'W; Nicht-Gauss: (X'WX)

  datamatrix help;

  datamatrix mu1;       // für MCMC: Gauss: y-eta; Nicht-Gauss: X'Wytilde

  ST::string utype;    // gauss oder nongauss

  datamatrix proposal;       // für MCMC - nongauss
  datamatrix weightiwls;
  datamatrix diff;
  datamatrix tildey;
  datamatrix mode;
  datamatrix linmode;

  // FUNCTION: compute_matrices
  // TASK: computes X1 = (X'WX)^-0.5
  //       computes X2 = (X'WX)^-1X'W

  void compute_matrices(void);

  FULLCOND fc_df;


  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_const_stepwise(void) : FULLCOND_const()
    {
    }


  // CONSTRUCTOR_1 linear effects

  FULLCOND_const_stepwise(MCMCoptions * o,DISTRIBUTION * dp,const datamatrix & d,
                          const ST::string & t, const int & constant,
                          const ST::string & fs,const ST::string & fr,
                          const unsigned & c=0);

  // CONSTRUCTOR_4 factor variable

  FULLCOND_const_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                          FULLCOND_const_stepwise * fcc, const datamatrix & d,
                          const ST::string & code, int & ref,
                          const ST::string & t,const ST::string & fs,
                          const ST::string & fr,const unsigned & c=0);


  void make_design(const datamatrix & d);  //, const ST::string & coding);

  void set_pointer_to_interaction(FULLCOND * inter);

  void get_interactionspointer(vector<FULLCOND*> & inter);

  void hierarchical(ST::string & possible);
  
  // COPY CONSTRUCTOR

  FULLCOND_const_stepwise(const FULLCOND_const_stepwise & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_const_stepwise & operator=(const FULLCOND_const_stepwise & m);

  void update_intercept(double & m);

  void update_linold(void); 

  void update_linold_vc(void); 

  void posteriormode_intercept(double & m);

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  void init_name(const ST::string & na);

  void init_names(const vector<ST::string> & na);

  void set_datanames_fixed_only(const vector<ST::string> & na);

  void outresults(void);

  void outoptions(void);

  void compute_lambdavec(vector<double> & lvec, int & number);

  const datamatrix & get_data_forfixedeffects(void)
    {
    assert(fctype==MCMC::factor);

    return data;
    }

  void update_stepwise(double la);

  ST::string get_effect(void);

  void include_effect(const vector<ST::string> & names, const datamatrix & newx);

  void posteriormode_single(const vector<ST::string> & names, datamatrix newx, const bool include);

  void safe_const(void);

  void set_const_old(void);

  void posteriormode_const(void);

  void update_fix_effect(const unsigned & pos, double & value, datamatrix fix);

  void posteriormode_const_varcoeff(datamatrix newx);

  void reset_effect(const unsigned & pos);

  void reset(void);

  void set_effect_zero(void);

  double compute_df(void);

  void update_bootstrap(const bool & uncond=false);

  void update_beta_average(unsigned & samplesize);

  void update_bootstrap_df(void);

  void save_betamean(void);
    
  void update_bootstrap_betamean(void);

  void outresults_df(unsigned & size);

  void update(void);

  void update_gauss(void);

  void update_nongauss(void);

  void compute_XWtildey(datamatrix * linb, double & invscale);

  void set_utype(void)
    {
    utype = "nongauss";
    }

  };


// noch nicht abgeleitet!!!!
//------------------------------------------------------------------------------
//------------------- CLASS: FULLCOND_const_gaussian_special -------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_const_gaussian_special : public FULLCOND_const
  {

  protected:

  datamatrix datatransformed;

  datamatrix mu;

  void compute_datatransformed(double lambda);

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_const_gaussian_special(void) : FULLCOND_const()
    {
    }

  //CONSTRUCTOR1

  FULLCOND_const_gaussian_special(MCMCoptions * o,DISTRIBUTION * dp,
                                  const datamatrix & d,const ST::string & t,
                                  const ST::string & fs,const ST::string & fr,
                                  const unsigned & c);

  // COPY CONSTRUCTOR

  FULLCOND_const_gaussian_special(const FULLCOND_const_gaussian_special & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_const_gaussian_special & operator=(
  const FULLCOND_const_gaussian_special & m);


  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  void outresults(void);

  void outoptions(void)
    {
    FULLCOND_const::outoptions();
    }

  const datamatrix & get_data_forfixedeffects(void)
    {
    return data;
    }

  ST::string  get_effect(void);

  void reset_effect(const unsigned & pos);

  void compute_lambdavec(vector<double> & lvec, int & number);

  double compute_df(void);

  void update_stepwise(double la);

  };


} // end: namespace MCMC

#endif
