
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#ifndef fullcond_nonp_gaussian_stepwiseH
#define fullcond_nonp_gaussian_stepwiseH

#include"fullcond_nonp_gaussian.h"


namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_nonp_gaussian_stepwise ---------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_nonp_gaussian_stepwise : public FULLCOND_nonp_gaussian
  {


  protected:

  double intercept;
  double df_lambdaold;
  double lambdaold;
  double lambdaold_unstr;
  double df_lambdaold_unstr;

  vector<FULLCOND*> interactions_pointer;
  //int lambda_nr;
  //datamatrix lambdas_local;

  datamatrix data_varcoeff_fix;
  datamatrix effmodi;
  datamatrix XVX;

  FULLCOND * fcunstruct;
  bool spatialtotal;
  //bool gleichwertig;

  vector<envmatdouble> all_precenv;      // vector of all possible (X'X + lambda_i P)
  vector<double> lambdavec;

  FULLCOND fc_df;
  bool isbootstrap;

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_nonp_gaussian_stepwise(void) : FULLCOND_nonp_gaussian()
    {
    }


  // additive Effekte, RW1 RW2 und season

  FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,
                      DISTRIBUTION * dp,
                      const datamatrix & d,
                      FULLCOND_const * fcc,
                      const unsigned & maxint,const fieldtype & ft,
                      const ST::string & ti,
                      const ST::string & fp, const ST::string & pres,
                      const unsigned & c,const double & l,
                      const unsigned & per=12);

  // varying coefficients , RW1 RW2 und season

  FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                       const datamatrix & d,
                       const datamatrix & intvar,
                       FULLCOND_const * fcc,
                       const unsigned & maxint,
                       const fieldtype & ft,const ST::string & ti,
                       const ST::string & fp, const ST::string & pres,
                       const unsigned & c,const double & l, const bool & vccent,
                       const unsigned & per=12);

  // spatial covariates

  FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,
                        DISTRIBUTION * dp,const datamatrix & d,
                        FULLCOND_const * fcc,
                        const MAP::map & m, const ST::string & mn,
                        const ST::string & ti,
                        const ST::string & fp, const ST::string & pres,
                        const unsigned & c,const double & l);

  // varying coefficients , spatial covariates as effect modifier

  FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,
                        DISTRIBUTION * dp,
                        FULLCOND_const * fcc,
                        const MAP::map & m,
                        const ST::string & mn,
                        const datamatrix & d,
                        const datamatrix & d2,
                        const ST::string & ti,
                        const ST::string & fp, const ST::string & pres,
                        const unsigned & c, const double & l, const bool & vccent);

  // COPY CONSTRUCTOR

  FULLCOND_nonp_gaussian_stepwise(const FULLCOND_nonp_gaussian_stepwise & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_nonp_gaussian_stepwise & operator=(const FULLCOND_nonp_gaussian_stepwise & fc);


// --------------------------- FOR STEPWISE ------------------------------------

  double compute_df(void);

  //double compute_df_andererteil(void);

  //void set_gleichwertig(const bool & gleich, bool weiter);  // f�r spatialtotal

  void update_stepwise(double la);
    /*{
    lambda=la;
    } */

  double get_lambda(void)
    {
    return lambda;
    }

  // FUNCTION: get_effect
  // TASK: returns a string of the estimated effect

  ST::string get_effect(void);

  ST::string get_befehl(void);

  void init_names(const vector<ST::string> & na);

  void reset_effect(const unsigned & pos);

  void hierarchie_rw1(vector<double> & untervector, int dfo);

  void compute_lambdavec(vector<double> & lvec, int & number);

  void update_fix_effect(double & intercept);

  void const_varcoeff(void);

  void hierarchical(ST::string & possible);

  void set_pointer_to_interaction(FULLCOND * inter);

  void get_interactionspointer(vector<FULLCOND*> & inter);

  void init_spatialtotal(FULLCOND * unstructp);

  const datamatrix & get_data_forfixedeffects(void);

  void create_weight(datamatrix & w);

  void update_bootstrap(const bool & uncond=false);

  void save_betamean(void);

  void update_bootstrap_betamean(void);
  
  void update(void);

  void update_IWLS(void);

  void update_bootstrap_df(void);

  void outresults_df(unsigned & size);

  void change_Korder(double lam);

  void undo_Korder(void);

  void outresults(void);

  //void save_betas(vector<double> & modell, int & anzahl);

  //void average_posteriormode(vector<double> & crit_weights);

  //void effect_sort(datamatrix & effect, const double & m, const unsigned & beg,
  //                 const unsigned & end,const statmatrix<int> & index);

// Vorschlag:
//  void effect_sort(datamatrix & effect, const double & m, unsigned & row);
  //void effect_sort(datamatrix & effect, double m, unsigned row);

  void createreml(datamatrix & X,datamatrix & Z, const unsigned & Xpos, const unsigned & Zpos);

// ------------------------- END: FOR STEPWISE ---------------------------------

  bool posteriormode(void);

  //void search_for_interaction(void);

  //void hierarchical(ST::string & possible);


  // DESTRUCTOR

  ~FULLCOND_nonp_gaussian_stepwise() {}

  };



} // end: namespace MCMC

#endif
