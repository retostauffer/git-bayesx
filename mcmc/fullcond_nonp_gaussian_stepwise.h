//---------------------------------------------------------------------------
#ifndef fullcond_nonp_gaussian_stepwiseH
#define fullcond_nonp_gaussian_stepwiseH

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif

#include<fullcond_nonp_gaussian.h>


namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_nonp_gaussian_stepwise ---------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_nonp_gaussian_stepwise : public FULLCOND_nonp_gaussian
  {


  protected:

  double intercept;

  vector< vector<double> > beta_average;
  vector<FULLCOND*> interactions_pointer;
  //int lambda_nr;
  //datamatrix lambdas_local;

  datamatrix data_varcoeff_fix;
  datamatrix effmodi;


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
                       const unsigned & c,const double & l, const bool & nofixed,
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
                        const unsigned & c, const double & l, const bool & nofixed);

  // COPY CONSTRUCTOR

  FULLCOND_nonp_gaussian_stepwise(const FULLCOND_nonp_gaussian_stepwise & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_nonp_gaussian_stepwise & operator=(const FULLCOND_nonp_gaussian_stepwise & fc);


// --------------------------- FOR STEPWISE ------------------------------------

  double compute_df(void);

  void update_stepwise(double la)
    {
    lambda=la;
    }

  // FUNCTION: get_effect
  // TASK: returns a string of the estimated effect

  ST::string get_effect(void);

  ST::string get_befehl(void);

  void reset_effect(const unsigned & pos);

  void hierarchie_rw1(vector<double> & untervector);

  void compute_lambdavec(vector<double> & lvec, int & number);

  void update_fix_effect(double & intercept);

  void const_varcoeff(void);

  void hierarchical(ST::string & possible);

  void set_pointer_to_interaction(FULLCOND * inter);

  const datamatrix & get_data_forfixedeffects(void);

  void save_betas(vector<double> & modell, unsigned & anzahl);

  void average_posteriormode(vector<double> & crit_weights);

  void effect_sort(datamatrix & effect, const double & m, const unsigned & beg,
                   const unsigned & end,const statmatrix<int> & index);

  void effect_sort(datamatrix & effect, const double & m, unsigned & row);

// ------------------------- END: FOR STEPWISE ---------------------------------

  bool posteriormode(void);

  //void search_for_interaction(void);

  //void hierarchical(ST::string & possible);


  // DESTRUCTOR

  ~FULLCOND_nonp_gaussian_stepwise() {}

  };



} // end: namespace MCMC

#endif
