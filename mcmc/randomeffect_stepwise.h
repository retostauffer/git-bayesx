#ifndef randomeffect_stepwiseH
#define randomeffect_stepwiseH

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif


#include<randomeffect.h>


namespace MCMC
{


//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_random_stepwise ----------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_random_stepwise : public FULLCOND_random
  {


  protected:

  FULLCOND_nonp_basis * fbasisp;

  double intercept;

  vector< vector<double> > beta_average;
  vector<FULLCOND*> interactions_pointer;
  //int lambda_nr;
  //datamatrix lambdas_local;

  datamatrix data_varcoeff_fix;
  datamatrix effmodi;


  public:

  // DEFAULT CONSTRUCTOR:

  FULLCOND_random_stepwise(void) : FULLCOND_random()
    {
    }

  // CONSTRUCTOR1
  // random intercept

  FULLCOND_random_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                  FULLCOND_const * fcc,
                  const datamatrix & d, const ST::string & t,
                  const ST::string & fp,const ST::string & pr,
                  const double & la, const unsigned & c=0);

  // CONSTRUCTOR2
  // random slope

  FULLCOND_random_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                  FULLCOND_const * fcc,
                  const datamatrix & intvar,const datamatrix & effmod,
                  const ST::string & t,
                  const ST::string & fp,const ST::string & pr,
                  const ST::string & prf,
                  const double & la, const bool & inclfix,
                  const unsigned & c=0);


  // COPY CONSTRUCTOR

  FULLCOND_random_stepwise(const FULLCOND_random_stepwise & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_random_stepwise & operator=(
                        const FULLCOND_random_stepwise & fc);

  // DESTRUCTOR

  ~FULLCOND_random_stepwise() {}



  bool posteriormode(void);

  // FUNCTION: get_effect
  // TASK: returns a string of the estimated effect
  ST::string get_effect(void);

  ST::string get_befehl(void);

  // FUNCTION: reset_effect
  // TASK: resets the effect, subtracts the current effect from linearpred
  void reset_effect(const unsigned & pos);

  void compute_lambdavec(vector<double> & lvec, int & number);

  void update_fix_effect(double & intercept);

  void const_varcoeff(void);

  void hierarchical(ST::string & possible);

  void set_pointer_to_interaction(FULLCOND * inter);

  const datamatrix & get_data_forfixedeffects(void);

  // FUNCTION: compute_df
  // TASK: returns the approximate degrees of freedom of a smoother
  double compute_df(void);

  // FUNCTION: update_stepwise
  // TASK: returns (usually) the current smoothing parameter
  void update_stepwise(double la);

  void save_betas(vector<double> & modell, unsigned & anzahl);

  void average_posteriormode(vector<double> & crit_weights);

  void effect_sort(datamatrix & effect, const double & m, const unsigned & beg,
                   const unsigned & end,const statmatrix<int> & index);

  void effect_sort(datamatrix & effect, const double & m, unsigned & row);

  void init_spatialtotal(FULLCOND_nonp_basis * sp,const ST::string & pnt,
                         const ST::string & prt);

  };     // end: class FULLCOND_random_stepwise

}   // end: namespace MCMC


#endif

