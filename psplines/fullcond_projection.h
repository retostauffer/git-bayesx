
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#ifndef fullcond_projectionH
#define fullcond_projectionH

#include<deque>
#include "mcmc_nonpbasis.h"
#include "spline_basis.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_pspline_gaussian ---------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_projection : public spline_basis
  {


  protected:

  unsigned nrvar;
  datamatrix pp_weights;
  datamatrix original_data;
  unsigned nrterms;
  vector<FULLCOND_projection*> pp_pointer;
  datamatrix Bderiv;

  double df_lambdaold;
  double lambdaold;

ofstream outw;

  datamatrix gesamt;

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_projection(void) : spline_basis()
    {
    }

  // CONSTRUCTOR 1  (for additive models)
  // o    : pointer to MCMCoptions object
  // dp   : pointer to DISTRIBUTION object
  // fcc  : pointer to FULLCOND_const object
  // d    : data
  // nrk  : number of knots
  // degr : degree of splines
  // kp   : position of knots (equidistant or quantiles)
  // ft   : field type (RW1, RW2)
  // monotone: increasing || decreasing || unrestricted
  // ti   : title of the object
  // fp   : file where sampled parameters are stored
  // pres : file where results are stored
  // deriv: should the first derivative be computed?
  // l    : starting value for lambda
  // gs   : gridsize
  // diag : should the diagonal transformation be performed?
  // c    : column of the linear predictor (ususally 0)

  FULLCOND_projection(MCMCoptions * o,DISTRIBUTION * dp,
                          FULLCOND_const * fcc, const datamatrix & d,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const fieldtype & ft, const ST::string & monotone, const ST::string & ti,
                         const ST::string & fp, const ST::string & pres, const bool & deriv,
                         const int & gs, vector<FULLCOND_projection*> & zeiger, const unsigned & nterms,
                         const unsigned & c=0);

  // COPY CONSTRUCTOR

  FULLCOND_projection(const FULLCOND_projection & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_projection & operator=(const FULLCOND_projection & fc);


  void compute_linear_combination(bool eins);

  bool posteriormode(void);

  void create_weight(datamatrix & w);

  void reset_effect(const unsigned & pos);

  void update_stepwise(double la);

  double compute_df(void);

  void hierarchie_rw1(vector<double> & untervector, int dfo);

  void compute_lambdavec(vector<double> & lvec, int & number);

  const datamatrix & get_data_forfixedeffects(void);

  ST::string get_effect(void);

  void init_names(const vector<ST::string> & na);

  void outresults(void);

  const datamatrix & get_gesamt(void)
    {
    return gesamt;
    }

  void make_Bspline(const datamatrix & md, const bool & deriv);

  datamatrix bspline_derivative(const double & x);

  void add_linearpred_multBS(const bool & deriv, const bool & current = true);


  // DESTRUCTOR

  ~FULLCOND_projection() {}

  };


} // end: namespace MCMC

#endif

