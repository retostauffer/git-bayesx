// Date: 4.12.99

#if !defined (MCMCpsplinegaussian_INCLUDED)
#define MCMCpsplinegaussian_INCLUDED

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif

#include<deque>
#include<mcmc_nonpbasis.h>
#include<spline_basis.h>

namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_pspline_gaussian ---------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_pspline_gaussian : public spline_basis
  {


  protected:

  bool samplecentered;
  bool diagtransform;
  bool hierarchical;

  double lineff;
  double lineffsum;
  vector<double> lineffsamples;
  datamatrix gamma;

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_pspline_gaussian(void) : spline_basis()
    {
    }

  // CONSTRUCTOR 1  (for additive models)
  // o    : pointer to MCMCoptions object
  // d    : data
  // nrk  : number of knots
  // degr : degree of splines
  // kp   : position of knots (equidistant or quantiles)
  // dp   : pointer to distribution object
  // fp   : file where sampled parameters are stored
  // pres : file where results are stored

  FULLCOND_pspline_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                          FULLCOND_const * fcc, const datamatrix & d,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const fieldtype & ft, const ST::string & monotone, const ST::string & ti,
                         const ST::string & fp, const ST::string & pres, const bool & deriv,
                         const double & l, const int & gs, const bool & diag, const unsigned & c=0);

  // CONSTRUCTOR 2  (for  varying coefficients term)
  // effmod: values of the effect modifier
  // intact: values of the interaction variable

  FULLCOND_pspline_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_const * fcc,
                         const datamatrix & effmod, const datamatrix & intact,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const fieldtype & ft, const ST::string & monotone, const ST::string & ti,
                         const ST::string & fp, const ST::string & pres, const bool & deriv,
                         const double & l, const int & gs, const unsigned & c=0);


  // COPY CONSTRUCTOR

  FULLCOND_pspline_gaussian(const FULLCOND_pspline_gaussian & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_pspline_gaussian & operator=(const FULLCOND_pspline_gaussian & fc);

  void update(void);

  void update_isotonic(void);

  void update_diagtransform(void);

  void outresults(void);

  // FUNCTION: updateK
  // TASK: updates the penalty matrix K

  void updateK(const datamatrix & q)
    {
    FULLCOND_nonp_basis::updateK(q);
    }

  void compute_u(datamatrix & u)
    {
    FULLCOND_nonp_basis::compute_u(u);
    }

  double compute_quadform(void);

  void outoptions(void);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    FULLCOND_nonp_basis::reset();
    }

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  // FUNCTION: predict
  // TASK: predicts the mean for a new observation Xnew

  void predict(const datamatrix & newX, datamatrix & linpred);

  // FUNCTION: write_contour  
  // writes the mean of the full conditional, 1/scale, 1/sigma2 etc. to an ascii file

  void write_contour();

  // FUNCTION: compute_contourprob
  // computes the contour probabilities for differences of order 'diff'

  void compute_contourprob(void);

  void compute_contourprob(const int & diff);

  void compute_pseudocontourprob(const int & diff);

  // DESTRUCTOR

  ~FULLCOND_pspline_gaussian() {}

  };


} // end: namespace MCMC

#endif

