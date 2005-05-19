
#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif


#if !defined (MCMCnonpbasis_INCLUDED)

#define MCMCnonpbasis_INCLUDED

#include<mcmc.h>
#include<fullcond.h>
#include<distribution.h>
#include<mcmc_const.h>
#include<bandmat.h>
#include<envmatrix.h>


namespace MCMC
{

enum fieldtype {
                RE,
                RW1,
                RW2,
                seasonal,
                mrf,
                mrfkronecker,
                mrflinear,
                mrflinearband,
                mrfquadratic8,
                mrfquadratic12,
                mrfkr1,
                mrfkr2,
                npspline,
                smoothspline,
                kriging
                };


enum knotpos {equidistant,quantiles};
enum model {additiv,varcoeff};
enum cmode {center_beta,center_spline,center_linpred};
enum rate {alpha_10=10,alpha_30=30,alpha_50=50,alpha_60=60,alpha_70=70,alpha_80=80};

class __EXPORT_TYPE FULLCOND_nonp_basis : public FULLCOND
  {


  protected:

  fieldtype type;
  unsigned period;

  double distance;

  statmatrix<int> index;
  vector<int> posbeg;
  vector<int> posend;

  vector<ST::string> effectvalues;
  vector<double> effectvdouble;

  bandmatdouble K;
  envmatdouble Kenv;
  envmatdouble invprec;
  SparseMatrix Ksp;
  unsigned rankK;

  bool adaptiv;
  datamatrix F1;
  datamatrix F2;
  datamatrix g;

  DISTRIBUTION * likep;

  double sigma2;                     // Varianze parameter /tau^2 in the paper

  bool polex; 

  ST::string pathresults;

  bool varcoeff;

  bool changingweight;

  datamatrix tildey;

  bool interaction;

  double lambda;

  double f;
  unsigned oldacceptance;
  unsigned oldnrtrials;

  FULLCOND fc_contour;         // speichert in jeder Iteration das mu der full
                            // conditional, sowie 1/scale und 1/sigma2
  int contourprob;          // contourprob == -1: keine contour probabilities
                            // contourprob == n : contour probabilities für
                            // n-te Differenzen (und kleiner)

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_nonp_basis(void) : FULLCOND()
    {
    }

  // CONSTRUCTOR
  // o    : pointer to MCMCoptions object
  // dp   : pointer to distribution object
  // minb : Minimum blocksize (minblock)
  // maxb : Maximum blocksize (maxblock)
  // fp   : file where sampled parameters are stored
  // pres : file where results are stored
  // cp   : maximal differences for contour probabilities

  FULLCOND_nonp_basis(MCMCoptions * o,DISTRIBUTION * dp, const fieldtype & ft,
                      const ST::string & ti, const ST::string & fp,
                      const ST::string & pres, const unsigned & c,
                      const unsigned & per=12);

  // CONSTRUCTOR für REML
  
  FULLCOND_nonp_basis(MCMCoptions * o, const ST::string & ti) : FULLCOND(o,ti){}

  // COPY CONSTRUCTOR

  FULLCOND_nonp_basis(const FULLCOND_nonp_basis & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_nonp_basis & operator=(const FULLCOND_nonp_basis & fc);

  void set_interaction(void)
    {
    interaction=true;
    }

  //----------------------------------------------------------------------------
  // ------------- FUNCTIONS FOR UPDATING THE VARIANCE PARAMETER ---------------
  //----------------------------------------------------------------------------

/*
  void update_sigma2samplep(double * p)
    {
    sigma2samplep = p;
    }
*/

  // FUNCTION: update_sigma2
  // TASK: updates sigma2

  void update_sigma2(const double & s)
    {
    sigma2 = s;
    }

  // FUNCTION: compute_quadform
  // TASK: computes beta' K beta

  virtual double compute_quadform(void);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // ------------- FUNCTIONS FOR UPDATING THE PENALTY MATRIX  ------------------
  //----------------------------------------------------------------------------

  // FUNCTION: updateK
  // TASK: updates the penalty matrix K

  virtual void updateK(const datamatrix & q);

  virtual void updateKenv(const datamatrix & q);

  virtual void compute_u(datamatrix & u);

  double compute_ui(unsigned i);

  void set_adaptiv(void);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: write_contour
  // writes the mean of the full conditional, 1/scale and 1/sigma2 to fc_contour

  void write_contour(const datamatrix & m, const double & scaleinv,
                         const double & sigma2inv,
                         const double & bXXb, const double & bKb, const double & mPm, const double & logDetP,
                         const envmatdouble * prec_envp);

  double centerbeta(void)
    {
    return FULLCOND::centerbeta();
    }


  // FUNCTION: get_type
  // TASK: returns prior type (RW1,RW2, etc.)

  const fieldtype & get_type(void)
    {
    return type;
    }

  // FUNCTION: get_rankK
  // TASK: returns the rank of the penalty matrix K

  const unsigned & get_rankK(void)
    {
    return rankK;
    }

  unsigned get_rankK2(void)  // neu!!!
    {                        // neu!!!
    return rankK;            //neu
    }                        //neu

  const double & getlambda(void)
    {
    return lambda;
    }

  const double & get_sigma2(void)
    {
    return sigma2;
    }

  void update(void);

  bool posteriormode(void);

  void outresults(void);

  void outoptions(void);

  void tune_updatetau(const rate & r);

  const vector<ST::string> & get_effectvalues(void)
    {
    return effectvalues;
    }

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void);

  // FUNCTION: predict (virtual)
  // TASK: predicts the mean for a new observation Xnew

  void predict(const datamatrix & newX, datamatrix & linpred);

  // DESTRUCTOR

  ~FULLCOND_nonp_basis() {}

  };


} // end: namespace MCMC

#endif
