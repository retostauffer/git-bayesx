//---------------------------------------------------------------------------
#ifndef spline_basisH
#define spline_basisH

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif

#include<deque>
#include<fullcond.h>
#include<mcmc_nonpbasis.h>
#include<bsplinemat.h>

namespace MCMC
{

//---------------------------------------------------------------------------
//----------------------- class: spline_basis -------------------------------
//---------------------------------------------------------------------------


class __EXPORT_TYPE spline_basis : public FULLCOND_nonp_basis
  {

  protected:

  FULLCOND_const * fcconst;

  bool pseudocontourprob;
  bool approx;
  int lengthstart;

  bool lambdaconst;
  bool outbsplines;
  double lambda_prec;

  bool predictright;
  unsigned nrparpredictright;
  bool predictleft;
  unsigned nrparpredictleft;

  bool derivative;
  vector<int> index2;

  bool increasing;
  bool decreasing;

  datamatrix W;
  datamatrix betaold;
  datamatrix betaprop;

  bandmatdouble XX;
  bandmatdouble prec;

  envmatdouble prec_env;
  envmatdouble XX_env;

  datamatrix mu;
  datamatrix muy;
  datamatrix standnormal;
  datamatrix betahelp;


  FULLCOND fchelp;
  FULLCOND fcderivative;

  bsplinemat Bderivative;
  datamatrix splinederivative;

  datamatrix Bcolmean;

  unsigned nrknots;
  unsigned degree;
  unsigned nrdiffobs;

  int gridsize;
  double intercept;
  knotpos knpos;

  vector<int> freq;
  vector<int> freqoutput;
  deque<int> firstnonzero;
  deque<int> lastnonzero;
  deque<double> knot;

  datamatrix xvalues;
  datamatrix spline;
  datamatrix splinehelp;
  datamatrix betaweight;

  datamatrix B;                         // 0 für nonp, X für varcoeff
  datamatrix BS;                        // X für nonp, XZ für varcoeff
  datamatrix G;                         // Transformationsmatrix für diagtransform

  datamatrix X_VCM;                     // für REML VCM
  datamatrix Z_VCM;                     // für REML VCM

  vector<int> begcol;

  datamatrix DG;
  vector<int> DGfirst;


  void make_index(const datamatrix & moddata);

  void make_index(const datamatrix & em,const datamatrix & ia);

  void make_Bspline(const datamatrix & md,const bool & minnull = false);

  void make_BS(const datamatrix & ia);

  void make_DG(void);

// initialisiert xvalues, fchelp und fcderivative

  void init_fchelp(const datamatrix & d);

  void compute_betaweight(void);

  void compute_intercept(void);

  void compute_intercept(const datamatrix & beta);

  void subtr_spline(void);

  // FUNCTION add_linearpred_multBS
  // TASK: multiplies the X-matrix with beta and adds the result to
  //       the current (proposed) linear predictor

  void add_linearpred_multBS(const bool & current = true);
  void add_linearpred_multBS(const datamatrix & beta,const bool & current = true);
  void add_linearpred_multBS(const datamatrix & beta1,const datamatrix & beta2, const bool & current = true);

  // condprior

  void add_linearpred_multBS_Block(const unsigned a, const unsigned e, const datamatrix & b);

  // für Cox-Modell, wenn baseline direkt modelliert wird anstatt log(baseline)  (condprior)

  void add_linearpred_multBS_Block2(const unsigned a, const unsigned e, const datamatrix & b);

  void compute_XWX(const datamatrix & weight);     // nur für diagtransform
  void compute_XWXenv(const datamatrix & weight, const unsigned & c=0);
  void compute_XWtildey(const datamatrix & weight, const double & scale);
  void compute_XWXenv_XWtildey(const datamatrix & weight, const double & scale, const unsigned & c=0);

  // für posteriormode

  void compute_XWtildey(const datamatrix & weight, const datamatrix & tildey, const double & scale, const unsigned & c=0);

  // FUNCTION: sample_centered
  // TASK: Sample under condition x|Ax=0
  //       x - Q^-1AT(AQ^-1AT)^-1(Ax), V=Q^-1AT

  void sample_centered(datamatrix & beta);
  void sample_centered_env(datamatrix & beta);

  void compute_Kweights(void);

  // FUNCTION: bspline_rek
  // TASK: needed in function 'predict' to compute B-Splines for a single observation

  double bspline_rek(unsigned l, unsigned knot, const datamatrix & X);

  void write_spline(void);

  void write_spline(const datamatrix & beta);

  void write_bsplinefunctions(const datamatrix & beta,datamatrix & bsplines);

  void write_derivative(void);

  void make_index2(void);

  void change_K(void);

  void update_prediction(void);


  public:


  // DEFAULT CONSTRUCTOR

  spline_basis(void) : FULLCOND_nonp_basis()
    {
    }

  // CONSTRUCTOR

  spline_basis(MCMCoptions * o, DISTRIBUTION * dp,
               FULLCOND_const * fcc, const fieldtype & ft,
                const ST::string & ti, const unsigned & nrk, const unsigned & degr,
                const MCMC::knotpos & kp, const int & gs, const ST::string & fp,
                const ST::string & pres, const bool & deriv, const unsigned & c);

  // CONSTRUCTOR für REML

  spline_basis(MCMCoptions * o, const datamatrix & d, const unsigned & nrk, const unsigned & degr,
               const knotpos & kp, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl);

  // CONSTRUCTOR für REML VCM

  spline_basis(MCMCoptions * o, const datamatrix & d1, const datamatrix & d2,
               const unsigned & nrk, const unsigned & degr,
               const knotpos & kp, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl);

  // COPY CONSTRUCTOR

  spline_basis(const spline_basis & sp);

  // OVERLOADED ASSIGNMENT OPERATOR

  const spline_basis & operator=(const spline_basis & sp);

  // zum ändern der Haupteffekte bei Interaktionen

  void change(const datamatrix & main, const double & inter);   // Gauss-Fall
  void change(const datamatrix & main);                         // Nicht-Gauss-Fall
  bool changeposterior(const datamatrix & main, const double & inter);   // Gauss-Fall
  bool changeposterior(const datamatrix & main);                         // Nicht-Gauss-Fall

  datamatrix bspline(const double & x);

  void multBS(datamatrix & res, const datamatrix & beta);

  void multBS_index(datamatrix & res, const datamatrix & beta);

  void multDG(datamatrix & res, const datamatrix & b);

  void outoptions(void);

  void outresults(void);

  unsigned & get_nrknots(void)
    {
    return nrknots;
    }

  deque<double> & get_knots(void)
    {
    return knot;
    }

  datamatrix & get_spline(void)
    {
    return spline;
    }

  datamatrix & get_splinehelp(void)
    {
    return splinehelp;
    }

  int * get_indexp(void)
    {
    return index.getV();
    }

  vector<int>::iterator get_freqit(void)
    {
    return freq.begin();
    }

  double get_intercept(void)
    {
    return intercept;
    }

  double * get_fchelpbetap(void)
    {
    return fchelp.getbetapointer();
    }

  void fchelpupdate(void)
    {
    fchelp.update();
    }

  void reset_effect(unsigned & pos);

  void hierarchie_rw1(vector<double> & untervector);  

  void compute_lambdavec(vector<double> & lvec,unsigned & number);  

  bandmatdouble & get_XX(void)
    {
    return XX;
    }

  int & get_gridsize(void)
    {
    return gridsize;
    }

  unsigned & get_degree(void)
    {
    return degree;
    }

  knotpos & get_knotpos(void)
    {
    return knpos;
    }

  void getX(datamatrix & X);

  void init_name(const ST::string & na);

  void init_names(const vector<ST::string> & na);

  void set_lambdaconst(double la);

  void set_contour(int cp, bool pseudocp, bool app, int ls,
                    const datamatrix & b = datamatrix(1,1,0.0));

  void set_outbsplines(void);

  // --------------------------- FOR STEPWISE ----------------------------------

  double compute_df(void);

  void update_stepwise(double la)
    {
    lambda=la;
    }

  // FUNCTION: get_effect
  // TASK: returns a string of the estimated effect

  ST::string  get_effect(void);

  const datamatrix & get_data_forfixedeffects(void);
      

  // REML

  void createreml(datamatrix & X,datamatrix & Z,const unsigned & Xpos,
                  const unsigned & Zpos);

  double outresultsreml(datamatrix & X,datamatrix & Z,
                                     datamatrix & betareml,
                                     datamatrix & betacov,
                                     datamatrix & thetareml,
                                     const unsigned & Xpos,
                                     const unsigned & Zpos,
                                     const unsigned & thetapos,
                                     const bool & dispers,
                                     const unsigned & betaXpos,
                                     const unsigned & betaZpos,
                                     const double & category,
                                     const bool & ismultinomial,
                                     const unsigned plotpos);

  void outoptionsreml();

  // DESTRUCTOR

  ~spline_basis(){}

  };


} // end: namespace MCMC

//---------------------------------------------------------------------------
#endif
