
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#ifndef fullcond_pspline_stepwiseH
#define fullcond_pspline_stepwiseH

#include"fullcond_pspline_gaussian.h"
#include "fullcond_nonp_gaussian.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_pspline_stepwise ---------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE FULLCOND_pspline_stepwise : public FULLCOND_pspline_gaussian
  {


  protected:

  vector< vector<double> > beta_average; 
  int lambda_nr;
  datamatrix lambdas_local;

  datamatrix data_varcoeff_fix;
  datamatrix effmodi;
  datamatrix XVX;

  double df_lambdaold;
  double lambdaold;

  vector<envmatdouble> all_precenv;      // vector of all possible (X'X + lambda_i P)
  vector<double> lambdavec;

  envmatdouble Menv;
  bool concave;
  bool convex;
  double lambdamono;

  FULLCOND fc_df;
  bool isbootstrap;
  updatetype utype;   // gaussian || iwls 


  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_pspline_stepwise(void) : FULLCOND_pspline_gaussian()
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

  FULLCOND_pspline_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                          FULLCOND_const * fcc, const datamatrix & d,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const fieldtype & ft, const ST::string & monotone, const ST::string & ti,
                         const ST::string & fp, const ST::string & pres, const bool & deriv,
                         const double & l, const int & gs, const bool & diag, const unsigned & c=0);

  // CONSTRUCTOR 2  (for  varying coefficients term)
  // effmod: values of the effect modifier
  // intact: values of the interaction variable

  FULLCOND_pspline_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_const * fcc,
                         const datamatrix & effmod, const datamatrix & intact,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const fieldtype & ft, const ST::string & monotone, const ST::string & ti,
                         const ST::string & fp, const ST::string & pres, const bool & deriv,
                         const double & l, const int & gs, const bool & vccent,
                         const unsigned & c=0);


  // COPY CONSTRUCTOR

  FULLCOND_pspline_stepwise(const FULLCOND_pspline_stepwise & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_pspline_stepwise & operator=(const FULLCOND_pspline_stepwise & fc);


  bool posteriormode(void);

  bool changeposterior3(const datamatrix & betamain, const datamatrix & main, const double & inter);

  bool changeposterior_varcoeff(const datamatrix & betamain, const datamatrix & main, const double & inter);

  /*void hilfeee(void)        // nur für Kontrolle!!!
    {
    ofstream out(("c:\\cprog\\test\\results\\spline_" + datanames[0] + ".txt").strtochar());
    spline.prettyPrint(out);
    } */

  void reset_effect(const unsigned & pos);

  void reset(void);

  void hierarchie_rw1(vector<double> & untervector, int dfo);

  void compute_lambdavec(vector<double> & lvec, int & number);

  double compute_df(void);

  //double compute_df_eigen(void);

  void update_stepwise(double la);         // neu!!!
/*    {
    if(smoothing == "global")
      {
      lambda=la;

if(likep->iwlsweights_constant() == true)
  {
  bool gefunden = false;
  unsigned i = 0;
  while(i<lambdavec.size() && gefunden == false)
    {
    if(lambda == lambdavec[i])
      gefunden = true;
    i++;
    }
  if(gefunden == true)
    {
    prec_env = all_precenv[i-1];
    lambda_prec = lambda;
    }
  }

      }
    else
      {
      lambda_prec = -1;
      lambdaold = -1;
      lambda = 1;
      lambdas_local(lambda_nr+nrpar-rankK,0) = 1/la;
      updateK(lambdas_local);
      }
    }*/

  void set_lambdas_vector(double & la)       // neu!!!
    {
    lambdas_local = datamatrix(nrpar,1,1/la);
    g = datamatrix(nrpar,1,1);
    if(type==RW2)
      {
      F1 = datamatrix(nrpar,1,-2);
      F2 = datamatrix(nrpar,1,1);
      }
    }

  void set_lambda_nr(void)
    {
    lambda_nr += 1;
    if(lambda_nr >= rankK)  
      lambda_nr = 0;
    }

  double get_lambda(void)
    {
    return lambda;
    }

  /*double compute_penal_lambda(void) //Versuch!!!
    {
    double penal = 0;
    unsigned i;
    double lambda1,lambda2;

    // Bestrafung von Differenzen 1. Ordnung:
    /*for(i=1;i<rankK;i++)
      {
      lambda1 = log10(lambdas_local(i+nrpar-rankK,0));
      lambda2 = log10(lambdas_local(i-1+nrpar-rankK,0));
      penal += (lambda1-lambda2) * (lambda1-lambda2);
      }*/
    /*
    // Bestrafung von Differenzen 2. Ordnung:
    double lambda3;
    for(i=2;i<rankK;i++)
      {
      lambda1 = log10(lambdas_local(i+nrpar-rankK,0));
      lambda2 = log10(lambdas_local(i-1+nrpar-rankK,0));
      lambda3 = log10(lambdas_local(i-2+nrpar-rankK,0));
      penal += (lambda1-2*lambda2+lambda3) * (lambda1-2*lambda2+lambda3);
      }

    return 0.25 * penal;
    }*/

  //void multBS_sort(datamatrix & res, const datamatrix & beta);

  void create_weight(datamatrix & w);

  // FUNCTION: get_effect
  // TASK: returns a string of the estimated effect

  ST::string get_effect(void);

  ST::string get_befehl(void);

  void init_names(const vector<ST::string> & na);

  const datamatrix & get_data_forfixedeffects(void);

  void update_fix_effect(void);

  void const_varcoeff(void);

  //void save_betas(vector<double> & modell, int & anzahl);

  //void average_posteriormode(vector<double> & crit_weights);

  void set_pointer_to_interaction(FULLCOND * inter);

  void get_interactionspointer(vector<FULLCOND*> & inter);

  bool search_for_interaction(void);

  void hierarchical(ST::string & possible);

  void createreml(datamatrix & X,datamatrix & Z,
                           const unsigned & Xpos, const unsigned & Zpos);

  void updateMenv(void);

  void set_spmonotone(double & spmono)
    {
    lambdamono = spmono;
    }

  void update_bootstrap(const bool & uncond=false);

  void update_beta_average(unsigned & samplesize);

  void save_betamean(void);

  void update_bootstrap_betamean(void);

  void update(void);

  void update_bootstrap_df(void);

  void outresults_df(unsigned & size);

  void change_Korder(double lam);

  void undo_Korder(void);

  void get_samples(const ST::string & filename,const unsigned & step) const;

  void change_varcoeff(const datamatrix & betamain,const datamatrix & main,const double & inter);

  void update_gauss(void);

  void update_IWLS(void);

  void set_utype(void)
    {
    utype = iwls;
    }

  void outresults(void);

  vector<int>::iterator get_freqoutputit(void)
    {
    return freqoutput.begin();
    }

  void set_spline(datamatrix & sp)
    {
    spline.assign(sp);
    }

  // DESTRUCTOR

  ~FULLCOND_pspline_stepwise() {}

  };


} // end: namespace MCMC


//---------------------------------------------------------------------------
#endif
