//---------------------------------------------------------------------------
#ifndef fullcond_pspline_stepwiseH
#define fullcond_pspline_stepwiseH

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif

#include<fullcond_pspline_gaussian.h> 


namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_pspline_stepwise ---------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_pspline_stepwise : public FULLCOND_pspline_gaussian
  {


  protected:

  vector< vector<double> > beta_average; 
  vector<FULLCOND*> interactions_pointer;
  int lambda_nr;
  datamatrix lambdas_local;


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
                         const double & l, const int & gs, const unsigned & c=0);


  // COPY CONSTRUCTOR

  FULLCOND_pspline_stepwise(const FULLCOND_pspline_stepwise & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_pspline_stepwise & operator=(const FULLCOND_pspline_stepwise & fc);


  bool posteriormode(void);

  bool changeposterior(const datamatrix & main, const double & inter);   

  bool changeposterior2(const datamatrix & main,const double & inter);    

  /*void hilfeee(void)        // nur für Kontrolle!!!
    {
    ofstream out(("c:\\cprog\\test\\results\\spline_" + datanames[0] + ".txt").strtochar());
    spline.prettyPrint(out);
    } */

  void reset_effect(const unsigned & pos);

  void hierarchie_rw1(vector<double> & untervector);

  void compute_lambdavec(vector<double> & lvec, int & number);

  //double compute_df(void);

  //double compute_df_eigen(void);

  void update_stepwise(double la)         // neu!!!
    {
    if(smoothing == "global")
      lambda=la;
    else
      {
      lambda = 1;
      lambdas_local(lambda_nr+nrpar-rankK,0) = 1/la;
      updateK(lambdas_local);
      }
    }

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

  void set_lambda_nr(void)               // neu!!!
    {
    lambda_nr += 1;
    if(lambda_nr == rankK)
      lambda_nr = 0;
    }

  void multBS_sort(datamatrix & res, const datamatrix & beta);

  // FUNCTION: get_effect
  // TASK: returns a string of the estimated effect

  ST::string get_effect(void);

  const datamatrix & get_data_forfixedeffects(void);

  void save_betas(vector<double> & modell, unsigned & anzahl);

  void average_posteriormode(vector<double> & crit_weights);

  void set_pointer_to_interaction(FULLCOND * inter);

  void search_for_interaction(void);

  void wiederholen(FULLCOND * haupt, bool konst);

  void wiederholen_fix(FULLCOND * haupt, int vorzeichen, bool inter);


  // DESTRUCTOR

  ~FULLCOND_pspline_stepwise() {}

  };


} // end: namespace MCMC


//---------------------------------------------------------------------------
#endif
