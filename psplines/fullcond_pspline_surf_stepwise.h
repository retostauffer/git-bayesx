
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#ifndef fullcond_pspline_surf_stepwiseH
#define fullcond_pspline_surf_stepwiseH

#include"FULLCOND_pspline_surf_gaussian.h"
#include"FULLCOND_pspline_stepwise.h"


namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_pspline_surf_stepwise ----------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_pspline_surf_stepwise : public FULLCOND_pspline_surf_gaussian
  {


  protected:

  FULLCOND_pspline_stepwise * mainpoi1;
  FULLCOND_pspline_stepwise * mainpoi2;

  datamatrix splineo1;
  datamatrix splineo2;

  datamatrix data_varcoeff_fix;
  datamatrix effmodi;
  datamatrix XVX;
  double centervalue;

  unsigned maineffectsexisting;      // gibt an, welche Haupteffekte (Fullcond-Obj.) zu Beginn angegeben sind
                                     // Kombinationen: 0 (kein HE), 11 (beide HE)
  double df_lambdaold;
  double lambdaold;
  double lambdaxold;
  double lambdayold;
  double lambdax_prec;
  double lambday_prec;

  envmatdouble KHenv;
  envmatdouble Kxenv;
  envmatdouble Kyenv;

  vector<envmatdouble> all_precenv;      // vector of all possible (X'X + lambda_i P)
  vector<double> lambdavec;

  void create(const datamatrix & v1, const datamatrix & v2, const datamatrix & intact=datamatrix(1,1));

  FULLCOND fc_df;

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_pspline_surf_stepwise(void) : FULLCOND_pspline_surf_gaussian()
    {
    }

  // CONSTRUCTOR

  // o    : pointer to MCMCoptions object
  // dp   : pointer to DISTRIBUTION object
  // fcc  : pointer to FULLCOND_const object
  // v1   : covariate vector 1
  // v2   : covariate vector 2
  // ti   : title
  // nrk  : number of knots
  // degr : degree of splines
  // kp   : position of knots (equidistant or quantiles)
  // l    : lambda
  // gs   : gridsize
  // ft   : fieldtype
  // fp   : file where sampled parameters are stored
  // pres : file where results are stored
  // of   : outfile (wird nicht verwendet)
  // sb   : singleblock

  FULLCOND_pspline_surf_stepwise(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                         const datamatrix & v1, const datamatrix & v2, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const int & gs,
                         const fieldtype & ft, const ST::string & fp, const ST::string & pres,
                         const ST::string & of, const bool & gauss, const unsigned & c=0);

  // CONSTRUCTOR 3: geosplines

  // mp   : map object
  // mn   : name of the map object

  FULLCOND_pspline_surf_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_const * fcc,const datamatrix & region,const MAP::map & mp,
                         const ST::string & mn,
                         const ST::string & ti, const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const int & gs, const fieldtype & ft, const ST::string & fp,
                         const ST::string & pres, const bool & gauss, const unsigned & c=0);


  // CONSTRUCTOR 5: varying coefficients

  // intact: Interaktionsvariable

  FULLCOND_pspline_surf_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_const * fcc,const datamatrix &  intact,
                         const datamatrix & v1, const datamatrix & v2, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const int & gs,
                         const fieldtype & ft, const ST::string & fp, const ST::string & pres,
                         const ST::string & of, const bool & gauss, const bool & vccent, const unsigned & c=0);


  // CONSTRUCTOR 7: geosplines varying coefficients

  FULLCOND_pspline_surf_stepwise(MCMCoptions * o,DISTRIBUTION * dp, FULLCOND_const * fcc,
                         const datamatrix &  intact,
                         const datamatrix & region,const MAP::map & mp, const ST::string & mn,
                         const ST::string & ti, const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const int & gs, const fieldtype & ft, const ST::string & fp,
                         const ST::string & pres, const bool & gauss, const bool & vccent, const unsigned & c=0);


  void init_maineffects(FULLCOND_pspline_stepwise * mp1,FULLCOND_pspline_stepwise * mp2,
                         const ST::string & pnt,const ST::string & prt);

  // COPY CONSTRUCTOR

  FULLCOND_pspline_surf_stepwise(const FULLCOND_pspline_surf_stepwise & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_pspline_surf_stepwise & operator=(const FULLCOND_pspline_surf_stepwise & fc);


  bool posteriormode(void);


// ----------------------- f�r stepwise ----------------------------------------

 /*void hilfeee(void)        // nur f�r Kontrolle!!!
    {
    ST::string test = datanames[0];
    test = test.replaceallsigns('*', '_');
    ofstream outi(("c:\\cprog\\test\\results\\inter_" + test + ".txt").strtochar());
    spline.prettyPrint(outi);
    }*/

  void reset_effect(const unsigned & pos);

  void reset(void);

  void remove_centering_fix(void);

  void hierarchical(ST::string & possible);

  void get_interactionspointer(vector<FULLCOND*> & inter);

  void hierarchie_rw1(vector<double> & untervector, int dfo);

  void compute_lambdavec(vector<double> & lvec, int & number);

  double compute_df(void);  

  void update_stepwise(double la);
    /*{
    lambda=la;
    }*/

  double get_lambda(void)
    {
    return lambda;
    }

  // FUNCTION: get_effect
  // TASK: returns a string of the estimated effect

  ST::string  get_effect(void);

  ST::string get_befehl(void);

  const datamatrix & get_data_forfixedeffects(void);

  void update_fix_effect(void);

  void const_varcoeff(void);

  void update_bootstrap(const bool & uncond=false);

  void update_bootstrap_betamean(void);
  
  void get_samples(const ST::string & filename,const unsigned & step) const;

  void compute_main(void);

  void safe_splines(bool & interact);

  void set_splines_old(void);

  void compute_main_varcoeff(void);

  void multBS_index(datamatrix & res, const datamatrix & b);

  // DESTRUCTOR

  ~FULLCOND_pspline_surf_stepwise() {}

  };


} // end: namespace MCMC



//---------------------------------------------------------------------------
#endif
