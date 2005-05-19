//---------------------------------------------------------------------------
#ifndef fullcond_pspline_surf_stepwiseH
#define fullcond_pspline_surf_stepwiseH

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif

#include<FULLCOND_pspline_surf_gaussian.h>
#include<FULLCOND_pspline_stepwise.h>


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

  double interhaupt;

  unsigned centerboth;               // gibt an, welche(r) Haupteffekt(e) aktuell im Modell ist
  unsigned maineffectsexisting;      // gibt an, welche Haupteffekte (Fullcond-Obj.) zu Beginn angegeben sind
                                     // Kombinationen: 0 (kein HE), 1 (HE Nr. 1), 10 (HE Nr. 2), 11 (beide HE)
  unsigned centerfix;                // gibt an, welche Haupteffekte linear modelliert werden
                                     // Kombinationen: s. "maineffectsexisting"
  

  void create(const datamatrix & v1, const datamatrix & v2, const datamatrix & intact=datamatrix(1,1));

  void compute_main(void);


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
                         const ST::string & of, const bool & gauss, const unsigned & c=0);


  // CONSTRUCTOR 7: geosplines varying coefficients

  FULLCOND_pspline_surf_stepwise(MCMCoptions * o,DISTRIBUTION * dp, FULLCOND_const * fcc,
                         const datamatrix &  intact,
                         const datamatrix & region,const MAP::map & mp, const ST::string & mn,
                         const ST::string & ti, const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const int & gs, const fieldtype & ft, const ST::string & fp,
                         const ST::string & pres, const bool & gauss, const unsigned & c=0);


  void init_maineffects(FULLCOND_pspline_stepwise * mp1,FULLCOND_pspline_stepwise * mp2,
                         const ST::string & pnt,const ST::string & prt);

  void init_maineffect(FULLCOND_pspline_stepwise * mp1,const ST::string & pnt,
                             const ST::string & prt, const unsigned & number);

  void search_maineffects(void);

  // COPY CONSTRUCTOR

  FULLCOND_pspline_surf_stepwise(const FULLCOND_pspline_surf_stepwise & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_pspline_surf_stepwise & operator=(const FULLCOND_pspline_surf_stepwise & fc);


  bool posteriormode(void);


// ----------------------- für stepwise ----------------------------------------

 /*void hilfeee(void)        // nur für Kontrolle!!!
    {
    ofstream out(("c:\\cprog\\test\\results\\inter_" + datanames[0] + ".txt").strtochar());
    spline.prettyPrint(out);
    }*/

  void reset_effect(const unsigned & pos);

  void remove_centering(void);

  void remove_centering_fix(void);

  void get_zentrierung(FULLCOND * haupt, bool & konst);

  void set_zentrierung(FULLCOND * haupt, int & vorzeichen, bool & inter);    

  void hierarchie_rw1(vector<double> & untervector);

  void compute_lambdavec(vector<double> & lvec, int & number);

  void update_stepwise(double la)
    {
    lambda=la;
    }

  // FUNCTION: get_effect
  // TASK: returns a string of the estimated effect

  ST::string  get_effect(void);

  const datamatrix & get_data_forfixedeffects(void);


  // DESTRUCTOR

  ~FULLCOND_pspline_surf_stepwise() {}

  };


} // end: namespace MCMC



//---------------------------------------------------------------------------
#endif
