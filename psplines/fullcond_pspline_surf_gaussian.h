//---------------------------------------------------------------------------
#ifndef fullcond_pspline_surf_gaussianH
#define fullcond_pspline_surf_gaussianH

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif

//#include<mcmc.h>
//#include<fullcond.h>
//#include<distribution.h>
//#include<mcmc_pspline.h>
#include<sparsemat.h>
#include<bandmat.h>
#include<deque>
#include<mcmc_nonpbasis.h>
#include<FULLCOND_pspline_gaussian.h>
#include <spline_basis_surf.h>
#include <fullcond_nonp_gaussian.h>

namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_pspline_surf_gaussian ----------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_pspline_surf_gaussian : public spline_basis_surf
  {


  protected:

  double lambda_prec;

  double f2;
  datamatrix kappaburnin;

  updatetype utype;

  double kappa;
  double kappaprop;
  double kappacurrent;
  double kappamean;
  double kappamode;
  double kappavar;

  bool samplecentered;

  unsigned updateW;

  double a_invgamma;
  double b_invgamma;

  datamatrix beta_ab;
  datamatrix prop_ab;
  datamatrix beta_mode_ab;

  datamatrix W;
  datamatrix proposal;
  bandmatdouble prec2;

  bandmatdouble XX;
  envmatdouble XX_env;
  envmatdouble prec_env;

  vector<datamatrix> Xblock;
  vector<envmatdouble> Kblock;

  bool singleblock;

  bandmatdouble prec;
  datamatrix mu;
  datamatrix muy;
  datamatrix muyhelp;
  datamatrix betahelp;
  datamatrix betahelp2;
  datamatrix betasample;
  datamatrix standnormal;

  void create(const datamatrix & v1, const datamatrix & v2, const datamatrix & intact=datamatrix(1,1));

  void compute_q(const datamatrix & beta, const unsigned & an, const unsigned & en,
                 const unsigned & beg, const unsigned & end, const double & sigma2,
                 const bool & current = true);

  // FUNCTION: add_linearpred_multBS
  // TASK: multipliziert 'X' (Desing-Matrix) mit 'b', speichert das Ergebnis in 'spline'
  //       und addiert das Ergebnis zum linearen Pr�diktor

  void add_linearpred_multBS(const datamatrix & b);

  // FUNCTION: add_linearpred_multBS2
  // TASK: zieht 'spline' vom linearen Pr�diktor ab,
  //       multipliziert 'X' (Desing-Matrix) mit 'b', speichert das Ergebnis in 'spline'
  //       und addiert das Ergebnis zum linearen Pr�diktor

  void add_linearpred_multBS2(const datamatrix & b);

  void compute_XWX(const datamatrix & W,const unsigned & col=0);
  void compute_XWXenv(const datamatrix & W,const unsigned & col=0);

  void compute_XWX_Block(const datamatrix & W,const unsigned a,const unsigned e,
                         const unsigned beg,const unsigned end,const unsigned & col=0);
  void compute_XWXenv_Block(const datamatrix & W,const unsigned a,const unsigned e,
                            const unsigned beg,const unsigned end,const unsigned & col=0);

  void compute_XWtildey(const datamatrix & W,const double & scale);

  void compute_XWtildey_Block(const datamatrix & W,const double & scale,const unsigned & beg,const unsigned & end);

  // f�r posteriormode

  void compute_XWtildey(const datamatrix & W,const datamatrix & tildey,const double & scale,const unsigned & col=0);

  void sample_centered(datamatrix & beta);

  void update_IWLS(void);

  void update_IWLS_mode(void);

  void update_IWLS_hyperblock(void);

  void update_IWLS_hyperblock_mode(void);


  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_pspline_surf_gaussian(void) : spline_basis_surf()
    {
    }

  // CONSTRUCTOR

  // o    : pointer to MCMCoptions object
  // dp   : pointer to DISTRIBUTION object
  // v1   : covariate vector 1
  // v2   : covariate vector 2
  // ti   : title
  // nrk  : number of knots
  // degr : degree of splines
  // kp   : position of knots (equidistant or quantiles)
  // dp   : pointer to distribution object
  // fp   : file where sampled parameters are stored
  // pres : file where results are stored
  // sb   :

  FULLCOND_pspline_surf_gaussian(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                         const datamatrix & v1, const datamatrix & v2, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const int & gs,
                         const fieldtype & ft, const ST::string & fp, const ST::string & pres,
                         const ST::string & of, const bool & sb, const unsigned & c=0);

  // CONSTRUCTOR 2: IWLS

  FULLCOND_pspline_surf_gaussian(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                         const datamatrix & v1, const datamatrix & v2,
                         const bool & mode, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const unsigned & upW, const bool & updatetau,
                         const double & fstart, const double & a, const double & b, const int & gs,
                         const fieldtype & ft, const ST::string & fp, const ST::string & pres,
                         const ST::string & of, const bool & iw, const bool & sb, const unsigned & c=0);

  // CONSTRUCTOR 3: geosplines

  FULLCOND_pspline_surf_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_const * fcc,const datamatrix & region,const MAP::map & mp,
                         const ST::string & mn,
                         const ST::string & ti, const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const int & gs, const fieldtype & ft, const ST::string & fp,
                         const ST::string & pres, const bool & sb, const unsigned & c=0);

  // CONSTRUCTOR 4: IWLS geosplines

  FULLCOND_pspline_surf_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_const * fcc, const datamatrix & region,const MAP::map & mp,
                         const ST::string & mn, const bool & mode, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const unsigned & upW, const bool & updatetau,
                         const double & fstart, const double & a, const double & b,
                         const int & gs, const fieldtype & ft, const ST::string & fp,
                         const ST::string & pres, const bool & iw, const bool & sb, const unsigned & c=0);


  // CONSTRUCTOR 5: varying coefficients

  FULLCOND_pspline_surf_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_const * fcc,const datamatrix &  intact,
                         const datamatrix & v1, const datamatrix & v2, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const int & gs,
                         const fieldtype & ft, const ST::string & fp, const ST::string & pres,
                         const ST::string & of, const bool & sb, const unsigned & c=0);

  // CONSTRUCTOR 6: IWLS varying coefficients

  FULLCOND_pspline_surf_gaussian(MCMCoptions * o,DISTRIBUTION * dp, FULLCOND_const * fcc,
                         const datamatrix &  intact,
                         const datamatrix & v1, const datamatrix & v2, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const unsigned & upW, const bool & updatetau,
                         const double & fstart, const double & a, const double & b, const int & gs,
                         const fieldtype & ft, const ST::string & fp, const ST::string & pres,
                         const ST::string & of, const bool & iw, const bool & sb, const unsigned & c=0);

  // CONSTRUCTOR 7: geosplines varying coefficients

  FULLCOND_pspline_surf_gaussian(MCMCoptions * o,DISTRIBUTION * dp, FULLCOND_const * fcc,
                         const datamatrix &  intact,
                         const datamatrix & region,const MAP::map & mp, const ST::string & mn,
                         const ST::string & ti, const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const int & gs, const fieldtype & ft, const ST::string & fp,
                         const ST::string & pres, const bool & sb, const unsigned & c=0);

  // CONSTRUCTOR 8: IWLS geosplines varying coefficients

  FULLCOND_pspline_surf_gaussian(MCMCoptions * o,DISTRIBUTION * dp, FULLCOND_const * fcc,
                          const datamatrix &  intact,
                         const datamatrix & region,const MAP::map & mp, const ST::string & mn,
                         const ST::string & ti, const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const unsigned & upW, const bool & updatetau,
                         const double & fstart, const double & a, const double & b,
                         const int & gs, const fieldtype & ft, const ST::string & fp,
                         const ST::string & pres, const bool & iw, const bool & sb, const unsigned & c=0);


   void init_maineffects(spline_basis * mp1,spline_basis * mp2,
                         const ST::string & pnt,const ST::string & prt);


  // COPY CONSTRUCTOR

  FULLCOND_pspline_surf_gaussian(const FULLCOND_pspline_surf_gaussian & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_pspline_surf_gaussian & operator=(const FULLCOND_pspline_surf_gaussian & fc);

  void update(void);

  void outresults(void);

  void outoptions(void);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    FULLCOND_nonp_basis::reset();
    }

  double compute_quadform(void)
    {
    if(centertotal)
      return FULLCOND_nonp_basis::compute_quadform();
    else
      return K.compute_quadform(beta_uncentered,0);
    }

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  double compute_squareddiff(unsigned i,unsigned j,unsigned k, unsigned l, unsigned nr);

  // NEU
  void compute_squareddiff(datamatrix & u);

  double * getdiagpointer(void)
    {
    return K.getdiagpointer();
    }

  double * getupperpointer(void)
    {
    return K.getupperpointer();
    }

  void setK(unsigned i, unsigned j, double t)
    {
    K.set(i,j,t);
    }

  // stepwise

  void reset_effect(const unsigned & pos);

  double compute_df(void);

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

  ~FULLCOND_pspline_surf_gaussian() {}

  };


} // end: namespace MCMC


//---------------------------------------------------------------------------
#endif
