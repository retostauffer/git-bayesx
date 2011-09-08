/* BayesX - Software for Bayesian Inference in
Structured Additive Regression Models.
Copyright (C) 2011  Christiane Belitz, Andreas Brezger,
Thomas Kneib, Stefan Lang, Nikolaus Umlauf

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */



#ifndef fullcond_pspline_surf_gaussianH
#define fullcond_pspline_surf_gaussianH

#include"../export_type.h"
#include "sparsemat.h"
#include "bandmat.h"
#include "bandmat_penalty.h"
#include<deque>
#include "mcmc_nonpbasis.h"
#include "fullcond_pspline_gaussian.h"
#include "spline_basis_surf.h"
#include "fullcond_nonp_gaussian.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_pspline_surf_gaussian ----------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_pspline_surf_gaussian : public spline_basis_surf
  {


  protected:

  double lambda_prec;           // lambda mit dem die aktuelle Pr�zisionsmatrix berechnet wurde

  double f2;                    // Skalierungsfaktor f�r die Parameter a,b bei Gamma-Proposal f�r kappa
                                // bei IWLS update basierend auf dem posteriori Modus
  datamatrix kappaburnin;       // Samples f�r kappa aus der ersten H�lfte der burnin Phase

  updatetype utype;             // iwls || iwlsmode || hyperblock || hyperblockmode
                                // increasing || decreasing || diagtransform

  double kappa;                 // 1/sigma2
  double kappaprop;             // vorgeschlagenes kappa
  double kappacurrent;          // f�r zeilen- und spaltenweises updaten
  double kappamean;             // Hilfvariable f�r Startwert von kappamode
  double kappamode;             // kappa f�r hyperblock
  double kappavar;              // Varianz aus kappaburnin

  bool samplecentered;          // Samplen unter der Nebenbedingung, dass die Zentrierungskonstante = 0 ist.

  unsigned updateW;             // jede wievielte Iteration soll IWLS-Gewicht W neu berechnet werden?

  double a_invgamma;            // Parameter a der IG(a,b) f�r sigma2
  double b_invgamma;            // Parameter a der IG(a,b) f�r sigma2

  datamatrix W;                 // IWLS Gewichtsmatrix
  datamatrix proposal;          // Vorgeschlagenes beta

  bandmatdouble XX;             // X'X (f�r Gauss und zeilenweises updaten)
  bandmatdouble prec;           // Pr�zisionsmatrix (X'WX + 1/sigma2*K) als Band-Matrix (f�r Gauss und zeilenweises updaten)

  envmatdouble XX_env;          // X'X als Envelope-Matrix (f�r IWLS)
  envmatdouble prec_env;        // Pr�zisionsmatrix (X'WX + 1/sigma2*K) als Envelope-Matrix (f�r IWLS)

  datamatrix mu;                // tildey
  datamatrix muy;               // X'W*tildey
  datamatrix standnormal;       // N(0,1)-verteilte ZV

// f�r zeilen- und spaltenweises updaten

  bool singleblock;             // zeilen- und spaltenweises updaten oder komplett?

  vector<datamatrix> Xblock;    // f�r zeilen- und spaltenweises updaten
  vector<envmatdouble> Kblock;  // f�r zeilen- und spaltenweises updaten

  datamatrix beta_ab;           // f�r zeilen- und spaltenweises updaten
  datamatrix prop_ab;           // f�r zeilen- und spaltenweises updaten
  datamatrix beta_mode_ab;      // f�r zeilen- und spaltenweises updaten

  datamatrix muyhelp;           // Hilfsvariable
  datamatrix betahelp;          // Hilfsvariable
  datamatrix betahelp2;         // Hilfsvariable
  bandmatdouble prec2;          // Hilfsmatrix

  void create(const datamatrix & v1, const datamatrix & v2, const datamatrix & intact=datamatrix(1,1));

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
  void compute_XWtildey(const datamatrix & W,const double & scale);

  // f�r zeilen und spalten weises updaten

  void compute_XWX_Block(const datamatrix & W,const unsigned a,const unsigned e,
                         const unsigned beg,const unsigned end,const unsigned & col=0);
  void compute_XWXenv_Block(const datamatrix & W,const unsigned a,const unsigned e,
                            const unsigned beg,const unsigned end,const unsigned & col=0);
  void compute_XWtildey_Block(const datamatrix & W,const double & scale,const unsigned & beg,const unsigned & end);

  // FUNCTION: compute_q
  // TASK: berechnet den Erwartungswertvektor f�r die IWLS proposal bei zeilenweisem updaten

  void compute_q(const datamatrix & beta, const unsigned & an, const unsigned & en,
                 const unsigned & beg, const unsigned & end, const double & sigma2,
                 const bool & current = true);

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

  FULLCOND_pspline_surf_gaussian(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                         const datamatrix & v1, const datamatrix & v2, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const int & gs,
                         const fieldtype & ft, const ST::string & fp, const ST::string & pres,
                         const ST::string & of, const bool & sb, const unsigned & c=0);

  // CONSTRUCTOR 2: IWLS

  // upW  : updateW
  // updatetau : updatetau
  // fstart: Startwert f�r Tuningparameter f
  // a,b  : Hyperparameter f�r IG-prior

  FULLCOND_pspline_surf_gaussian(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                         const datamatrix & v1, const datamatrix & v2,
                         const bool & mode, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const unsigned & upW, const bool & updatetau,
                         const double & fstart, const double & a, const double & b, const int & gs,
                         const fieldtype & ft, const ST::string & fp, const ST::string & pres,
                         const ST::string & of, const bool & iw, const bool & sb, const unsigned & c=0);

  // CONSTRUCTOR 3: geosplines

  // mp   : map object
  // mn   : name of the map object

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

  // intact: Interaktionsvariable

  FULLCOND_pspline_surf_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_const * fcc,const datamatrix &  intact,
                         const datamatrix & v1, const datamatrix & v2, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const int & gs,
                         const fieldtype & ft, const ST::string & fp, const ST::string & pres,
                         const ST::string & of, const bool & sb,
                         const bool & ce, const unsigned & c=0);

  // CONSTRUCTOR 6: IWLS varying coefficients

  FULLCOND_pspline_surf_gaussian(MCMCoptions * o,DISTRIBUTION * dp, FULLCOND_const * fcc,
                         const datamatrix &  intact,
                         const datamatrix & v1, const datamatrix & v2, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const unsigned & upW, const bool & updatetau,
                         const double & fstart, const double & a, const double & b, const int & gs,
                         const fieldtype & ft, const ST::string & fp, const ST::string & pres,
                         const ST::string & of, const bool & iw, const bool & sb,
                         const bool & ce, const unsigned & c=0);

  // CONSTRUCTOR 7: geosplines varying coefficients

  FULLCOND_pspline_surf_gaussian(MCMCoptions * o,DISTRIBUTION * dp, FULLCOND_const * fcc,
                         const datamatrix &  intact,
                         const datamatrix & region,const MAP::map & mp, const ST::string & mn,
                         const ST::string & ti, const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const int & gs, const fieldtype & ft, const ST::string & fp,
                         const ST::string & pres, const bool & sb,const bool & ce, const unsigned & c=0);

  // CONSTRUCTOR 8: IWLS geosplines varying coefficients

  FULLCOND_pspline_surf_gaussian(MCMCoptions * o,DISTRIBUTION * dp, FULLCOND_const * fcc,
                          const datamatrix &  intact,
                         const datamatrix & region,const MAP::map & mp, const ST::string & mn,
                         const ST::string & ti, const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const double & l, const unsigned & upW, const bool & updatetau,
                         const double & fstart, const double & a, const double & b,
                         const int & gs, const fieldtype & ft, const ST::string & fp,
                         const ST::string & pres, const bool & iw, const bool & sb,
                         const bool & ce, const unsigned & c=0);


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

  double compute_df(void);



  // DESTRUCTOR

  ~FULLCOND_pspline_surf_gaussian() {}

  };


} // end: namespace MCMC


//---------------------------------------------------------------------------
#endif
