
// DATE: today

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif


#if !defined (remlreg_INCLUDED)

#define remlreg_INCLUDED

#include<statobj.h>
#include<dataobj.h>
#include<map.h>
#include<mapobject.h>
#include<remlest.h>
#include<remlest_multi.h>
#include<remlest_multi2.h>

#include<mcmc.h>


#include<mcmc_const.h>

#include<fullcond_nonp_gaussian.h>

#include<spline_basis.h>
#include<spline_basis_surf.h>

#include<randomeffect.h>

#include<kriging.h>
#include<baseline_reml.h>

#include<model_remlreg.h>

//#if defined(JAVA_OUTPUT_WINDOW)
//#include<adminparse_basic.h>
//#endif

using MCMC::MCMCoptions;
using MCMC::FULLCOND;
using MCMC::FULLCOND_const;

using MCMC::FULLCOND_nonp_gaussian;

using MCMC::spline_basis;
using MCMC::spline_basis_surf;

using MCMC::FULLCOND_random;

using MCMC::FULLCOND_kriging;
using MCMC::baseline_reml;


class __EXPORT_TYPE remlreg : public statobject
  {


  private :

  void make_paths(unsigned collinpred,ST::string & pathnonp,
                          ST::string & pathres,ST::string & title,
                          ST::string varname1,ST::string varname2,
                          ST::string endingraw,
                          ST::string endingres,ST::string endingtitle) ;

  void clear(void);
  void initpointers(void);

  //------------------------- PRIVATE VARIABLES --------------------------------

  // vector of pointers to current statobjects

  vector<statobject*> * statobj;

  // pointer to functions

  typedef void (* runpointer )(remlreg & b);

  runpointer functions[10];

  datamatrix D;
  vector<ST::string> modelvarnamesv;

  // global options

  fileoption outfile;

  optionlist globaloptions;

  // for method drawmap

  modelStandard mdrawmap;

  simpleoption replace;
  simpleoption swapcolors;
  simpleoption nolegend;
  simpleoption color;
  stroption title2;
  stroption outfile4;
  doubleoption upperlimit;
  doubleoption lowerlimit;
  intoption nrcolors;
  stroption plotvar;
  simpleoption pcat;
  simpleoption drawnames;

  optionlist drawmapoptions;

  use udrawmap;

  friend void drawmaprun(remlreg & b);


  // for method plotnonp

  vector<ST::string> resultfiles;

  modelStandard mplotnonp;

  stroption xlab;
  stroption ylab;
  intoption height;
  intoption width;
  doubleoption ylimtop;
  doubleoption ylimbottom;
  doubleoption ystep;
  stroption levels;
  simpleoption median;
  stroption outfile2;
  stroption title;
  simpleoption replace2;
  doubleoption xlimtop;
  doubleoption xlimbottom;
  doubleoption xstep;


  optionlist plotnonpoptions;

  use uplotnonp;

  friend void plotnonprun(remlreg & b);


  // ---------------------  for method 'regress'  ------------------------------

  remlest RE;
  remlest_multinomial RE_M;
  remlest_ordinal RE_O;
  datamatrix cats;
  bool ismultinomial;
  doubleoption reference;

  friend void remlrun(remlreg & b);

  vector<FULLCOND*> fullcond;       // Vector of pointers to full conditionals

  MCMCoptions generaloptions;

  bool create_data(datamatrix & weight);

  bool create_response(datamatrix & response, datamatrix & weight);

  // OPTIONS for method regress

  ST::string add_name;

  doubleoption level1;                 // Nominal level 1 of credible intervals
  doubleoption level2;                 // Nominal level 2 of credible intervals

  intoption maxint;
  vector<ST::string> families;          // Response families
  stroption family;                     // specifies the response distribution
  stroption knots;                      // equidistant knots or non equidistant
                                        // knots (P-splines)
  intoption maxit;
  doubleoption lowerlim;
  doubleoption eps;                                       

  vector<ST::string> knotsdef;

  optionlist regressoptions;

  // end: OPTIONS for method regress

  vector<basic_termtype*> termtypes;
  modelterm modreg;
  vector<term> terms;

  use udata;



  bool resultsyesno;


//--------------------------- for the offset -----------------------------------

  term_offset offset;
  bool create_offset(datamatrix & o);

//-------------------------  for fixed effects ---------------------------------


  FULLCOND_const * interceptpointer;
  vector<FULLCOND_const> fcconst;

  basic_termtype fixedeffects;

  bool create_const(const unsigned & colllinpred=0);

// ----------------------- end: for fixed effects ------------------------------

//------------------------ for nonparametric terms -----------------------------

  vector<FULLCOND_nonp_gaussian> fcnonpgaussian;

  term_autoreg_remlreg nonprw1rw2;
  term_season_remlreg nonpseason;

  bool create_nonprw1rw2(const unsigned & collinpred=0);
  bool create_nonpseason(const unsigned & collinpred=0);

  term_autoreg_varcoef_remlreg nonprw1rw2_varcoef;
  term_season_varcoef_remlreg nonpseason_varcoef;

  bool create_nonprw1rw2_varcoef(const unsigned & collinpred=0);
  bool create_nonpseason_varcoef(const unsigned & collinpred=0);

  term_spatial_remlreg nonpspatial;
  term_spatialxy nonpspatialxy;

  bool create_spatial(const unsigned & collinpred=0);
  bool create_spatialxy(const unsigned & collinpred=0);

  term_spatial_varcoef_remlreg nonpspatial_varcoef;
  bool create_spatial_varcoef(const unsigned & collinpred=0);

  vector<FULLCOND_kriging> fckriging;
  term_kriging_remlreg nonpspatial_kriging;
  bool create_kriging(const unsigned & collinpred=0);
  term_geokriging_remlreg nonpspatial_geokriging;
  bool create_geokriging(const unsigned & collinpred=0);

  vector<baseline_reml> fcbaseline;
  term_baseline_remlreg nonp_baseline;
  bool create_baseline(const unsigned & collinpred=0);
  vector<baseline_reml> fcbaseline_varcoeff;
  term_baseline_varcoeff_remlreg nonp_baseline_varcoeff;
  bool create_baseline_varcoeff(const unsigned & collinpred=0);

  vector<spline_basis> fcpspline;
  vector<spline_basis_surf> fcpsplinesurf;
  term_varcoeff_pspline_remlreg nonpvarcoeffpspline;
  term_pspline_remlreg nonppspline;
  term_interactpspline_remlreg nonpinteractpspline;
  term_geospline_remlreg nonpgeospline;
  term_varcoeff_geospline nonpvarcoeffgeospline;
  bool create_pspline(const unsigned & collinpred=0);
  bool create_varcoeffpspline(const unsigned & collinpred=0);
  bool create_interactionspspline(const unsigned & collinpred=0);
  bool create_geospline(const unsigned & collinpred=0);
  bool create_varcoeff_geospline(const unsigned & collinpred=0);

//------------------------ for nonparametric terms -----------------------------


  //------------------------- for random effects -------------------------------

  term_random_remlreg randomeff;
  term_randomslope_remlreg randomeffslope;

  vector<FULLCOND_random> fcrandom;

  bool create_random(const unsigned & collinpred=0);
  bool create_randomslope(const unsigned & collinpred=0);

  //-------------------- end: for random effects -------------------------------


  //------------------------ end: for method regress ---------------------------

  void create(void);


  public:


  //------------------------- PUBLIC FUNCTIONS ---------------------------------

  // DEFAULT CONSTRUCTOR

  remlreg (void) : statobject()
         {
         type = "remlreg";
         resultsyesno = false;
         }

  // CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - name = n

  remlreg (
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adb,
  #endif
  const ST::string & n,ofstream * lo,istream * i,
                                ST::string p,vector<statobject*> * st);

  // COPY CONSTRUCTOR

  remlreg (const remlreg & b);

  // DESTRUCTOR

  ~remlreg()
         {
         }


  // OVERLOADED ASSIGNMENT OPERATOR

  const remlreg & operator=(const remlreg & b);


  int parse(const ST::string & c);

  void describe(optionlist & globaloptions = optionlist());

//void plotnonprun(remlreg & b);


  };


#endif





