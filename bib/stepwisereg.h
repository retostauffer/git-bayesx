
// DATE: today


#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif


#if !defined (stepwisereg_INCLUDED)

#define stepwisereg_INCLUDED

#include<statobj.h>
#include<dataobj.h>
#include<map.h>
#include<mapobject.h>

#include<mcmc.h>

#include<distribution.h>
#include<nbinomial.h>

#include<mcmc_const.h>

#include<fullcond_nonp_gaussian.h>
#include<variance_nonp.h>
#include<fullcond_surf_gaussian.h>

#include<fullcond_pspline_gaussian.h>
#include<IWLS_pspline.h>
#include<mcmc_pspline_surf.h>
#include<fullcond_pspline_surf_gaussian.h>
#include<mcmc_pspline.h>

#include<randomeffect.h>

#include<mcmcsimul2.h>

#include<model_stepwise.h>

using MCMC::MCMCoptions;
using MCMC::DISTRIBUTION;
using MCMC::DISTRIBUTION_gaussian;
using MCMC::DISTRIBUTION_binomial;
using MCMC::DISTRIBUTION_poisson;
using MCMC::DISTRIBUTION_gamma;
using MCMC::DISTRIBUTION_nbinomial;
using MCMC::FULLCOND;
using MCMC::FULLCOND_const;
using MCMC::FULLCOND_const_gaussian;
using MCMC::FULLCOND_const_gamma;
using MCMC::FULLCOND_const_nongaussian;
using MCMC::FULLCOND_const_gaussian_special;
using MCMC::FULLCOND_nonp_gaussian;
using MCMC::FULLCOND_variance_nonp;
using MCMC::FULLCOND_pspline;
using MCMC::FULLCOND_pspline_gaussian;
using MCMC::IWLS_pspline;
using MCMC::FULLCOND_pspline_surf;
using MCMC::FULLCOND_pspline_surf_gaussian;
using MCMC::FULLCOND_random_nongaussian;
using MCMC::FULLCOND_random_gaussian;
using MCMC::STEPWISErun;


class __EXPORT_TYPE stepwisereg : public statobject
  {


  private :

  ST::string pathres;
  ST::string title;
  ST::string pathnonp;

  void make_paths(unsigned collinpred,ST::string & pathnonp,
                          ST::string & pathres,ST::string & title,
                          ST::string  varname1,ST::string  varname2,
                          ST::string  endingraw,
                          ST::string  endingres,ST::string  endingtitle) ;

  bool check_gaussian(void);

  bool check_nongaussian(void);

  void clear(void);
  void initpointers(void);

  //------------------------- PRIVATE VARIABLES --------------------------------

  // vector of pointers to current statobjects

  vector<statobject*> * statobj;

  // pointer to functions

  typedef void (* runpointer )(stepwisereg & b);

  runpointer functions[10];

  datamatrix D;
  vector<ST::string> modelvarnamesv;

  // global options

  fileoption outfile;
  optionlist globaloptions;

  // for stepwise

  stroption criterion;

  intoption steps;

  stroption trace;

  intoption number;

  stroption startmodel;

  intoption increment;

  simpleoption fine_tuning;


  intoption maxint;

  vector<ST::string> outfiles;

  vector<FULLCOND*> fullcond;       // Vector of pointers to full conditionals
  STEPWISErun runobj;              

  vector<MCMCoptions> generaloptions;
  bool create_generaloptions(void);

  // OPTIONS for method regress

  ST::string add_name;

  simpleoption constscale;

  // options gamma distributed response
  stroption scalegamma;
  doubleoption scalevalue;
  doubleoption gamvar;
  intoption cit;
  // options gamma distributed response

  // options negative binomial distributed response
  doubleoption propvar;
  stroption distopt;
  stroption propopt;
  simpleoption hierarchical;
  // options negative binomial distributed resposne

  simpleoption predict;                 // indicates that predicted values,
                                        // deviances, etc. should be computed
  simpleoption predictmu;
  intoption predictuntil;

  vector<ST::string> families;          // Response families
  stroption family;                     // specifies the response distribution

  stroption knots;                      // equidistant knots or non equidistant
                                        // knots (P-splines)
  vector<ST::string> knotsdef;


  optionlist regressoptions;


  vector<unsigned> begin_fc;
  vector<unsigned> end_fc;


  friend void regressrun(stepwisereg & b);

  // end for method stepwise

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

  friend void drawmaprun(stepwisereg & b);

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
  stroption title0;
  simpleoption replace2;
  doubleoption xlimtop;
  doubleoption xlimbottom;
  doubleoption xstep;

  optionlist plotnonpoptions;

  use uplotnonp;

  friend void plotnonprun(stepwisereg & b);

//------------------------------ DISTRIBUTION ----------------------------------

  vector<ST::string> distrstring;
  vector<unsigned> distrposition;

  unsigned nrcategories;           // number of categories of the response

  vector<DISTRIBUTION_gaussian> distr_gaussian;
  DISTRIBUTION_binomial distr_binomial;
  DISTRIBUTION_binomial_latent distr_binomlat;
  DISTRIBUTION_poisson distr_poisson;
  DISTRIBUTION_gamma distr_gamma;
  DISTRIBUTION_nbinomial distr_nbinomial;

  vector<DISTRIBUTION *> distr;              // Pointer to distribution objects


  bool create_distribution(void);

  vector<basic_termtype*> termtypes;
  modelterm modreg;
  vector<term> terms;

  use udata;

  bool resultsyesno;

//--------------------------- for the offset -----------------------------------

  term_offset offset;
  bool create_offset(datamatrix & o);

//-------------------------  for fixed effects ---------------------------------

  vector<FULLCOND_const_gaussian> factorgaussian;
  vector<FULLCOND_const_nongaussian> factornongaussian;
  vector<FULLCOND_const_gaussian> normalconst;
  vector<FULLCOND_const_gaussian_special> normalconst_special;
  vector<FULLCOND_const_nongaussian> nongaussianconst;
  vector<FULLCOND_const_gamma> gammaconst;
  FULLCOND_const * fcconst_intercept;

  basic_termtype fixedeffects;

  term_factor_stepwise termfactor;
  term_nonlinearf_stepwise termnonlinearf;

  bool create_const(const unsigned & collinpred=0);
  bool create_factor(const unsigned & collinpred=0);
  bool create_nonlinearf(const unsigned & collinpred=0);

// ----------------------- end: for fixed effects ------------------------------

//------------------------ for nonparametric terms -----------------------------

  // vector<FULLCOND_variance_nonp> fcvarnonp;

  vector<FULLCOND_nonp_gaussian> fcnonpgaussian;
  term_autoreg_stepwise nonprw1rw2;
  term_season_stepwise nonpseason;

  bool create_nonprw1rw2(const unsigned & collinpred=0);
  bool create_nonpseason(const unsigned & collinpred=0);

  term_spatial_stepwise nonpspatial;

  bool create_spatial(const unsigned & collinpred=0);

  vector<FULLCOND_pspline> fcpspline;
  vector<FULLCOND_pspline_gaussian> fcpsplinegaussian;
  vector<IWLS_pspline> fciwlspspline;
  vector<FULLCOND_pspline_surf> fcpsplinesurf;
  vector<FULLCOND_pspline_surf_gaussian> fcpsplinesurfgaussian;
  term_pspline_stepwise nonppspline;
  bool create_pspline(const unsigned & collinpred=0);

//------------------------ for nonparametric terms -----------------------------

//------------------------- for random effects ---------------------------------

  term_random_stepwise randomeff;
  term_randomslope_stepwise randomeffslope;

  vector<FULLCOND_random_nongaussian> fcrandom;
  vector<FULLCOND_random_gaussian> fcrandomgaussian;

  bool create_random(const unsigned & collinpred=0);
  bool create_randomslope(const unsigned & collinpred=0);


  //-------------------- end: for random effects -------------------------------

  void create(void);

  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_pointer * adminp_p;
  #endif

  public:

  //------------------------- PUBLIC FUNCTIONS ---------------------------------

  // DEFAULT CONSTRUCTOR

  stepwisereg (void) : statobject()
         {
         type = "stepwisereg";
         resultsyesno = false;
         }

  // CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - name = n


  #if defined(JAVA_OUTPUT_WINDOW)
  stepwisereg (administrator_basic * adb, administrator_pointer * adp,
               const ST::string & n,ofstream * lo,istream * i,
                                ST::string p,vector<statobject*> * st);
  #else
  stepwisereg (const ST::string & n,ofstream * lo,istream * i,
                                ST::string p,vector<statobject*> * st);
  #endif

  // COPY CONSTRUCTOR

  stepwisereg (const stepwisereg & b);

  // DESTRUCTOR

  ~stepwisereg()
         {
         }


  // OVERLOADED ASSIGNMENT OPERATOR

  const stepwisereg & operator=(const stepwisereg & b);

  int parse(const ST::string & c);

  void describe(optionlist & globaloptions = optionlist());


  };


#endif





