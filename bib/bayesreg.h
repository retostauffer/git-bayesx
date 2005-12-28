
// DATE: today

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif


#if !defined (BAYESREG_INCLUDED)

#define BAYESREG_INCLUDED

#include<statobj.h>
#include<dataobj.h>
#include<map.h>
#include<mapobject.h>

#include<mcmc.h>

#include<distribution.h>
#include<multgaussian.h>
#include<nbinomial.h>
#include<zip.h>
#include<cox.h>
#include<multistate.h>

#include<mcmc_const.h>

#include<mcmc_nonp.h>
#include<fullcond_nonp_gaussian.h>
#include<tvariance.h>
#include<tvariance2dim.h>
#include<fullcond_adaptiv.h>
#include<variance_nonp.h>
#include<fullcond_surf_gaussian.h>

#include<fullcond_pspline_gaussian.h>
#include<IWLS_pspline.h>
#include<mcmc_pspline_surf.h>
#include<fullcond_pspline_surf_gaussian.h>
#include<mcmc_pspline.h>
//#include<isotonic.h>

#include<randomeffect.h>
#include<mixture.h>

#include<kriging2.h>
#include<baseline.h>
#include<multibaseline.h>
//#include<IWLS_baseline.h>

#include<mcmcsimul.h>


using randnumbers::uniform;
using randnumbers::rand_normvek;
using randnumbers::rand_normal;
using randnumbers::rand_invgamma;
using MCMC::MCMCoptions;
using MCMC::DISTRIBUTION;
using MCMC::DISTRIBUTION_gaussian;
using MCMC::DISTRIBUTION_lognormal;
using MCMC::DISTRIBUTION_multgaussian;
using MCMC::DISTRIBUTION_binomial;
using MCMC::DISTRIBUTION_binomial_latent;
using MCMC::DISTRIBUTION_binomial_logit_latent;
using MCMC::DISTRIBUTION_poisson;
using MCMC::DISTRIBUTION_gamma;
using MCMC::DISTRIBUTION_vargaussian;
using MCMC::DISTRIBUTION_nbinomial;
using MCMC::DISTRIBUTION_zip;
using MCMC::DISTRIBUTION_multinom;
using MCMC::DISTRIBUTION_multinomial_latent;
using MCMC::DISTRIBUTION_cumulative_latent3;
using MCMC::DISTRIBUTION_coxmodel;
using MCMC::DISTRIBUTION_multistatemodel;
using MCMC::FULLCOND;
using MCMC::FULLCOND_const;
using MCMC::FULLCOND_const_gaussian;
using MCMC::FULLCOND_const_nongaussian;
using MCMC::FULLCOND_const_nbinomial;
using MCMC::PenaltyMatrix;
using MCMC::FULLCOND_nonp;
using MCMC::FULLCOND_nonp_gaussian;
using MCMC::FULLCOND_tvariance;
using MCMC::FULLCOND_tvariance2dim;
using MCMC::FULLCOND_adaptiv;
using MCMC::FULLCOND_variance_nonp;
using MCMC::FULLCOND_pspline;
using MCMC::FULLCOND_pspline_gaussian;
using MCMC::IWLS_pspline;
//using MCMC::ISOTONIC;
//using MCMC::ISOTONIC_IWLS;
using MCMC::FULLCOND_pspline_surf;
using MCMC::FULLCOND_pspline_surf_gaussian;
using MCMC::FULLCOND_random_nongaussian;
using MCMC::FULLCOND_random_gaussian;
using MCMC::FULLCOND_mixture;
using MCMC::FULLCOND_kriging2;
using MCMC::pspline_baseline;
//using MCMC::IWLS_baseline;
using MCMC::pspline_multibaseline;
using MCMC::MCMCsimulate;


class __EXPORT_TYPE bayesreg : public statobject
  {

  private :

  ST::string pathres;
  ST::string title;
  ST::string pathnonp;

  void make_paths(unsigned collinpred,ST::string & pathnonp,
                          ST::string & pathres,ST::string & title,
                          ST::string  varname1,ST::string  varname2,
                          ST::string  endingraw,ST::string  endingres,
                          ST::string  endingtitle) ;

  bool check_gaussian(void);

  bool check_nongaussian(void);

  bool check_iwls(bool & iwls);

  void clear(void);
  void initpointers(void);

  //------------------------- PRIVATE VARIABLES --------------------------------

  // vector of pointers to current statobjects

  vector<statobject*> * statobj;

  // pointer to functions

  typedef void (* runpointer )(bayesreg & b);

  runpointer functions[10];

  datamatrix D;
  unsigned statepos;
  unsigned begpos;
  vector<ST::string> modelvarnamesv;

  // global options

  fileoption outfile;
  intoption iterationsprint;

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
  simpleoption hclcolors;

  optionlist drawmapoptions;

  use udrawmap;

  friend void __EXPORT_TYPE drawmaprun(bayesreg & b);


  // for method plotnonp

  vector<ST::string> resultfiles;

  modelStandard mplotnonp;

  stroption xlab;
  stroption ylab;
  stroption connect;  
  intoption height;
  intoption width;
  doubleoption ylimtop;
  doubleoption ylimbottom;
  doubleoption ystep;
  doubleoption ystart;
  stroption levels;
  simpleoption median;
  stroption outfile2;
  stroption title0;
  simpleoption replace2;
  doubleoption xlimtop;
  doubleoption xlimbottom;
  doubleoption xstep;
  doubleoption xstart;
  intoption linewidth;
  intoption fontsize;
  intoption pointsize;
  stroption linecolor;
  doubleoption titlescale;

  optionlist plotnonpoptions;

  use uplotnonp;

  friend void __EXPORT_TYPE plotnonprun(bayesreg & b);

  // for method plotautocor

  stroption outfile3;
  intoption maxlag2;
  simpleoption meanonly;
  simpleoption replaceautocor;

  optionlist plotautocoroptions;

  modelStandard mplotautocor;
  use uplotautocor;

  friend void __EXPORT_TYPE plotautocorrun(bayesreg & b);


  // for method 'autocor'

  intoption maxlag;

  optionlist autocorroptions;

  modelStandard ma;

  useDataset ad;

  friend void __EXPORT_TYPE autocorrrun(bayesreg & b);

  // ---------------------  for method 'getsample' -----------------------------

  optionlist getsampleoptions;

  modelStandard mgetsample;

  useDataset usegetsample;

  friend void __EXPORT_TYPE getsamplerun(bayesreg & b);

  // --------------------- for method outresults -------------------------------

  stroption transformtype;
  optionlist outresultsoptions;
  modelStandard moutresults;
  use useoutresults;

  friend void __EXPORT_TYPE outresultsrun(bayesreg & b);

    // for method texsummary

  modelStandard mtexsummary;

  optionlist texsummaryoptions;

  use utexsummary;

  friend void __EXPORT_TYPE texsummaryrun(bayesreg & b);

  // ---------------------  for method 'regress'  ------------------------------

  vector<ST::string> outfiles;

  vector<FULLCOND*> fullcond;       // Vector of pointers to full conditionals
  MCMCsimulate simobj;              // Simulation object;

  vector<MCMCoptions> generaloptions;
  bool create_generaloptions(void);

  // OPTIONS for method regress

  ST::string add_name;

  bool varianceest;
  unsigned varianceend_fc;

  simpleoption modeonly;               // Computes the posterior mode only
  simpleoption noposteriormode;
  intoption setseed;
  simpleoption nographs;

  simpleoption pseudocontourprob;
  simpleoption uniformprior;
  simpleoption approx;
  intoption lengthstart;

  stroption predictind;                // 0/1 indikator variable that
                                       // indicates for which observations
                                       // predictions are necessary

  // general MCMC options
  intoption iterations;                // Number of iterations
  intoption burnin;                    // Number of burnin iterations
  intoption step;                      // Thinning parameter
  doubleoption level1;                 // Nominal level 1 of credible intervals
  doubleoption level2;                 // Nominal level 2 of credible intervals
  simpleoption nosamples;              // Sampled parameters are not stored
                                       // i.e. computation of credible
                                       // interavals etc. not possible
  // general MCMC options

  intoption nutlink;

  simpleoption constscale;

  // options gamma/vargaussian distributed response
  stroption scalegamma;
  doubleoption scalevalue;
  doubleoption gamvar;
  intoption cit;
  // options gamma distributed response

  // options for cumprobit

  simpleoption nosort;

  // end: options for cumprobit

  // options negative binomial distributed response
  doubleoption propvar;
  stroption distopt;
  stroption propopt;
  simpleoption hierarchical;
  // options negative binomial distributed response

  // options zero inflated distributed response
  stroption zipdistopt;
  // options zero inflated distributed response

  // options cox model
  stroption begin;
  stroption state;
  // options multistate model

  simpleoption predict;                 // indicates that predicted values,
                                        // deviances, etc. should be computed
  simpleoption predictmu;
  intoption predictuntil;

  simpleoption constlambda;

  simpleoption noconst;                 // supresses the constant, works only
                                        // for one covariate

  intoption maxint;
  vector<ST::string> families;          // Response families
  stroption family;                     // specifies the response distribution
  doubleoption aresp;                   // Hyperparameter a of overal variance
                                        // (Gaussian response)
  doubleoption bresp;                   // Hyperparameter b of overal variance
                                        // (Gaussian response)

  stroption knots;                      // equidistant knots or non equidistant
                                        // knots (P-splines)
  vector<ST::string> knotsdef;

  intoption nu;

  optionlist regressoptions;

  simpleoption missingreg;
  stroption missingind;
  bool missingest;
  unsigned missingend_fc;
  vector<datamatrix> mind;


  vector<unsigned> begin_fc;
  vector<unsigned> end_fc;


  // end: OPTIONS for method regress

//------------------------------ DISTRIBUTION ----------------------------------

  vector<ST::string> distrstring;
  vector<unsigned> distrposition;

  unsigned nrcategories;           // number of categories of the response

  vector<DISTRIBUTION_gaussian> distr_gaussian;
  DISTRIBUTION_multgaussian distr_multgaussian;
  DISTRIBUTION_lognormal distr_lognormal;  
  DISTRIBUTION_binomial distr_binomial;
  DISTRIBUTION_binomial_latent distr_binomlat;
  DISTRIBUTION_binomial_logit_latent distr_binomlogitlat;
  DISTRIBUTION_poisson distr_poisson;
  DISTRIBUTION_gamma distr_gamma;
  DISTRIBUTION_vargaussian distr_vargaussian;  
  DISTRIBUTION_cumulative_latent3 distr_cumlat3;
  DISTRIBUTION_nbinomial distr_nbinomial;
  DISTRIBUTION_zip distr_zip;  
  DISTRIBUTION_multinom distr_multinom;
  DISTRIBUTION_multinomial_latent distr_multinom_latent;
  DISTRIBUTION_coxmodel distr_cox;
  DISTRIBUTION_multistatemodel distr_multistatemodel;

  doubleoption reference;

  vector<DISTRIBUTION *> distr;              // Pointer to distribution objects


  bool create_distribution(void);

  vector<basic_termtype*> termtypes;
  modelterm modreg;
  vector<term> terms;

  // for multivariate regression

  modeltermmult modregmult;
  vector < vector <term> > termsmult;

  use udata;



  bool resultsyesno;
  bool posteriormode;

//--------------------------- for the offset -----------------------------------

  term_offset offset;
  bool create_offset(datamatrix & o);

//-------------------------  for fixed effects ---------------------------------

  vector<FULLCOND_const_gaussian> normalconst;
  vector<FULLCOND_const_nongaussian> nongaussianconst;
  vector<FULLCOND_const_nbinomial> nbinomialconst;
  FULLCOND_const * fcconst_intercept;


  basic_termtype fixedeffects;

  bool create_const(const unsigned & colllinpred=0);

// ----------------------- end: for fixed effects ------------------------------

//------------------------ for nonparametric terms -----------------------------

  vector<FULLCOND_variance_nonp> fcvarnonp;

  vector<PenaltyMatrix> Pmatrices;
  vector<FULLCOND_nonp> fcnonp;
  vector<FULLCOND_nonp_gaussian> fcnonpgaussian;
  vector<FULLCOND_tvariance> fctvariance;
  vector<FULLCOND_tvariance2dim> fctvariance2dim;
  vector<FULLCOND_adaptiv> fcadaptiv;
  term_autoreg nonprw1rw2;
  term_season nonpseason;

  bool create_nonprw1rw2(const unsigned & collinpred=0);
  bool create_nonpseason(const unsigned & collinpred=0);

  term_spatial nonpspatial;
  term_spatialxy nonpspatialxy;

  bool create_spatial(const unsigned & collinpred=0);
  bool create_spatialxy(const unsigned & collinpred=0);


  vector<FULLCOND_pspline> fcpspline;
  vector<FULLCOND_pspline_gaussian> fcpsplinegaussian;
  vector<IWLS_pspline> fciwlspspline;
  vector<FULLCOND_pspline_surf> fcpsplinesurf;
  vector<FULLCOND_pspline_surf_gaussian> fcpsplinesurfgaussian;
//  vector<ISOTONIC> fcisotonic;
//  vector<ISOTONIC_IWLS> fcisotoniciwls;
  term_varcoeff_pspline nonpvarcoeffpspline;
  term_pspline nonppspline;
  term_interactpspline nonpinteractpspline;
  term_geospline nonpgeospline;
  term_varcoeff_geospline nonpvarcoeffgeospline;
  bool create_pspline(const unsigned & collinpred=0);
  bool create_varcoeffpspline(const unsigned & collinpred=0);
  bool create_interactionspspline(const unsigned & collinpred=0);
  bool create_geospline(const unsigned & collinpred=0);
  bool create_varcoeff_geospline(const unsigned & collinpred=0);

  vector<FULLCOND_kriging2> fckriging;
  term_geokriging nonpspatial_geokriging;
  bool create_geokriging(const unsigned & collinpred=0);

//------------------------ for nonparametric terms -----------------------------

//--------------------------- for baseline terms -------------------------------
  vector<pspline_baseline> fcbaseline;
//  vector<IWLS_baseline> fcbaselineiwls;
  vector<pspline_multibaseline> fcmultibaseline;
  term_baseline baseline;
  term_varcoeff_baseline varcoeffbaseline;
  bool create_baseline(const unsigned & collinpred=0);
  bool create_varcoeffbaseline(const unsigned & collinpred=0);
  bool create_multibaseline(const unsigned & collinpred=0);
  bool create_varcoeffmultibaseline(const unsigned & collinpred=0);
//------------------------ for baseline terms ----------------------------------

  //------------------------- for random effects -------------------------------

  term_random randomeff;
  term_randomslope randomeffslope;
  term_mixture mixtureeff;

  vector<FULLCOND_random_nongaussian> fcrandom;
  vector<FULLCOND_random_gaussian> fcrandomgaussian;
  vector<FULLCOND_mixture> fcmixture;

  bool create_random(const unsigned & collinpred=0);
  bool create_randomslope(const unsigned & collinpred=0);
  bool create_mixture(const unsigned & collinpred=0);


  //-------------------- end: for random effects -------------------------------

  friend void __EXPORT_TYPE regressrun(bayesreg & b);

  friend void __EXPORT_TYPE mregressrun(bayesreg & b);

  //------------------------ end: for method regress ---------------------------

  void create(void);

  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_pointer * adminp_p;
  #endif

  public:


  //------------------------- PUBLIC FUNCTIONS ---------------------------------

  // DEFAULT CONSTRUCTOR

  bayesreg (void) : statobject()
         {
         type = "bayesreg";
         resultsyesno = false;
         }

  // CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - name = n

  #if defined(JAVA_OUTPUT_WINDOW)
  bayesreg (administrator_basic * adb, administrator_pointer * adp,
            const ST::string & n,ofstream * lo,istream * i,
                                ST::string p,vector<statobject*> * st);
  #else
  bayesreg (const ST::string & n,ofstream * lo,istream * i,
                                ST::string p,vector<statobject*> * st);
  #endif

  // COPY CONSTRUCTOR

  bayesreg (const bayesreg & b);

  // DESTRUCTOR

  ~bayesreg()
         {
         }


  // OVERLOADED ASSIGNMENT OPERATOR

  const bayesreg & operator=(const bayesreg & b);


  int parse(const ST::string & c);

  void describe(optionlist & globaloptions = optionlist());


  };

#endif

