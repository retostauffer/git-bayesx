
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (STEPMULTIrun_INCLUDED)

#define STEPMULTIrun_INCLUDED


#include"mcmcsimul.h"


namespace MCMC
{

class __EXPORT_TYPE STEPMULTIrun : public MCMCsimulate
  {

  protected:

// -----------------------------------------------------------------------------
// ------------ Variablen, die in der ganzen Datei bekannt sind: ---------------
// -----------------------------------------------------------------------------

  datamatrix D;
  vector<ST::string> modelv;
  vector<FULLCOND*> fullcond_alle;
  ST::string algorithm;
  ST::string minim;
  ST::string minim2;
  ST::string criterion;
  int increment;
  int steps;
  ST::string startmodel;
  ST::string trace;
  double kriterium_tex;
  ofstream outmodels;
  ofstream outcriterium;
  ofstream outtex;
  bool hierarchical;
  int bootstrap;
  bool isboot;
  bool unconditional;

  unsigned kategorien;
  unsigned katje;
  unsigned anz_fullcond;

  vector<vector<double> > lambdavec;
  vector<ST::string> names_fixed;
  vector<vector<ST::string> > names_nonp;
  vector<double> modell_neu;
  vector<double> modell_alt;
  double kriterium_alt;
  double kriterium_neu;
  ST::string text_alt;
  vector<vector<vector<double> > > modellematrix;
  bool fertig;
  int steps_aktuell;
  //int window;
  vector<ST::string> posttitle;

  void schaetzen(int z, double & kriterium, bool neu, ST::string variante);

// -----------------------------------------------------------------------------
// -------------- Funktionen, f�r Stepwise / Stepmin ---------------------------
// -----------------------------------------------------------------------------

  bool stepfunctions(void);

  unsigned stepwise_fixfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration);

  void stepwise_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z);

  void stepmin_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z);

  void minexact_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z);

// -----------------------------------------------------------------------------
// ------------------ Funktionen f�r Stepmin -----------------------------------
// -----------------------------------------------------------------------------

  void step_minfix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration);

  void stepmin_fix_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration, unsigned & i);

  void stepmin_leer_fix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration, unsigned & i);

  unsigned step_minfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration);

  void stepmin_factor_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z);

  void stepmin_leer_factor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z);

  void stepmin_nonp_nonp(unsigned & z, vector<double> & krit_fkt, double & kriterium);

  void stepmin_nonp_fix(unsigned & z, vector<double> & krit_fkt, double & kriterium);

  void stepmin_nonp_leer(unsigned & z, vector<double> & krit_fkt, double & kriterium);

  void minexact_nonp_nonp(unsigned & z, vector<double> & krit_fkt, double & kriterium);

  void minexact_nonp_fix(unsigned & z, vector<double> & krit_fkt, double & kriterium);

  void minexact_nonp_leer(unsigned & z, vector<double> & krit_fkt, double & kriterium);

  double criterion_min(const double & df);

  double criterion_min(const double & df, const ST::string & auswahl);

// -----------------------------------------------------------------------------
// ------------------ Funktionen f�r Koordinatenmethode ------------------------
// -----------------------------------------------------------------------------

  bool koordabstieg(void);

  void koord_minfix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell);

  void koord_fix_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & i);

  void koord_leer_fix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & i);

  unsigned koord_minfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell);

  void koord_factor_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & z);

  void koord_leer_factor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & z);

  void koord_minnonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z, double & kriterium_aktuell);

  unsigned koordexact_fixfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell);

  void koordexact_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z, double & kriterium_aktuell);

// -----------------------------------------------------------------------------
// ------- Funktionen f�r die Erstellung des Startmodels -----------------------
// -----------------------------------------------------------------------------

  bool vcm_doppelt(void);

  void initialise_lambdas(vector<vector<ST::string> > & namen_nonp,
       vector<ST::string> & namen_fix, vector<vector<double> > & lambdavector,
       const int & number, const bool & gewichte);

  void initialise_weights(double prop);

  unsigned search_lambdaindex(const double & m, const vector<double> lam,
                                            bool & b) const;

  unsigned search_lambdastartindex(const double & start,
                           const vector<double> & lambdas) const;

  void startwerte(const ST::string & startmodel,
       vector<vector<unsigned> > & startindex, vector<vector<double> > & startfix);

// -----------------------------------------------------------------------------
// ------- Funktionen f�r die Berechnung neuer Modelle bei Stepwise ------------
// -----------------------------------------------------------------------------

  void newmodel(vector<double> & krit, vector<vector<double> > & mi,
    vector<ST::string> & textit);

  void newmodel_fix(const double & mo, vector<double> & krit,
      vector<vector<double> > & mi, vector<ST::string> & textit,
      const ST::string & name);

  void newmodel_factor(const double & mo,  const unsigned & index,
      vector<double> & krit, vector<vector<double> > & mi,
      vector<ST::string> & textit, const vector<ST::string> & name);

  void newmodel_nonp(const unsigned & index, vector<double> & krit,
     vector<vector<double> > & mi, vector<ST::string> & textit);

  bool modelcomparison(const vector<double> & m,
       const vector<vector<vector<double> > > & mmatrix);

// -----------------------------------------------------------------------------
// ------- Funktionen f�r die Erstellung des fullcondp-Vektors -----------------
// -----------------------------------------------------------------------------

  void fullcond_einzeln(const vector<double> & modell1,
         const vector<double> & modell2, const unsigned & index);

  void fullcond_komplett(const vector<double> & m);

  void fix_ganz_komplett(const vector<double> &  modell);
  
  void fix_komplett(const vector<double> & modell);

  void reset_fix(const ST::string & name);

  void include_fix(const ST::string & name);

  int column_for_fix(const ST::string & name);

// -----------------------------------------------------------------------------
// ------- Funktionen f�r die Ausgabe im Output-Fenster ------------------------
// -----------------------------------------------------------------------------

  bool make_pause(void);

  void maketext(const ST::string & h, const vector<double> & m,
                const double & a, ST::string & text, const bool & neutext,
                const ST::string & tr,const bool & datei);

  void options_text(const int & number, const vector<vector<double> > & startfix,
      const vector<vector<unsigned> > & startindex, const ST::string & name);

// -----------------------------------------------------------------------------
// ------- Funktionen f�r die Ausgabe im Tex-File ------------------------------
// -----------------------------------------------------------------------------

  void make_graphics(const ST::string & name, vector<vector<unsigned> > & startindex);

  void make_tex_end(ST::string & path, const vector<double> & modell);

  void make_options(void);

  void make_predictor(void);

  void make_model(void);

  void make_prior(vector<vector<unsigned> > & startindex);

  void make_fixed_table(void);

  void make_plots(ST::string & path_batch,ST::string & path_splus);
                 // ST::string & path_stata);


// -----------------------------------------------------------------------------
// ------------- Model Averaging -----------------------------------------------
// -----------------------------------------------------------------------------

  void update_bootstrap(unsigned & zaehler);

  bool confidence_intervals(const ST::string & CI,
          const vector<double> & modell_final,const double & kriterium_final,
          vector<FULLCOND*> & fullcond_z);

  bool confidence_bootstrap(const vector<double> & modell_final,const double & kriterium_final,
                                          vector<FULLCOND*> & fullcond_z);

  bool confidence_MCMCbootstrap(const vector<double> & modell_final,const double & kriterium_final,
                                          vector<FULLCOND*> & fullcond_z);

  bool confidence_MCMCselect(const vector<double> & modell_final,const double & kriterium_final,
                                          vector<FULLCOND*> & fullcond_z);

  bool simulate(const vector<ST::string> & header, const int & seed,
                           const unsigned & startit, const unsigned & endit);



  public:

  // DEFAULT CONSTRUCTOR

  STEPMULTIrun(void)
    {
    }

  // CONSTRUCTOR1
  // TASK: initializes the MCMC simulation object with general MCMC options 'go'
  //       distribuiton object dp and a vector of full conditionals 'fc'

  STEPMULTIrun(MCMCoptions * go,DISTRIBUTION * dp,vector<FULLCOND*> & fc);


  // COPY CONSTRUCTOR

  STEPMULTIrun(const STEPMULTIrun & s);

  // OVERLOADED ASSIGNMENT CONSTRUCTOR

  const STEPMULTIrun & operator=(const STEPMULTIrun & s);

  bool posteriormode(const vector<ST::string> & header, const bool & presim);

  bool single_stepwise(const vector<unsigned> & start,
                         const vector<double> & startfix, const bool & tex);

  bool stepwise(const ST::string & procedure, const ST::string & minimum,
         const ST::string & crit, const int & stp, const ST::string & trac,
         const int & number, const ST::string & stam, const int & inc,
         const int & boot, const bool & uncond, 
         const datamatrix & D,const vector<ST::string> & modelv,
         const ST::string & name, vector<FULLCOND*> & fullcond_z, ST::string & path,
         const ST::string & CI, bool & hier, const double & prop);

  double compute_criterion(void);

  // DESTRUCTOR

  ~STEPMULTIrun() {}

  };


} // end: namespace MCMC

#endif










 