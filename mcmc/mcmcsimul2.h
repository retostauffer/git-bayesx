
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (STEPWISErun_INCLUDED)

#define STEPWISErun_INCLUDED

#include"mcmcsimul.h"


namespace MCMC
{

class __EXPORT_TYPE STEPWISErun : public MCMCsimulate
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
  bool fine_tuning;
  bool finelocal;
  ST::string trace;
  double kriterium_tex;
  ofstream outmodels;
  ofstream outcriterium;
  ofstream outtex;
  ST::string smoothing;    // für Unterscheidung globaler / lokaler Glättungsparameter
  bool hierarchical;
  bool miniback_off;
  int bootstrap;
  bool isboot;
  bool unconditional;

  bool df_exact;
  bool ganze_matrix;
  vector<unsigned> xcut;
  vector<unsigned> zcut;
  datamatrix X;
  datamatrix Z;

  //unsigned BIC_min;

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

  bool finetuning(vector<double> & modell);
  bool fine_local(vector<double> & modell);  // für Wahl lokaler Glättungsparameter

  void schaetzen(int z, double & kriterium, bool neu, ST::string variante);

// -----------------------------------------------------------------------------
// -------------- Funktionen, für Stepwise / Stepmin ---------------------------
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
// ------------------ Funktionen für Stepmin -----------------------------------
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
// ------------------ Funktionen für Koordinatenmethode ------------------------
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
// --------------------------- Mini-Backfitting --------------------------------
// -----------------------------------------------------------------------------

  bool blockbilden(vector<FULLCOND*> & fullcond_block, unsigned & z, unsigned & pos);

  void minifullcond_aendern(FULLCOND* & fullcondz, vector<FULLCOND*> & fullcond, unsigned & pos);

  void minibackfitting(vector<FULLCOND*> & fullcond);

// -----------------------------------------------------------------------------
// ------- Funktionen für die Erstellung des Startmodels -----------------------
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
// ------- Funktionen für die Berechnung neuer Modelle bei Stepwise ------------
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

  double df_ganzehatmatrix(void);

// -----------------------------------------------------------------------------
// ------- Funktionen für die Erstellung des fullcondp-Vektors -----------------
// -----------------------------------------------------------------------------

  void fullcond_einzeln(const vector<double> & modell1,
         const vector<double> & modell2, const unsigned & index);

  void fullcond_komplett(const vector<double> & m);

  void fix_komplett(const vector<double> & modell);

  void fix_ganz_komplett(const vector<double> &  modell);

  void reset_fix(const ST::string & name);

  void include_fix(const ST::string & name);

  int column_for_fix(const ST::string & name);

  void korrektur(void);

// -----------------------------------------------------------------------------
// ------- Funktionen für die Ausgabe im Output-Fenster ------------------------
// -----------------------------------------------------------------------------

  bool make_pause(void);

  void maketext(const ST::string & h, const vector<double> & m,
                const double & a, ST::string & text, const bool & neutext,
                const ST::string & tr,const bool & datei);

  void options_text(const int & number, const vector<vector<double> > & startfix,
      const vector<vector<unsigned> > & startindex, const ST::string & name);

// -----------------------------------------------------------------------------
// ------- Funktionen für die Ausgabe im Tex-File ------------------------------
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
// ------- Funktionen für Golden Section Search --------------------------------
// -----------------------------------------------------------------------------

  unsigned golden_section(unsigned & z, double & kriterium);

//  double startbedingungen(unsigned & z, double & kriterium);

//  void approx_zurueck(unsigned & z);

//  void exact_zurueck(unsigned & z);

//  double wert_einzeln(unsigned & z, unsigned i, double & df);

//  double approx_einzeln(unsigned & z, unsigned & i, double & df);

//  double exact_einzeln(unsigned & z, unsigned & i, double & df);

//  int index_suchen(const unsigned & index, const vector<unsigned> & index_vec);

//  int start_b_suchen(vector<unsigned> & index_vec, vector<double> & krit_vec,
//                   double & df, unsigned & b, unsigned & z);

//  double compute_findex(vector<unsigned> & index_vec, vector<double> & krit_vec,
//                  unsigned & index, unsigned & z, double & df);

// -----------------------------------------------------------------------------
// ------------- Model Averaging -----------------------------------------------
// -----------------------------------------------------------------------------

  void update_bootstrap(void);

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

//  void compute_average(ofstream & outweight);

//  void save_alle_betas(vector<double> modell);

/*  void alle_modelle_berechnen(double z, vector<double> & hilf,
                      const vector<double> & best,
                      const vector<double> & unten, const vector<double> & oben,
                      vector<double> & kriterien_alle, vector<double> & priori,
                      vector<vector<double> > & modelle, vector<ST::string> & ausgabe);   */

/*  void occam(vector<double> & kriterien_alle, vector<double> & priori,
           vector<vector<double> > & modelle, vector<ST::string> & ausgabe,
           ofstream & outweight, vector<int> & vorgekommen);   */

/*  void mc3(vector<double> & kriterien_alle, vector<double> & priori,
           vector<vector<double> > & modelle, vector<ST::string> & ausgabe,
           ofstream & outweight, vector<int> & vorgekommen);     */


  public:

  // DEFAULT CONSTRUCTOR

  STEPWISErun(void)
    {
    }

  // CONSTRUCTOR1
  // TASK: initializes the MCMC simulation object with general MCMC options 'go'
  //       distribuiton object dp and a vector of full conditionals 'fc'

  STEPWISErun(MCMCoptions * go,DISTRIBUTION * dp,vector<FULLCOND*> & fc);


  // COPY CONSTRUCTOR

  STEPWISErun(const STEPWISErun & s);

  // OVERLOADED ASSIGNMENT CONSTRUCTOR

  const STEPWISErun & operator=(const STEPWISErun & s);

  bool posteriormode(const vector<ST::string> & header, const bool & presim);

  bool single_stepwise(const vector<unsigned> & start,
                         const vector<double> & startfix, const bool & tex);

  bool stepwise(const ST::string & procedure, const ST::string & minimum,
         const ST::string & crit, const int & stp, const ST::string & trac,
         const int & number, const ST::string & stam, const int & inc, const bool & finet,
         const bool & fineloc, const int & boot, const bool & uncond, 
         const datamatrix & D,const vector<ST::string> & modelv,
         const ST::string & name, vector<FULLCOND*> & fullcond_z, ST::string & path,
         const ST::string & CI, bool & hier, bool & gm, const double & prop, const bool & minib);

  double compute_criterion(void);

  // DESTRUCTOR

  ~STEPWISErun() {}

  };


} // end: namespace MCMC

#endif










