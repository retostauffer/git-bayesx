// DATE 14.10.99

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif


#if !defined (STEPWISErun_INCLUDED)

#define STEPWISErun_INCLUDED


#include<mcmcsimul.h>


namespace MCMC
{

class __EXPORT_TYPE STEPWISErun : public MCMCsimulate
  {

  protected:

     //  Variablen, die in der ganzen Datei bekannt sind:
  vector<FULLCOND*> fullcond_alle;
  ST::string criterion;
  int increment;
  int steps;
  ST::string startmodel;
  bool fine_tuning;
  ST::string trace;
  double kriterium_tex;
  ofstream outmodels;
  ofstream outcriterium;
  ofstream outtex;

  bool finetuning(vector<double> & modell, vector<vector<double> > & lambdavec,
                 vector<ST::string> & names_fixed,
                 const datamatrix & D, const vector<ST::string> & modelv);

// -----------------------------------------------------------------------------
// ------- Funktionen für die Erstellung des Startmodels -----------------------
// -----------------------------------------------------------------------------

  bool STEPWISErun::vcm_doppelt(const vector<ST::string> & names_fixed,
      const vector<vector<ST::string> > & names_nonp);

  void initialise_lambdas(vector<vector<ST::string> > & names_nonp,
                 vector<ST::string> & names_fixed, vector<vector<double> > & lambdavec,
                 const int & number, const bool & gewichte);

  unsigned search_lambdaindex(const double & m, const vector<double> lam,
                                            bool & b) const;

  unsigned search_lambdastartindex(const double & start,
                           const vector<double> & lambdas) const;

  void startwerte(const ST::string & startmodel,
       const vector<vector<double> > & lambdavec, const vector<ST::string> & names_fixed,
       vector<vector<unsigned> > & startindex, vector<vector<double> > & startfix);

// -----------------------------------------------------------------------------
// ------- Funktionen für die Berechnung neuer Modelle -------------------------
// -----------------------------------------------------------------------------

  void newmodel(bool & fertig, const vector<double> & modell,
    vector<double> & krit, vector<vector<double> > & mi,
    vector<ST::string> & textit);

  void newmodel_fix(bool & fertig, const double & mo, vector<double> & krit,
      vector<vector<double> > & mi, const vector<double> & modell,
      vector<ST::string> & textit, const ST::string & name, const datamatrix & D,
      const vector<ST::string> & modelv);

  void newmodel_factor(bool & fertig, const double & mo,  const unsigned & index,
      vector<double> & krit, vector<vector<double> > & mi, const vector<double> & modell,
      vector<ST::string> & textit, const vector<ST::string> & name,
      const datamatrix & D, const vector<ST::string> & modelv);

  void newmodel_nonp(bool & f, const unsigned & index, const vector<double> & modell,
     const vector<double> & modell_alt, vector<double> & krit, vector<vector<double> > & mi,
     vector<ST::string> & textit, const vector<ST::string> & names_fixed,
     const vector<vector<ST::string> > & names_nonp,
     const datamatrix & D, const vector<ST::string> & modelv);

  bool modelcomparison(const vector<double> & m,
       const vector<vector<vector<double> > > & mmatrix);

// -----------------------------------------------------------------------------
// ------- Funktionen für die Erstellung des fullcondp-Vektors -----------------
// -----------------------------------------------------------------------------

  void fullcond_einzeln(const vector<double> & modell_neu,
         const vector<double> & modell_alt, const unsigned & index,
         const unsigned & nf_size, const vector<vector<ST::string> > & names_nonp,
         const datamatrix & D, const vector<ST::string> & modelv);

  void fullcond_komplett(const vector<double> & m, const unsigned & nf_size,
                     const vector<vector<ST::string> > & names_nonp,
                     const datamatrix & D, const vector<ST::string> & modelv);

  void fix_komplett(const vector<double> & modell,
            const vector<ST::string> & names_fixed,
            const vector<vector<ST::string> > & names_nonp, const datamatrix & D,
            const vector<ST::string> & modelv);

  void reset_fix(const ST::string & name);

  void include_fix(const ST::string & name, const datamatrix & D,
         const vector<ST::string> & modelv);

// -----------------------------------------------------------------------------
// ------- Funktionen für die Ausgabe im Output-Fenster ------------------------
// -----------------------------------------------------------------------------

  bool STEPWISErun::make_pause(void);

  void maketext(const ST::string & h, const vector<double> & m,
                const double & a, ST::string & text, const bool & neutext,
                const ST::string & tr,const bool & datei);

  void options_text(const int & number, const vector<vector<double> > & lambdavec,
      const vector<ST::string> & names_fixed,
      const vector<vector<ST::string> > & names_nonp,
      const vector<vector<double> > & startfix,
      const vector<vector<unsigned> > & startindex, const ST::string & name);

// -----------------------------------------------------------------------------
// ------- Funktionen für die Ausgabe im Tex-File ------------------------------
// -----------------------------------------------------------------------------

  void make_graphics(const ST::string & name, vector<vector<double> > & lambdavec,
                 vector<vector<unsigned> > & startindex);

  void make_tex_end(ST::string & path, const vector<double> & modell,
                    const vector<ST::string> & names_f);

  void make_options(void);

  void make_predictor(void);

  void make_model(void);

  void make_prior(vector<vector<double> > & lambdavec,
                  vector<vector<unsigned> > & startindex);

  void make_fixed_table(void);

  void make_plots(ST::string & path_batch,ST::string & path_splus);
                 // ST::string & path_stata);

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


  bool single_stepwise(double & kriterium_alt, vector<double> & modell_alt, ST::string & text_alt,
                       const vector<vector<double> > & lambdavec, const vector<unsigned> & start,
                       const vector<double> & startfix, const vector<ST::string> & names_fixed,
                       const vector<vector<ST::string> > & names_nonp,
                       const datamatrix & D, const vector<ST::string> & modelv, const bool & tex);

  bool stepwise(const ST::string & crit, const int & stp, const ST::string & trac,
                const int & number, const ST::string & stam, const int & inc,
                const bool & finet, const datamatrix & D,const vector<ST::string> & modelv,
                const ST::string & name, vector<FULLCOND*> & fullcond_z,
                ST::string & path);

  double compute_criterion(void);

  // DESTRUCTOR

  ~STEPWISErun() {}

  };


} // end: namespace MCMC

#endif










