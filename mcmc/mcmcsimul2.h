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

  void modelcomparison(bool & s, const vector<double> & m,
                   const vector<vector<vector<double> > > & mmatrix);

  void maketext(const ST::string & h, const vector<double> & m,
                const double & a, ST::string & text, const bool & neutext,
                const ST::string & tr,const bool & datei);

  void lambdas_update(const vector<double> & m);

  void newmodel_rechnen(bool & f, const vector<double> & m,
     vector<double> & a, vector<vector<double> > & mi,
     const vector<vector<vector<double> > > & mmatrix, const vector<FULLCOND*> & alle,
     const vector<unsigned> & beg, const vector<unsigned> & ende, vector<ST::string> & textit,
     const vector<ST::string> & names_fixed, const vector<vector<ST::string> > & names_nonp,
     const datamatrix & D, const vector<ST::string> & modelv);

  unsigned search_lambdaindex(const double & m, const vector<double> lam,
                                            bool & b) const;
  unsigned search_lambdastartindex(const double & start, const vector<double> & lambdas) const;

  void startwerte(const ST::string & startmodel, const vector<vector<double> > & lambdavec,
                     const vector<ST::string> & names_fixed, vector<vector<unsigned> > & startindex,
                     vector<vector<double> > & startfix);

  void fullcond_entfernen(vector<FULLCOND*> & fullc, const vector<double> & m,
                     const vector<unsigned> & beg, const vector<unsigned> & ende,
                     const vector<ST::string> & names_fixed, const vector<vector<ST::string> > & names_nonp,
                     const datamatrix & D, const vector<ST::string> & modelv);

  void fixed_entfernen(FULLCOND* & fullc, const vector<double> & modell,
                     const vector<ST::string> & names_fixed,
                     const vector<vector<ST::string> > & names_nonp, const datamatrix & D,
                     const vector<ST::string> & modelv);

  void options_text(const int & number, const vector<vector<double> > & lambdavec,
                  const vector<ST::string> & names_fixed,
                  const vector<vector<ST::string> > & names_nonp, const vector<vector<double> > & startfix,
                  const vector<vector<unsigned> > & startindex, const ST::string & name);

  bool finetuning(vector<double> & modell, vector<vector<double> > & lambdavec,
                 vector<ST::string> & names_fixed,
                 const vector<unsigned> & beg, const vector<unsigned> & end,
                 const datamatrix & D, const vector<ST::string> & modelv);

  void initialise_lambdas(vector<vector<ST::string> > & names_nonp,
                 vector<ST::string> & names_fixed, vector<vector<double> > & lambdavec,
                 vector<unsigned> & anfang, vector<unsigned> & ende, const int & number);

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
                       const vector<unsigned> & anfang, const vector<unsigned> & ende,
                       const vector<vector<double> > & lambdavec, const vector<unsigned> & start,
                       const vector<double> & startfix, const vector<FULLCOND*> & fullcond_fest,
                       const vector<ST::string> & names_fixed, const vector<vector<ST::string> > & names_nonp,
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










