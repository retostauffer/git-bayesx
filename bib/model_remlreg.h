

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif


#if !defined (MODELREMLREG_INCLUDED)
#define MODELREMLREG_INCLUDED

#include<model.h>

//------------------------------------------------------------------------------
//-------------------------- class term_autoreg_remlreg ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_autoreg_remlreg : public basic_termtype
  {

  protected:

  doubleoption lambda;
  doubleoption lambdastart;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_autoreg_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a first or second order random walk

  bool checkvector(const vector<term>  & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_autoreg_remlreg() {}

  };

//------------------------------------------------------------------------------
//-------------------- class term_autoreg_varcoef_remlreg ----------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_autoreg_varcoef_remlreg : public basic_termtype
  {

  protected:

  doubleoption lambda;
  doubleoption lambdastart;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_autoreg_varcoef_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a first or second order random walk

  bool checkvector(const vector<term>  & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_autoreg_varcoef_remlreg() {}

  };

//------------------------------------------------------------------------------
//-------------------------- class term_season_remlreg -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_season_remlreg : public basic_termtype
  {

  protected:

  intoption period;
  doubleoption lambda;
  doubleoption lambdastart;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_season_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a seasonal component

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_season_remlreg() {}

  };

//------------------------------------------------------------------------------
//------------------- class term_season_varcoef_remlreg ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_season_varcoef_remlreg : public basic_termtype
  {

  protected:

  intoption period;
  doubleoption lambda;
  doubleoption lambdastart;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_season_varcoef_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a seasonal component

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_season_varcoef_remlreg() {}

  };

//------------------------------------------------------------------------------
//------------------------- class term_pspline_remlreg -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_pspline_remlreg : public basic_termtype
  {

  protected:

  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  intoption gridsize;
  simpleoption diagtransform;
  simpleoption derivative;
  doubleoption lambdastart;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_pspline_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a psline with rw1 oder rw2 Penalty

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_pspline_remlreg() {}

  };

//------------------------------------------------------------------------------
//-------------------------- class term_spatial_remlreg ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_spatial_remlreg : public basic_termtype
  {

  protected:

  stroption map;
  doubleoption lambda;
  doubleoption lambdastart;


  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_spatial_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a spatial component

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_spatial_remlreg() {}

  };

//------------------------------------------------------------------------------
//--------------------- class term_spatial_varcoef_remlreg ---------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_spatial_varcoef_remlreg : public basic_termtype
  {

  protected:

  stroption map;
  doubleoption lambda;
  doubleoption lambdastart;


  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_spatial_varcoef_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a spatial component

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_spatial_varcoef_remlreg() {}

  };

//------------------------------------------------------------------------------
//----------------------- class term_randomslope_remlreg -----------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_randomslope_remlreg : public basic_termtype
  {

  protected:

  doubleoption lambda;
  doubleoption lambdastart;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_randomslope_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a random effect

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_randomslope_remlreg() {}

  };

//------------------------------------------------------------------------------
//---------------------- class term_random_remlreg -----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_random_remlreg : public basic_termtype
  {

  protected:

  doubleoption lambda;
  doubleoption lambdastart;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_random_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a random effect

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_random_remlreg() {}

  };

//------------------------------------------------------------------------------
//-------------------- class term_interactpspline_remlreg ----------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_interactpspline_remlreg : public basic_termtype
  {

  protected:

  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  doubleoption lambdastart;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_interactpspline_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an interaction term

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_interactpspline_remlreg() {}

  };


//------------------------------------------------------------------------------
//--------------------------- class term_geospline_remlreg ---------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_geospline_remlreg : public basic_termtype
  {

  protected:

  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  stroption map;
  doubleoption lambdastart;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_geospline_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an interaction term

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_geospline_remlreg() {}

  };

//------------------------------------------------------------------------------
//------------------- class term_varcoeff_pspline_remlreg ----------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE term_varcoeff_pspline_remlreg : public basic_termtype
  {


  protected:

  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  doubleoption lambdastart;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_varcoeff_pspline_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an interaction term

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ((terms[i].type == "varpsplinerw1") || (terms[i].type == "varpsplinerw2"))
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_varcoeff_pspline_remlreg() {}

  };

//------------------------------------------------------------------------------
//----------------------- class term _kriging_remlreg --------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE term_kriging_remlreg : public basic_termtype
  {


  protected:

  intoption numberknots;
  doubleoption nu;
  doubleoption maxdist;
  simpleoption full;
  stroption knotdata;
  doubleoption p;
  doubleoption q;
  intoption maxsteps;
  doubleoption lambda;
  doubleoption lambdastart;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_kriging_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an interaction term

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if (terms[i].type == "kriging")
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_kriging_remlreg() {}

  };

//------------------------------------------------------------------------------
//----------------------- class term_geokriging_remlreg ------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE term_geokriging_remlreg : public basic_termtype
  {


  protected:

  intoption numberknots;
  doubleoption nu;
  doubleoption maxdist;
  simpleoption full;
  stroption knotdata;
  doubleoption p;
  doubleoption q;
  intoption maxsteps;
  doubleoption lambda;
  doubleoption lambdastart;
  stroption map;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_geokriging_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an interaction term

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if (terms[i].type == "geokriging")
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_geokriging_remlreg() {}

  };

//------------------------------------------------------------------------------
//--------------------- class term_baseline_remlreg ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_baseline_remlreg : public basic_termtype
  {


  protected:

  intoption degree;
  intoption numberknots;
  intoption tgrid;
  stroption gridchoice;
  doubleoption lambda;
  doubleoption lambdastart;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_baseline_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an interaction term

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_baseline_remlreg() {}

  };

//------------------------------------------------------------------------------
//------------------ class term_baseline_varcoef_remlreg -----------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_baseline_varcoeff_remlreg : public basic_termtype
  {

  protected:

  intoption degree;
  intoption numberknots;
  intoption tgrid;
  doubleoption lambda;
  doubleoption lambdastart;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_baseline_varcoeff_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an interaction term

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_baseline_varcoeff_remlreg() {}

  };

#endif
