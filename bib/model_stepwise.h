

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif


#if !defined (MODELSTEPWISE_INCLUDED)
#define MODELSTEPWISE_INCLUDED

#include<model.h>

class __EXPORT_TYPE term_nonlinearf_stepwise : public basic_termtype
  {

  protected:

  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption forced_into;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_nonlinearf_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a nonlinearf component

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( terms[i].type == "nonlinearf" )
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_nonlinearf_stepwise() {}

  };


class __EXPORT_TYPE term_autoreg_stepwise : public basic_termtype
  {

  protected:

  doubleoption lambda;
  doubleoption lambdamin;
  doubleoption lambdamax;
  doubleoption lambdastart;
  simpleoption forced_into;
  doubleoption df_for_lambdamax;
  doubleoption df_for_lambdamin;
  simpleoption lambdamax_opt;
  simpleoption lambdamin_opt;
  intoption number;
  simpleoption df_equidist;
  doubleoption df_accuracy;
  simpleoption nofixed;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_autoreg_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a first or second order random walk

  bool checkvector(const vector<term>  & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( (terms[i].type == "rw1") || (terms[i].type == "rw2") ||
         (terms[i].type == "varcoeffrw1") || (terms[i].type == "varcoeffrw2")
       )
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_autoreg_stepwise() {}

  };


class __EXPORT_TYPE term_season_stepwise : public basic_termtype
  {

  protected:

  intoption period;
  doubleoption lambda;
  doubleoption lambdamin;
  doubleoption lambdamax;
  doubleoption lambdastart;
  simpleoption forced_into;
  doubleoption df_for_lambdamax;
  doubleoption df_for_lambdamin;
  simpleoption lambdamax_opt;
  simpleoption lambdamin_opt;
  intoption number;
  simpleoption df_equidist;
  doubleoption df_accuracy;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_season_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a seasonal component

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( ( terms[i].type == "season" ) ||
         ( terms[i].type == "varcoeffseason")
       )
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_season_stepwise() {}

  };


class __EXPORT_TYPE term_pspline_stepwise : public basic_termtype
  {

  protected:

  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  intoption gridsize;
  simpleoption diagtransform;
  simpleoption derivative;
  doubleoption lambdamin;
  doubleoption lambdamax;
  doubleoption lambdastart;
  simpleoption forced_into;
  doubleoption df_for_lambdamax;
  doubleoption df_for_lambdamin;
  simpleoption lambdamax_opt;
  simpleoption lambdamin_opt;
  intoption number;
  simpleoption df_equidist;
  doubleoption df_accuracy;
  stroption monotone;
  simpleoption nofixed;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_pspline_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a psline with rw1 oder rw2 Penalty

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( (terms[i].type == "psplinerw1") || (terms[i].type == "psplinerw2") ||
         (terms[i].type == "varpsplinerw1") || (terms[i].type == "varpsplinerw2") )
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_pspline_stepwise() {}

  };


class __EXPORT_TYPE term_spatial_stepwise : public basic_termtype
  {

  protected:

  stroption map;
  doubleoption lambda;
  doubleoption lambdamin;
  doubleoption lambdamax;
  doubleoption lambdastart;
  simpleoption forced_into;
  doubleoption df_for_lambdamax;
  doubleoption df_for_lambdamin;
  simpleoption lambdamax_opt;
  simpleoption lambdamin_opt;
  intoption number;
  simpleoption df_equidist;
  doubleoption df_accuracy;
  simpleoption nofixed;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_spatial_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a spatial component

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( (terms[i].type == "spatial") ||  (terms[i].type == "varcoeffspatial"))
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_spatial_stepwise() {}

  };


class __EXPORT_TYPE term_randomslope_stepwise : public basic_termtype
  {

  protected:

  simpleoption nofixed;
  doubleoption lambda;
  doubleoption lambdamin;
  doubleoption lambdamax;
  doubleoption lambdastart;
  simpleoption forced_into;
  doubleoption df_for_lambdamax;
  doubleoption df_for_lambdamin;
  simpleoption lambdamax_opt;
  simpleoption lambdamin_opt;
  intoption number;
  simpleoption df_equidist;
  doubleoption df_accuracy;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_randomslope_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a random effect

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( terms[i].type == "randomslope" )
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_randomslope_stepwise() {}

  };


class __EXPORT_TYPE term_random_stepwise : public basic_termtype
  {

  protected:

  doubleoption lambda;
  doubleoption lambdamin;
  doubleoption lambdamax;
  doubleoption lambdastart;
  simpleoption forced_into;
  doubleoption df_for_lambdamax;
  doubleoption df_for_lambdamin;
  simpleoption lambdamax_opt;
  simpleoption lambdamin_opt;
  intoption number;
  simpleoption df_equidist;
  doubleoption df_accuracy;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_random_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a random effect

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( terms[i].type == "random" )
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_random_stepwise() {}

  };


class __EXPORT_TYPE term_factor_stepwise : public basic_termtype
  {

  protected:

  stroption coding;
  doubleoption reference;
  intoption lambdastart;
  simpleoption forced_into;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_factor_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a first or second order random walk

  bool checkvector(const vector<term>  & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if (terms[i].type == "factor")
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_factor_stepwise() {}

  };


//------------------------------------------------------------------------------
//---------------------- class term_interactpspline  ---------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_interactpspline_stepwise : public basic_termtype
  {

  protected:

  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  intoption gridsize;
  doubleoption lambdamin;
  doubleoption lambdamax;
  doubleoption lambdastart;
  simpleoption forced_into;
  doubleoption df_for_lambdamax;
  doubleoption df_for_lambdamin;
  simpleoption lambdamax_opt;
  simpleoption lambdamin_opt;
  intoption number;
  simpleoption df_equidist;
  doubleoption df_accuracy;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_interactpspline_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an interaction term

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_interactpspline_stepwise() {}

  };


//------------------------------------------------------------------------------
//--------------------------- class term_geospline  ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_geospline_stepwise : public basic_termtype
  {

  protected:

  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  stroption map;
  doubleoption lambdamin;
  doubleoption lambdamax;
  doubleoption lambdastart;
  simpleoption forced_into;
  doubleoption df_for_lambdamax;
  doubleoption df_for_lambdamin;
  simpleoption lambdamax_opt;
  simpleoption lambdamin_opt;
  intoption number;
  simpleoption df_equidist;
  doubleoption df_accuracy;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_geospline_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an interaction term

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_geospline_stepwise() {}

  };




#endif
