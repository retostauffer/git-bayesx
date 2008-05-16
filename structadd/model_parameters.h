#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (MODELparameters_INCLUDED)
#define MODELparameters_INCLUDED

#include"clstring.h"
#include<vector>
#include"data.h"
#include"option.h"
#include"model.h"



//------------------------------------------------------------------------------
//--------------------------- class term_nonp ----------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE term_nonp : public basic_termtype
  {
  protected:

  intoption degree;
  intoption numberknots;
  intoption difforder;
  doubleoption lambda;
  doubleoption a;
  doubleoption b;
  simpleoption nocenter;
  stroption map;
  doubleoption lambda_re;
  doubleoption a_re;
  doubleoption b_re;
  simpleoption internal_mult;


  vector<ST::string> termnames;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_nonp(void) {}

  // CONSTRUCTOR

  term_nonp(vector<ST::string> & na);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a nonparametric term

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_nonp() {}

  };

#endif
