// DATE: 24.10.1999


#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif


#if !defined (MODEL_INCLUDED)
#define MODEL_INCLUDED

#include<clstring.h>
#include<vector>
#include<data.h>
#include<option.h>



//------------------------------------------------------------------------------
//--------------------------- CLASS: model -------------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE model

  {

  protected:


  // ---------------------- PROTECTED VARIABLES --------------------------------

  // curent model

  ST::string modeltext;

  // true, if a correct model has been specified

  bool modelexisting;

  // contains current errormessages, if modelspecification is not correct

  vector<ST::string> errormessages;
//  errorm::messages errormessages;

  // Variablenames of the corresponding model designmatrix
  // (the designmatrix can be obtained with method 'makeModelMatrix'

  list<ST::string> modelVarnames;


  // --------------------  PROTECTED FUNCTIONS ---------------------------------

  // VIRTUAL FUNCTION: clear
  // TASK: clears 'modelVarnames',
  //       modelexisting = false

  virtual void clear(void)
	 {
     modeltext = "";
	 modelexisting = false;
	 if (! modelVarnames.empty())
		modelVarnames.erase(modelVarnames.begin(),modelVarnames.end());
	 }


  public:


  // ------------------------- PUBLIC FUNCTIONS --------------------------------

  // DEFAULT CONSTRUCTOR

  model(void)
	 {
	 modelexisting = false;
	 }

  // DESTRUCTOR

  ~model() {}

  // COPY CONSTRUCTOR

  model(const model & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const model & operator=(const model & m);

  // FUNCTION: geterrormessages
  // TASK: returns current errormessages

  const vector<ST::string> & geterrormessages()
//  const errorm::messages & geterrormessages()
	 {
	 return errormessages;
	 }

  // VIRTUAL FUNCTION: parse
  // TASK: basis function for inherited classes
  //       parses the string m
  //       before parsing modeltext,modelVarnames,modelextisting and
  //       errormessages are cleared

  virtual void parse(const ST::string & m)
    {
    clear();
    errormessages.erase(errormessages.begin(),errormessages.end());
    }

  // FUNCTION: makeModelMatrix
  // TASK: creates the designmatrix of the model (stored in datamatrix d)
  // POSSIBLE ERRORS:
  // - variables specified in the model are not existing in dataset ds
  // ADDITIONAL INFORMATION:
  // - if modelVarnames is empty, the datametrix consists of all variables in
  //   the dataset ds

  void makeModelMatrix(dataset & ds,datamatrix & d)
	 {
	 if (modelexisting == true)
		{
		ds.makematrix(modelVarnames,d);
		errormessages = ds.geterrormessages();
		if (! errormessages.empty())
		  clear();
		}
	 }


  // FUNCTION: makeModelMatrix_j
  // TASK: creates j th column of the the designmatrix of the model
  //       (stored in datamatrix d)
  // POSSIBLE ERRORS:
  // - variables specified in the model are not existing in dataset ds
  // ADDITIONAL INFORMATION:
  // - if modelVarnames is empty, the datametrix consists of all variables in
  //   the dataset ds

  void makeModelMatrix_j(dataset & ds,datamatrix & d,const unsigned & j);

  // FUNCTION: getModelVarnames
  // TASK: returns modelVarnames (names of the variables in the model
  //       design matrix)

  const list<ST::string> & getModelVarnames()
	 {
	 return modelVarnames;
	 }

  // FUNCTION: getModelVarnamesAsVector
  // TASK: returns modelVarnames as a vector (names of the variables in the
  //       model design matrix)

  vector<ST::string> getModelVarnamesAsVector();

  // FUNCTION: getModelText
  // TASK: returns current model (modeltext)

  const ST::string & getModelText(void)
	 {
	 return modeltext;
	 }

  };


//------------------------------------------------------------------------------
//------------------------ CLASS: modelStandard --------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE modelStandard : public model

  {

  protected:


  // ------------------------ PROTECTED FUNCTIONS ------------------------------

  // VIRTUAL FUNCTION: clear
  // TASK: clears 'modelVarnames',
  //       modelexisting = false
  //       modeltext = ""

  void clear(void)
	 {
	 model::clear();
	 }


  public:


  // -------------------------- PUBLIC FUNCTIONS -------------------------------

  // DEFAULT CONSTRUCTOR

  modelStandard(void) : model()  {}

  // COPY CONSTRUCTOR

  modelStandard(const modelStandard & m) : model(model(m)) {}

  // OVERLOADED ASSIGNMENT OPERATOR

  const modelStandard & operator=(const modelStandard & m);

  // FUNCTION: parse
  // TASK: parses model 'm'
  //       m wird zerlegt in einzelne token und speichert diese in
  //       modelVarnames. Es findet keine Überprüfung statt, ob die einzelnen
  //       token gültige Variablennamen sind. Handelt es sich bei 'm' um einen
  //       leeren string, so bleibt auch modelVarnames leer. Auch in diesem
  //       Fall wird keine Fehlermeldung zurückgegeben.
  //       modeltexisting wird auf 'true' gesetzt.

  void parse(const ST::string & m);


  };


//------------------------------------------------------------------------------
//--------------------------- CLASS expression ---------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE expression : public model

  {

  protected:


  //------------------------ PROTECTED VARIABLES -------------------------------

  ST::string varname;

  ST::string expr;


  void clear(void)
	 {
     model::clear();
	 varname = "";
	 expr = "";
	 }


  public:


  // ------------------------ PUBLIC FUNCTIONS ---------------------------------

  // DEFAULT CONSTRUCTOR

  expression(void)
	 {
	 modelexisting = false;
	 }

  // COPY CONSTRUCTOR

  expression(const expression & e);

  // OVERLOADED ASSIGNMENT OPERATOR

  const expression & operator=(const expression & e);

  // DESTRUCTOR

  ~expression(void)
	 {
	 }

  // FUNCTION: parse
  // TASK: parses string 'e'

  void parse(const ST::string & e);

  // FUNCTION: getvarname
  // TASK: returns 'varname'

  const ST::string & getvarname(void)
	 {
	 return varname;
	 }

  // FUNCTION: getexpression
  // TASK: returns 'expr'

  const ST::string & getexpression(void)
	 {
	 return expr;
	 }

  // FUNCTION: addvariable
  // TASK: adds a variable with variable name 'varname' and expression 'expr'
  //       to the dataset
  // ADDITIONAL INFORMATION:
  //

  void addvariable(dataset & d)
	 {
	 if (modelexisting == true)
		d.addvariable(varname,expr);
//	 errormessages.insert_back(d.geterrormessages());

	 if (! (d.geterrormessages()).empty())
		errormessages.insert(errormessages.end(),(d.geterrormessages()).begin(),
        (d.geterrormessages()).end());


	 }


  };


//------------------------------------------------------------------------------
//-------------------------------- CLASS: term ---------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE term
  {

  protected:

  // contains current errormessages, if modelspecification is not correct

//  errorm::messages errormessages;
  vector<ST::string> errormessages;

  public:


  ST::string type;

  vector<ST::string> varnames;

  vector<ST::string> options;

  void clear(void)
    {
    if (errormessages.size() > 0)
      errormessages.erase(errormessages.begin(),errormessages.end());
    type="";
    if (varnames.size() > 0)
      varnames.erase(varnames.begin(),varnames.end());
    if (options.size() > 0)
      options.erase(options.begin(),options.end());
    }

  // FUNKTION: parse
  // AUFGABE: parsen des strings 't'
  //          2 Möglichkeiten:
  //            1. t = Variablenname, dann wird 'varname' = t gesetzt.
  //            2. t = Funktion, d.h. funktionsname(option1,option2,...)
  //               dann wird 'varname' = funktionsname gesetzt und
  //               'options' = option1,option2 etc.

  void parse(const ST::string & t);


  const vector<ST::string> & geterrormessages(void)
    {
    return errormessages;
    }

  void print(ostream & out)
    {
    unsigned i;
    out << "Varnames: " << endl;
    out << endl;
    for (i=0;i<varnames.size();i++)
      out << varnames[i] << endl;

    out << endl;
    out << "Options: " << endl;
    out << endl;
    for (i=0;i<options.size();i++)
      out << options[i] << endl;


    }


  term(void) {}

  ~term() {}

  };


//------------------------------------------------------------------------------
//-------------------------- CLASS: basic_termtype -----------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE basic_termtype
  {

  friend term;

  protected:

  ST::string type;

  public:

  // DEFAULT CONSTRUCTOR

  basic_termtype(void)
    {
    type = "basic_termtype";
    }


  virtual bool check(term & t)
    {
    if ( (t.varnames.size() == 1) &&  (t.varnames[0].isvarname()==0)
          && (t.options.size()==0) )
      {
      t.type = "basic_termtype";
      return true;
      }
    else
      return false;
    }

  // FUNKTION: get_constvariables
  // TASK: returns all variables belonging to the basic_termtype

  vector<ST::string> get_constvariables(vector<term> & terms);

  // DESTRUCTOR

  ~basic_termtype() {}

  };


//------------------------------------------------------------------------------
//-------------------------- class term_autoreg  -------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE term_varcoeff : public basic_termtype
  {


  };


class __EXPORT_TYPE term_autoreg : public basic_termtype
  {

  protected:

  intoption min;
  intoption max;
  intoption minvar;
  intoption maxvar;
  doubleoption startv;
  doubleoption lambda;
  doubleoption a;
  doubleoption b;
  stroption proposal;
  intoption updateW;
  simpleoption updatetau;
  doubleoption f;
  doubleoption lambdamin;
  doubleoption lambdamax;
  doubleoption lambdastart;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_autoreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a first or second order random walk

  bool checkvector(const vector<term>  & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( (terms[i].type == "rw1") || (terms[i].type == "rw2") ||
         (terms[i].type == "trw1") || (terms[i].type == "trw2") ||
         (terms[i].type == "rw1vrw1") || (terms[i].type == "rw2vrw1") ||
         (terms[i].type == "rw1vrw2") || (terms[i].type == "rw2vrw2") ||
         (terms[i].type == "varcoeffrw1") || (terms[i].type == "varcoeffrw2")
       )
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_autoreg() {}

  };


//------------------------------------------------------------------------------
//--------------------------- class term_season --------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_varcoeffseason : public basic_termtype
  {

  };

class __EXPORT_TYPE term_season : public basic_termtype
  {


  protected:

  intoption min;
  intoption max;
  intoption period;
  doubleoption lambda;
  doubleoption a;
  doubleoption b;
  stroption proposal;
  intoption updateW;
  simpleoption updatetau;
  doubleoption f;
  doubleoption lambdamin;
  doubleoption lambdamax;
  doubleoption lambdastart;
  simpleoption uniformprior;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_season(void);

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

  ~term_season() {}

  };


//------------------------------------------------------------------------------
//--------------------------- class term_pspline -------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE term_pspline : public basic_termtype
  {


  protected:

  intoption min;
  intoption max;
  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  doubleoption a;
  doubleoption b;
  simpleoption uniformb;
  intoption gridsize;
  intoption minvar;
  intoption maxvar;
  doubleoption startv;
  stroption proposal;
  stroption monotone;
  intoption updateW;
  simpleoption updatetau;
  doubleoption f;
  simpleoption diagtransform;
  simpleoption derivative;
  simpleoption bsplinebasis;
  intoption contourprob;
  simpleoption uniformprior;
  stroption beta_0;
  simpleoption discrete;
  intoption df;
//  doubleoption lambdamin;
//  doubleoption lambdamax;
//  doubleoption lambdastart;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_pspline(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a psline with rw1 oder rw2 Penalty

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( (terms[i].type == "psplinerw1") || (terms[i].type == "psplinerw2") ||
         (terms[i].type == "tpsplinerw1") || (terms[i].type == "tpsplinerw2") ||
         (terms[i].type == "psplinerw1vrw1") || (terms[i].type == "psplinerw1vrw2") ||
         (terms[i].type == "psplinerw2vrw1") || (terms[i].type == "psplinerw2vrw2") )
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_pspline() {}

  };


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//--------------------------- class term_spatial -------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_varcoeff_spatial : public basic_termtype
  {

  };

class __EXPORT_TYPE term_spatial : public basic_termtype
  {


  protected:

  stroption map;
  intoption min;
  intoption max;
  doubleoption lambda;
  doubleoption a;
  doubleoption b;
  stroption proposal;
  intoption updateW;
  simpleoption updatetau;
  doubleoption f;
  doubleoption lambdamin;
  doubleoption lambdamax;
  doubleoption lambdastart;
  simpleoption uniformprior;


  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_spatial(void);

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

  ~term_spatial() {}

  };


//------------------------------------------------------------------------------
//--------------------------- class term_spatial -------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE term_spatialxy : public basic_termtype
  {


  protected:

  intoption min;
  intoption max;
  doubleoption lambda;
  doubleoption a;
  doubleoption b;
  doubleoption maxdist;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_spatialxy(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a spatial component

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( terms[i].type == "spatialxy")
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_spatialxy() {}

  };


//------------------------------------------------------------------------------
//----------------------- class term_geokriging --------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE term_geokriging : public basic_termtype
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

  doubleoption a;
  doubleoption b;
  stroption proposal;
  intoption updateW;
  simpleoption updatetau;
  doubleoption f;
  simpleoption uniformprior;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_geokriging(void);

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

  ~term_geokriging() {}

  };


//------------------------------------------------------------------------------
//------------------------ class term_randomslope  -----------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE term_randomslope : public basic_termtype
  {

  protected:

  simpleoption nofixed;
  doubleoption lambda;
  doubleoption a;
  doubleoption b;
  stroption proposal;
  simpleoption updatetau;
  simpleoption uniformprior;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_randomslope(void);

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

  ~term_randomslope() {}

  };



//------------------------------------------------------------------------------
//--------------------------- class term_random  -------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE term_random : public basic_termtype
  {


  protected:

  doubleoption lambda;
  doubleoption a;
  doubleoption b;
  stroption proposal;
  simpleoption updatetau;
  simpleoption uniformprior;


  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_random(void);

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

  ~term_random() {}

  };


//------------------------------------------------------------------------------
//--------------------------- class term_offset  -------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE term_offset : public basic_termtype
  {


  protected:

  public:

  // DEFAULT CONSTRUCTOR

  term_offset(void)
    {
    type = "term_offset";
    }

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an offset

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( terms[i].type == "offset" )
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_offset() {}

  };


//------------------------------------------------------------------------------
//---------------------- class term_interactpspline  ---------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_interactpspline : public basic_termtype
  {

  protected:

  intoption min;
  intoption max;
  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  simpleoption reduced;
  doubleoption a;
  doubleoption b;
  simpleoption singleblock;
  intoption gridsize;
  stroption proposal;
  intoption updateW;
  simpleoption updatetau;
  doubleoption f;
  simpleoption uniformprior;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_interactpspline(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an interaction term

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_interactpspline() {}

  };


//------------------------------------------------------------------------------
//--------------------------- class term_geospline  ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_geospline : public basic_termtype
  {

  protected:

  intoption min;
  intoption max;
  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  simpleoption reduced;
  stroption map;
  simpleoption singleblock;
  doubleoption a;
  doubleoption b;
  stroption proposal;
  intoption updateW;
  simpleoption updatetau;
  doubleoption f;
  simpleoption uniformprior;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_geospline(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an interaction term

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_geospline() {}

  };


//------------------------------------------------------------------------------
//---------------------- class term_varcoeff_interactpspline  ------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_varcoeff_interactpspline : public basic_termtype
  {

  protected:

  intoption min;
  intoption max;
  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  simpleoption reduced;
  doubleoption a;
  doubleoption b;
  simpleoption singleblock;
  intoption gridsize;
  stroption proposal;
  intoption updateW;
  simpleoption updatetau;
  doubleoption f;
  simpleoption uniformprior;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_varcoeff_interactpspline(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an interaction term

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_varcoeff_interactpspline() {}

  };


//------------------------------------------------------------------------------
//--------------------------- class term_varcoeff_geospline  -------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_varcoeff_geospline : public basic_termtype
  {

  protected:

  intoption min;
  intoption max;
  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  simpleoption reduced;
  stroption map;
  simpleoption singleblock;
  doubleoption a;
  doubleoption b;
  stroption proposal;
  intoption updateW;
  simpleoption updatetau;
  doubleoption f;
  simpleoption uniformprior;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_varcoeff_geospline(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an interaction term

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_varcoeff_geospline() {}

  };


//------------------------------------------------------------------------------
//------------------------ class term_varcoeff_pspline  ------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE term_varcoeff_pspline : public basic_termtype
  {


  protected:

  intoption min;
  intoption max;
  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  intoption gridsize;
  doubleoption a;
  doubleoption b;
  stroption proposal;
  stroption monotone;
  intoption updateW;
  simpleoption updatetau;
  doubleoption f;
  simpleoption diagtransform;
  simpleoption derivative;
  intoption contourprob;
  simpleoption uniformprior;
  stroption beta_0;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_varcoeff_pspline(void);

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

  ~term_varcoeff_pspline() {}

  };

//------------------------------------------------------------------------------
//------------------------- class term_baseline -------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE term_baseline : public basic_termtype
  {


  protected:

  intoption min;
  intoption max;
  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  doubleoption a;
  doubleoption b;
  simpleoption uniformb;
  intoption gridsize;
  simpleoption uniformprior;


  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_baseline(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a psline with rw1 oder rw2 Penalty

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if (terms[i].type == "baseline")

      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_baseline() {}

  };


//------------------------------------------------------------------------------
//------------------------ class term_varcoeff_baseline  ------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE term_varcoeff_baseline : public basic_termtype
  {


  protected:

  intoption min;
  intoption max;
  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  intoption gridsize;
  doubleoption a;
  doubleoption b;
  simpleoption uniformprior;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_varcoeff_baseline(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an interaction term

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if (terms[i].type == "varbaseline")
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_varcoeff_baseline() {}

  };



//------------------------------------------------------------------------------
//--------------------------- CLASS: modelterm ---------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE modelterm : public model
  {

  protected:

  vector<basic_termtype*> * termtypes;

  vector <term> terms;

  ST::string responsevar;

  void clear(void)
    {
    model::clear();
    responsevar="";
    }

  public:

  // DEFAULT CONSTRUCTOR

  modelterm(void)
    {
    modelexisting = false;
    }

  // CONSTRUCTOR

  modelterm(vector<basic_termtype*> * t) : model()
    {
    termtypes = t;
    }

  // COPY CONSTRUCTOR

  modelterm(const modelterm & m) : model(model(m))
    {
    responsevar = m.responsevar;
    termtypes = m.termtypes;
    terms = m.terms;
    }

  const modelterm & operator=(const modelterm & m);

  // FUNKTION: parse
  // AUFGABE:  parsen des Modells m
  //           Syntax von m: responsevariable = term1 + term2 + ... + termn
  //           Syntax eines Terms: unterterm1 * unterterm2 * unterterm3
  //           Syntax eines Unterterms funktions/variablenname(option1,option2,...)
  //           alle möglichen Unterterme sind in 'termtypes' gespeichert
  //           Responsevariable wird in 'responsevar' gespeichert
  //           Terme mit ihren Untertermen werden in 'terms' gespeichert

  void parse(const ST::string & m);

  vector<term> getterms(void)
    {
    return terms;
    }

  const ST::string & getresponsevar(void)
    {
    return responsevar;
    }

  ~modelterm() {}

  };


//------------------------------------------------------------------------------
//--------------------------- CLASS: modelterm ---------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE modeltermmult : public model
  {

  protected:

  vector<basic_termtype*> * termtypes;

  vector < vector <term> > terms;

  vector<ST::string> responsevar;

  vector<unsigned> responsecol;

  void clear(void);

  public:

  // DEFAULT CONSTRUCTOR

  modeltermmult(void)
    {
    modelexisting = false;
    }

  // CONSTRUCTOR

  modeltermmult(vector<basic_termtype*> * t) : model()
    {
    termtypes = t;
    }

  // COPY CONSTRUCTOR

  modeltermmult(const modeltermmult & m) : model(model(m))
    {
    responsevar = m.responsevar;
    responsecol = m.responsecol;
    termtypes = m.termtypes;
    terms = m.terms;
    }

  const modeltermmult & operator=(const modeltermmult & m);

  // FUNKTION: parse
  // AUFGABE:  parsen des Modells m
  //           Syntax von m: responsevariable1 = term1 + term2 + ... + termn :
  //                         responsevariable2 = term1 + term2 + ... + termn :
  //                                 .
  //                                 .
  //                                 .
  //           Syntax eines Terms: unterterm1 * unterterm2 * unterterm3
  //           Syntax eines Unterterms funktions/variablenname(option1,option2,...)
  //           alle möglichen Unterterme sind in 'termtypes' gespeichert
  //           Responsevariable wird in 'responsevar' gespeichert
  //           Terme mit ihren Untertermen werden in 'terms' gespeichert

  void parse(const ST::string & m);

  vector<term> getterms(const unsigned & i)
    {
    return terms[i];
    }

  vector< vector<term> > getterms(void)
    {
    return terms;
    }

  const ST::string & getresponsevar(const unsigned & i)
    {
    return responsevar[i];
    }

  const vector<unsigned> & getresponsecol(void)
    {
    return responsecol;
    }

  ~modeltermmult() {}

  };


#endif
