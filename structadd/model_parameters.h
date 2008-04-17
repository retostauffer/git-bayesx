#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (MODEL_INCLUDED)
#define MODEL_INCLUDED

#include"clstring.h"
#include<vector>
#include"data.h"
#include"option.h"



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

  friend class term;

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

  vector<ST::string> termnames;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_nonp(vector<ST::string> na);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a nonparametric term

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_nonp() {}

  };

#endif
