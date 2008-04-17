


#include"model_parameters.h"
#include<algorithm>

//------------------------------------------------------------------------------
//-------------- CLASS model: implementation of member functions ---------------
//------------------------------------------------------------------------------


model::model(const model & m)
  {                                                  
  modelexisting = m.modelexisting;
  modeltext = m.modeltext;
  errormessages = m.errormessages;
  modelVarnames = m.modelVarnames;
  }


const model & model::operator=(const model & m)
  {
  if (this == &m)
	 return *this;
  modelexisting = m.modelexisting;
  modeltext = m.modeltext;
  errormessages = m.errormessages;
  modelVarnames = m.modelVarnames;
  return *this;
  }


vector<ST::string> model::getModelVarnamesAsVector()
  {
  vector<ST::string> help;
  if (! modelVarnames.empty())
    {
    list<ST::string>::iterator it;
    for (it = modelVarnames.begin();it != modelVarnames.end();++it)
      help.push_back(*it);
    }
  return help;
  }

void model::makeModelMatrix_j(dataset & ds,datamatrix & d,const unsigned & j)
  {
  if (modelexisting == true)
    {
//    int j2 = j;
//    list<ST::string>::iterator it = modelVarnames.begin()+j2;
//    ST::string var_j = *it;
    vector<ST::string> varn = getModelVarnamesAsVector();
    ST::string var_j = varn[j];
    ds.makematrix(var_j,d);
    errormessages = ds.geterrormessages();
    if (! errormessages.empty())
      clear();
    }
  }



//------------------------------------------------------------------------------
//---------- CLASS modelStandard: implementation of member functions -----------
//------------------------------------------------------------------------------


const modelStandard & modelStandard::operator=(const modelStandard & m)
  {
  if (this == &m)
	 return *this;
  model::operator=(model(m));
  return *this;
  }


void modelStandard::parse (const ST::string & m)

  {
  model::parse(m);
  modelVarnames = m.strtokenlist(" ");
  modelexisting =  true;
  modeltext = m;
  }


//------------------------------------------------------------------------------
//------------ CLASS expression: implementation of member functions ------------
//------------------------------------------------------------------------------


expression::expression(const expression & e) : model(model(e))
  {
  varname = e.varname;
  expr = e.expr;
  }


const expression & expression::operator=(const expression & e)
  {
  if (this == &e)
	 return *this;
  model::operator=(model(e));
  varname = e.varname;
  expr = e.expr;
  return *this;
  }


void expression::parse(const ST::string & e)
  {

  model::parse(e);

  int equalsignpos = e.checksign('=');
  if (equalsignpos == -1)
	 errormessages.push_back("ERROR: \"=\" expected\n");
  else if (e.length() <= equalsignpos+1)
	 errormessages.push_back("ERROR: expression expected\n");
  else
	 {
     if (equalsignpos > 0)
       {
	   varname = e.substr(0,equalsignpos);
	   varname = varname.eatwhitespace();
	   if (varname.isvarname() == 1)
		 errormessages.push_back("ERROR: " + varname + " invalid varname\n");
       }
     else
       {
       errormessages.push_back("ERROR: new varname expected\n");
       }

     if (e.length()-equalsignpos-1>0)
       {
	   expr = e.substr(equalsignpos+1,e.length()-equalsignpos-1);
	   expr = expr.eatwhitespace();
       }
     else
       errormessages.push_back("ERROR: expression expected\n");  

	 }

  if (errormessages.empty())
	 {
	 modelexisting = true;
	 modeltext = e;
	 }
  else
	 clear();

  }


//------------------------------------------------------------------------------
//--------------- CLASS term: implementation of member functions ---------------
//------------------------------------------------------------------------------


void term::parse(const ST::string & t)
  {

  clear();

  ST::string te;
  te = t.eatallwhitespace();
  if (te.length() == 0)
    {
    errormessages.push_back("ERROR: invalid term specification");
    }
  else
    {
    ST::string functionname;
    ST::string argument;
    int isfunc = te.isfunction(functionname,argument);
    if (isfunc == 0)
      {
      vector<ST::string> token = te.strtoken(" *",false);
      unsigned i;
      for(i=0;i<token.size();i++)
        {
        if (token[i].isvarname() == 0)
          varnames.push_back(token[i]);
        else
          errormessages.push_back("ERROR: " + token[i] +
                                  " is not a valid varname\n");
        }
      }
    else if (isfunc==-1)
      errormessages.push_back("ERROR: missing bracket(s) in " + te + "\n");
    else
      {
      vector<ST::string> token = functionname.strtoken(" *",false);
      unsigned i;
      for(i=0;i<token.size();i++)
        {
        if (token[i].isvarname() == 0)
          varnames.push_back(token[i]);
        else
          errormessages.push_back("ERROR: " + token[i] +
                                  " is not a valid varname\n");
        }

      if (argument.length() > 0)
        options = argument.strtoken(",",false);

      }
    }
  }


//------------------------------------------------------------------------------
//---------- CLASS basic_termtype: implementation of member functions ----------
//------------------------------------------------------------------------------

vector<ST::string> basic_termtype::get_constvariables(vector<term> & terms)
  {
  vector<ST::string> res;
  unsigned i;
  for(i=0;i<terms.size();i++)
    {
    if (terms[i].type == "basic_termtype")
      res.push_back(terms[i].varnames[0]);

    }
  return res;
  }

//------------------------------------------------------------------------------
//----------- class term_nonp: implementation of member functions --------------
//------------------------------------------------------------------------------

term_nonp::term_nonp(void)
  {
  type = "term_psline";
  degree=intoption("degree",3,0,5);
  numberknots=intoption("nrknots",20,5,500);
  difforder =  intoption("difforder",2,1,3);
  lambda = doubleoption("lambda",0.1,0,10000000);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  }

void term_nonp::setdefault(void)
  {
  degree.setdefault();
  numberknots.setdefault();
  difforder.setdefault();
  lambda.setdefault();
  a.setdefault();
  b.setdefault();
  }


bool term_nonp::checkvector(const vector<term> & terms,const unsigned & i)
  {
  assert(i< terms.size());

  bool f = false;
  unsigned j;
  while ( (j<termnames.size()) && (f == false) )
    {
    if terms[i].type == ternames[j]
      {
      f = true;
      }
    j ++;
    }

  return f;
  }


bool term_nonp::check(term & t)
  {

  if ( (t.varnames.size()<=2)  && (t.options.size() >= 1)
        && (t.options.size() <= 50) )
    {

    bool f = false;
    unsigned j;

    while ( (j<termnames.size()) && (f == false) )
      {
      if (t.options[0] == ternames[j])
        {
        f = true;
        }
      j ++;
      }

    if (f==false)
      {
      setdefault();
      return false;
      }

    long minim,maxim;

    optionlist optlist;
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&difforder);
    optlist.push_back(&lambda);
    optlist.push_back(&a);
    optlist.push_back(&b);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(50);
    t.options[0] = t.type;
    t.options[1] = ST::inttostring(degree.getvalue());
    t.options[2] = ST::inttostring(numberknots.getvalue());
    t.options[3] = ST::inttostring(difforder.getvalue());
    t.options[4] = ST::doubletostring(lambda.getvalue());
    t.options[5] = ST::doubletostring(a.getvalue());
    t.options[6] = ST::doubletostring(b.getvalue());

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }





