
#include"model_parameters.h"
#include<algorithm>


//------------------------------------------------------------------------------
//----------- class term_nonp: implementation of member functions --------------
//------------------------------------------------------------------------------

term_nonp::term_nonp(vector<ST::string> & na)
  {
  termnames = na;
  degree=intoption("degree",3,0,5);
  numberknots=intoption("nrknots",20,5,500);
  difforder =  intoption("difforder",2,1,3);
  lambda = doubleoption("lambda",0.1,0,10000000);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  nocenter = simpleoption("nocenter",false);
  map=stroption("map");
  lambda_re = doubleoption("lambda_re",0.1,0,10000000);
  a_re = doubleoption("a_re",0.001,-1.0,500);
  b_re = doubleoption("b_re",0.001,0,500);
  internal_mult = simpleoption("internal_mult",false);
  samplemult = simpleoption("samplemult",false);
  vector<ST::string> ctypes;
  ctypes.push_back("unconstrained");
  ctypes.push_back("increasing");
  ctypes.push_back("decreasing");
  constraints = stroption("constraints",ctypes,"unconstrained");
  round = doubleoption("round",-1,0,500);
  vector<ST::string> centermethods;
  centermethods.push_back("mean");
  centermethods.push_back("meanintegral");
  centermethods.push_back("meaninvvar");
  centermethods.push_back("nullspace");
  centermethods.push_back("meansimple");
  centermethod = stroption("centermethod",centermethods,"mean");
  internal_multexp = simpleoption("internal_multexp",false);
  pvalue = simpleoption("pvalue",false);
  meaneffect = simpleoption("meaneffect",false);
  }

void term_nonp::setdefault(void)
  {
  degree.setdefault();
  numberknots.setdefault();
  difforder.setdefault();
  lambda.setdefault();
  a.setdefault();
  b.setdefault();
  nocenter.setdefault();
  map.setdefault();
  lambda_re.setdefault();
  a_re.setdefault();
  b_re.setdefault();
  internal_mult.setdefault();
  samplemult.setdefault();
  constraints.setdefault();
  round.setdefault();
  centermethod.setdefault();
  internal_multexp.setdefault();
  pvalue.setdefault();
  meaneffect.setdefault();
  }


bool term_nonp::checkvector(const vector<term> & terms,const unsigned & i)
  {
  assert(i< terms.size());

  bool f = false;
  unsigned j=0;
  while ( (j<termnames.size()) && (f == false) )
    {
    if (terms[i].type == termnames[j])
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
    unsigned j=0;
    unsigned namespos;

    while ( (j<termnames.size()) && (f == false) )
      {
      if (t.options[0] == termnames[j])
        {
        namespos=j;
        f = true;
        }
      j ++;
      }

    if (f==false)
      {
      setdefault();
      return false;
      }

    optionlist optlist;
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&difforder);
    optlist.push_back(&lambda);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&nocenter);
    optlist.push_back(&map);
    optlist.push_back(&lambda_re);
    optlist.push_back(&a_re);
    optlist.push_back(&b_re);
    optlist.push_back(&internal_mult);
    optlist.push_back(&samplemult);
    optlist.push_back(&constraints);
    optlist.push_back(&round);
    optlist.push_back(&centermethod);
    optlist.push_back(&internal_multexp);
    optlist.push_back(&pvalue);
    optlist.push_back(&meaneffect);

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
    t.options[0] = termnames[namespos];
    t.options[1] = ST::inttostring(degree.getvalue());
    t.options[2] = ST::inttostring(numberknots.getvalue());
    t.options[3] = ST::inttostring(difforder.getvalue());
    t.options[4] = ST::doubletostring(lambda.getvalue());
    t.options[5] = ST::doubletostring(a.getvalue());
    t.options[6] = ST::doubletostring(b.getvalue());

    if(nocenter.getvalue() == false)
      t.options[7] = "false";
    else
      t.options[7] = "true";

    t.options[8] = map.getvalue();

    t.options[9] = ST::doubletostring(lambda_re.getvalue());

    t.options[10] = ST::doubletostring(a_re.getvalue());

    t.options[11] = ST::doubletostring(b_re.getvalue());

    if(internal_mult.getvalue() == false)
      t.options[12] = "false";
    else
      t.options[12] = "true";

    if(samplemult.getvalue() == false)
      t.options[13] = "false";
    else
      t.options[13] = "true";

    t.options[14] = constraints.getvalue();

    t.options[15] = ST::doubletostring(round.getvalue());

    t.options[16] = centermethod.getvalue();

    if(internal_mult.getvalue() == false)
      t.options[17] = "false";
    else
      t.options[17] = "true";

    if(pvalue.getvalue() == false)
      t.options[18] = "false";
    else
      t.options[18] = "true";

    if(meaneffect.getvalue() == false)
      t.options[19] = "false";
    else
      t.options[19] = "true";

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }



