
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

    while ( (j<termnames.size()) && (f == false) )
      {
      if (t.options[0] == termnames[j])
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





