
#include "model_stepwise.h"


//------------------------------------------------------------------------------
//---------- class term_nonlinearf: implementation of member functions --------
//------------------------------------------------------------------------------

term_nonlinearf_stepwise::term_nonlinearf_stepwise(void)
  {
  type = "term_nonlinearf";
  lambda = doubleoption("lambda",-1,-1,10);
  lambdastart = doubleoption("lambdastart",-1,-1,10);
  forced_into = simpleoption("forced_into",false);
  }


void term_nonlinearf_stepwise::setdefault(void)
  {
  lambda.setdefault();
  lambdastart.setdefault();
  forced_into.setdefault();
  }


bool term_nonlinearf_stepwise::check(term & t)
  {
  if ( (t.varnames.size()<=1)  && (t.options.size()<=4) &&
       (t.options.size() >= 1) )
    {

    if (t.options[0] == "nonlinearf")
      t.type = "nonlinearf";
    else
      {
      setdefault();
      return false;
      }


    vector<ST::string> opt;
    optionlist optlist;

    optlist.push_back(&lambda);
    optlist.push_back(&lambdastart);
    optlist.push_back(&forced_into);

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
    t.options = vector<ST::string>(4);
    t.options[0] = t.type;
    t.options[1] = ST::doubletostring(lambda.getvalue());
    t.options[2] = ST::doubletostring(lambdastart.getvalue());
    if (forced_into.getvalue()==false)
      t.options[3] = "false";
    else
      t.options[3] = "true";

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


//------------------------------------------------------------------------------
//---------- class term_autoreg: implementation of member functions ------------
//------------------------------------------------------------------------------


term_autoreg_stepwise::term_autoreg_stepwise(void)
  {
  type = "term_autoreg";
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdamin = doubleoption("lambdamin",0.0001,0.000001,100000000);
  lambdamax = doubleoption("lambdamax",10000,0.000001,100000000);
  lambdastart = doubleoption("lambdastart",-1,-1,100000000);
  forced_into = simpleoption("forced_into",false);
  df_for_lambdamax = doubleoption("df_for_lambdamax",2,0,200);         // Unterscheidung rw1/rw2 bei default-Wert!!!
  df_for_lambdamin = doubleoption("df_for_lambdamin",10,0,200);
  lambdamax_opt = simpleoption("lambdamax_opt",false);
  lambdamin_opt = simpleoption("lambdamin_opt",false);
  number = intoption("number",0,0,50);
  df_equidist = simpleoption("df_equidist",false);
  df_accuracy = doubleoption("df_accuracy",0.05,0.01,0.5);
  // df_start =  = doubleoption("df_start",-1,-1,100);
  }

void term_autoreg_stepwise::setdefault(void)
  {
  lambda.setdefault();
  lambdamin.setdefault();
  lambdamax.setdefault();
  lambdastart.setdefault();
  forced_into.setdefault();
  df_for_lambdamax.setdefault();
  df_for_lambdamin.setdefault();
  lambdamax_opt.setdefault();
  lambdamin_opt.setdefault();
  number.setdefault();
  df_equidist.setdefault();
  df_accuracy.setdefault();
  // df_start.setdefault();
  }


bool term_autoreg_stepwise::check(term & t)
  {

  if ( (t.varnames.size() <= 2) && (t.varnames.size() >= 1) &&
       (t.options.size()<=13) && (t.options.size() >= 1) )
    {

    if (t.options[0] == "rw1" && t.varnames.size() == 1)
      t.type = "rw1";
    else if (t.options[0] == "rw2" && t.varnames.size() == 1)
      t.type = "rw2";
    else if (t.options[0] == "rw1" && t.varnames.size() == 2)
      t.type = "varcoeffrw1";
    else if (t.options[0] == "rw2" && t.varnames.size() == 2)
      t.type = "varcoeffrw2";
    else
      {
      setdefault();
      return false;
      }

    double minl, maxl, startl, df_max, df_min;

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&lambda);
    optlist.push_back(&lambdamin);
    optlist.push_back(&lambdamax);
    optlist.push_back(&lambdastart);
    optlist.push_back(&forced_into);
    optlist.push_back(&df_for_lambdamax);
    optlist.push_back(&df_for_lambdamin);
    optlist.push_back(&lambdamax_opt);
    optlist.push_back(&lambdamin_opt);
    optlist.push_back(&number);
    optlist.push_back(&df_equidist);
    optlist.push_back(&df_accuracy);
    // optlist.push_back(&df_start);

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

    optlist.push_back(&lambda);
    optlist.push_back(&lambdamin);
    optlist.push_back(&lambdamax);
    optlist.push_back(&lambdastart);
    optlist.push_back(&forced_into);
    optlist.push_back(&df_for_lambdamax);
    optlist.push_back(&df_for_lambdamin);
    optlist.push_back(&lambdamax_opt);
    optlist.push_back(&lambdamin_opt);
    optlist.push_back(&number);
    optlist.push_back(&df_equidist);
    optlist.push_back(&df_accuracy);
    // optlist.push_back(&df_start);


    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(13);
    t.options[0] = t.type;
    t.options[1] = ST::doubletostring(lambda.getvalue());
    t.options[2] = ST::doubletostring(lambdamin.getvalue());
    t.options[3] = ST::doubletostring(lambdamax.getvalue());
    t.options[4] = ST::doubletostring(lambdastart.getvalue());
    if (forced_into.getvalue()==false)
        t.options[5] = "false";
    else
        t.options[5] = "true";
    t.options[6] = ST::doubletostring(df_for_lambdamax.getvalue());
    t.options[7] = ST::doubletostring(df_for_lambdamin.getvalue());
    if (lambdamax_opt.getvalue()==false)
        t.options[8] = "false";
    else
        t.options[8] = "true";
    if (lambdamin_opt.getvalue()==false)
        t.options[9] = "false";
    else
        t.options[9] = "true";
    t.options[10] = ST::inttostring(number.getvalue());
    if (df_equidist.getvalue()==false)
        t.options[11] = "false";
    else
        t.options[11] = "true";
    t.options[12] = ST::doubletostring(df_accuracy.getvalue());


    int b = t.options[2].strtodouble(minl);
    b = t.options[3].strtodouble(maxl);
    b = t.options[4].strtodouble(startl);
    b = t.options[6].strtodouble(df_max);
    b = t.options[7].strtodouble(df_min);

    if (b==1)
      {
      setdefault();
      return false;
      }

    if (minl >= maxl)
      {
      setdefault();
      return false;
      }

    if (maxl < startl)
      {
      setdefault();
      return false;
      }

    if(df_max>=df_min)
      {
      setdefault();
      return false;
      }

    if ((df_min <=1 && t.options[0] == "rw2") ||            // kann man auf Anzahl der Parameter zugreifen???
            (df_max <=1 && t.options[0] == "rw2"))
      {
      setdefault();
      return false;
      }

    if ((df_min ==1 && t.options[0] == "rw1") ||            // kann man auf Anzahl der Parameter zugreifen???
            (df_max ==1 && t.options[0] == "rw1"))
      {
      setdefault();
      return false;
      }

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


//------------------------------------------------------------------------------
//---------- class term_season: implementation of member functions ------------
//------------------------------------------------------------------------------

term_season_stepwise::term_season_stepwise(void)
  {
  type = "term_season";
  period = intoption("period",12,2,72);
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdamin = doubleoption("lambdamin",0.0001,0.000001,10000000);
  lambdamax = doubleoption("lambdamax",10000,0.000001,10000000);
  lambdastart = doubleoption("lambdastart",10000,0,10000000);
  forced_into = simpleoption("forced_into",false);
  df_for_lambdamax = doubleoption("df_for_lambdamax",11,0,200);         // Unterscheidung rw1/rw2 bei default-Wert!!!
  df_for_lambdamin = doubleoption("df_for_lambdamin",15,0,200);
  lambdamax_opt = simpleoption("lambdamax_opt",false);
  lambdamin_opt = simpleoption("lambdamin_opt",false);
  number = intoption("number",0,0,50);
  df_equidist = simpleoption("df_equidist",false);
  df_accuracy = doubleoption("df_accuracy",0.05,0.01,0.5);
  }


void term_season_stepwise::setdefault(void)
  {
  period.setdefault();
  lambda.setdefault();
  lambdamin.setdefault();
  lambdamax.setdefault();
  lambdastart.setdefault();
  forced_into.setdefault();
  df_for_lambdamax.setdefault();
  df_for_lambdamin.setdefault();
  lambdamax_opt.setdefault();
  lambdamin_opt.setdefault();
  number.setdefault();
  df_equidist.setdefault();
  df_accuracy.setdefault();
  }


bool term_season_stepwise::check(term & t)
  {
  if ( (t.varnames.size()<=2)  && (t.options.size()<=14) &&
       (t.options.size() >= 1) )
    {

    if ( (t.options[0] == "season") && (t.varnames.size() == 1) )
      t.type = "season";
    else if( (t.options[0] == "season") && (t.varnames.size() == 2) )
      t.type = "varcoeffseason";
    else
      {
      setdefault();
      return false;
      }

    long per;
    double minl, maxl, df_max, df_min, pe;

    vector<ST::string> opt;
    optionlist optlist;

    optlist.push_back(&period);
    optlist.push_back(&lambda);
    optlist.push_back(&lambdamin);
    optlist.push_back(&lambdamax);
    optlist.push_back(&lambdastart);
    optlist.push_back(&forced_into);
    optlist.push_back(&df_for_lambdamax);
    optlist.push_back(&df_for_lambdamin);
    optlist.push_back(&lambdamax_opt);
    optlist.push_back(&lambdamin_opt);
    optlist.push_back(&number);
    optlist.push_back(&df_equidist);
    optlist.push_back(&df_accuracy);

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
    t.options = vector<ST::string>(14);
    t.options[0] = t.type;
    t.options[1] = ST::inttostring(period.getvalue());
    t.options[2] = ST::doubletostring(lambda.getvalue());
    t.options[3] = ST::doubletostring(lambdamin.getvalue());
    t.options[4] = ST::doubletostring(lambdamax.getvalue());
    t.options[5] = ST::doubletostring(lambdastart.getvalue());
    if (forced_into.getvalue()==false)
      t.options[6] = "false";
    else
      t.options[6] = "true";
    t.options[7] = ST::doubletostring(df_for_lambdamax.getvalue());
    t.options[8] = ST::doubletostring(df_for_lambdamin.getvalue());
    if (lambdamax_opt.getvalue()==false)
        t.options[9] = "false";
    else
        t.options[9] = "true";
    if (lambdamin_opt.getvalue()==false)
        t.options[10] = "false";
    else
        t.options[10] = "true";
    t.options[11] = ST::inttostring(number.getvalue());
    if(df_equidist.getvalue()==false)
        t.options[12] = "false";
    else
        t.options[12] = "true";
    t.options[13] = ST::doubletostring(df_accuracy.getvalue());


    if (t.options[1].strtolong(per) == 1)
      {
      setdefault();
      return false;
      }

    int b = t.options[3].strtodouble(minl);
    b = t.options[4].strtodouble(maxl);
    b = t.options[7].strtodouble(df_max);
    b = t.options[8].strtodouble(df_min);
    b = t.options[1].strtodouble(pe);

    if (b==1)
      {
      setdefault();
      return false;
      }

    if (minl >= maxl)
      {
      setdefault();
      return false;
      }

    if (df_min < (pe-2) ||            // kann man auf Anzahl der Parameter zugreifen???
            df_max < (pe-2))
      {
      setdefault();
      return false;
      }

    if(df_max>=df_min)
      {
      setdefault();
      return false;
      }

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


//------------------------------------------------------------------------------
//----------- class term_pspline: implementation of member functions -----------
//------------------------------------------------------------------------------

term_pspline_stepwise::term_pspline_stepwise(void)
  {
  type = "term_psline";
  degree=intoption("degree",3,0,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  gridsize = intoption("gridsize",-1,10,500);
  diagtransform = simpleoption("diagtransform",false);
  derivative = simpleoption("derivative",false);
  lambdamin = doubleoption("lambdamin",0.0001,0.000001,100000000);
  lambdamax = doubleoption("lambdamax",10000,0.000001,100000000);
  lambdastart = doubleoption("lambdastart",-1,-1,100000000);
  forced_into = simpleoption("forced_into",false);
  df_for_lambdamax = doubleoption("df_for_lambdamax",2,0,200);         // Unterscheidung rw1/rw2 bei default-Wert!!!
  df_for_lambdamin = doubleoption("df_for_lambdamin",10,0,200);
  lambdamax_opt = simpleoption("lambdamax_opt",false);
  lambdamin_opt = simpleoption("lambdamin_opt",false);
  number = intoption("number",0,-1,50);
  df_equidist = simpleoption("df_equidist",false);
  df_accuracy = doubleoption("df_accuracy",0.05,0.01,0.5);
  }

void term_pspline_stepwise::setdefault(void)
  {
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  gridsize.setdefault();
  diagtransform.setdefault();
  derivative.setdefault();
  lambdamin.setdefault();
  lambdamax.setdefault();
  lambdastart.setdefault();
  forced_into.setdefault();
  df_for_lambdamax.setdefault();
  df_for_lambdamin.setdefault();
  lambdamax_opt.setdefault();
  lambdamin_opt.setdefault();
  number.setdefault();
  df_equidist.setdefault();
  df_accuracy.setdefault();
  }

bool term_pspline_stepwise::check(term & t)
  {

  if ( (t.varnames.size()<=2)  && (t.varnames.size() >= 1)
       && (t.options.size() >= 1) && (t.options.size() <= 18) )
    {

    if (t.options[0] == "psplinerw1" && t.varnames.size() == 1)
      t.type = "psplinerw1";
    else if (t.options[0] == "psplinerw2" && t.varnames.size() == 1)
      t.type = "psplinerw2";
    else if (t.options[0] == "psplinerw1" && t.varnames.size() == 2)
      t.type = "varpsplinerw1";
    else if (t.options[0] == "psplinerw2" && t.varnames.size() == 2)
      t.type = "varpsplinerw2";
    else
      {
      setdefault();
      return false;
      }

    double maxl,minl,startl,df_max,df_min;

    optionlist optlist;
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&lambda);
    optlist.push_back(&gridsize);
    optlist.push_back(&diagtransform);
    optlist.push_back(&derivative);
    optlist.push_back(&lambdamin);
    optlist.push_back(&lambdamax);
    optlist.push_back(&lambdastart);
    optlist.push_back(&forced_into);
    optlist.push_back(&df_for_lambdamax);
    optlist.push_back(&df_for_lambdamin);
    optlist.push_back(&lambdamax_opt);
    optlist.push_back(&lambdamin_opt);
    optlist.push_back(&number);
    optlist.push_back(&df_equidist);
    optlist.push_back(&df_accuracy);

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
   t.options = vector<ST::string>(18);
   t.options[0] = t.type;
   t.options[1] = ST::inttostring(degree.getvalue());
   t.options[2] = ST::inttostring(numberknots.getvalue());
   t.options[3] = ST::doubletostring(lambda.getvalue());
   t.options[4] = ST::inttostring(gridsize.getvalue());
   if (diagtransform.getvalue() == false)
     t.options[5] = "false";
   else
     t.options[5] = "true";
   if (derivative.getvalue() == false)
     t.options[6] = "false";
   else
     t.options[6] = "true";
    t.options[7] = ST::doubletostring(lambdamin.getvalue());
    t.options[8] = ST::doubletostring(lambdamax.getvalue());
    t.options[9] = ST::doubletostring(lambdastart.getvalue());
    if (forced_into.getvalue()==false)
      t.options[10] = "false";
    else
      t.options[10] = "true";
    t.options[11] = ST::doubletostring(df_for_lambdamax.getvalue());
    t.options[12] = ST::doubletostring(df_for_lambdamin.getvalue());
    if (lambdamax_opt.getvalue()==false)
        t.options[13] = "false";
    else
        t.options[13] = "true";
    if (lambdamin_opt.getvalue()==false)
        t.options[14] = "false";
    else
        t.options[14] = "true";
    t.options[15] = ST::inttostring(number.getvalue());
    if (df_equidist.getvalue()==false)
       t.options[16] = "false";
    else
       t.options[16] = "true";
    t.options[17] = ST::doubletostring(df_accuracy.getvalue());

      
   if (lambda.getvalue() < 0)
     {
     setdefault();
     return false;
     }

    int b = t.options[7].strtodouble(minl);
    b = t.options[8].strtodouble(maxl);
    b = t.options[9].strtodouble(startl);
    b = t.options[11].strtodouble(df_max);
    b = t.options[12].strtodouble(df_min);

    if (b==1)
      {
      setdefault();
      return false;
      }


    if (minl >= maxl)
      {
      setdefault();
      return false;
      }

    if (maxl < startl)
      {
      setdefault();
      return false;
      }

    if ((df_min <= 1 && t.options[0] == "psplinerw2") ||            // kann man auf Anzahl der Parameter zugreifen???
            (df_max <= 1 && t.options[0] == "psplinerw2"))
      {
      setdefault();
      return false;
      }

    if ((df_min == 1 && t.options[0] == "psplinerw1") ||            // kann man auf Anzahl der Parameter zugreifen???
            (df_max == 1 && t.options[0] == "psplinerw1"))
      {
      setdefault();
      return false;
      }

    if(df_max>=df_min)
      {
      setdefault();
      return false;
      }

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


//------------------------------------------------------------------------------
//----------- class term_spatial: implementation of member functions -----------
//------------------------------------------------------------------------------

term_spatial_stepwise::term_spatial_stepwise(void)
  {
  type = "term_spatial";
  map=stroption("map");
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdamin = doubleoption("lambdamin",0.0001,0.000001,10000000);
  lambdamax = doubleoption("lambdamax",10000,0.000001,10000000);
  lambdastart = doubleoption("lambdastart",10000,-1,10000000);
  forced_into = simpleoption("forced_into",false);
  df_for_lambdamax = doubleoption("df_for_lambdamax",1,0,500);         // Unterscheidung rw1/rw2 bei default-Wert!!!
  df_for_lambdamin = doubleoption("df_for_lambdamin",10,0,500);
  lambdamax_opt = simpleoption("lambdamax_opt",false);
  lambdamin_opt = simpleoption("lambdamin_opt",false);
  number = intoption("number",0,0,50);
  df_equidist = simpleoption("df_equidist",false);
  df_accuracy = doubleoption("df_accuracy",0.05,0.01,0.5);
  }

void term_spatial_stepwise::setdefault(void)
  {
  map.setdefault();
  lambda.setdefault();
  lambdamin.setdefault();
  lambdamax.setdefault();
  lambdastart.setdefault();
  forced_into.setdefault();
  df_for_lambdamax.setdefault();
  df_for_lambdamin.setdefault();
  lambdamax_opt.setdefault();
  lambdamin_opt.setdefault();
  number.setdefault();
  df_equidist.setdefault();
  df_accuracy.setdefault();
  }


bool term_spatial_stepwise::check(term & t)
  {

  if ( (t.varnames.size()<=2)  && (t.varnames.size()>=1) &&
       (t.options.size()<=14) && (t.options.size() >= 1) )
    {

    if (t.options[0] == "spatial" && t.varnames.size()==1)
      t.type = "spatial";
    else  if (t.options[0] == "spatial" && t.varnames.size()==2)
      t.type = "varcoeffspatial";
    else
      {
      setdefault();
      return false;
      }

    double minl,maxl,startl,df_max,df_min;

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&map);
    optlist.push_back(&lambda);
    optlist.push_back(&lambdamin);
    optlist.push_back(&lambdamax);
    optlist.push_back(&lambdastart);
    optlist.push_back(&forced_into);
    optlist.push_back(&df_for_lambdamax);
    optlist.push_back(&df_for_lambdamin);
    optlist.push_back(&lambdamax_opt);
    optlist.push_back(&lambdamin_opt);
    optlist.push_back(&number);
    optlist.push_back(&df_equidist);
    optlist.push_back(&df_accuracy);

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
    t.options = vector<ST::string>(14);
    t.options[0] = t.type;
    t.options[1] = map.getvalue();
    t.options[2] = ST::doubletostring(lambda.getvalue());
    t.options[3] = ST::doubletostring(lambdamin.getvalue());
    t.options[4] = ST::doubletostring(lambdamax.getvalue());
    t.options[5] = ST::doubletostring(lambdastart.getvalue());
    if (forced_into.getvalue()==false)
      t.options[6] = "false";
    else
      t.options[6] = "true";
    t.options[7] = ST::doubletostring(df_for_lambdamax.getvalue());
    t.options[8] = ST::doubletostring(df_for_lambdamin.getvalue());
    if (lambdamax_opt.getvalue()==false)
        t.options[9] = "false";
    else
        t.options[9] = "true";
    if (lambdamin_opt.getvalue()==false)
        t.options[10] = "false";
    else
        t.options[10] = "true";
    t.options[11] = ST::inttostring(number.getvalue());
    if (df_equidist.getvalue()==false)
        t.options[12] = "false";
    else
        t.options[12] = "true";
    t.options[13] = ST::doubletostring(df_accuracy.getvalue());


    int b = t.options[3].strtodouble(minl);
    b = t.options[4].strtodouble(maxl);
    b = t.options[5].strtodouble(startl);
    b = t.options[7].strtodouble(df_max);
    b = t.options[8].strtodouble(df_min);

    if (b==1)
      {
      setdefault();
      return false;
      }

    if (minl >= maxl)
      {
      setdefault();
      return false;
      }

    if (startl == -1 && t.options[0] == "spatial")
      {
      setdefault();
      return false;
      }

    if(df_max>=df_min)
      {
      setdefault();
      return false;
      }

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


//------------------------------------------------------------------------------
//--------- class term_randomslope: implementation of member functions ---------
//------------------------------------------------------------------------------

term_randomslope_stepwise::term_randomslope_stepwise(void)
  {
  type = "term_randomslope";
  nofixed = simpleoption("nofixed",false);
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdamin = doubleoption("lambdamin",0.0001,0.000001,10000000);
  lambdamax = doubleoption("lambdamax",10000,0.000001,10000000);
  lambdastart = doubleoption("lambdastart",-1,-1,10000000);
  forced_into = simpleoption("forced_into",false);
  df_for_lambdamax = doubleoption("df_for_lambdamax",2,1.01,500);
  df_for_lambdamin = doubleoption("df_for_lambdamin",10,1.01,500);
  lambdamax_opt = simpleoption("lambdamax_opt",false);
  lambdamin_opt = simpleoption("lambdamin_opt",false);
  number = intoption("number",0,0,50);
  df_equidist = simpleoption("df_equidist",false);
  df_accuracy = doubleoption("df_accuracy",0.05,0.01,0.5);
  }

void term_randomslope_stepwise::setdefault(void)
  {
  nofixed.setdefault();
  lambda.setdefault();
  lambdamin.setdefault();
  lambdamax.setdefault();
  lambdastart.setdefault();
  forced_into.setdefault();
  df_for_lambdamax.setdefault();
  df_for_lambdamin.setdefault();
  lambdamax_opt.setdefault();
  lambdamin_opt.setdefault();
  number.setdefault();
  df_equidist.setdefault();
  df_accuracy.setdefault();
  }


bool term_randomslope_stepwise::check(term & t)
  {

  if ( (t.varnames.size()==2)  && (t.options.size()<=14) )
    {

    if (t.options[0] == "random")
      t.type = "randomslope";
    else
      {
      setdefault();
      return false;
      }

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&nofixed);
    optlist.push_back(&lambda);
    optlist.push_back(&lambdamin);
    optlist.push_back(&lambdamax);
    optlist.push_back(&lambdastart);
    optlist.push_back(&forced_into);
    optlist.push_back(&df_for_lambdamax);
    optlist.push_back(&df_for_lambdamin);
    optlist.push_back(&lambdamax_opt);
    optlist.push_back(&lambdamin_opt);
    optlist.push_back(&number);
    optlist.push_back(&df_equidist);
    optlist.push_back(&df_accuracy);

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
    t.options = vector<ST::string>(14);
    t.options[0] = t.type;
    if (nofixed.getvalue() == true)
      t.options[1] = "true";
    else
      t.options[1] = "false";
    t.options[2] = ST::doubletostring(lambda.getvalue());
    t.options[3] = ST::doubletostring(lambdamin.getvalue());
    t.options[4] = ST::doubletostring(lambdamax.getvalue());
    t.options[5] = ST::doubletostring(lambdastart.getvalue());
    if (forced_into.getvalue()==false)
      t.options[6] = "false";
    else
      t.options[6] = "true";
    t.options[7] = ST::doubletostring(df_for_lambdamax.getvalue());
    t.options[8] = ST::doubletostring(df_for_lambdamin.getvalue());
    if (lambdamax_opt.getvalue()==false)
        t.options[9] = "false";
    else
        t.options[9] = "true";
    if (lambdamin_opt.getvalue()==false)
        t.options[10] = "false";
    else
        t.options[10] = "true";
    t.options[11] = ST::inttostring(number.getvalue());
    if (df_equidist.getvalue()==false)
        t.options[12] = "false";
    else
        t.options[12] = "true";
    t.options[13] = ST::doubletostring(df_accuracy.getvalue());


    double minl, maxl, startl,df_max,df_min;

    int b = t.options[3].strtodouble(minl);
    b = t.options[4].strtodouble(maxl);
    b = t.options[5].strtodouble(startl);
    b = t.options[7].strtodouble(df_max);
    b = t.options[8].strtodouble(df_min);

    if (b==1)
      {
      setdefault();
      return false;
      }


    if (minl >= maxl)
      {
      setdefault();
      return false;
      }

    if (maxl < startl)
      {
      setdefault();
      return false;
      }

    if(df_max>=df_min)
      {
      setdefault();
      return false;
      }

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


//------------------------------------------------------------------------------
//----------- class term_random: implementation of member functions ------------
//------------------------------------------------------------------------------

term_random_stepwise::term_random_stepwise(void)
  {
  type = "term_random";
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdamin = doubleoption("lambdamin",0.0001,0.000001,10000000);
  lambdamax = doubleoption("lambdamax",10000,0.000001,10000000);
  lambdastart = doubleoption("lambdastart",10000,0,10000000);
  forced_into = simpleoption("forced_into",false);
  df_for_lambdamax = doubleoption("df_for_lambdamax",1,0,500);         // Unterscheidung rw1/rw2 bei default-Wert!!!
  df_for_lambdamin = doubleoption("df_for_lambdamin",10,0,500);
  lambdamax_opt = simpleoption("lambdamax_opt",false);
  lambdamin_opt = simpleoption("lambdamin_opt",false);
  number = intoption("number",0,0,50);
  df_equidist = simpleoption("df_equidist",false);
  df_accuracy = doubleoption("df_accuracy",0.05,0.01,0.5);
  }


void term_random_stepwise::setdefault(void)
  {
  lambda.setdefault();
  lambdamin.setdefault();
  lambdamax.setdefault();
  lambdastart.setdefault();
  forced_into.setdefault();
  df_for_lambdamax.setdefault();
  df_for_lambdamin.setdefault();
  lambdamax_opt.setdefault();
  lambdamin_opt.setdefault();
  number.setdefault();
  df_equidist.setdefault();
  df_accuracy.setdefault();
  }


bool term_random_stepwise::check(term & t)
  {

  if ( (t.varnames.size()==1)  && (t.options.size()<=13) )
    {


    if (t.options[0] == "random")
      t.type = "random";

    else
      {
      setdefault();
      return false;
      }

    double minl, maxl,df_max,df_min;

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&lambda);
    optlist.push_back(&lambdamin);
    optlist.push_back(&lambdamax);
    optlist.push_back(&lambdastart);
    optlist.push_back(&forced_into);
    optlist.push_back(&df_for_lambdamax);
    optlist.push_back(&df_for_lambdamin);
    optlist.push_back(&lambdamax_opt);
    optlist.push_back(&lambdamin_opt);
    optlist.push_back(&number);
    optlist.push_back(&df_equidist);
    optlist.push_back(&df_accuracy);

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
    t.options = vector<ST::string>(13);
    t.options[0] = t.type;
    t.options[1] = ST::doubletostring(lambda.getvalue());
    t.options[2] = ST::doubletostring(lambdamin.getvalue());
    t.options[3] = ST::doubletostring(lambdamax.getvalue());
    t.options[4] = ST::doubletostring(lambdastart.getvalue());
    if (forced_into.getvalue()==false)
      t.options[5] = "false";
    else
      t.options[5] = "true";
    t.options[6] = ST::doubletostring(df_for_lambdamax.getvalue());
    t.options[7] = ST::doubletostring(df_for_lambdamin.getvalue());
    if (lambdamax_opt.getvalue()==false)
        t.options[8] = "false";
    else
        t.options[8] = "true";
    if (lambdamin_opt.getvalue()==false)
        t.options[9] = "false";
    else
        t.options[9] = "true";
    t.options[10] = ST::inttostring(number.getvalue());
    if (df_equidist.getvalue()==false)
        t.options[11] = "false";
    else
        t.options[11] = "true";
    t.options[12] = ST::doubletostring(df_accuracy.getvalue());


    int b = t.options[2].strtodouble(minl);
    b = t.options[3].strtodouble(maxl);
    b = t.options[6].strtodouble(df_max);
    b = t.options[7].strtodouble(df_min);

    if (b==1)
      {
      setdefault();
      return false;
      }


    if (minl >= maxl)
      {
      setdefault();
      return false;
      }

    if(df_max>=df_min)
      {
      setdefault();
      return false;
      }

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


//------------------------------------------------------------------------------
//---------- class term_factor: implementation of member functions ------------
//------------------------------------------------------------------------------

term_factor_stepwise::term_factor_stepwise(void)
  {
  type = "term_factor";
  vector<ST::string> code;
  code.push_back("dummy");
  code.push_back("effect");
  coding = stroption("coding",code,"dummy");
  reference = doubleoption("reference",1,-100,100);
  lambdastart = intoption("lambdastart",-1,-1,0);
  forced_into = simpleoption("forced_into",false);
  }

void term_factor_stepwise::setdefault(void)
  {
  coding.setdefault();
  reference.setdefault();
  lambdastart.setdefault();
  forced_into.setdefault();
  }


bool term_factor_stepwise::check(term & t)
  {

  if ( (t.varnames.size() == 1) &&
       (t.options.size()<=5) && (t.options.size() >= 1) )
    {

    if ( (t.options[0] == "factor") && (t.varnames.size() == 1) )
      t.type = "factor";
    else
      {
      setdefault();
      return false;
      }

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&coding);
    optlist.push_back(&reference);
    optlist.push_back(&lambdastart);
    optlist.push_back(&forced_into);

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

    optlist.push_back(&coding);
    optlist.push_back(&reference);
    optlist.push_back(&lambdastart);
    optlist.push_back(&forced_into);

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(5);
    t.options[0] = t.type;
    t.options[1] = coding.getvalue();
    t.options[2] = ST::doubletostring(reference.getvalue());
    t.options[3] = ST::inttostring(lambdastart.getvalue());
    if (forced_into.getvalue()==false)
      t.options[4] = "false";
    else
      t.options[4] = "true";


    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }

//------------------------------------------------------------------------------
//------ class term_interactpspline: implementation of member functions --------
//------------------------------------------------------------------------------

term_interactpspline_stepwise::term_interactpspline_stepwise(void)
  {
  type = "term_interactpspline";
  degree=intoption("degree",3,0,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  gridsize = intoption("gridsize",-1,10,35);
  lambdamin = doubleoption("lambdamin",0.000001,0.000001,100000000);
  lambdamax = doubleoption("lambdamax",10,0.000001,100000000);
  lambdastart = doubleoption("lambdastart",-1,-1,100000000);
  forced_into = simpleoption("forced_into",false);
  df_for_lambdamax = doubleoption("df_for_lambdamax",2,0,400);         // Unterscheidung rw1/rw2 bei default-Wert!!!
  df_for_lambdamin = doubleoption("df_for_lambdamin",10,0,400);
  lambdamax_opt = simpleoption("lambdamax_opt",false);
  lambdamin_opt = simpleoption("lambdamin_opt",false);
  number = intoption("number",0,-1,50);
  df_equidist = simpleoption("df_equidist",false);
  df_accuracy = doubleoption("df_accuracy",0.05,0.01,1);
  }


void term_interactpspline_stepwise::setdefault(void)
  {
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  gridsize.setdefault();
  lambdamin.setdefault();
  lambdamax.setdefault();
  lambdastart.setdefault();
  forced_into.setdefault();
  df_for_lambdamax.setdefault();
  df_for_lambdamin.setdefault();
  lambdamax_opt.setdefault();
  lambdamin_opt.setdefault();
  number.setdefault();
  df_equidist.setdefault();
  df_accuracy.setdefault();
  }

bool term_interactpspline_stepwise::check(term & t)
  {

  optionlist optlist;
  optlist.push_back(&degree);
  optlist.push_back(&numberknots);
  optlist.push_back(&lambda);
  optlist.push_back(&gridsize);
  optlist.push_back(&lambdamin);
  optlist.push_back(&lambdamax);
  optlist.push_back(&lambdastart);
  optlist.push_back(&forced_into);
  optlist.push_back(&df_for_lambdamax);
  optlist.push_back(&df_for_lambdamin);
  optlist.push_back(&lambdamax_opt);
  optlist.push_back(&lambdamin_opt);
  optlist.push_back(&number);
  optlist.push_back(&df_equidist);
  optlist.push_back(&df_accuracy);

  if ( (t.varnames.size()==2)  && (t.options.size() >= 1)
        && (t.options.size() <= 16) )
    {

    if (t.options[0] == "pspline2dimrw1")
      t.type = "pspline2dimrw1";
    else if (t.options[0] == "pspline2dimrw2")
      t.type = "pspline2dimrw2";
    //else if (t.options[0] == "tpspline2dimrw1")
    //  t.type = "tpspline2dimrw1";
    //else if (t.options[0] == "pspline2dimband")
    //  t.type = "pspline2dimband";
    //else if (t.options[0] == "tpspline2dimband")
    //  t.type = "tpspline2dimband";
    //else if (t.options[0] == "psplinekrrw1")
    //  t.type = "psplinekrrw1";
    //else if (t.options[0] == "psplinekrrw2")
    //  t.type = "psplinekrrw2";
    else
      {
      setdefault();
      return false;
      }

    double maxl,minl,startl,df_max,df_min;

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
   t.options = vector<ST::string>(16);
   t.options[0] = t.type;
   t.options[1] = ST::inttostring(degree.getvalue());
   t.options[2] = ST::inttostring(numberknots.getvalue());
   t.options[3] = ST::doubletostring(lambda.getvalue());
   t.options[4] = ST::inttostring(gridsize.getvalue());
   t.options[5] = ST::doubletostring(lambdamin.getvalue());
   t.options[6] = ST::doubletostring(lambdamax.getvalue());
   t.options[7] = ST::doubletostring(lambdastart.getvalue());
   if (forced_into.getvalue()==false)
     t.options[8] = "false";
   else
     t.options[8] = "true";
   t.options[9] = ST::doubletostring(df_for_lambdamax.getvalue());
   t.options[10] = ST::doubletostring(df_for_lambdamin.getvalue());
   if (lambdamax_opt.getvalue()==false)
       t.options[11] = "false";
   else
       t.options[11] = "true";
   if (lambdamin_opt.getvalue()==false)
       t.options[12] = "false";
   else
       t.options[12] = "true";
   t.options[13] = ST::inttostring(number.getvalue());
   if (df_equidist.getvalue()==false)
      t.options[14] = "false";
   else
      t.options[14] = "true";
   t.options[15] = ST::doubletostring(df_accuracy.getvalue());


   if (lambda.getvalue() < 0)
     {
     setdefault();
     return false;
     }

    int b = t.options[5].strtodouble(minl);
    b = t.options[6].strtodouble(maxl);
    b = t.options[7].strtodouble(startl);
    b = t.options[9].strtodouble(df_max);
    b = t.options[10].strtodouble(df_min);

    if (b==1)
      {
      setdefault();
      return false;
      }


    if (minl >= maxl)
      {
      setdefault();
      return false;
      }

    if (maxl < startl)
      {
      setdefault();
      return false;
      }

    if ((df_min == 1 && t.options[0] == "pspline2dimrw1") ||            // kann man auf Anzahl der Parameter zugreifen???
            (df_max == 1 && t.options[0] == "pspline2dimrw1"))
      {
      setdefault();
      return false;
      }

    if(df_max>=df_min)
      {
      setdefault();
      return false;
      }

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


bool term_interactpspline_stepwise::checkvector(const vector<term> & terms,
                                       const unsigned & i)
  {

  assert(i< terms.size());

  if ((terms[i].type == "pspline2dimrw1") || (terms[i].type == "pspline2dimrw2"))
     //|| (terms[i].type == "psplinekrrw1") || (terms[i].type == "psplinekrrw2")
     //|| (terms[i].type == "tpspline2dimrw1")
     //|| (terms[i].type == "pspline2dimband") || (terms[i].type == "tpspline2dimband") )
    return true;

  return false;
  }


//------------------------------------------------------------------------------
//---------- class term_geospline: implementation of member functions ----------
//------------------------------------------------------------------------------

term_geospline_stepwise::term_geospline_stepwise(void)
  {
  type = "term_geospline";
  map=stroption("map");
  degree=intoption("degree",3,1,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdamin = doubleoption("lambdamin",0.0001,0.000001,10000000);
  lambdamax = doubleoption("lambdamax",10000,0.000001,10000000);
  lambdastart = doubleoption("lambdastart",10000,-1,10000000);
  forced_into = simpleoption("forced_into",false);
  df_for_lambdamax = doubleoption("df_for_lambdamax",1,0,500);         // Unterscheidung rw1/rw2 bei default-Wert!!!
  df_for_lambdamin = doubleoption("df_for_lambdamin",10,0,500);
  lambdamax_opt = simpleoption("lambdamax_opt",false);
  lambdamin_opt = simpleoption("lambdamin_opt",false);
  number = intoption("number",0,0,50);
  df_equidist = simpleoption("df_equidist",false);
  df_accuracy = doubleoption("df_accuracy",0.05,0.01,0.5);
  }


void term_geospline_stepwise::setdefault(void)
  {
  map.setdefault();
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  lambdamin.setdefault();
  lambdamax.setdefault();
  lambdastart.setdefault();
  forced_into.setdefault();
  df_for_lambdamax.setdefault();
  df_for_lambdamin.setdefault();
  lambdamax_opt.setdefault();
  lambdamin_opt.setdefault();
  number.setdefault();
  df_equidist.setdefault();
  df_accuracy.setdefault();
  }

bool term_geospline_stepwise::check(term & t)
  {

  if ( (t.varnames.size()==1)  && (t.options.size() >= 1)
        && (t.options.size() <= 16) )
    {

    if (t.options[0] == "geospline")
      t.type = "geospline";
    else if (t.options[0] == "geosplinerw1")
      t.type = "geospline";
    else if (t.options[0] == "geosplinerw2")
      t.type = "geosplinerw2";
    else
      {
      setdefault();
      return false;
      }

    double maxl,minl,startl,df_max,df_min;

    optionlist optlist;
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&lambda);
    optlist.push_back(&map);
    optlist.push_back(&lambdamin);
    optlist.push_back(&lambdamax);
    optlist.push_back(&lambdastart);
    optlist.push_back(&forced_into);
    optlist.push_back(&df_for_lambdamax);
    optlist.push_back(&df_for_lambdamin);
    optlist.push_back(&lambdamax_opt);
    optlist.push_back(&lambdamin_opt);
    optlist.push_back(&number);
    optlist.push_back(&df_equidist);
    optlist.push_back(&df_accuracy);

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
    t.options = vector<ST::string>(16);
    t.options[0] = t.type;
    t.options[1] = ST::inttostring(degree.getvalue());
    t.options[2] = ST::inttostring(numberknots.getvalue());
    t.options[3] = ST::doubletostring(lambda.getvalue());
    t.options[4] = map.getvalue();
    t.options[5] = ST::doubletostring(lambdamin.getvalue());
    t.options[6] = ST::doubletostring(lambdamax.getvalue());
    t.options[7] = ST::doubletostring(lambdastart.getvalue());
    if (forced_into.getvalue()==false)
      t.options[8] = "false";
    else
      t.options[8] = "true";
    t.options[9] = ST::doubletostring(df_for_lambdamax.getvalue());
    t.options[10] = ST::doubletostring(df_for_lambdamin.getvalue());
    if (lambdamax_opt.getvalue()==false)
        t.options[11] = "false";
    else
        t.options[11] = "true";
    if (lambdamin_opt.getvalue()==false)
        t.options[12] = "false";
    else
        t.options[12] = "true";
    t.options[13] = ST::inttostring(number.getvalue());
    if (df_equidist.getvalue()==false)
       t.options[14] = "false";
    else
       t.options[14] = "true";
    t.options[15] = ST::doubletostring(df_accuracy.getvalue());


    int b = t.options[5].strtodouble(minl);
    b = t.options[6].strtodouble(maxl);
    b = t.options[7].strtodouble(startl);
    b = t.options[9].strtodouble(df_max);
    b = t.options[10].strtodouble(df_min);

    if (b==1)
      {
      setdefault();
      return false;
      }


    if (minl >= maxl)
      {
      setdefault();
      return false;
      }

    if (maxl < startl)
      {
      setdefault();
      return false;
      }

    if(df_max>=df_min)
      {
      setdefault();
      return false;
      }

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


bool term_geospline_stepwise::checkvector(const vector<term> & terms, const unsigned & i)
  {

  assert(i< terms.size());

  if (terms[i].type == "geospline")
    return true;

  return false;
  }


