//---------------------------------------------------------------------------
#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include<StatwinFrame.h>
#include<statwin_haupt.h>

#endif


#include<remlreg.h>
#include<fullcond.h>
#include<typeinfo.h>
#include<stddef.h>


//------------------------------------------------------------------------------
//------------- CLASS remlreg: implementation of member functions -------------
//------------------------------------------------------------------------------

void remlreg::make_paths(unsigned  collinpred,ST::string & pathnonp,
                          ST::string & pathres,ST::string & title,
                          ST::string varname1,ST::string varname2,
                          ST::string endingraw,ST::string endingres,
                          ST::string endingtitle)
  {
  if (collinpred == 0)
    {
    if (varname2=="")
      {
      pathnonp = defaultpath + "\\temp\\" + name +  add_name + "_f_" +
                 varname1 + endingraw;

      pathres = outfile.getvalue() + add_name + "_f_" + varname1 + endingres;

      title = "f_"+ varname1 + endingtitle +  add_name;
      }
    else
      {
      pathnonp = defaultpath + "\\temp\\" + name +  add_name + "_" +
                 varname2 +  "_f_" + varname1 + endingraw;

      pathres = outfile.getvalue() + add_name + "_" + varname2 +
                "_f_" + varname1 + endingres;

      title = varname2 + "_f_"+ varname1 + endingtitle +  add_name;
      }
    }
  else
    {
    if (varname2=="")
      {
      pathnonp = defaultpath + "\\temp\\" + name +  add_name + "_f_" +
                 ST::inttostring(collinpred+1) + "_" +
                     varname1 + endingraw;

      pathres = outfile.getvalue() + add_name + "_f_" +
                ST::inttostring(collinpred+1) + "_" +
                varname1 + endingres;

      title = "f_"+ ST::inttostring(collinpred+1) + "_" +
               varname1 + endingtitle + add_name;
      }
    else
      {
      pathnonp = defaultpath + "\\temp\\" + name +  add_name + "_" + varname2
                 + "_f_" + ST::inttostring(collinpred+1) + "_" +
                     varname1 + endingraw;

      pathres = outfile.getvalue() + add_name + "_" + varname2 + "_f_" +
                ST::inttostring(collinpred+1) + "_" +
                varname1 + endingres;

      title = varname2 + "_f_"+ ST::inttostring(collinpred+1) + "_" +
               varname1 + endingtitle + add_name;
      }
    }
  }

void remlreg::initpointers(void)
  {
  unsigned i;

  for(i=0;i<fcpsplinesurf.size();i++)
    fullcond.push_back(&fcpsplinesurf[i]);
  for(i=0;i<fcconst.size();i++)
    fullcond.push_back(&fcconst[i]);
  for(i=0;i<fcnonpgaussian.size();i++)
    fullcond.push_back(&fcnonpgaussian[i]);
  for(i=0;i<fcpspline.size();i++)
    fullcond.push_back(&fcpspline[i]);
  for(i=0;i<fcrandom.size();i++)
    fullcond.push_back(&fcrandom[i]);
  for(i=0; i<fckriging.size();i++)
    fullcond.push_back(&fckriging[i]);
  for(i=0; i<fcbaseline.size(); i++)
    fullcond.push_back(&fcbaseline[i]);
  for(i=0; i<fcbaseline_varcoeff.size(); i++)
    fullcond.push_back(&fcbaseline_varcoeff[i]);
  }

void remlreg::clear(void)
  {
  fullcond.erase(fullcond.begin(),fullcond.end());
  fullcond.reserve(50);
  fcconst.erase(fcconst.begin(),fcconst.end());
  fcconst.reserve(20);
  fcnonpgaussian.erase(fcnonpgaussian.begin(),fcnonpgaussian.end());
  fcnonpgaussian.reserve(20);
  fcpspline.erase(fcpspline.begin(),fcpspline.end());
  fcpspline.reserve(20);
  fcpsplinesurf.erase(fcpsplinesurf.begin(),fcpsplinesurf.end());
  fcpsplinesurf.reserve(20);
  fcrandom.erase(fcrandom.begin(),fcrandom.end());
  fcrandom.reserve(20);
  fckriging.erase(fckriging.begin(),fckriging.end());
  fckriging.reserve(20);
  fcbaseline.erase(fcbaseline.begin(),fcbaseline.end());
  fcbaseline.reserve(20);
  fcbaseline_varcoeff.erase(fcbaseline_varcoeff.begin(),fcbaseline_varcoeff.end());
  fcbaseline_varcoeff.reserve(20);
  }

void remlreg::create(void)
  {
  // SYNTAX OF COMMANDS:
  // name [model] [weight varname] [by varname] [if expression]
  //      [, options] [using usingtext]

  add_name="";
  ST::string h = defaultpath+"\\output\\"+name;

  outfile = fileoption("outfile",h,false);

  globaloptions.push_back(&outfile);

//------------------------------------------------------------------------------
// ----------------------------- method regress --------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------- for the offset -----------------------------------
//------------------------------------------------------------------------------

  offset = term_offset();
  termtypes.push_back(&offset);

//------------------------------------------------------------------------------
//-------------------------  for fixed effects ---------------------------------
//------------------------------------------------------------------------------

  fixedeffects = basic_termtype();
  termtypes.push_back(&fixedeffects);

//------------------------------------------------------------------------------
//------------------------ for nonparametric terms -----------------------------
//------------------------------------------------------------------------------

  nonprw1rw2 = term_autoreg_remlreg();
  termtypes.push_back(&nonprw1rw2);

  nonprw1rw2_varcoef = term_autoreg_varcoef_remlreg();
  termtypes.push_back(&nonprw1rw2_varcoef);

  nonpseason = term_season_remlreg();
  termtypes.push_back(&nonpseason);

  nonpseason_varcoef = term_season_varcoef_remlreg();
  termtypes.push_back(&nonpseason_varcoef);

  nonpspatial = term_spatial_remlreg();
  termtypes.push_back(&nonpspatial);

  nonpspatial_varcoef = term_spatial_varcoef_remlreg();
  termtypes.push_back(&nonpspatial_varcoef);

//------------------------------------------------------------------------------
//-------------------------- for p-spline terms --------------------------------
//------------------------------------------------------------------------------

  nonppspline = term_pspline_remlreg();
  termtypes.push_back(&nonppspline);

  nonpvarcoeffpspline = term_varcoeff_pspline_remlreg();
  termtypes.push_back(&nonpvarcoeffpspline);

  nonpinteractpspline = term_interactpspline_remlreg();
  termtypes.push_back(&nonpinteractpspline);

  nonpvarcoeffinteractpspline = term_interactpspline_varcoeff_remlreg();
  termtypes.push_back(&nonpvarcoeffinteractpspline);

  nonpgeospline = term_geospline_remlreg();
  termtypes.push_back(&nonpgeospline);

  nonpvarcoeffgeospline = term_geospline_varcoeff_remlreg();
  termtypes.push_back(&nonpvarcoeffgeospline);

//------------------------------------------------------------------------------
//-------------------------- for kriging terms ---------------------------------
//------------------------------------------------------------------------------

  nonpspatial_kriging = term_kriging_remlreg();
  termtypes.push_back(&nonpspatial_kriging);

  nonp_kriging = term_kriging_1dim_remlreg();
  termtypes.push_back(&nonp_kriging);

  nonpspatial_kriging_varcoeff = term_kriging_varcoeff_remlreg();
  termtypes.push_back(&nonpspatial_kriging_varcoeff);

  nonpspatial_geokriging = term_geokriging_remlreg();
  termtypes.push_back(&nonpspatial_geokriging);

  nonpspatial_geokriging_varcoeff = term_geokriging_varcoeff_remlreg();
  termtypes.push_back(&nonpspatial_geokriging_varcoeff);

//------------------------------------------------------------------------------
//-------------------------- for baseline terms --------------------------------
//------------------------------------------------------------------------------

  nonp_baseline = term_baseline_remlreg();
  termtypes.push_back(&nonp_baseline);

  nonp_baseline_varcoeff = term_baseline_varcoeff_remlreg();
  termtypes.push_back(&nonp_baseline_varcoeff);

//------------------------------------------------------------------------------
//------------------------ for random effect terms -----------------------------
//------------------------------------------------------------------------------

  randomeff = term_random_remlreg();
  termtypes.push_back(&randomeff);

  randomeffslope = term_randomslope_remlreg();
  termtypes.push_back(&randomeffslope);

//------------------------------------------------------------------------------
//--------------------------- General options ----------------------------------
//------------------------------------------------------------------------------

  modreg = modelterm(&termtypes);

  udata = use();

  knotsdef.push_back("equidistant");
  knotsdef.push_back("quantiles");
  knots = stroption("knots",knotsdef,"equidistant");

  level1 = doubleoption("level1",95,40,99);
  level2 = doubleoption("level2",80,40,99);
  maxint = intoption("maxint",150,0,20000);

  families.reserve(15);
  families.push_back("gaussian");
  families.push_back("binomial");
  families.push_back("binomialprobit");
  families.push_back("binomialcomploglog");
  families.push_back("poisson");
  families.push_back("gamma");
  families.push_back("poissondispers");
  families.push_back("binomialdispers");
  families.push_back("binomialprobitdispers");
  families.push_back("multinomial");
  families.push_back("cumlogit");
  families.push_back("cumprobit");
  families.push_back("seqlogit");
  families.push_back("seqprobit");
  families.push_back("cox");
  families.push_back("coxinterval");
  families.push_back("coxinterval2");
  family = stroption("family",families,"binomial");

  maxit = intoption("maxit",400,1,100000);
  lowerlim = doubleoption("lowerlim",0.001,0,1);
  eps = doubleoption("eps",0.00001,0,1);

  reference = doubleoption("reference",0,-10000,10000);

  noconst = simpleoption("noconst",false);

  leftint = stroption("leftint");
  lefttrunc = stroption("lefttrunc");

  regressoptions.reserve(100);

  regressoptions.push_back(&level1);
  regressoptions.push_back(&level2);
  regressoptions.push_back(&maxint);
  regressoptions.push_back(&family);

  regressoptions.push_back(&knots);
  regressoptions.push_back(&maxit);
  regressoptions.push_back(&lowerlim);
  regressoptions.push_back(&eps);
  regressoptions.push_back(&reference);

  regressoptions.push_back(&noconst);

  regressoptions.push_back(&leftint);
  regressoptions.push_back(&lefttrunc);

//------------------------------------------------------------------------------
//-------------------------- methods[0]: remlrun -------------------------------
//------------------------------------------------------------------------------

  methods.push_back(command("regress",&modreg,&regressoptions,&udata,required,
			 optional,optional,optional,optional,required));

  functions[0] = remlrun;

//------------------------------------------------------------------------------
// ----------------------------- method plotnonp -------------------------------
//------------------------------------------------------------------------------

  uplotnonp = use();
  mplotnonp = modelStandard();

  xlab = stroption("xlab");
  ylab = stroption("ylab");
  connect = stroption("connect");
  height = intoption("height",210,0,500);
  width = intoption("width",356,0,500);
  ylimtop = doubleoption("ylimtop",-MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  ylimbottom = doubleoption("ylimbottom",MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  xlimtop = doubleoption("xlimtop",-MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  xlimbottom = doubleoption("xlimbottom",MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  xstep = doubleoption("xstep",0.0,-MAXDOUBLE,MAXDOUBLE);
  ystep = doubleoption("ystep",0.0,-MAXDOUBLE,MAXDOUBLE);
  xstart = doubleoption("xstart",MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  ystart = doubleoption("ystart",MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  linewidth = intoption("linewidth",5,0,100);
  fontsize = intoption("fontsize",12,0,100);
  pointsize = intoption("pointsize",20,0,100);
  linecolor = stroption("linecolor");
  titlescale = doubleoption("titlesize",1.5,0.0,MAXDOUBLE);

  vector<ST::string> levelchoice;
  levelchoice.reserve(4);
  levelchoice.push_back("all");
  levelchoice.push_back("1");
  levelchoice.push_back("2");
  levelchoice.push_back("none");
  outfile2=stroption("outfile");
  title = stroption("title");
  replace2 = simpleoption("replace",false);

  levels = stroption("levels",levelchoice,"all");
  median = simpleoption("median",false);

  plotnonpoptions.push_back(&xlab);
  plotnonpoptions.push_back(&ylab);
  plotnonpoptions.push_back(&connect);
  plotnonpoptions.push_back(&height);
  plotnonpoptions.push_back(&width);
  plotnonpoptions.push_back(&ylimtop);
  plotnonpoptions.push_back(&ylimbottom);
  plotnonpoptions.push_back(&ystep);
  plotnonpoptions.push_back(&ystart);
  plotnonpoptions.push_back(&xlimtop);
  plotnonpoptions.push_back(&xlimbottom);
  plotnonpoptions.push_back(&xstep);
  plotnonpoptions.push_back(&xstart);
  plotnonpoptions.push_back(&levels);
  plotnonpoptions.push_back(&median);
  plotnonpoptions.push_back(&title);
  plotnonpoptions.push_back(&outfile2);
  plotnonpoptions.push_back(&replace2);
  plotnonpoptions.push_back(&linewidth);
  plotnonpoptions.push_back(&fontsize);
  plotnonpoptions.push_back(&pointsize);
  plotnonpoptions.push_back(&linecolor);
  plotnonpoptions.push_back(&titlescale);

  // methods[1]:

  methods.push_back(command("plotnonp",&mplotnonp,&plotnonpoptions,&uplotnonp,
                   required,notallowed,notallowed,notallowed,optional,
                   notallowed));

  functions[1] = plotnonprun;

//------------------------------------------------------------------------------
// ---------------------------- method drawmaprun ------------------------------
//------------------------------------------------------------------------------

  udrawmap = use();

  mdrawmap = modelStandard();

  outfile4 = stroption("outfile");
  title2 = stroption("title");
  upperlimit = doubleoption("upperlimit",1,-MAXDOUBLE,MAXDOUBLE);
  lowerlimit = doubleoption("lowerlimit",0,-MAXDOUBLE,MAXDOUBLE);
  nrcolors = intoption("nrcolors",256,1,256);
  color = simpleoption("color",false);
  nolegend = simpleoption("nolegend",false);
  swapcolors = simpleoption("swapcolors",false);
  replace = simpleoption("replace",false);
  plotvar = stroption("plotvar","pmode");
  pcat = simpleoption("pcat",false);
  drawnames = simpleoption("drawnames",false);

  drawmapoptions.push_back(&outfile4);
  drawmapoptions.push_back(&title2);
  drawmapoptions.push_back(&upperlimit);
  drawmapoptions.push_back(&lowerlimit);
  drawmapoptions.push_back(&nrcolors);
  drawmapoptions.push_back(&color);
  drawmapoptions.push_back(&nolegend);
  drawmapoptions.push_back(&swapcolors);
  drawmapoptions.push_back(&replace);
  drawmapoptions.push_back(&plotvar);
  drawmapoptions.push_back(&pcat);
  drawmapoptions.push_back(&drawnames);

  // methods[2]:
  methods.push_back(command("drawmap",&mdrawmap,&drawmapoptions,&udrawmap,
                   required,notallowed,notallowed,notallowed,optional,
                   notallowed));

  functions[2] = drawmaprun;

  }

//------------------------------------------------------------------------------
// ------------------------------- Constructor ---------------------------------
//------------------------------------------------------------------------------

  // CONSTRUCTOR 1:
  // ADDITIONAL INFORMATION: name = n

remlreg::remlreg(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adb,
  #endif
  const ST::string & n,ofstream * lo,istream * in,
						 ST::string p,vector<statobject*> * st)
						 : statobject(
                         #if defined(JAVA_OUTPUT_WINDOW)
                         adb,
                         #endif
                         n,"remlreg",lo,in,p)
  {
  statobj = st;
  create();
  resultsyesno = false;
  describetext.push_back("CURRENT REGRESSION RESULTS: none\n");
  }

  // COPY CONSTRUCTOR

remlreg::remlreg(const remlreg & b) : statobject(statobject(b))
  {
  create();
  statobj = b.statobj;
  D = b.D;
  modelvarnamesv = b.modelvarnamesv;
  terms = b.terms;
  fcpspline = b.fcpspline;
  resultsyesno = b.resultsyesno;
  initpointers();
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const remlreg & remlreg::operator=(const remlreg & b)
  {
  if (this == & b)
	 return *this;
  statobject::operator=(statobject(b));
  create();
  statobj = b.statobj;
  D = b.D;
  modelvarnamesv = b.modelvarnamesv;
  terms = b.terms;
  fcpspline = b.fcpspline;
  resultsyesno = b.resultsyesno;
  initpointers();
  return *this;
  }

int remlreg::parse(const ST::string & c)
  {
  int u = statobject::parse(c);
  int pos = statobject::parsecom(c,methods,globaloptions);

  if (pos >= 0)
	 (*functions[pos])(*this);
  }

void remlreg::describe(optionlist & globaloptions)
  {
  statobject::describe(globaloptions);
  }

// -----------------------------------------------------------------------------
// ---------------------------- Create data ------------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_data(datamatrix & weight)
  {
  unsigned i;
  bool failure=false;

  // find the dataset object

  dataobject * datap;               // pointer to dataset
  int objpos = findstatobject(*statobj,udata.getusingtext(),"dataset");
  statobject * s;
  if (objpos >= 0)
    {
    s = statobj->at(objpos);
    datap = dynamic_cast<dataobject*>(s);
    }
  else
    {
    if (objpos == -1)
      outerror("ERROR: " + udata.getusingtext() + " is not existing\n");
    else
      outerror("ERROR: " + udata.getusingtext() + " is not a dataset object\n");
    return true;
    }

  // create data matrix and waits

  ST::string rname;
  ST::string wn;
  ST::string ifexpression;
  ST::string predictindicator;
  ST::string missingindicator;
  int wpos=-1;

  modelvarnamesv = modreg.getModelVarnamesAsVector();
  rname = modelvarnamesv[0].to_bstr();
  wn = methods[0].get_weight_variable().to_bstr();

  if (wn.length() != 0)
    {
    modelvarnamesv.push_back(wn);
    wpos = modelvarnamesv.size()-1;
    }

  ifexpression = methods[0].getexpression();

  // F�r Cox model: variable 'leftint' an 2. Stelle setzen
  vector<ST::string>::iterator modelit = modelvarnamesv.begin()+1;
  if(leftint.getvalue() != "")
    modelvarnamesv.insert(modelit,1,leftint.getvalue());
  // F�r Cox model: variable 'lefttrunc' an 3. Stelle setzen
  modelit++;
  if(lefttrunc.getvalue() != "")
    modelvarnamesv.insert(modelit,2,lefttrunc.getvalue());

  // testing, wether all variables specified are already existing
  vector<ST::string> notex;
  if ((datap->allexisting(modelvarnamesv,notex)) == false)
    {
    for (i=0;i<notex.size();i++)
      if (notex[i] != "const")
        {
        outerror("ERROR: variable " + notex[i] + " is not existing\n");
        failure = true;
        }
    if (failure)
      return true;
    }

  datap->makematrix(modelvarnamesv,D,ifexpression);

  // Check for error messages

  errormessages = datap->geterrormessages();

  if (!errormessages.empty())
    return true;

  // extract weights

  if (wpos==-1)
    {
    // No weights specified
    weight = datamatrix(D.rows(),1,1);
    }
  else
    {
    // check for correct distributions
    if(family.getvalue()=="multinomial" ||
       family.getvalue()=="cumlogit" ||
       family.getvalue()=="cumprobit" ||
       family.getvalue()=="seqlogit" ||
       family.getvalue()=="seqprobit")
      {
      outerror("ERROR: weight not allowed for multicategorical response\n");
      return true;
      }
    if(family.getvalue()=="cox")
      {
      outerror("ERROR: weight not allowed for family=cox\n");
      return true;
      }

    weight = D.getCol(wpos);

    // check for negative weights
    if(weight.min(0)<0)
      {
      outerror("ERROR: negative weights encountered\n");
      return true;
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// -------------------------- Create Response ----------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_response(datamatrix & response, datamatrix & weight)
  {

  // extract response variable
  response = D.getCol(0);

  // family=binomial
  if(family.getvalue()=="binomial" ||
     family.getvalue()=="binomialprobit" ||
     family.getvalue()=="binomialcomploglog" ||
     family.getvalue()=="binomialdispers" ||
     family.getvalue()=="binomialprobitdispers")
    {
    unsigned i;
    for(i=0; i<response.rows(); i++)
      {
      if (response(i,0) != int(response(i,0)))
        {
        outerror("ERROR: response cannot be binomial\; values must be integer numbers\n");
        return true;
        }

      if (response(i,0) < 0)
        {
        outerror("ERROR: response cannot be binomial\; some values are negative\n");
        return true;
        }

      if (response(i,0) > weight(i,0))
        {
        outerror("ERROR: response cannot be binomial\;\n");
        outerror("       number of successes larger than number of trials for some values\n");
        return true;
        }
      // Transform response for binomial families  (to match usual GLM definition)
      response(i,0) = response(i,0)/weight(i,0);
      }
    }

  // family=gamma
  if(family.getvalue()=="gamma")
    {
    if(response.min(0)<0)
      {
      outerror("ERROR: response cannot be gamma distributed\; some values are negative\n");
      return true;
      }
    }

  // family=poisson
  if(family.getvalue()=="poisson")
    {
    if(response.min(0)<0)
      {
      outerror("ERROR: response cannot be poisson distributed\; some values are negative\n");
      return true;
      }
    }

  // family=cox
  if(family.getvalue()=="cox")
    {
    unsigned i;
    for(i=0; i<response.rows(); i++)
      {
      if (response(i,0)!=0 && response(i,0)!=1)
        {
        outerror("ERROR: response must be either zero or one\n");
        return true;
        }
      }

    // check whether there is a baseline effect
    bool baselineexisting = false;
    for(i=0;i<terms.size();i++)
      {
      if(nonp_baseline.checkvector(terms,i) == true )
        {
        baselineexisting = true;
        }
      }
    if(baselineexisting == false)
      {
      outerror("ERROR: no baseline specified\n");
      return true;
      }
    }

  // check whether a baseline effect is specified if family != cox
  if(family.getvalue()!="cox" && family.getvalue()!="coxinterval" &&
     family.getvalue()!="coxinterval2")
    {
    unsigned i;
    bool baselineexisting = false;
    for(i=0;i<terms.size();i++)
      {
      if(nonp_baseline.checkvector(terms,i) == true )
        {
        baselineexisting = true;
        }
      }
    if(baselineexisting == true)
      {
      outerror("ERROR: term baseline can only be used with family=cox\n");
      return true;
      }
    }

  // family=multinomial / family=cumlogit / family=cumprobit
  if (family.getvalue()=="multinomial")
    {
    ismultinomial=true;
    }
  else
    {
    ismultinomial=false;
    }

  if (family.getvalue()=="multinomial" || family.getvalue()=="cumlogit" ||
      family.getvalue()=="cumprobit" || family.getvalue()=="seqlogit" ||
      family.getvalue()=="seqprobit")
    {
    // extract categories
    datamatrix resphelp=response;
    resphelp.sort(0,resphelp.rows()-1,0);
    vector<double> categories;
    categories.push_back(resphelp(0,0));
    unsigned i;
    unsigned refpos;
    for(i=1; i<resphelp.rows(); i++)
      {
      if(resphelp(i,0)!=resphelp(i-1,0))
        {
        categories.push_back(resphelp(i,0));
        }
      }

    // check for proper categories
    if (categories.size() == 1)
      {
      outerror("ERROR: response variable does not vary\n");
      return true;
      }
    if (categories.size() > 12)
      {
      outerror("ERROR: too many values for the response variable\n");
      return true;
      }

    // Define / extract the reference category
    if(family.getvalue()=="multinomial")
      {
      refpos=0; //categories.size()-1;
      if(reference.changed())
        {
        bool existing = false;
        double ref=reference.getvalue();
        i=0;
        while(!existing && i<categories.size())
          {
          if(categories[i] == ref)
            {
            existing=true;
            refpos=i;
            }
          i++;
          }
        if(!existing)
          {
          outerror("ERROR: reference category is not existing\n");
          return true;
          }
        }
      }
    else
      {
      if(reference.changed())
        {
        outerror("ERROR: Option reference is not allowed in cumulative models\n");
        return true;
        }
      else
        {
        refpos=categories.size()-1;
        }
      }
    cats=datamatrix(categories.size()-1,1,0);
    unsigned j=0;
    for(i=0; i<categories.size(); i++)
      {
      if(i!=refpos)
        {
        cats(j,0)=categories[i];
        j++;
        }
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// --------------------------- Create_offset -----------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_offset(datamatrix & o)
  {
  unsigned i;
  for(i=0;i<terms.size();i++)
    {
    if ( offset.checkvector(terms,i) == true)
      {
      // check for right distributions
      if(family.getvalue()=="multinomial" ||
         family.getvalue()=="cumlogit" ||
         family.getvalue()=="cumprobit" ||
         family.getvalue()=="seqlogit" ||
         family.getvalue()=="seqprobit")
        {
        outerror("ERROR: offset not allowed for multicategorical response\n");
        return true;
        }
      if(family.getvalue()=="cox")
        {
        outerror("ERROR: offset not allowed for family=cox\n");
        return true;
        }

      // Extract offset if specified
      unsigned j = terms[i].varnames[0].isinlist(modelvarnamesv);
      if (o.rows() < D.rows())
        o = datamatrix(D.rows(),1,0);

      register unsigned k;
      double * worko = o.getV();
      double * workD = D.getV()+j;
      unsigned size = D.cols();
      for(k=0;k<D.rows();k++,worko++,workD+=size)
        *worko += *workD;
      }
    else
      {
      // Default offset (=0)
      o = datamatrix(D.rows(),1,0);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// --------------------------- Create_const -----------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_const(const unsigned & collinpred)
  {
  unsigned i;
  int j;

  vector<ST::string> varnames;
  vector<ST::string> varnamesh =  fixedeffects.get_constvariables(terms);

  if(noconst.getvalue()==false)
    {
    varnames.push_back("const");
    }

  for(i=0;i<varnamesh.size();i++)
    varnames.push_back(varnamesh[i]);

  unsigned nr = varnames.size();

  ST::string title;
  ST::string pathconst;
  ST::string pathconstres;

  if (collinpred == 0)
    {
    title = "FixedEffects" + add_name;
    pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                           + add_name + "_FixedEffects" + ".raw";

    pathconstres = outfile.getvalue() + add_name +
                     "_FixedEffects"  + ".res";
    }
  else
    {
    title = "FixedEffects"  "_" +
                            ST::inttostring(collinpred+1) + add_name;
    pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                           + add_name + "_FixedEffects"  +
                           "_" + ST::inttostring(collinpred+1) + ".raw";

    pathconstres = outfile.getvalue() + add_name + "_FixedEffects" + "_" +
                   ST::inttostring(collinpred+1) + ".res";
    }

  if (pathconst.isvalidfile() == 1)
    {
    errormessages.push_back("ERROR: unable to open file " + pathconst +
                                 " for writing\n");
    return true;
    }

  datamatrix X(D.rows(),nr,1);

  for(i=0;i<varnames.size();i++)
    {
    if (varnames[i] != "const")
      {
      j = varnames[i].isinlist(modelvarnamesv);

      if (j != -1)
        {
        unsigned l;
        double * workX=X.getV()+i;
        double * workD=D.getV()+j;
        for (l=0;l<X.rows();l++,workX+=X.cols(),workD+=D.cols())
          *workX = *workD;
        }
      }
    }

  fcconst.push_back(FULLCOND_const(&generaloptions,X,title,0,
                                   pathconst,pathconstres));

  fcconst[fcconst.size()-1].init_names(varnames);
  fcconst[fcconst.size()-1].set_fcnumber(fullcond.size());
  fullcond.push_back(&fcconst[fcconst.size()-1]);

  return false;
  }

// -----------------------------------------------------------------------------
// ------------------------- create_nonprw1rw2 ---------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_nonprw1rw2(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string title;

  double hd;
  double lambda, startlambda;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonprw1rw2.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "rw1")
        type = MCMC::RW1;
      else
        type = MCMC::RW2;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtodouble(hd);
      lambda = hd;
      f = (terms[i].options[2]).strtodouble(hd);
      startlambda = hd;

      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_rw.raw","_rw.res","_rw");

      fcnonpgaussian.push_back(FULLCOND_nonp_gaussian(&generaloptions,D.getCol(j),
                         unsigned(maxint.getvalue()),type,
                         title,pathres,lambda,startlambda));

      fcnonpgaussian[fcnonpgaussian.size()-1].init_name(terms[i].varnames[0]);
      fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// --------------------- create_nonprw1rw2_varcoef -----------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_nonprw1rw2_varcoef(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string title;

  double hd;
  double lambda, startlambda;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonprw1rw2_varcoef.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "varcoeffrw1")
        type = MCMC::RW1;
      else
        type = MCMC::RW2;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interacting var
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv); // effectmod

      f = (terms[i].options[1]).strtodouble(hd);
      lambda = hd;
      f = (terms[i].options[2]).strtodouble(hd);
      startlambda = hd;

      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],"_rw.raw","_rw.res","_rw");

      fcnonpgaussian.push_back(FULLCOND_nonp_gaussian(&generaloptions,
                         D.getCol(j2),D.getCol(j1),
                         unsigned(maxint.getvalue()),type,
                         title,pathres,lambda,startlambda));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fcnonpgaussian[fcnonpgaussian.size()-1].init_names(na);
      fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ------------------------- create_nonpseason ---------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_nonpseason(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;

  long h;
  double hd;
  double lambda, startlambda;
  unsigned per;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpseason.checkvector(terms,i) == true )
      {
      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      per = unsigned(h);
      f = (terms[i].options[2]).strtodouble(hd);
      lambda = hd;
      f = (terms[i].options[3]).strtodouble(hd);
      startlambda = hd;

      if (f==1)
        return true;

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_season.raw","_season.res","_season");

      fcnonpgaussian.push_back(
      FULLCOND_nonp_gaussian(&generaloptions,D.getCol(j),
                                         unsigned(maxint.getvalue()),
                                         MCMC::seasonal,title,pathres,
                                         lambda,startlambda,per));

      fcnonpgaussian[fcnonpgaussian.size()-1].init_name(terms[i].varnames[0]);
      fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_nonpseason_varcoef ---------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_nonpseason_varcoef(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double hd;
  double lambda, startlambda;
  unsigned per;
  int f;

  unsigned i;
  int j1, j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpseason_varcoef.checkvector(terms,i) == true )
      {
      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interacting var
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv); // effectmod

      f = (terms[i].options[1]).strtolong(h);
      per = unsigned(h);
      f = (terms[i].options[2]).strtodouble(hd);
      lambda = hd;
      f = (terms[i].options[3]).strtodouble(hd);
      startlambda = hd;

      if (f==1)
        return true;

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],"_season.raw","_season.res",
                 "_season");

      fcnonpgaussian.push_back(
      FULLCOND_nonp_gaussian(&generaloptions,D.getCol(j2),D.getCol(j1),
                                         unsigned(maxint.getvalue()),
                                         MCMC::seasonal,title,pathres,
                                         lambda,startlambda,per));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fcnonpgaussian[fcnonpgaussian.size()-1].init_names(na);
      fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// --------------------------- create_spatial ---------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_spatial(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  double hd;
  int f;
  double lambda, startlambda;
  unsigned i;
  int j;

  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial.checkvector(terms,i) == true )
      {
      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      mapobject * mapp;                           // pointer to mapobject
      int objpos = findstatobject(*statobj,terms[i].options[1],"map");
      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          {
          if ((terms[i].options[1] == "") || (terms[i].options[1] == " "))
            outerror("ERROR: map object must be specified to estimate a spatial effect\n");
          else
            outerror("ERROR: map object " + terms[i].options[1] + " is not existing\n");
          }
        else
          outerror("ERROR: " + terms[i].options[1] + " is not a map object\n");
        return true;
        }
      MAP::map m = mapp->getmap();

      if (m.isconnected()==false)
        {
        outerror("ERROR: map is disconnected, spatial effect cannot be estimated\n");
        return true;
        }

      f = (terms[i].options[2]).strtodouble(hd);
      lambda = hd;
      f = (terms[i].options[3]).strtodouble(hd);
      startlambda = hd;

      if (f==1)
        return true;

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_spatial.raw","_spatial.res","_spatial");

      fcnonpgaussian.push_back(
      FULLCOND_nonp_gaussian(&generaloptions,D.getCol(j),m,terms[i].options[1],
                             title,pathnonp,pathres,lambda,startlambda)
                           );

      if (fcnonpgaussian[fcnonpgaussian.size()-1].get_errors().size() > 0)
        {
        unsigned i;
        for(i=0;i<fcnonpgaussian[fcnonpgaussian.size()-1].get_errors().size();i++)
          errormessages.push_back(fcnonpgaussian[fcnonpgaussian.size()-1].get_errors()[i]);
        return true;
        }

      fcnonpgaussian[fcnonpgaussian.size()-1].init_name(terms[i].varnames[0]);
      fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ------------------------- create_spatial_varcoef ----------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_spatial_varcoef(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  double hd;
  int f;
  double lambda, startlambda;
  unsigned i;
  int j1, j2;

  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial_varcoef.checkvector(terms,i) == true )
      {
      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interacting var
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv); // effectmod

      mapobject * mapp;                           // pointer to mapobject
      int objpos = findstatobject(*statobj,terms[i].options[1],"map");
      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          {
          if ((terms[i].options[1] == "") || (terms[i].options[1] == " "))
            outerror("ERROR: map object must be specified to estimate a spatial effect\n");
          else
            outerror("ERROR: map object " + terms[i].options[1] + " is not existing\n");
          }
        else
          outerror("ERROR: " + terms[i].options[1] + " is not a map object\n");
        return true;
        }
      MAP::map m = mapp->getmap();

      f = (terms[i].options[2]).strtodouble(hd);
      lambda = hd;
      f = (terms[i].options[3]).strtodouble(hd);
      startlambda = hd;

      if (f==1)
        return true;

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],"_spatial.raw","_spatial.res","_spatial");

      fcnonpgaussian.push_back(
      FULLCOND_nonp_gaussian(&generaloptions,D.getCol(j2),D.getCol(j1),m,
                             terms[i].options[1],title,pathnonp,pathres,lambda,
                             startlambda)
                           );

      if (fcnonpgaussian[fcnonpgaussian.size()-1].get_errors().size() > 0)
        {
        unsigned i;
        for(i=0;i<fcnonpgaussian[fcnonpgaussian.size()-1].get_errors().size();i++)
          errormessages.push_back(fcnonpgaussian[fcnonpgaussian.size()-1].get_errors()[i]);
        return true;
        }

      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fcnonpgaussian[fcnonpgaussian.size()-1].init_names(na);
      fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// -------------------------- create_pspline -----------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_pspline(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string title;
  long h;
  unsigned degree,nrknots;
  double lambda, startlambda;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonppspline.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "psplinerw1")
        type = MCMC::RW1;
      else
        type = MCMC::RW2;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      degree = unsigned(h);
      f = (terms[i].options[2]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[3]).strtodouble(lambda);
      f = (terms[i].options[7]).strtodouble(startlambda);

      if (f==1)
        return true;

      MCMC::knotpos po;
      if (knots.getvalue() == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline.raw","_pspline.res","_pspline");

      fcpspline.push_back( spline_basis(&generaloptions,
                                              D.getCol(j),
                                              nrknots,
                                              degree,
                                              po,
                                              type,
                                              title,
                                              pathnonp,
                                              pathres,
                                              lambda,
                                              startlambda
                                             )
                           );
      fcpspline[fcpspline.size()-1].init_name(terms[i].varnames[0]);
      fcpspline[fcpspline.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcpspline[fcpspline.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ------------------------- create_varcoeffpspline ----------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_varcoeffpspline(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string monotone;
  long h;
  unsigned degree,nrknots;
  double lambda,startlambda;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpvarcoeffpspline.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "varpsplinerw1")
        type = MCMC::RW1;
      else
        type = MCMC::RW2;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interacting var
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv); // effectmod

      f = (terms[i].options[1]).strtolong(h);
      degree = unsigned(h);
      f = (terms[i].options[2]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[3]).strtodouble(lambda);
      f = (terms[i].options[4]).strtodouble(startlambda);

      if (f==1)
        return true;

      MCMC::knotpos po;
      if (knots.getvalue() == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_pspline.raw","_pspline.res","_pspline");

      fcpspline.push_back( spline_basis(&generaloptions,
                                              D.getCol(j2),
                                              D.getCol(j1),
                                              nrknots,
                                              degree,
                                              po,
                                              type,
                                              title,
                                              pathnonp,
                                              pathres,
                                              lambda,
                                              startlambda
                                             )
                           );
      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fcpspline[fcpspline.size()-1].init_names(na);
      fcpspline[fcpspline.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcpspline[fcpspline.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ------------------- create_interactionpspline -------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_interactionspspline(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double lambda,startlambda;
  unsigned nrknots,degree;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpinteractpspline.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "pspline2dimrw1")
        type = MCMC::mrflinear;
      else if (terms[i].options[0] == "pspline2dimrw2")
        type = MCMC::mrfquadratic8;
      else if( terms[i].options[0] =="pspline2dimbiharmonic")
        type = MCMC::mrfquadratic12;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      degree = unsigned(h);
      f = (terms[i].options[2]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[3]).strtodouble(lambda);
      f = (terms[i].options[4]).strtodouble(startlambda);

      if (f==1)
        return true;

      ST::string title;

      ST::string help  = terms[i].varnames[0] + "_" + terms[i].varnames[1];

      make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_pspline.raw","_pspline.res","_pspline");

      fcpsplinesurf.push_back(
      spline_basis_surf(&generaloptions,
                                      D.getCol(j1),
                                      D.getCol(j2),
                                      nrknots,degree,type,
                                      title,
                                      pathnonp,
                                      pathres,
                                      lambda,
                                      startlambda
                                      ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fcpsplinesurf[fcpsplinesurf.size()-1].init_names(na);
      fcpsplinesurf[fcpsplinesurf.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcpsplinesurf[fcpsplinesurf.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ---------------------- create_varcoeffinteractionpspline --------------------
// -----------------------------------------------------------------------------

bool remlreg::create_varcoeffinteractionspspline(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double lambda,startlambda;
  unsigned nrknots,degree;
  int f;

  unsigned i;
  int j1,j2,j3;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpvarcoeffinteractpspline.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "varpspline2dimrw1")
        type = MCMC::mrflinear;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);  //interaction variable
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);
      j3 = terms[i].varnames[2].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      degree = unsigned(h);
      f = (terms[i].options[2]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[3]).strtodouble(lambda);
      f = (terms[i].options[4]).strtodouble(startlambda);

      if (f==1)
        return true;

      ST::string title;

      ST::string help  = terms[i].varnames[1] + "_" + terms[i].varnames[2];

      make_paths(collinpred,pathnonp,pathres,title,help,terms[i].varnames[0],
                 "_pspline.raw","_pspline.res","_pspline");

      fcpsplinesurf.push_back(
      spline_basis_surf(&generaloptions,
                                      D.getCol(j1),
                                      D.getCol(j2),
                                      D.getCol(j3),
                                      nrknots,degree,type,
                                      title,
                                      pathnonp,
                                      pathres,
                                      lambda,
                                      startlambda
                                      ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[2]);
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fcpsplinesurf[fcpsplinesurf.size()-1].init_names(na);
      fcpsplinesurf[fcpsplinesurf.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcpsplinesurf[fcpsplinesurf.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ------------------------- create_geospline ----------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_geospline(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double lambda,startlambda;
  unsigned nrknots,degree;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpgeospline.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "geosplinerw1")
        type = MCMC::mrflinear;
      else if (terms[i].options[0] == "geosplinerw2")
        type = MCMC::mrfquadratic8;
      else if (terms[i].options[0] == "geosplinebiharmonic")
        type = MCMC::mrfquadratic12;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      degree = unsigned(h);
      f = (terms[i].options[2]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[3]).strtodouble(lambda);
      f = (terms[i].options[5]).strtodouble(startlambda);

      if (f==1)
        return true;

      mapobject * mapp;                           // pointer to mapobject
      int objpos = findstatobject(*statobj,terms[i].options[4],"map");
      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          outerror("ERROR: map object " + terms[i].options[4] + " is not existing\n");
        else
          outerror("ERROR: " + terms[i].options[4] + " is not a map object\n");
        return true;
        }
      MAP::map m = mapp->getmap();

      if(!m.centroids_existing())
        {
        outerror("ERROR: map object doesn�t contain centroids\n");
        return true;
        }

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_geospline.raw","_geospline.res","_geospline");

      fcpsplinesurf.push_back(
      spline_basis_surf(&generaloptions,
                                      D.getCol(j),m,terms[i].options[4],
                                      nrknots,degree,
                                      type,
                                      title,
                                      pathnonp,
                                      pathres,
                                      lambda,
                                      startlambda
                                      ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[0]);
      fcpsplinesurf[fcpsplinesurf.size()-1].init_names(na);
      fcpsplinesurf[fcpsplinesurf.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcpsplinesurf[fcpsplinesurf.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ------------------------- create_geospline_varcoef --------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_geospline_varcoeff(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double lambda,startlambda;
  unsigned nrknots,degree;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpvarcoeffgeospline.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type = MCMC::mrflinear;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interacting var
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv); // effectmod

      f = (terms[i].options[1]).strtolong(h);
      degree = unsigned(h);
      f = (terms[i].options[2]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[3]).strtodouble(lambda);
      f = (terms[i].options[5]).strtodouble(startlambda);

      if (f==1)
        return true;

      mapobject * mapp;                           // pointer to mapobject
      int objpos = findstatobject(*statobj,terms[i].options[4],"map");
      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          outerror("ERROR: map object " + terms[i].options[4] + " is not existing\n");
        else
          outerror("ERROR: " + terms[i].options[4] + " is not a map object\n");
        return true;
        }
      MAP::map m = mapp->getmap();

      if(!m.centroids_existing())
        {
        outerror("ERROR: map object doesn�t contain centroids\n");
        return true;
        }

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_geospline.raw","_geospline.res","_geospline");

      fcpsplinesurf.push_back(
      spline_basis_surf(&generaloptions,
                                      D.getCol(j1),D.getCol(j2),
                                      m,terms[i].options[4],
                                      nrknots,degree,
                                      type,
                                      title,
                                      pathnonp,
                                      pathres,
                                      lambda,
                                      startlambda
                                      ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fcpsplinesurf[fcpsplinesurf.size()-1].init_names(na);
      fcpsplinesurf[fcpsplinesurf.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcpsplinesurf[fcpsplinesurf.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_kriging --------------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_kriging(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double nu,maxdist,p,q,lambda,startlambda;
  unsigned nrknots, maxsteps;
  bool full;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial_kriging.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "kriging")
        type = MCMC::kriging;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[2]).strtodouble(nu);
      if(nu!=0.5 && nu!=1.5&& nu!=2.5 && nu!=3.5)
        {
        outerror("ERROR: Invalid value for nu\n");
        return true;
        }
      f = (terms[i].options[3]).strtodouble(maxdist);
      if(maxdist<=0) // w�hle maxdist so, dass Korrelation f�r Punkte mit maximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;
          }
        }
      if(terms[i].options[4] == "true")
        {
        full=true;
        }
      else
        {
        full=false;
        }
      f = (terms[i].options[6]).strtodouble(p);
      f = (terms[i].options[7]).strtodouble(q);
      f = (terms[i].options[8]).strtolong(h);
      maxsteps = unsigned(h);
      f = (terms[i].options[9]).strtodouble(lambda);
      f = (terms[i].options[10]).strtodouble(startlambda);

      if (f==1)
        return true;

      // read knots from a specified dataset object
      datamatrix knotdata;
      if(terms[i].options[5]!="" && !full)
        {
        dataobject * datap;                           // pointer to datsetobject
        int objpos = findstatobject(*statobj,terms[i].options[5],"dataset");
        if (objpos >= 0)
          {
          statobject * s = statobj->at(objpos);
          datap = dynamic_cast<dataobject*>(s);
          if (datap->obs()==0 || datap->getVarnames().size()==0)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " does not contain any data\n");
            return true;
            }
          else if (datap->getVarnames().size()>2)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " contains more than two variables\n");
            return true;
            }
          }
        else
          {
          outerror("ERROR: dataset object " + terms[i].options[5] + " is not existing\n");
          return true;
          }
        list<ST::string> knotnames = datap->getVarnames();
        ST::string expr = "";
        datap->makematrix(knotnames,knotdata,expr);
        }
      else
        {
        knotdata = datamatrix(1,1,0);
        }

      ST::string title;

      ST::string help  = terms[i].varnames[0] + "_" + terms[i].varnames[1];

      make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_kriging.raw","_kriging.res","_kriging");

      fckriging.push_back(
      FULLCOND_kriging(&generaloptions,
              D.getCol(j1),
              D.getCol(j2),
              knotdata,
              nrknots,nu,maxdist,p,q,maxsteps,full,type,
              title,
              pathnonp,
              pathres,
              lambda,
              startlambda
              ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fckriging[fckriging.size()-1].init_names(na);
      fckriging[fckriging.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fckriging[fckriging.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_kriging_1dim ---------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_kriging_1dim(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  double nu,maxdist,lambda,startlambda;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonp_kriging.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "1dimkriging")
        type = MCMC::kriging;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtodouble(nu);
      if(nu!=0.5 && nu!=1.5&& nu!=2.5 && nu!=3.5)
        {
        outerror("ERROR: Invalid value for nu\n");
        return true;
        }
      f = (terms[i].options[2]).strtodouble(maxdist);
      if(maxdist<=0) // w�hle maxdist so, dass Korrelation f�r Punkte mitmaximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;
          }
        }
      f = (terms[i].options[3]).strtodouble(lambda);
      f = (terms[i].options[4]).strtodouble(startlambda);

      if (f==1)
        return true;

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_kriging.raw","_kriging.res","_kriging");

      fckriging.push_back(
      FULLCOND_kriging(&generaloptions,
              D.getCol(j),
              nu,maxdist,type,
              title,
              pathnonp,
              pathres,
              lambda,
              startlambda
              ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[0]);
      fckriging[fckriging.size()-1].init_names(na);
      fckriging[fckriging.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fckriging[fckriging.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_kriging_varcoeff -----------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_kriging_varcoeff(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double nu,maxdist,p,q,lambda,startlambda;
  unsigned nrknots, maxsteps;
  bool full;
  int f;

  unsigned i;
  int j1,j2,j3;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial_kriging_varcoeff.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "varkriging")
        type = MCMC::kriging;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); //interaction variable
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);
      j3 = terms[i].varnames[2].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      nrknots = unsigned(h);

      f = (terms[i].options[2]).strtodouble(nu);
      if(nu!=0.5 && nu!=1.5&& nu!=2.5 && nu!=3.5)
        {
        outerror("ERROR: Invalid value for nu\n");
        return true;
        }
      f = (terms[i].options[3]).strtodouble(maxdist);
      if(maxdist<=0) // w�hle maxdist so, dass Korrelation f�r Punkte mitmaximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;
          }
        }
      if(terms[i].options[4] == "true")
        {
        full=true;
        }
      else
        {
        full=false;
        }

      f = (terms[i].options[6]).strtodouble(p);
      f = (terms[i].options[7]).strtodouble(q);
      f = (terms[i].options[8]).strtolong(h);
      maxsteps = unsigned(h);

      f = (terms[i].options[9]).strtodouble(lambda);
      f = (terms[i].options[10]).strtodouble(startlambda);

      if (f==1)
        return true;

      datamatrix knotdata;
      if(terms[i].options[5]!="" && !full)
        {
        dataobject * datap;                           // pointer to datsetobject
        int objpos = findstatobject(*statobj,terms[i].options[5],"dataset");
        if (objpos >= 0)
          {
          statobject * s = statobj->at(objpos);
          datap = dynamic_cast<dataobject*>(s);
          if (datap->obs()==0 || datap->getVarnames().size()==0)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " does not contain any data\n");
            return true;
            }
          else if (datap->getVarnames().size()>2)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " contains more than two variables\n");
            return true;
            }
          }
        else
          {
          outerror("ERROR: dataset object " + terms[i].options[5] + " is not existing\n");
          return true;
          }
        list<ST::string> knotnames = datap->getVarnames();
        ST::string expr = "";
        datap->makematrix(knotnames,knotdata,expr);
        }
      else
        {
        knotdata = datamatrix(1,1,0);
        }

      ST::string title;

      ST::string help  = terms[i].varnames[1] + "_" + terms[i].varnames[2];

      make_paths(collinpred,pathnonp,pathres,title,help,terms[i].varnames[0],
                 "_kriging.raw","_kriging.res","_kriging");

      fckriging.push_back(
      FULLCOND_kriging(&generaloptions,
              D.getCol(j1),
              D.getCol(j2),
              D.getCol(j3),
              knotdata,
              nrknots,nu,maxdist,p,q,maxsteps,full,type,
              title,
              pathnonp,
              pathres,
              lambda,
              startlambda
              ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[2]);
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fckriging[fckriging.size()-1].init_names(na);
      fckriging[fckriging.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fckriging[fckriging.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_geokriging -----------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_geokriging(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double nu,maxdist,p,q,lambda,startlambda;
  unsigned nrknots, maxsteps;
  bool full;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial_geokriging.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "geokriging")
        type = MCMC::kriging;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[2]).strtodouble(nu);
      if(nu!=0.5 && nu!=1.5&& nu!=2.5 && nu!=3.5)
        {
        outerror("ERROR: Invalid value for nu\n");
        return true;
        }
      f = (terms[i].options[3]).strtodouble(maxdist);
      if(maxdist<=0) // w�hle maxdist so, dass Korrelation f�r Punkte mitmaximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;
          }
        }
      if(terms[i].options[4] == "true")
        {
        full=true;
        }
      else
        {
        full=false;
        }
      f = (terms[i].options[6]).strtodouble(p);
      f = (terms[i].options[7]).strtodouble(q);
      f = (terms[i].options[8]).strtolong(h);
      maxsteps = unsigned(h);
      f = (terms[i].options[9]).strtodouble(lambda);
      f = (terms[i].options[10]).strtodouble(startlambda);

      if (f==1)
        return true;

      mapobject * mapp;                           // pointer to mapobject
      int objpos = findstatobject(*statobj,terms[i].options[11],"map");
      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          outerror("ERROR: map object " + terms[i].options[11] + " is not existing\n");
        else
          outerror("ERROR: " + terms[i].options[11] + " is not a map object\n");
        return true;
        }
      MAP::map m = mapp->getmap();
      if(!m.centroids_existing())
        {
        outerror("ERROR: map object doesn�t contain centroids\n");
        return true;
        }

      datamatrix knotdata;
      if(terms[i].options[5]!="" && !full)
        {
        dataobject * datap;                           // pointer to datasetobject
        int objpos = findstatobject(*statobj,terms[i].options[5],"dataset");
        if (objpos >= 0)
          {
          statobject * s = statobj->at(objpos);
          datap = dynamic_cast<dataobject*>(s);
          if (datap->obs()==0 || datap->getVarnames().size()==0)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " does not contain any data\n");
            return true;
            }
          else if (datap->getVarnames().size()>2)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " contains more than two variables\n");
            return true;
            }
          }
        else
          {
          outerror("ERROR: dataset object " + terms[i].options[5] + " is not existing\n");
          return true;
          }
        list<ST::string> knotnames = datap->getVarnames();
        ST::string expr = "";
        datap->makematrix(knotnames,knotdata,expr);
        }
      else
        {
        knotdata = datamatrix(1,1,0);
        }

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_geokriging.raw","_geokriging.res","_geokriging");

      fckriging.push_back(
      FULLCOND_kriging(&generaloptions,
              D.getCol(j),m,terms[i].options[11],
              knotdata,
              nrknots,nu,maxdist,p,q,maxsteps,full,type,
              title,
              pathnonp,
              pathres,
              lambda,
              startlambda
              ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[0]);
      fckriging[fckriging.size()-1].init_names(na);
      fckriging[fckriging.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fckriging[fckriging.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_geokriging_varcoeff --------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_geokriging_varcoeff(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double nu,maxdist,p,q,lambda,startlambda;
  unsigned nrknots, maxsteps;
  bool full;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial_geokriging_varcoeff.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "vargeokriging")
        type = MCMC::kriging;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interacting var
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv); // effectmod

      f = (terms[i].options[1]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[2]).strtodouble(nu);
      if(nu!=0.5 && nu!=1.5&& nu!=2.5 && nu!=3.5)
        {
        outerror("ERROR: Invalid value for nu\n");
        return true;
        }
      f = (terms[i].options[3]).strtodouble(maxdist);
      if(maxdist<=0) // w�hle maxdist so, dass Korrelation f�r Punkte mitmaximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;
          }
        }
      if(terms[i].options[4] == "true")
        {
        full=true;
        }
      else
        {
        full=false;
        }
      f = (terms[i].options[6]).strtodouble(p);
      f = (terms[i].options[7]).strtodouble(q);
      f = (terms[i].options[8]).strtolong(h);
      maxsteps = unsigned(h);
      f = (terms[i].options[9]).strtodouble(lambda);
      f = (terms[i].options[10]).strtodouble(startlambda);

      if (f==1)
        return true;

      mapobject * mapp;                           // pointer to mapobject
      int objpos = findstatobject(*statobj,terms[i].options[11],"map");
      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          outerror("ERROR: map object " + terms[i].options[11] + " is not existing\n");
        else
          outerror("ERROR: " + terms[i].options[11] + " is not a map object\n");
        return true;
        }
      MAP::map m = mapp->getmap();
      if(!m.centroids_existing())
        {
        outerror("ERROR: map object doesn�t contain centroids\n");
        return true;
        }

      datamatrix knotdata;
      if(terms[i].options[5]!="" && !full)
        {
        dataobject * datap;                           // pointer to datasetobject
        int objpos = findstatobject(*statobj,terms[i].options[5],"dataset");
        if (objpos >= 0)
          {
          statobject * s = statobj->at(objpos);
          datap = dynamic_cast<dataobject*>(s);
          if (datap->obs()==0 || datap->getVarnames().size()==0)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " does not contain any data\n");
            return true;
            }
          else if (datap->getVarnames().size()>2)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " contains more than two variables\n");
            return true;
            }
          }
        else
          {
          outerror("ERROR: dataset object " + terms[i].options[5] + " is not existing\n");
          return true;
          }
        list<ST::string> knotnames = datap->getVarnames();
        ST::string expr = "";
        datap->makematrix(knotnames,knotdata,expr);
        }
      else
        {
        knotdata = datamatrix(1,1,0);
        }

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_geokriging.raw","_geokriging.res","_geokriging");

      fckriging.push_back(
      FULLCOND_kriging(&generaloptions,
              D.getCol(j1),D.getCol(j2),m,terms[i].options[11],
              knotdata,
              nrknots,nu,maxdist,p,q,maxsteps,full,type,
              title,
              pathnonp,
              pathres,
              lambda,
              startlambda
              ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fckriging[fckriging.size()-1].init_names(na);
      fckriging[fckriging.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fckriging[fckriging.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_baseline -------------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_baseline(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string title;
  long h;
  unsigned degree,nrknots,tgrid,nrquant,nrbetween;
  double lambda, startlambda;
  int f;
  MCMC::knotpos gridpo;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonp_baseline.checkvector(terms,i) == true )
      {
      if(fcbaseline.size()>0)
        {
        outerror("ERROR: More than one baseline term specified!\n");
        return true;
        }

      MCMC::fieldtype type;
      if (terms[i].options[0] == "baseline")
        type = MCMC::RW2;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      // read lower interval boundary
      datamatrix lowerint;
      if(leftint.getvalue()!="")
        {
        lowerint = D.getCol(1);
        }
      else
        {
        lowerint = datamatrix(1,1,0);
        }

      // read left truncation time
      datamatrix lowertrunc;
      if(lefttrunc.getvalue()!="")
        {
        lowertrunc = D.getCol(2);
        }
      else
        {
        lowertrunc = datamatrix(1,1,0);
        }

/*      datamatrix lower;
      if(terms[i].options[9]!="")
        {
        dataobject * datap;                           // pointer to datsetobject
        int objpos = findstatobject(*statobj,terms[i].options[9],"dataset");
        if (objpos >= 0)
          {
          statobject * s = statobj->at(objpos);
          datap = dynamic_cast<dataobject*>(s);
          if (datap->obs()==0 || datap->getVarnames().size()==0)
            {
            outerror("ERROR: dataset object " + terms[i].options[9] + " does not contain any data\n");
            return true;
            }
          else if (datap->getVarnames().size()>1)
            {
            outerror("ERROR: dataset object " + terms[i].options[9] + " contains more than one variables\n");
            return true;
            }
          }
        else
          {
          outerror("ERROR: dataset object " + terms[i].options[9] + " is not existing\n");
          return true;
          }
        list<ST::string> lowname = datap->getVarnames();
        ST::string expr = "";
        datap->makematrix(lowname,lower,expr);
        }
      else
        {
        lower = datamatrix(1,1,0);
        }*/

      f = (terms[i].options[1]).strtolong(h);
      degree = unsigned(h);
      f = (terms[i].options[2]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[3]).strtolong(h);
      tgrid = unsigned(h);
      if(terms[i].options[4] == "equidistant")
        {
        gridpo = MCMC::equidistant;
        }
      else
        {
        gridpo = MCMC::quantiles;
        }
      f = (terms[i].options[5]).strtolong(h);
      nrquant = unsigned(h);
      f = (terms[i].options[6]).strtolong(h);
      nrbetween = unsigned(h);
      f = (terms[i].options[7]).strtodouble(lambda);
      f = (terms[i].options[8]).strtodouble(startlambda);

      if (f==1)
        return true;

      MCMC::knotpos po = MCMC::equidistant;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_logbaseline.raw","_logbaseline.res","_logbaseline");

      fcbaseline.push_back( baseline_reml(&generaloptions,
                                              D.getCol(j),
                                              lowerint,
                                              lowertrunc,
                                              nrknots,
                                              degree,
                                              tgrid,
                                              nrquant,
                                              nrbetween,
                                              po,
                                              type,
                                              title,
                                              pathnonp,
                                              pathres,
                                              lambda,
                                              startlambda,
                                              gridpo
                                             )
                           );

      fcbaseline[fcbaseline.size()-1].init_name(terms[i].varnames[0]);
      fcbaseline[fcbaseline.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcbaseline[fcbaseline.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_baseline_varcoeff-----------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_baseline_varcoeff(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string title;
  double lambda, startlambda;
  unsigned degree,nrknots,tgrid;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonp_baseline_varcoeff.checkvector(terms,i) == true )
      {
      if(fcbaseline.size()<1)
        {
        outerror("ERROR: Time-varying effects without baseline effect!\n");
        return true;
        }

      MCMC::fieldtype type = MCMC::RW2;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interacting var
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv); // effectmod

      f = (terms[i].options[1]).strtodouble(lambda);
      f = (terms[i].options[2]).strtodouble(startlambda);

      if (f==1)
        return true;

      tgrid = fcbaseline[0].get_tgrid();
      degree = fcbaseline[0].get_degree();
      nrknots = fcbaseline[0].get_nrknots();

      MCMC::knotpos po = MCMC::equidistant;

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_logbaseline.raw","_logbaseline.res","_logbaseline");

      fcbaseline_varcoeff.push_back( baseline_reml(&generaloptions,
                                              D.getCol(j2),
                                              D.getCol(j1),
                                              nrknots,
                                              degree,
                                              tgrid,
                                              po,
                                              type,
                                              title,
                                              pathnonp,
                                              pathres,
                                              lambda,
                                              startlambda
                                             )
                           );

      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fcbaseline_varcoeff[fcbaseline_varcoeff.size()-1].init_names(na);
      fcbaseline_varcoeff[fcbaseline_varcoeff.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcbaseline_varcoeff[fcbaseline_varcoeff.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ------------------------ create_random --------------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_random(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string title;
  double lambda, startlambda;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( randomeff.checkvector(terms,i) == true )
      {
      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtodouble(lambda);
      f = (terms[i].options[2]).strtodouble(startlambda);

      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_random.raw","_random.res","_random");

      fcrandom.push_back(FULLCOND_random(&generaloptions,
                         D.getCol(j),title,pathnonp,
                         pathres,lambda,startlambda));

      fcrandom[fcrandom.size()-1].init_name(terms[i].varnames[0]);
      fcrandom[fcrandom.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcrandom[fcrandom.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_randomslope ----------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_randomslope(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string title;
  unsigned i;
  int j1,j2;
  double lambda, startlambda;
  int f;
  for(i=0;i<terms.size();i++)
    {
    if ( randomeffslope.checkvector(terms,i) == true )
      {
      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtodouble(lambda);
      f = (terms[i].options[2]).strtodouble(startlambda);

      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_random.raw","_random.res","_random");

      fcrandom.push_back(FULLCOND_random(&generaloptions,D.getCol(j1),
                          D.getCol(j2),title,pathnonp,pathres,lambda,
                          startlambda));

      vector<ST::string> varnameshelp;
      varnameshelp.push_back(terms[i].varnames[1]);
      varnameshelp.push_back(terms[i].varnames[0]);
      fcrandom[fcrandom.size()-1].init_names(varnameshelp);
      fcrandom[fcrandom.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcrandom[fcrandom.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------------- remlrun ---------------------------------------
// -----------------------------------------------------------------------------

void remlrun(remlreg & b)
  {
  b.resultsyesno = false;
  b.terms = b.modreg.getterms();

  b.describetext.erase(b.describetext.begin(),b.describetext.end());
  b.describetext.push_back("LAST ESTIMATED MODEL: \n");
  b.describetext.push_back("\n");
  b.describetext.push_back(b.modreg.getModelText());
  b.describetext.push_back("\n");

  b.clear();

  #if defined(BORLAND_OUTPUT_WINDOW)
    bool failure = false;
  #elif defined(JAVA_OUTPUT_WINDOW)
    bool failure = b.adminb_p->breakcommand();
  #endif

  b.generaloptions = MCMCoptions(
  #if defined(JAVA_OUTPUT_WINDOW)
  b.adminb_p,
  #endif
  12000,2000,100,b.logout,b.level1.getvalue(),
                               b.level2.getvalue());

  ST::string header;
  bool dispers;

// Read design matrix, compute weights, etc.
  datamatrix weight;
  if (!failure)
    failure = b.create_data(weight);

// Define and check response
  datamatrix response;
  if(!failure)
    failure = b.create_response(response,weight);

// Compute offset
  datamatrix offset;
  if(!failure)
    failure = b.create_offset(offset);

// Compute different model terms
  if (!failure)
    failure = b.create_const(0);
  if( !failure)
    failure = b.create_baseline(0);
  if( !failure)
    failure = b.create_baseline_varcoeff(0);
  if (!failure)
    failure = b.create_nonprw1rw2(0);
  if (!failure)
    failure = b.create_nonprw1rw2_varcoef(0);
  if (!failure)
    failure = b.create_pspline(0);
  if (!failure)
    failure = b.create_nonpseason(0);
  if (!failure)
    failure = b.create_nonpseason_varcoef(0);
  if (!failure)
    failure = b.create_spatial(0);
  if (!failure)
    failure = b.create_spatial_varcoef(0);
  if (!failure)
    failure = b.create_geospline(0);
  if (!failure)
    failure = b.create_geospline_varcoeff(0);
  if (!failure)
    failure = b.create_varcoeffpspline(0);
  if (!failure)
    failure = b.create_random(0);
  if (!failure)
    failure = b.create_randomslope(0);
  if (!failure)
    failure = b.create_interactionspspline(0);
  if (!failure)
    failure = b.create_varcoeffinteractionspspline(0);
  if (!failure)
    failure = b.create_kriging(0);
  if (!failure)
    failure = b.create_kriging_1dim(0);
  if (!failure)
    failure = b.create_kriging_varcoeff(0);
  if (!failure)
    failure = b.create_geokriging(0);
  if (!failure)
    failure = b.create_geokriging_varcoeff(0);

  if (!failure)
    {
    header= "remlreg object " + b.name.to_bstr() + ": reml procedure" ;

// Nominale Modelle
    if (b.family.getvalue()=="multinomial")
      {
      b.RE_M = remlest_multinomial(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond,response,b.family.getvalue(),b.outfile.getvalue(),
      b.maxit.getvalue(),b.lowerlim.getvalue(),b.eps.getvalue(),b.cats,b.logout);
      if (b.fullcond.size() == 1)    // fixed effects only
        failure = b.RE_M.estimate_glm(response,offset,weight);
      else
        failure = b.RE_M.estimate(response,offset,weight);
      }
// Ordinale Modelle
    else if (b.family.getvalue()=="cumlogit" ||
        b.family.getvalue()=="cumprobit" ||
        b.family.getvalue()=="seqlogit" ||
        b.family.getvalue()=="seqprobit")
      {
      b.RE_O = remlest_ordinal(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond,response,b.family.getvalue(),b.outfile.getvalue(),
      b.maxit.getvalue(),b.lowerlim.getvalue(),b.eps.getvalue(),b.cats,b.logout);
      if (b.fullcond.size() == 1)    // fixed effects only
        failure = b.RE_O.estimate_glm(response,offset,weight);
      else
        failure = b.RE_O.estimate(response,offset,weight);
      }
// Univariate Modelle ohne Dispersionsparameter
    else if (b.family.getvalue()=="binomial" ||
        b.family.getvalue()=="binomialprobit" ||
        b.family.getvalue()=="binomialcomploglog" ||
        b.family.getvalue()=="poisson")
      {
      dispers=false;
      b.RE = remlest(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond,response,dispers,b.family.getvalue(),b.outfile.getvalue(),
      b.maxit.getvalue(),b.lowerlim.getvalue(),b.eps.getvalue(),b.logout);
      if (b.fullcond.size() == 1)    // fixed effects only
        failure = b.RE.estimate_glm(response,offset,weight);
      else
        failure = b.RE.estimate(response,offset,weight);
      }
// Cox-Modell
    else if (b.family.getvalue()=="cox")
      {
      dispers=false;
      b.RE = remlest(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond,response,dispers,b.family.getvalue(),b.outfile.getvalue(),
      b.maxit.getvalue(),b.lowerlim.getvalue(),b.eps.getvalue(),b.logout);
      failure = b.RE.estimate_survival(response,offset,weight);
      }
// Cox-Modell mit Intervallzensierung
    else if (b.family.getvalue()=="coxinterval")
      {
      dispers=false;
      b.RE = remlest(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond,response,dispers,b.family.getvalue(),b.outfile.getvalue(),
      b.maxit.getvalue(),b.lowerlim.getvalue(),b.eps.getvalue(),b.logout);
      failure = b.RE.estimate_survival_interval(response,offset,weight);
      }
// Cox-Modell mit Intervallzensierung & Linkstrunkierung
    else if (b.family.getvalue()=="coxinterval2")
      {
      dispers=false;
      b.RE = remlest(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond,response,dispers,b.family.getvalue(),b.outfile.getvalue(),
      b.maxit.getvalue(),b.lowerlim.getvalue(),b.eps.getvalue(),b.logout);
      failure = b.RE.estimate_survival_interval2(response,offset,weight);
      }
// Univariate Modelle mit Dispersionsparameter
    else
      {
      dispers=true;
      b.RE = remlest(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond,response,dispers,b.family.getvalue(),b.outfile.getvalue(),
      b.maxit.getvalue(),b.lowerlim.getvalue(),b.eps.getvalue(),b.logout);
      if (b.fullcond.size() == 1)    // fixed effects only
        failure = b.RE.estimate_glm_dispers(response,offset,weight);
      else
        failure = b.RE.estimate_dispers(response,offset,weight);
      }
    }

// Compute graphics
#if defined(JAVA_OUTPUT_WINDOW)
  if(!failure)
    {
    for(unsigned j=0;j<b.fullcond.size();j++)
      {
      MCMC::plotstyles plst = b.fullcond[j]->get_plotstyle();
      if(plst != MCMC::noplot)
        {
        vector<ST::string> varnames = b.fullcond[j]->get_datanames();
        ST::string xvar = varnames[0];
        ST::string effect = xvar;
        if(varnames.size()>1)
          {
          effect = varnames[1] + "*" + effect;
          }
        if(b.family.getvalue()=="multinomial")
          {
          for(unsigned i=0; i<b.cats.rows(); i++)
            {
            ST::string pathresult = b.fullcond[j]->get_pathresult();
            pathresult = pathresult.insert_after_string(ST::doubletostring(b.cats(i,0),6)+"_","_f_");
            ST::string pathps = pathresult.substr(0, pathresult.length()-4);
            if(plst == MCMC::plotnonp)
              {
              b.newcommands.push_back(b.name + ".plotnonp " + ST::inttostring(i*b.fullcond.size()+j)
              + ", title = \"Effect of " + effect +"\" xlab = " + xvar
              + " ylab = \" \" outfile = " + pathps + ".ps replace");
              }
            else if(plst==MCMC::drawmap)
              {
              double u = b.fullcond[j]->get_level1();
              double o = b.fullcond[j]->get_level2();
              ST::string u_str = ST::doubletostring(u,0);
              ST::string o_str = ST::doubletostring(o,0);
              b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(i*b.fullcond.size()+j)
              + ", color outfile = " + pathps + "_pmode.ps replace");
              b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(i*b.fullcond.size()+j)
              + ", plotvar = pcat" + u_str + " nolegend  pcat outfile = " + pathps
              + "_pcat" + u_str + ".ps replace");
              b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(i*b.fullcond.size()+j)
              + ", plotvar = pcat" + o_str + " nolegend  pcat outfile = " + pathps
              + "_pcat" + o_str + ".ps replace");
              }
            }
          }
        else
          {
          ST::string pathresult = b.fullcond[j]->get_pathresult();
          ST::string pathps = pathresult.substr(0, pathresult.length()-4);
          if(plst == MCMC::plotnonp)
            {
            b.newcommands.push_back(b.name + ".plotnonp " + ST::inttostring(j)
            + ", title = \"Effect of " + effect +"\" xlab = " + xvar
            + " ylab = \" \" outfile = " + pathps + ".ps replace");
            }
          else if(plst==MCMC::drawmap)
            {
            double u = b.fullcond[j]->get_level1();
            double o = b.fullcond[j]->get_level2();
            ST::string u_str = ST::doubletostring(u,0);
            ST::string o_str = ST::doubletostring(o,0);
            b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
            + ", color outfile = " + pathps + "_pmode.ps replace");
            b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
            + ", plotvar = pcat" + u_str + " nolegend  pcat outfile = " + pathps
            + "_pcat" + u_str + ".ps replace");
            b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
            + ", plotvar = pcat" + o_str + " nolegend  pcat outfile = " + pathps
            + "_pcat" + o_str + ".ps replace");
            }
          }
        }
      }
    }
#endif

// Produce batch-file for graphics and model summary in tex
  if(!failure)
    {
    ST::string path = b.outfile.getvalue() + "_graphics.prg";
    ST::string path2 = b.outfile.getvalue() + "_model_summary.tex";
    ST::string path3 = b.outfile.getvalue() +  "_splus.txt";

    if(b.family.getvalue()=="multinomial")
      {
      b.RE_M.make_graphics(header,path,path2,path3,
                         b.modreg.getModelVarnamesAsVector()[0].to_bstr());
      }
    else if(b.family.getvalue()=="cumlogit" || b.family.getvalue()=="cumprobit" ||
            b.family.getvalue()=="seqlogit" || b.family.getvalue()=="seqprobit")
      {
      b.RE_O.make_graphics(header,path,path2,path3,
                         b.modreg.getModelVarnamesAsVector()[0].to_bstr());
      }
    else
      {
      b.RE.make_graphics(header,path,path2,path3,
                         b.modreg.getModelVarnamesAsVector()[0].to_bstr(),
                         dispers);
      }
    }

  if (!failure)
    {
    b.resultsyesno = true;
    }
  else
    {
    b.describetext.erase(b.describetext.begin(),b.describetext.end());
    b.describetext.push_back("CURRENT REGRESSION RESULTS: none\n");
    b.resultsyesno = false;
    }
  }

// -----------------------------------------------------------------------------
// ----------------------------- drawmaprun ------------------------------------
// -----------------------------------------------------------------------------

void drawmaprun(remlreg & b)
  {
#if defined(BORLAND_OUTPUT_WINDOW)
  b.outerror("ERROR: method drawmap is not available in this version\n");
#elif defined(JAVA_OUTPUT_WINDOW)
  bool error = false;

  vector<ST::string> varnames = b.mdrawmap.getModelVarnamesAsVector();
  if (varnames.size() != 1)
    {
    b.outerror("ERROR: syntax error for method drawmap\n");
    error = true;
    }

  long nr;
  if (varnames[0].strtolong(nr) != 0)
    {
    b.outerror("ERROR: syntax error for method drawmap\n");
    error = true;
    }

  unsigned catnum;
  if(b.ismultinomial)
    {
    catnum = nr / b.fullcond.size();
    nr = nr % b.fullcond.size();
    if(catnum >= b.cats.rows())
      {
      b.outerror("ERROR: syntax error for method plotnonp\n");
      error = true;
      }
    }

  if (nr < 0 || nr >= b.fullcond.size())
    {
    b.outerror("ERROR: syntax error for method drawmap\n");
    error = true;
    }

  if (error == false)
    {
    if (b.fullcond[nr]->get_plotstyle() != MCMC::drawmap)
      {
      error = true;
      b.outerror("ERROR: results cannot be visualized with method drawmap\n");
      }
    }

  if (error==false)
    {

    ST::string path = b.fullcond[nr]->get_pathresult();
    if(b.ismultinomial)
      {
      path = path.insert_after_string(ST::doubletostring(b.cats(catnum,0),6)+"_","_f_");
      }

    vector<ST::string> vnames;
    ifstream in(path.strtochar());
    ST::string h;
    ST::getline(in,10000,h);
    vnames = 	h.strtoken(" ");

    ST::string graphname = "_" + b.name + "_graph";
    b.newcommands.push_back("graph " + graphname);

    ST::string datasetname = "_" + b.name + "_r0";
    b.newcommands.push_back("dataset " + datasetname);
    b.newcommands.push_back(datasetname + ".infile , nonote using " + path);

    ST::string plotvar;
      plotvar = b.plotvar.getvalue() + " " + vnames[1] + " ";

    ST::string ot="map=" + b.fullcond[nr]->getinfo() + " ";

    ot = ot + "nrcolors="+b.nrcolors.getValueAsString()+" ";
    ot = ot + "title=\""+b.title2.getvalue() + "\" ";
    if (b.outfile4.getvalue().length() > 0)
      ot = ot + "outfile=\""+b.outfile4.getvalue() + "\" ";
    if (b.nolegend.getvalue() == true)
      ot = ot + "nolegend ";
    if (b.color.getvalue() == true)
      ot = ot + "color ";
    if (b.swapcolors.getvalue() == true)
      ot = ot + "swapcolors ";
    if (b.replace.getvalue() == true)
      ot = ot + "replace ";
    if (b.lowerlimit.changed() == true)
      ot = ot + "lowerlimit="  + b.lowerlimit.getValueAsString() + " ";
    if (b.upperlimit.changed() == true)
      ot = ot + "upperlimit=" + b.upperlimit.getValueAsString() + " ";
    if (b.pcat.getvalue() == true)
      ot = ot + "pcat ";
    if (b.drawnames.getvalue() == true)
      ot = ot + "drawnames ";

    if (ot.length() == 0)
      b.newcommands.push_back(graphname + ".drawmap " + plotvar + " using " + datasetname);
    else
      b.newcommands.push_back(graphname + ".drawmap " + plotvar + "," + ot + " using "
      + datasetname);

    b.newcommands.push_back("drop " + graphname + " " + datasetname);
    }
#endif
  }

// -----------------------------------------------------------------------------
// -------------------------- plotnonprun --------------------------------------
// -----------------------------------------------------------------------------

void plotnonprun(remlreg & b)
  {
#if defined(BORLAND_OUTPUT_WINDOW)
  b.outerror("ERROR: method plotnonp is not available in this version\n");
#elif defined(JAVA_OUTPUT_WINDOW)
  bool error = false;

  vector<ST::string> varnames = b.mplotnonp.getModelVarnamesAsVector();
  if (varnames.size() != 1)
    {
    b.outerror("ERROR: syntax error for method plotnonp\n");
    error = true;
    }

  long nr;
  if (varnames[0].strtolong(nr) != 0)
    {
    b.outerror("ERROR: syntax error for method plotnonp\n");
    error = true;
    }

  unsigned catnum;
  if(b.ismultinomial)
    {
    catnum = nr / b.fullcond.size();
    nr = nr % b.fullcond.size();
    if(catnum >= b.cats.rows())
      {
      b.outerror("ERROR: syntax error for method plotnonp\n");
      error = true;
      }
    }

  if (nr < 0 || nr >= b.fullcond.size())
    {
    b.outerror("ERROR: syntax error for method plotnonp\n");
    error = true;
    }

  if (error == false)
    {
    if (b.fullcond[nr]->get_plotstyle() != MCMC::plotnonp)
      {
      error = true;
      b.outerror("ERROR: results cannot be visualized with method plotnonp\n");
      }
    }

  if (error==false)
    {

    ST::string path = b.fullcond[nr]->get_pathresult();

    if(b.ismultinomial)
      {
      path = path.insert_after_string(ST::doubletostring(b.cats(catnum,0),6)+"_","_f_");
      }

    vector<ST::string> vnames;
    ifstream in(path.strtochar());
    ST::string h;
    ST::getline(in,10000,h);
    vnames = 	h.strtoken(" ");

    ST::string graphname = "_" + b.name + "_graph";
    b.newcommands.push_back("graph " + graphname);

    ST::string datasetname = "_" + b.name + "_r0";
    b.newcommands.push_back("dataset " + datasetname);
    b.newcommands.push_back(datasetname + ".infile , nonote using " + path);

    ST::string plotvar;
    if (b.levels.getvalue()=="all")
      {
      plotvar = vnames[1] + " ";
      if (b.median.getvalue() == true)
        plotvar = plotvar + vnames[5] + " ";
      else
        plotvar = plotvar + vnames[2] + " ";
      plotvar = plotvar + vnames[3] + " " +
                         vnames[4] + " " +
                         vnames[6] + " " +
                         vnames[7] + " ";
      }
    else if (b.levels.getvalue()=="none")
      {
      if (b.median.getvalue() == true)
        plotvar = vnames[1] + " " + vnames[5];
      else
        plotvar = vnames[1] + " " + vnames[2];

      }
    else if (b.levels.getvalue()=="1")
      {
      if (b.median.getvalue() == true)
        plotvar = vnames[1] + " " + vnames[5] + " " + vnames[3] + " " +
                  vnames[7] + " ";
      else
        plotvar = vnames[1] + " " + vnames[2] + " " + vnames[3] + " " +
                  vnames[7] + " ";
      }
    else
      {
      if (b.median.getvalue() == true)
        plotvar = vnames[1] + " " + vnames[5] + " " + vnames[4] + " " +
                  vnames[6] + " ";
      else
        plotvar = vnames[1] + " " + vnames[2] + " " + vnames[4] + " " +
                  vnames[6] + " ";
      }


    ST::string ot;
    ot = "xlab=\""+b.xlab.getvalue() + "\" ";
    ot = ot + "ylab=\""+b.ylab.getvalue() + "\" ";
    ot = ot + "title=\""+b.title.getvalue() + "\" ";
    if (b.outfile2.getvalue().length() > 0)
      ot = ot + "outfile=\""+b.outfile2.getvalue() + "\" ";
    ot = ot + "height="+b.height.getValueAsString() + " ";
    ot = ot + "width="+b.width.getValueAsString() + " ";
    if (b.replace2.getvalue() == true)
      ot = ot + " replace ";
    if (b.connect.changed() == true)
      ot = ot + "connect="+b.connect.getvalue() + " ";
    if (b.ylimbottom.changed() == true)
      ot = ot + "ylimbottom="+b.ylimbottom.getValueAsString() + " ";
    if (b.ylimtop.changed() == true)
      ot = ot + "ylimtop="+b.ylimtop.getValueAsString() + " ";
    if (b.ystep.changed() == true)
      ot = ot + "ystep="+b.ystep.getValueAsString() + " ";
    if (b.ystart.changed() == true)
      ot = ot + "ystart="+b.ystart.getValueAsString() + " ";
    if (b.xlimbottom.changed() == true)
      ot = ot + "xlimbottom="+b.xlimbottom.getValueAsString() + " ";
    if (b.xlimtop.changed() == true)
      ot = ot + "xlimtop="+b.xlimtop.getValueAsString() + " ";
    if (b.xstep.changed() == true)
      ot = ot + "xstep="+b.xstep.getValueAsString() + " ";
    if (b.xstart.changed() == true)
      ot = ot + "xstart="+b.xstart.getValueAsString() + " ";
    if (b.linewidth.changed() == true)
      ot = ot + "linewidth="+b.linewidth.getValueAsString() + " ";
    if (b.fontsize.changed() == true)
      ot = ot + "fontsize="+b.fontsize.getValueAsString() + " ";
    if (b.pointsize.changed() == true)
      ot = ot + "pointsize="+b.pointsize.getValueAsString() + " ";
    if (b.linecolor.changed() == true)
      ot = ot + "linecolor="+b.linecolor.getValueAsString() + " ";
    if (b.titlescale.changed() == true)
      ot = ot + "titlesize="+b.titlescale.getValueAsString() + " ";

    if (ot.length() == 0)
      b.newcommands.push_back(graphname + ".plot " + plotvar + " using " + datasetname);
    else
      b.newcommands.push_back(graphname + ".plot " + plotvar + "," + ot + " using "
      + datasetname);

    b.newcommands.push_back("drop " + graphname + " " + datasetname);
    }
#endif
  b.plotnonpoptions.setdefault();
  }

#if defined(BORLAND_OUTPUT_WINDOW)
#pragma package(smart_init)
#endif














