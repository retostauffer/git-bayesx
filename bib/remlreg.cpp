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


void remlreg::create(void)
  {

  add_name="";

  ST::string h = defaultpath+"\\output\\"+name;

  outfile = fileoption("outfile",h,false);

  globaloptions.push_back(&outfile);

  // SYNTAX OF COMMANDS:
  // name [model] [weight varname] [by varname] [if expression]
  //      [, options] [using usingtext]

  // --------------------------- method regress --------------------------------

  offset = term_offset();
  fixedeffects = basic_termtype();
  nonprw1rw2 = term_autoreg_remlreg();
  nonprw1rw2_varcoef = term_autoreg_varcoef_remlreg();
  nonpseason = term_season_remlreg();
  nonpseason_varcoef = term_season_varcoef_remlreg();
  nonppspline = term_pspline_remlreg();
  nonpspatial = term_spatial_remlreg();
  nonpspatial_varcoef = term_spatial_varcoef_remlreg();
  nonpspatialxy = term_spatialxy();
  randomeff = term_random_remlreg();
  randomeffslope = term_randomslope_remlreg();
  nonpvarcoeffpspline = term_varcoeff_pspline_remlreg();
  nonpinteractpspline = term_interactpspline_remlreg();
  nonpvarcoeffinteractpspline = term_interactpspline_varcoeff_remlreg();
  nonpgeospline = term_geospline_remlreg();
  nonpvarcoeffgeospline = term_geospline_varcoeff_remlreg();
  nonpspatial_kriging = term_kriging_remlreg();
  nonp_kriging = term_kriging_1dim_remlreg();
  nonpspatial_kriging_varcoeff = term_kriging_varcoeff_remlreg();
  nonpspatial_geokriging = term_geokriging_remlreg();
  nonpspatial_geokriging_varcoeff = term_geokriging_varcoeff_remlreg();
  nonp_baseline = term_baseline_remlreg();
  nonp_baseline_varcoeff = term_baseline_varcoeff_remlreg();

  termtypes.push_back(&offset);
  termtypes.push_back(&fixedeffects);
  termtypes.push_back(&nonprw1rw2);
  termtypes.push_back(&nonpseason);
  termtypes.push_back(&nonprw1rw2_varcoef);
  termtypes.push_back(&nonpseason_varcoef);
  termtypes.push_back(&nonppspline);
  termtypes.push_back(&nonpspatial);
  termtypes.push_back(&nonpspatial_varcoef);
  termtypes.push_back(&randomeff);
  termtypes.push_back(&randomeffslope);
  termtypes.push_back(&nonpvarcoeffpspline);
  termtypes.push_back(&nonpinteractpspline);
  termtypes.push_back(&nonpvarcoeffinteractpspline);
  termtypes.push_back(&nonpspatialxy);
  termtypes.push_back(&nonpgeospline);
  termtypes.push_back(&nonpvarcoeffgeospline);
  termtypes.push_back(&nonpspatial_kriging);
  termtypes.push_back(&nonp_kriging);
  termtypes.push_back(&nonpspatial_kriging_varcoeff);
  termtypes.push_back(&nonpspatial_geokriging);
  termtypes.push_back(&nonpspatial_geokriging_varcoeff);
  termtypes.push_back(&nonp_baseline);
  termtypes.push_back(&nonp_baseline_varcoeff);

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
  families.push_back("poisson");
  families.push_back("gamma");
  families.push_back("poissondispers");
  families.push_back("binomialdispers");
  families.push_back("binomialprobitdispers");

  families.push_back("multinomial");
//  families.push_back("multinomialprobit");
  families.push_back("cumlogit");
  families.push_back("cumprobit");
  families.push_back("cox");
  family = stroption("family",families,"binomial");
  maxit = intoption("maxit",400,1,100000);
  lowerlim = doubleoption("lowerlim",0.001,0,1);
  eps = doubleoption("eps",0.00001,0,1);
  reference = doubleoption("reference",0,-10000,10000);

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
  // methods 0

  methods.push_back(command("regress",&modreg,&regressoptions,&udata,required,
			 optional,optional,optional,optional,required));

  functions[0] = remlrun;


  // --------------------------- method plotnonp -------------------------------

  uplotnonp = use();
  mplotnonp = modelStandard();

  xlab = stroption("xlab");
  ylab = stroption("ylab");
  height = intoption("height",210,0,500);
  width = intoption("width",356,0,500);
  ylimtop = doubleoption("ylimtop",-MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  ylimbottom = doubleoption("ylimbottom",MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  xlimtop = doubleoption("xlimtop",-MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  xlimbottom = doubleoption("xlimbottom",MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  xstep = doubleoption("xstep",0.0,-MAXDOUBLE,MAXDOUBLE);
  ystep = doubleoption("ystep",0.0,-MAXDOUBLE,MAXDOUBLE);

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
  plotnonpoptions.push_back(&height);
  plotnonpoptions.push_back(&width);
  plotnonpoptions.push_back(&ylimtop);
  plotnonpoptions.push_back(&ylimbottom);
  plotnonpoptions.push_back(&ystep);
  plotnonpoptions.push_back(&xlimtop);
  plotnonpoptions.push_back(&xlimbottom);
  plotnonpoptions.push_back(&xstep);
  plotnonpoptions.push_back(&levels);
  plotnonpoptions.push_back(&median);
  plotnonpoptions.push_back(&title);
  plotnonpoptions.push_back(&outfile2);
  plotnonpoptions.push_back(&replace2);


  // SYNTAX OF COMMANDS:
  // name [model] [weight varname] [by varname] [if expression]
  //      [, options] [using usingtext]

  // methods 1
  methods.push_back(command("plotnonp",&mplotnonp,&plotnonpoptions,&uplotnonp,
                   required,notallowed,notallowed,notallowed,optional,
                   notallowed));

  functions[1] = plotnonprun;

  // -------------------------- method drawmaprun ------------------------------

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

  // SYNTAX OF COMMANDS:
  // name [model] [weight varname] [by varname] [if expression]
  //      [, options] [using usingtext]

  // methods 2
  methods.push_back(command("drawmap",&mdrawmap,&drawmapoptions,&udrawmap,
                   required,notallowed,notallowed,notallowed,optional,
                   notallowed));

  functions[2] = drawmaprun;

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


bool remlreg::create_offset(datamatrix & o)
  {
  unsigned i;
  for(i=0;i<terms.size();i++)
    {
    if ( offset.checkvector(terms,i) == true)
      {
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
      o = datamatrix(D.rows(),1,0);
      }
    }
  return false;
  }


bool remlreg::create_data(datamatrix & weight)
  {
  unsigned i;
  bool failure=false;


  //--------------------- reading dataset information --------------------------

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

  //------------------ end: reading dataset information ------------------------


  //---------------- reading data, creating designmatrices ---------------------

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

    } // end: if ((datap->allexisting(modelvarnamesv,notex)) == false)


  datap->makematrix(modelvarnamesv,D,ifexpression);

  errormessages = datap->geterrormessages();

  if (!errormessages.empty())
    return true;

  if (wpos==-1)
    {
    weight = datamatrix(D.rows(),1,1);
    }
  else
    {
    weight = D.getCol(wpos);
    }

  return false;
  }

bool remlreg::create_response(datamatrix & response, datamatrix & weight)
  {
  response = D.getCol(0);

// Transform response for binomial families
  if(family.getvalue()=="binomial" ||
     family.getvalue()=="binomialprobit" ||
     family.getvalue()=="binomialdispers" ||
     family.getvalue()=="binomialprobitdispers")
    {
    unsigned i;
    for(i=0; i<response.rows(); i++)
      {
      response(i,0) = response(i,0)/weight(i,0);
      }
    }

// Read categories for multicategorial response
  if (family.getvalue()=="multinomial")
    {
    ismultinomial=true;
    }
  else
    {
    ismultinomial=false;
    }
  if (family.getvalue()=="multinomial" || family.getvalue()=="cumlogit" ||
      family.getvalue()=="cumprobit")
    {
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
    if(family.getvalue()=="multinomial")
      {
      refpos=0;//categories.size()-1;
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

bool remlreg::create_const(const unsigned & collinpred)
  {

  unsigned i;
  int j;

  vector<ST::string> varnames;
  vector<ST::string> varnamesh =  fixedeffects.get_constvariables(terms);

  varnames.push_back("const");

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

      }   // end: if ( nonpspatial.checkvector(terms,i) == true )

    } //  end:  for(i=0;i<terms.size();i++)

  return false;

  }

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

      }   // end: if ( nonpspatial.checkvector(terms,i) == true )

    } //  end:  for(i=0;i<terms.size();i++)

  return false;

  }

bool remlreg::create_spatialxy(const unsigned & collinpred)
  {

  /*

  ST::string pathnonp;
  ST::string pathres;
  ST::string pathmap;
  ST::string mapname;

  long h;
  unsigned min,max;
  double a1,b1;
  double lambda;
  double maxdist;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatialxy.checkvector(terms,i) == true )
      {

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);


      f = (terms[i].options[1]).strtolong(h);
      min = unsigned(h);

      f = (terms[i].options[2]).strtolong(h);
      max = unsigned(h);

      f = (terms[i].options[3]).strtodouble(lambda);

      f = (terms[i].options[4]).strtodouble(a1);

      f = (terms[i].options[5]).strtodouble(b1);

      f = (terms[i].options[6]).strtodouble(maxdist);

      ST::string title;

      pathmap = outfile.getvalue()+ add_name + "_" + terms[i].varnames[0] + "_" +
                terms[i].varnames[1] + "_dist" + ST::doubletostring(maxdist)
                + ".bnd";

      mapname = terms[i].varnames[0] + "_" + terms[i].varnames[1] + "_dist"
                + ST::doubletostring(maxdist) + add_name;

      if (collinpred == 0)
        {
        pathnonp = defaultpath + "\\temp\\" + name + add_name + "_f_" +
                   terms[i].varnames[0] + "_" + terms[i].varnames[1] +
                   "_spatial.raw";

        pathres = outfile.getvalue() + add_name + "_f_" + terms[i].varnames[0] + "_" +
                  terms[i].varnames[1] + "_spatial.res";

        title = "f_" + terms[i].varnames[0] + "_" + terms[i].varnames[1] + add_name;

        }
      else
        {
        pathnonp = defaultpath + "\\temp\\" + name + add_name + terms[i].varnames[0] +
                   terms[i].varnames[1] + "_spatial_" +
                   ST::inttostring(collinpred+1) + ".raw";

        pathres = outfile.getvalue() + add_name + terms[i].varnames[0] +
                  terms[i].varnames[1] + "_spatial_" +
                  ST::inttostring(collinpred+1) + ".res";

        title = "f_" + ST::inttostring(collinpred+1) +  "_" +
                terms[i].varnames[0] + "_" + terms[i].varnames[1]+add_name;

        }


      if (
          (family.getvalue() == "gaussian") ||
          (family.getvalue() == "cumprobit") ||
          (family.getvalue() == "binomialprobit")  ||
          (family.getvalue() == "bernoullilogit") ||
          (family.getvalue() == "binomialtlink") ||
          (family.getvalue() == "multinomialprobit")
         )

        {

        fcnonpgaussian.push_back( FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
                                        distr[distr.size()-1],
                                        D.getCol(j1),D.getCol(j2),
                                        lambda,maxdist,mapname,title,
                                        pathnonp, pathres,pathmap,collinpred
                                        )
                           );

        fcnonpgaussian[fcnonpgaussian.size()-1].init_name("regionnr");

        if (collinpred==0)
          {
          pathnonp = defaultpath + "\\temp\\" + name + add_name + "_f_" +
                       terms[i].varnames[0] + "_" + terms[i].varnames[1] +
                       "_spatial_var.raw";

          pathres =  outfile.getvalue() + add_name + "_f_" +
                       terms[i].varnames[0] + "_" + terms[i].varnames[1] +
                       "_spatial_var.res";

          title = "f_" + terms[i].varnames[0] + "_" + terms[i].varnames[1] +
                  "_variance"+add_name;
          }
        else
          {
          pathnonp = defaultpath + "\\temp\\" + name + add_name + "_f_" +
                     ST::inttostring(collinpred+1) + "_" +
                     terms[i].varnames[0] + "_" + terms[i].varnames[1] +
                     "_spatial_var.raw";

          pathres =  outfile.getvalue() + add_name + "_f_" +
                     ST::inttostring(collinpred+1) + "_" +
                     terms[i].varnames[0] + "_" + terms[i].varnames[1] +
                     "_spatial_var.res";

          title = "f_" + ST::inttostring(collinpred+1) + "_" +
                  terms[i].varnames[0] + "_" + terms[i].varnames[1] +
                  "_variance"+add_name;
          }

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                       &fcnonpgaussian[fcnonpgaussian.size()-1],
                                       distr[distr.size()-1],a1,b1,title,
                                       pathnonp,pathres,false,collinpred)
                           );

        if (constlambda.getvalue() == true)
          fcvarnonp[fcvarnonp.size()-1].set_constlambda();



        fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

        }
      else
        {


        }


      }   // end: if ( nonpspatial.checkvector(terms,i) == true )

    } //  end:  for(i=0;i<terms.size();i++)

  */

  return false;

  }


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
      if(maxdist<=0) // wähle maxdist so, dass Korrelation für Punkte mitmaximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;//4.605170186;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;//6.638352068;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;//8.022007057;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;//9.158140446;
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
      if(maxdist<=0) // wähle maxdist so, dass Korrelation für Punkte mitmaximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;//4.605170186;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;//6.638352068;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;//8.022007057;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;//9.158140446;
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
      if(maxdist<=0) // wähle maxdist so, dass Korrelation für Punkte mitmaximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;//4.605170186;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;//6.638352068;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;//8.022007057;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;//9.158140446;
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
      if(maxdist<=0) // wähle maxdist so, dass Korrelation für Punkte mitmaximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;//4.605170186;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;//6.638352068;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;//8.022007057;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;//9.158140446;
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
        outerror("ERROR: map object doesn´t contain centroids\n");
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
      if(maxdist<=0) // wähle maxdist so, dass Korrelation für Punkte mitmaximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;//4.605170186;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;//6.638352068;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;//8.022007057;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;//9.158140446;
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
        outerror("ERROR: map object doesn´t contain centroids\n");
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

      // --------------- reading options, term information ---------------------
      MCMC::fieldtype type;
      if (terms[i].options[0] == "baseline")
        type = MCMC::RW2;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

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

      // -------------end: reading options, term information -------------------


      //--------- creating path for samples and and results, creating title ----

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_logbaseline.raw","_logbaseline.res","_logbaseline");

      //----- end: creating path for samples and and results, creating title ---

      fcbaseline.push_back( baseline_reml(&generaloptions,
                                              D.getCol(j),
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

      // -------------end: reading options, term information -------------------


      //--------- creating path for samples and and results, creating title ----


      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_logbaseline.raw","_logbaseline.res","_logbaseline");


      //----- end: creating path for samples and and results, creating title ---


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

      MCMC::fieldtype type = MCMC::mrflinear;

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
        outerror("ERROR: map object doesn´t contain centroids\n");
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
        outerror("ERROR: map object doesn´t contain centroids\n");
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

      // --------------- reading options, term information ---------------------
      MCMC::fieldtype type;
      if ( (terms[i].options[0] == "psplinerw1") ||
           (terms[i].options[0] == "tpsplinerw1") ||
           (terms[i].options[0] == "psplinerw1vrw1") ||
           (terms[i].options[0] == "psplinerw1vrw2") )
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

      // -------------end: reading options, term information -------------------


      //--------- creating path for samples and and results, creating title ----


      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline.raw","_pspline.res","_pspline");


      //----- end: creating path for samples and and results, creating title ---

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

      // -------------end: reading options, term information -------------------


      //--------- creating path for samples and and results, creating title ----


      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_pspline.raw","_pspline.res","_pspline");


      //----- end: creating path for samples and and results, creating title ---


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
    if (b.ylimbottom.changed() == true)
      ot = ot + "ylimbottom="+b.ylimbottom.getValueAsString() + " ";
    if (b.ylimtop.changed() == true)
      ot = ot + "ylimtop="+b.ylimtop.getValueAsString() + " ";
    if (b.ystep.changed() == true)
      ot = ot + "ystep="+b.ystep.getValueAsString() + " ";
    if (b.xlimbottom.changed() == true)
      ot = ot + "xlimbottom="+b.xlimbottom.getValueAsString() + " ";
    if (b.xlimtop.changed() == true)
      ot = ot + "xlimtop="+b.xlimtop.getValueAsString() + " ";
    if (b.xstep.changed() == true)
      ot = ot + "xstep="+b.xstep.getValueAsString() + " ";

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

// Multinomiale Modelle
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
        b.family.getvalue()=="cumprobit")
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

    } // end: if (!failure)

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
              + ", title = \"Effect of " + xvar +"\" xlab = " + xvar
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
            + ", title = \"Effect of " + xvar +"\" xlab = " + xvar
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
    else if(b.family.getvalue()=="cumlogit" || b.family.getvalue()=="cumprobit")
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

      // -------------- reading options, term information ----------------------
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

      // -------------- reading options, term information ----------------------

      // -------- creating paths for samples and results, titles ---------------

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_rw.raw","_rw.res","_rw");

      // -------- end: creating paths for samples and results, titles ----------

      fcnonpgaussian.push_back(FULLCOND_nonp_gaussian(&generaloptions,D.getCol(j),
                         unsigned(maxint.getvalue()),type,
                         title,pathres,lambda,startlambda));

      fcnonpgaussian[fcnonpgaussian.size()-1].init_name(terms[i].varnames[0]);
      fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);

      } // end: if ( nonprw1rw2.checkvector(terms,i) == true )

    }

  return false;
  }

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

      // -------------- reading options, term information ----------------------
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

      // -------------- reading options, term information ----------------------

      // -------- creating paths for samples and results, titles ---------------

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],"_rw.raw","_rw.res","_rw");

      // -------- end: creating paths for samples and results, titles ----------

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
      } // end: if ( nonprw1rw2.checkvector(terms,i) == true )
    }
  return false;
  }

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

  void remlreg::describe(optionlist & globaloptions)
  {
  statobject::describe(globaloptions);
  }


#if defined(BORLAND_OUTPUT_WINDOW)
//------------------------------------------------------------------------------
#pragma package(smart_init)
#endif














