//---------------------------------------------------------------------------
#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include<StatwinFrame.h>
#include<statwin_haupt.h>
#endif


#include<stepwisereg.h>
#include<fullcond.h>
#include<typeinfo.h>
#include<stddef.h>


//------------------------------------------------------------------------------
//------------- CLASS stepwisereg: implementation of member functions -------------
//------------------------------------------------------------------------------

bool stepwisereg::check_gaussian(void)
  {

  if ( (family.getvalue() == "gaussian") ||
       (family.getvalue() == "multgaussian") ||
     (family.getvalue() == "binomialprobit") ||
     (family.getvalue() == "bernoullilogit") ||
     (family.getvalue() == "binomialtlink") ||
     (family.getvalue() == "multinomialprobit") ||
     (family.getvalue() == "cumprobit")
     )
     return true;
  else
    return false;
  }


bool  stepwisereg::check_nongaussian(void)
  {
  if ( (family.getvalue() == "binomial") || (family.getvalue() == "poisson") ||
       (family.getvalue() == "gamma") || (family.getvalue() == "nbinomial") ||
       (family.getvalue() == "multinomial") || (family.getvalue() == "cox") )
     return true;
  else
    return false;
  }



void stepwisereg::make_paths(unsigned  collinpred,ST::string & pathnonp,
                          ST::string & pathres,
                          ST::string & title,ST::string  varname1,
                          ST::string varname2,
                          ST::string  endingraw,ST::string endingres,
                          ST::string  endingtitle)
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



void stepwisereg::create(void)
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
  nonprw1rw2 = term_autoreg_stepwise();
  nonpseason = term_season_stepwise();
  nonppspline = term_pspline_stepwise();
  nonpspatial = term_spatial_stepwise();
  randomeff = term_random_stepwise();
  randomeffslope = term_randomslope_stepwise();
  termfactor = term_factor_stepwise();
  termnonlinearf = term_nonlinearf_stepwise();

  termtypes.push_back(&offset);
  termtypes.push_back(&fixedeffects);
  termtypes.push_back(&nonprw1rw2);
  termtypes.push_back(&nonpseason);
  termtypes.push_back(&nonppspline);
  termtypes.push_back(&nonpspatial);
  termtypes.push_back(&randomeff);
  termtypes.push_back(&randomeffslope);
  termtypes.push_back(&termfactor);
  termtypes.push_back(&termnonlinearf);

  modreg = modelterm(&termtypes);

  udata = use();

  //------------STEPWISE ---------------------------------------
  vector<ST::string> cr;
  cr.push_back("AIC");
  cr.push_back("AIC_imp");
  cr.push_back("GCV");
  cr.push_back("BIC");

  criterion = stroption("criterion",cr,"GCV");

  steps = intoption("steps",1000,1,10000);

  vector<ST::string> tr;
  tr.push_back("trace_on");
  tr.push_back("trace_off");
  trace = stroption("trace",tr,"trace_on");

  number = intoption("number",20,1,50);

  vector<ST::string> stmodel;
  stmodel.push_back("empty");
  stmodel.push_back("full");
  stmodel.push_back("both");
  stmodel.push_back("userdefined");
  startmodel = stroption("startmodel",stmodel,"empty");

  increment = intoption("increment",1,1,5);

  fine_tuning = simpleoption("fine_tuning",false);

  //----------END: STEPWISE --------------------------------------

  vector<ST::string> scalega;
  scalega.push_back("fixed");
  scalega.push_back("phi");
  scalega.push_back("random");
  scalegamma = stroption("scalegamma",scalega,"phi");
  gamvar = doubleoption("gammavar",0.001,0,1000);
  cit = intoption("cit",500,0,10000000);
  scalevalue = doubleoption("scale",1,0,1000000);
  constscale = simpleoption("constscale",false);

  knotsdef.push_back("equidistant");
  knotsdef.push_back("quantiles");
  knots = stroption("knots",knotsdef,"equidistant");

  maxint = intoption("maxint",150,0,20000);
  families.reserve(10);
  families.push_back("gaussian");
  families.push_back("binomial");
  families.push_back("poisson");
  families.push_back("gamma");
  families.push_back("nbinomial");
  family = stroption("family",families,"binomial");

  vector<ST::string> dop;
  dop.push_back("nb");
  dop.push_back("poga");
  dop.push_back("poig");
  distopt = stroption("distopt",dop,"nb");

  vector<ST::string> pro;
  pro.push_back("uniform");
  pro.push_back("gamma");
  propopt = stroption("propopt",pro,"uniform");

  propvar = doubleoption("propvar",0.1,0,500);

  predict = simpleoption("predict",false);
  predictmu = simpleoption("predictmu",false);
  predictuntil=intoption("predictuntil",0,1,100000000000);

  regressoptions.reserve(100);

  regressoptions.push_back(&maxint);
  regressoptions.push_back(&family);

  regressoptions.push_back(&gamvar);
  regressoptions.push_back(&cit);
  regressoptions.push_back(&scalevalue);
  regressoptions.push_back(&scalegamma);
  regressoptions.push_back(&constscale);

  regressoptions.push_back(&knots);

  regressoptions.push_back(&predict);
  regressoptions.push_back(&predictmu);
  regressoptions.push_back(&predictuntil);

  regressoptions.push_back(&propvar);
  regressoptions.push_back(&propopt);
  regressoptions.push_back(&distopt);

  regressoptions.push_back(&criterion);
  regressoptions.push_back(&steps);
  regressoptions.push_back(&trace);
  regressoptions.push_back(&number);
  regressoptions.push_back(&startmodel);
  regressoptions.push_back(&increment);
  regressoptions.push_back(&fine_tuning);

  // method 0

  methods.push_back(command("regress",&modreg,&regressoptions,&udata,required,
			 optional,optional,optional,optional,required));

  functions[0] = regressrun;


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
  title0 = stroption("title");
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
  plotnonpoptions.push_back(&title0);
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
  plotvar = stroption("plotvar","pmean");
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


void stepwisereg::initpointers(void)
  {

  unsigned i;

  for(i=0;i<distrstring.size();i++)
    {
    if (distrstring[i] == "gaussian")
    distr.push_back(&distr_gaussian[distrposition[i]]);
    else if (distrstring[i] == "binomial")
      distr.push_back(&distr_binomial);
    else if (distrstring[i] == "poisson")
      distr.push_back(&distr_poisson);
    else if (distrstring[i] == "gamma")
      distr.push_back(&distr_gamma);
    else if (distrstring[i] == "nbinomial")
      distr.push_back(&distr_nbinomial);
    }


  for(i=0;i<fcpsplinesurfgaussian.size();i++)
    fullcond.push_back(&fcpsplinesurfgaussian[i]);

  for(i=0;i<normalconst.size();i++)
    fullcond.push_back(&normalconst[i]);

  for(i=0;i<gammaconst.size();i++)
    fullcond.push_back(&gammaconst[i]);

  for(i=0;i<nongaussianconst.size();i++)
    fullcond.push_back(&nongaussianconst[i]);

  for(i=0;i<factorconst.size();i++)
    fullcond.push_back(&factorconst[i]);

  for(i=0;i<normalconst_special.size();i++)
    fullcond.push_back(&normalconst_special[i]);

  for(i=0;i<fcnonpgaussian.size();i++)
    fullcond.push_back(&fcnonpgaussian[i]);

  for(i=0;i<fcpspline.size();i++)
    fullcond.push_back(&fcpspline[i]);

  for(i=0;i<fciwlspspline.size();i++)
    fullcond.push_back(&fciwlspspline[i]);

  for(i=0;i<fcpsplinegaussian.size();i++)
    fullcond.push_back(&fcpsplinegaussian[i]);

  for(i=0;i<fcpsplinesurf.size();i++)
    fullcond.push_back(&fcpsplinesurf[i]);

  for(i=0;i<fcrandom.size();i++)
    fullcond.push_back(&fcrandom[i]);

  for(i=0;i<fcrandomgaussian.size();i++)
    fullcond.push_back(&fcrandomgaussian[i]);

  for(i=0;i<fcvarnonp.size();i++)
    fullcond.push_back(&fcvarnonp[i]);

  }


void stepwisereg::clear(void)
  {

  outfiles.erase(outfiles.begin(),outfiles.end());
  outfiles.reserve(10);

  generaloptions.erase(generaloptions.begin(),generaloptions.end());
  generaloptions.reserve(10);

  distrstring.erase(distrstring.begin(),distrstring.end());
  distrstring.reserve(10);

  distrposition.erase(distrposition.begin(),distrposition.end());
  distrposition.reserve(10);

  distr.erase(distr.begin(),distr.end());
  distr.reserve(10);

  distr_gaussian.erase(distr_gaussian.begin(),distr_gaussian.end());
  distr_gaussian.reserve(5);

  fullcond.erase(fullcond.begin(),fullcond.end());
  fullcond.reserve(50);

  normalconst.erase(normalconst.begin(),normalconst.end());
  normalconst.reserve(20);

  normalconst_special.erase(normalconst_special.begin(),
  normalconst_special.end());
  normalconst_special.reserve(20);

  factorconst.erase(factorconst.begin(),factorconst.end());
  factorconst.reserve(20);

  gammaconst.erase(gammaconst.begin(),gammaconst.end());
  gammaconst.reserve(10);

  nongaussianconst.erase(nongaussianconst.begin(),nongaussianconst.end());
  nongaussianconst.reserve(10);

  fcvarnonp.erase(fcvarnonp.begin(),fcvarnonp.end());
  fcvarnonp.reserve(40);

  fcnonpgaussian.erase(fcnonpgaussian.begin(),fcnonpgaussian.end());
  fcnonpgaussian.reserve(40);

  fcpspline.erase(fcpspline.begin(),fcpspline.end());
  fcpspline.reserve(40);

  fciwlspspline.erase(fciwlspspline.begin(),fciwlspspline.end());
  fciwlspspline.reserve(20);

  fcpsplinegaussian.erase(fcpsplinegaussian.begin(),fcpsplinegaussian.end());
  fcpsplinegaussian.reserve(40);

  fcpsplinesurf.erase(fcpsplinesurf.begin(),fcpsplinesurf.end());
  fcpsplinesurf.reserve(20);

  fcpsplinesurfgaussian.erase(fcpsplinesurfgaussian.begin(),
                              fcpsplinesurfgaussian.end());
  fcpsplinesurfgaussian.reserve(20);

  fcrandom.erase(fcrandom.begin(),fcrandom.end());
  fcrandom.reserve(40);

  fcrandomgaussian.erase(fcrandomgaussian.begin(),fcrandomgaussian.end());
  fcrandomgaussian.reserve(40);

  }


stepwisereg::stepwisereg(
#if defined(JAVA_OUTPUT_WINDOW)
administrator_basic * adb, administrator_pointer * adp,
#endif
const ST::string & n,ofstream * lo,istream * in,ST::string p,
vector<statobject*> * st)
: statobject(
#if defined(JAVA_OUTPUT_WINDOW)
adb,
#endif
n,"stepwisereg",lo,in,p)
  {
  statobj = st;
  create();
  resultsyesno = false;
  describetext.push_back("CURRENT REGRESSION RESULTS: none\n");
  }


stepwisereg::stepwisereg(const stepwisereg & b) : statobject(statobject(b))
  {
  create();
  #if defined(JAVA_OUTPUT_WINDOW)
  adminp_p = b.adminp_p;
  #endif
  statobj = b.statobj;
  D = b.D;
  distrstring = b.distrstring;
  modelvarnamesv = b.modelvarnamesv;
  runobj = b.runobj;
  distr_gaussian = b.distr_gaussian;
  distr_binomial = b.distr_binomial;
  distr_poisson = b.distr_poisson;
  distr_gamma = b.distr_gamma;
  distr_nbinomial = b.distr_nbinomial;
  terms = b.terms;
  normalconst = b.normalconst;
  nongaussianconst = b.nongaussianconst;
  fcpspline = b.fcpspline;
  resultsyesno = b.resultsyesno;
  initpointers();
  }


const stepwisereg & stepwisereg::operator=(const stepwisereg & b)
  {
  if (this == & b)
	 return *this;
  statobject::operator=(statobject(b));
  create();
  #if defined(JAVA_OUTPUT_WINDOW)
  adminp_p = b.adminp_p;
  #endif
  statobj = b.statobj;
  D = b.D;
  distrstring = b.distrstring;
  modelvarnamesv = b.modelvarnamesv;
  runobj = b.runobj;
  distr_gaussian = b.distr_gaussian;
  distr_binomial = b.distr_binomial;
  distr_poisson = b.distr_poisson;
  distr_gamma = b.distr_gamma;
  distr_nbinomial = b.distr_nbinomial;
  terms = b.terms;
  normalconst = b.normalconst;
  nongaussianconst = b.nongaussianconst;
  fcpspline = b.fcpspline;
  resultsyesno = b.resultsyesno;
  initpointers();
  return *this;
  }


int stepwisereg::parse(const ST::string & c)
  {

  int u = statobject::parse(c);  

  int pos = statobject::parsecom(c,methods,globaloptions);

  if (pos >= 0)
	 (*functions[pos])(*this);

  }


bool stepwisereg::create_generaloptions(void)
  {

  generaloptions.push_back(MCMCoptions(
  #if defined(JAVA_OUTPUT_WINDOW)
  adminb_p,
  #endif
  1000,100,10,logout,95,80));

  return false;

  }



bool stepwisereg::create_distribution(void)
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
  unsigned weightpos;

  modelvarnamesv = modreg.getModelVarnamesAsVector();
  rname = modelvarnamesv[0].to_bstr();
  wn = methods[0].get_weight_variable().to_bstr();
  if (wn.length() != 0)
    {
    modelvarnamesv.push_back(wn);
    weightpos = modelvarnamesv.size()-1;
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


  datamatrix offs(1,1,0);
  failure = create_offset(offs);


  datamatrix w;
  datamatrix pind;

  if (wn.length() > 0)
    {
    w = D.getCol(weightpos);
    }
  else
    w = datamatrix(1,1);


  ST::string path = outfile.getvalue() + add_name + "_predictmean.raw";
  ST::string pathfull = outfile.getvalue() + add_name + "_predictmu.raw";
  ST::string pathfullsample = defaultpath + "\\temp\\" + name + add_name +
                              "_predictmu.raw";
  ST::string pathdev = outfile.getvalue() + add_name + "_deviance_sample.raw";


//---------------------------- Gaussian response -------------------------------
  if (family.getvalue() == "gaussian")
    {

    ST::string path2 = outfile.getvalue() + add_name + "_scale.res";
    ST::string path3 = defaultpath + "\\temp\\" + name + add_name + "_scale.raw";

    if (offs.rows() == 1)
      distr_gaussian.push_back(DISTRIBUTION_gaussian(1.0,0.005,
      &generaloptions[generaloptions.size()-1],
      D.getCol(0),path2,path3,w));
    else
      distr_gaussian.push_back(DISTRIBUTION_gaussian(offs,1.0,0.005,
      &generaloptions[generaloptions.size()-1],
      D.getCol(0),path2,path3,w));


    distr_gaussian[distr_gaussian.size()-1].init_names(rname,wn);


    // prediction stuff

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_gaussian[distr_gaussian.size()-1].set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_gaussian[distr_gaussian.size()-1].set_predictfull(pathfullsample,pathfull,n);
      }

    if (pind.rows() > 1)
      distr_gaussian[distr_gaussian.size()-1].set_predictresponse(pind);

    // end: prediction stuff

    distr.push_back(&distr_gaussian[distr_gaussian.size()-1]);
    distrstring.push_back("gaussian");
    distrposition.push_back(distr_gaussian.size()-1);
    nrcategories = 1;
    }
//-------------------------- END: Gaussian response ----------------------------
//----------------------- binomial response, logit link ------------------------
  else if (family.getvalue() == "binomial")
    {

    if (offs.rows() == 1)
      distr_binomial = DISTRIBUTION_binomial(
      &generaloptions[generaloptions.size()-1],D.getCol(0),w);
    else
      distr_binomial = DISTRIBUTION_binomial(offs,
      &generaloptions[generaloptions.size()-1],D.getCol(0),w);

    if (distr_binomial.geterrors().size() > 0)
      {
      outerror(distr_binomial.geterrors());
      return true;
      }

    distr_binomial.init_names(rname,wn);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_binomial.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_binomial.set_predictfull(pathfullsample,pathfull,n);
      }

    distr.push_back(&distr_binomial);
    distrstring.push_back("binomial");
    distrposition.push_back(0);
    nrcategories = 1;
    }
//--------------------- END: binomial response, logit link ---------------------
//------------------------------ gamma response --------------------------------
  else if (family.getvalue() == "gamma")
    {

    ST::string path2 = outfile.getvalue() + add_name + "_scale.res";
    ST::string path3 = defaultpath + "\\temp\\" + name + add_name + "_scale.raw";

    int st = cit.getvalue();
    double v1 = gamvar.getvalue();

    double sv = scalevalue.getvalue();

    ST::string type = scalegamma.getvalue();

    if (type=="fixed")
      {
      distr_gamma = DISTRIBUTION_gamma(sv,&generaloptions[generaloptions.size()-1],D.getCol(0),path2,
                                       path3,w);
      }
    else if (type=="phi")
      {
      distr_gamma = DISTRIBUTION_gamma(1.0,0.005,st,
               &generaloptions[generaloptions.size()-1],D.getCol(0),path2,path3,
                                     w);
      }
    else
      {
      distr_gamma = DISTRIBUTION_gamma(1.0,0.005,v1,st,
      &generaloptions[generaloptions.size()-1],D.getCol(0),path2,path3,w);
      }

    if (offs.rows() > 1)
      {
      distr_gamma.init_offset(offs);
      }


    distr_gamma.init_names(rname,wn);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_gamma.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
        n = D.rows();
      distr_gamma.set_predictfull(pathfullsample,pathfull,n);
      }

    distr.push_back(&distr_gamma);
    distrstring.push_back("gamma");
    distrposition.push_back(0);
    nrcategories = 1;
    }
//----------------------------- END: gamma response ----------------------------
//--------------------------- poisson response ---------------------------------
  else if (family.getvalue() == "poisson")
    {

    if (offs.rows() == 1)
      distr_poisson = DISTRIBUTION_poisson(
      &generaloptions[generaloptions.size()-1],D.getCol(0),w);
    else
      distr_poisson = DISTRIBUTION_poisson(offs,
      &generaloptions[generaloptions.size()-1],D.getCol(0),w);
    distr_poisson.init_names(rname,wn);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_poisson.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_poisson.set_predictfull(pathfullsample,pathfull,n);
      }


    distr.push_back(&distr_poisson);
    distrstring.push_back("poisson");
    distrposition.push_back(0);
    nrcategories = 1;

    }
//------------------------ END: poisson response -------------------------------
//----------------------- negative binomial response ---------------------------
  else
    {

    ST::string path2 = outfile.getvalue() + add_name + "_scale.res";
    ST::string path3 = defaultpath + "\\temp\\" + name + add_name + "_scale.raw";

    MCMC::vertopt vo;

    if (distopt.getvalue() == "nb")
      vo = MCMC::nb;
    else if (distopt.getvalue() == "poga")
      vo = MCMC::poga;
    else
      vo = MCMC::poig;

    MCMC::propscale po;

    if (propopt.getvalue() == "uniform")
      po = MCMC::unif;
    else
      po = MCMC::gam;

    if (offs.rows() == 1)       // without offset
      distr_nbinomial = DISTRIBUTION_nbinomial(1.0,0.005,
      propvar.getvalue(),vo,po,hierarchical.getvalue(),
      &generaloptions[generaloptions.size()-1],D.getCol(0),path2,path3, w);
    else
      distr_nbinomial = DISTRIBUTION_nbinomial(1.0,0.005,
      propvar.getvalue(),vo,po,hierarchical.getvalue(),
      offs,&generaloptions[generaloptions.size()-1],D.getCol(0),path2,path3,w);
    distr_nbinomial.init_names(rname,wn);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_nbinomial.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_nbinomial.set_predictfull(pathfullsample,pathfull,n);
      }

    distr.push_back(&distr_nbinomial);
    distrstring.push_back("nbinomial");
    distrposition.push_back(0);
    nrcategories = 1;

    }
//--------------------- END: negative binomial response ------------------------


//----------------- end: reading data, creating designmatrices -----------------

  return false;

  }


bool stepwisereg::create_offset(datamatrix & o)
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

    }

  return false;

  }


bool stepwisereg::create_factor(const unsigned & collinpred)
  {

  int j;
  int f;
  double reference;
  double lambdastart;
  bool forced_into;
  unsigned i;

  for(i=0;i<terms.size();i++)
    {
    if ( termfactor.checkvector(terms,i) == true )
      {

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[2]).strtodouble(reference);

      f = (terms[i].options[3]).strtodouble(lambdastart);

      if (terms[i].options[4] == "true")
         forced_into = true;
      else
         forced_into = false;


      if (f==1)
        return true;

      ST::string title;
      ST::string pathconst;
      ST::string pathconstres;

      make_paths(collinpred,pathconst,pathres,title,terms[i].varnames[0],"",
                 "_factor.raw","_factor.res","");


      if ( check_gaussian())
        {

        factorconst.push_back(
        FULLCOND_const(&generaloptions[generaloptions.size()-1],
                                distr[distr.size()-1],D.getCol(j),
                                terms[i].options[1],
                                reference,title,pathconst,pathconstres,
                                collinpred));

        factorconst[factorconst.size()-1].init_name(terms[i].varnames[0]);

        factorconst[factorconst.size()-1].set_stepwise_options(
        lambdastart,0,-1,forced_into,0,0,false,false,0,false);

        factorconst[factorconst.size()-1].set_fcnumber(fullcond.size());

        fullcond.push_back(&factorconst[factorconst.size()-1]);

        }
      else
        {

        // NONGAUSSISAN fehlt

        }

      }

    }

   return false;

  }


bool stepwisereg::create_nonlinearf(const unsigned & collinpred)
  {

  int j;
  int f;
  double lambdastart;
  bool forced_into;
  unsigned i;

  for(i=0;i<terms.size();i++)
    {
    if ( termnonlinearf.checkvector(terms,i) == true )
      {

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[2]).strtodouble(lambdastart);

      if (terms[i].options[3] == "true")
         forced_into = true;
      else
         forced_into = false;


      if (f==1)
        return true;

      ST::string title;
      ST::string pathconst;
      ST::string pathconstres;

      make_paths(collinpred,pathconst,pathres,title,terms[i].varnames[0],"",
                 "_nonlinearf.raw","_nonlinearf.res","");


      if ( check_gaussian())
        {

        normalconst_special.push_back(
        FULLCOND_const_gaussian_special(&generaloptions[generaloptions.size()-1],
                                distr[distr.size()-1],D.getCol(j),
                                title,pathconst,pathconstres,
                                collinpred));

        normalconst_special[normalconst_special.size()-1].init_name(terms[i].varnames[0]);

        normalconst_special[normalconst_special.size()-1].set_stepwise_options(
        lambdastart,2,1,forced_into,0,0,false,false,0,false);

        normalconst_special[normalconst_special.size()-1].set_fcnumber(fullcond.size());

        fullcond.push_back(&normalconst_special[normalconst_special.size()-1]);

        }
      else
        {

        // NONGAUSSISAN fehlt

        }

      }

    }

   return false;

  }


bool stepwisereg::create_const(const unsigned & collinpred)
  {
  unsigned i;
  int j;
  vector<ST::string> varnames;
  vector<ST::string> varnamesh =  fixedeffects.get_constvariables(terms);

  varnames.push_back("const");

  for(i=0;i<varnamesh.size();i++)
    varnames.push_back(varnamesh[i]);

  unsigned nr = varnames.size();

  if (nr > 0)
    {
    vector< vector<ST::string> > varnamesvec;

    varnamesvec.push_back(varnames);

    for(unsigned k=0;k<varnamesvec.size();k++)
      {

      bool constincl = false;
      int constpos = -1;

      nr = varnamesvec[k].size();

      ST::string title;
      ST::string pathconst;
      ST::string pathconstres;

      if (collinpred == 0)
        {
        title = "FixedEffects" + ST::inttostring(k+1) + add_name;
        pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                           + add_name + "_FixedEffects" + ST::inttostring(k+1) + ".raw";

        pathconstres = outfile.getvalue() + add_name + "_FixedEffects" + ST::inttostring(k+1)
                       + ".res";

        }
      else
        {
        title = "FixedEffects" + ST::inttostring(k+1) + "_" +
                            ST::inttostring(collinpred+1) + add_name;
        pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                           + add_name + "_FixedEffects" + ST::inttostring(k+1) +
                           "_" + ST::inttostring(collinpred+1) + ".raw";

        pathconstres = outfile.getvalue() + add_name + "_FixedEffects" + ST::inttostring(k+1)
                        + "_" + ST::inttostring(collinpred+1) + ".res";

        }

      if (pathconst.isvalidfile() == 1)
        {
        errormessages.push_back("ERROR: unable to open file " + pathconst +
                                 " for writing\n");
        return true;
        }

      datamatrix X(D.rows(),nr,1);

      if (varnamesvec[k].size() != 0)
        {

        for(i=0;i<varnamesvec[k].size();i++)
          {

          if (varnamesvec[k][i] != "const")
            {

            j = varnamesvec[k][i].isinlist(modelvarnamesv);

            if (j != -1)
              {
              unsigned l;
              double * workX=X.getV()+i;
              double * workD=D.getV()+j;
              for (l=0;l<X.rows();l++,workX+=X.cols(),workD+=D.cols())
                *workX = *workD;
              }

            }
          else
            {
            constincl = true;
            constpos = 0;
            }

          }

        if ( check_gaussian())
          {

          normalconst.push_back(FULLCOND_const_gaussian(&generaloptions[generaloptions.size()-1],distr[distr.size()-1],X,
                                title,constpos,pathconst,pathconstres,
                                collinpred));
          normalconst[normalconst.size()-1].init_names(varnamesvec[k]);

          normalconst[normalconst.size()-1].set_fcnumber(fullcond.size());

          if (constincl == true)
            fcconst_intercept = &normalconst[normalconst.size()-1];
          fullcond.push_back(&normalconst[normalconst.size()-1]);

          } // end: if (family.getvalue() == "gaussian")
/*        else if (family.getvalue() == "gamma")
          {

          gammaconst.push_back(FULLCOND_const_gamma(&generaloptions[generaloptions.size()-1],
                                   distr[distr.size()-1],X,title,constpos,pathconst,
                                   pathconstres));
          gammaconst[gammaconst.size()-1].init_names(varnamesvec[k]);

          gammaconst[gammaconst.size()-1].set_fcnumber(fullcond.size());

          if (constincl == true)
            fcconst_intercept = &gammaconst[gammaconst.size()-1];
          fullcond.push_back(&gammaconst[gammaconst.size()-1]);

          }
*/
        else
          {

          nongaussianconst.push_back(FULLCOND_const_nongaussian(&generaloptions[generaloptions.size()-1],
                                   distr[distr.size()-1],X,title,constpos,pathconst,pathconstres,
                                   collinpred));
          nongaussianconst[nongaussianconst.size()-1].init_names(varnamesvec[k]);

          nongaussianconst[nongaussianconst.size()-1].set_fcnumber(fullcond.size());

          if (constincl == true)
            fcconst_intercept = &nongaussianconst[nongaussianconst.size()-1];
          fullcond.push_back(&nongaussianconst[nongaussianconst.size()-1]);


          }

        } // end: if (varnamesvec[k].size() != 0)

      }  // end:   for(k=0;k<varnamesvec.size();k++)

    } // end: if (nr > 0)

  return false;
  }

bool stepwisereg::create_nonprw1rw2(const unsigned & collinpred)
  {

  double hd;
  double lambda;
  int f;
  double lambdamin;
  double lambdamax;
  double lambdastart;
  bool forced_into;
  double df_lambdamax;
  double df_lambdamin;
  bool lambdamax_opt;
  bool lambdamin_opt;
  double numb;
  bool df_equidist;
  double df_accuracy;

  bool varcoeff;
  MCMC::fieldtype type;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonprw1rw2.checkvector(terms,i) == true )
      {

      // -------------- reading options, term information ----------------------

      if ( (terms[i].options[0] == "rw1") ||
           (terms[i].options[0] == "varcoeffrw1") )
        type = MCMC::RW1;
      else
        type = MCMC::RW2;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);

      if (terms[i].varnames.size()==1)
        {
        varcoeff=false;
//        if (terms[i].options[0] == "rw1")
//          type = MCMC::RW1;
//        else
//          type = MCMC::RW2;
        }
      else
        {
        varcoeff=true;
//        if (terms[i].options[0] == "rw1")
//          type = MCMC::RW1;
//        else
//          type = MCMC::RW2;

        j2 = terms[i].varnames[1].isinlist(modelvarnamesv);
        }




      f = (terms[i].options[1]).strtodouble(hd);
      lambda = hd;

      f = (terms[i].options[2]).strtodouble(lambdamin);
      f = (terms[i].options[3]).strtodouble(lambdamax);
      f = (terms[i].options[4]).strtodouble(lambdastart);

      if (terms[i].options[5] == "true")
         forced_into = true;
      else
         forced_into = false;

      f = (terms[i].options[6]).strtodouble(df_lambdamax);
      f = (terms[i].options[7]).strtodouble(df_lambdamin);

      if (terms[i].options[8] == "true")
         lambdamax_opt = true;
      else
         lambdamax_opt = false;

      if (terms[i].options[9] == "true")
         lambdamin_opt = true;
      else
         lambdamin_opt = false;

      f = (terms[i].options[10]).strtodouble(numb);
      if (terms[i].options[11] == "true")
         df_equidist = true;
      else
         df_equidist = false;
      f = (terms[i].options[12]).strtodouble(df_accuracy);


      if (f==1)
        return true;

      // -------------- reading options, term information ----------------------

      // -------- creating paths for samples and results, titles ---------------

      if (varcoeff==false)
        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_rw.raw","_rw.res","");
     else
        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                terms[i].varnames[0],
                 "_rw.raw","_rw.res","_rw");

      // -------- end: creating paths for samples and results, titles ----------

      if (varcoeff==false)
        {

        fcnonpgaussian.push_back(
        FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
        distr[distr.size()-1],D.getCol(j1),fcconst_intercept,
        unsigned(maxint.getvalue()),type,title,pathnonp,pathres,
        collinpred,lambda));

        fcnonpgaussian[fcnonpgaussian.size()-1].init_name(terms[i].varnames[0]);
        }
      else
        {

        fcnonpgaussian.push_back(
        FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
        distr[distr.size()-1],D.getCol(j2),D.getCol(j1),fcconst_intercept,
        unsigned(maxint.getvalue()),type,title,pathnonp,pathres,collinpred,
        lambda));

        vector<ST::string> na;
        na.push_back(terms[i].varnames[1]);
        na.push_back(terms[i].varnames[0]);
        fcnonpgaussian[fcnonpgaussian.size()-1].init_names(na);

        }

      fcnonpgaussian[fcnonpgaussian.size()-1].set_stepwise_options(
      lambdastart,lambdamax,lambdamin,forced_into,df_lambdamax,df_lambdamin,lambdamax_opt,lambdamin_opt,
      numb,df_equidist);
      fcnonpgaussian[fcnonpgaussian.size()-1].set_stepwise_accuracy(df_accuracy);


      if (check_nongaussian())
        fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS(1);

      fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());

      fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);


      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_rw_var.raw","_rw_var.res","_variance");

      fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                         &fcnonpgaussian[fcnonpgaussian.size()-1],
                                         distr[distr.size()-1],
                                         1.0,
                                         0.005,
                                         title,pathnonp,pathres,false,collinpred)
                             );


      fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

      } // end: if ( nonprw1rw2.checkvector(terms,i) == true )

    }

  return false;

  }


bool stepwisereg::create_pspline(const unsigned & collinpred)
  {

  long h;
  unsigned degree,nrknots;
  double lambda;
  int gridsize;
  int f;
  double lambdamin;
  double lambdamax;
  double lambdastart;
  bool forced_into;
  bool varcoeff;
  double df_lambdamax;
  double df_lambdamin;
  bool lambdamax_opt;
  bool lambdamin_opt;
  double numb;
  bool df_equidist;
  double df_accuracy;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonppspline.checkvector(terms,i) == true )
      {

      // --------------- reading options, term information ---------------------
      MCMC::fieldtype type;

      if ( (terms[i].options[0] == "psplinerw1") ||
           (terms[i].options[0] == "varpsplinerw1") )
        type = MCMC::RW1;
      else
        type = MCMC::RW2;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);

      if (terms[i].varnames.size()==1)
        {
        varcoeff=false;
        }
      else
        {
        varcoeff=true;
        j2 = terms[i].varnames[1].isinlist(modelvarnamesv);
        }

      f = (terms[i].options[1]).strtolong(h);
      degree = unsigned(h);

      f = (terms[i].options[2]).strtolong(h);
      nrknots = unsigned(h);

      f = (terms[i].options[3]).strtodouble(lambda);

      f = (terms[i].options[4]).strtolong(h);
      gridsize = unsigned(h);

      f = (terms[i].options[7]).strtodouble(lambdamin);
      f = (terms[i].options[8]).strtodouble(lambdamax);
      f = (terms[i].options[9]).strtodouble(lambdastart);

      if (terms[i].options[10] == "true")
         forced_into = true;
      else
         forced_into = false;

      f = (terms[i].options[11]).strtodouble(df_lambdamax);
      f = (terms[i].options[12]).strtodouble(df_lambdamin);

      if (terms[i].options[13] == "true")
         lambdamax_opt = true;
      else
         lambdamax_opt = false;

      if (terms[i].options[14] == "true")
         lambdamin_opt = true;
      else
         lambdamin_opt = false;

      f = (terms[i].options[15]).strtodouble(numb);

      if (terms[i].options[16] == "true")
         df_equidist = true;
      else
         df_equidist = false;
      f = (terms[i].options[17]).strtodouble(df_accuracy);

      MCMC::knotpos po;

      if (knots.getvalue() == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      if (f==1)
        return true;


      // -------------end: reading options, term information -------------------


      //--------- creating path for samples and and results, creating title ----

      if(varcoeff==false)
        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline.raw","_pspline.res","");
      else
        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],terms[i].varnames[0],
                 "_pspline.raw","_pspline.res","_pspline");

      //----- end: creating path for samples and and results, creating title ---

      // ---------------------- gaussian response, etc. ------------------------
      if (check_gaussian())
        {

        if(varcoeff==false)
          {
          fcpsplinegaussian.push_back(
          FULLCOND_pspline_gaussian(&generaloptions[generaloptions.size()-1],
                                              distr[distr.size()-1],
                                              fcconst_intercept,
                                              D.getCol(j1),
                                              nrknots,
                                              degree,
                                              po,
                                              type,
                                              "unrestricted",
                                              title,
                                              pathnonp,
                                              pathres,
                                              false,
                                              lambda,
                                              gridsize,
                                              false,
                                              collinpred
                                             )
                           );

          fcpsplinegaussian[fcpsplinegaussian.size()-1].init_name(terms[i].varnames[0]);
          }
        else
          {
          fcpsplinegaussian.push_back(
          FULLCOND_pspline_gaussian(&generaloptions[generaloptions.size()-1],
                                              distr[distr.size()-1],
                                              fcconst_intercept,
                                              D.getCol(j2),
                                              D.getCol(j1),
                                              nrknots,
                                              degree,
                                              po,
                                              type,
                                              "unrestricted",
                                              title,
                                              pathnonp,
                                              pathres,
                                              false,
                                              lambda,
                                              gridsize,
                                              collinpred
                                             )
                           );

          vector<ST::string> na;
          na.push_back(terms[i].varnames[1]);
          na.push_back(terms[i].varnames[0]);
          fcpsplinegaussian[fcpsplinegaussian.size()-1].init_names(na);
          }

        fcpsplinegaussian[fcpsplinegaussian.size()-1].set_stepwise_options(
               lambdastart,lambdamax,lambdamin,forced_into,df_lambdamax,df_lambdamin,lambdamax_opt,lambdamin_opt,
               numb,df_equidist);
        fcpsplinegaussian[fcpsplinegaussian.size()-1].set_stepwise_accuracy(df_accuracy);


        fcpsplinegaussian[fcpsplinegaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcpsplinegaussian[fcpsplinegaussian.size()-1]);

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                     "_pspline_var.raw","_pspline_var.res","_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(
                              &generaloptions[generaloptions.size()-1],
                              &fcpsplinegaussian[fcpsplinegaussian.size()-1],
                              distr[distr.size()-1],1.0,0.005,title,
                              pathnonp,pathres,false,collinpred)
                              );

        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());



        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);


        // ---------------------- end: gaussian response, etc. -----------------
        }
      else
        {
        // -------------------- non-gaussian response, etc. --------------------

        fciwlspspline.push_back( IWLS_pspline(&generaloptions[generaloptions.size()-1],
                                                distr[distr.size()-1],
                                                fcconst_intercept,
                                                D.getCol(j1),
                                                false,
                                                nrknots,
                                                degree,
                                                po,
                                                lambda,
                                                type,
                                                "unrestricted",
                                                1,
                                                false,
                                                2,
                                                1.0,0.005,
                                                title,
                                                pathnonp,
                                                pathres,
                                                false,
                                                gridsize,
                                                false,
                                                collinpred
                                               )
                             );

          fciwlspspline[fciwlspspline.size()-1].init_name(terms[i].varnames[0]);

          fciwlspspline[fciwlspspline.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fciwlspspline[fciwlspspline.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline_var.raw","_pspline_var.res","_variance");

          fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                &fciwlspspline[fciwlspspline.size()-1],
                                distr[distr.size()-1],1.0,0.005,title,
                                pathnonp,pathres,false,collinpred)
                                );

          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());


          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

        // ----------------- end:  non-gaussian response, etc. -----------------
        }

      }

    }

  return false;
  }


bool stepwisereg::create_spatial(const unsigned & collinpred)
  {

  ST::string pathnonpv;
  ST::string pathresv;

  double hd;
  int f;
  double lambda;
  unsigned i;
  int j1,j2;
  double lambdamin;
  double lambdamax;
  double lambdastart;
  bool forced_into;
  double df_lambdamax;
  double df_lambdamin;
  bool lambdamax_opt;
  bool lambdamin_opt;
  double numb;
  bool df_equidist;
  double df_accuracy;

  bool varcoeff;

  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial.checkvector(terms,i) == true )
      {

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);

      if (terms[i].varnames.size()==1)
        {
        varcoeff=false;
        }
      else
        {
        varcoeff=true;
        j2 = terms[i].varnames[1].isinlist(modelvarnamesv);
        }

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

      f = (terms[i].options[3]).strtodouble(lambdamin);
      f = (terms[i].options[4]).strtodouble(lambdamax);
      f = (terms[i].options[5]).strtodouble(lambdastart);

      if (terms[i].options[6] == "true")
         forced_into = true;
      else
         forced_into = false;

      f = (terms[i].options[7]).strtodouble(df_lambdamax);
      f = (terms[i].options[8]).strtodouble(df_lambdamin);

      if (terms[i].options[9] == "true")
         lambdamax_opt = true;
      else
         lambdamax_opt = false;

      if (terms[i].options[10] == "true")
         lambdamin_opt = true;
      else
         lambdamin_opt = false;

      f = (terms[i].options[11]).strtodouble(numb);

      if (terms[i].options[12] == "true")
         df_equidist = true;
      else
         df_equidist = false;
      f = (terms[i].options[13]).strtodouble(df_accuracy);


      ST::string start = startmodel.getvalue();
      if(lambdamax<lambdastart && start=="userdefined")
        {
        lambdastart = (lambdamax + lambdamin) / 2;                //nur bei startmodel=userdefined!!!
        outerror("ATTENTION: You forgot to specify a correct starting value for the smoothing parameter!");
        }


      if (f==1)
        return true;

      ST::string title, titlev;

      if (varcoeff==false)
        {
        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_spatial.raw","_spatial.res","");

        make_paths(collinpred,pathnonpv,pathresv,titlev,terms[i].varnames[0],"",
                   "_spatial_var.raw","_spatial_var.res","_variance");
        }
      else
        {

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                   terms[i].varnames[0],"_spatial.raw","_spatial.res","_spatial");

        make_paths(collinpred,pathnonpv,pathresv,titlev,terms[i].varnames[1],
                   terms[i].varnames[0],"_spatial_var.raw","_spatial_var.res",
                   "_spatial_variance");
        }

      if ( varcoeff==false)
        {
        fcnonpgaussian.push_back(
        FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
                                          distr[distr.size()-1],
                                          D.getCol(j1),
                                          fcconst_intercept,
                                          m,terms[i].options[1],
                                          title,
                                          pathnonp,
                                          pathres,collinpred,lambda
                                          )
                             );

        fcnonpgaussian[fcnonpgaussian.size()-1].init_name(terms[i].varnames[0]);
        }
      else
        {
        fcnonpgaussian.push_back(
        FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
        distr[distr.size()-1],fcconst_intercept,m,terms[i].options[1],
        D.getCol(j2),D.getCol(j1),title,pathnonp,pathres,collinpred,lambda));

        vector<ST::string> na;        
        na.push_back(terms[i].varnames[1]);
        na.push_back(terms[i].varnames[0]);
        fcnonpgaussian[fcnonpgaussian.size()-1].init_names(na);

        }

      fcnonpgaussian[fcnonpgaussian.size()-1].set_stepwise_options(
      lambdastart,lambdamax,lambdamin,forced_into,df_lambdamax,df_lambdamin,lambdamax_opt,lambdamin_opt,
      numb,df_equidist);
      fcnonpgaussian[fcnonpgaussian.size()-1].set_stepwise_accuracy(df_accuracy);


      if (fcnonpgaussian[fcnonpgaussian.size()-1].get_errors().size() > 0)
        {
        unsigned i;
        for(i=0;i<fcnonpgaussian[fcnonpgaussian.size()-1].get_errors().size();i++)
          errormessages.push_back(fcnonpgaussian[fcnonpgaussian.size()-1].get_errors()[i]);
        return true;
        }



      if (check_nongaussian())
        fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS(1);


      fcvarnonp.push_back(FULLCOND_variance_nonp(
      &generaloptions[generaloptions.size()-1],
      &fcnonpgaussian[fcnonpgaussian.size()-1],distr[distr.size()-1],1,0.005,
      titlev,pathnonpv,pathresv,false,collinpred));

      fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
      fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

      }   // end: if ( nonpspatial.checkvector(terms,i) == true )

    } //  end:  for(i=0;i<terms.size();i++)

  return false;
  }

bool stepwisereg::create_randomslope(const unsigned & collinpred)
  {


  ST::string pathnonp2;
  ST::string pathres2;
  ST::string pathresfixed;
  ST::string title2;

  unsigned i;
  int j1,j2;
  bool inclf;
  double lambda;
  double lambdamin,lambdamax,lambdastart;
  bool forced_into;
  double df_lambdamax;
  double df_lambdamin;
  bool lambdamax_opt;
  bool lambdamin_opt;
  double numb;
  bool df_equidist;
  double df_accuracy;

  int f;
  for(i=0;i<terms.size();i++)
    {
    if ( randomeffslope.checkvector(terms,i) == true )
      {

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);

      if (terms[i].options[1] == "true")
        inclf = false;
      else
        inclf = true;

      f = (terms[i].options[2]).strtodouble(lambda);

      f = (terms[i].options[3]).strtodouble(lambdamin);
      f = (terms[i].options[4]).strtodouble(lambdamax);
      f = (terms[i].options[5]).strtodouble(lambdastart);

      if (terms[i].options[6] == "true")
         forced_into = true;
      else
         forced_into = false;

      f = (terms[i].options[7]).strtodouble(df_lambdamax);
      f = (terms[i].options[8]).strtodouble(df_lambdamin);

      if (terms[i].options[9] == "true")
         lambdamax_opt = true;
      else
         lambdamax_opt = false;

      if (terms[i].options[10] == "true")
         lambdamin_opt = true;
      else
         lambdamin_opt = false;

      f = (terms[i].options[11]).strtodouble(numb);

      if (terms[i].options[12] == "true")
         df_equidist = true;
      else
         df_equidist = false;
      f = (terms[i].options[13]).strtodouble(df_accuracy);

      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_random.raw","_random.res","");

      make_paths(collinpred,pathnonp2,pathres2,title2,terms[i].varnames[1],
                 terms[i].varnames[0],"_random_var.raw",
                 "_random_var.res","_variance");

      if (check_gaussian())
        {

        fcrandomgaussian.push_back(
        FULLCOND_random_gaussian(&generaloptions[generaloptions.size()-1],
                                                        distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        D.getCol(j1),
                                                        D.getCol(j2),
                                                        title,
                                                        pathnonp,
                                                        pathres,pathresfixed,
                                                        lambda,
                                                        inclf,collinpred
                                                        )
                          );

        vector<ST::string> na;
        na.push_back(terms[i].varnames[1]);
        na.push_back(terms[i].varnames[0]);
        fcrandomgaussian[fcrandomgaussian.size()-1].init_names(na);
//        fcrandomgaussian[fcrandomgaussian.size()-1].init_name(terms[i].varnames[1]);

        fcrandomgaussian[fcrandomgaussian.size()-1].set_stepwise_options(
        lambdastart,lambdamax,lambdamin,forced_into,df_lambdamax,df_lambdamin,lambdamax_opt,lambdamin_opt,
        numb,df_equidist);
        fcrandomgaussian[fcrandomgaussian.size()-1].set_stepwise_accuracy(df_accuracy);


        fcrandomgaussian[fcrandomgaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcrandomgaussian[fcrandomgaussian.size()-1]);

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fcrandomgaussian[fcrandomgaussian.size()-1],
                            distr[distr.size()-1],1,0.005,title2,pathnonp2,pathres2,
                            false,collinpred)
                            );

        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);


        }

      else
        {

        fcrandom.push_back(
        FULLCOND_random_nongaussian(&generaloptions[generaloptions.size()-1],
                                                        distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        D.getCol(j1),
                                                        D.getCol(j2),
                                                        title,
                                                        pathnonp,
                                                        pathres,pathresfixed,
                                                        lambda,false,
                                                        collinpred
                                                        )
                          );

        fcrandom[fcrandom.size()-1].init_name(terms[i].varnames[1]);

        fcrandom[fcrandom.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcrandom[fcrandom.size()-1]);

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fcrandom[fcrandom.size()-1],
                            distr[distr.size()-1],1,0.005,
                            title2,pathnonp2,pathres2,false,collinpred)
                            );

        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

        }

      }

    }

  return false;

  }

bool stepwisereg::create_random(const unsigned & collinpred)
  {

  ST::string pathnonp2;
  ST::string pathres2;
  ST::string title2;
  double lambda;
  int f;
  double lambdamin,lambdamax,lambdastart;
  bool forced_into;
  double df_lambdamax;
  double df_lambdamin;
  bool lambdamax_opt;
  bool lambdamin_opt;
  double numb;
  bool df_equidist;
  double df_accuracy;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( randomeff.checkvector(terms,i) == true )
      {

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtodouble(lambda);

      f = (terms[i].options[2]).strtodouble(lambdamin);
      f = (terms[i].options[3]).strtodouble(lambdamax);
      f = (terms[i].options[4]).strtodouble(lambdastart);

      if (terms[i].options[5] == "true")
         forced_into = true;
      else
         forced_into = false;

      f = (terms[i].options[6]).strtodouble(df_lambdamax);
      f = (terms[i].options[7]).strtodouble(df_lambdamin);

      if (terms[i].options[8] == "true")
         lambdamax_opt = true;
      else
         lambdamax_opt = false;

      if (terms[i].options[9] == "true")
         lambdamin_opt = true;
      else
         lambdamin_opt = false;

      f = (terms[i].options[10]).strtodouble(numb);

      if (terms[i].options[11] == "true")
         df_equidist = true;
      else
         df_equidist = false;
      f = (terms[i].options[12]).strtodouble(df_accuracy);

      ST::string start = startmodel.getvalue();
      if(lambdamax<lambdastart && start=="userdefined")
        {
        lambdastart = (lambdamax + lambdamin) / 2;
        outerror("ATTENTION: You forgot to specify a correct starting value for the smoothing parameter!");
        }

      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_random.raw","_random.res","");

      make_paths(collinpred,pathnonp2,pathres2,title2,terms[i].varnames[0],"",
                 "_random_var.raw","_random_var.res","_variance");

      if (check_gaussian())
        {

        FULLCOND_nonp_gaussian * structuredp;
        unsigned structured=0;

        unsigned j1;
        for (j1=0;j1<fcnonpgaussian.size();j1++)
          {
          if  ( ((fcnonpgaussian[j1].get_datanames()).size() == 1) &&
                (fcnonpgaussian[j1].get_datanames()[0]==terms[i].varnames[0]) &&
                (fcnonpgaussian[j1].get_col() == collinpred) &&
                (fcnonpgaussian[j1].get_type() == MCMC::mrf)
              )
            {
            structuredp = &fcnonpgaussian[j1];
            structured ++;
            }
          }

        fcrandomgaussian.push_back(
        FULLCOND_random_gaussian(&generaloptions[generaloptions.size()-1],
                                                        distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        D.getCol(j),
                                                        title,
                                                        pathnonp,pathres,
                                                        lambda,collinpred
                                                        )
                          );

        fcrandomgaussian[fcrandomgaussian.size()-1].set_stepwise_options(
        lambdastart,lambdamax,lambdamin,forced_into,df_lambdamax,df_lambdamin,lambdamax_opt,lambdamin_opt,
        numb,df_equidist);
        fcrandomgaussian[fcrandomgaussian.size()-1].set_stepwise_accuracy(df_accuracy);

        if (structured==1)
          {

          ST::string pathnonpt = defaultpath + "\\temp\\" + name + add_name +
                 terms[i].varnames[0] +
                 "_spatialtotal.raw";

          ST::string pathrest =
          outfile.getvalue() + add_name + "_" + terms[i].varnames[0] +
                                "_spatialtotal.res";

          fcrandomgaussian[fcrandomgaussian.size()-1].init_spatialtotal(
          structuredp,pathnonpt,pathrest);
          }
        else if (structured==0)
          {
          }
        else
          {
          // FEHLERMELDUNG
          }

        fcrandomgaussian[fcrandomgaussian.size()-1].init_name(terms[i].varnames[0]);

        fcrandomgaussian[fcrandomgaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcrandomgaussian[fcrandomgaussian.size()-1]);

        fcvarnonp.push_back(
        FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fcrandomgaussian[fcrandomgaussian.size()-1],
                            distr[distr.size()-1],1.0,0.005,
                            title2,pathnonp2,pathres2,false,collinpred)
                            );

        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

        } // end: gaussian, probit, etc. ...

      else
        {

        FULLCOND_nonp_gaussian * structuredp;
        unsigned structured=0;

        unsigned j1;
        for (j1=0;j1<fcnonpgaussian.size();j1++)
          {
          if  ( ((fcnonpgaussian[j1].get_datanames()).size() == 1) &&
                (fcnonpgaussian[j1].get_datanames()[0]==terms[i].varnames[0]) &&
                (fcnonpgaussian[j1].get_col() == collinpred) //&&
                //(fcnonp[j1].get_type() == MCMC::mrf)
              )
            {
            structuredp = &fcnonpgaussian[j1];
            structured ++;
            }
          }

        fcrandom.push_back(
        FULLCOND_random_nongaussian(&generaloptions[generaloptions.size()-1],
                                                        distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        D.getCol(j),
                                                        title,
                                                        pathnonp,
                                                        pathres,
                                                        lambda,
                                                        collinpred
                                                        )
                          );

        fcrandomgaussian[fcrandomgaussian.size()-1].set_stepwise_options(
        lambdastart,lambdamax,lambdamin,forced_into,df_lambdamax,df_lambdamin,lambdamax_opt,lambdamin_opt,
        numb,df_equidist);
        fcrandomgaussian[fcrandomgaussian.size()-1].set_stepwise_accuracy(df_accuracy);

       if (structured==1)
          {

          ST::string pathnonpt = defaultpath + "\\temp\\" + name + add_name +
                 terms[i].varnames[0] +
                 "_spatialtotal.raw";

          ST::string pathrest = outfile.getvalue() + add_name + "_" + terms[i].varnames[0] +
                                "_spatialtotal.res";

//          fcrandom[fcrandom.size()-1].init_spatialtotal(
//          structuredp,pathnonpt,pathrest);
          }
        else if (structured==0)
          {
          }
        else
          {
          // FEHLERMELDUNG
          }


        fcrandom[fcrandom.size()-1].init_name(terms[i].varnames[0]);

        fcrandom[fcrandom.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcrandom[fcrandom.size()-1]);

        fcvarnonp.push_back(
        FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fcrandom[fcrandom.size()-1],
                            distr[distr.size()-1],1.0,0.005,title2,
                            pathnonp2,pathres2,false,collinpred)
                            );

        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

        }

      }

    }

  return false;

  }


void drawmaprun(stepwisereg & b)
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


void plotnonprun(stepwisereg & b)
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
    ot = ot + "title=\""+b.title0.getvalue() + "\" ";
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



void regressrun(stepwisereg & b)
  {

  b.resultsyesno = false;

  b.terms = b.modreg.getterms();

  b.describetext.erase(b.describetext.begin(),b.describetext.end());
  b.describetext.push_back("LAST ESTIMATED MODEL: \n");
  b.describetext.push_back("\n");
  b.describetext.push_back(b.modreg.getModelText());
  b.describetext.push_back("\n");

  b.clear();

  b.outfiles.push_back(b.outfile.getvalue());

  bool failure = false;

  if (!failure)
    failure = b.create_generaloptions();

  if (!failure)
    failure = b.create_distribution();



  if (!failure)
    {

    if (!failure)
      failure = b.create_const(0);

    if (!failure)
      failure = b.create_factor(0);

    if (!failure)
      failure = b.create_nonlinearf(0);

    if (!failure)
        failure = b.create_nonprw1rw2(0);

      if (!failure)
        failure = b.create_pspline(0);

      if (!failure)
        failure = b.create_nonpseason(0);

      if (!failure)
        failure = b.create_spatial(0);

      if (!failure)
        failure = b.create_random(0);

      if (!failure)
        failure = b.create_randomslope(0);

    } // end: if (!failure)


  if (!failure)
    {

    //vector<ST::string> header;
    //header.push_back("STEPWISEREG OBJECT " + b.name.to_bstr() +
    //                   ": stepwise procedure" );

    ST::string name = b.name.to_bstr();
    ST::string cr = b.criterion.getvalue();
    int steps = b.steps.getvalue();
    ST::string tr = b.trace.getvalue();
    int number = b.number.getvalue();
    ST::string stmodel = b.startmodel.getvalue();
    int increment = b.increment.getvalue();
    bool fine_tuning = b.fine_tuning.getvalue();
    vector<FULLCOND*> fullcond_z;

    b.runobj = STEPWISErun(&b.generaloptions[0],b.distr[0],b.fullcond);

    ST::string path = b.outfiles[0];

    failure = b.runobj.stepwise(cr,steps,tr,number,stmodel,increment,
    fine_tuning,b.D,b.modelvarnamesv,name,fullcond_z,path);

    if(!failure)
      {
      /*
      vector<ST::string> path;
      vector<ST::string> path2;
      vector<ST::string> path3;
      vector<ST::string> path4;

      unsigned i;
      for(i=0;i<b.outfiles.size();i++)
         {
         path.push_back(b.outfiles[i] + "_graphics.prg");
         path2.push_back(b.outfiles[i] + "_model_summary.tex");
         path3.push_back(b.outfiles[i] +  "_splus.txt");
         path4.push_back(b.outfiles[i] +  "_stata.do");
         }

      b.runobj.make_graphics(header,path,path2,path3,path4);
      */

#if defined(JAVA_OUTPUT_WINDOW)

      b.fullcond = fullcond_z;
      for(unsigned j=0;j<b.fullcond.size();j++)
         {
         MCMC::plotstyles plst = b.fullcond[j]->get_plotstyle();
         if(plst != MCMC::noplot)
           {
           vector<ST::string> varnames = b.fullcond[j]->get_datanames();
           ST::string xvar = varnames[0];
           ST::string pathresult = b.fullcond[j]->get_pathresult();
           ST::string pathps = pathresult.substr(0, pathresult.length()-4);
           if(plst == MCMC::plotnonp)
                   {
                   b.newcommands.push_back(b.name + ".plotnonp " + ST::inttostring(j)
                   + ", title = \"Effect of " + xvar +"\" xlab = " + xvar
                   + " ylab = \" \" outfile = " + pathps + ".ps replace");
                   }

           if(plst==MCMC::drawmap)
                   {
                  // double u = b.fullcond[j]->get_level1();
                  // double o = b.fullcond[j]->get_level2();
                  // ST::string u_str = ST::doubletostring(u,0);
                  // ST::string o_str = ST::doubletostring(o,0);
                   b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
                   + ", color outfile = " + pathps + "_pmean.ps replace");
                  // b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
                  // + ", plotvar = pcat" + u_str + " nolegend  pcat outfile = " + pathps
                  // + "_pcat" + u_str + ".ps replace");
                  // b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
                  // + ", plotvar = pcat" + o_str + " nolegend  pcat outfile = " + pathps
                  // + "_pcat" + o_str + ".ps replace");

                   }
           }
         }

#endif
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

bool stepwisereg::create_nonpseason(const unsigned & collinpred)
  {

  ST::string pathnonpv;
  ST::string pathresv;

  long h;
  double hd;
  unsigned per;
  double lambda;
  int f;
  double lambdamin;
  double lambdamax;
  double lambdastart;
  bool forced_into;
  double df_lambdamax;
  double df_lambdamin;
  bool lambdamax_opt;
  bool lambdamin_opt;
  double numb;
  bool df_equidist;
  double df_accuracy;

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

      f = (terms[i].options[3]).strtodouble(lambdamin);
      f = (terms[i].options[4]).strtodouble(lambdamax);
      f = (terms[i].options[5]).strtodouble(lambdastart);

      if (terms[i].options[6] == "true")
         forced_into = true;
      else
         forced_into = false;

      f = (terms[i].options[7]).strtodouble(df_lambdamax);
      f = (terms[i].options[8]).strtodouble(df_lambdamin);

      if (terms[i].options[9] == "true")
         lambdamax_opt = true;
      else
         lambdamax_opt = false;

      if (terms[i].options[10] == "true")
         lambdamin_opt = true;
      else
         lambdamin_opt = false;

      f = (terms[i].options[11]).strtodouble(numb);

      if (terms[i].options[12] == "true")
         df_equidist = true;
      else
         df_equidist = false;
      f = (terms[i].options[13]).strtodouble(df_accuracy);

      if (f==1)
        return true;

      ST::string start = startmodel.getvalue();
      if(lambdamax<lambdastart && start=="userdefined")
        {
        lambdastart = (lambdamax + lambdamin) / 2;
        outerror("ATTENTION: You forgot to specify a correct starting value for the smoothing parameter!");
        }

      ST::string titlev;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_rw_season.raw","_rw_season.res","");

      make_paths(collinpred,pathnonpv,pathresv,titlev,terms[i].varnames[0],"",
                 "_rw_season.raw","_season_var.res","_variance");


      fcnonpgaussian.push_back(
      FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
                                         distr[distr.size()-1],
                                         D.getCol(j),
                                         fcconst_intercept,
                                         unsigned(maxint.getvalue()),
                                         MCMC::seasonal,
                                         title,
                                         pathnonp,
                                         pathres,
                                         collinpred,lambda,per
                                        )
                        );

      fcnonpgaussian[fcnonpgaussian.size()-1].init_name(terms[i].varnames[0]);

      fcnonpgaussian[fcnonpgaussian.size()-1].set_stepwise_options(
      lambdastart,lambdamax,lambdamin,forced_into,df_lambdamax,df_lambdamin,lambdamax_opt,lambdamin_opt,
      numb,df_equidist);
      fcnonpgaussian[fcnonpgaussian.size()-1].set_stepwise_accuracy(df_accuracy);


      if (check_nongaussian())
        fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS(1);


      fcvarnonp.push_back(
      FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
      &fcnonpgaussian[fcnonpgaussian.size()-1],distr[distr.size()-1],
                                       1.0,0.005,titlev,pathnonpv,pathresv,
                                       false,collinpred));


      fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
      fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

      }

    }

  return false;
  }






void stepwisereg::describe(optionlist & globaloptions)
  {
  statobject::describe(globaloptions);
  }


#if defined(BORLAND_OUTPUT_WINDOW)
//------------------------------------------------------------------------------
#pragma package(smart_init)
#endif













