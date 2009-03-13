



#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include<StatwinFrame.h>
#include<statwin_haupt.h>

#endif

#include"superbayesreg.h"
#include <mapobject.h>

// Vorschlag:
//#include<typeinfo.h>

#include<stddef.h>


//------------------------------------------------------------------------------
//------------- CLASS superbayesreg: implementation of member functions --------
//------------------------------------------------------------------------------


void superbayesreg::make_paths(ST::string & pathnonp,
                          ST::string & pathres,
                          ST::string & title,vector<ST::string> vn,
                          ST::string  endingraw,ST::string endingres,
                          ST::string  endingtitle)
  {

  unsigned modnr = equations.size()-1;

  ST::string h = equations[modnr].paths;

  ST::string varname1 = vn[0];
  ST::string varname2;

  if (vn.size() >=2)
    varname2 = vn[1];
  else
    varname2 = "";

  if (varname2=="")
    {
    pathnonp = defaultpath + "\\temp\\" + name + "_" + h + "_f_"
                                  + varname1 + endingraw;

    pathres = outfile.getvalue() + "_" + h + "_" +
                     endingres + "_" + varname1 + ".res";

    title = h + ": " + endingtitle + varname1;

    }
  else
    {
    pathnonp = defaultpath + "\\temp\\" + name + "_" + h + "_" + varname2
                 + "_f_"  + varname1 + endingraw;

    pathres = outfile.getvalue() + "_" + h + "_" +
                     endingres + "_" + varname1 + "_" + varname2 + ".res";


    title = h + ": " + endingtitle + varname1 + " and " + varname2;

    }

  }

// forward declaration
void hregressrun(superbayesreg & b);

void superbayesreg::create_hregress(void)
  {

  // SYNTAX OF COMMANDS:
  // name [model] [weight varname] [by varname] [if expression]
  //      [, options] [using usingtext]

  // --------------------------- method hregress -------------------------------

  vector<ST::string> tnames;
  tnames.push_back("pspline");
  tnames.push_back("hrandom");
  tnames.push_back("spatial");
  tnames.push_back("hrandom_pspline");
  tnames.push_back("hrandomexp_pspline");


  tnonp = term_nonp(tnames);
  lineareffects = basic_termtype();
  termtypes.push_back(&tnonp);
  termtypes.push_back(&lineareffects);

  modreg = modelterm(&termtypes);

  udata = use();

  modeonly = simpleoption("modeonly",false);
  setseed = intoption("setseed",-1,0,MAXINT);

  iterations = intoption("iterations",52000,1,10000000);
  burnin = intoption("burnin",2000,0,500000);
  step = intoption("step",50,1,1000);
  level1 = doubleoption("level1",95,40,99);
  level2 = doubleoption("level2",80,40,99);
  storesample = simpleoption("storesample",false);

  families.reserve(20);
  families.push_back("gaussian");
  families.push_back("loggaussian");  
  families.push_back("gaussian_re");
  families.push_back("gaussian_exp");
  families.push_back("gaussian_mult");
  families.push_back("binomial_logit");
  family = stroption("family",families,"gaussian");
  aresp = doubleoption("aresp",0.001,-1.0,500);
  bresp = doubleoption("bresp",0.001,0.0,500);

  hlevel = intoption("hlevel",1,1,3);
  equationnr = intoption("equation",1,1,50);
  equationtypes.reserve(20);
  equationtypes.push_back("mean");
  equationtypes.push_back("variance");
  equationtype = stroption("equationtype",equationtypes,"mean");

  predictop.reserve(20);
  predictop.push_back("no");
  predictop.push_back("nosamples");
  predictop.push_back("full");
  predict = stroption("predict",predictop,"no");

  MSEop.reserve(20);
  MSEop.push_back("no")
  MSEop.push_back("all")
  MSEop.push_back("weightszero")
  mse = stroption("MSE",MSEop,"no");


  regressoptions.reserve(100);

  regressoptions.push_back(&modeonly);
  regressoptions.push_back(&setseed);
  regressoptions.push_back(&iterations);
  regressoptions.push_back(&burnin);
  regressoptions.push_back(&step);
  regressoptions.push_back(&level1);
  regressoptions.push_back(&level2);
  regressoptions.push_back(&storesample);
  regressoptions.push_back(&family);
  regressoptions.push_back(&aresp);
  regressoptions.push_back(&bresp);
  regressoptions.push_back(&hlevel);
  regressoptions.push_back(&equationnr);
  regressoptions.push_back(&equationtype);
  regressoptions.push_back(&predict);
  regressoptions.push_back(&mse);  


  // methods 0
  methods.push_back(command("hregress",&modreg,&regressoptions,&udata,required,
			 optional,optional,optional,optional,required));

  functions[0] = hregressrun;

  }

// forward declaration
void autocorrrun(superbayesreg & b);

void superbayesreg::create_autocorr(void)
  {

  // --------------------------- method autocor --------------------------------

  maxlag = intoption("maxlag",250,1,500);

  autocorroptions.push_back(&maxlag);

  // methods 1
  methods.push_back(command("autocor",&ma,&autocorroptions,&ad,notallowed,
						  notallowed,notallowed,notallowed,optional,
                          notallowed));

  functions[1] = autocorrrun;

  }

// forward declaration
void getsamplerun(superbayesreg & b);

void superbayesreg::create_getsample(void)
  {

  // methods 2
  methods.push_back(command("getsample",&mgetsample,&getsampleoptions,
						  &usegetsample,notallowed,notallowed,notallowed,
                          notallowed,notallowed,notallowed));

  functions[2] = getsamplerun;

  }

void superbayesreg::create(void)
  {

  generaloptions_yes = false;
  run_yes=false;

  ST::string h = defaultpath+"\\output\\"+name;

  outfile = fileoption("outfile",h,false);

  globaloptions.push_back(&outfile);

  create_hregress();

  create_autocorr();

  create_getsample();

  }

/*
void superbayesreg::initpointers(void)
  {

  unsigned i;

  for(i=0;i<FC_nonps.size();i++)
    {
    equations[0].FCpointer.push_back(&FC_nonps[i]);
    equations[0].FCpointer.push_back(&FC_nonp_variances[i]);
    }

  }
*/

void superbayesreg::clear(void)
  {

  equations.erase(equations.begin(),equations.end());
  equations.reserve(10);

  distr_gaussians.erase(distr_gaussians.begin(),distr_gaussians.end());
  distr_gaussians.reserve(10);

  distr_loggaussians.erase(distr_loggaussians.begin(),distr_loggaussians.end());
  distr_loggaussians.reserve(10);

  distr_gaussian_res.erase(distr_gaussian_res.begin(),distr_gaussian_res.end());
  distr_gaussian_res.reserve(10);

  distr_gaussian_exps.erase(distr_gaussian_exps.begin(),distr_gaussian_exps.end());
  distr_gaussian_exps.reserve(10);

  distr_gaussian_mults.erase(distr_gaussian_mults.begin(),
  distr_gaussian_mults.end());
  distr_gaussian_mults.reserve(10);

  distr_binomials.erase(distr_binomials.begin(),distr_binomials.end());
  distr_binomials.reserve(10);

  FC_linears.erase(FC_linears.begin(),FC_linears.end());
  FC_linears.reserve(30);

  design_psplines.erase(design_psplines.begin(),design_psplines.end());
  design_psplines.reserve(30);

  design_mrfs.erase(design_mrfs.begin(),design_mrfs.end());
  design_mrfs.reserve(30);

  design_hrandoms.erase(design_hrandoms.begin(),design_hrandoms.end());
  design_hrandoms.reserve(30);

  FC_nonps.erase(FC_nonps.begin(),FC_nonps.end());
  FC_nonps.reserve(30);

  FC_hrandoms.erase(FC_hrandoms.begin(),FC_hrandoms.end());
  FC_hrandoms.reserve(30);

  FC_mults.erase(FC_mults.begin(),FC_mults.end());
  FC_mults.reserve(30);

  FC_nonp_variances.erase(FC_nonp_variances.begin(),FC_nonp_variances.end());
  FC_nonp_variances.reserve(30);

  FC_hrandom_variances.erase(FC_hrandom_variances.begin(),
  FC_hrandom_variances.end());
  FC_hrandom_variances.reserve(30);

  FC_predicts.erase(FC_predicts.begin(),FC_predicts.end());
  FC_predicts.reserve(10);

  }


#if defined(JAVA_OUTPUT_WINDOW)
superbayesreg::superbayesreg(
administrator_basic * adb, administrator_pointer * adp,
const ST::string & n,ofstream * lo,istream * in,
						 ST::string p,vector<statobject*> * st)
						 : statobject(adb,n,"mcmcreg",lo,in,p)
  {
  clear();
  adminp_p = adp;
  statobj = st;
  master = MASTER_OBJ();
  create();
  resultsyesno = false;
  run_yes = false;
  posteriormode = false;
  describetext.push_back("CURRENT REGRESSION RESULTS: none\n");
  }
#else
superbayesreg::superbayesreg(const ST::string & n,ofstream * lo,istream * in,
						 ST::string p,vector<statobject*> * st)
						 : statobject(n,"mcmcreg",lo,in,p)
  {
  clear();
  statobj = st;
  master =MASTER_OBJ();
  create();
  resultsyesno = false;
  run_yes = false;
  posteriormode = false;
  describetext.push_back("CURRENT REGRESSION RESULTS: none\n");
  }
#endif


superbayesreg::superbayesreg(const superbayesreg & b) : statobject(statobject(b))
  {
  create();
  #if defined(JAVA_OUTPUT_WINDOW)
  adminp_p = b.adminp_p;
  #endif

  pathres = b.pathres;
  title = b.title;
  pathnonp = b.pathnonp;

  statobj = b.statobj;

  D = b.D;

  modelvarnamesv = b.modelvarnamesv;

  equations=b.equations;
  simobj = b.simobj;

  master = b.master;

  generaloptions = b.generaloptions;

  distr_gaussians = b.distr_gaussians;
  distr_loggaussians = b.distr_loggaussians;  
  distr_gaussian_res = b.distr_gaussian_res;
  distr_gaussian_exps = b.distr_gaussian_exps;
  distr_gaussian_mults = b.distr_gaussian_mults;
  distr_binomials = b.distr_binomials;

  resultsyesno = b.resultsyesno;
  run_yes = b.run_yes;
  posteriormode = b.posteriormode;

  FC_linears = b.FC_linears;
  design_psplines = b.design_psplines;
  FC_nonps = b.FC_nonps;
  FC_nonp_variances = b.FC_nonp_variances;
  FC_predicts = b.FC_predicts;

  design_hrandoms = b.design_hrandoms;
  FC_hrandom_variances = b.FC_hrandom_variances;

  }


const superbayesreg & superbayesreg::operator=(const superbayesreg & b)
  {
  if (this == & b)
	 return *this;
  statobject::operator=(statobject(b));
  create();
  #if defined(JAVA_OUTPUT_WINDOW)
  adminp_p = b.adminp_p;
  #endif

  pathres = b.pathres;
  title = b.title;
  pathnonp = b.pathnonp;

  statobj = b.statobj;

  D = b.D;

  modelvarnamesv = b.modelvarnamesv;

  equations=b.equations;
  simobj = b.simobj;

  master = b.master;

  generaloptions = b.generaloptions;

  distr_gaussians = b.distr_gaussians;
  distr_loggaussians = b.distr_loggaussians;
  distr_gaussian_res = b.distr_gaussian_res;
  distr_gaussian_exps = b.distr_gaussian_exps;
  distr_gaussian_mults = b.distr_gaussian_mults;
  distr_binomials = b.distr_binomials;

  resultsyesno = b.resultsyesno;
  run_yes = b.run_yes;
  posteriormode = b.posteriormode;

  FC_linears = b.FC_linears;
  design_psplines = b.design_psplines;
  FC_nonps = b.FC_nonps;
  FC_nonp_variances = b.FC_nonp_variances;
  FC_predicts = b.FC_predicts;

  design_hrandoms = b.design_hrandoms;
  FC_hrandom_variances = b.FC_hrandom_variances;

  return *this;
  }


int superbayesreg::parse(const ST::string & c)
  {

  int u = statobject::parse(c);

  int pos = statobject::parsecom(c,methods,globaloptions);

  if (pos >= 0)
	 (*functions[pos])(*this);

  }


bool superbayesreg::create_generaloptions(void)
  {


  if (iterations.getvalue()- burnin.getvalue() < 100)
    {
    outerror("ERROR: number of iterations must exceed number of burnin iterations about 100\n");
    return true;
    }

  if (step.getvalue() >= iterations.getvalue() - burnin.getvalue())
    {
    outerror("ERROR: thinning parameter too large\n");
    return true;
    }


  generaloptions = MCMC::GENERAL_OPTIONS(
  #if defined(JAVA_OUTPUT_WINDOW)
  adminb_p,
  #endif
  iterations.getvalue(),burnin.getvalue(),
                               step.getvalue(),logout,
                               level1.getvalue(),level2.getvalue());

  describetext.push_back("ESTIMATION OPTIONS:\n");
  describetext.push_back("\n");
  describetext.push_back("Number of Iterations: "
                            + ST::inttostring(iterations.getvalue()) + "\n");
  describetext.push_back("Burnin: " + ST::inttostring(burnin.getvalue()) + "\n");
  describetext.push_back("Thinning parameter: " +
                            ST::inttostring(step.getvalue()) + "\n");

//  generaloptions.nrout(iterationsprint.getvalue());

  return false;

  }


void superbayesreg::make_header(unsigned & modnr)
  {
  if (equations[modnr].hlevel == 1)
    {
    if (equations[modnr].equationtype == "mean")
      {
      equations[modnr].header = "MCMCREG OBJECT " + name.to_bstr() +
                              ": MAIN REGRESSION";
      equations[modnr].paths = "MAIN_REGRESSION";
      }
    else if (equations[modnr].equationtype == "variance")
      {
      equations[modnr].header = "MCMCREG OBJECT " + name.to_bstr() +
                              ": MAIN VARIANCE REGRESSION";
      equations[modnr].paths = "MAIN_VARIANCE_REGRESSION";
      }
    }
  else if (equations[modnr].hlevel == 2)
    {
    if (equations[modnr].equationtype == "mean")
      {
      equations[modnr].header = "MCMCREG OBJECT " + name.to_bstr() +
                              ": RANDOM EFFECTS REGRESSION";
      equations[modnr].paths = "RANDOM_EFFECTS";
      }
    else if (equations[modnr].equationtype == "variance")
      {
      equations[modnr].header = "MCMCREG OBJECT " + name.to_bstr() +
                              ": RANDOM EFFECTS VARIANCE REGRESSION";
      equations[modnr].paths = "RANDOM_EFFECTS_VARIANCE";
      }
    }

  }


void hregressrun(superbayesreg & b)
  {

  if ((b.equations.size()==0) || (b.run_yes==true))
    {
    b.clear();
    b.run_yes=false;
    b.generaloptions_yes = false;
    }

  b.resultsyesno = false;
  if (b.modeonly.getvalue() == true)
    b.posteriormode = true;
  else
    b.posteriormode = false;

  b.terms = b.modreg.getterms();

  b.describetext.erase(b.describetext.begin(),b.describetext.end());
  b.describetext.push_back("LAST ESTIMATED MODEL: \n");
  b.describetext.push_back("\n");
  b.describetext.push_back(b.modreg.getModelText());
  b.describetext.push_back("\n");

  b.equations.push_back(equation(b.equationnr.getvalue(),b.hlevel.getvalue(),
  b.equationtype.getvalue()));
  unsigned modnr = b.equations.size()-1;

  b.make_header(modnr);

//  b.outfiles.push_back(b.outfile.getvalue()+b.add_name);

  bool failure = false;

  if (!failure && b.generaloptions_yes==false)
    {
    failure = b.create_generaloptions();
    b.generaloptions_yes = true;
    }

  if (!failure)
    failure = b.create_distribution();
  if (!failure)
    failure = b.create_linear();
  if (!failure && b.terms.size() >= 1)
    failure = b.create_nonp();

  if (!failure)
    b.create_predict();


  if ((! failure) && (b.hlevel.getvalue() == 1) &&
      (b.equationtype.getvalue()=="mean"))
    {


    b.run_yes=true;

    b.simobj = MCMCsim(&b.generaloptions,b.equations);

    ST::string pathgraphs = b.outfile.getvalue();

    if (b.modeonly.getvalue())
      {
      failure = b.simobj.posteriormode(pathgraphs,false);
      }
    else
      failure = b.simobj.simulate(pathgraphs,b.setseed.getvalue(),true);

    if (!failure)
      b.resultsyesno = true;
    }

  }


bool superbayesreg::create_distribution(void)
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

  vector<ST::string> modelvarnamesh;
  modelvarnamesh = modreg.getModelVarnamesAsVector();
  modelvarnamesv.erase(modelvarnamesv.begin(),modelvarnamesv.end());
  for(i=0;i<modelvarnamesh.size();i++)
    if (modelvarnamesh[i] != "noconst")
      modelvarnamesv.push_back(modelvarnamesh[i]);

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
      if (notex[i] != "const" && notex[i] != "noconst")
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

//  datamatrix offs(1,1,0);
//  failure = create_offset(offs);


  datamatrix w;

  if (wn.length() > 0)
    {
    w = D.getCol(weightpos);
    }
  else
    w = datamatrix(1,1);



  describetext.push_back("Response distribution: "
                           + family.getvalue() + "\n");


  unsigned modnr = equations.size()-1;

//---------------------------- Gaussian response -------------------------------
  if (family.getvalue() == "gaussian")
    {

    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";

    distr_gaussians.push_back(DISTR_gaussian(aresp.getvalue(),bresp.getvalue(),
                                      &generaloptions,D.getCol(0),path,w) );

    equations[modnr].distrp = &distr_gaussians[distr_gaussians.size()-1];
    equations[modnr].pathd = defaultpath + "\\temp\\" + name  + "_scale.res";

    }
//-------------------------- END: Gaussian response ----------------------------
//-------------------------- log-Gaussian response -----------------------------
  else if (family.getvalue() == "loggaussian")
    {

    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";

    distr_loggaussians.push_back(DISTR_loggaussian(aresp.getvalue(),
                                       bresp.getvalue(),
                                      &generaloptions,D.getCol(0),path,w) );

    equations[modnr].distrp = &distr_loggaussians[distr_loggaussians.size()-1];
    equations[modnr].pathd = defaultpath + "\\temp\\" + name  + "_scale.res";

    }
//------------------------- END: log-Gaussian response -------------------------
//--------------- Gaussian response, exponential response function -------------
  else if (family.getvalue() == "gaussian_exp")
    {

    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";

    distr_gaussian_exps.push_back(DISTR_gaussian_exp(
                                  aresp.getvalue(),bresp.getvalue(),
                                  &generaloptions,D.getCol(0),path,w) );

    equations[modnr].distrp = &distr_gaussian_exps[distr_gaussian_exps.size()-1];
    equations[modnr].pathd = defaultpath + "\\temp\\" + name  + "_scale.res";

    }
//------------- END: Gaussian response, exponential response function ----------
//------------- Gaussian response, multiplicative random effects allowed -------
  else if (family.getvalue() == "gaussian_mult")
    {

    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";

    distr_gaussian_mults.push_back(DISTR_gaussian_mult(
                                  aresp.getvalue(),bresp.getvalue(),
                                  &generaloptions,D.getCol(0),path,w) );

    equations[modnr].distrp = &distr_gaussian_mults[distr_gaussian_mults.size()-1];
    equations[modnr].pathd = defaultpath + "\\temp\\" + name  + "_scale.res";

    }
//---------- END: Gaussian response, multiplicative random effects allowed -----
//-------------------------- Gaussian random effect ----------------------------
  else if (family.getvalue() == "gaussian_re")
    {

    distr_gaussian_res.push_back(DISTR_gaussian_re(
                                &generaloptions,D.getCol(0),w) );

    equations[modnr].distrp = &distr_gaussian_res[distr_gaussian_res.size()-1];
    equations[modnr].pathd = "";
    }
//------------------------- END: Gaussian random effect ------------------------
//---------------------------- Binomial response -------------------------------
  else if (family.getvalue() == "binomial_logit")
    {

    distr_binomials.push_back(DISTR_binomial(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_binomials[distr_binomials.size()-1];
    equations[modnr].pathd = "";

    }
//-------------------------- END: Binomial response ----------------------------


  equations[modnr].distrp->responsename=rname;
  equations[modnr].distrp->weightname=wn;

  if (equations[modnr].hlevel==1)
    master.level1_likep = equations[modnr].distrp;

  return false;

  }


void superbayesreg::extract_data(unsigned i, datamatrix & d,datamatrix & iv)
  {
  int j1,j2;
  unsigned p = terms[i].varnames.size()-1;
  j1 = terms[i].varnames[p].isinlist(modelvarnamesv);
  d = D.getCol(j1);

  if (terms[i].varnames.size() > 1)
    {
    j2 = terms[i].varnames[0].isinlist(modelvarnamesv);
    iv = D.getCol(j2);
    }
  }


void superbayesreg::create_predict(void)
  {
  if (predict.getvalue() != "no")
    {

    unsigned modnr = equations.size()-1;

    ST::string h = equations[modnr].paths;

    ST::string pathnonp = defaultpath + "\\temp\\" + name + "_" + h +
                            "_predict.raw";

    ST::string pathnonp2 = defaultpath + "\\temp\\" + name + "_" + h +
                            "_deviance.raw";


    ST::string pathres = outfile.getvalue() +  "_" + h + "_predict.res";

    FC_predicts.push_back(FC_predict(&generaloptions,
                         equations[modnr].distrp,"",pathnonp,pathnonp2,D,modelvarnamesv));

    equations[modnr].add_FC(&FC_predicts[FC_predicts.size()-1],pathres);


    }


  }


bool superbayesreg::create_linear(void)
  {

  unsigned modnr = equations.size()-1;

  unsigned i;
  int j;

  vector<ST::string> varnames;
  vector<ST::string> varnamesh =  lineareffects.get_constvariables(terms);

  bool noconst=false;

  for(i=0;i<varnamesh.size();i++)
    if (varnamesh[i] == "noconst")
      noconst=true;

  if (equations[modnr].hlevel == 1 && (noconst==false))
    varnames.push_back("const");

  for(i=0;i<varnamesh.size();i++)
    if (varnamesh[i] != "noconst")
      varnames.push_back(varnamesh[i]);

  unsigned nr = varnames.size();


  ST::string title;
  ST::string pathconst;
  ST::string pathconstres;

  ST::string h = equations[modnr].paths;

  title = h + ": linear effects" ;

  pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                           + "_LinearEffects"  +
                           "_" + h + ".raw";

  pathconstres = outfile.getvalue() + "_" + h + "_LinearEffects.res";

  if (pathconst.isvalidfile() == 1)
    {
    errormessages.push_back("ERROR: unable to open file " + pathconst +
                                 " for writing\n");
    return true;
    }

  datamatrix X;
  if (nr > 0)
    X = datamatrix(D.rows(),nr,1);

  for (i=0;i<varnames.size();i++)
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

  bool store=false ;
  if (storesample.getvalue() == true)
    store = true;

  FC_linears.push_back(FC_linear(&master,&generaloptions,equations[modnr].distrp,X,
                         varnames,title,pathconst,store));

  equations[modnr].add_FC(&FC_linears[FC_linears.size()-1],pathconstres);

  return false;

  }


void superbayesreg::create_pspline(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_pspline.raw","nonlinear_pspline_effect_of"," Nonlinear effect (P-spline) of ");

  datamatrix d,iv;
  extract_data(i,d,iv);

  // TEST
  //  ofstream out("c:\\bayesx\\testh\\results\\d.res");
  //  d.prettyPrint(out);
  // TEST

  design_psplines.push_back(DESIGN_pspline(d,iv,equations[modnr].distrp,
                            &FC_linears[FC_linears.size()-1],
                            terms[i].options,terms[i].varnames));

  bool store=false ;
  if (storesample.getvalue() == true)
    store = true;

  FC_nonps.push_back(FC_nonp(&master,&generaloptions,equations[modnr].distrp,title,
                     pathnonp,&design_psplines[design_psplines.size()-1],
                     terms[i].options,terms[i].varnames,store));


  equations[modnr].add_FC(&FC_nonps[FC_nonps.size()-1],pathres);


  // variances

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_pspline_var.raw","variance_of_nonlinear_pspline_effect_of","Variance of nonlinear effect of ");

  FC_nonp_variances.push_back(FC_nonp_variance(&generaloptions,equations[modnr].distrp,
                                title,pathnonp,&design_psplines[design_psplines.size()-1],
                                &FC_nonps[FC_nonps.size()-1],terms[i].options,
                                terms[i].varnames,store));

  equations[modnr].add_FC(&FC_nonp_variances[FC_nonp_variances.size()-1],pathres);

  }


bool superbayesreg::findREdistr(ST::string & na,equation & maine,unsigned & fnr)
  {
  bool found = false;
  unsigned i;
  for (i=0;i<equations.size()-1;i++)
    {

/*
    TEST
    int hl = equations[i].hlevel;
    ST::string typem = maine.equationtype;
    ST::string type = equations[i].equationtype;
    ST::string name = equations[i].distrp->responsename;
*/

    if ((equations[i].distrp->responsename==na) &&
        (equations[i].hlevel==2) &&
        (equations[i].equationtype==maine.equationtype)
        )
      {
      found = true;
      fnr = i;
      }
    }
  return found;
  }


bool superbayesreg::create_hrandom(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_hrandom.raw","random_effect_of","Random effect of ");

  ST::string pathnonp2 = pathnonp.substr(0,pathnonp.length()-4)  + "_2.res";


  datamatrix d,iv;
  extract_data(i,d,iv);

  unsigned fnr;
  ST::string na = terms[i].varnames[terms[i].varnames.size()-1];
  bool found = findREdistr(na,equations[modnr],fnr);

  if (found==false)
    {
    outerror("ERROR: level 2 equation for variable " + na + " not found \n");
    return true;
    }


  design_hrandoms.push_back(DESIGN_hrandom(d,iv,equations[modnr].distrp,
                            &FC_linears[FC_linears.size()-1],
                             equations[fnr].distrp,
                            terms[i].options,terms[i].varnames));

  bool store=false ;
  if (storesample.getvalue() == true)
    store = true;

  FC_hrandoms.push_back(FC_hrandom(&master,&generaloptions,equations[modnr].distrp,
                        equations[fnr].distrp, title,pathnonp,pathnonp2,
                        &design_hrandoms[design_hrandoms.size()-1],
                        terms[i].options,terms[i].varnames,store));

  equations[modnr].add_FC(&FC_hrandoms[FC_hrandoms.size()-1],pathres);

  // variances

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_hrandom_var.raw","variance_of_random_effect_of","Variance of random effect of ");

  FC_hrandom_variances.push_back(FC_hrandom_variance(&generaloptions,equations[modnr].distrp,
                                 equations[fnr].distrp,
                                title,pathnonp,&design_hrandoms[design_hrandoms.size()-1],
                                &FC_hrandoms[FC_hrandoms.size()-1],terms[i].options,
                                terms[i].varnames,store));

  equations[modnr].add_FC(&FC_hrandom_variances[FC_hrandom_variances.size()-1],pathres);

  return false;

  }


bool  superbayesreg::create_random_pspline(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  terms[i].options[12] = "true";
  bool multexp = false;
  if (terms[i].options[0] == "hrandomexp_pspline")
    {
    terms[i].options[17] = "true";
    multexp=true;
    }

  create_pspline(i);
  FC_nonp * fcnp_pspline = &FC_nonps[FC_nonps.size()-1];
  MCMC::DESIGN * dp_pspline = &design_psplines[design_psplines.size()-1];
  dp_pspline->changingdesign=true;
  datamatrix effect(D.rows(),1,1);
  dp_pspline->set_intvar(effect,0);

  FC_mults.push_back(FC_mult(true,multexp));
  equations[modnr].add_FC(&FC_mults[FC_mults.size()-1],"");

  term helpt = terms[i];
  terms[i].varnames.erase(terms[i].varnames.begin(),terms[i].varnames.end());
  terms[i].varnames.push_back(helpt.varnames[1]);
  terms[i].varnames.push_back(helpt.varnames[0]);

  terms[i].options[12] = "true";
  create_hrandom(i);
  FC_nonp * fcnp_hrandom = &FC_hrandoms[FC_hrandoms.size()-1];
  MCMC::DESIGN * dp_hrandom = &design_hrandoms[design_hrandoms.size()-1];
  dp_hrandom->changingdesign=true;

  FC_mults.push_back(FC_mult(false,multexp));

  FC_mults[FC_mults.size()-2].set_effectp(dp_pspline,fcnp_pspline);
  FC_mults[FC_mults.size()-2].set_intp(dp_hrandom,fcnp_hrandom);

  FC_mults[FC_mults.size()-1].set_effectp(dp_hrandom,fcnp_hrandom);
  FC_mults[FC_mults.size()-1].set_intp(dp_pspline,fcnp_pspline);

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_mult.raw","_multiplicative_effect_of",
             "Multiplicative effect of ");

  bool samplem;
  if (terms[i].options[13] == "false")
    samplem = false;
  else
    samplem = true;

  bool store=false ;
  if (storesample.getvalue() == true)
    store = true;

  FC_mults[FC_mults.size()-1].set_multeffects(&generaloptions,title,pathnonp,
         samplem,store);

  equations[modnr].add_FC(&FC_mults[FC_mults.size()-1],pathres);

  return false;

  }


bool superbayesreg::create_mrf(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_spatial.raw","spatial_MRF_effect_of","Spatial effect (MRF) of ");

  datamatrix d,iv;
  extract_data(i,d,iv);

  mapobject * mapp;                           // pointer to mapobject

  int objpos = findstatobject(*statobj,terms[i].options[8],"map");

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
  bool isconnected = m.isconnected();
  if (isconnected==false)
    {
    outerror("ERROR: map is disconnected, spatial effect cannot be estimated\n");
    return true;
    }

  design_mrfs.push_back(DESIGN_mrf(d,iv,equations[modnr].distrp,
                           &FC_linears[FC_linears.size()-1],m,
                            terms[i].options,terms[i].varnames));

  bool store=false ;
  if (storesample.getvalue() == true)
    store = true;

  FC_nonps.push_back(FC_nonp(&master,&generaloptions,equations[modnr].distrp,title,
                     pathnonp,&design_mrfs[design_mrfs.size()-1],
                     terms[i].options,terms[i].varnames,store));

  equations[modnr].add_FC(&FC_nonps[FC_nonps.size()-1],pathres);

  // variances

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_spatial_var.raw","variance_of_spatial_MRF_effect_of","Variance of spatial effect of ");

  FC_nonp_variances.push_back(FC_nonp_variance(&generaloptions,equations[modnr].distrp,
                                title,pathnonp,&design_mrfs[design_mrfs.size()-1],
                                &FC_nonps[FC_nonps.size()-1],terms[i].options,
                                terms[i].varnames,store));

  equations[modnr].add_FC(&FC_nonp_variances[FC_nonp_variances.size()-1],pathres);

  return false;
  }


bool superbayesreg::create_nonp(void)
  {

  unsigned i;
  bool error=false;

  for(i=0;i<terms.size();i++)
    {
    if (terms[i].options.size() > 0)
      {
      if (terms[i].options[0] == "pspline")
        create_pspline(i);
      if (terms[i].options[0] == "hrandom")
        error = create_hrandom(i);
      if (terms[i].options[0] == "spatial")
        error = create_mrf(i);
      if ((terms[i].options[0] == "hrandom_pspline") ||
          (terms[i].options[0] == "hrandomexp_pspline"))
        error = create_random_pspline(i);
      }

    if (error)
      return error;
    }

  return false;
  }

void autocorrrun(superbayesreg & b)
  {

  if (b.resultsyesno==true)
	 {
     if (b.posteriormode == false)
       {
       ST::string path = b.outfile.getvalue()+  "_autocor" + ".raw";
       if (b.generaloptions.samplesize < unsigned(b.maxlag.getvalue()*4))
         b.outerror("ERROR: samplesize too small\n");
       else
         b.simobj.autocorr(b.maxlag.getvalue());
       }
     else
       b.outerror("ERROR: no MCMC simulation results\n");
	 }
  else
	 b.outerror("ERROR: no regression results\n");

  }


void getsamplerun(superbayesreg & b)
  {

  if (b.resultsyesno == true)
    {
    if (b.posteriormode == false)
      {
      #if defined(JAVA_OUTPUT_WINDOW)

      b.simobj.get_samples(b.newcommands,b.outfile.getvalue() + "_");
      #else
      b.simobj.get_samples();
      #endif
      }
    else
      b.outerror("ERROR: no MCMC simulation results\n");
    }
  else
    b.outerror("ERROR: no regression results\n");

  }


void superbayesreg::describe(const optionlist & globaloptions)
  {
  statobject::describe(globaloptions);
  }


#if defined(BORLAND_OUTPUT_WINDOW)
//------------------------------------------------------------------------------
#pragma package(smart_init)
#endif











