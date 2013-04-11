/* BayesX - Software for Bayesian Inference in
Structured Additive Regression Models.
Copyright (C) 2011  Christiane Belitz, Andreas Brezger,
Thomas Kneib, Stefan Lang, Nikolaus Umlauf

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */






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
#if defined(__BUILDING_LINUX)
    pathnonp = defaultpath + "/temp/" + name + "_" + h + "_f_"
                                  + varname1 + endingraw;
#else
    pathnonp = defaultpath + "\\temp\\" + name + "_" + h + "_f_"
                                  + varname1 + endingraw;
#endif

    pathres = outfile.getvalue() + "_" + h + "_" +
                     endingres + "_" + varname1 + ".res";

    title = h + ": " + endingtitle + varname1;

    }
  else
    {
#if defined(__BUILDING_LINUX)
    pathnonp = defaultpath + "/temp/" + name + "_" + h + "_" + varname2
                 + "_f_"  + varname1 + endingraw;
#else
    pathnonp = defaultpath + "\\temp\\" + name + "_" + h + "_" + varname2
                 + "_f_"  + varname1 + endingraw;
#endif

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
  tnames.push_back("kriging");
  tnames.push_back("geokriging");
  tnames.push_back("hrandom_pspline");
  tnames.push_back("hrandomexp_pspline");
  tnames.push_back("ridge");
  tnames.push_back("lasso");
  tnames.push_back("offset");


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

  standardize = simpleoption("standardize",false);

  families.reserve(80);
  families.push_back("gaussian");
  families.push_back("hetgaussian");
  families.push_back("vargaussian");
  families.push_back("gaussian_mixture");
  families.push_back("loggaussian");
  families.push_back("quantreg");
  families.push_back("gaussian_re");
  families.push_back("gaussian_exp");
  families.push_back("gaussian_mult");
  families.push_back("binomial_logit");
  families.push_back("poisson");
  families.push_back("poisson_ext");
  families.push_back("poisson_extlin");
  families.push_back("binomial_probit");
  families.push_back("binomial_svm");
  families.push_back("multinom_probit");
  families.push_back("multgaussian");
  families.push_back("binomial_logit_l1");
  families.push_back("multinom_logit");
  families.push_back("zip_lambda");
  families.push_back("zip_pi");
  families.push_back("zinb_mu");
  families.push_back("zinb_pi");
  families.push_back("zinb_delta");
  families.push_back("zip_pi_cloglog");
  families.push_back("zip_lambda_cloglog");
  families.push_back("negbin_mu");
  families.push_back("negbin_delta");
  families.push_back("beta_mu");
  families.push_back("beta_sigma2");
  family = stroption("family",families,"gaussian");
  aresp = doubleoption("aresp",0.001,-1.0,500);
  bresp = doubleoption("bresp",0.001,0.0,500);
  H = intoption("H",6,2,6);

  hlevel = intoption("hlevel",1,1,3);
  equationnr = intoption("equation",1,1,50);
  equationtypes.reserve(20);
  equationtypes.push_back("mean");
  equationtypes.push_back("meanservant");
  equationtypes.push_back("variance");
  equationtypes.push_back("pi");
  equationtypes.push_back("delta");
  equationtypes.push_back("sigma2");
  equationtype = stroption("equationtype",equationtypes,"mean");

  predictop.reserve(20);
  predictop.push_back("no");
  predictop.push_back("nosamples");
  predictop.push_back("full");
  predictop.push_back("fulls");
  predictop.push_back("predictor");

  predict = stroption("predict",predictop,"no");

  pred_check = simpleoption("pred_check",false);

  cv = simpleoption("cv",false);

  MSEop.reserve(20);
  MSEop.push_back("no");
  MSEop.push_back("yes");
  MSEop.push_back("quadratic");
  MSEop.push_back("check");
  mse = stroption("MSE",MSEop,"no");

  mseparam = doubleoption("mseparam",0.5,-999999999,999999999);

  centerlinear = simpleoption("centerlinear",false);

  quantile = doubleoption("quantile",0.5,0.001,0.999);

  includescale = simpleoption("includescale",false);

  modemaxit = intoption("modemaxit",1000,0,5000);

  stopsum = doubleoption("stopsum",1,0.9,1);
  stoprmax = intoption("stoprmax",5,1,100000);
  fraclimit = doubleoption("fraclimit",0.001,0,1);

  scaleconst = simpleoption("scaleconst",false);

  utilities = simpleoption("utilities",false);

  aexp = doubleoption("aexp",0,0.00000000001,100000000);
  bexp = doubleoption("bexp",1,0.00000000001,100000000);

  adaptexp = simpleoption("adaptexp",false);

  slow = simpleoption("slow",false);

  regressoptions.reserve(100);

  regressoptions.push_back(&modeonly);
  regressoptions.push_back(&setseed);
  regressoptions.push_back(&iterations);
  regressoptions.push_back(&burnin);
  regressoptions.push_back(&step);
  regressoptions.push_back(&level1);
  regressoptions.push_back(&level2);
  regressoptions.push_back(&family);
  regressoptions.push_back(&aresp);
  regressoptions.push_back(&bresp);
  regressoptions.push_back(&H);
  regressoptions.push_back(&hlevel);
  regressoptions.push_back(&equationnr);
  regressoptions.push_back(&equationtype);
  regressoptions.push_back(&predict);
  regressoptions.push_back(&pred_check);
  regressoptions.push_back(&mse);
  regressoptions.push_back(&mseparam);
  regressoptions.push_back(&centerlinear);
  regressoptions.push_back(&quantile);
  regressoptions.push_back(&cv);
  regressoptions.push_back(&includescale);
  regressoptions.push_back(&standardize);
  regressoptions.push_back(&modemaxit);
  regressoptions.push_back(&stopsum);
  regressoptions.push_back(&stoprmax);
  regressoptions.push_back(&fraclimit);
  regressoptions.push_back(&scaleconst);
  regressoptions.push_back(&utilities);
  regressoptions.push_back(&aexp);
  regressoptions.push_back(&bexp);
  regressoptions.push_back(&adaptexp);
  regressoptions.push_back(&slow);  

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

#if defined(__BUILDING_LINUX)
  ST::string h = defaultpath+"/output/"+name;
#else
  ST::string h = defaultpath+"\\output\\"+name;
#endif

  outfile = fileoption("outfile",h,false);

  globaloptions.push_back(&outfile);

  create_hregress();

  create_autocorr();

  create_getsample();

  }


void superbayesreg::clear(void)
  {

  equations.erase(equations.begin(),equations.end());
  equations.reserve(20);

  distr_gaussians.erase(distr_gaussians.begin(),distr_gaussians.end());
  distr_gaussians.reserve(20);

  distr_hetgaussians.erase(distr_hetgaussians.begin(),distr_hetgaussians.end());
  distr_hetgaussians.reserve(20);

  distr_vargaussians.erase(distr_vargaussians.begin(),distr_vargaussians.end());
  distr_vargaussians.reserve(20);

  distr_gaussianmixtures.erase(distr_gaussianmixtures.begin(),
  distr_gaussianmixtures.end());
  distr_gaussianmixtures.reserve(20);

  distr_quantregs.erase(distr_quantregs.begin(),distr_quantregs.end());
  distr_quantregs.reserve(20);

  distr_loggaussians.erase(distr_loggaussians.begin(),distr_loggaussians.end());
  distr_loggaussians.reserve(20);

  distr_gaussian_res.erase(distr_gaussian_res.begin(),distr_gaussian_res.end());
  distr_gaussian_res.reserve(20);

  distr_gaussian_exps.erase(distr_gaussian_exps.begin(),distr_gaussian_exps.end());
  distr_gaussian_exps.reserve(20);

  distr_gaussian_mults.erase(distr_gaussian_mults.begin(),
  distr_gaussian_mults.end());
  distr_gaussian_mults.reserve(20);

  distr_binomials.erase(distr_binomials.begin(),distr_binomials.end());
  distr_binomials.reserve(20);

  distr_poissons.erase(distr_poissons.begin(),distr_poissons.end());
  distr_poissons.reserve(20);

  distr_poisson_exts.erase(distr_poisson_exts.begin(),distr_poisson_exts.end());
  distr_poisson_exts.reserve(20);

  distr_poisson_extlins.erase(distr_poisson_extlins.begin(),
                              distr_poisson_extlins.end());
  distr_poisson_extlins.reserve(20);

  distr_binomialprobits.erase(distr_binomialprobits.begin(),
                              distr_binomialprobits.end());
  distr_binomialprobits.reserve(20);

  distr_binomialsvms.erase(distr_binomialsvms.begin(),
                              distr_binomialsvms.end());
  distr_binomialsvms.reserve(20);

  distr_multinomprobits.erase(distr_multinomprobits.begin(),
                              distr_multinomprobits.end());
  distr_multinomprobits.reserve(20);

  distr_multgaussians.erase(distr_multgaussians.begin(),
                              distr_multgaussians.end());
  distr_multgaussians.reserve(20);

  distr_multinomlogits.erase(distr_multinomlogits.begin(),
                              distr_multinomlogits.end());
  distr_multinomlogits.reserve(20);

  distr_logit_fruehwirths.erase(distr_logit_fruehwirths.begin(),
                               distr_logit_fruehwirths.end());
  distr_logit_fruehwirths.reserve(20);

  distr_ziplambdas.erase(distr_ziplambdas.begin(),distr_ziplambdas.end());
  distr_ziplambdas.reserve(20);

  distr_zippis.erase(distr_zippis.begin(),distr_zippis.end());
  distr_zippis.reserve(20);

  distr_negbinzip_mus.erase(distr_negbinzip_mus.begin(),distr_negbinzip_mus.end());
  distr_negbinzip_mus.reserve(20);

  distr_negbinzip_pis.erase(distr_negbinzip_pis.begin(),distr_negbinzip_pis.end());
  distr_negbinzip_pis.reserve(20);

  distr_negbinzip_deltas.erase(distr_negbinzip_deltas.begin(),distr_negbinzip_deltas.end());
  distr_negbinzip_deltas.reserve(20);

  distr_zip_cloglog_mus.erase(distr_zip_cloglog_mus.begin(),distr_zip_cloglog_mus.end());
  distr_zip_cloglog_mus.reserve(20);

  distr_zip_cloglog_pis.erase(distr_zip_cloglog_pis.begin(),distr_zip_cloglog_pis.end());
  distr_zip_cloglog_pis.reserve(20);

  distr_negbin_mus.erase(distr_negbin_mus.begin(),distr_negbin_mus.end());
  distr_negbin_mus.reserve(20);

  distr_negbin_deltas.erase(distr_negbin_deltas.begin(),distr_negbin_deltas.end());
  distr_negbin_deltas.reserve(20);

  distr_beta_mus.erase(distr_beta_mus.begin(),distr_beta_mus.end());
  distr_beta_mus.reserve(20);

  distr_beta_sigma2s.erase(distr_beta_sigma2s.begin(),distr_beta_sigma2s.end());
  distr_beta_sigma2s.reserve(20);

  FC_linears.erase(FC_linears.begin(),FC_linears.end());
  FC_linears.reserve(50);

  design_psplines.erase(design_psplines.begin(),design_psplines.end());
  design_psplines.reserve(100);

  design_mrfs.erase(design_mrfs.begin(),design_mrfs.end());
  design_mrfs.reserve(30);

  design_krigings.erase(design_krigings.begin(),design_krigings.end());
  design_krigings.reserve(30);

  design_hrandoms.erase(design_hrandoms.begin(),design_hrandoms.end());
  design_hrandoms.reserve(50);

  FC_nonps.erase(FC_nonps.begin(),FC_nonps.end());
  FC_nonps.reserve(200);

  FC_hrandoms.erase(FC_hrandoms.begin(),FC_hrandoms.end());
  FC_hrandoms.reserve(50);

  FC_mults.erase(FC_mults.begin(),FC_mults.end());
  FC_mults.reserve(30);

  FC_nonp_variances.erase(FC_nonp_variances.begin(),FC_nonp_variances.end());
  FC_nonp_variances.reserve(200);


  FC_linear_pens.erase(FC_linear_pens.begin(),FC_linear_pens.end());
  FC_linear_pens.reserve(50);

  FC_variance_pen_vectors.erase(FC_variance_pen_vectors.begin(),
                                FC_variance_pen_vectors.end());
  FC_variance_pen_vectors.reserve(50);


  FC_hrandom_variances.erase(FC_hrandom_variances.begin(),
  FC_hrandom_variances.end());
  FC_hrandom_variances.reserve(50);

  FC_hrandom_variance_vecs.erase(FC_hrandom_variance_vecs.begin(),
  FC_hrandom_variance_vecs.end());
  FC_hrandom_variance_vecs.reserve(50);

  FC_hrandom_variance_vec_nmigs.erase(FC_hrandom_variance_vec_nmigs.begin(),
  FC_hrandom_variance_vec_nmigs.end());
  FC_hrandom_variance_vec_nmigs.reserve(50);

  FC_hrandom_variance_ssvss.erase(FC_hrandom_variance_ssvss.begin(),
  FC_hrandom_variance_ssvss.end());
  FC_hrandom_variance_ssvss.reserve(50);

  FC_predicts.erase(FC_predicts.begin(),FC_predicts.end());
  FC_predicts.reserve(30);

  FC_predicts_mult.erase(FC_predicts_mult.begin(),FC_predicts_mult.end());
  FC_predicts_mult.reserve(30);

  predict_mult_distrs.erase(predict_mult_distrs.begin(),predict_mult_distrs.end());
  predict_mult_distrs.reserve(30);

  FC_predict_predictors.erase(FC_predict_predictors.begin(),
                              FC_predict_predictors.end());
  FC_predict_predictors.reserve(30);


  FC_predictive_checks.erase(FC_predictive_checks.begin(),
  FC_predictive_checks.end());
  FC_predictive_checks.reserve(30);

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
  computemodeforstartingvalues = true;
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
  computemodeforstartingvalues = true;
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
  nrlevel1 = b.nrlevel1;
  simobj = b.simobj;

  master = b.master;

  generaloptions = b.generaloptions;

  distr_gaussians = b.distr_gaussians;
  distr_hetgaussians = b.distr_hetgaussians;
  distr_vargaussians = b.distr_vargaussians;
  distr_quantregs = b.distr_quantregs;
  distr_gaussianmixtures = b.distr_gaussianmixtures;
  distr_loggaussians = b.distr_loggaussians;
  distr_gaussian_res = b.distr_gaussian_res;
  distr_gaussian_exps = b.distr_gaussian_exps;
  distr_gaussian_mults = b.distr_gaussian_mults;
  distr_binomials = b.distr_binomials;
  distr_poissons = b.distr_poissons;
  distr_poisson_exts = b.distr_poisson_exts;
  distr_poisson_extlins = b.distr_poisson_extlins;
  distr_binomialprobits = b.distr_binomialprobits;
  distr_binomialsvms = b.distr_binomialsvms;
  distr_multinomprobits = b.distr_multinomprobits;
  distr_multgaussians = b.distr_multgaussians;
  distr_multinomlogits = b.distr_multinomlogits;
  distr_logit_fruehwirths = b.distr_logit_fruehwirths;
  distr_ziplambdas = b.distr_ziplambdas;
  distr_zippis = b.distr_zippis;
  distr_negbinzip_mus = b.distr_negbinzip_mus;
  distr_negbinzip_pis = b.distr_negbinzip_pis;
  distr_negbinzip_deltas = b.distr_negbinzip_deltas;
  distr_zip_cloglog_mus = b.distr_zip_cloglog_mus;
  distr_zip_cloglog_pis = b.distr_zip_cloglog_pis;
  distr_negbin_mus = b.distr_negbin_mus;
  distr_negbin_deltas = b.distr_negbin_deltas;
  distr_beta_mus = b.distr_beta_mus;
  distr_beta_sigma2s = b.distr_beta_sigma2s;

  resultsyesno = b.resultsyesno;
  run_yes = b.run_yes;
  posteriormode = b.posteriormode;
  computemodeforstartingvalues = b.computemodeforstartingvalues;

  FC_linears = b.FC_linears;
  design_psplines = b.design_psplines;
  FC_nonps = b.FC_nonps;
  FC_nonp_variances = b.FC_nonp_variances;
  FC_predicts = b.FC_predicts;
  FC_predicts_mult = b.FC_predicts_mult;
  FC_predict_predictors = b.FC_predict_predictors;
  FC_predictive_checks = b.FC_predictive_checks;

  ridge = b.ridge;
  lasso = b.lasso;
  ridge_linear = b.ridge_linear;
  lasso_linear = b.lasso_linear;
  FC_variance_pen_vectors = b.FC_variance_pen_vectors;
  FC_linear_pens = b.FC_linear_pens;

  design_hrandoms = b.design_hrandoms;
  FC_hrandom_variances = b.FC_hrandom_variances;
  FC_hrandom_variance_vecs = b.FC_hrandom_variance_vecs;
  FC_hrandom_variance_vec_nmigs = b.FC_hrandom_variance_vec_nmigs;
  FC_hrandom_variance_ssvss = b.FC_hrandom_variance_ssvss;

  design_krigings = b.design_krigings;
  design_mrfs = b.design_mrfs;

  errors = b.errors;

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
  nrlevel1 = b.nrlevel1;
  simobj = b.simobj;

  master = b.master;

  generaloptions = b.generaloptions;

  distr_gaussians = b.distr_gaussians;
  distr_hetgaussians = b.distr_hetgaussians;
  distr_vargaussians = b.distr_vargaussians;
  distr_quantregs = b.distr_quantregs;
  distr_gaussianmixtures = b.distr_gaussianmixtures;
  distr_loggaussians = b.distr_loggaussians;
  distr_gaussian_res = b.distr_gaussian_res;
  distr_gaussian_exps = b.distr_gaussian_exps;
  distr_gaussian_mults = b.distr_gaussian_mults;
  distr_binomials = b.distr_binomials;
  distr_poissons = b.distr_poissons;
  distr_poisson_exts = b.distr_poisson_exts;
  distr_poisson_extlins = b.distr_poisson_extlins;
  distr_binomialprobits = b.distr_binomialprobits;
  distr_binomialsvms = b.distr_binomialsvms;
  distr_multinomprobits = b.distr_multinomprobits;
  distr_multgaussians = b.distr_multgaussians;
  distr_multinomlogits = b.distr_multinomlogits;
  distr_logit_fruehwirths = b.distr_logit_fruehwirths;
  distr_ziplambdas = b.distr_ziplambdas;
  distr_zippis = b.distr_zippis;
  distr_negbinzip_mus = b.distr_negbinzip_mus;
  distr_negbinzip_pis = b.distr_negbinzip_pis;
  distr_negbinzip_deltas = b.distr_negbinzip_deltas;
  distr_zip_cloglog_mus = b.distr_zip_cloglog_mus;
  distr_zip_cloglog_pis = b.distr_zip_cloglog_pis;
  distr_negbin_mus = b.distr_negbin_mus;
  distr_negbin_deltas = b.distr_negbin_deltas;
  distr_beta_mus = b.distr_beta_mus;
  distr_beta_sigma2s = b.distr_beta_sigma2s;


  resultsyesno = b.resultsyesno;
  run_yes = b.run_yes;
  posteriormode = b.posteriormode;
  computemodeforstartingvalues = b.computemodeforstartingvalues;

  FC_linears = b.FC_linears;
  design_psplines = b.design_psplines;
  FC_nonps = b.FC_nonps;
  FC_nonp_variances = b.FC_nonp_variances;
  FC_predicts = b.FC_predicts;
  FC_predicts_mult = b.FC_predicts_mult;
  FC_predict_predictors = b.FC_predict_predictors;
  FC_predictive_checks = b.FC_predictive_checks;

  ridge = b.ridge;
  lasso = b.lasso;
  ridge_linear = b.ridge_linear;
  lasso_linear = b.lasso_linear;
  FC_variance_pen_vectors = b.FC_variance_pen_vectors;
  FC_linear_pens = b.FC_linear_pens;

  design_hrandoms = b.design_hrandoms;
  FC_hrandom_variances = b.FC_hrandom_variances;
  FC_hrandom_variance_vecs = b.FC_hrandom_variance_vecs;
  FC_hrandom_variance_vec_nmigs = b.FC_hrandom_variance_vec_nmigs;
  FC_hrandom_variance_ssvss = b.FC_hrandom_variance_ssvss;

  design_krigings = b.design_krigings;
  design_mrfs = b.design_mrfs;

  errors = b.errors;

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

    ST::string rn = equations[modnr].distrp->responsename;
    if (equations[modnr].equationtype == "mean")
      {
      equations[modnr].header = "MCMCREG OBJECT " + name.to_bstr() +
                              ": MAIN REGRESSION_"+ rn;
      equations[modnr].paths = "MAIN_REGRESSION_"+ rn;
      }
    else if (equations[modnr].equationtype == "variance")
      {
      equations[modnr].header = "MCMCREG OBJECT " + name.to_bstr() +
                              ": MAIN VARIANCE REGRESSION_"+ rn;
      equations[modnr].paths = "MAIN_VARIANCE_REGRESSION_"+ rn;
      }
    else if (equations[modnr].equationtype == "pi")
      {
      equations[modnr].header = "MCMCREG OBJECT " + name.to_bstr() +
                              ": MAIN ZERO INFLATION REGRESSION_"+ rn;
      equations[modnr].paths = "MAIN_ZERO_INFLATION_REGRESSION_"+ rn;
      }
    else if (equations[modnr].equationtype == "delta")
      {
      equations[modnr].header = "MCMCREG OBJECT " + name.to_bstr() +
                              ": MAIN DELTA REGRESSION_"+ rn;
      equations[modnr].paths = "MAIN_DELTA_REGRESSION_"+ rn;
      }
    else if (equations[modnr].equationtype == "sigma2")
      {
      equations[modnr].header = "MCMCREG OBJECT " + name.to_bstr() +
                              ": MAIN SIGMA2 REGRESSION_"+ rn;
      equations[modnr].paths = "MAIN_SIGMA2_REGRESSION_"+ rn;
      }
    else if (equations[modnr].equationtype == "meanservant")
      {
      equations[modnr].header = "MCMCREG OBJECT " + name.to_bstr() +
                              ": MAIN REGRESSION_" + rn;
      equations[modnr].paths = "MAIN_REGRESSION_" + rn;
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
    else if (equations[modnr].equationtype == "pi")
      {
      equations[modnr].header = "MCMCREG OBJECT " + name.to_bstr() +
                              ": RANDOM EFFECTS ZERO INFLATION REGRESSION";
      equations[modnr].paths = "RANDOM_EFFECTS_ZERO_INFLATION";
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
    b.errors=false;
    b.nrlevel1 = 0;
    }

  if (b.errors==false)    // errors from level 2 equations
    {
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

//    b.make_header(modnr);

    bool failure = false;

    if (!failure && b.generaloptions_yes==false)
      {
      failure = b.create_generaloptions();
      b.generaloptions_yes = true;
      }

    if (!failure)
      failure = b.create_distribution();

    if (!failure)
      b.make_header(modnr);


    if (!failure)
      failure = b.create_linear();
    if (!failure && b.terms.size() >= 1)
      failure = b.create_nonp();

    if (!failure)
      failure = b.create_predict();

    if (!failure)
      b.create_predictive_check();

    if (!failure)
      b.create_cv();

    if (b.hlevel.getvalue() == 1)
      b.nrlevel1++;

    if ((! failure) && (b.hlevel.getvalue() == 1) &&
        (b.equationtype.getvalue()=="mean"))
      {

      if (!failure)
        failure = b.check_errors();

      if (!failure)
        {
        b.run_yes=true;

        unsigned mit = b.modemaxit.getvalue();
        b.simobj = MCMCsim(&b.generaloptions,b.equations,mit);

        ST::string pathgraphs = b.outfile.getvalue();

        if (b.modeonly.getvalue())
          {
          failure = b.simobj.posteriormode(pathgraphs,false);
          }
        else
          failure = b.simobj.simulate(pathgraphs,b.setseed.getvalue(),
          b.computemodeforstartingvalues);
        }

      if (!failure)
        b.resultsyesno = true;
      }


    if (failure)
      b.errors=true;

    } // end: if errors == false

  }


bool superbayesreg::check_errors(void)
  {

  bool err=false;


  unsigned j,k;

  for (k=0;k<equations.size();k++)
    {

    if (equations[k].distrp->errors == true)
      {
      err = true;
      for (j=0;j<equations[k].distrp->errormessages.size();j++)
        outerror(equations[k].distrp->errormessages[j]);
      }

    }


/*
  for (k=0;k<equations.size();k++)
    {
    for (i=0;i<equations[k].distrp.size();i++)
      {
      if (equations[k].distrp[i]->errors == true)
        {
        err = true;
        for (j=0;j<equations[k].distrp[i]->errormessages.size();j++)
          outerror(equations[k].distrp[i]->errormessages[j]);
        }

      }

    }
*/

  return err;
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
  unsigned weightpos=0;

  vector<ST::string> modelvarnamesh;
  modelvarnamesh = modreg.getModelVarnamesAsVector();
  modelvarnamesv.erase(modelvarnamesv.begin(),modelvarnamesv.end());
  for(i=0;i<modelvarnamesh.size();i++)
    if (modelvarnamesh[i] != "const")
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
      {
      ST::string test = notex[i];
      if (notex[i] != "const" && notex[i] != "noconst")
        {
        outerror("ERROR: variable " + notex[i] + " is not existing\n");
        failure = true;
        }
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
    computemodeforstartingvalues = true;

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif

    distr_gaussians.push_back(DISTR_gaussian(aresp.getvalue(),bresp.getvalue(),
                                      &generaloptions,D.getCol(0),path,w) );

    equations[modnr].distrp = &distr_gaussians[distr_gaussians.size()-1];
    equations[modnr].pathd = outfile.getvalue() + "_scale.res";

    }
//-------------------------- END: Gaussian response ----------------------------

//------------- variance equation of heteroscedastic Gaussian ------------------
  if (family.getvalue() == "hetgaussian" && equationtype.getvalue()=="variance")
    {

    distr_vargaussians.push_back(DISTR_vargaussian(&generaloptions,D.getCol(0)));

    equations[modnr].distrp = &distr_vargaussians[distr_vargaussians.size()-1];
    equations[modnr].pathd = "";

    }
//---------- END: variance equation of heteroscedastic Gaussian ----------------

//---------------------- heteroscedastic Gaussian response ---------------------
  else if (family.getvalue() == "hetgaussian" && equationtype.getvalue()=="mean")
    {

    computemodeforstartingvalues = true;

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif

    distr_hetgaussians.push_back(DISTR_hetgaussian(aresp.getvalue(),
                                 bresp.getvalue(), &generaloptions,
                                 D.getCol(0),path,scaleconst.getvalue(),w) );

    equations[modnr].distrp = &distr_hetgaussians[distr_hetgaussians.size()-1];

#if defined(__BUILDING_LINUX)
    equations[modnr].pathd = defaultpath + "/temp/" + name  + "_scale.res";
#else
    equations[modnr].pathd = defaultpath + "\\temp\\" + name  + "_scale.res";
#endif

    if (distr_vargaussians.size() != 1)
      {
      outerror("ERROR: Variance equation for heteroscedastic Gausian responses is missing");
      return true;
      }
    else
      {
      distr_vargaussians[0].dgaussian = &distr_hetgaussians[distr_hetgaussians.size()-1];
      }

    }
//----------------- END: heteroscedastic Gaussian response ---------------------

//---------------------------------- ZINB pi -----------------------------------
  else if (family.getvalue() == "zinb_pi" && equationtype.getvalue()=="pi")
    {

    computemodeforstartingvalues = true;

    distr_negbinzip_pis.push_back(DISTR_negbinzip_pi(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_negbinzip_pis[distr_negbinzip_pis.size()-1];
    equations[modnr].pathd = "";

    }
//------------------------------- END: ZINB pi ---------------------------------

//---------------------------------- ZINB delta --------------------------------
  else if (family.getvalue() == "zinb_delta" && equationtype.getvalue()=="delta")
    {

    computemodeforstartingvalues = true;

    int strmax = stoprmax.getvalue();
    double flimit = fraclimit.getvalue();

    distr_negbinzip_deltas.push_back(DISTR_negbinzip_delta(&generaloptions,
                                     D.getCol(0),flimit,
                                     strmax,w));

    equations[modnr].distrp = &distr_negbinzip_deltas[distr_negbinzip_deltas.size()-1];
    equations[modnr].pathd = "";

    }
//------------------------------ END: ZINB delta -------------------------------


//----------------------------------- ZINB mu ----------------------------------
  else if (family.getvalue() == "zinb_mu" && equationtype.getvalue()=="mean")
    {

    computemodeforstartingvalues = true;

    distr_negbinzip_mus.push_back(DISTR_negbinzip_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_negbinzip_mus[distr_negbinzip_mus.size()-1];
    equations[modnr].pathd = "";

    if (distr_negbinzip_pis.size() != 1)
      {
      outerror("ERROR: Equation for pi is missing");
      return true;
      }

    if (distr_negbinzip_deltas.size() != 1)
      {
      outerror("ERROR: Equation for delta is missing");
      return true;
      }

    predict_mult_distrs.push_back(&distr_negbinzip_deltas[distr_negbinzip_deltas.size()-1]);
    predict_mult_distrs.push_back(&distr_negbinzip_pis[distr_negbinzip_pis.size()-1]);
    predict_mult_distrs.push_back(&distr_negbinzip_mus[distr_negbinzip_mus.size()-1]);

    distr_negbinzip_pis[distr_negbinzip_pis.size()-1].distrmu =
    &distr_negbinzip_mus[distr_negbinzip_mus.size()-1];

    distr_negbinzip_pis[distr_negbinzip_pis.size()-1].distrdelta =
    &distr_negbinzip_deltas[distr_negbinzip_deltas.size()-1];


    distr_negbinzip_mus[distr_negbinzip_mus.size()-1].distrpi =
    &distr_negbinzip_pis[distr_negbinzip_pis.size()-1];

    distr_negbinzip_mus[distr_negbinzip_mus.size()-1].distrdelta =
    &distr_negbinzip_deltas[distr_negbinzip_deltas.size()-1];


    distr_negbinzip_deltas[distr_negbinzip_deltas.size()-1].distrpi =
    &distr_negbinzip_pis[distr_negbinzip_pis.size()-1];

    distr_negbinzip_deltas[distr_negbinzip_deltas.size()-1].distrmu =
    &distr_negbinzip_mus[distr_negbinzip_mus.size()-1];

    }
//------------------------------ END: ZINB mu ----------------------------------


//-------------------------------- negbin delta --------------------------------
  else if (family.getvalue() == "negbin_delta" && equationtype.getvalue()=="delta")
    {

    computemodeforstartingvalues = true;

    double flimit = fraclimit.getvalue();
    int strmax = stoprmax.getvalue();
    bool sl = slow.getvalue();

    distr_negbin_deltas.push_back(DISTR_negbin_delta(&generaloptions,D.getCol(0),
    flimit,strmax,sl,w));

    equations[modnr].distrp = &distr_negbin_deltas[distr_negbin_deltas.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_negbin_deltas[distr_negbin_deltas.size()-1]);

    }
//---------------------------- END: negbin delta -------------------------------

//------------------------------- negbin mu ------------------------------------
  else if (family.getvalue() == "negbin_mu" && equationtype.getvalue()=="mean")
    {

    computemodeforstartingvalues = true;

    distr_negbin_mus.push_back(DISTR_negbin_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_negbin_mus[distr_negbin_mus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_negbin_mus[distr_negbin_mus.size()-1]);

    if (distr_negbin_deltas.size() != 1)
      {
      outerror("ERROR: Equation for delta is missing");
      return true;
      }
    else
      {
      distr_negbin_deltas[distr_negbin_deltas.size()-1].distrp.push_back
      (&distr_negbin_mus[distr_negbin_mus.size()-1]);

      distr_negbin_mus[distr_negbin_mus.size()-1].distrp.push_back
      (&distr_negbin_deltas[distr_negbin_deltas.size()-1]);
      }

    }
//------------------------------- END: negbin mu -------------------------------

//-------------------------------- beta sigma2 ---------------------------------
  else if (family.getvalue() == "beta_sigma2" && equationtype.getvalue()=="sigma2")
    {

    computemodeforstartingvalues = true;

    distr_beta_sigma2s.push_back(DISTR_beta_sigma2(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_beta_sigma2s[distr_beta_sigma2s.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_beta_sigma2s[distr_beta_sigma2s.size()-1]);

    }
//---------------------------- END: beta sigma2 -------------------------------

//------------------------------- beta mu ------------------------------------
  else if (family.getvalue() == "beta_mu" && equationtype.getvalue()=="mean")
    {

    computemodeforstartingvalues = true;

    distr_beta_mus.push_back(DISTR_beta_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_beta_mus[distr_beta_mus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_beta_mus[distr_beta_mus.size()-1]);

    if (distr_beta_sigma2s.size() != 1)
      {
      outerror("ERROR: Equation for sigma2 is missing");
      return true;
      }
    else
      {
      distr_beta_sigma2s[distr_beta_sigma2s.size()-1].distrp.push_back
      (&distr_beta_mus[distr_beta_mus.size()-1]);

      distr_beta_mus[distr_beta_mus.size()-1].distrp.push_back
      (&distr_beta_sigma2s[distr_beta_sigma2s.size()-1]);
      }

    }
//------------------------------- END: beta mu -------------------------------


//---------------------------------- ZIP pi ------------------------------------
  else if (family.getvalue() == "zip_pi" && equationtype.getvalue()=="pi")
    {

    computemodeforstartingvalues = true;

    distr_zippis.push_back(DISTR_zippi(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_zippis[distr_zippis.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_zippis[distr_zippis.size()-1]);

    }
//------------------------------- END: ZIP pi ----------------------------------

//-------------------------------- ZIP lambda ----------------------------------
  else if (family.getvalue() == "zip_lambda" && equationtype.getvalue()=="mean")
    {

    computemodeforstartingvalues = true;

    distr_ziplambdas.push_back(DISTR_ziplambda(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_ziplambdas[distr_ziplambdas.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_ziplambdas[distr_ziplambdas.size()-1]);

    if (distr_zippis.size() != 1)
      {
      outerror("ERROR: Equation for pi is missing");
      return true;
      }
    else
      {
      distr_zippis[distr_zippis.size()-1].distrlambda =
      &distr_ziplambdas[distr_ziplambdas.size()-1];

      distr_ziplambdas[distr_ziplambdas.size()-1].distrpi =
      &distr_zippis[distr_zippis.size()-1];

      }

    }
//------------------------------- END: ZIP lambda ------------------------------


//-------------------------------- ZIP pi cloglog ------------------------------
  else if (family.getvalue() == "zip_pi_cloglog" && equationtype.getvalue()=="pi")
    {

    computemodeforstartingvalues = true;

    distr_zip_cloglog_pis.push_back(DISTR_zip_cloglog_pi(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_zip_cloglog_pis[distr_zip_cloglog_pis.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_zip_cloglog_pis[distr_zip_cloglog_pis.size()-1]);

    }
//------------------------------ END: ZIP pi cloglog----------------------------


//------------------------------ ZIP lambda cloglog ----------------------------
  else if (family.getvalue() == "zip_lambda_cloglog" && equationtype.getvalue()=="mean")
    {

    computemodeforstartingvalues = true;

    distr_zip_cloglog_mus.push_back(DISTR_zip_cloglog_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_zip_cloglog_mus[distr_zip_cloglog_mus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_zip_cloglog_mus[distr_zip_cloglog_mus.size()-1]);

    if (distr_zip_cloglog_pis.size() != 1)
      {
      outerror("ERROR: Equation for pi is missing");
      return true;
      }
    else
      {
      distr_zip_cloglog_pis[distr_zip_cloglog_pis.size()-1].distrp.push_back
      (&distr_zip_cloglog_mus[distr_zip_cloglog_mus.size()-1]);

      distr_zip_cloglog_mus[distr_zip_cloglog_mus.size()-1].distrp.push_back
      (&distr_zip_cloglog_pis[distr_zip_cloglog_pis.size()-1]);
      }

    }
//---------------------------- END: ZIP lambda cloglog -------------------------


//----------------------- multinomial probit response --------------------------
  else if (family.getvalue() == "multinom_probit" && equationtype.getvalue()=="meanservant")
    {

    computemodeforstartingvalues = true;
    unsigned cat = distr_multinomprobits.size();

    distr_multinomprobits.push_back(DISTR_multinomprobit(&generaloptions,cat,false,D.getCol(0)));

    equations[modnr].distrp = &distr_multinomprobits[distr_multinomprobits.size()-1];
    equations[modnr].pathd = "";

    }
  else if (family.getvalue() == "multinom_probit" && equationtype.getvalue()=="mean")
    {

    computemodeforstartingvalues = true;
    unsigned cat = distr_multinomprobits.size();

    distr_multinomprobits.push_back(DISTR_multinomprobit(&generaloptions,cat,true,D.getCol(0),w) );

    equations[modnr].distrp = &distr_multinomprobits[distr_multinomprobits.size()-1];
    equations[modnr].pathd = "";

    if (distr_multinomprobits.size() > 1)
      {
      unsigned i;
      for (i=0;i<distr_multinomprobits.size()-1;i++)
        {
        distr_multinomprobits[distr_multinomprobits.size()-1].assign_othercat(&distr_multinomprobits[i]);
        }

      }

    distr_multinomprobits[distr_multinomprobits.size()-1].create_responsecat();

    }
//--------------------- END: multinomial probit response -----------------------

//----------------------- multivariate gaussian response -----------------------
  else if (family.getvalue() == "multgaussian" && equationtype.getvalue()=="meanservant")
    {

    computemodeforstartingvalues = true;
    unsigned cat = distr_multgaussians.size();

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif

    distr_multgaussians.push_back(DISTR_multgaussian(aresp.getvalue(),
    bresp.getvalue(),cat,&generaloptions,false,D.getCol(0),path) );

    equations[modnr].distrp = &distr_multgaussians[distr_multgaussians.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_multgaussians[distr_multgaussians.size()-1]);

    }
  else if (family.getvalue() == "multgaussian" && equationtype.getvalue()=="mean")
    {

    computemodeforstartingvalues = true;
    unsigned cat = distr_multgaussians.size();

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif


    distr_multgaussians.push_back(DISTR_multgaussian(aresp.getvalue(),
    bresp.getvalue(),cat,&generaloptions,true,D.getCol(0),path,w) );

    equations[modnr].distrp = &distr_multgaussians[distr_multgaussians.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_multgaussians[distr_multgaussians.size()-1]);

    unsigned i;
    vector<DISTR *> dp;
    for (i=0;i<distr_multgaussians.size();i++)
      {
      dp.push_back(&distr_multgaussians[i]);
      }

    for (i=0;i<distr_multgaussians.size();i++)
      distr_multgaussians[i].assign_distributions(dp);

    }
//--------------------- END: multivariate gaussian response --------------------


//----------------------- multinomial logit response ---------------------------
  else if (family.getvalue() == "multinom_logit" && equationtype.getvalue()=="meanservant")
    {

    computemodeforstartingvalues = true;
    unsigned cat = distr_multinomlogits.size();

    distr_multinomlogits.push_back(DISTR_multinomlogit(&generaloptions,cat,false,D.getCol(0)));

    equations[modnr].distrp = &distr_multinomlogits[distr_multinomlogits.size()-1];
    equations[modnr].pathd = "";

    }
  else if (family.getvalue() == "multinom_logit" && equationtype.getvalue()=="mean")
    {

    computemodeforstartingvalues = true;
    unsigned cat = distr_multinomlogits.size();

    distr_multinomlogits.push_back(DISTR_multinomlogit(&generaloptions,cat,true,D.getCol(0),w) );

    equations[modnr].distrp = &distr_multinomlogits[distr_multinomlogits.size()-1];
    equations[modnr].pathd = "";

    if (distr_multinomlogits.size() > 1)
      {
      unsigned i;
      for (i=0;i<distr_multinomlogits.size()-1;i++)
        {
        distr_multinomlogits[distr_multinomlogits.size()-1].assign_othercat(&distr_multinomlogits[i]);
        }

      distr_multinomlogits[distr_multinomlogits.size()-1].create_responsecat();

      }

    }
//--------------------- END: multinomial logit response ------------------------


//----------------------------- Poisson response -------------------------------
  else if (family.getvalue() == "poisson")
    {
    computemodeforstartingvalues = true;

    distr_poissons.push_back(DISTR_poisson(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_poissons[distr_poissons.size()-1];
    equations[modnr].pathd = "";

    }
//-------------------------- END: poisson response -----------------------------

//-------------- Poisson response with extended response function---------------
  else if (family.getvalue() == "poisson_ext")
    {
    computemodeforstartingvalues = true;

    distr_poisson_exts.push_back(DISTR_poisson_ext(
    &generaloptions,D.getCol(0),aexp.getvalue(),
    bexp.getvalue(),adaptexp.getvalue(),w));

    equations[modnr].distrp = &distr_poisson_exts[distr_poisson_exts.size()-1];
    equations[modnr].pathd = "";

    }
//--------------End: Poisson response with extended response function-----------

//---- Poisson response with extended response function, linear for eta > 0-----
  else if (family.getvalue() == "poisson_extlin")
    {
    computemodeforstartingvalues = true;

    distr_poisson_extlins.push_back(DISTR_poisson_extlin(
    &generaloptions,D.getCol(0),w));

    equations[modnr].distrp =
    &distr_poisson_extlins[distr_poisson_extlins.size()-1];
    equations[modnr].pathd = "";

    }
//-- End: Poisson response with extended response function, linear for eta > 0--

//-------------------------- log-Gaussian response -----------------------------
  else if (family.getvalue() == "loggaussian")
    {

    computemodeforstartingvalues = true;

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif

    distr_loggaussians.push_back(DISTR_loggaussian(aresp.getvalue(),
                                       bresp.getvalue(),
                                      &generaloptions,D.getCol(0),path,w) );

    equations[modnr].distrp = &distr_loggaussians[distr_loggaussians.size()-1];
    equations[modnr].pathd =  outfile.getvalue() + "_scale.res";

    }
//------------------------- END: log-Gaussian response -------------------------

//----------------------------- quantreg response ------------------------------
  else if (family.getvalue() == "quantreg")
    {

    computemodeforstartingvalues = true;

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif

    double quant = quantile.getvalue();
    distr_quantregs.push_back(DISTR_quantreg(aresp.getvalue(),
                                              bresp.getvalue(),
                                     &generaloptions,D.getCol(0),path,quant, w) );

    equations[modnr].distrp = &distr_quantregs[distr_quantregs.size()-1];
    equations[modnr].pathd = outfile.getvalue() + "_scale.res";

    }
//------------------------- END: quantreg response -----------------------------

//-------------------------- gaussian_mixture response -------------------------
  else if (family.getvalue() == "gaussian_mixture")
    {

    computemodeforstartingvalues = true;

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif

    distr_gaussianmixtures.push_back(DISTR_gaussianmixture(aresp.getvalue(),
                                              bresp.getvalue(),
                                     &generaloptions,D.getCol(0),path,w) );

    equations[modnr].distrp = &distr_gaussianmixtures[distr_gaussianmixtures.size()-1];
    equations[modnr].pathd = outfile.getvalue()  + "_scale.res";

    }
//------------------------- END: gaussian mixture response ---------------------

//--------------- Gaussian response, exponential response function -------------
  else if (family.getvalue() == "gaussian_exp")
    {

    computemodeforstartingvalues = true;

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif

    distr_gaussian_exps.push_back(DISTR_gaussian_exp(
                                  aresp.getvalue(),bresp.getvalue(),
                                  &generaloptions,D.getCol(0),path,w) );

    equations[modnr].distrp = &distr_gaussian_exps[distr_gaussian_exps.size()-1];
    equations[modnr].pathd =  outfile.getvalue() +  "_scale.res";

    }
//------------- END: Gaussian response, exponential response function ----------
//------------- Gaussian response, multiplicative random effects allowed -------
  else if (family.getvalue() == "gaussian_mult")
    {

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif

    distr_gaussian_mults.push_back(DISTR_gaussian_mult(
                                  aresp.getvalue(),bresp.getvalue(),
                                  &generaloptions,D.getCol(0),path,w) );

    equations[modnr].distrp = &distr_gaussian_mults[distr_gaussian_mults.size()-1];

#if defined(__BUILDING_LINUX)
    equations[modnr].pathd = defaultpath + "/temp/" + name  + "_scale.res";
#else
    equations[modnr].pathd = defaultpath + "\\temp\\" + name  + "_scale.res";
#endif

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

    computemodeforstartingvalues = true;

    distr_binomials.push_back(DISTR_binomial(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_binomials[distr_binomials.size()-1];
    equations[modnr].pathd = "";

    }
//-------------------------- END: Binomial response ----------------------------
//---------------------------- Binomial response probit-------------------------
  else if (family.getvalue() == "binomial_probit")
    {

    computemodeforstartingvalues = true;

    #if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_latentutilities.raw";
    #else
    ST::string path = defaultpath + "\\temp\\" + name  + "_latentutilities.raw";
    #endif


    distr_binomialprobits.push_back(DISTR_binomialprobit(
    &generaloptions,D.getCol(0),utilities.getvalue(),path,w));

    equations[modnr].distrp = &distr_binomialprobits[distr_binomialprobits.size()-1];
    equations[modnr].pathd = outfile.getvalue() + "_latentvariables.res";
    }
//-------------------------- END: Binomial response probit ---------------------
//---------------------------- Binomial response SVM ---------------------------
  else if (family.getvalue() == "binomial_svm")
    {
    computemodeforstartingvalues = false;

    distr_binomialsvms.push_back(DISTR_binomialsvm(
    &generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_binomialsvms[distr_binomialsvms.size()-1];
    equations[modnr].pathd = "";

    }
//-------------------------- END: Binomial response SVM ------------------------
//----------------------- Binomial response logit fruehwirth--------------------
  else if (family.getvalue() == "binomial_logit_l1")
    {

    computemodeforstartingvalues = true;

    distr_logit_fruehwirths.push_back(DISTR_logit_fruehwirth(H.getvalue(),
    &generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_logit_fruehwirths[distr_logit_fruehwirths.size()-1];
    equations[modnr].pathd = "";

    }
//-------------------------- END: Binomial response probit ---------------------

  equations[modnr].distrp->responsename=rname;
  equations[modnr].distrp->weightname=wn;

  if (equations[modnr].hlevel==1)
    master.level1_likep.push_back(equations[modnr].distrp);

  return false;

  }


void superbayesreg::extract_data(unsigned i, datamatrix & d,datamatrix & iv,
                                 unsigned dim_dm)
  {

  d = datamatrix(D.rows(),dim_dm);

  int j1,j2;
  unsigned p;
  for (p=0;p<dim_dm;p++)
    {
    j1 = terms[i].varnames[terms[i].varnames.size()-1-p].isinlist(modelvarnamesv);
    d.putCol(dim_dm-1-p,D.getCol(j1));
    }

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\d.res");
  // d.prettyPrint(out);
  // TEST

  if (terms[i].varnames.size() > dim_dm)
    {
    j2 = terms[i].varnames[0].isinlist(modelvarnamesv);
    iv = D.getCol(j2);
    }
  }


bool superbayesreg::create_predict(void)
  {
  if (predict.getvalue() != "no")
    {

    unsigned modnr = equations.size()-1;

    ST::string h = equations[modnr].paths;

#if defined(__BUILDING_LINUX)
    ST::string pathnonp = defaultpath + "/temp/" + name + "_" + h +
                            "_predict.raw";
#else
    ST::string pathnonp = defaultpath + "\\temp\\" + name + "_" + h +
                            "_predict.raw";
#endif

#if defined(__BUILDING_LINUX)
    ST::string pathnonp2 = defaultpath + "/temp/" + name + "_" + h +
                            "_deviance.raw";
#else
    ST::string pathnonp2 = defaultpath + "\\temp\\" + name + "_" + h +
                            "_deviance.raw";
#endif


    ST::string pathres = outfile.getvalue() +  "_" + h + "_predict.res";


    if ((predict.getvalue() == "full") || (predict.getvalue() == "fulls"))
      {

      if (equations[modnr].distrp->maindistribution == false)
        {
        outerror("ERROR: predict=full is allowed only for the main regression equation\n");
        outerror("HINT: Use predict=predictor to get predictions of the predictor for this equation\n");
        return true;
        }
      else
        {
        if (equations[modnr].distrp->predict_mult == false)
          {
          FC_predicts.push_back(FC_predict(&generaloptions,
                             equations[modnr].distrp,"",pathnonp,
                             pathnonp2,D,modelvarnamesv));

          if (predict.getvalue() == "fulls")
            FC_predicts[FC_predicts.size()-1].nosamples=false;

          if (mse.getvalue() ==  "yes")
            FC_predicts[FC_predicts.size()-1].MSE = MCMC::quadraticMSE;

          if (mse.getvalue() ==  "quadratic")
            FC_predicts[FC_predicts.size()-1].MSE = MCMC::quadraticMSE;

          if (mse.getvalue() ==  "check")
            {
            FC_predicts[FC_predicts.size()-1].MSE = MCMC::checkMSE;
            FC_predicts[FC_predicts.size()-1].MSEparam = mseparam.getvalue();
            }

          equations[modnr].add_FC(&FC_predicts[FC_predicts.size()-1],pathres);
          }
        else // predict_mult==true
          {

          FC_predicts_mult.push_back(FC_predict_mult(&generaloptions,
                             predict_mult_distrs,"",pathnonp,
                             pathnonp2,D,modelvarnamesv));

          if (predict.getvalue() == "fulls")
            FC_predicts_mult[FC_predicts_mult.size()-1].nosamples=false;

          equations[modnr].add_FC(&FC_predicts_mult[FC_predicts_mult.size()-1],pathres);


          }

        }
      }
    else if (predict.getvalue() == "predictor")
      {

      FC_predict_predictors.push_back(FC_predict_predictor(&generaloptions,
                           equations[modnr].distrp,"",pathnonp,pathnonp2,D,
                           modelvarnamesv));

      equations[modnr].add_FC(&FC_predict_predictors[FC_predict_predictors.size()-1],pathres);

      }


    }

  return false;
  }


void superbayesreg::create_predictive_check(void)
  {

  if (pred_check.getvalue() == true)
    {

    unsigned modnr = equations.size()-1;

    ST::string h = equations[modnr].paths;

    ST::string pathres = outfile.getvalue() +  "_" + h + "_predictive_check.res";

    FC_predictive_checks.push_back(FC_predictive_check(&generaloptions,
                         equations[modnr].distrp,"","",D,modelvarnamesv));

    equations[modnr].add_FC(&FC_predictive_checks[FC_predictive_checks.size()-1]
                           ,pathres);

    }

  }


void superbayesreg::create_cv(void)
  {

  if (cv.getvalue() == true)
    {

    unsigned modnr = equations.size()-1;

    ST::string h = equations[modnr].paths;


#if defined(__BUILDING_LINUX)
    ST::string pathnonp = defaultpath + "/temp/" + name + "_" + h +
                            "_cv_responses.raw";
#else
    ST::string pathnonp = defaultpath + "\\temp\\" + name + "_" + h +
                            "_cv_responses.raw";
#endif


    ST::string pathres = outfile.getvalue() +  "_" + h + "_cv.res";

    FCcv = FC_cv(&generaloptions,equations[modnr].distrp,"",pathnonp,&FC_hrandoms);

    equations[modnr].add_FC(&FCcv,pathres);

    }

  }

void superbayesreg::create_offset(unsigned i)
  {

  datamatrix d,iv;
  extract_data(i,d,iv,1);

  unsigned modnr = equations.size()-1;

  if (equations[modnr].distrp->linpred_current==1)
    equations[modnr].distrp->linearpred1.plus(d);
  else
    equations[modnr].distrp->linearpred1.plus(d);

  equations[modnr].distrp->offsetname = terms[i].varnames[0];

  }

bool superbayesreg::create_linear(void)
  {

  unsigned modnr = equations.size()-1;

  unsigned i;
  int j;

  vector<ST::string> varnames;
  vector<ST::string> varnamesh =  lineareffects.get_constvariables(terms);



  for(i=0;i<varnamesh.size();i++)
    varnames.push_back(varnamesh[i]);

  unsigned nr = varnames.size();


  ST::string title;
  ST::string pathconst;
  ST::string pathconstres;

  ST::string h = equations[modnr].paths;

  title = h + ": linear effects" ;

#if defined(__BUILDING_LINUX)
  pathconst = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                           + "_LinearEffects"  +
                           "_" + h + ".raw";
#else
  pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                           + "_LinearEffects"  +
                           "_" + h + ".raw";
#endif

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

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\X.res");
  // X.prettyPrint(out);
  // TEST


  FC_linears.push_back(FC_linear(&master,nrlevel1,&generaloptions,equations[modnr].distrp,X,
                         varnames,title,pathconst,
                         centerlinear.getvalue()));

  equations[modnr].add_FC(&FC_linears[FC_linears.size()-1],pathconstres);

  return false;

  }


void superbayesreg::create_pspline(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_pspline.raw","nonlinear_pspline_effect_of"," Nonlinear effect (P-spline) of ");

  datamatrix d,iv;
  extract_data(i,d,iv,1);

  // TEST
  //  ofstream out("c:\\bayesx\\testh\\results\\d.res");
  //  d.prettyPrint(out);
  // TEST

  design_psplines.push_back(DESIGN_pspline(d,iv,&generaloptions,equations[modnr].distrp,
                            &FC_linears[FC_linears.size()-1],
                            terms[i].options,terms[i].varnames));

  FC_nonps.push_back(FC_nonp(&master,nrlevel1,&generaloptions,equations[modnr].distrp,title,
                     pathnonp,&design_psplines[design_psplines.size()-1],
                     terms[i].options,terms[i].varnames));


  equations[modnr].add_FC(&FC_nonps[FC_nonps.size()-1],pathres);


  // variances

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_pspline_var.raw","variance_of_nonlinear_pspline_effect_of","Variance of nonlinear effect of ");

  FC_nonp_variances.push_back(FC_nonp_variance(&master,nrlevel1,&generaloptions,equations[modnr].distrp,
                                title,pathnonp,&design_psplines[design_psplines.size()-1],
                                &FC_nonps[FC_nonps.size()-1],terms[i].options,
                                terms[i].varnames));

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
  extract_data(i,d,iv,1);

  unsigned fnr;
  ST::string na = terms[i].varnames[terms[i].varnames.size()-1];
  bool found = findREdistr(na,equations[modnr],fnr);

  if (found==false)
    {
    outerror("ERROR: level 2 equation for variable " + na + " not found \n");
    return true;
    }


  design_hrandoms.push_back(DESIGN_hrandom(d,iv,&generaloptions,
                            equations[modnr].distrp,
                            &FC_linears[FC_linears.size()-1],
                             equations[fnr].distrp,
                            terms[i].options,terms[i].varnames));

  FC_hrandoms.push_back(FC_hrandom(&master,nrlevel1,&generaloptions,equations[modnr].distrp,
                        equations[fnr].distrp, title,pathnonp,pathnonp2,
                        &design_hrandoms[design_hrandoms.size()-1],
                        terms[i].options,terms[i].varnames));

  equations[modnr].add_FC(&FC_hrandoms[FC_hrandoms.size()-1],pathres);

  // variances

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_hrandom_var.raw","variance_of_random_effect_of","Variance of random effect of ");

  if (terms[i].options[35] == "iid")
    {
    FC_hrandom_variances.push_back(FC_hrandom_variance(&master,nrlevel1,
                                   &generaloptions,equations[modnr].distrp,
                                   equations[fnr].distrp,
                                  title,pathnonp,
                                  &design_hrandoms[design_hrandoms.size()-1],
                                  &FC_hrandoms[FC_hrandoms.size()-1],
                                  terms[i].options, terms[i].varnames));

    equations[modnr].add_FC(&FC_hrandom_variances[FC_hrandom_variances.size()-1],pathres);

    }
  else if (terms[i].options[35] == "lasso")
    {

    FC_hrandom_variance_vecs.push_back(FC_hrandom_variance_vec(&master,nrlevel1,
                                   &generaloptions,equations[modnr].distrp,
                                   equations[fnr].distrp,
                                  title,pathnonp,
                                  &design_hrandoms[design_hrandoms.size()-1],
                                  &FC_hrandoms[FC_hrandoms.size()-1],
                                  terms[i].options, terms[i].varnames));

    equations[modnr].add_FC(&FC_hrandom_variance_vecs[FC_hrandom_variance_vecs.size()-1],pathres);

    }
  else if (terms[i].options[35] == "nmig")
    {

    FC_hrandom_variance_vec_nmigs.push_back(FC_hrandom_variance_vec_nmig(&master,
                                        nrlevel1,&generaloptions,equations[modnr].distrp,
                                        equations[fnr].distrp,
                                        title,pathnonp,
                                        &design_hrandoms[design_hrandoms.size()-1],
                                        &FC_hrandoms[FC_hrandoms.size()-1],
                                        terms[i].options, terms[i].varnames));

    equations[modnr].add_FC(&FC_hrandom_variance_vec_nmigs[FC_hrandom_variance_vec_nmigs.size()-1],pathres);

    }
  else if ((terms[i].options[35] == "ssvs"))
    {

    FC_hrandom_variance_ssvss.push_back(FC_hrandom_variance_ssvs(&master,nrlevel1,
                                        &generaloptions,equations[modnr].distrp,
                                        equations[fnr].distrp,
                                        title,pathnonp,
                                        &design_hrandoms[design_hrandoms.size()-1],
                                        &FC_hrandoms[FC_hrandoms.size()-1],
                                        terms[i].options, terms[i].varnames));

    equations[modnr].add_FC(&FC_hrandom_variance_ssvss[FC_hrandom_variance_ssvss.size()-1],pathres);

    }

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
  dp_pspline->meaneffectintvar = 1;

  FC_mults.push_back(FC_mult(true,multexp));
  equations[modnr].add_FC(&FC_mults[FC_mults.size()-1],"");

  term helpt = terms[i];
  terms[i].varnames.erase(terms[i].varnames.begin(),terms[i].varnames.end());
  terms[i].varnames.push_back(helpt.varnames[1]);
  terms[i].varnames.push_back(helpt.varnames[0]);

  terms[i].options[12] = "true";
  bool error=false;
  error = create_hrandom(i);
  if (error==true)
    return true;
  FC_nonp * fcnp_hrandom = &FC_hrandoms[FC_hrandoms.size()-1];
  MCMC::DESIGN * dp_hrandom = &design_hrandoms[design_hrandoms.size()-1];
  dp_hrandom->changingdesign=true;
  dp_hrandom->meaneffectintvar=0;

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

  bool me=false ;
  if (terms[i].options[19] == "true")
    me = true;

  double mec;
  int f;
   f = (terms[i].options[34]).strtodouble(mec);


  FC_mults[FC_mults.size()-1].set_multeffects(&master,nrlevel1,&generaloptions,
                                              title,pathnonp,samplem,me,mec);

  equations[modnr].add_FC(&FC_mults[FC_mults.size()-1],pathres);

  return false;

  }


bool  superbayesreg::find_map(unsigned i,MAP::map & m)
    {

    mapobject * mapp;

    int objpos = findstatobject(*statobj,terms[i].options[8],"map");

    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      mapp = dynamic_cast<mapobject*>(s);
      m = mapp->getmap();
      return false;
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

    }


bool superbayesreg::create_mrf(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_spatial.raw","spatial_MRF_effect_of","Spatial effect (MRF) of ");

  datamatrix d,iv;
  extract_data(i,d,iv,1);

  /*
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
  */

  MAP::map m;
  if (find_map(i,m) == true)
    return true;

  bool isconnected = m.isconnected();
  if (isconnected==false)
    {
    outerror("ERROR: map is disconnected, spatial effect cannot be estimated\n");
    return true;
    }

  design_mrfs.push_back(DESIGN_mrf(d,iv,&generaloptions,equations[modnr].distrp,
                           &FC_linears[FC_linears.size()-1],m,
                            terms[i].options,terms[i].varnames));

  FC_nonps.push_back(FC_nonp(&master,nrlevel1,&generaloptions,equations[modnr].distrp,title,
                     pathnonp,&design_mrfs[design_mrfs.size()-1],
                     terms[i].options,terms[i].varnames));

  if (FC_nonps[FC_nonps.size()-1].errors==true)
    return true;


  equations[modnr].add_FC(&FC_nonps[FC_nonps.size()-1],pathres);

  // variances

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_spatial_var.raw","variance_of_spatial_MRF_effect_of","Variance of spatial effect of ");

  FC_nonp_variances.push_back(FC_nonp_variance(&master,nrlevel1,
                                &generaloptions,equations[modnr].distrp,
                                title,pathnonp,&design_mrfs[design_mrfs.size()-1],
                                &FC_nonps[FC_nonps.size()-1],terms[i].options,
                                terms[i].varnames));

  equations[modnr].add_FC(&FC_nonp_variances[FC_nonp_variances.size()-1],pathres);

  return false;
  }


bool superbayesreg::create_kriging(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_kriging.raw","2dim_kriging_effect_of","2dim effect (kriging) of ");

  datamatrix d,iv;
  extract_data(i,d,iv,2);

  //----------------------------------------------------------------------------

  datamatrix knotdata;
  if (terms[i].options[37] != "")
    {
    dataobject * datap;                           // pointer to datasetobject
    int objpos = findstatobject(*statobj,terms[i].options[37],"dataset");
    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      if (datap->obs()==0 || datap->getVarnames().size()==0)
        {
        outerror("ERROR: dataset object " + terms[i].options[37] + " does not contain any data\n");
        return true;
        }
      else if (datap->getVarnames().size()>2)
        {
        outerror("ERROR: dataset object " + terms[i].options[37] + " contains more than two variables\n");
        return true;
        }
      }
    else
      {
      outerror("ERROR: dataset object " + terms[i].options[37] + " is not existing\n");
      return true;
      }
    list<ST::string> knotnames = datap->getVarnames();
    ST::string expr = "";
    datap->makematrix(knotnames,knotdata,expr);

    }

  //----------------------------------------------------------------------------

  design_krigings.push_back(DESIGN_kriging(d,iv,&generaloptions,
                            equations[modnr].distrp,
                           &FC_linears[FC_linears.size()-1],
                            terms[i].options,terms[i].varnames,knotdata));

  FC_nonps.push_back(FC_nonp(&master,nrlevel1,&generaloptions,equations[modnr].distrp,
                     title, pathnonp,&design_krigings[design_krigings.size()-1],
                     terms[i].options,terms[i].varnames));

  equations[modnr].add_FC(&FC_nonps[FC_nonps.size()-1],pathres);

  // variances

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_kriging_var.raw","variance_of_2dim_kriging_effect_of",
  "Variance of 2dim effect of ");

  FC_nonp_variances.push_back(FC_nonp_variance(&master,nrlevel1,
                                &generaloptions,equations[modnr].distrp,
                                title,pathnonp,&design_krigings[design_krigings.size()-1],
                                &FC_nonps[FC_nonps.size()-1],terms[i].options,
                                terms[i].varnames));

  equations[modnr].add_FC(&FC_nonp_variances[FC_nonp_variances.size()-1],pathres);

  return false;
  }


bool superbayesreg::create_geokriging(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_geokriging.raw","geokriging_effect_of",
             "spatial effect (kriging) of ");

  datamatrix d,iv;
  extract_data(i,d,iv,1);

  /*
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
        outerror("ERROR: map object must be specified to estimate a geokriging effect\n");
      else
        outerror("ERROR: map object " + terms[i].options[1] + " is not existing\n");
      }
    else
      outerror("ERROR: " + terms[i].options[1] + " is not a map object\n");
    return true;
    }

  MAP::map m = mapp->getmap();
  */

  MAP::map m;
  if (find_map(i,m) == true)
    return true;

  design_krigings.push_back(DESIGN_kriging(d,iv,m,&generaloptions,
                            equations[modnr].distrp,
                           &FC_linears[FC_linears.size()-1],
                            terms[i].options,terms[i].varnames));

  FC_nonps.push_back(FC_nonp(&master,nrlevel1,&generaloptions,equations[modnr].distrp,
                     title, pathnonp,&design_krigings[design_krigings.size()-1],
                     terms[i].options,terms[i].varnames));

  equations[modnr].add_FC(&FC_nonps[FC_nonps.size()-1],pathres);

  // variances

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_geokriging_var.raw","variance_of_geokriging_effect_of",
  "Variance of geokriging effect of ");

  FC_nonp_variances.push_back(FC_nonp_variance(&master,nrlevel1,
                                &generaloptions,equations[modnr].distrp,
                                title,pathnonp,&design_krigings[design_krigings.size()-1],
                                &FC_nonps[FC_nonps.size()-1],terms[i].options,
                                terms[i].varnames));

  equations[modnr].add_FC(&FC_nonp_variances[FC_nonp_variances.size()-1],pathres);

  return false;
  }



bool superbayesreg::create_ridge_lasso(unsigned i)
  {
  bool isridge;
  if (terms[i].options[0] == "ridge")
    isridge = true;
  else
    isridge = false;

  if ( ( (ridge == -1) && (terms[i].options[0] == "ridge")  ) ||
       ( (lasso == -1) && (terms[i].options[0] == "lasso")  )
     )
    {

    unsigned modnr = equations.size()-1;

    ST::string title;
    ST::string pathpen;
    ST::string pathpenres;

    ST::string titlevar;
    ST::string pathpenvar;
    ST::string pathpenresvar;


    ST::string h = equations[modnr].paths;

    if (isridge)
      {

      title = h + ": linear effects with ridge penalty";

#if defined(__BUILDING_LINUX)
      pathpen = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                             + "_LinearEffects_ridgepenalty"  +
                             "_" + h + ".raw";
#else
      pathpen = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                             + "_LinearEffects_ridgepenalty"  +
                             "_" + h + ".raw";
#endif

      pathpenres = outfile.getvalue() + "_" + h +
                   "_LinearEffects_ridgepenalty.res";

      titlevar = h + ": linear effects with ridge penalty (var)";

#if defined(__BUILDING_LINUX)
      pathpenvar = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                             + "_LinearEffects_ridgepenalty_var"  +
                             "_" + h + ".raw";
#else
      pathpenvar = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                             + "_LinearEffects_ridgepenalty_var"  +
                             "_" + h + ".raw";
#endif

      pathpenresvar = outfile.getvalue() + "_" + h +
                   "_LinearEffects_ridgepenalty_var.res";

      }
    else
      {
      title = h + ": linear effects with lasso penalty";

#if defined(__BUILDING_LINUX)
      pathpen = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                             + "_LinearEffects_lassopenalty"  +
                             "_" + h + ".raw";
#else
      pathpen = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                             + "_LinearEffects_lassopenalty"  +
                             "_" + h + ".raw";
#endif

      pathpenres = outfile.getvalue() + "_" + h +
                   "_LinearEffects_lassopenalty.res";

      titlevar = h + ": linear effects with lasso penalty (var)";

#if defined(__BUILDING_LINUX)
      pathpenvar = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                             + "_LinearEffects_lassopenalty_var"  +
                             "_" + h + ".raw";
#else
      pathpenvar = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                             + "_LinearEffects_lassopenalty_var"  +
                             "_" + h + ".raw";
#endif

      pathpenresvar = outfile.getvalue() + "_" + h +
                   "_LinearEffects_lassopenalty_var.res";

      }


    if (pathpen.isvalidfile() == 1)
      {
      errormessages.push_back("ERROR: unable to open file " + pathpen +
                                 " for writing\n");
      return true;
      }


    datamatrix d,iv;
    extract_data(i,d,iv,1);

    FC_linear_pens.push_back(FC_linear_pen(&master,nrlevel1,&generaloptions,
                         equations[modnr].distrp,d,terms[i].varnames,title,
                         pathpen, centerlinear.getvalue()));

    equations[modnr].add_FC(&FC_linear_pens[FC_linear_pens.size()-1],pathpenres);

    FC_variance_pen_vectors.push_back(
    FC_variance_pen_vector(&master,&generaloptions,
                            &(FC_linear_pens[FC_linear_pens.size()-1])  ,
                            equations[modnr].distrp,titlevar, pathpenvar,
                            isridge));

    if (isridge)
      {
      ridge_linear = FC_linear_pens.size()-1;
      ridge = FC_variance_pen_vectors.size()-1;
      FC_variance_pen_vectors[ridge].add_variable(d,terms[i].options,
                                                  terms[i].varnames);
      }
    else
      {
      lasso_linear = FC_linear_pens.size()-1;
      lasso = FC_variance_pen_vectors.size()-1;
      FC_variance_pen_vectors[lasso].add_variable(d,terms[i].options,
                                                  terms[i].varnames);
      }

    equations[modnr].add_FC(&FC_variance_pen_vectors[
    FC_variance_pen_vectors.size()-1],pathpenresvar);

    }
  else
    {
    datamatrix d,iv;
    extract_data(i,d,iv,1);

    if (isridge)
      {
      FC_linear_pens[ridge_linear].add_variable(d,terms[i].varnames[0]);
      FC_variance_pen_vectors[ridge].add_variable(d,terms[i].options,
                                                  terms[i].varnames);
      }
    else
      {
      FC_linear_pens[lasso_linear].add_variable(d,terms[i].varnames[0]);
      FC_variance_pen_vectors[lasso].add_variable(d,terms[i].options,
                                                  terms[i].varnames);
      }

    }

  return false;
  }



bool superbayesreg::create_nonp(void)
  {

  unsigned i;
  bool error=false;

  lasso = -1;
  ridge = -1;
  lasso_linear = -1;
  ridge_linear = -1;

  for(i=0;i<terms.size();i++)
    {
    if (terms[i].options.size() > 0)
      {
      if (terms[i].options[0] == "offset")
        create_offset(i);
      if (terms[i].options[0] == "pspline")
        create_pspline(i);
      if (terms[i].options[0] == "hrandom")
        error = create_hrandom(i);
      if (terms[i].options[0] == "spatial")
        error = create_mrf(i);
      if (terms[i].options[0] == "kriging")
        error = create_kriging(i);
      if (terms[i].options[0] == "geokriging")
        error = create_geokriging(i);
      if ((terms[i].options[0] == "hrandom_pspline") ||
          (terms[i].options[0] == "hrandomexp_pspline"))
        error = create_random_pspline(i);
      if (terms[i].options[0] == "ridge")
        error = create_ridge_lasso(i);
      if (terms[i].options[0] == "lasso")
        error = create_ridge_lasso(i);
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
         {
         ST::string pathgraphs = b.outfile.getvalue();
         b.simobj.autocorr(b.maxlag.getvalue(),pathgraphs);
         }
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
    ST::string pathgraphs = b.outfile.getvalue();
    if (b.posteriormode == false)
      {
      #if defined(JAVA_OUTPUT_WINDOW)

// STEFAN: CHECKEN
// zweites Argument sollte ein Vektor sein.
      // b.simobj.get_samples(b.newcommands,b.outfile.getvalue() + "_");
      ST::string aString( b.outfile.getvalue() + "_" );
      b.simobj.get_samples(aString, b.newcommands);
      #else
      b.simobj.get_samples(pathgraphs);
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














