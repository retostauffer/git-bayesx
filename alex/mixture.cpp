#include "mixture.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//----------------- class FULLCOND_mixture --------------------------------------
//------------------------------------------------------------------------------

void FULLCOND_mixture::init_name(const ST::string & na)
    {
    char charh = '_';
    ST::string stringh = "\\_";

    FULLCOND::init_name(na);

    ST::string helpname = na.insert_string_char(charh,stringh);
    term_symbolic = "f_{" +  helpname + "}("+helpname+")";

    if (column > 0)
      priorassumptions.push_back("$" + term_symbolic + "$" +
       " (" + ST::inttostring(column+1) + ". \\mbox{ } response \\mbox{ } category)");
    else
      priorassumptions.push_back("$" + term_symbolic + "$");

    priorassumptions.push_back("i.i.d. Gaussian random effects");
    }


void FULLCOND_mixture::init_names(const vector<ST::string> & na)
    {
    char charh = '_';
    ST::string stringh = "\\_";

    FULLCOND::init_names(na);
    if (na.size()==1)
      {
      ST::string helpname = na[0].insert_string_char(charh,stringh);
      term_symbolic = "f_{" +  helpname + "}("+helpname+")";
      }
    else
      {
      ST::string helpname1 = na[0].insert_string_char(charh,stringh);
      ST::string helpname2 = na[1].insert_string_char(charh,stringh);
      term_symbolic = "f_{" +  helpname1 + "}("+helpname1+") \\cdot "
                        + helpname2;
      }

    if (column > 0)
      priorassumptions.push_back("$" + term_symbolic + "$" +
       " (" + ST::inttostring(column+1) + ". \\mbox{ } response \\mbox{ } category)");
    else
      priorassumptions.push_back("$" + term_symbolic + "$");

    priorassumptions.push_back("i.i.d. Gaussian random effects");
    }



double FULLCOND_mixture::centerbeta(void)
  {
  unsigned i;

  double sum=0;
  double * workbeta = beta.getV();

  for (i=0;i<nrpar;i++,workbeta++)
    sum+= *workbeta;

  sum /= nrpar;

  double v = sigma2/double(nrpar);

  sum = sum+sqrt(v)*rand_normal();
 
  workbeta = beta.getV();

  for (i=0;i<nrpar;i++,workbeta++)
    *workbeta-= sum;

  return sum;

  }


void FULLCOND_mixture::compute_XWX(const datamatrix & weightmat,
                                  const unsigned & col)
  {

  register unsigned j,i;

  double * workXX = XX.getV();
  int *  workindex = index.getV();
  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();
  unsigned n = posbeg.size();

  if (!randomslope)
    {
    for(j=0;j<n;j++,workXX++,++itbeg,++itend)
      {
      *workXX = 0;
      for (i=*itbeg;i<=*itend;i++,workindex++)
        *workXX += weightmat(*workindex,col);
      }
    }
  else
    {

    double * datap = data.getV();
    for(j=0;j<n;j++,workXX++,++itbeg,++itend)
      {
      *workXX = 0;
      for (i=*itbeg;i<=*itend;i++,workindex++,datap++)
        {
        *workXX += weightmat(*workindex,col) * (*datap) * (*datap);
        }
      }
    }

  }



void FULLCOND_mixture::set_lambdaconst(double la)
  {
  lambda=la;
  lambdaconst = true;
  }


FULLCOND_mixture::FULLCOND_mixture(MCMCoptions * o,DISTRIBUTION * dp,
                              FULLCOND_const * fcc,
                              const datamatrix & d, const ST::string & t,
                              const ST::string & fp,const ST::string & pr,const int & nrc,
                              const double & la, const unsigned & c)
                            : FULLCOND(o,datamatrix(1,1),t,1,1,fp)
  {

// Mixture stuff
  nrcomp = nrc;
  compweight = datamatrix(nrcomp,1,1.0/nrcomp);
  cwprior = datamatrix(nrcomp,1,1.0);
  csize = statmatrix<unsigned>(nrcomp,1,33);
    compmean = datamatrix(nrcomp,1,0);
    compvar = datamatrix(nrcomp,1,1);
    cmpriorm = 0;
    cmpriorv = 10;
    cvpriorsh = 2;
    cvpriorsc = 0.5;
// End of Mixture stuff


  fcconst = fcc;

  fctype = randomeffects;

  spatialtotal = false;
  randomslope = false;
  includefixed = false;

  changingweight = dp->get_changingweight();

  column = c;

  likep = dp;
  pathresult = pr;
  pathcurrent = pr;

  lambda = la;
  lambdaold1 = -1;
  lambdaold2 = -1;
  lambdaconst=false;

  index = statmatrix<int>(d.rows(),1);
  index2 = statmatrix<int>(d.rows(),1);
  index.indexinit();
  d.indexsort(index,0,d.rows()-1,0,0);

  unsigned j;
  int * workindex = index.getV();
  int * workindex2 = index2.getV();
  *workindex2 = *workindex;
  int help = *workindex;
  workindex++;
  workindex2++;

  for(j=1;j<d.rows();j++,workindex++,workindex2++)
    {
    *workindex2 = *workindex-help;
    help = *workindex;
    }

  posbeg = vector<unsigned>();
  posend = vector<unsigned>();

  posbeg.push_back(0);
  workindex=index.getV()+1;
  help = index(0,0);
  for(j=1;j<d.rows();j++,workindex++)
    {
    if ( d(*workindex,0) != d(help,0) )
      {
      posbeg.push_back(j);
      posend.push_back(j-1);
      }

    help = *workindex;

    }

  posend.push_back(d.rows()-1);

  effvalues = datamatrix(posbeg.size(),1);
  double * workeffvalues = effvalues.getV();
  for(j=0;j<posbeg.size();j++,workeffvalues++)
    *workeffvalues = d(index(posbeg[j],0),0);

  XX = datamatrix(posbeg.size());
  compute_XWX(likep->get_weight(),0);


  setbeta(posbeg.size(),1,0);

  identifiable = false;

  compind = statmatrix<unsigned>(nrpar,1,1);
  muy = datamatrix(nrpar,1);
  mu = datamatrix(index.rows(),1);

//  identifiable =true;

  }


FULLCOND_mixture::FULLCOND_mixture(const FULLCOND_mixture & fc)
                            : FULLCOND(FULLCOND(fc))
  {
  nrcomp = fc.nrcomp;
  compweight=fc.compweight;
  cwprior=fc.cwprior;
  csize=fc.csize;
  compind=fc.compind;
    compmean=fc.compmean;
    compvar=fc.compvar;
    cmpriorm=fc.cmpriorm;
    cmpriorv=fc.cmpriorv;
    cvpriorsh=fc.cvpriorsh;
    cvpriorsc=fc.cvpriorsc;

  muy = fc.muy;
    mu = fc.mu;
  fcconst = fc.fcconst;
  randomslope = fc.randomslope;
  includefixed = fc.includefixed;
  XX = fc.XX;
  likep = fc.likep;
  index = fc.index;
  index2 = fc.index2;
  posbeg = fc.posbeg;
  posend = fc.posend;
  effvalues = fc.effvalues;
  sigma2 = fc.sigma2;
  pathsample_total = fc.pathsample_total;
  ftotal = fc.ftotal;
  spatialtotal = fc.spatialtotal;
  lambda = fc.lambda;
  lambdaold1 = fc.lambdaold1;
  lambdaold2 = fc.lambdaold2;
  df_lambdaold1 = fc.df_lambdaold1;
  df_lambdaold2 = fc.df_lambdaold2;
  lambdaconst=fc.lambdaconst;
  data2 = fc.data2;
  }


const FULLCOND_mixture & FULLCOND_mixture::
         operator=(const FULLCOND_mixture & fc)
  {
  if (this==&fc)
    return *this;

  FULLCOND::operator=(FULLCOND(fc));

  nrcomp = fc.nrcomp;
  compweight=fc.compweight;
  compind=fc.compind;
    compmean=fc.compmean;
    compvar=fc.compvar;
    cmpriorm=fc.cmpriorm;
    cmpriorv=fc.cmpriorv;
    cvpriorsh=fc.cvpriorsh;
    cvpriorsc=fc.cvpriorsc;

  muy = fc.muy;
      mu = fc.mu;
  fcconst = fc.fcconst;
  randomslope = fc.randomslope;
  includefixed = fc.includefixed;
  XX = fc.XX;
  likep = fc.likep;
  index = fc.index;
  index2 = fc.index2;
  posbeg = fc.posbeg;
  posend = fc.posend;
  effvalues = fc.effvalues;
  sigma2 = fc.sigma2;
  pathsample_total = fc.pathsample_total;
  ftotal = fc.ftotal;
  spatialtotal = fc.spatialtotal;
  lambda = fc.lambda;
  lambdaold1 = fc.lambdaold1;
  lambdaold2 = fc.lambdaold2;
  df_lambdaold1 = fc.df_lambdaold1;
  df_lambdaold2 = fc.df_lambdaold2;
  lambdaconst=fc.lambdaconst;
  data2 = fc.data2;

  return *this;
  }



void FULLCOND_mixture::update(void)
  {
  unsigned i,k;

// Update component indicators
  for(i=0;i<compind.rows();i++)
  {
  datamatrix cprob(nrcomp,1,1.0/nrcomp); // probabilities psi_{ik}
  datamatrix cptemp(nrcomp,1,1.0/nrcomp);

  // calculate psi_{ik}
  //  double cptempsum=1.0;
  for(k=0;k<nrcomp;k++)
    {
    double cmean=compmean(k,0);
    double cvar=compvar(k,0);
    cptemp(k,0)=compweight(k,0) * (1.0/(sqrt(cvar))) * exp(-0.5*(beta(i,0)-cmean)*(1.0/cvar)*(beta(i,0)-cmean));
  //    cptemp(k,0)=compweight(k,0) * exp(-0.5*(beta(i,0)-cmean)*(beta(i,0)-cmean));
  //    cptempsum+=cptemp(k,0);
    }
  double cptempsum=cptemp.sum(0);
  cptemp=(1.0/cptempsum)*cptemp;
  cprob.assign(cptemp);

  // sample component indicator
  double u;
  u=uniform();

  double cprobsum=0.0;
  for(k=0;k<nrcomp;k++)
      {
      if ( (cprobsum<u) && (u<=cprobsum+cprob(k,0)))
        compind(i,0) = k+1;
      cprobsum+=cprob(k,0);
      }
/*
  // Sampling
  //  if(i==0)
  //  {
  ST::string pathtemp = pathcurrent.substr(0,pathcurrent.length()-4)+"_comptemp.res";
  ofstream outrest(pathtemp.strtochar(),ios::app);
  for(unsigned k=0;k<nrcomp;k++)
    {
    outrest << cprob(k,0) << "  " ;
    }
  outrest  << cptempsum << endl;
  //  }
*/
  }


// Update component sizes
statmatrix<unsigned> cstemp(nrcomp,1,0);

for(k=0;k<nrcomp;k++)
{
  for(i=0;i<nrpar;i++)
    {
    if(compind(i,0)==k+1) cstemp(k,0)+=1;
    }
}
csize.assign(cstemp);


// Update random effects
  if (optionsp->get_nriter()==1 || changingweight)
    compute_XWX(likep->get_weight(),0);
  update_linpred(false);

  likep->compute_respminuslinpred(mu,column);
  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();

  double * workbeta = beta.getV();
  int * workindex2 = index2.getV();
  double * workmuy = muy.getV();
  double * mup = mu.getV();
  likep->set_weightp();
  for(i=0;i<nrpar;i++,workmuy++,++itbeg,++itend)
      {
      *workmuy = 0;
      for(unsigned j=*itbeg;j<=*itend;j++,workindex2++)
        {
        mup += *workindex2;
        *workmuy+= likep->get_weight(*workindex2)* *mup;
        }
      }


  // update of random effects
  datamatrix retemp(beta.rows(),1,0);
  double scaletemp=likep->get_scale(column); // scale parameter sigma_i^2
  for(i=0;i<beta.rows();i++)
  {
  double remtemp,revtemp;
  double indobs=20; // number of observations for individual i=X_i'X_i
  double sumworkres=100; // sum of X_i'y_i
  unsigned comp=compind(i,0);
  revtemp = 1.0/( indobs*(1.0/scaletemp)+(1.0/compvar(comp-1,0)) );
  remtemp = revtemp*( sumworkres*(1.0/scaletemp)+(1.0/compvar(comp-1,0))*compmean(comp-1,0) );
  retemp(i,0)=remtemp+revtemp*rand_normal();
  }
  beta.assign(retemp);


  update_linpred(true);
  if (center)
    {
    double m = centerbeta();
    fcconst->update_intercept(m);
    }
  acceptance++;
  transform = likep->get_trmult(column);



// Update component means
  datamatrix cmeantemp(nrcomp,1,0);

  for(k=0;k<nrcomp;k++)
  {
  double mtemp,vtemp;
  vtemp = 1.0 / ( csize(k,0)*(1.0/compvar(k,0))+(1.0/cmpriorv) ) ;
  mtemp = vtemp*( csize(k,0)*(1.0/compvar(k,0))+(1.0/cmpriorv)*cmpriorm );

  cmeantemp(k,0)=mtemp+vtemp*rand_normal();
  }
  compmean.assign(cmeantemp);


// Update component variances
  datamatrix cvartemp(nrcomp,1,0);

  for(k=0;k<nrcomp;k++)
  {
  double shtemp,sctemp,resum;
  resum=0.0;
  for(i=0;i<nrpar;i++)
  {
    if(compind(i,0)==k+1)
      resum+=(beta(i,0)-compmean(k,0))*(beta(i,0)-compmean(k,0));
  }

  shtemp = cvpriorsh+csize(k,0);
  sctemp = cvpriorsc+resum;
  cvartemp(k,0)=rand_gamma(shtemp,sctemp);
  }
  compvar.assign(cvartemp);


// Update component weights
  datamatrix cwtemp(nrcomp,1,0);
  double cwtempsum;
  for(k=0;k<nrcomp;k++)
  {
  cwtemp(k,0)=rand_gamma(cwprior(k,0)+csize(k,0),1.0);
  }
  cwtempsum=cwtemp.sum(0);
  cwtemp=(1.0/cwtempsum)*cwtemp;
  compweight.assign(cwtemp);


  transform = likep->get_trmult(column);

  FULLCOND::update();

  }


void FULLCOND_mixture::outresults(void)
  {
  FULLCOND::outresults();

  ST::string vstr;

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);

  ST::string nl1 = l1;
  ST::string nl2 = l2;
  ST::string nu1 = u1;
  ST::string nu2 = u2;
  nl1 = nl1.replaceallsigns('.','p');
  nl2 = nl2.replaceallsigns('.','p');
  nu1 = nu1.replaceallsigns('.','p');
  nu2 = nu2.replaceallsigns('.','p');

  ST::string pathre = pathcurrent.substr(0,pathcurrent.length()-4)+"_re.res";
  ST::string pathcompind = pathcurrent.substr(0,pathcurrent.length()-4)+"_compind.res";

    optionsp->out("  Results for random effects are stored in file\n");
    optionsp->out("  " + pathre + "\n");
    optionsp->out("\n");
    optionsp->out("  Results for random effects mixture component indicators are stored in file\n");
    optionsp->out("  " + pathcompind + "\n");

/*
    if (lambdaconst==true)
    {
    optionsp->out("\n");
    optionsp->out("  Constant smoothing parameter: " +
    ST::doubletostring(lambda,6) + "\n");
    optionsp->out("\n");
    }
*/
    optionsp->out("\n");

    optionsp->out("\n");

//  Variable  mean           Std. Dev.      2.5% quant.    median         97.5% quant.
//  const     0.728983       0.0271208      0.67481        0.729393       0.782481

    optionsp->out("Results for means:\n");
    for(int k=0;k<nrcomp;k++)
       {
       optionsp->out("  Component " + ST::inttostring(k+1) + "  "
                     + ST::doubletostring(compmean(k,0),6) + "   ("
                     + ST::doubletostring(compmean(k,0),4) + ")"
                     );
       }
    optionsp->out("\n");

    optionsp->out("Results for variances:\n");
    for(int k=0;k<nrcomp;k++)
       {
       optionsp->out("  Component " + ST::inttostring(k+1) + "  "
                     + ST::doubletostring(compvar(k,0),6) + "   ("
                     + ST::doubletostring(compvar(k,0),4) + ")"
                     );
       }
    optionsp->out("\n");

    optionsp->out("Results for weights:\n");
    for(int k=0;k<nrcomp;k++)
       {
       optionsp->out("  Component " + ST::inttostring(k+1) + "  "
                     + ST::doubletostring(compweight(k,0),6) + "   ("
                     + ST::doubletostring(compweight(k,0),4) + ")"
                     );
       }
    optionsp->out("\n");

  unsigned i;
  ST::string name = datanames[0];

  // Ausgabe mixture random effects
  ofstream outres(pathre.strtochar());
  assert(!outres.fail());

  outres << "intnr" << "   ";
  outres << name << "   ";
  outres << "pmean   ";
  outres << "pqu"  << nl1  << "   ";
  outres << "pqu"  << nl2  << "   ";
  outres << "pmed   ";
  outres << "pqu"  << nu1  << "   ";
  outres << "pqu"  << nu2  << "   ";
  outres << "pcat" << level1 << "   ";
  outres << "pcat" << level2 << "   ";
  outres << endl;

     for(i=0;i<beta.rows();i++)
       {
       outres << (i+1) << "   ";
       outres << effvalues(i,0) << "   ";
       outres << betamean(i,0) << "   ";
       outres << betaqu_l1_lower(i,0) << "   ";
       outres << betaqu_l2_lower(i,0) << "   ";
       outres << betaqu50(i,0) << "   ";
       outres << betaqu_l2_upper(i,0) << "   ";
       outres << betaqu_l1_upper(i,0) << "   ";
       if (betaqu_l1_lower(i,0) > 0)
       outres << "1   ";
       else if (betaqu_l1_upper(i,0) < 0)
       outres << "-1   ";
       else
       outres << "0   ";
       if (betaqu_l2_lower(i,0) > 0)
       outres << "1   ";
       else if (betaqu_l2_upper(i,0) < 0)
         outres << "-1   ";
       else
         outres << "0   ";
      outres << endl;
      }

  // Ausgabe mixture component indicators
/*  ofstream outres2(pathcompind.strtochar());
  outres2 << "intnr" << "   ";
  outres2 << name << "   ";
  outres2 << "pmean   ";
  outres2 << "pmed   ";
  outres2 << endl;

     for(i=0;i<beta.rows();i++)
      {
      outres2 << (i+1) << "   ";
      outres2 << effvalues(i,0) << "   ";
      outres2 << betamean(i,0) << "   ";
      outres2 << betaqu50(i,0) << "   ";
      outres2 << endl;
      }
*/
  ofstream outres2(pathcompind.strtochar());
     for(i=0;i<compind.rows();i++)
      {
      outres2 << (i+1) << "   ";
      outres2 << compind(i,0);
      outres2 << endl;
      }

  ST::string pathcompvar = pathcurrent.substr(0,pathcurrent.length()-4)+"_compvar.res";
  ofstream outres3(pathcompvar.strtochar());
   for(unsigned i=0;i<compvar.rows();i++)
    {
    outres3 << (i+1) << "   ";
    outres3 << compvar(i,0);
    outres3 << endl;
    }

  }


void FULLCOND_mixture::outoptions(void)
  {

  optionsp->out("  OPTIONS FOR RANDOM EFFECT MIXTURE: " + title + "\n",true);
  optionsp->out("\n");

  optionsp->out("  Type of Mixture: Normal\n",true);
  optionsp->out("  Number of components: " + ST::inttostring(nrcomp) + "\n",true);

  optionsp->out("  Hyperparameter for component means:\n",true);
  optionsp->out("    Prior Means: " + ST::doubletostring(cmpriorm,2) + "\n",true);
  optionsp->out("    Prior Variances: " + ST::doubletostring(cmpriorv,2) + "\n",true);

  optionsp->out("  Hyperparameter for component variances:\n",true);
  optionsp->out("    Prior shape: " + ST::doubletostring(cvpriorsh,2) + "\n",true);
  optionsp->out("    Prior scale: " + ST::doubletostring(cvpriorsc,2) + "\n",true);

//  optionsp->out("  Prior variance of component means: " + ST::inttostring(cmpriorv) + "\n",true);
//  optionsp->out("  Prior shape of component variances: " + ST::inttostring(cvpriorsh) + "\n",true);
//  optionsp->out("  Prior scale of component variances: " + ST::inttostring(cvpriorsc) + "\n",true);
  optionsp->out("\n");
  }




void FULLCOND_mixture::update_linpred(const bool & add)
  {
  unsigned i,j;

  unsigned n = nrpar;

  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();
  double * workbeta = beta.getV();

  int * workindex2;

  if (add==false)
    {
//      likep->set_linpredp_current(column);
      for (i=0;i<nrpar;i++,++itbeg,++itend,workbeta++)
        if (*itbeg != -1)
          likep->add_linearpred2(-(*workbeta),*itbeg,*itend,index,index2,column);
    } // end: if (add == false)
  else
    {
//      likep->set_linpredp_current(column);
      for (i=0;i<nrpar;i++,++itbeg,++itend,workbeta++)
        if (*itbeg != -1)
          likep->add_linearpred2(*workbeta,*itbeg,*itend,index,index2,column);
    }
  }




bool FULLCOND_mixture::posteriormode(void)
  {

  unsigned n = nrpar;
  if (includefixed)
    n = nrpar-1;

  update_linpred(false);

  compute_XWX(likep->get_weightiwls(),column);

  likep->compute_weightiwls_workingresiduals(column);

  unsigned i,j;
  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();
  int * workindex2 = index2.getV();
  itbeg = posbeg.begin();
  itend = posend.begin();
  double * workmuy = muy.getV();
  likep->set_workingresp();

  for(i=0;i<nrpar;i++,workmuy++,++itbeg,++itend)
      {
      *workmuy = 0;
      for(j=*itbeg;j<=*itend;j++,workindex2++)
        {
        *workmuy+= likep->get_workingres(*workindex2);
        }
      }

  itbeg = posbeg.begin();
  itend = posend.begin();
  workmuy = muy.getV();
  double * workbeta = beta.getV();
  double * workXX = XX.getV();

  for(i=0;i<n;i++,workmuy++,++itbeg,++itend,workbeta++,workXX++)
    {
    *workbeta = (*workmuy)/(*workXX+lambda);
    }

  update_linpred(true);

  transform = likep->get_trmult(column);

  return FULLCOND::posteriormode();

  }


//------------------------------------------------------------------------------
//------------------- class FULLCOND_mixture_gaussian ---------------------------
//------------------------------------------------------------------------------



void FULLCOND_mixture_gaussian::update(void)
  {
/*
  double var;
  double m;
  unsigned i,j;
  unsigned n = nrpar;


  if (optionsp->get_nriter()==1 || changingweight)
    compute_XWX(likep->get_weight(),0);

  if (lambdaconst == false)
    lambda = likep->get_scale(column)/sigma2;
  else
    sigma2 = likep->get_scale(column)/lambda;

  double sqrtscale = sqrt(likep->get_scale(column));

  update_linpred(false);

  // nicht verändern wegen SUR-Modellen
  likep->compute_respminuslinpred(mu,column);

  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();
  double * workbeta = beta.getV();

  int * workindex2 = index2.getV();
  double * workmuy = muy.getV();
  double * mup = mu.getV();
  likep->set_weightp();

  for(i=0;i<nrpar;i++,workmuy++,++itbeg,++itend)
      {
      *workmuy = 0;
      for(j=*itbeg;j<=*itend;j++,workindex2++)
        {
        mup += *workindex2;
        *workmuy+= likep->get_weight(*workindex2)* *mup;
        }
      }

  workbeta = beta.getV();
  workmuy = muy.getV();
  double * workXX = XX.getV();
  for (i=0;i<n;i++,workbeta++,workmuy++,workXX++)
    {
    var = 1.0/(*workXX  + lambda);
    m = var * *workmuy;
    *workbeta = m + sqrtscale*sqrt(var)*rand_normal();
    }

  update_linpred(true);


  if (center)
    {
    double m = centerbeta();
    fcconst->update_intercept(m);
    }

  acceptance++;

  transform = likep->get_trmult(column);
*/
  FULLCOND_mixture::update();
  }


} // end: namespace MCMC




