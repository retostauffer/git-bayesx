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

    priorassumptions.push_back("Gaussian mixture random effect");
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

    priorassumptions.push_back("Gaussian mixture random effect");
    }



FULLCOND_mixture::FULLCOND_mixture(MCMCoptions * o,DISTRIBUTION * dp,
                              FULLCOND_const * fcc,
                              const datamatrix & d, const ST::string & t,
                              const ST::string & fp,const ST::string & pr,
                              const int & nrc,const double & pw,
                              const double & pmm,const double & pmv,
                              const double & pva,const double & pvb,
                              const bool & s,
                              const unsigned & c)
                            : FULLCOND(o,datamatrix(1,1),t,1,1,fp)
  {
  nrcomp = nrc;
  compweight = datamatrix(nrcomp,1,1.0/nrcomp);
  compmean = datamatrix(nrcomp,1,0);
  compvar = datamatrix(nrcomp,1,1);

  cwprior = datamatrix(nrcomp,1,pw);
  cmpriorm = pmm;
  cmpriorv = pmv;
  cvpriora = pva;
  cvpriorb = pvb;
  nosamples = s;

  csize = statmatrix<unsigned>(nrcomp,1,1);
  temp = datamatrix(nrcomp,1,0);
  checkorder = false;

  fcconst = fcc;
  fctype = randomeffects;

  column = c;

  likep = dp;
  pathresult = pr;
  pathcurrent = pr;

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

  setbeta(posbeg.size(),1,0);

  compind = statmatrix<unsigned>(beta.rows(),1,1);

  ST::string path = samplepath.substr(0,samplepath.length()-4)+"_compparsample.raw";
  cpar_fc = FULLCOND(o,datamatrix(1,1),t+"_cpar_fc",nrcomp,3,path);
  cpar_fc.setflags(MCMC::norelchange | MCMC::nooutput);

  ST::string path2 = samplepath.substr(0,samplepath.length()-4)+"_compindsample.raw";
  cind_fc = FULLCOND(o,datamatrix(1,1),t+"_cind_fc",beta.rows(),1,path2);
  cind_fc.setflags(MCMC::norelchange | MCMC::nooutput);

  transform = likep->get_trmult(column);

  identifiable = false;
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
  cvpriora=fc.cvpriora;
  cvpriorb=fc.cvpriorb;
  nosamples=fc.nosamples;
  temp=fc.temp;
  checkorder = fc.checkorder;

  fcconst = fc.fcconst;
  cpar_fc = fc.cpar_fc;
  cind_fc = fc.cind_fc;

  likep = fc.likep;
  index = fc.index;
  index2 = fc.index2;
  posbeg = fc.posbeg;
  posend = fc.posend;
  effvalues = fc.effvalues;
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
  cvpriora=fc.cvpriora;
  cvpriorb=fc.cvpriorb;
  nosamples=fc.nosamples;
  temp=fc.temp;
  checkorder = fc.checkorder;

  fcconst = fc.fcconst;
  cpar_fc = fc.cpar_fc;
  cind_fc = fc.cind_fc;

  likep = fc.likep;
  index = fc.index;
  index2 = fc.index2;
  posbeg = fc.posbeg;
  posend = fc.posend;
  effvalues = fc.effvalues;

  return *this;
  }


double FULLCOND_mixture::centerbeta(void)
  {
  unsigned i;
  double sum=0;
  sum = beta.sum(0);
  sum/=beta.rows();
  for (i=0;i<beta.rows();i++)
    beta(i,0) -= sum;
  return sum;
 }


void FULLCOND_mixture::update_weights(void)
  {
  double cwtempsum;
  unsigned k;
  for(k=0;k<nrcomp;k++)
    {
    temp(k,0)=rand_gamma(cwprior(k,0)+csize(k,0),1.0);
    }
  cwtempsum=temp.sum(0);
  temp=(1.0/cwtempsum)*temp;
  compweight.assign(temp);

  checkorder=true;
  for(k=0;k<nrcomp-1;k++)
    {
    if(compweight(k,0)>compweight(k+1,0)) checkorder=false;
    }
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
  for(k=0;k<nrcomp;k++)
    {
    double cmean=compmean(k,0);
    double cvar=compvar(k,0);
    cptemp(k,0)=compweight(k,0) * (1.0/(sqrt(cvar))) * exp(-0.5*(beta(i,0)-cmean)*(1.0/cvar)*(beta(i,0)-cmean));
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
  }


// Update component sizes
  statmatrix<unsigned> cstemp(nrcomp,1,0);
  for(k=0;k<nrcomp;k++)
  {
    for(i=0;i<beta.rows();i++)
      {
      if(compind(i,0)==k+1) cstemp(k,0)+=1;
      }
  }
  csize.assign(cstemp);



// Update random effects
  double scaletemp=likep->get_scale(column); // scale parameter sigma_i^2
  double remtemp,revtemp,indobs,sumworkres;
  unsigned comp,j;

  for(i=0;i<beta.rows();i++)
  {
  likep->add_linearpred(-1.0*beta(i,0),posbeg[i],posend[i],index,0,true);

  indobs=posend[i]-posbeg[i]+1; // number of observations for individual i=X_i'X_i
  sumworkres=0;
  for(j=posbeg[i];j<=posend[i];j++)
    {
    sumworkres += likep->get_response(j,0)-likep->get_linearpred(j,0);  // sum of X_i'(y_i-eta_i)
    }

  comp=compind(i,0);
  revtemp = 1.0/( indobs*(1.0/scaletemp)+(1.0/compvar(comp-1,0)) );
  remtemp = revtemp*( sumworkres*(1.0/scaletemp)+(1.0/compvar(comp-1,0))*compmean(comp-1,0) );
  beta(i,0)=remtemp+sqrt(revtemp)*rand_normal();

  likep->add_linearpred(beta(i,0),posbeg[i],posend[i],index,0,true);
  }



// Update component means
  double mtemp,vtemp,remean,resum;

  for(k=0;k<nrcomp;k++)
  {
  resum=0.0;
  for(i=0;i<beta.rows();i++)
    {
    if(compind(i,0)==k+1)
      resum+=beta(i,0);
    }
  if(csize(k,0)==0)
    remean=0.0;
  else
    remean=resum/csize(k,0);

  vtemp = 1.0 / ( csize(k,0)*(1.0/compvar(k,0))+(1.0/cmpriorv) ) ;
  mtemp = vtemp*( csize(k,0)*(1.0/compvar(k,0))*remean+(1.0/cmpriorv)*cmpriorm );

  compmean(k,0)=mtemp+sqrt(vtemp)*rand_normal();
  }



// Update component variances
  for(k=0;k<nrcomp;k++)
  {
  double shtemp,sctemp,resum;
  resum=0.0;
  for(i=0;i<beta.rows();i++)
  {
    if(compind(i,0)==k+1)
      resum+=(beta(i,0)-compmean(k,0))*(beta(i,0)-compmean(k,0));
  }

  shtemp = cvpriora+0.5*csize(k,0);
  sctemp = cvpriorb+0.5*resum;
  compvar(k,0)=rand_invgamma(shtemp,sctemp);
  }



// Update component weights
checkorder=false;
while(checkorder==false)
  update_weights();
/*
  double cwtempsum;
  for(unsigned k=0;k<nrcomp;k++)
    {
    temp(k,0)=rand_gamma(cwprior(k,0)+csize(k,0),1.0);
    }
  cwtempsum=temp.sum(0);
  temp=(1.0/cwtempsum)*temp;
  compweight.assign(temp);
*/

  double m = centerbeta();
  fcconst->update_intercept(m);


  double * cp_p = cpar_fc.getbetapointer();
  for(k=0;k<nrcomp;k++)
  {
  *cp_p = compweight(k,0);
  cp_p++;
  *cp_p = compmean(k,0)*transform;
  cp_p++;
  *cp_p = compvar(k,0)*transform*transform;
  cp_p++;
  }
  cpar_fc.update();


  double * ci_p = cind_fc.getbetapointer();
  for(i=0;i<beta.rows();i++)
  {
  *ci_p = compind(i,0);
  ci_p++;
  }
  cind_fc.update();

  acceptance++;

  FULLCOND::update();

/*
 // Sampling
  ST::string pt1 = pathcurrent.substr(0,pathcurrent.length()-4)+"_cmeantemp.res";
  ofstream ot1(pt1.strtochar(),ios::app);
  for(unsigned k=0;k<nrcomp;k++)
    {
    ot1 << compmean(k,0)*likep->get_trmult(column) << "  " ;
    }
  ot1 << endl;

  ST::string pt2 = pathcurrent.substr(0,pathcurrent.length()-4)+"_cvartemp.res";
  ofstream ot2(pt2.strtochar(),ios::app);
  for(unsigned k=0;k<nrcomp;k++)
    {
    ot2 << compvar(k,0)*likep->get_trmult(column) << "  " ;
    }
  ot2 << endl;
*/
 }


void FULLCOND_mixture::outresults(void)
  {
  FULLCOND::outresults();

  unsigned i,k;

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
  ST::string pathcpar = pathcurrent.substr(0,pathcurrent.length()-4)+"_comppar.res";
  ST::string name = datanames[0];

  // Ausgabe Schätzungen random effects
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

// Ausgabe Schätzungen mixture component parameters
  cpar_fc.outresults();
  ofstream outres2(pathcpar.strtochar());
  outres2 << "paramnr   ";
  outres2 << "varname   ";
  outres2 << "pmean    ";
  outres2 << "pstddev    ";
  outres2 << "pqu"  << nl1  << "   ";
  outres2 << "pqu"  << nl2  << "   ";
  outres2 << "pmed   ";
  outres2 << "pqu"  << nu1  << "   ";
  outres2 << "pqu"  << nu2  << "   ";
  outres2 << endl;

  for(k=0;k<nrcomp;k++)
    {
    outres2 << (k+1) << "   ";
    outres2 << "weight" << (k+1) << "   ";
    outres2 << cpar_fc.get_betamean(k,0) << "   ";
    if(nrcomp==1)
      outres2 << 0 << "   ";
    else
      outres2 << sqrt(cpar_fc.get_betavar(k,0)) << "   ";
    outres2 << cpar_fc.get_beta_lower1(k,0) << "   ";
    outres2 << cpar_fc.get_beta_lower2(k,0) << "   ";
    outres2 << cpar_fc.get_betaqu50(k,0) << "   ";
    outres2 << cpar_fc.get_beta_upper2(k,0) << "   ";
    outres2 << cpar_fc.get_beta_upper1(k,0) << "   ";
    outres2 << endl;
    }
  for(k=0;k<nrcomp;k++)
    {
    outres2 << (k+nrcomp+1) << "   ";
    outres2 << "mean" << (k+1) << "   ";
    outres2 << cpar_fc.get_betamean(k,1) << "   ";
    outres2 << sqrt(cpar_fc.get_betavar(k,1)) << "   ";
    outres2 << cpar_fc.get_beta_lower1(k,1) << "   ";
    outres2 << cpar_fc.get_beta_lower2(k,1) << "   ";
    outres2 << cpar_fc.get_betaqu50(k,1) << "   ";
    outres2 << cpar_fc.get_beta_upper2(k,1) << "   ";
    outres2 << cpar_fc.get_beta_upper1(k,1) << "   ";
    outres2 << endl;
    }
  for(k=0;k<nrcomp;k++)
    {
    outres2 << (k+2*nrcomp+1) << "   ";    
    outres2 << "var" << (k+1) << "   ";
    outres2 << cpar_fc.get_betamean(k,2) << "   ";
    outres2 << sqrt(cpar_fc.get_betavar(k,2)) << "   ";
    outres2 << cpar_fc.get_beta_lower1(k,2) << "   ";
    outres2 << cpar_fc.get_beta_lower2(k,2) << "   ";
    outres2 << cpar_fc.get_betaqu50(k,2) << "   ";
    outres2 << cpar_fc.get_beta_upper2(k,2) << "   ";
    outres2 << cpar_fc.get_beta_upper1(k,2) << "   ";
    outres2 << endl;
    }


  ST::string help =  ST::doubletostring(lower1,4) + "% quant.";
  ST::string levell = help + ST::string(' ',15-help.length());
  help = ST::doubletostring(upper2,4) + "% quant.";
  ST::string levelu = help + ST::string(' ',15-help.length());

  optionsp->out("  Results for component weights:\n");
  optionsp->out("   Component   Post. Mean     Std. Dev.      " +
                  levell + "Post. Median   " + levelu + "\n");
  for(k=0;k<nrcomp;k++)
     {
     if(nrcomp==1)
       {
       optionsp->out("     " + ST::outresults(7,ST::inttostring(k+1),cpar_fc.get_betamean(k,0),
                      0,cpar_fc.get_beta_lower1(k,0),
                      cpar_fc.get_betaqu50(k,0),cpar_fc.get_beta_upper1(k,0)) + "\n");
       }
     else
       {
       optionsp->out("     " + ST::outresults(7,ST::inttostring(k+1),cpar_fc.get_betamean(k,0),
                      sqrt(cpar_fc.get_betavar(k,0)),cpar_fc.get_beta_lower1(k,0),
                      cpar_fc.get_betaqu50(k,0),cpar_fc.get_beta_upper1(k,0)) + "\n");
       }
     }
  optionsp->out("\n");

  optionsp->out("  Results for component means:\n");
  optionsp->out("   Component   Post. Mean     Std. Dev.      " +
                  levell + "Post. Median   " + levelu + "\n");
  for(k=0;k<nrcomp;k++)
     {
     optionsp->out("     " + ST::outresults(7,ST::inttostring(k+1),cpar_fc.get_betamean(k,1),
                      sqrt(cpar_fc.get_betavar(k,1)),cpar_fc.get_beta_lower1(k,1),
                      cpar_fc.get_betaqu50(k,1),cpar_fc.get_beta_upper1(k,1)) + "\n");
     }
  optionsp->out("\n");

  optionsp->out("  Results for component variances:\n");
  optionsp->out("   Component   Post. Mean     Std. Dev.      " +
                  levell + "Post. Median   " + levelu + "\n");
  for(k=0;k<nrcomp;k++)
     {
     optionsp->out("     " + ST::outresults(7,ST::inttostring(k+1),cpar_fc.get_betamean(k,2),
                      sqrt(cpar_fc.get_betavar(k,2)),cpar_fc.get_beta_lower1(k,2),
                      cpar_fc.get_betaqu50(k,2),cpar_fc.get_beta_upper1(k,2)) + "\n");
     }
  optionsp->out("\n\n");

  optionsp->out("  Results for estimated mixture component parameters are also stored in file\n");
  optionsp->out("  " + pathcpar + "\n");
  optionsp->out("\n");


  if (nosamples==false)
    {
    ST::string file = pathcurrent.substr(0,pathcurrent.length()-4) + "_comppar_sample.raw";
    cpar_fc.get_samples(file);
    optionsp->out("  Sampling paths for mixture component parameters are stored in file\n");
    optionsp->out("  " + file + "\n");
    optionsp->out("\n");

    cind_fc.outresults();
    ST::string file2 = pathcurrent.substr(0,pathcurrent.length()-4) + "_compind_sample.raw";
    cind_fc.get_samples(file2);
    optionsp->out("  Sampling paths for mixture component indicators are stored in file\n");
    optionsp->out("  " + file2 + "\n");
    optionsp->out("\n");
    }

  optionsp->out("  Results for estimated random effects are stored in file\n");
  optionsp->out("  " + pathre + "\n");
  optionsp->out("\n");
 }



void FULLCOND_mixture::outoptions(void)
  {

  optionsp->out("  OPTIONS FOR MIXTURE RANDOM EFFECT: " + title + "\n",true);
  optionsp->out("\n");

  optionsp->out("  Number of components: " + ST::inttostring(nrcomp) + "\n",true);
  optionsp->out("  Type of Mixture: Normal\n",true);
  optionsp->out("  Prior parameter for component weights: " + ST::doubletostring(cwprior(0,0)) + "\n",true);
  optionsp->out("  Prior parameters for component means:\n",true);
  optionsp->out("    Hyperprior for means: " + ST::doubletostring(cmpriorm,4) + "\n",true);
  optionsp->out("    Hyperprior for variances: " + ST::doubletostring(cmpriorv,4) + "\n",true);
  optionsp->out("  Prior parameters for component variances:\n",true);
  optionsp->out("    Hyperprior for a: " + ST::doubletostring(cvpriora,4) + "\n",true);
  optionsp->out("    Hyperprior for b: " + ST::doubletostring(cvpriorb,4) + "\n",true);
  optionsp->out("\n");
  }




bool FULLCOND_mixture::posteriormode(void)
  {
  return true;
  }




} // end: namespace MCMC




