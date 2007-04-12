
#include "first.h"

#include "randomeffect_stepwise.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//----------------- class FULLCOND_random_stepwise -----------------------------
//------------------------------------------------------------------------------

FULLCOND_random_stepwise::FULLCOND_random_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                              FULLCOND_const * fcc,
                              const datamatrix & d, const ST::string & t,
                              const ST::string & fp,const ST::string & pr,
                              const double & la, const unsigned & c)
                            : FULLCOND_random(o,dp,fcc,d,t,fp,pr,la,c)
  {
  identifiable = false;
  grenzfall = 0;
  intercept = 0.0;

  dimX = 0;
  dimZ = nrpar;

  //gleichwertig = true;
  }

// randomslope
FULLCOND_random_stepwise::FULLCOND_random_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                  FULLCOND_const * fcc,
                  const datamatrix & intvar,const datamatrix & effmod,
                  const ST::string & t,
                  const ST::string & fp,const ST::string & pr,
                  const ST::string & prf,
                  const double & la,
                  const bool & inclfixed,const unsigned & c)
                  : FULLCOND_random(o,dp,fcc,intvar,effmod,t,fp,pr,prf,la,false,c)
  {
  grenzfall = 0;
  dimX = 0;
  dimZ = nrpar;

  if(inclfixed == false)
    {
    identifiable = false;
    }

  includefixed = false;
  intercept = 0.0;
  //gleichwertig = true;
  }


FULLCOND_random_stepwise::FULLCOND_random_stepwise(const FULLCOND_random_stepwise & fc)
                            : FULLCOND_random(FULLCOND_random(fc))
  {
  intercept = fc.intercept;
  beta_average = fc.beta_average;
  data_varcoeff_fix = fc.data_varcoeff_fix;
  effmodi = fc.effmodi;
  fbasisp = fc.fbasisp;
  df_unstruct = fc.df_unstruct;
  fc_df = fc.fc_df;
  //gleichwertig = fc.gleichwertig;
  }


const FULLCOND_random_stepwise & FULLCOND_random_stepwise::
         operator=(const FULLCOND_random_stepwise & fc)
  {
  if (this==&fc)
    return *this;

  FULLCOND_random::operator=(FULLCOND_random(fc));

  intercept = fc.intercept;
  beta_average = fc.beta_average;
  data_varcoeff_fix = fc.data_varcoeff_fix;
  effmodi = fc.effmodi;
  fbasisp = fc.fbasisp;
  df_unstruct = fc.df_unstruct;
  fc_df = fc.fc_df;  
  //gleichwertig = fc.gleichwertig;

  return *this;
  }


bool FULLCOND_random_stepwise::posteriormode(void)
  {
  unsigned n = nrpar;

  update_linpred(false);

  if(calculate_xwx == true)
    {
    calculate_xwx = false;
    compute_XWX(likep->get_weightiwls(),column);
    }

  likep->compute_weightiwls_workingresiduals(column);

  unsigned i,j;
  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();
  int * workindex2 = index2.getV();
  itbeg = posbeg.begin();
  itend = posend.begin();
  double * workmuy = muy.getV();
  likep->set_workingresp();

  if (!randomslope)
    {
    for(i=0;i<nrpar;i++,workmuy++,++itbeg,++itend)
      {
      *workmuy = 0;
      for(j=*itbeg;j<=*itend;j++,workindex2++)
        {
        *workmuy+= likep->get_workingres(*workindex2);
        }
      }
    }
  else
    {
    double * datap = data.getV();
    for(i=0;i<n;i++,workmuy++,++itbeg,++itend)
      {
      *workmuy = 0;
      for(j=*itbeg;j<=*itend;j++,workindex2++,datap++)
        {
        *workmuy+= likep->get_workingres(*workindex2)* (*datap);
        }
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

  if (randomslope && center)
    {
    double * workbeta = beta.getV();
    double sum=0;
    for (i=0;i<n;i++,workbeta++)
      {
      sum += *workbeta;
      }

    intercept = sum/double(n);

    workbeta = beta.getV();
    for (i=0;i<n;i++,workbeta++)
      *workbeta -= intercept;

    update_fix_effect(intercept);
    intercept = 0.0;
    }

  update_linpred(true);

  transform = likep->get_trmult(column);

  return FULLCOND::posteriormode();

  }


// BEGIN: For Varying Coefficients ---------------------------------------------

void FULLCOND_random_stepwise::update_fix_effect(double & intercept)
  {
  bool raus = false;
  unsigned j = 1;
  ST::string name_richtig = datanames[0];
  while(j<fcconst->get_datanames().size() && raus==false)
     {
     if(fcconst->get_datanames()[j] == datanames[0])
        {
        raus = true;
        }
     if(fcconst->get_datanames()[j] == (datanames[0]+"_1"))
        {
        raus = true;
        name_richtig = datanames[0] + "_1";
        }
     j = j + 1;
     }
  if(raus == true)
    {
    fcconst->update_fix_effect(j-1,intercept,data_forfixed);
    }
  else
    {
    vector<ST::string> names;
    names.push_back(name_richtig);
    fcconst->include_effect(names,data_forfixed);
    interactions_pointer[0]->set_inthemodel(-1);
    fcconst->update_fix_effect(j,intercept,data_forfixed);
    }
  }


void FULLCOND_random_stepwise::set_pointer_to_interaction(FULLCOND * inter)
  {
  interactions_pointer.push_back(inter);
  }

void FULLCOND_random_stepwise::get_interactionspointer(vector<FULLCOND*> & inter)
  {
  inter = interactions_pointer;
  }


void FULLCOND_random_stepwise::hierarchical(ST::string & possible)
  {
  if(!randomslope)
    {
    possible = "alles";
    }
  else
    {
    possible = "valles";
    }
  }


void FULLCOND_random_stepwise::const_varcoeff(void)
  {
  if(randomslope)
    fcconst->posteriormode_const_varcoeff(data_forfixed);
  }

// END: For Varying Coefficients -----------------------------------------------


void FULLCOND_random_stepwise::create_weight(datamatrix & w)
  {
  if(!spatialtotal)
    {
    vector<unsigned>::iterator itbeg = posbeg.begin();
    vector<unsigned>::iterator itend = posend.begin();
    int * workindex = index.getV();
    unsigned i;

    unsigned j;
    for(i=0;i<nrpar;i++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        {
        w(*workindex,0) = 1;
        for(j=*itbeg;j<=*itend;j++)
          workindex++;
        }
      }
    }
  }


void FULLCOND_random_stepwise::compute_lambdavec(vector<double> & lvec, int & number)
  {
  if (df_equidist==true && spfromdf==true)
     FULLCOND::compute_lambdavec_equi(lvec,number);
  else
     FULLCOND::compute_lambdavec(lvec,number);
  if (!nofixed && randomslope && identifiable)
    hierarchie_fix(lvec,1);
  if(forced_into==false)
     lvec.push_back(0);

  // Startwert für lambda aus df:
  if(spfromdf==true)
    {
    double lambdavorg = 1000;
    if(!randomslope)
      {
      if(dfstart==0)
        lambdastart = 0;
      else
        lambdastart = lambda_from_df(dfstart,lambdavorg);
      }
    else
      {
      if(!nofixed && dfstart==1 && identifiable)
        lambdastart = -1;
      else if(dfstart==0)
        lambdastart = 0;
      else
        lambdastart = lambda_from_df(dfstart,lambdavorg);
      }
    if(lambdastart==-9 || lambdastart==1000000000)    // falls dfstart nicht erreicht werden kann
      lambdastart = 0;
    }

  }


void FULLCOND_random_stepwise::hierarchie_fix(vector<double> & untervector, int dfo)
  {

  unsigned number = untervector.size()-1;

  update_stepwise(untervector[0]);
  double df_max = compute_df();

  update_stepwise(untervector[number]);
  double df_min = compute_df();

  if(df_max > dfo && df_min < dfo)
     {
     bool geordnet = false;
     unsigned stelle_oben = number;
     unsigned stelle_unten = 0;
     while(geordnet==false)
        {
        unsigned stelle = stelle_oben + stelle_unten;
        update_stepwise(untervector[stelle/2]);
        double df_mitteunten = compute_df();
        update_stepwise(untervector[stelle/2 + 1]);
        double df_mitteoben = compute_df();

        if(df_mitteunten > dfo && df_mitteoben > dfo)
          stelle_unten = stelle/2;
        else if(df_mitteunten < dfo && df_mitteoben < dfo)
          stelle_oben = stelle/2 + 1;
        else
          {
          geordnet = true;
          vector<double> hilf;
          unsigned i;
          stelle_unten = stelle/2;
          stelle_oben = stelle/2 + 1;
          for(i=0;i<=stelle_unten;i++)
             hilf.push_back(untervector[i]);
          hilf.push_back(-1);
          for(i=stelle_oben;i<untervector.size();i++)
            hilf.push_back(untervector[i]);
          untervector = hilf;
          }
        }
     }
  else if(df_min >= dfo)
     untervector.push_back(-1);
  else
     {
     vector<double> hilf;
     hilf.push_back(-1);
     unsigned i;
     for(i=0;i<untervector.size();i++)
        hilf.push_back(untervector[i]);
     untervector = hilf;
     }
  }


const datamatrix & FULLCOND_random_stepwise::get_data_forfixedeffects(void)
  {
  // useful for randomslopes only
  if ( (data_forfixed.rows() < data.rows()) && (randomslope==true) )
    {
    data_forfixed=datamatrix(data.rows(),1);
    unsigned i;
    int * workindex = index.getV();
    double * workdata = data.getV();
    for (i=0;i<data.rows();i++,workindex++,workdata++)
      {
      data_forfixed(*workindex,0) = *workdata;
      }
    }

  return data_forfixed;
  }


double FULLCOND_random_stepwise::compute_df(void)
  {
  double df = 0;
  if(inthemodel == true)
    {
    unsigned i;
    bool struct_included = false;
    if(spatialtotal)
      {
      bool fix;
      fbasisp->get_inthemodel(struct_included,fix);
      if(struct_included == true)
        {
        df = df_unstruct;
        }
      }

    if(!spatialtotal || !struct_included)
      {
      if (lambdaold1==lambda && likep->get_iwlsweights_notchanged() == true && !spatialtotal)
        {
        df = df_lambdaold1;
        }
      else
        {
        if(calculate_xwx == true)
          {
          calculate_xwx = false;
          compute_XWX(likep->get_weightiwls(),column);
          }
        double * workXX=XX.getV();
        unsigned n = nrpar;

        if(!identifiable)
          {
          double c = 0;
          double w = 0;
          for(i=0;i<n;i++,workXX++)
            {
            c += (*workXX * *workXX)/(*workXX+lambda);
            w += *workXX;
            }
         c = 1/(w - c);

         workXX = XX.getV();
         for(i=0;i<n;i++,workXX++)
           {
           df += (*workXX * (lambda + *workXX * (-c * (*workXX + 2*lambda) + 1)))/((*workXX+lambda) * (*workXX+lambda));
           }
         df += w*c - 1;
         }
        else
          {
          for(i=0;i<n;i++,workXX++)
            df += (*workXX)/(*workXX+lambda);
          }

        df_lambdaold1 = df;
        lambdaold1 = lambda;
        }
      }
    }

  return df;
  }


void FULLCOND_random_stepwise::set_dfunstruct(const double & df_unstr)
  {
  df_unstruct = df_unstr;
  }


void FULLCOND_random_stepwise::update_stepwise(double la)
  {
  lambda = la;
  }

double FULLCOND_random_stepwise::get_lambda(void)
  {
  return lambda;
  }


ST::string FULLCOND_random_stepwise::get_effect(void)
  {

  ST::string h;

  if(randomslope)
    h = datanames[0] + "*" + datanames[1];
  else
    h = datanames[0];

  h = h + "(random,df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";

  return h;

  }


void FULLCOND_random_stepwise::init_names(const vector<ST::string> & na)
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
      ST::string helpname1 = na[1].insert_string_char(charh,stringh);
      ST::string helpname2 = na[0].insert_string_char(charh,stringh);
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


void FULLCOND_random_stepwise::reset_effect(const unsigned & pos)
  {

  update_linpred(false);

  unsigned i;
  double * work;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;

  intercept = 0.0;
  }


// BEGIN: MODEL-AVERAGING ------------------------------------------------------

void FULLCOND_random_stepwise::update_bootstrap(const bool & uncond)
  {
  if(optionsp->get_samplesize()==1)
    {
    ST::string path = samplepath.substr(0,samplepath.length()-4)+"_df.raw";
    fc_df = FULLCOND(optionsp,datamatrix(1,1),"title?",1,1,path);
    fc_df.setflags(MCMC::norelchange | MCMC::nooutput);
    }

  datamatrix betaold = beta;

  if(fixornot==true)
    {
    bool raus = false;
    unsigned j = 1;
    ST::string name_richtig = datanames[0];
    while(j<fcconst->get_datanames().size() && raus==false)
      {
      if(fcconst->get_datanames()[j] == datanames[0])
        raus = true;
      j = j + 1;
      }
    unsigned index_fix = j-1;
    double fix = fcconst->getbeta(index_fix,0);
    unsigned i;
    double * workbeta = beta.getV();
    vector<unsigned>::iterator itbeg = posbeg.begin();
    vector<unsigned>::iterator itend = posend.begin();
    //int * workindex = index2.getV();
    //int k;
    for(i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        {
        *workbeta = fix;
        //for(k=*itbeg;k<=*itend;k++)
        //  workindex++;
        }
      }
    workbeta = beta.getV();

    FULLCOND::update_bootstrap();
    fc_df.setbetavalue(0,0,-1);
    fc_df.update_bootstrap();
    }
  else if(inthemodel==false && fixornot==false)
    {
    beta = datamatrix(nrpar,1,0);
    FULLCOND::update_bootstrap();
    fc_df.setbetavalue(0,0,0);
    fc_df.update_bootstrap();
    }
  else
    {
    FULLCOND::update_bootstrap();
    fc_df.setbetavalue(0,0,lambda);
    fc_df.update_bootstrap();
    }
  beta = betaold;
  }


void FULLCOND_random_stepwise::update_bootstrap_betamean(void)
  {
  FULLCOND::update_bootstrap_betamean();
  FULLCOND::setflags(MCMC::norelchange);

  fc_df.update_bootstrap_betamean();
  //fc_df.outresults();
  double * workmean = fc_df.get_betameanp();

  ST::string pathdf = pathcurrent.substr(0,pathcurrent.length()-4)+"_df.res";
  ofstream outres(pathdf.strtochar());

  outres << "value   ";
  outres << "frequency  ";
  outres << "pmean   " << endl;

// Häufigkeitstabelle:

  //samplestream.close();
  datamatrix sample(optionsp->get_samplesize(),1);
  fc_df.readsample(sample,0);
  unsigned i;

  vector<unsigned> number;
  vector<unsigned> number1;
  vector<unsigned> number2;
  vector<unsigned> cumnumber1;
  vector<unsigned> cumnumber;

  statmatrix<int> index(sample.rows(),1);
  index.indexinit();
  sample.indexsort(index,0,sample.rows()-1,0,0);

  i = 0;
  unsigned j,anz;
  while(i<index.rows())
     {
     anz=0;
     int* p = index.getV() + i;
     int* q = index.getV() + i;
     j=i;
     while(j<index.rows() && (sample.get(*p,0) == sample.get(*q,0)))
        {
        anz = anz+1;
        j++;
        p++;  
        }
     if(sample.get(*q,0) <= 0)
       number1.push_back(anz);
     else if(sample.get(*q,0) > 0)
       number2.push_back(anz);
     if(cumnumber1.size()>0)
       cumnumber1.push_back(cumnumber1[cumnumber1.size()-1]+anz);
     else
       cumnumber1.push_back(anz);
     i = i + anz;
     }

  int k;
  for(k=number1.size()-1;k>=0;k--)
    {
    cumnumber.push_back(cumnumber1[k]);
    number.push_back(number1[k]);
    }
  for(k=number2.size()-1;k>=0;k--)
    {
    cumnumber.push_back(cumnumber1[k+number1.size()]);
    number.push_back(number2[k]);
    }

  for(i=0;i<number.size();i++)
    {
    double help = sample.get(index(cumnumber[i]-1,0),0);
    double dfs = -1*help;
    if(help>0)
      {
      update_stepwise(help);
      set_inthemodel(help);
      dfs = compute_df();
      }
    outres << ST::doubletostring(dfs,6) << "   " << ST::inttostring(number[i]) << "   ";
    if(*workmean == help)
      outres << "selected"; // ST::doubletostring(*workmean,6);
    else
      outres << "-";
    outres << endl;
    }
  }


void FULLCOND_random_stepwise::init_spatialtotal(FULLCOND_nonp_basis * sp,
                                        const ST::string & pnt,
                                        const ST::string & prt)
  {
  df_unstruct = 0;
  fbasisp = sp;
  vector<ST::string> ev = sp->get_effectvalues();

  FULLCOND_random::init_spatialtotal(ev,pnt,prt);

  }


ST::string FULLCOND_random_stepwise::get_befehl(void)
  {
  ST::string h;

  if(randomslope)
    h = datanames[0] + "*" + datanames[1];
  else
    h = datanames[0];

  h = h + "(random,lambda=" + ST::doubletostring(lambda,6) + ")";
  return h;
  }


} // end: namespace MCMC




