
#include "mcmc_const.h"

namespace MCMC
{

ST::string FULLCOND_const::get_effect(void)
  {
  ST::string h="";
  unsigned i;
  if (fctype==MCMC::factor)
    {
    for(i=0;i<datanames.size();i++)
      h = h + " + " + datanames[i];
    }
  else
    {
    h = datanames[0];
    for(i=1;i<datanames.size();i++)
      h = h + " + " + datanames[i];
    }
  return h;
  }

double FULLCOND_const::compute_df(void)
  {
  if (lambda==-1)
    return data.cols();
  else
    return 0;
  }


void FULLCOND_const::update_stepwise(double la)
  {
  lambda=la;
  }


void FULLCOND_const::transfer_interceptsample(void)
  {

  if (identifiable==false)
    {
    datamatrix s(optionsp->compute_samplesize(),1);
    readsample(s,0);

    likep->set_interceptsample(s,column);
    }

  }


bool FULLCOND_const::is_missing(const ST::string & na)
  {
  unsigned i=0;
  bool found=false;
  while ( (i < datanames.size()) && (found==false) )
    {
    if (datanames[i] == na)
      found=true;
    i++;
    }

  return found;
  }


void FULLCOND_const::init_name(const ST::string & na)
    {
    if (fctype == MCMC::factor)
      {
      vector<ST::string> nam;
      unsigned i;
      for (i=0;i<diff_categories.size();i++)
        {
        if (diff_categories[i] != reference)
          nam.push_back(na+"_" + ST::doubletostring(diff_categories[i]));
        }
      FULLCOND::init_names(nam);

      char charh ='_';
      ST::string stringh = "\\_";

      ST::string helpname;

      for(i=0;i<nam.size();i++)
        {
        helpname = nam[i].insert_string_char(charh,stringh);
        term_symbolic = term_symbolic + "\\gamma_{"+helpname+"}"+helpname;
        if (i+1<nam.size())
          term_symbolic = term_symbolic + " + ";
        }
      int c = column;
      if(c==0)
         {
         priorassumptions.push_back("Factor $" + na.insert_string_char(charh,stringh) + "$:");
         ST::string liste = "";
         for(i=0;i<nam.size()-1;i++)
            liste = liste + "$" + nam[i].insert_string_char(charh,stringh) + "$, ";
         liste = liste +  + "$" + nam[nam.size()-1].insert_string_char(charh,stringh) + "$";
         priorassumptions.push_back("Resulting variables: " + liste);
         priorassumptions.push_back("diffuse priors");
         priorassumptions.push_back("Coding: " + coding);
         priorassumptions.push_back("\\\\");
         }
      if(c>0)
         {
         priorassumptions.push_back(
         "Factor " + na + " (" + ST::inttostring(c+1) + ". response category):");
         priorassumptions.push_back("diffuse priors");
         priorassumptions.push_back("\\\\");
         }

      }
    else
      {
      FULLCOND::init_name(na);
      char charh ='_';
      ST::string stringh = "\\_";
      ST::string helpname = na.insert_string_char(charh,stringh);
      term_symbolic = "\\gamma_{"+helpname+"}"+helpname;
      int c = column;
      if(c==0)
        {
        priorassumptions.push_back("Fixed effects:");
        priorassumptions.push_back("diffuse priors");
        priorassumptions.push_back("\\\\");
        }
      if(c>0)
        {
        priorassumptions.push_back("Fixed effects (" + ST::inttostring(c+1) + ". response category):");
        priorassumptions.push_back("diffuse priors");
        priorassumptions.push_back("\\\\");
        }
      }

    }


void FULLCOND_const::init_names(const vector<ST::string> & na)
    {
    FULLCOND::init_names(na);
    unsigned i;

    char charh ='_';
    ST::string stringh = "\\_";

    ST::string helpname;
    term_symbolic = "";    //"\\gamma_{"+helpname+"}"+helpname + " + ";
    for(i=0;i<na.size();i++)
      {
      helpname = na[i].insert_string_char(charh,stringh);
      term_symbolic = term_symbolic + "\\gamma_{"+helpname+"}"+helpname;
      if (i+1<na.size())
        term_symbolic = term_symbolic + " + ";
      }
    int c = column;
    if(c==0)
       {
       priorassumptions.push_back("Fixed effects:");
       priorassumptions.push_back("diffuse priors");
       priorassumptions.push_back("\\\\");
       }
    if(c>0)
       {
       priorassumptions.push_back("Fixed effects (" + ST::inttostring(c+1) + ". response category):");
       priorassumptions.push_back("diffuse priors");
       priorassumptions.push_back("\\\\");
       }
    }


  // CONSTRUCTOR_3 (REML)

FULLCOND_const::FULLCOND_const(MCMCoptions * op,const datamatrix & d,
                               const ST::string & t,
                               const int & constant,const ST::string & fs,
                               const ST::string & fr)
                  : FULLCOND(op,t)
  {

  data = d;

  nrconst = data.cols();

  dimX = nrconst;
  dimZ = 0;

  plotstyle=MCMC::noplot;

  pathresult = fr;

  pathcurrent = fr;

  results_type="fixed";

  }


void FULLCOND_const::createreml(datamatrix & X,datamatrix & Z,
                                const unsigned & Xpos, const unsigned & Zpos)
  {
  unsigned i,j;

  double * workdata= data.getV();
  double * workX = X.getV()+Xpos;
  unsigned s = X.cols()-nrconst;
  for (i=0;i<data.rows();i++,workX+=s)
    for(j=0;j<nrconst;j++,workdata++,workX++)
      {
      *workX = *workdata;
      }

  }


double FULLCOND_const::outresultsreml(datamatrix & X,datamatrix & Z,
                                    datamatrix & betareml, datamatrix & betacov,
                                    datamatrix & thetareml,
                                    const unsigned & Xpos,const unsigned & Zpos,
                                    const unsigned & thetapos,
                                    const bool & dispers,
                                    const unsigned & betaXpos,
                                    const unsigned & betaZpos,
                                    const double & category,
                                    const bool & ismultinomial,
                                    const unsigned plotpos)
  {
  double meanhelp=0;
  betamean=datamatrix(nrconst,1,0);
  datamatrix betastd=datamatrix(nrconst,1,0);
  betaqu_l1_lower=datamatrix(nrconst,1,0);
  betaqu_l1_upper=datamatrix(nrconst,1,0);
  betaqu_l2_lower=datamatrix(nrconst,1,0);
  betaqu_l2_upper=datamatrix(nrconst,1,0);
  datamatrix betapval=datamatrix(nrconst,1,0);

  unsigned i;
  for(i=0;i<nrconst;i++)
    {
//    betamean(i,0) = betareml(i,0);
//    betastd(i,0) = sqrt(betacov(i,i));
    betamean(i,0) = betareml(betaXpos+i,0);
    betastd(i,0) = sqrt(betacov(betaXpos+i,betaXpos+i));
    betaqu_l1_lower(i,0) = betamean(i,0)+randnumbers::invPhi2(lower1/100)*betastd(i,0);
    betaqu_l1_upper(i,0) = betamean(i,0)+randnumbers::invPhi2(upper2/100)*betastd(i,0);
    betaqu_l2_lower(i,0) = betamean(i,0)+randnumbers::invPhi2(lower2/100)*betastd(i,0);
    betaqu_l2_upper(i,0) = betamean(i,0)+randnumbers::invPhi2(upper1/100)*betastd(i,0);
    if (betamean(i,0)>0)
      {
      betapval(i,0) = 2*(1-randnumbers::Phi2(betamean(i,0)/betastd(i,0)));
      }
    else
      {
      betapval(i,0) = 2*(randnumbers::Phi2(betamean(i,0)/betastd(i,0)));
      }
    }

  FULLCOND::outresults();

  ST::string outest=pathcurrent;
  if(ismultinomial)
    {
    outest = pathcurrent.insert_after_string("_"+ST::doubletostring(category,6),"FixedEffects");
    }
  ofstream outp(outest.strtochar());

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  vector<ST::string> resultstable(6);

  outp << "paramnr varname pmode ci" << level1 << "lower ci" << level2 <<
             "lower std ci" << level2 << "upper ci" << level1 << "upper pcat"
             << level1 << " pcat" << level2 << " pvalue" << endl;

  ST::string l;
  int maxvarnamelength = 0;
  int len;

  for(i=0;i<nrconst;i++)
    {
    len = datanames[i].length();
    if (len > maxvarnamelength)
      maxvarnamelength = len;
    }

  if (maxvarnamelength>10)
    l = ST::string(' ',maxvarnamelength-6);
  else
    l = "  ";

  ST::string help =  ST::doubletostring(level1,4) + "% Confidence Interval";
  ST::string ci = help + ST::string(' ',30-help.length());

  optionsp->out("  Variable" + l +
                "Post. Mode     " +
                "Std. Dev.      " +
                "p-value        " +
                ci + "\n");

  ST::string mean;
  ST::string std;
  ST::string qu10;
  ST::string qu50;
  ST::string qu90;

  unsigned nsp;

  for (i=0;i<nrconst;i++)
    {

    if (maxvarnamelength  > 10)
      nsp = 2+maxvarnamelength-datanames[i].length();
    else
      nsp = 10-datanames[i].length();

    outp << (i+1) << "   ";
    outp << datanames[i] << "   ";
    outp << betamean(i,0) << "   ";
    outp << betaqu_l1_lower(i,0) << "   ";
    outp << betaqu_l2_lower(i,0) << "   ";
    outp << betastd(i,0) << "   ";
    outp << betaqu_l2_upper(i,0) << "   ";
    outp << betaqu_l1_upper(i,0) << "   ";
    if (betaqu_l1_lower(i,0) > 0)
      outp << "1   ";
    else if (betaqu_l1_upper(i,0) < 0)
      outp << "-1   ";
    else
      outp << "0   ";

    if (betaqu_l2_lower(i,0) > 0)
      outp << "1   ";
    else if (betaqu_l2_upper(i,0) < 0)
      outp << "-1   ";
    else
      outp << "0   ";
    outp << betapval(i,0) << "   ";
    outp << endl;

    optionsp->out(ST::outresults(nsp,datanames[i],betamean(i,0),
                    betastd(i,0),betapval(i,0),betaqu_l1_lower(i,0),
                    betaqu_l1_upper(i,0)) + "\n");

    char hchar = '_';
    ST::string hstring = "\\_";
    resultstable[0] = datanames[i].insert_string_char(hchar,hstring);

    resultstable[1] = ST::doubletostring(betamean(i,0),6);
    resultstable[2] = ST::doubletostring(betastd(i,0),6);
    resultstable[3] = ST::doubletostring(betapval(i,0),6);
    resultstable[4] = ST::doubletostring(betaqu_l1_lower(i,0),6);
    resultstable[5] = ST::doubletostring(betaqu_l1_upper(i,0),6);

    results_latex.push_back(ST::make_latextable(resultstable));

    }

  optionsp->out("\n");

  optionsp->out("  Results for fixed effects are also stored in file\n");
  optionsp->out("  " + outest + "\n");

  optionsp->out("\n");
  return meanhelp;
  }

void FULLCOND_const::outresultsreml_ordinal(datamatrix & X,datamatrix & Z,
                     datamatrix & betareml, datamatrix & betacov,
                     unsigned nrcat2)
  {
  if(term_symbolic.length()>20)
    {
    term_symbolic = "\\theta "+term_symbolic.substr(20,term_symbolic.length()-20);
    }
  else
    {
    term_symbolic="\\theta ";
    }
  nrconst=dimX+nrcat2-1;
  datanames[0]="theta_"+ST::inttostring(nrcat2);
  int j;
  for(j=nrcat2-1; j>0; j--)
    {
    datanames.insert(datanames.begin(),"theta_"+ST::inttostring(j));
    }
  outresultsreml(X,Z,betareml,betacov,datamatrix(1,1,0),0,0,0,false,0,0,0,false,0);
  }


void FULLCOND_const::outresults(void)
  {
  FULLCOND::outresults();

  ofstream outp(pathcurrent.strtochar());

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  vector<ST::string> resultstable(6);

  outp << "paramnr varname pmean pstd pqu" << l1 << " pqu" << l2 <<
            " pmed pqu" << u1 << " pqu" << u2 << " pcat" << level1
            << " pcat" << level2 << endl;

  unsigned i;

  optionsp->out("\n");

  ST::string l;
  int maxvarnamelength = 0;
  int len;

  for(i=0;i<nrconst;i++)
    {
    len = datanames[i].length();
    if (len > maxvarnamelength)
      maxvarnamelength = len;
    }

  if (maxvarnamelength>10)
    l = ST::string(' ',maxvarnamelength-6);
  else
    l = "  ";

    ST::string help =  ST::doubletostring(lower1,4) + "% quant.";
    ST::string levell = help + ST::string(' ',15-help.length());
    help = ST::doubletostring(upper2,4) + "% quant.";
    ST::string levelu = help + ST::string(' ',15-help.length());

    optionsp->out("  Variable" + l +
                  "mean           " +
                  "Std. Dev.      " +
                  levell +
                  "median         " +
                  levelu + "\n");

    ST::string mean;
    ST::string std;
    ST::string qu10;
    ST::string qu50;
    ST::string qu90;

    double m,stddouble;

    unsigned nsp;

    for (i=0;i<nrconst;i++)
      {

      if (maxvarnamelength  > 10)
        nsp = 2+maxvarnamelength-datanames[i].length();
      else
        nsp = 10-datanames[i].length();

      if ((optionsp->get_samplesize() == 0) &&
          (interceptyes) &&
          (i==interceptpos)
          )
        {
        m = betamean(i,0)+likep->get_addinterceptsample();
        }
      else
        {
        m= betamean(i,0);
        }

      if (betavar(i,0) == 0)
        stddouble = 0;
      else
        stddouble = sqrt(betavar(i,0));

      outp << (i+1) << "   ";
      outp << datanames[i] << "   ";
      outp << m << "   ";
      outp << stddouble << "   ";
      outp << betaqu_l1_lower(i,0) << "   ";
      outp << betaqu_l2_lower(i,0) << "   ";
      outp << betaqu50(i,0) << "   ";
      outp << betaqu_l2_upper(i,0) << "   ";
      outp << betaqu_l1_upper(i,0) << "   ";
      if (betaqu_l1_lower(i,0) > 0)
        outp << "1   ";
      else if (betaqu_l1_upper(i,0) < 0)
        outp << "-1   ";
      else
        outp << "0   ";

      if (betaqu_l2_lower(i,0) > 0)
        outp << "1   ";
      else if (betaqu_l2_upper(i,0) < 0)
        outp << "-1   ";
      else
        outp << "0   ";

      outp << endl;

      optionsp->out(ST::outresults(nsp,datanames[i],m,
                      stddouble,betaqu_l1_lower(i,0),
                      betaqu50(i,0),betaqu_l1_upper(i,0)) + "\n");


      char hchar = '_';
      ST::string hstring = "\\_";
      resultstable[0] = datanames[i].insert_string_char(hchar,hstring);
      resultstable[1] = ST::doubletostring(m,6);
      resultstable[2] = ST::doubletostring(stddouble,6);
      resultstable[3] = ST::doubletostring(betaqu_l1_lower(i,0),6);
      resultstable[4] = ST::doubletostring(betaqu50(i,0),6);
      resultstable[5] = ST::doubletostring(betaqu_l1_upper(i,0),6);

      results_latex.push_back(ST::make_latextable(resultstable));

      }

    optionsp->out("\n");

    optionsp->out("  Results for fixed effects are also stored in file\n");
    optionsp->out("  " + pathcurrent + "\n");

    optionsp->out("\n");

  }


void FULLCOND_const::make_design(datamatrix & d)   //, const ST::string & code)
  {
  //coding = code;
  vector<unsigned> zaehlen;

  statmatrix<int> index(d.rows(),1);
  index.indexinit();
  d.indexsort(index,0,d.rows()-1,0,0);

  unsigned i = 0;
  unsigned j;
  while(i<index.rows())
     {
     double anz=0;
     int* p = index.getV() + i;
     int* q = index.getV() + i;
     for(j=i;j<index.rows();j++,p++)
        {
        if (d.get(*p,0) == d.get(*q,0))
                  anz = anz+1;
               }
               zaehlen.push_back(anz);
               diff_categories.push_back(d.get(*q,0));
               i = i + anz;
          }

  bool gefunden = false;
  for(i=0;i<diff_categories.size();i++)
     {
     if(diff_categories[i]==reference)
        gefunden = true;
     }
  if(gefunden==false)
     {
     optionsp->out("WARNING: The value for the reference category doesn't exist!\n");
     optionsp->out("Category " + ST::doubletostring(diff_categories[0]) + " used instead\n");
     reference = diff_categories[0];
     }

  data = datamatrix(d.rows(),diff_categories.size()-1);
  unsigned spalte = 0;
  unsigned zeile = 0;
  unsigned ref;
  for(i=0;i<diff_categories.size();i++)
     {
     if(diff_categories[i]!=reference)
       {
       int* q = index.getV();
       for(j=0;j<zeile;j++,q++)
          data(*q,spalte) = 0;
       q = index.getV() + zeile;
       for(j=zeile;j<zeile+zaehlen[i];j++,q++)
          data(*q,spalte) = 1;
       q = index.getV() + zeile + zaehlen[i];
       for(j=zeile+zaehlen[i];j<d.rows();j++,q++)
          data(*q,spalte) = 0;
       spalte = spalte + 1;
       }

     else if(diff_categories[i]==reference && coding=="effect")
       ref = i;

     zeile = zeile + zaehlen[i];
     }

  if(coding=="effect")
     {
     zeile = 0;
     for(i=0;i<ref;i++)
        zeile = zeile + zaehlen[i];
     int* q = index.getV() + zeile;
     for(i=zeile;i<zeile+zaehlen[ref];i++,q++)
        {
        for(j=0;j<data.cols();j++)
           data(*q,j) = -1;
        }
     }
  }


void FULLCOND_const::compute_lambdavec(vector<double> & lvec,unsigned & number)
  {

  assert(fctype == MCMC::factor);
  lvec.push_back(-1);
  get_forced();
  if(forced_into==false)
     lvec.push_back(0);

  }

// For stepwise factor
FULLCOND_const::FULLCOND_const(MCMCoptions * o,DISTRIBUTION * dp,
                 const datamatrix & d,const ST::string & code, int
                 & ref,
                 const ST::string & t,
                 const ST::string & fs,const ST::string & fr,
                 const unsigned & c)
  {

  lambda=-1;

  fctype = MCMC::factor;

  reference = ref;

  coding = code;

  make_design(d);

  interceptadd=0;

  likep = dp;

  sumold = 0;

  datamatrix w = likep->get_weight();

  column = c;

  nrconst = data.cols();

  interceptyes = false;

  pathresult = fr;
  pathcurrent = fr;

  results_type="fixed";

//  negbin=false;


  }



FULLCOND_const::FULLCOND_const (MCMCoptions * o,DISTRIBUTION * dp,
                                const datamatrix & d, const ST::string & t,
                                const int & constant, const ST::string & fs,
                                const ST::string & fr,
                                const unsigned & c)
          : FULLCOND(o,d,t,d.cols(),1,fs)
  {

  lambda=-1;

  interceptadd=0;

  fctype = MCMC::fixed;

  likep = dp;

  sumold = 0;

  datamatrix w = likep->get_weight();

  column = c;

  nrconst = data.cols();

  linold = datamatrix(likep->get_nrobs(),1,0);
  linnew = linold;
  linoldp = &linold;
  linnewp = &linnew;

  if (constant>=0)
    {
    identifiable = false;
    interceptpos = constant;
    interceptyes = true;
    }
  else
    {
    interceptyes = false;
    interceptpos = constant;
    }

  pathresult = fr;
  pathcurrent = fr;

  results_type="fixed";

//  negbin=false;

  }


FULLCOND_const::FULLCOND_const(const FULLCOND_const & m) : FULLCOND(FULLCOND(m))
  {
  lambda=m.lambda;
  reference = m.reference;
  coding = m.coding;
  diff_categories = m.diff_categories;
  interceptadd = m.interceptadd;
//  negbin=m.negbin;
  likep = m.likep;
  nrconst = m.nrconst;
  linold = m.linold;
  linnew = m.linnew;
  sumold = m.sumold;
  pathresult = m.pathresult;
  table_results = m.table_results;
  interceptpos = m.interceptpos;
  interceptyes = m.interceptyes;
  }


const FULLCOND_const & FULLCOND_const::operator=(const FULLCOND_const & m)
  {
  if (this==&m)
	 return *this;
  FULLCOND::operator=(FULLCOND(m));
  lambda=m.lambda;
  reference = m.reference;
  diff_categories = m.diff_categories;
  interceptadd = m.interceptadd;
//  negbin=m.negbin;
  likep = m.likep;
  nrconst = m.nrconst;
  linold = m.linold;
  linnew = m.linnew;
  sumold = m.sumold;
  pathresult = m.pathresult;
  table_results = m.table_results;
  interceptpos = m.interceptpos;
  interceptyes = m.interceptyes;
  return *this;
  }


void FULLCOND_const::update(void)
  {

  if (interceptyes)
    {
    beta(interceptpos,0) += likep->get_addinterceptsample();
    }

  FULLCOND::update();

  if (interceptyes)
    {
    beta(interceptpos,0) -= likep->get_addinterceptsample();
    }


  if  (optionsp->get_nriter() == optionsp->get_iterations())
    {
    FULLCOND::outresults();
    transfer_interceptsample();
    }

  }


bool FULLCOND_const::posteriormode(void)
  {

  return FULLCOND::posteriormode();

  }



void FULLCOND_const::outoptions(void)
  {

  optionsp->out("  OPTIONS FOR FIXED EFFECTS:\n",true);
  optionsp->out("\n");
  optionsp->out("  Priors:\n");
  optionsp->out("\n");
  optionsp->out("  diffuse priors\n");
  optionsp->out("\n");

  }


//------------------------------------------------------------------------------
//------------------ CLASS: FULLCOND_const_gaussian_special --------------------
//------------------------------------------------------------------------------

double FULLCOND_const_gaussian_special::compute_df(void)
  {
  if (lambda==0)
    return 0;
  else
    return 1.0;
  }


FULLCOND_const_gaussian_special::FULLCOND_const_gaussian_special(
                                  MCMCoptions * o,DISTRIBUTION * dp,
                                  const datamatrix & d,const ST::string & t,
                                  const ST::string & fs,const ST::string & fr,
                                  const unsigned & c)
  : FULLCOND_const(o,dp,d,t,false,fs,fr,c)
  {

  fctype = MCMC::nonlinearf;
  lambda = -1;
  transform = likep->get_trmult(c);
  datatransformed=data;
  mu = datamatrix(data.rows(),1);
  }


FULLCOND_const_gaussian_special::FULLCOND_const_gaussian_special(
const FULLCOND_const_gaussian_special & m)  : FULLCOND_const(FULLCOND_const(m))
  {
  datatransformed = m.datatransformed;
  mu = m.mu;
  }



const FULLCOND_const_gaussian_special &
FULLCOND_const_gaussian_special::operator=(
const FULLCOND_const_gaussian_special & m)
  {
  if (this==&m)
	 return *this;
  FULLCOND_const::operator=(FULLCOND_const(m));
  datatransformed = m.datatransformed;
  mu = m.mu;
  return *this;
  }


void FULLCOND_const_gaussian_special::update_stepwise(double la)
  {
  lambda=la;
  unsigned i=0;
  reset_effect(i);
  compute_datatransformed(lambda);
  }



bool FULLCOND_const_gaussian_special::posteriormode(void)
  {
  /*
  double xwx =0;
  unsigned i;
  double * workdatatransformed = datatransformed.getV();
  likep->set_weightp(0);
  double * workweight = likep->get_weightp();
  double * worklinold=linold.getV();
  for(i=0;i<data.rows();i++,workdatatransformed++,workweight++,worklinold++)
    {
    *worklinold = beta(0,0) * *workdatatransformed;
    xwx += *workdatatransformed * *workdatatransformed * *workweight;

    }

  likep->substr_linearpred_m(linold,column);

  mu.minus(likep->get_response(),likep->get_linearpred());

  double sumy = 0;
  likep->set_weightp(0);
  workweight = likep->get_weightp();
  double * workmu = mu.getV();
  workdatatransformed = datatransformed.getV();
  for(i=0;i<data.rows();i++,workdatatransformed++,workweight++,workmu++)
    {
    sumy += *workweight * *workdatatransformed * *workmu;
    }

  beta(0,0) = sumy/xwx;

  double * worklinnew=linnew.getV();
  workdatatransformed = datatransformed.getV();
  for(i=0;i<data.rows();i++,worklinnew++,workdatatransformed++)
    *worklinnew = beta(0,0) * *workdatatransformed;

  likep->add_linearpred_m(linnew,column);
  */
  return FULLCOND_const::posteriormode();
  }


bool FULLCOND_const_gaussian_special::posteriormode_converged(
const unsigned & itnr)
  {
  return likep->posteriormode_converged_fc(beta,beta_mode,itnr);
  }



void FULLCOND_const_gaussian_special::outresults(void)
  {
  FULLCOND_const::outresults();
  }


void FULLCOND_const_gaussian_special::compute_datatransformed(double lambda)
  {
  double * workdata=data.getV();
  double * workdatatransformed = datatransformed.getV();
  unsigned i;
  for (i=0;i<data.rows();i++,workdata++,workdatatransformed++)
    {
    if (lambda == -1)
      *workdatatransformed = *workdata;
    else if (lambda == 1)
      *workdatatransformed = log(*workdata);
    else if ( lambda == 2)
      *workdatatransformed = 1.0/(*workdata+1);
    }

  }


void FULLCOND_const_gaussian_special::compute_lambdavec(
  vector<double> & lvec,unsigned & number)
  {

  lvec.push_back(2);       // 1/x+1
  lvec.push_back(1);       // ln
  lvec.push_back(-1);
  get_forced();
  if(forced_into==false)
     lvec.push_back(0);

  }


ST::string  FULLCOND_const_gaussian_special::get_effect(void)
  {
  ST::string h;
  if (lambda==-1)
    h = datanames[0];
  else if (lambda==1)
    h = "log(" + datanames[0] + ")";
  else if (lambda == 2)
    h = "1/(" + datanames[0] + "+1)";

  return h;
  }


void FULLCOND_const_gaussian_special::reset_effect(unsigned & pos)
  {
  double * worklinnew=linnew.getV();
  double * workdatatransformed = datatransformed.getV();
  unsigned i;
  for(i=0;i<data.rows();i++,worklinnew++,workdatatransformed++)
    *worklinnew = -beta(0,0) * *workdatatransformed;

  likep->add_linearpred_m(linnew,column);

  beta(0,0) = 0;
  }

//------------------------------------------------------------------------------
//------------------ CLASS: FULLCOND_const_gaussian ----------------------------
//------------------------------------------------------------------------------

void FULLCOND_const_gaussian::update_missings(datamatrix & x,
                                              const datamatrix & linpredx,
                                              const statmatrix<unsigned> & w,
                                              const ST::string & na,
                                              double & scalex)
  {

  unsigned i=0;
  bool found=false;
  unsigned nr;
  while ( (i < datanames.size()) && (found==false) )
    {
    if (datanames[i] == na)
      nr = i;
    i++;
    }

  double mux;
  double sigma2x;

  unsigned * workw = w.getV();

  double * workdata = data.getV() + (nr + data.cols()*(*workw));

  double dataold;
  double * workx = x.getV()+(*workw);
  double * worklinpredx = linpredx.getV()+(*workw);
  unsigned respos = *workw;
  double r;
  for(i=0;i<w.rows();i++,workw++,workx+=(*workw),worklinpredx+=(*workw),
      workdata+=data.cols()*(*workw),respos+=(*workw))
    {

    sigma2x = (beta(nr,0)*beta(nr,0))/likep->get_scale(column) +
               1.0/scalex;
    sigma2x = 1.0/sigma2x;

    r = likep->get_response(respos,0)-likep->get_linearpred(respos,0)
        +(*workdata)*beta(nr,0);

    mux = sigma2x*(beta(nr,0)*r/likep->get_scale(column)
          + (*worklinpredx)/scalex);
    *workx = mux + sqrt(sigma2x)*rand_normal();
    dataold = *workdata;
    *workdata = *workx;

    likep->add_linearpred(beta(nr,0)*(*workdata-dataold),respos,column);
    linold(respos,0) += beta(nr,0)*(*workdata-dataold);

    }

  compute_matrices();
  }

/*
void FULLCOND_const_gaussian::include_effect(vector<ST::string> & names,
datamatrix & newx)
  {

  ST::string name = names[0];
  unsigned i,j;

  nrconst++;
  datamatrix dataold = data;

  data = datamatrix(data.rows(),nrconst);

  double * workold = dataold.getV();
  double * workdata = data.getV();
  double * worknew = newx.getV();

  for(i=0;i<dataold.rows();i++,worknew++,workdata++)
    {

    for (j=0;j<dataold.cols();j++,workold++,workdata++)
      {
      *workdata = *workold;
      }

    *workdata = *worknew;

    }

  datanames.push_back(name);

  X1 = datamatrix(nrconst,nrconst,0);

  X2 = datamatrix(nrconst,likep->get_nrobs());


  datamatrix betao = beta;

  setbeta(nrconst,1,0);

  double * workbeta = beta.getV();
  double * workbetao = betao.getV();
  double * workbetameanold = betameanold.getV();
  for(i=0;i<betao.rows();i++,workbetao++,workbeta++,workbetameanold++)
    {
    *workbeta=*workbetao;
    *workbetameanold = *workbeta;
    }

  }
*/


void FULLCOND_const_gaussian::include_effect(vector<ST::string> & names,
datamatrix & newx)
  {
  unsigned i,j;

  nrconst+=names.size();
  datamatrix dataold = data;

  data = datamatrix(data.rows(),nrconst);

  double * workold = dataold.getV();
  double * workdata = data.getV();
  double * worknew = newx.getV();

  for(i=0;i<dataold.rows();i++)
    {

    for (j=0;j<dataold.cols();j++,workold++,workdata++)
      {
      *workdata = *workold;
      }

    for (j=0;j<newx.cols();j++,workdata++,worknew++)
      *workdata = *worknew;

    }


  for (j=0;j<names.size();j++)
    {
    datanames.push_back(names[j]);
    }

  X1 = datamatrix(nrconst,nrconst,0);

  X2 = datamatrix(nrconst,likep->get_nrobs());

  datamatrix betao = beta;

  setbeta(nrconst,1,0);

  double * workbeta = beta.getV();
  double * workbetao = betao.getV();
  double * workbetameanold = betameanold.getV();
  for(i=0;i<betao.rows();i++,workbetao++,workbeta++,workbetameanold++)
    {
    *workbeta=*workbetao;
    *workbetameanold = *workbeta;
    }

  }


void FULLCOND_const_gaussian::reset_effect(unsigned & pos)
  {

  unsigned i,j;

  nrconst--;
  datamatrix dataold = data;

  data = datamatrix(data.rows(),nrconst);

  double * workold = dataold.getV();
  double * workdata = data.getV();

  for(i=0;i<dataold.rows();i++)
    {

    for (j=0;j<dataold.cols();j++,workold++)
      {
      if (j!=pos)
        {
        *workdata = *workold;
        workdata++;
        }
      }

    }

  vector<ST::string> dn = datanames;

  datanames.erase(datanames.begin(),datanames.end());
//  datanames.reserve(dn.size()-1);
  for(i=0;i<dn.size();i++)
    {
    if (i!=pos)
      datanames.push_back(dn[i]);
    }

  X1 = datamatrix(nrconst,nrconst,0);

  X2 = datamatrix(nrconst,likep->get_nrobs());


  datamatrix betao = beta;

  setbeta(nrconst,1,0);

  double * workbeta = beta.getV();
  double * workbetao = betao.getV();
  double * workbetameanold = betameanold.getV();
  for(i=0;i<betao.rows();i++,workbetao++)
    {
    if (i!=pos)
      {
      *workbeta=*workbetao;
      *workbetameanold = *workbeta;
      workbeta++;
      workbetameanold++;
      }
    }


  likep->substr_linearpred_m(linold,column);
  linold.mult(data,beta);
  likep->add_linearpred_m(linold,column);


  }


FULLCOND_const_gaussian::FULLCOND_const_gaussian(MCMCoptions * o,
                         DISTRIBUTION* dp, const datamatrix & d,
                         const ST::string & t, const int & constant,
                         const ST::string & fs,const ST::string & fr,
                         const unsigned & c)
                         : FULLCOND_const(o,dp,d,t,constant,fs,fr,c)
  {

  transform = likep->get_trmult(c);

  changingweight = likep->get_changingweight();

  mu1 = datamatrix(likep->get_nrobs(),1);

  X1 = datamatrix(nrconst,nrconst,0);

  help = datamatrix(nrconst,likep->get_nrobs(),0);

  X2 = datamatrix(nrconst,likep->get_nrobs());

  compute_matrices();

  if (X1.rows() < nrconst)
    errors.push_back("ERROR: design matrix for fixed effects is rank deficient\n");

  }


FULLCOND_const_gaussian::FULLCOND_const_gaussian(
  const FULLCOND_const_gaussian & m)
  : FULLCOND_const(FULLCOND_const(m))
  {
  X1 = m.X1;
  X2 = m.X2;
  mu1 = m.mu1;
  help = m.help;
  changingweight = m.changingweight;
  }


const FULLCOND_const_gaussian & FULLCOND_const_gaussian::operator=
(const FULLCOND_const_gaussian & m)
  {
  if (this==&m)
	 return *this;
  FULLCOND_const::operator=(FULLCOND_const(m));
  changingweight = m.changingweight;
  X1 = m.X1;
  X2 = m.X2;
  mu1 = m.mu1;
  help = m.help;
  return *this;
  }


void FULLCOND_const_gaussian::compute_matrices(void)
  {

  // computing X1

  unsigned i,j,k;
  double * workweight;
  double * workdata_i_j;
  double * workdata_i_k;


  for (j=0;j<nrconst;j++)
    for(k=0;k<=j;k++)
      {
      X1(j,k) = 0;
      workweight = likep->get_weightp();
      workdata_i_j = data.getV()+j;
      workdata_i_k = data.getV()+k;
      for (i=0;i<likep->get_nrobs();i++,workweight++,workdata_i_j+=nrconst
                         ,workdata_i_k+=nrconst)
        X1(j,k) += *workweight * (*workdata_i_j) * (*workdata_i_k);
//        X1(j,k) += (likep->get_weight)(i,0)*data(i,j)*data(i,k);
      if (k!=j)
        X1(k,j) = X1(j,k);
      }

   X1 = X1.cinverse();

  // computing X2

  double * workhelp = help.getV();
  for (j=0;j<nrconst;j++)
    {
    workweight = likep->get_weightp();
    workdata_i_j = data.getV()+j;
    for(i=0;i<likep->get_nrobs();i++,workweight++,workhelp++,
                      workdata_i_j+=nrconst)
      *workhelp = *workweight * (*workdata_i_j);
//      help(j,i) = (likep->get_weight)(i,0)*data(i,j);

    }

  if (X1.rows() == nrconst)
    {
    X2.mult(X1,help);

    X1 = X1.root();
    }

  }


void FULLCOND_const_gaussian::update_intercept(double & m)
  {
  interceptadd+=m;
  beta(interceptpos,0) +=m;
  }


void FULLCOND_const_gaussian::update(void)
  {

  FULLCOND_const::update();  

  if (changingweight || optionsp->get_nriter()==1)
    {
    compute_matrices();
    }

  unsigned i;
  double * worklinold=linold.getV();
  for(i=0;i<linold.rows();i++,worklinold++)
    *worklinold += interceptadd;
  interceptadd=0;

  likep->substr_linearpred_m(linold,column);

  likep->compute_respminuslinpred(mu1,column);

  beta.mult(X2,mu1);
  beta+= sqrt(likep->get_scale(column))*X1*rand_normvek(nrconst);

  linold.mult(data,beta);
  likep->add_linearpred_m(linold,column);

  acceptance++;

  transform = likep->get_trmult(column);

//  FULLCOND_const::update();

  }


void FULLCOND_const_gaussian::posteriormode_intercept(double & m)
  {
  interceptadd+=m;
  beta(interceptpos,0) +=m;
  betameanold(interceptpos,0) +=m;
  }


bool FULLCOND_const_gaussian::posteriormode(void)
  {
  unsigned i;
  double * worklinold=linold.getV();
  for(i=0;i<linold.rows();i++,worklinold++)
    *worklinold += interceptadd;
  interceptadd=0;
  likep->fisher(X1,data,column);
  X1.assign((X1.cinverse()));
  likep->substr_linearpred_m(linold,column);
  likep->compute_weightiwls_workingresiduals(column);
  beta = X1*data.transposed()*likep->get_workingresiduals();
  linold.mult(data,beta);
  likep->add_linearpred_m(linold,column);
  return FULLCOND_const::posteriormode();
  }


bool FULLCOND_const_gaussian::posteriormode_converged(const unsigned & itnr)
  {
  return likep->posteriormode_converged_fc(beta,beta_mode,itnr);
  }


void FULLCOND_const_gaussian::outresults(void)
  {
  acceptance=optionsp->get_nriter();
  FULLCOND_const::outresults();
  }


//------------------------------------------------------------------------------
//---- CLASS FULLCOND_const_nongaussian: implementation of member functions ----
//------------------------------------------------------------------------------


FULLCOND_const_nongaussian::FULLCOND_const_nongaussian(MCMCoptions* o,
                            DISTRIBUTION * dp, const datamatrix & d,
                            const ST::string & t, const int & constant,
                            const ST::string & fs,const ST::string & fr,
                            const unsigned & c)
                            : FULLCOND_const(o,dp,d,t,constant,fs,fr,c)
  {

  step = o->get_step();
  diff = linnew;
  weightiwls = datamatrix(likep->get_nrobs(),1,1);
  tildey = weightiwls;
  proposal = beta;
  XWX = datamatrix(nrconst,nrconst);
  XWXold = XWX;
  help = beta;
  muy=datamatrix(nrconst,1);

  compute_XWX(XWXold);
  datamatrix test = XWXold.cinverse();
  if (test.rows() < nrconst)
    errors.push_back("ERROR: design matrix for fixed effects is rank deficient\n");

  }

FULLCOND_const_nongaussian::FULLCOND_const_nongaussian(
                          const FULLCOND_const_nongaussian & m)
  : FULLCOND_const(FULLCOND_const(m))
  {
  step = m.step;
  XWX = m.XWX;
  XWXold = m.XWXold;
  diff = m.diff;
  tildey = m.tildey;
  weightiwls = m.weightiwls;
  proposal = m.proposal;
  help = m.help;
  muy = m.muy;
  }


const FULLCOND_const_nongaussian & FULLCOND_const_nongaussian::operator=(
                           const FULLCOND_const_nongaussian & m)
  {
  if (this == &m)
    return *this;
  FULLCOND_const::operator=(FULLCOND_const(m));
  step = m.step;
  XWX = m.XWX;
  XWXold = m.XWXold;
  diff = m.diff;
  tildey = m.tildey;
  weightiwls = m.weightiwls;
  proposal = m.proposal;
  help = m.help;
  muy = m.muy;
  return *this;
  }


void FULLCOND_const_nongaussian::posteriormode_intercept(double & m)
  {
  interceptadd+=m;
  beta(interceptpos,0) +=m;
  betameanold(interceptpos,0) +=m;
  }


bool FULLCOND_const_nongaussian::posteriormode(void)
  {

  likep->fisher(XWX,data,column);
  XWX.assign(XWX.cinverse());

  unsigned i;
  double * worklinold=linold.getV();
  for(i=0;i<linold.rows();i++,worklinold++)
    *worklinold += interceptadd;
  interceptadd=0;


  likep->substr_linearpred_m(linold,column);

  likep->compute_weightiwls_workingresiduals(column);

  beta = XWX*data.transposed()*likep->get_workingresiduals();

  linnew.mult(data,beta);

  likep->add_linearpred_m(linnew,column);
  linold.assign(linnew);

  return FULLCOND_const::posteriormode();

  }


bool FULLCOND_const_nongaussian::posteriormode_converged(const unsigned & itnr)
  {
  return likep->posteriormode_converged_fc(beta,beta_mode,itnr);
  }


void FULLCOND_const_nongaussian::update_intercept(double & m)
  {
  interceptadd+=m;
  beta(interceptpos,0) +=m;
  }


void FULLCOND_const_nongaussian::compute_XWX(datamatrix & XWXw)
  {

  unsigned p,k,i;

  double * workw = weightiwls.getV();
  double*workXp;
  double* workXk;

  for (p=0;p<nrconst;p++)
    for (k=p;k<nrconst;k++)
      {
      XWXw(p,k)=0;
      workw = weightiwls.getV();
      workXp = data.getV()+p;
      workXk = data.getV()+k;
      for(i=0;i<weightiwls.rows();i++,workw++,workXp+=nrconst,workXk+=nrconst)
        XWXw(p,k)+= *workw  *  *workXp * *workXk;
      XWXw(k,p) = XWXw(p,k);
      }

  }


void FULLCOND_const_nongaussian::compute_XWtildey(datamatrix * linb)
  {
  unsigned i,j;

  double * worktildey=tildey.getV();
  double * worklinb = linb->getV();
  double * workw = weightiwls.getV();
  double * workdata = data.getV();
  double h;

  muy = datamatrix(nrconst,1,0);

  for(i=0;i<tildey.rows();i++,worktildey++,worklinb++,workw++)
    {
    h = *workw * (*worktildey + *worklinb);
    for (j=0;j<nrconst;j++,workdata++)
      muy(j,0) += *workdata * h;
    }

  }


void  FULLCOND_const_nongaussian::update(void)
  {

//  update_iwls();


  FULLCOND_const::update();

  double qoldbeta;
  double qnewbeta;

  if (optionsp->get_nriter() == 1)
    {

    linoldp = &linold;
    linnewp = &linnew;
    linold.mult(data,beta);
    linmode = datamatrix(data.rows(),1);
    mode =beta;
    }


  if (interceptyes && interceptadd!=0)
    {
    unsigned i;

    double * work = linoldp->getV();
    for (i=0;i<linoldp->rows();i++,work++)
      *work += interceptadd;
    interceptadd=0;
    }


  double logold = likep->loglikelihood();

  linmode.mult(data,mode);
  diff.minus(linmode,*linoldp);
  likep->add_linearpred_m(diff,column);

  likep->compute_IWLS_weight_tildey(weightiwls,tildey,column);
  compute_XWX(XWXold);
  XWX.assign((XWXold.cinverse()).root());

  compute_XWtildey(&linmode);

  mode = XWXold.solve(muy);

  help.minus(beta,mode);
  qoldbeta = -0.5*XWXold.compute_quadform(help);

  proposal.plus(mode,XWX*rand_normvek(nrconst));

  help.minus(proposal,mode);

  qnewbeta = -0.5*XWXold.compute_quadform(help);

  linnewp->mult(data,proposal);

  diff.minus(*linnewp,linmode);

  likep->add_linearpred_m(diff,column);          // (mit proposed)

  double logprop = likep->loglikelihood();     // mit proposed


  double u = log(uniform());

  if (u <= (logprop + qoldbeta - logold - qnewbeta) )
    {
    datamatrix * mp = linoldp;
    linoldp = linnewp;
    linnewp = mp;
//    linold.assign(linnew);
    beta.assign(proposal);

//    likep->swap_linearpred();

    acceptance++;
    }
  else
    {
    diff.minus(*linoldp,*linnewp);
    likep->add_linearpred_m(diff,column);
    }



  }


void  FULLCOND_const_nongaussian::update_iwls(void)
  {

  FULLCOND_const::update();

  double qold;
  double qnew;

  if (optionsp->get_nriter() == 1)
    {

    linoldp = &linold;
    linnewp = &linnew;
    linold.mult(data,beta);
    mode =beta;
    }


  if (interceptyes && interceptadd!=0)
    {
    unsigned i;

    double * work = linoldp->getV();
    for (i=0;i<linoldp->rows();i++,work++)
      *work += interceptadd;
    interceptadd=0;
    }


  double logold = likep->compute_IWLS(weightiwls,tildey,true,column);
  compute_XWX(XWXold);
  XWX.assign((XWXold.cinverse()).root());

// eigene Version!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  linoldp->mult(data,beta);

  compute_XWtildey(linoldp);



  mode = XWXold.solve(muy);

  proposal.plus(mode,XWX*rand_normvek(nrconst));

// eigene Version!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  proposal(interceptpos, 0) = beta(interceptpos, 0);

  help.minus(proposal,mode);

  qold = XWXold.compute_quadform(help);
  qold += log(1.0/XWXold.det());
  qold *= -0.5;


  linnewp->mult(data,proposal);

  diff.minus(*linnewp,*linoldp);

  likep->addtocurrentcol(diff,column);                // (mit proposed)

  // mir (proposed)
  double logprop = likep->compute_IWLS(weightiwls,tildey,true,column,false);
  compute_XWX(XWXold);
//  XWX.assign((XWXold.cinverse()).root());

  compute_XWtildey(linnewp);

  mode = XWXold.solve(muy);

  help.minus(beta,mode);

  qnew = XWXold.compute_quadform(help);
  qnew += log(1.0/XWXold.det());
  qnew *= -0.5;

  double u = log(uniform());

  if (u <= (logprop + qnew - logold - qold) )
    {
    datamatrix * mp = linoldp;
    beta.assign(proposal);
    linoldp = linnewp;
    linnewp = mp;

    likep->swap_linearpred();

    acceptance++;
    }


  }


//------------------------------------------------------------------------------
//---- CLASS FULLCOND_const_nbinomial: implementation of member functions ----
//------------------------------------------------------------------------------


FULLCOND_const_nbinomial::FULLCOND_const_nbinomial(MCMCoptions* o,
                            DISTRIBUTION * dp, DISTRIBUTION_nbinomial * nb,const datamatrix & d,
                            const ST::string & t, const int & constant,
                            const ST::string & fs,const ST::string & fr,
                            const unsigned & c)
                            : FULLCOND_const_nongaussian(o,dp,d,t,constant,fs,fr,c)
  {
  nblikep = nb;
  }

FULLCOND_const_nbinomial::FULLCOND_const_nbinomial(
                          const FULLCOND_const_nbinomial & m)
  : FULLCOND_const_nongaussian(FULLCOND_const_nongaussian(m))
  {
  nblikep = m.nblikep;
  }


const FULLCOND_const_nbinomial & FULLCOND_const_nbinomial::operator=(
                           const FULLCOND_const_nbinomial & m)
  {
  if (this == &m)
    return *this;
  FULLCOND_const_nongaussian::operator=(FULLCOND_const_nongaussian(m));

  nblikep = m.nblikep;

  return *this;
  }


bool FULLCOND_const_nbinomial::posteriormode_converged(const unsigned & itnr)
  {
  if(FULLCOND_const_nongaussian::posteriormode_converged(itnr))
  nblikep->initialize_hierint(beta(interceptpos,0));
  return FULLCOND_const_nongaussian::posteriormode_converged(itnr);
  }


void FULLCOND_const_nbinomial::update_intercept(double & m)
  {
  FULLCOND_const_nongaussian::update_intercept(m);
  nblikep->exchange_hierint(m);
  nblikep->add_nu(exp(m));
  }


void FULLCOND_const_nbinomial::update(void)
  {
/*
    Mit der Funktion geht es nicht! die Acceptance-Raten für FixedEffect
    sind viel zu klein!

  FULLCOND_const_nongaussian::update();

*/

//  FULLCOND_const_nongaussian::update_iwls();


  FULLCOND_const::update();

  double qold;
  double qnew;

  if (optionsp->get_nriter() == 1)
    {

    linoldp = &linold;
    linnewp = &linnew;
    linold.mult(data,beta);
    mode =beta;
    }


  if (interceptyes && interceptadd!=0)
    {
    unsigned i;

    double * work = linoldp->getV();
    for (i=0;i<linoldp->rows();i++,work++)
      *work += interceptadd;
    interceptadd=0;
    }


  double logold = likep->compute_IWLS(weightiwls,tildey,true,column);
  compute_XWX(XWXold);
  XWX.assign((XWXold.cinverse()).root());

// eigene Version!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  linoldp->mult(data,beta);




  compute_XWtildey(linoldp);

// eigene Version!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  mode = XWXold.solve(muy);

  proposal.plus(mode,XWX*rand_normvek(nrconst));

  proposal(interceptpos, 0) = beta(interceptpos, 0);

  help.minus(proposal,mode);


//    proposal.plus(beta,XWX*rand_normvek(nrconst));

//    proposal(interceptpos, 0) = beta(interceptpos, 0);

//    help.minus(proposal,beta);

  qold = XWXold.compute_quadform(help);
  qold += log(1.0/XWXold.det());
  qold *= -0.5;


  linnewp->mult(data,proposal);

  diff.minus(*linnewp,*linoldp);

  likep->addtocurrentcol(diff,column);                // (mit proposed)

  // mir (proposed)
  double logprop = likep->compute_IWLS(weightiwls,tildey,true,column,false);
  compute_XWX(XWXold);
//  XWX.assign((XWXold.cinverse()).root());

  compute_XWtildey(linnewp);

// eigene Version!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  mode = XWXold.solve(muy);

  help.minus(beta,mode);

//    help.minus(beta,proposal);


  qnew = XWXold.compute_quadform(help);
  qnew += log(1.0/XWXold.det());
  qnew *= -0.5;

  double u = log(uniform());

  if (u <= (logprop + qnew - logold - qold) )
    {
    datamatrix * mp = linoldp;
    beta.assign(proposal);
    linoldp = linnewp;
    linnewp = mp;

    likep->swap_linearpred();

    acceptance++;
    }


    beta(interceptpos, 0) = update_hierint();

  }


double FULLCOND_const_nbinomial::update_hierint(void) const
{
//    unsigned i;
    double intaux;
    double intwork = nblikep->get_hierint();
    double scalework = likep->get_scale();
    double sum = nblikep->get_sum_nu();
    double sum2 = nblikep->get_sum2_nu();
    int ver = nblikep->get_distopt();
    unsigned nrobs = likep->get_nrobs();


    if(ver==1)
    {
        intaux = intwork;
        intwork = -log(randnumbers::rand_gamma( scalework*nrobs,
                        scalework*sum));
        double diff = intwork - intaux;
        likep->add_linearpred_m(diff,0);
        nblikep->exchange_hierint(diff);
        diff = exp(diff);
        nblikep->add_nu(diff);


    }
    else  // if poig
    {

        double intpost;
        double priori_ratio=0.0; // Weil diffuse Prior
        double likelihood_ratio;
        double proposal_ratio=0.0; // Weil unabhängig von aktuellem Zustand
        double pwork=nblikep->get_pvar();

        intpost = (nrobs*0.5-1+
               sqrt((nrobs*0.5-1)*(nrobs*0.5-1)+ scalework* scalework * sum *sum2))/
               (scalework* sum2);
        intpost = log(intpost);
        if(intpost > pwork) intaux = intpost-pwork + 2*pwork*randnumbers::uniform();
        else intaux = (intpost + pwork)*(randnumbers::uniform());

        likelihood_ratio = nrobs*(intaux - intwork)*0.5+
                       scalework*0.5*(1/exp(intwork)-1/exp(intaux))*sum+
                       scalework*0.5*(exp(intwork)-exp(intaux))*sum2;

        double alpha = likelihood_ratio + priori_ratio + proposal_ratio;
        double h =log(randnumbers::uniform( ));
        if(h<=alpha)
        {
            double diff = intaux - intwork;
            likep->add_linearpred_m(diff,0);
            nblikep->exchange_hierint(diff);
            diff = exp(diff);
            nblikep->add_nu(diff);
            intwork = intaux;
            nblikep->exchange_accept();

        }
        
    }

    return intwork;
}





} // end: namespace MCMC






