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

/*
unsigned FULLCOND_mixture::get_rankK(void)
  {
  if (randomslope && includefixed)
    return nrpar-1;
  else
    return nrpar;
  }
*/

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
  compind = statmatrix<int>(d.rows(),1,1);
  compmean = datamatrix(nrcomp,1,0);
  compvar = datamatrix(nrcomp,1,1);
  compweights = datamatrix(nrcomp,1,1.0/nrcomp);
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

  muy = datamatrix(nrpar,1);

//  identifiable =true;

  }


FULLCOND_mixture::FULLCOND_mixture(const FULLCOND_mixture & fc)
                            : FULLCOND(FULLCOND(fc))
  {
  nrcomp = fc.nrcomp;
  compind=fc.compind;
  compmean=fc.compmean;
  compvar=fc.compvar;
  compweights=fc.compweights;

  muy = fc.muy;
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
  compind=fc.compind;
  compmean=fc.compmean;
  compvar=fc.compvar;
  compweights=fc.compweights;

  muy = fc.muy;
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

  if (randomslope && includefixed)
    {
    optionsp->out("  Fixed effect:\n");
    optionsp->out("\n");

    ST::string help =  ST::doubletostring(lower1,4) + "% quant.";
    ST::string levell = help + ST::string(' ',15-help.length());
    help = ST::doubletostring(upper2,4) + "% quant.";
    ST::string levelu = help + ST::string(' ',15-help.length());
    help = ST::string(' ',+2);

    optionsp->out(help +
                    "mean           " +
                    "Std. Dev.      " +
                    levell +
                    "median         " +
                    levelu + "\n");

    optionsp->out(ST::outresults(0,"",betamean(nrpar-1,0),
                  sqrt(betavar(nrpar-1,0)),betaqu_l1_lower(nrpar-1,0),
                  betaqu50(nrpar-1,0),betaqu_l1_upper(nrpar-1,0)));

    optionsp->out("\n");

    optionsp->out("  Results for the fixed effect are also stored in file \n");
    optionsp->out("  " + pathcurrent2 + "\n");

    optionsp->out("\n");

    ofstream outfixed(pathcurrent2.strtochar());

    outfixed << "pmean   ";
    outfixed << "pqu"  << nl1  << "   ";
    outfixed << "pqu"  << nl2  << "   ";
    outfixed << "pmed   ";
    outfixed << "pqu"  << nu1  << "   ";
    outfixed << "pqu"  << nu2  << "   ";
    outfixed << "pcat" << level1 << "   ";
    outfixed << "pcat" << level2 << "   ";
    outfixed << endl;

    outfixed << betamean(nrpar-1,0) << "   ";
    outfixed << betaqu_l1_lower(nrpar-1,0) << "   ";
    outfixed << betaqu_l2_lower(nrpar-1,0) << "   ";
    outfixed << betaqu50(nrpar-1,0) << "   ";
    outfixed << betaqu_l2_upper(nrpar-1,0) << "   ";
    outfixed << betaqu_l1_upper(nrpar-1,0) << "   ";
    if (betaqu_l1_lower(nrpar-1,0) > 0)
      outfixed << "1   ";
    else if (betaqu_l1_upper(nrpar-1,0) < 0)
      outfixed << "-1   ";
    else
      outfixed << "0   ";
    if (betaqu_l2_lower(nrpar-1,0) > 0)
      outfixed << "1   ";
    else if (betaqu_l2_upper(nrpar-1,0) < 0)
      outfixed << "-1   ";
    else
      outfixed << "0   ";

    outfixed << endl;
    }

  if (randomslope)
    {
    optionsp->out("  Results for random slopes are stored in file\n");
    optionsp->out("  " + pathcurrent + "\n");
    }
  else
    {
    optionsp->out("  Results for random effects are stored in file\n");
    optionsp->out("  " + pathcurrent + "\n");

  if (lambdaconst==true)
    {
    optionsp->out("\n");
    optionsp->out("  Constant smoothing parameter: " +
    ST::doubletostring(lambda,6) + "\n");
    optionsp->out("\n");
    }

    optionsp->out("\n");

    for(unsigned i=0;i<nrcomp;i++)
       {
       optionsp->out("  Component means: " + ST::doubletostring(compmean(i,0),6) + "\n");
       optionsp->out("  Component variances: " + ST::doubletostring(compvar(i,0),6) + "\n");
       optionsp->out("  Component weights: " + ST::doubletostring(compweights(i,0),6) + "\n");
       }
    }

  if (optionsp->get_samplesize() == 0)
    {
    optionsp->out("\n");
    double df = compute_df();
    optionsp->out("  Approximate degrees of freedom: "
                    + ST::doubletostring(df,6) + "\n");
    }

  optionsp->out("\n");

  unsigned i;

  ofstream outres(pathcurrent.strtochar());
  assert(!outres.fail());

  ST::string name = datanames[0];

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

  double * workmean = betamean.getV();
  double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
  double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
  double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
  double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
  double * workbetaqu50 = betaqu50.getV();

  for(i=0;i<nrpar;i++,workmean++,workbetaqu_l1_lower_p++,
      workbetaqu_l2_lower_p++,workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++,
      workbetaqu50++)
    {
    if (randomslope && includefixed && i == nrpar-1)
      {
      }
    else
      {
      outres << (i+1) << "   ";
      outres << effvalues(i,0) << "   ";
      outres << *workmean << "   ";
      outres << *workbetaqu_l1_lower_p << "   ";
      outres << *workbetaqu_l2_lower_p << "   ";
      outres << *workbetaqu50 << "   ";
      outres << *workbetaqu_l2_upper_p << "   ";
      outres << *workbetaqu_l1_upper_p << "   ";

      if (*workbetaqu_l1_lower_p > 0)
        outres << "1   ";
      else if (*workbetaqu_l1_upper_p < 0)
        outres << "-1   ";
      else
        outres << "0   ";

      if (*workbetaqu_l2_lower_p > 0)
        outres << "1   ";
      else if (*workbetaqu_l2_upper_p < 0)
        outres << "-1   ";
      else
        outres << "0   ";

      outres << endl;
      }
    }
  }


void FULLCOND_mixture::outoptions(void)
  {

  optionsp->out("  OPTIONS FOR RANDOM EFFECT MIXTURE: " + title + "\n",true);
  optionsp->out("\n");

  optionsp->out("  Type of Mixture: Normal\n",true);
  optionsp->out("  Number of components: " + ST::inttostring(nrcomp) + "\n",true);
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
    if (!randomslope)
      {
//      likep->set_linpredp_current(column);
      for (i=0;i<nrpar;i++,++itbeg,++itend,workbeta++)
        if (*itbeg != -1)
          likep->add_linearpred2(-(*workbeta),*itbeg,*itend,index,index2,column);
      }
    else
      {
      workindex2 = index2.getV();
      double * datap = data.getV();
      if (includefixed)
        {
        n = nrpar-1;
        double ms = beta(nrpar-1,0);
        double h;
        likep->set_linpredp_current(column);
        for (i=0;i<n;i++,++itbeg,++itend,workbeta++)
          {
          if (*itbeg != -1)
            {
            h = *workbeta+ms;
            for(j=*itbeg;j<=*itend;j++,workindex2++,datap++)
              likep->add_linearpred2(-h*(*datap),*workindex2);
            }
          }
        }
      else
        {
        n = nrpar;
        likep->set_linpredp_current(column);
        for (i=0;i<n;i++,++itbeg,++itend,workbeta++)
          {
          if (*itbeg != -1)
            {
            for(j=*itbeg;j<=*itend;j++,workindex2++,datap++)
              {
              likep->add_linearpred2(-*workbeta*(*datap),*workindex2);
              }
            }

          }
        }

      }


    } // end: if (add == false)
  else
    {

    if (!randomslope)
      {
//      likep->set_linpredp_current(column);
      for (i=0;i<nrpar;i++,++itbeg,++itend,workbeta++)
        if (*itbeg != -1)
          likep->add_linearpred2(*workbeta,*itbeg,*itend,index,index2,column);
      }
    else
      {
      workindex2 = index2.getV();
      double * datap = data.getV();
      if (includefixed)
        {
        n = nrpar-1;
        double ms = beta(nrpar-1,0);
        double h;
        likep->set_linpredp_current(column);
        for (i=0;i<n;i++,++itbeg,++itend,workbeta++)
          {
          if (*itbeg != -1)
            {
            h = *workbeta+ms;
            for(j=*itbeg;j<=*itend;j++,workindex2++,datap++)
              likep->add_linearpred2(h*(*datap),*workindex2);
            }
          }
        }
      else
        {
        n = nrpar;
        likep->set_linpredp_current(column);
        for (i=0;i<n;i++,++itbeg,++itend,workbeta++)
          {
          if (*itbeg != -1)
            {
            for(j=*itbeg;j<=*itend;j++,workindex2++,datap++)
              {
              likep->add_linearpred2(*workbeta*(*datap),*workindex2);
              }
            }
          }
        }

      }

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
    if (includefixed)
      {
      double ms = beta(nrpar-1,0);
      likep->set_linpredp_current(column);
      for (i=0;i<n;i++,workmuy++,++itbeg,++itend)
        {
        *workmuy = 0;
        for(j=*itbeg;j<=*itend;j++,workindex2++,datap++)
          *workmuy += likep->get_workingres(*workindex2)* (*datap);

        *workmuy+= lambda* ms;
        }
      }
    else
      {
      for(i=0;i<n;i++,workmuy++,++itbeg,++itend)
        {
        *workmuy = 0;
        for(j=*itbeg;j<=*itend;j++,workindex2++,datap++)
          {
          *workmuy+= likep->get_workingres(*workindex2)* (*datap);
          }

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


  if (randomslope && includefixed)
    {
    double * workbeta = beta.getV();
    double sum=0;
    for (i=0;i<n;i++,workbeta++)
      {
      sum += *workbeta;
      }

    beta(nrpar-1,0) = sum/double(n);

    workbeta = beta.getV();
    double ms = beta(nrpar-1,0);
    for (i=0;i<n;i++,workbeta++)
      *workbeta -= ms;

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

  double var;
  double m;
  unsigned i,j;
  unsigned n = nrpar;

  if (randomslope && includefixed)
  n = nrpar-1;


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

  if (!randomslope)
    {
    for(i=0;i<nrpar;i++,workmuy++,++itbeg,++itend)
      {
      *workmuy = 0;
      for(j=*itbeg;j<=*itend;j++,workindex2++)
        {
        mup += *workindex2;
        *workmuy+= likep->get_weight(*workindex2)* *mup;
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
        mup += *workindex2;
        *workmuy+= likep->get_weight(*workindex2)* (*mup) * (*datap);
        }

      if (includefixed)
        *workmuy += beta(n,0)*lambda;

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


  if (randomslope && includefixed)
    {

    workbeta = beta.getV();
    double s=0;
    for (i=0;i<nrpar-1;i++,workbeta++)
      s += *workbeta;
    s /= double(nrpar-1);

    double v = sigma2/double(nrpar-1);

    beta(nrpar-1,0) = s+sqrt(v)*rand_normal();

    workbeta = beta.getV();
    double ms = beta(nrpar-1,0);
    for (i=0;i<nrpar-1;i++,workbeta++)
      *workbeta -= ms;

    }


  update_linpred(true);


  if (center)
    {
    double m = centerbeta();
    fcconst->update_intercept(m);
    }


  acceptance++;

  transform = likep->get_trmult(column);


  FULLCOND_mixture::update();

  if (spatialtotal)
    {
    double * ftotal_bp = ftotal.getbetapointer();
    workbeta=beta.getV();
    double * workbetaspat = fbasisp->getbetapointer();
    int * indexp = indextotal.getV();
    for (i=0;i<nrpar;i++,workbeta++,ftotal_bp++,indexp++)
      {
      workbetaspat+= *indexp;
      *ftotal_bp = *workbeta + *workbetaspat;
      }

    ftotal.set_transform(likep->get_trmult(column));

    ftotal.update();
    }

  }


} // end: namespace MCMC
