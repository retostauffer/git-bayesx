// DATE: 11.8.99

#include "fullcond.h"
#include <clstring.h>

//------------------------------------------------------------------------------
//--------------- CLASS: FULLCOND implementation of member functions -----------
//------------------------------------------------------------------------------


namespace MCMC
{


FULLCOND::FULLCOND(void)
  {
  }


FULLCOND::FULLCOND(MCMCoptions * o,const datamatrix & d,
                   const ST::string & t,const unsigned & rows,
                   const unsigned & cols,const ST::string & fp)
  {

  df_accuracy = 0.05;

  transformnonlinear = false;
  transformed =  false;

  optionsp = o;

  plotstyle = noplot;

  term_symbolic="";

  results_type = "";

  title = t;

  samplepath = fp;

  flags = bitset<flagnr>("000");

  data = d;

  for (unsigned i=0;i<data.cols();i++)
    datanames.push_back("X" + ST::inttostring(i));

  level1 = o->get_level1();
  level2 = o->get_level2();

  lower1 = (100.0-level1)/2;
  lower2 = (100.0-level2)/2;

  upper1 = 100.0 - lower2;
  upper2 = 100.0 - lower1;

  setbeta(rows,cols,0);

  acceptance = 0;
  nrtrials = 0;

  identifiable = true;
  center = false;
  baseline = false;

  column = 0;

  transform = 1;
  transformmult = datamatrix(1,1,1);
  addon = 0;

  }

// for REML
FULLCOND::FULLCOND(MCMCoptions * o,const ST::string & t)
  {

  title = t;

  column=0;

  optionsp = o;

  level1 = o->get_level1();
  level2 = o->get_level2();

  lower1 = (100.0-level1)/2;
  lower2 = (100.0-level2)/2;

  upper1 = 100.0 - lower2;
  upper2 = 100.0 - lower1;

  }

FULLCOND::FULLCOND(const FULLCOND & m)
  {

  //------------------------------ for stepwise --------------------------------

  lambdastart=m.lambdastart;
  lambdamin = m.lambdamin;
  lambdamax = m.lambdamax;
  data_forfixed = m.data_forfixed;

  //---------------------------- end: for stepwise -----------------------------

  fctype = m.fctype;

  optionsp = m.optionsp;

  title = m.title;

  samplepath = m.samplepath;
  pathresult = m.pathresult;
  pathresult2 = m.pathresult2;
  pathresult3 = m.pathresult3;
  pathcurrent = m.pathcurrent;
  pathcurrent2 = m.pathcurrent2;
  pathcurrent3 = m.pathcurrent3;

  fcnumber = m.fcnumber;
  plotstyle = m.plotstyle;

  term_symbolic = m.term_symbolic;
  priorassumptions = m.priorassumptions;

  flags = m.flags;

  nrpar = m.nrpar;

  data = m.data;
  datanames = m.datanames;

  errors = m.errors;

  beta = m.beta;
  beta_mode = m.beta_mode;
  betamean = m.betamean;
  betas2 = m.betas2;
  betavar = m.betavar;
  betamin = m.betamin;
  betamax = m.betamax;

  betaqu_l1_lower = m.betaqu_l1_lower;
  betaqu_l2_lower = m.betaqu_l2_lower;
  betaqu_l1_upper = m.betaqu_l1_upper;
  betaqu_l2_upper = m.betaqu_l2_upper;
  betaqu50 = m.betaqu50;
  level1 = m.level1;
  level2 = m.level2;
  lower1 = m.lower1;
  lower2 = m.lower2;
  upper1 = m.upper1;
  upper2 = m.upper2;

  betameanold = m.betameanold;
  betavarold = m.betavarold;
  betaminold = m.betaminold;
  betamaxold = m.betamaxold;

  transform = m.transform;
  transformmult = m.transformmult;
  addon = m.addon;

  results_latex = m.results_latex;
  results_type = m.results_type;

  transformnonlinear = m.transformnonlinear;
  transformed = m.transformed;
  transformtype = m.transformtype;

  acceptance = m.acceptance;
  nrtrials = m.nrtrials;
  identifiable = m.identifiable;
  center = m.center;
  baseline = m.baseline;

  column = m.column;

  weight = m.weight;

  //------------------------------ for REML ------------------------------------

  dimX = m.dimX;
  dimZ = m.dimZ;
  startlambda = m.startlambda;
  isnonparametric = m.isnonparametric;

  }


const FULLCOND & FULLCOND::operator=(const FULLCOND & m)
  {

  if (this==&m)
	 return *this;

  //------------------------------ for stepwise --------------------------------

  lambdastart=m.lambdastart;
  lambdamin = m.lambdamin;
  lambdamax = m.lambdamax;
  data_forfixed = m.data_forfixed;

  //------------------------------ end: stepwise -------------------------------

  fctype = m.fctype;

  optionsp = m.optionsp;

  title = m.title;

  samplepath = m.samplepath;

  pathresult = m.pathresult;
  pathresult2 = m.pathresult2;
  pathresult3 = m.pathresult3;

  pathcurrent = m.pathcurrent;
  pathcurrent2 = m.pathcurrent2;
  pathcurrent3 = m.pathcurrent3;

  fcnumber = m.fcnumber;
  plotstyle = m.plotstyle;

  term_symbolic = m.term_symbolic;
  priorassumptions = m.priorassumptions;

  flags = m.flags;

  nrpar = m.nrpar;

  data = m.data;
  datanames = m.datanames;

  errors = m.errors;

  beta = m.beta;
  beta_mode = m.beta_mode;
  betamean = m.betamean;
  betas2 = m.betas2;
  betavar = m.betavar;
  betamin = m.betamin;
  betamax = m.betamax;

  betaqu_l1_lower = m.betaqu_l1_lower;
  betaqu_l2_lower = m.betaqu_l2_lower;
  betaqu_l1_upper = m.betaqu_l1_upper;
  betaqu_l2_upper = m.betaqu_l2_upper;
  betaqu50 = m.betaqu50;
  level1 = m.level1;
  level2 = m.level2;
  lower1 = m.lower1;
  lower2 = m.lower2;
  upper1 = m.upper1;
  upper2 = m.upper2;

  betameanold = m.betameanold;
  betavarold = m.betavarold;
  betaminold = m.betaminold;
  betamaxold = m.betamaxold;

  transform = m.transform;
  transformmult = m.transformmult;
  addon = m.addon;

  results_latex = m.results_latex;
  results_type = m.results_type;

  transformnonlinear = m.transformnonlinear;
  transformed = m.transformed;
  transformtype = m.transformtype;

  acceptance = m.acceptance;
  nrtrials = m.nrtrials;

  identifiable = m.identifiable;
  center = m.center;
  baseline = m.baseline;

  column = m.column;

  weight = m.weight;

  //------------------------------ for REML ------------------------------------

  dimX = m.dimX;
  dimZ = m.dimZ;
  startlambda = m.startlambda; 
  isnonparametric = m.isnonparametric;

  return *this;

  }


void FULLCOND::set_transform(ST::string & suffix,ST::string & trtype)
  {
  transformnonlinear = true;
  transformtype=trtype;

  if (pathresult.length()>0)
    {

    int i=pathresult.length();
    bool found=false;
    while ( (i>0) && (found==false) )
      {
      i--;
      if (pathresult[i] == '.')
        found=true;
      }

    if (found==true)
      {
      pathcurrent = pathresult.substr(0,i) + suffix;
      }

    }

  if (pathresult2.length()>0)
    {

    int i=pathresult2.length();
    bool found=false;
    while ( (i>0) && (found==false) )
      {
      i--;
      if (pathresult2[i] == '.')
        found=true;
      }

    if (found==true)
      {
      pathcurrent2 = pathresult2.substr(0,i) + suffix;
      }

    }

  if (pathresult3.length()>0)
    {

    int i=pathresult3.length();
    bool found=false;
    while ( (i>0) && (found==false) )
      {
      i--;
      if (pathresult3[i] == '.')
        found=true;
      }

    if (found==true)
      {
      pathcurrent3 = pathresult3.substr(0,i) + suffix;
      }

    }


  }


void FULLCOND::outerrors(void)
  {
  unsigned i;
  for (i=0;i<errors.size();i++)
    optionsp->outerror(errors[i]);
  }


void FULLCOND::setbeta(const unsigned & rows,const unsigned & cols,
                       const double & v)
  {
  assert(rows > 0);
  assert(cols > 0);

  nrpar = rows*cols;
  beta = datamatrix(rows,cols,v);
  betamean = datamatrix(rows,cols,0);
  beta_mode = datamatrix(rows,cols,0);
  betas2 = datamatrix(rows,cols,0);
  betavar = datamatrix(rows,cols,0);
  betamin = datamatrix(rows,cols,0);
  betamax = datamatrix(rows,cols,0);
  betaqu_l1_lower = datamatrix(rows,cols,0);
  betaqu_l2_lower = datamatrix(rows,cols,0);
  betaqu_l1_upper = datamatrix(rows,cols,0);
  betaqu_l2_upper = datamatrix(rows,cols,0);
  betaqu50 = datamatrix(rows,cols,0);
  betameanold = datamatrix(rows,cols,0);
  betavarold = datamatrix(rows,cols,0);
  betaminold = datamatrix(rows,cols,0);
  betamaxold = datamatrix(rows,cols,0);
  }


void FULLCOND::setbetavalue(const unsigned & row,const unsigned & col,
                            const double & v)
  {
  beta(row,col) = v;
  }


void FULLCOND::setbeta(const datamatrix & betanew)
  {
  nrpar = betanew.rows()*betanew.cols();
  beta = betanew;
  betamean = beta;
  beta_mode = beta;
  betas2 = beta;
  betavar = beta;
  betamin = beta;
  betamax = beta;
  betaqu_l1_lower = beta;
  betaqu_l2_lower = beta;
  betaqu_l1_upper = beta;
  betaqu_l2_upper = beta;
  betaqu50 = beta;
  betameanold = beta;
  betavarold = beta;
  betaminold = beta;
  betamaxold = beta;
  }


double FULLCOND::centerbeta(void)
  {

  unsigned i;
  double sum=0;
  double * workbeta = beta.getV();

  #if defined(MICROSOFT_VISUAL)
  for (i=0;i<nrpar;i++,workbeta++)
    sum+= *workbeta;
  #else

  for (i=0;i<nrpar;i++,workbeta++)
    {
    sum+= *workbeta;
    }
  #endif

  workbeta = beta.getV();

  sum /= double(nrpar);

  for (i=0;i<nrpar;i++,workbeta++)
    *workbeta-= sum;

  return sum;

  }



double FULLCOND::centerbeta2(datamatrix & sumk1,datamatrix & sumk2)
  {
  unsigned i,j;


  double * workbeta = beta.getV();
  unsigned rK1 = sumk1.rows();
  unsigned rK2 = sumk2.rows();

  double * worksumk1=sumk1.getV();
  double * worksumk2;
  double sumtotal=0;

  for (j=0;j<rK1;j++,worksumk1++)
    {
    *worksumk1 = 0;
    worksumk2 = sumk2.getV();
    for (i=0;i<rK2;i++,workbeta++,worksumk2++)
      {
      if (j==0)
        *worksumk2 = 0;
      sumtotal+= *workbeta;
      *worksumk1 += *workbeta;
      *worksumk2 += (1.0/double(rK1))* *workbeta;
      }

    *worksumk1 = (*worksumk1)/double(rK2);

    }

  sumtotal/= double(rK1*rK2);

  workbeta = beta.getV();
  worksumk1 = sumk1.getV();

  for (j=0;j<rK1;j++,worksumk1++)
    {
    worksumk2 = sumk2.getV();
    for (i=0;i<rK2;i++,workbeta++,worksumk2++)
      *workbeta-= *worksumk1 + *worksumk2 - sumtotal;
    }

  return sumtotal;

  }



void FULLCOND::reset(void)
  {

  setbeta(beta.rows(),beta.cols(),0);
  acceptance = 0;
  nrtrials = 0;
  if (flags[0] == 0)
    {
    samplestream.close();
    remove(samplepath.strtochar());
    }
  errors.erase(errors.begin(),errors.end());
  }



void FULLCOND::readsample(datamatrix & sample,const unsigned & nr,
                          const unsigned & col) const
  {

  assert(nr < nrpar);
  assert(sample.rows() == optionsp->get_samplesize());



  unsigned size = sizeof(double);
  ifstream in;
  in.open(samplepath.strtochar(),ios::binary);
  unsigned i;
  double* work=sample.getV()+col;
  unsigned s = sample.cols();
  in.seekg(size*nr);
  for (i=0;i<optionsp->get_samplesize();i++,work+=s)
    {
//    in.seekg(size*nr+i*nrpar*size);
    in.read((char*) work,size);
    in.seekg(size*(nrpar-1),ios::cur);
    }
  }


void FULLCOND::readsample2(datamatrix & b,const unsigned & nr) const
  {

  assert(b.rows() == beta.rows());
  assert(b.cols() == beta.cols());

  unsigned size = sizeof(double);
  ifstream in;
  in.open(samplepath.strtochar(),ios::binary);
  double * work = b.getV();
  unsigned i;
  for (i=0;i<nrpar;i++,work++)
    {
    in.seekg(size*nrpar*nr + i*size);
    in.read((char*) work,size);
    }

  }


void FULLCOND::readsample3(datamatrix & b) const
  {

  assert(b.rows() == optionsp->get_samplesize());
  assert(b.cols() == nrpar);

  unsigned size = sizeof(double);
  ifstream in;
  in.open(samplepath.strtochar(),ios::binary);
  double * work = b.getV();
  unsigned i,j;
  for(i=0;i<b.rows();i++)
    for(j=0;j<b.cols();j++,work++)
      in.read((char *) work,size);

  }


void FULLCOND::setflags(const bitset<flagnr> & newflags)
  {
  flags = newflags;
  }


datamatrix FULLCOND::compute_autocorr(const unsigned & lag,const unsigned & row,
                                      const unsigned & col) const
  {
  assert(optionsp->get_samplesize() > 0);
  datamatrix sample(optionsp->get_samplesize(),1);
  unsigned nr = row*(beta.cols())+col;
  readsample(sample,nr);
  return sample.autocorr (1,lag,0);
  }


void FULLCOND::get_samples(const ST::string & filename,const unsigned & step) const
  {
  ofstream out(filename.strtochar());
  assert(!out.fail());

  unsigned i,j,k;

  datamatrix betac(beta.rows(),beta.cols());

  out << "intnr " << " ";
  if (beta.cols() > 1)
    {
    for (j=0;j<beta.rows();j+=step)
      for(k=0;k<beta.cols();k++)
        out << "b_" << (j+1) << "_" << (k+1) << " ";
    }
  else
    {
    for (j=0;j<nrpar;j+=step)
      out << "b_" << (j+1) << " ";
    }

  out << endl;

  for(i=0;i<optionsp->get_samplesize();i++)
    {
    readsample2(betac,i);
    out << (i+1) << " ";
    for (j=0;j<beta.rows();j+=step)
      for(k=0;k<betac.cols();k++)
        out << betac(j,k) << " ";

    out << endl;
    }

  out.close();

  }


void FULLCOND::get_covmatrix(datamatrix & r)
  {

  unsigned size = optionsp->get_samplesize();
  datamatrix b1(size,1);
  datamatrix b2(size,1);
  double m1,m2;
  double sum2;

  r = datamatrix(nrpar,nrpar,0);

  double * workmean = betamean.getV();
  double * workvar = betavar.getV();

  unsigned i,j,k;
  for(i=0;i<nrpar;i++,workvar++,workmean++)
    {
    readsample(b1,i);
    m1 = *workmean;
    r(i,i) = *workvar;
    for(j=i+1;j<nrpar;j++)
      {
      readsample(b2,j);
      m2 = b2.mean(0);

      sum2=0;
      for(k=0;k<size;k++)
        sum2 += b1(k,0)*b2(k,0);

      r(i,j) = (1.0/size)*sum2 - m1*m2;
      r(j,i) = r(i,j);

      }



    }

  }


void FULLCOND::get_covmatrix(const ST::string & file,const covstyle & st)
  {

  ofstream out(file.strtochar());
  assert(!out.fail());

  datamatrix r;
  get_covmatrix(r);

  if (st == precision)
    r = r.inverse();

  unsigned i,j;
  for (i=0;i<nrpar;i++)
    {
    for(j=0;j<nrpar;j++)
      {
      if (st == correlation)
        out << (r(i,j)/sqrt(r(i,i)*r(j,j))) << " ";
      else
        out << r(i,j) << " ";
      }
    out << endl;
    }

  }

void FULLCOND::update(void)
  {
  double diffmean;
  double diffvar;
  double diffmin;
  double diffmax;
  double normold;
  double rate;

  if(
        (optionsp->get_nriter() > optionsp->get_burnin())
     &&
        ((optionsp->get_nriter()-optionsp->get_burnin()-1) % (optionsp->get_step()) == 0)
    )
    {

    register unsigned i;
    double* workbeta = beta.getV();
    double* workbetamean = betamean.getV();
    double* workbetas2 = betas2.getV();
    double* workbetavar = betavar.getV();
    double* workbetamin = betamin.getV();
    double* workbetamax = betamax.getV();

    unsigned samplesize = optionsp->get_samplesize();
    if ((flags[0] != 1) && (samplesize == 1))
      {
      samplestream.open(samplepath.strtochar(),ios::binary);
      if (samplestream.fail())
        flags[0] = 1;
      }


    double betatransform;

    for(i=0;i<nrpar;i++,workbeta++,workbetamean++,workbetas2++,workbetavar++,
        workbetamin++,workbetamax++)

      {

      betatransform = transform * (*workbeta)+addon;


      // storing sampled parameters in binary mode

      if (flags[0] != 1)
        samplestream.write((char *) &betatransform,sizeof betatransform);


      // updating betamean
      if (samplesize==1)
        *workbetamean = betatransform;
      else
        *workbetamean = (1.0/(samplesize))*
                     ((samplesize-1)*(*workbetamean) + betatransform);

      // updating betavar
      *workbetas2 += betatransform*betatransform;
      *workbetavar = (1.0/samplesize)*
                     (*workbetas2)-(*workbetamean)*(*workbetamean);

      // updating betamin, betamax
      if (samplesize==1)
        {
        *workbetamin = betatransform;
        *workbetamax = betatransform;
        }
      else
        {
        if (betatransform < *workbetamin)
          *workbetamin = betatransform;
        if (betatransform > *workbetamax)
          *workbetamax = betatransform;
        }


      // initializing betameanold,betavarold,betaminold and betamaxold
      if (samplesize==1)
        {
        betameanold = betamean;
        betavarold = betavar;
        betaminold = betamin;
        betamaxold = betamax;
        }

      }  // end: for i=0; ...


    } // end: if ( (nriter > optionsp->burnin) && ((nriter-burnin-1) % optionsp->step == 0) )

  if(
        (flags[1]!= 1)
     &&
        (optionsp->get_nriter() > optionsp->get_burnin())
     &&
        ((optionsp->get_nriter()-optionsp->get_burnin()) %
          optionsp->get_nrbetween() == 0)
    )
    {

    optionsp->out("\n");
    optionsp->out("  " + title + "\n");
    optionsp->out("\n");

    if (nrtrials == 0)
      rate = (double(acceptance)/double(optionsp->get_nriter()) )*100;
    else
      rate = (double(acceptance)/double(nrtrials) )*100;
    optionsp->out("  Acceptance rate:    "  + ST::doubletostring(rate,4)
                  + " %\n");
    optionsp->out("\n");


    normold = norm(betameanold);
    if (normold==0)
#if defined(MICROSOFT_VISUAL)
      diffmean = DBL_MAX;
#else
	  diffmean = MAXDOUBLE;
#endif
    else
      diffmean = norm(betamean-betameanold)/normold;

    normold = norm(betavarold);
    if (normold==0)
#if defined(MICROSOFT_VISUAL)
      diffvar = DBL_MAX;
#else
	  diffvar = MAXDOUBLE;
#endif
	else
      diffvar = norm(betavar-betavarold)/normold;

    normold = norm(betaminold);
    if (normold==0)
#if defined(MICROSOFT_VISUAL)
      diffmin = DBL_MAX;
#else
	  diffmin = MAXDOUBLE;
#endif
    else
      diffmin = norm(betamin-betaminold)/normold;

    normold = norm(betamaxold);
    if (normold==0)
#if defined(MICROSOFT_VISUAL)
      diffmax = DBL_MAX;
#else
	  diffmax = MAXDOUBLE;
#endif
    else
      diffmax = norm(betamax-betamaxold)/normold;

    optionsp->out("  Relative Changes in  \n");
    optionsp->out("\n");

    optionsp->out("  Mean:               " + ST::doubletostring(diffmean,6) + "\n");
    optionsp->out("  Variance:           " + ST::doubletostring(diffvar,6) + "\n");
    optionsp->out("  Minimum:            " + ST::doubletostring(diffmin,6) + "\n");
    optionsp->out("  Maximum:            " + ST::doubletostring(diffmax,6) + "\n");

    optionsp->out("\n");
    optionsp->out("\n");

    betameanold.assign(betamean);
    betavarold.assign(betavar);
    betaminold.assign(betamin);
    betamaxold.assign(betamax);

    }  // end: if ((it > burnin) && ((it-burnin) % nrbetween == 0))

  }


void FULLCOND::updatemult(void)
  {
  double diffmean;
  double diffvar;
  double diffmin;
  double diffmax;
  double normold;
  double rate;

  if(
        (optionsp->get_nriter() > optionsp->get_burnin())
     &&
        ((optionsp->get_nriter()-optionsp->get_burnin()-1) % (optionsp->get_step()) == 0)
    )
    {

    register unsigned i;
    double* workbeta = beta.getV();
    double* workbetamean = betamean.getV();
    double* workbetas2 = betas2.getV();
    double* workbetavar = betavar.getV();
    double* workbetamin = betamin.getV();
    double* workbetamax = betamax.getV();

    unsigned samplesize = optionsp->get_samplesize();
    if ((flags[0] != 1) && (samplesize == 1))
      {
      samplestream.open(samplepath.strtochar(),ios::binary);
      if (samplestream.fail())
        flags[0] = 1;
      }


    double betatransform;
    double * worktransformmult = transformmult.getV();

    for(i=0;i<nrpar;i++,workbeta++,workbetamean++,workbetas2++,workbetavar++,
        workbetamin++,workbetamax++,worktransformmult++)

      {

      betatransform = *worktransformmult * (*workbeta);


      // storing sampled parameters in binary mode

      if (flags[0] != 1)
        samplestream.write((char *) &betatransform,sizeof betatransform);


      // updating betamean
      if (samplesize==1)
        *workbetamean = betatransform;
      else
        *workbetamean = (1.0/(samplesize))*
                     ((samplesize-1)*(*workbetamean) + betatransform);

      // updating betavar
      *workbetas2 += betatransform*betatransform;
      *workbetavar = (1.0/samplesize)*
                     (*workbetas2)-(*workbetamean)*(*workbetamean);

      // updating betamin, betamax
      if (samplesize==1)
        {
        *workbetamin = betatransform;
        *workbetamax = betatransform;
        }
      else
        {
        if (betatransform < *workbetamin)
          *workbetamin = betatransform;
        if (betatransform > *workbetamax)
          *workbetamax = betatransform;
        }


      // initializing betameanold,betavarold,betaminold and betamaxold
      if (samplesize==1)
        {
        betameanold = betamean;
        betavarold = betavar;
        betaminold = betamin;
        betamaxold = betamax;
        }

      }  // end: for i=0; ...


    } // end: if ( (nriter > optionsp->burnin) && ((nriter-burnin-1) % optionsp->step == 0) )

  if(
        (flags[1]!= 1)
     &&
        (optionsp->get_nriter() > optionsp->get_burnin())
     &&
        ((optionsp->get_nriter()-optionsp->get_burnin()) %
          optionsp->get_nrbetween() == 0)
    )
    {

    optionsp->out("\n");
    optionsp->out("  " + title + "\n");
    optionsp->out("\n");

    if (nrtrials == 0)
      rate = (double(acceptance)/double(optionsp->get_nriter()) )*100;
    else
      rate = (double(acceptance)/double(nrtrials) )*100;
    optionsp->out("  Acceptance rate:    "  + ST::doubletostring(rate,4)
                  + " %\n");
    optionsp->out("\n");


    normold = norm(betameanold);
    if (normold==0)
#if defined(MICROSOFT_VISUAL)
      diffmean = DBL_MAX;
#else
	  diffmean = MAXDOUBLE;
#endif
    else
      diffmean = norm(betamean-betameanold)/normold;

    normold = norm(betavarold);
    if (normold==0)
#if defined(MICROSOFT_VISUAL)
      diffvar = DBL_MAX;
#else
	  diffvar = MAXDOUBLE;
#endif
	else
      diffvar = norm(betavar-betavarold)/normold;

    normold = norm(betaminold);
    if (normold==0)
#if defined(MICROSOFT_VISUAL)
      diffmin = DBL_MAX;
#else
	  diffmin = MAXDOUBLE;
#endif
    else
      diffmin = norm(betamin-betaminold)/normold;

    normold = norm(betamaxold);
    if (normold==0)
#if defined(MICROSOFT_VISUAL)
      diffmax = DBL_MAX;
#else
	  diffmax = MAXDOUBLE;
#endif
    else
      diffmax = norm(betamax-betamaxold)/normold;

    optionsp->out("  Relative Changes in  \n");
    optionsp->out("\n");

    optionsp->out("  Mean:               " + ST::doubletostring(diffmean,6) + "\n");
    optionsp->out("  Variance:           " + ST::doubletostring(diffvar,6) + "\n");
    optionsp->out("  Minimum:            " + ST::doubletostring(diffmin,6) + "\n");
    optionsp->out("  Maximum:            " + ST::doubletostring(diffmax,6) + "\n");

    optionsp->out("\n");
    optionsp->out("\n");

    betameanold.assign(betamean);
    betavarold.assign(betavar);
    betaminold.assign(betamin);
    betamaxold.assign(betamax);

    }  // end: if ((it > burnin) && ((it-burnin) % nrbetween == 0))

  }



bool FULLCOND::posteriormode(void)
  {

  double diffmean;
  double normold;

  normold = norm(betameanold);

  if (normold==0)
    #if defined(MICROSOFT_VISUAL)
      diffmean = DBL_MAX;
    #else
      diffmean = MAXDOUBLE;
    #endif
  else
    diffmean = norm(beta-betameanold)/normold;

  betameanold.assign(beta);

  unsigned i;
  double* workbetamean = betamean.getV();
  double* workbeta = beta.getV();

  for(i=0;i<nrpar;i++,workbetamean++,workbeta++)
    {
    *workbetamean = transform * (*workbeta)+addon;
    }

  if (diffmean <= 0.00001)
    {
    return true;
    }
  else
    return false;

  }


void FULLCOND::posteriormode_set_beta_mode(void)
  {
  beta_mode.assign(beta);
  }


bool FULLCOND::posteriormode_converged(const unsigned & itnr)
  {
  return true;
  }


void FULLCOND::outresults(void)
  {

  double rate;

  if (flags[2] != 1)
    {
    optionsp->out("\n");
    optionsp->out("  " + title + "\n",true);
    optionsp->out("\n");
    optionsp->out("\n");
    }

  // computing acceptance rate
  if (optionsp->get_samplesize() > 0)
    {

    if (flags[2] != 1)
      {
      if (nrtrials == 0)
        rate = (double(acceptance)/double(optionsp->get_nriter()) )*100;
      else
        rate = (double(acceptance)/double(nrtrials))*100;
      optionsp->out("  Acceptance rate:    "  + ST::doubletostring(rate,4) + " %\n");
      optionsp->out("\n");
      }

    if (( (transformnonlinear== true) || (transformed==true) ) && (flags[0] != 1))
      {

      transformed = true;

      }
    else if (flags[0] != 1)
      {

      samplestream.close();
      datamatrix sample(optionsp->get_samplesize(),1);
      unsigned i;

      double* wqu1l = betaqu_l1_lower.getV();
      double* wqu2l = betaqu_l2_lower.getV();
      double* wqu50 = betaqu50.getV();
      double* wqu1u = betaqu_l2_upper.getV();
      double* wqu2u = betaqu_l1_upper.getV();

      for(i=0;i<nrpar;i++,wqu1l++,wqu2l++,wqu50++,wqu1u++,wqu2u++)
        {
        readsample(sample,i);
        *wqu1l = sample.quantile(lower1,0);
        *wqu2l = sample.quantile(lower2,0);
        *wqu50 = sample.quantile(50,0);
        *wqu1u = sample.quantile(upper1,0);
        *wqu2u = sample.quantile(upper2,0);
        }

      } // end: if (flags[0] != 1)

    if (flags[0] == 1)
      {
      optionsp->out("  NOTE: Sampled parameters have not been stored, i.e. posterior quantiles are not available!\n");
      optionsp->out("\n");
      }

    }  // if (optionsp->get_samplesize() > 0)
  }

// ---------------------------------------------------------------------------
// ------------------------FOR STEPWISE---------------------------------------
// ---------------------------------------------------------------------------

  // FUNCTION: compute_lambdavec
  // TASK: returns the values for the smoothing parameter (logarithmic scale)

void FULLCOND::compute_lambdavec(vector<double> & lvec,unsigned & number)
  {
  if(get_lambdamax_opt()==true)
    lambdamax = lambda_from_df(get_df_lambdamax(),lambdamax);
  if(get_lambdamin_opt()==true)
    lambdamin = lambda_from_df(get_df_lambdamin(),lambdamin);

  double l = log10(lambdamin);
  double u = log10(lambdamax);
  if(number==1)
    lvec.push_back(lambdamax);
  else
    {
    unsigned j;
    for(j=0;j<number;j++)
      lvec.push_back(pow(10,l+double(j)*((u-l)/(double(number)-1))));
    }
  }


  // FUNCTION: compute_lambdavec
  // TASK: returns the values for the smoothing parameter (the resulting df's are equidistant)

void FULLCOND::compute_lambdavec_equi(vector<double> & lvec,unsigned & number)
  {
  double diff = (get_df_lambdamin()-get_df_lambdamax())/(number-1);
  lambdamax = lambda_from_df(get_df_lambdamax(),lambdamax);
  lambdamin = lambda_from_df(get_df_lambdamin(),lambdamin);
  lvec.push_back(lambdamin);
  int i;
  for(i=number-2;i>=1;i--)
     lvec.push_back(lambda_from_df(get_df_lambdamax() + i*diff,
                                     lambdamin + i*(lambdamax-lambdamin)/(number-1)));
  lvec.push_back(lambdamax);
  }

  // FUNCTION: compute_lambdavec
  // TASK: returns the value for the smoothing parameter of a given df

double FULLCOND::lambda_from_df(double & df_wunsch, double & lambda_vorg)
  {
  update_stepwise(lambda_vorg);
  double df_vorg = compute_df();
  double lambda_unten, lambda_oben, df_mitte;
  if( (df_wunsch-df_vorg) < get_accuracy() && (df_wunsch-df_vorg) > -1*get_accuracy() )
     return lambda_vorg;
  else if((df_wunsch-df_vorg) >= get_accuracy())
     {
     double lambda_vers = lambda_vorg;
     double df_vers = df_vorg;
     while(df_vers < df_wunsch)
        {
        lambda_vers = lambda_vers*0.75;
        update_stepwise(lambda_vers);
        df_vers = compute_df();
        if( (df_wunsch-df_vers) < get_accuracy() && (df_wunsch-df_vers) > -1*get_accuracy() )
           return lambda_vers;
        }
     lambda_unten = lambda_vorg;
     lambda_oben = lambda_vers;
     df_mitte = df_vers;
     }
  else
     {
     double lambda_vers = lambda_vorg;
     double df_vers = df_vorg;
     while(df_vers > df_wunsch)
        {
        lambda_vers = lambda_vers*2;
        update_stepwise(lambda_vers);
        df_vers = compute_df();
        if( (df_wunsch-df_vers) < get_accuracy() && (df_wunsch-df_vers) > -1*get_accuracy())
          return lambda_vers;
        }
     lambda_unten = lambda_vers;
     lambda_oben = lambda_vorg;
     df_mitte = df_vers;
     }
  double lambda_mitte;
  while( (df_mitte-df_wunsch) >= get_accuracy() || (df_mitte-df_wunsch) <= -1*get_accuracy() )
     {
     lambda_mitte = lambda_oben + (lambda_unten - lambda_oben) / 2;
     update_stepwise(lambda_mitte);
     df_mitte = compute_df();
     if(df_mitte < df_wunsch)
       lambda_unten = lambda_mitte;
     else
       lambda_oben = lambda_mitte;
     }
  return lambda_mitte;
  }

// ---------------------------------------------------------------------------
// ------------------------END: FOR STEPWISE----------------------------------
// ---------------------------------------------------------------------------


} // end: namespace MCMC



