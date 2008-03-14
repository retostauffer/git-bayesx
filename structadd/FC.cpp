
#include "FC.h"
#include "clstring.h"

//------------------------------------------------------------------------------
//----------------- CLASS: FC implementation of member functions ---------------
//------------------------------------------------------------------------------


namespace MCMC
{


FC::FC(void)
  {
  }


FC::FC(GENERAL_OPTIONS * o,const ST::string & t,const unsigned & rows,
                   const unsigned & cols,const ST::string & fp)
  {

  changingdesign = false;

  optionsp = o;

  title = t;

  samplepath = fp;

  setbeta(rows,cols,0);

  acceptance = 0;
  nrtrials = 0;

  column = 0;

  transform = 1;

  addon = 0;

  }


FC::FC(const FC & m)
  {

  changingdesign = m.changingdesign;
  optionsp = m.optionsp;

  title = m.title;

  samplepath = m.samplepath;

  priorassumptions = m.priorassumptions;

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

  betameanold = m.betameanold;
  betavarold = m.betavarold;
  betaminold = m.betaminold;
  betamaxold = m.betamaxold;

  transform = m.transform;
  addon = m.addon;


  acceptance = m.acceptance;
  nrtrials = m.nrtrials;

  column = m.column;

  }



const FC & FC::operator=(const FC & m)
  {

  if (this==&m)
	 return *this;

  changingdesign = m.changingdesign;

  optionsp = m.optionsp;

  title = m.title;

  samplepath = m.samplepath;

  priorassumptions = m.priorassumptions;


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

  betameanold = m.betameanold;
  betavarold = m.betavarold;
  betaminold = m.betaminold;
  betamaxold = m.betamaxold;

  transform = m.transform;
  addon = m.addon;


  acceptance = m.acceptance;
  nrtrials = m.nrtrials;

  column = m.column;

  return *this;

  }



void FC::setbeta(const unsigned & rows,const unsigned & cols,
                       const double & v)
  {

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


void FC::setbetavalue(const unsigned & row,const unsigned & col,
                            const double & v)
  {
  beta(row,col) = v;
  }


void FC::setbeta(const datamatrix & betanew)
  {
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



void FC::reset(void)
  {

  setbeta(beta.rows(),beta.cols(),0);
  acceptance = 0;
  nrtrials = 0;
  }



void FC::readsample(datamatrix & sample,const unsigned & nr,
                          const unsigned & col) const
  {

  unsigned nrpar = beta.cols()*beta.rows();
  unsigned size = sizeof(double);
  ifstream in;
  in.open(samplepath.strtochar(),ios::binary);
  unsigned i;
  double* work=sample.getV()+col;
  unsigned s = sample.cols();
  in.seekg(size*nr);
  for (i=0;i<optionsp->samplesize;i++,work+=s)
    {
    in.read((char*) work,size);
    in.seekg(size*(nrpar-1),ios::cur);
    }
  }


void FC::readsample2(datamatrix & b,const unsigned & nr) const
  {

  unsigned nrpar = beta.cols()*beta.rows();
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


void FC::readsample3(datamatrix & b) const
  {

  unsigned size = sizeof(double);
  ifstream in;
  in.open(samplepath.strtochar(),ios::binary);
  double * work = b.getV();
  unsigned i,j;
  for(i=0;i<b.rows();i++)
    for(j=0;j<b.cols();j++,work++)
      in.read((char *) work,size);

  }




datamatrix FC::compute_autocorr(const unsigned & lag,const unsigned & row,
                                      const unsigned & col) const
  {
  datamatrix sample(optionsp->samplesize,1);
  unsigned nr = row*(beta.cols())+col;
  readsample(sample,nr);
  return sample.autocorr (1,lag,0);
  }


void FC::get_samples(const ST::string & filename,const unsigned & step) const
  {
  ofstream out(filename.strtochar());
  assert(!out.fail());

  unsigned i,j,k;
  unsigned nrpar = beta.cols()*beta.rows();

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

  for(i=0;i<optionsp->samplesize;i++)
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



void FC::update(void)
  {
  double diffmean;
  double diffvar;
  double diffmin;
  double diffmax;
  double normold;
  double rate;

  if(
     (optionsp->nriter > optionsp->burnin)
     &&
     ((optionsp->nriter-optionsp->burnin-1) % (optionsp->step) == 0)
    )
    {

    register unsigned i;
    double* workbeta = beta.getV();
    double* workbetamean = betamean.getV();
    double* workbetas2 = betas2.getV();
    double* workbetavar = betavar.getV();
    double* workbetamin = betamin.getV();
    double* workbetamax = betamax.getV();

    unsigned samplesize = optionsp->samplesize;

    samplestream.open(samplepath.strtochar(),ios::binary);

    double betatransform;

    unsigned nrpar=beta.cols()*beta.rows();


    for(i=0;i<nrpar;i++,workbeta++,workbetamean++,workbetas2++,workbetavar++,
        workbetamin++,workbetamax++)

      {

      betatransform = transform * (*workbeta)+addon;


      // storing sampled parameters in binary mode


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

      }  // end: for i=0; ...


    } // end: if ( (nriter > optionsp->burnin) && ((nriter-burnin-1) % optionsp->step == 0) )

  if(
        (optionsp->nriter > optionsp->burnin)
     &&
        ((optionsp->nriter-optionsp->burnin) %
          optionsp->nrbetween == 0)
    )
    {

    optionsp->out("\n");
    optionsp->out("  " + title + "\n");
    optionsp->out("\n");

    if (nrtrials == 0)
      rate = (double(acceptance)/double(optionsp->nriter) )*100;
    else
      rate = (double(acceptance)/double(nrtrials) )*100;
    optionsp->out("  Acceptance rate:    "  + ST::doubletostring(rate,4)
                  + " %\n");
    optionsp->out("\n");


    normold = norm(betameanold);
    if (normold==0)
	  diffmean = MAXDOUBLE;
    else
      diffmean = norm(betamean-betameanold)/normold;

    normold = norm(betavarold);
    if (normold==0)
	  diffvar = MAXDOUBLE;
	else
      diffvar = norm(betavar-betavarold)/normold;

    normold = norm(betaminold);
    if (normold==0)
	  diffmin = MAXDOUBLE;
    else
      diffmin = norm(betamin-betaminold)/normold;

    normold = norm(betamaxold);
    if (normold==0)
	  diffmax = MAXDOUBLE;
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


bool FC::posteriormode(void)
  {

  unsigned nrpar= beta.cols()*beta.rows();

  double diffmean;
  double normold;

  normold = norm(betameanold);

  if (normold==0)
      diffmean = MAXDOUBLE;
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



void FC::outresults(void)
  {

  unsigned nrpar=beta.rows()*beta.cols();

  double rate;

  optionsp->out("\n");
  optionsp->out("  " + title + "\n",true);
  optionsp->out("\n");
  optionsp->out("\n");


  // computing acceptance rate
  if (optionsp->samplesize > 0)
    {


    if (nrtrials == 0)
      rate = (double(acceptance)/double(optionsp->nriter))*100;
    else
      rate = (double(acceptance)/double(nrtrials))*100;
    optionsp->out("  Acceptance rate:    "  + ST::doubletostring(rate,4) + " %\n");
    optionsp->out("\n");


    samplestream.close();
    datamatrix sample(optionsp->samplesize,1);
    unsigned i;

    double* wqu1l = betaqu_l1_lower.getV();
    double* wqu2l = betaqu_l2_lower.getV();
    double* wqu50 = betaqu50.getV();
    double* wqu1u = betaqu_l2_upper.getV();
    double* wqu2u = betaqu_l1_upper.getV();

    for(i=0;i<nrpar;i++,wqu1l++,wqu2l++,wqu50++,wqu1u++,wqu2u++)
      {
      readsample(sample,i);
      *wqu1l = sample.quantile(optionsp->lower1,0);
      *wqu2l = sample.quantile(optionsp->lower2,0);
      *wqu50 = sample.quantile(50,0);
      *wqu1u = sample.quantile(optionsp->upper1,0);
      *wqu2u = sample.quantile(optionsp->upper2,0);
      }


    }  // if (optionsp->get_samplesize() > 0)
  }



} // end: namespace MCMC



