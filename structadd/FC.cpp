
#include "FC.h"
#include "clstring.h"

using std::ofstream;
using std::ifstream;
using std::ios;

//------------------------------------------------------------------------------
//----------------- CLASS: FC implementation of member functions ---------------
//------------------------------------------------------------------------------


namespace MCMC
{

void FC::read_options(vector<ST::string> & op,vector<ST::string> & vn)
  {

  }

FC::FC(void)
  {
  }


FC::FC(GENERAL_OPTIONS * o,const ST::string & t,const unsigned & rows,
                   const unsigned & cols,const ST::string & fp,bool sstore)
  {

  samplestore = sstore;
  nosamples = false;
  optionsp = o;

  title = t;

  samplepath = fp;

  setbeta(rows,cols,0);

  acceptance = 0;
  nrtrials = 0;

  column = 0;

  addon = 0;

  meaneffect = 0;

  }


FC::FC(const FC & m)
  {

  samplestore = m.samplestore;
  nosamples = m.nosamples;

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

  sampled_beta = m.sampled_beta;

  transform = m.transform;
  addon = m.addon;


  acceptance = m.acceptance;
  nrtrials = m.nrtrials;

  column = m.column;

  meaneffect = m.meaneffect;

  }



const FC & FC::operator=(const FC & m)
  {

  if (this==&m)
	 return *this;

  samplestore = m.samplestore;
  nosamples = m.nosamples;

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

  sampled_beta = m.sampled_beta;

  transform = m.transform;
  addon = m.addon;


  acceptance = m.acceptance;
  nrtrials = m.nrtrials;

  column = m.column;

  meaneffect = m.meaneffect;

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
  transform = datamatrix(cols,1,1);

  }


void FC::setbetavalue(const unsigned & row,const unsigned & col,
                            const double & v)
  {
  beta(row,col) = v;
  }


void FC::setbeta(const datamatrix & betanew)
  {
  beta = betanew;
  betamean = datamatrix(beta.rows(),beta.cols(),0);
  beta_mode = datamatrix(beta.rows(),beta.cols(),0);
  betas2 = datamatrix(beta.rows(),beta.cols(),0);
  betavar = datamatrix(beta.rows(),beta.cols(),0);
  betamin = datamatrix(beta.rows(),beta.cols(),0);
  betamax = datamatrix(beta.rows(),beta.cols(),0);
  betaqu_l1_lower = datamatrix(beta.rows(),beta.cols(),0);
  betaqu_l2_lower = datamatrix(beta.rows(),beta.cols(),0);
  betaqu_l1_upper = datamatrix(beta.rows(),beta.cols(),0);
  betaqu_l2_upper = datamatrix(beta.rows(),beta.cols(),0);
  betaqu50 = datamatrix(beta.rows(),beta.cols(),0);
  betameanold = datamatrix(beta.rows(),beta.cols(),0);
  betavarold = datamatrix(beta.rows(),beta.cols(),0);
  betaminold = datamatrix(beta.rows(),beta.cols(),0);
  betamaxold = datamatrix(beta.rows(),beta.cols(),0);
  transform = datamatrix(beta.cols(),1,1);

  }


void FC::readsample(datamatrix & sample,const unsigned & nr,
                          const unsigned & col) const
  {
  unsigned nrpar = beta.cols()*beta.rows();


  unsigned i;

  unsigned s = sample.cols();

  double* work=sample.getV()+col;

  if (samplestore==true)
    {
    unsigned size = sizeof(double);

    ifstream in;
    in.open(samplepath.strtochar(),ios::binary);

    in.seekg(size*nr);
    for (i=0;i<optionsp->samplesize;i++,work+=s)
      {
      in.read((char*) work,size);
      in.seekg(size*(nrpar-1),ios::cur);
      }
    }
  else
    {

    double * sampled_betap = sampled_beta.getV()+nr;
    for (i=0;i<optionsp->samplesize;i++,work+=s,sampled_betap+=nrpar)
      {
      *work = *sampled_betap;
      }

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

  if (nosamples == false)
    {
    if (samplestore)
      {
      unsigned size = sizeof(double);
      ifstream in;
      in.open(samplepath.strtochar(),ios::binary);
      double * work = b.getV();
      unsigned i,j;
      for(i=0;i<b.rows();i++)
        for(j=0;j<b.cols();j++,work++)
          {
          in.read((char *) work,size);
          }
      }
    else
      {
      b.assign(sampled_beta);
      }

    } // end: if (nosamples == false)


  }


datamatrix FC::compute_autocorr(const unsigned & lag,const unsigned & row,
                                      const unsigned & col) const
  {
  if (nosamples==false)
    {

    unsigned nr = row*(beta.cols())+col;

    if (samplestore)
      {
      unsigned o = optionsp->samplesize;
      datamatrix sample(o,1);
      readsample(sample,nr);
      return sample.autocorr (1,lag,0);
      }
    else
      {
      return sampled_beta.autocorr(1,lag,nr);
      }
    }
  else
    return datamatrix(1,1);
  }


double FC::compute_autocorr_single(const unsigned & lag,const unsigned & row,
                                      const unsigned & col) const
  {
  if (nosamples == false)
    {

    unsigned nr = row*(beta.cols())+col;

    if (samplestore==true)
      {
      unsigned o = optionsp->samplesize;
      datamatrix sample(o,1);
      readsample(sample,nr);
      return sample.autocorr (lag,0);
      }
    else
      {
      return sampled_beta.autocorr(lag,nr);
      }
    }
  else
    return 0;
  }


void FC::compute_autocorr(const ST::string & path, unsigned lag) const
  {
  if (nosamples == false)
    {
    unsigned k,i,c,l,r;
    unsigned nrpar;
    bool misstot;
    double min,max,mean,autoc;

    ofstream out(path.strtochar());

    optionsp->out(path);
    optionsp->out("\n");

    out << "lag  ";

    for (k=0;k<beta.cols();k++)
      for(i=0;i<beta.rows();i++)
        {
        if (beta.cols() == 1)
          out << "b_" << (i+1) << " ";
        else
          out << "b_" << (i+1) << "_" << (k+1) << " ";
        }

    out  << "b_min " << "b_mean " << "b_max " << endl;

    misstot = false;

    for(l=1;l<=lag;l++)
      {

      nrpar = 0;

      out << l << "  ";

      min = 1;
      max = -1;
      mean = 0;

      for(c=0;c<beta.cols();c++)
        for (r=0;r<beta.rows();r++)
          {
          autoc = compute_autocorr_single(l,r,c);
          if ( (autoc <= 1) && (autoc >= -1) )
            {
            nrpar++;
            if (autoc < min)
              min = autoc;
            if (autoc > max)
              max = autoc;
            mean += autoc;
            out << autoc << "  ";
            }
          else
            {
            out << "NA  " << endl;
            misstot = true;
            }
          }

        out << min << "  ";
        out << max << "  ";
        out << mean/nrpar << "  ";
        out << endl;
      }  // end: for(l=0;l<lag;l++)

    if (misstot==true)
      {
      optionsp->out("WARNING: There were undefined autocorrelations\n",true,true);
      optionsp->out("\n");
      }
    }
  }


void FC::compute_autocorr_all(const ST::string & path, unsigned lag) const
  {
  compute_autocorr(path,lag);
  }


void FC::get_samples(const ST::string & filename) const
  {
  if (nosamples == false)
    {

    unsigned i,j,k;

    unsigned nrpar = beta.rows()*beta.cols();

    ofstream out(filename.strtochar());
    assert(!out.fail());

    out << "intnr " << " ";
    if (beta.cols() > 1)
      {
      for (j=0;j<beta.rows();j++)
        for(k=0;k<beta.cols();k++)
          out << "b_" << (j+1) << "_" << (k+1) << " ";
      }
    else
      {
      for (j=0;j<nrpar;j++)
        out << "b_" << (j+1) << " ";
      }

    out << endl;


    if (samplestore)
      {

      datamatrix betac(beta.rows(),beta.cols());


      for(i=0;i<optionsp->samplesize;i++)
        {
        readsample2(betac,i);
        out << (i+1) << " ";
        for (j=0;j<beta.rows();j++)
          for(k=0;k<betac.cols();k++)
            out << betac(j,k) << " ";

        out << endl;
        }


      }
    else
      {

      double * sampled_betap = sampled_beta.getV();
      for(i=0;i<optionsp->samplesize;i++)
        {
        out << (i+1) << " ";
        for (j=0;j<nrpar;j++,sampled_betap++)
          out << *sampled_betap << " ";

        out << endl;
        }

      }

    out.close();

    } // end: if (nosamples == false)
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

    register unsigned i,j;
    double* workbeta = beta.getV();
    double* workbetamean = betamean.getV();
    double* workbetas2 = betas2.getV();
    double* workbetavar = betavar.getV();
    double* workbetamin = betamin.getV();
    double* workbetamax = betamax.getV();


    unsigned samplesize = optionsp->samplesize;

    if (samplesize==1)
     {
     if (samplestore==true)
       samplestream.open(samplepath.strtochar(),ios::binary);
     else
       {
       unsigned ssize = optionsp->compute_samplesize();
       unsigned npar = beta.rows()*beta.cols();
       sampled_beta = datamatrix(ssize,npar,0);
       }
     }


    double betatransform;

    double * sbetap = sampled_beta.getV() +  (samplesize-1)*sampled_beta.cols();

    for(i=0;i<beta.rows();i++)
      {
      for (j=0;j<beta.cols();j++,workbeta++,workbetamean++,workbetas2++,
           workbetavar++,workbetamin++,workbetamax++,sbetap++)

        {

        betatransform = transform(j,0) * (*workbeta)+addon;


        // storing sampled parameters in binary mode

        if (samplestore==true)
          samplestream.write((char *) &betatransform,sizeof betatransform);
        else
          *sbetap = betatransform;

        // updating betamean
        if (samplesize==1)
          *workbetamean = betatransform;
        else
          *workbetamean = (1.0/(samplesize))*
                       ((samplesize-1)*(*workbetamean) + betatransform);

        // updating betavar
        if (samplesize==1)
          *workbetas2 = betatransform*betatransform;
        else
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

        }

      }  // end: for i=0; ...


    } // end: if ( (nriter > optionsp->burnin) && ((nriter-burnin-1) % optionsp->step == 0) )

  if(
        (optionsp->nriter > optionsp->burnin)
     &&
        ((optionsp->nriter-optionsp->burnin) %
          optionsp->nrbetween == 0)
     && (title!="")
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

  if ( (optionsp->nriter==optionsp->iterations) && (samplestore==true))
    {
    samplestream.close();
    }

  }


bool FC::posteriormode(void)
  {

  double diffmean;
  double normold;

  normold = norm(betameanold);

  if (normold==0)
      diffmean = MAXDOUBLE;
  else
    diffmean = norm(beta-betameanold)/normold;

  betameanold.assign(beta);

  posteriormode_betamean();

  if (diffmean <= 0.00001)
    {
    return true;
    }
  else
    return false;

  }


void FC::posteriormode_betamean(void)
  {

  unsigned i,j;
  double* workbetamean = betamean.getV();
  double* workbeta = beta.getV();

  for(i=0;i<beta.rows();i++)
    for (j=0;j<beta.cols();j++,workbetamean++,workbeta++)
      {
      *workbetamean = transform(j,0) * (*workbeta)+addon;
      }
  }


void FC::outresults_acceptance(void)
  {
  if (optionsp->samplesize > 0)
    {
    double rate;
    if (nrtrials == 0)
      rate = (double(acceptance)/double(optionsp->nriter))*100;
    else
      rate = (double(acceptance)/double(nrtrials))*100;
    optionsp->out("    Acceptance rate:    "  + ST::doubletostring(rate,4) + " %\n");
    optionsp->out("\n");
    }
  }


void FC::outresults(ofstream & out_stata, ofstream & out_R,
                    const ST::string & pathresults)
  {

  unsigned nrpar=beta.rows()*beta.cols();

  if (title != "")
    {
    optionsp->out("\n");
    optionsp->out("  " + title + "\n",true);
    optionsp->out("\n");
    }

  if (optionsp->samplesize > 0)
    {

    if (samplestore==true)
      {
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
      }
    else
      {

      unsigned i;


      statmatrix<int> index(sampled_beta.rows(),1);

      double* wqu1l = betaqu_l1_lower.getV();
      double* wqu2l = betaqu_l2_lower.getV();
      double* wqu50 = betaqu50.getV();
      double* wqu1u = betaqu_l2_upper.getV();
      double* wqu2u = betaqu_l1_upper.getV();

      for(i=0;i<nrpar;i++,wqu1l++,wqu2l++,wqu50++,wqu1u++,wqu2u++)
        {

        index.indexinit();

        sampled_beta.indexsort(index,0,sampled_beta.rows()-1,i,0);


        *wqu1l = sampled_beta.quantile(optionsp->lower1,i,index);
        *wqu2l = sampled_beta.quantile(optionsp->lower2,i,index);
        *wqu50 = sampled_beta.quantile(50,i,index);
        *wqu1u = sampled_beta.quantile(optionsp->upper1,i,index);
        *wqu2u = sampled_beta.quantile(optionsp->upper2,i,index);
        }

      }

    }  // if (optionsp->get_samplesize() > 0)
  }


void FC::reset(void)
  {

  setbeta(beta.rows(),beta.cols(),0);
  acceptance = 0;
  nrtrials = 0;
  if (nosamples==false && samplestore == true)
    {
    samplestream.close();
    remove(samplepath.strtochar());
    }

  }


} // end: namespace MCMC



