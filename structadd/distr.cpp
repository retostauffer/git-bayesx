
#include "distr.h"

namespace MCMC
{




DISTR::DISTR(GENERAL_OPTIONS * o, const datamatrix & r,
             const datamatrix & w)
  {


  optionsp = o;

  family = "unknown";

  response = r;
  responsename = "Y";

  nrobs = response.rows();

  if (w.rows() == 1)
    {
    weight = datamatrix(r.rows(),1,1);
    }
  else
    {
    weight = w;
    }

  workingweight = weight;
  changingweight = false;

  weightname = "W";

  linearpred = datamatrix(nrobs,1,0);
  linearpredprop = datamatrix(nrobs,1,0);

  linpred_current = &linearpred;
  linpred_proposed = &linearpredprop;

  trmult=1;

  }



DISTR::DISTR(const DISTR & d)
  {

  optionsp = d.optionsp;
  nrobs = d.nrobs;

  response = d.response;
  responsename = d.responsename;

  weight = d.weight;
  weightname = d.weightname;

  workingweight = d.workingweight;
  changingweight = d.changingweight;

  linearpred = d.linearpred;
  linearpredprop = d.linearpredprop;
  linpred_current = d.linpred_current;
  linpred_proposed = d.linpred_proposed;

  family = d.family;

  trmult=d.trmult;
  }


const DISTR & DISTR::operator=(const DISTR & d)
  {
  if (this == &d)
    return *this;

  optionsp = d.optionsp;
  nrobs = d.nrobs;

  response = d.response;
  responsename = d.responsename;

  weight = d.weight;
  weightname = d.weightname;

  workingweight = d.workingweight;
  changingweight = d.changingweight;

  linearpred = d.linearpred;
  linearpredprop = d.linearpredprop;
  linpred_current = d.linpred_current;
  linpred_proposed = d.linpred_proposed;

  family = d.family;

  trmult=d.trmult;  

  return *this;
  }



void DISTR::outoptions(void)
  {
  optionsp->out("RESPONSE DISTRIBUTION:\n",true);
  optionsp->out("\n");
  optionsp->out("  Family: " + family + "\n");
  optionsp->out("  Number of observations: " + ST::inttostring(nrobs) + "\n");
  }



double DISTR::loglikelihood(const bool & current) const
  {

  register unsigned  i;
  double* workweight = weight.getV();
  double* workres = response.getV();
  double help = 0;

  double* worklin;
  if (current)
    worklin = (*linpred_current).getV();
  else
    worklin = (*linpred_proposed).getV();

  for (i=0;i<nrobs;i++,workweight++,worklin++,workres++)
    help += loglikelihood(workres,worklin,workweight);

  return help;

  }


double DISTR::loglikelihood(const unsigned & beg,const unsigned & end,
                                   const statmatrix<int> & index,
                                   const bool & current)
  {

  unsigned i;
  double help = 0;

  int* workind = index.getV()+beg;

  if (current)
    {
    for (i=beg;i<=end;i++,workind++)
      help+=loglikelihood(&response(*workind,0),&((*linpred_current)(*workind,0)),
                          &weight(*workind,0));

    }
  else
    {
    for (i=beg;i<=end;i++,workind++)
      help+=loglikelihood(&response(*workind,0),&((*linpred_proposed)(*workind,0)),
                          &weight(*workind,0));
    }

  return help;

  }


void DISTR::swap_linearpred(void)
  {
  datamatrix * help = linpred_current;
  linpred_current = linpred_proposed;
  linpred_proposed = help;
  }


void DISTR::update(void)
  {
  } // end: update



void DISTR::outresults(ST::string pathresults)
  {
  optionsp->out("\n");
  }


void DISTR::reset(void)
  {
  linearpred = datamatrix(nrobs,linearpred.cols(),0);
  linearpredprop = linearpred;
  linpred_current = &linearpred;
  linpred_proposed = &linearpredprop;
  }


double DISTR::get_scale(void)
  {
  return 1;
  }

//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_gaussian --------------------------
//------------------------------------------------------------------------------


DISTR_gaussian::DISTR_gaussian(const double & a,
                                             const double & b,
                                             GENERAL_OPTIONS * o,
                                             const datamatrix & r,
                                             const ST::string & ps,
                                             const datamatrix & w)
  : DISTR(o,r,w)

  {

  a_invgamma = a;
  b_invgamma = b;
  family = "Gaussian";

  standardise();

  FCsigma2 = FC(o,"Gaussian variance parameter",1,1,ps);
  FCsigma2.transform(0,0) = pow(trmult,2);

  sigma2 = 1;

  }


const DISTR_gaussian & DISTR_gaussian::operator=(
                                      const DISTR_gaussian & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  a_invgamma = nd.a_invgamma;
  b_invgamma = nd.b_invgamma;
  FCsigma2 = nd.FCsigma2;
  sigma2=nd.sigma2;
  return *this;
  }



DISTR_gaussian::DISTR_gaussian(const DISTR_gaussian & nd)
   : DISTR(DISTR(nd))
  {
  a_invgamma = nd.a_invgamma;
  b_invgamma = nd.b_invgamma;
  FCsigma2 = nd.FCsigma2;
  sigma2 = nd.sigma2;
  }



void DISTR_gaussian::standardise(void)
  {

  trmult = sqrt(response.var(0,weight));

  unsigned i;
  double * workresp = response.getV();
  double * worklin = (*linpred_current).getV();
  for (i=0;i<nrobs;i++,workresp++,worklin++)
   {
   *workresp = *workresp/trmult;
   *worklin = *worklin/trmult;
   }

  }

void DISTR_gaussian::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: identity\n");

  optionsp->out("  Hyperparameter a: " + ST::doubletostring(a_invgamma,6) + "\n");
  optionsp->out("  Hyperparameter b: " + ST::doubletostring(b_invgamma,6) + "\n");

  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussian::update(void)
  {

  register unsigned i;

  double help;

  double * worklin;
  double * workresp;
  double * workweight;


  // scaleparameter

  double sum = 0;

  worklin = (*linpred_current).getV();
  workresp = response.getV();
  workweight = weight.getV();

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {
    help = *workresp - *worklin;
    sum += *workweight*pow(help,2);
    }

  sigma2  = rand_invgamma(a_invgamma+0.5*nrobs,
                          b_invgamma+0.5*sum);

  FCsigma2.beta(0,0) = sigma2;
  FCsigma2.update();

  DISTR::update();

  }



double DISTR_gaussian::loglikelihood(double * res, double * lin,
                                     double * w) const
  {
  double help = *res-*lin;
  return  - *w * (pow(help,2))/(2* sigma2);
  }


double DISTR_gaussian::get_scale(void)
  {
  return sigma2;
  }



bool DISTR_gaussian::posteriormode(void)
  {

  unsigned i;

  double * worklin = (*linpred_current).getV();
  double * workresp = response.getV();
  double * workweight = weight.getV();

  double sum = 0;
  double help;

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {
    help = *workresp - *worklin;
    sum += *workweight*pow(help,2);
    }

  sigma2 = (1.0/nrobs)*sum;

  return true;

  }


void DISTR_gaussian::outresults(ST::string pathresults)
  {
  DISTR::outresults();

  FCsigma2.outresults("");


  ST::string l1 = ST::doubletostring(optionsp->lower1,4);
  ST::string l2 = ST::doubletostring(optionsp->lower2,4);
  ST::string u1 = ST::doubletostring(optionsp->upper1,4);
  ST::string u2 = ST::doubletostring(optionsp->upper2,4);

  ST::string nl1 = ST::doubletostring(optionsp->lower1,4);
  ST::string nl2 = ST::doubletostring(optionsp->lower2,4);
  ST::string nu1 = ST::doubletostring(optionsp->upper1,4);
  ST::string nu2 = ST::doubletostring(optionsp->upper2,4);
  nl1 = nl1.replaceallsigns('.','p');
  nl2 = nl2.replaceallsigns('.','p');
  nu1 = nu1.replaceallsigns('.','p');
  nu2 = nu2.replaceallsigns('.','p');

  double help;

  optionsp->out("  Estimation results for the scale parameter:\n",true);
  optionsp->out("\n");


  ST::string vstr;

  vstr = "  Mean:         ";
  optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(FCsigma2.betamean(0,0),6) + "\n");

  vstr = "  Std. dev.:    ";
  if (FCsigma2.betavar(0,0) < 0)
    help = 0;
  else
    help = sqrt(FCsigma2.betavar(0,0));
  optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
  ST::doubletostring(help,6) + "\n");

  vstr = "  " + l1 + "% Quantile: ";
  optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
  ST::doubletostring(FCsigma2.betaqu_l1_lower(0,0),6) + "\n");

  vstr = "  " + l2 + "% Quantile: ";
  optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
  ST::doubletostring(FCsigma2.betaqu_l2_lower(0,0),6) + "\n");

  vstr = "  50% Quantile: ";
  optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
  ST::doubletostring(FCsigma2.betaqu50(0,0),6) + "\n");

  vstr = "  " + u1 + "% Quantile: ";
  optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
  ST::doubletostring(FCsigma2.betaqu_l2_upper(0,0),6) + "\n");

  vstr = "  " + u2 + "% Quantile: ";
  optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
  ST::doubletostring(FCsigma2.betaqu_l1_upper(0,0),6) + "\n");





/*
  outscale << "pmean   pstddev   pqu" << nl1 << "   pqu" << nl2 <<
              "   pqu50   pqu" << nu1 << "   pqu" << nu2 << endl;

        ST::string vstr;

        vstr = "  Mean:         ";
        help = Scalesave.get_betamean(0,0);
        optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(help,6) + "\n");
        outscale << help << "  ";
        results_latex.push_back("Mean  \\> " +
                            ST::doubletostring(help,6) + " \\\\");


        vstr = "  Std. dev.:    ";
        if (Scalesave.get_betavar(0,0) < 0)
          help = 0;
        else
          help = sqrt(Scalesave.get_betavar(0,0));
        optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(help,6) + "\n");
        outscale << help << "  ";
        results_latex.push_back("Std. dev.:  \\> " +
                            ST::doubletostring(help,6) + " \\\\");

        vstr = "  " + l1 + "% Quantile: ";
        help = Scalesave.get_beta_lower1(0,0);
        optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(help,6) + "\n");
        outscale << help << "  ";
        results_latex.push_back(vstr.insert_string_char(hchar,helpstring) + " \\>" +
                            ST::doubletostring(help,6) + " \\\\");


        vstr = "  " + l2 + "% Quantile: ";
        help = Scalesave.get_beta_lower2(0,0);
        optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(help,6) + "\n");
        outscale << help << "  ";
        results_latex.push_back(vstr.insert_string_char(hchar,helpstring) + " \\>" +
                            ST::doubletostring(help,6) + " \\\\");

        vstr = "  50% Quantile: ";
        help = Scalesave.get_betaqu50(0,0);
        optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(help,6) + "\n");
        outscale << help << "  ";
        results_latex.push_back(vstr.insert_string_char(hchar,helpstring) + " \\>" +
                            ST::doubletostring(help,6) + " \\\\");

        vstr = "  " + u1 + "% Quantile: ";
        help = Scalesave.get_beta_upper2(0,0);
        optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(help,6) + "\n");
        outscale << help << "  ";
        results_latex.push_back(vstr.insert_string_char(hchar,helpstring) + " \\>" +
                            ST::doubletostring(help,6) + " \\\\");

        vstr = "  " + u2 + "% Quantile: ";
        help = Scalesave.get_beta_upper1(0,0);
        optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(help,6) + "\n");
        outscale << help << "  ";
        results_latex.push_back(vstr.insert_string_char(hchar,helpstring) + " \\>" +
                            ST::doubletostring(help,6) + " \\\\");

        outscale << endl;

        optionsp->out("\n");
*/

  optionsp->out("\n");



  }



//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_gaussian_re -----------------------
//------------------------------------------------------------------------------
/*
void DISTRIBUTION_gaussian_re::set_constscale(double s)
  {
  scale(0,0) = s/(trmult(0,0)*trmult(0,0));
  constscale=true;
  }

void DISTRIBUTION_gaussian_re::undo_constscale(void)
  {
  constscale = false;
  }

void DISTRIBUTION_gaussian_re::set_uniformprior(void)
  {
  uniformprior=true;
  }



void DISTRIBUTION_gaussian_re::outoptions(void)
  {
  DISTRIBUTION::outoptions();
  optionsp->out("  Response function: identity\n");
  if(uniformprior)
    {
    optionsp->out("  Uniform prior on sigma\n");
    }
  else
    {
    optionsp->out("  Hyperparameter a: " + ST::doubletostring(a_invgamma,6) + "\n");
    optionsp->out("  Hyperparameter b: " + ST::doubletostring(b_invgamma,6) + "\n");
    }
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTRIBUTION_gaussian_re::update(void)
  {

  register unsigned i;

  double help;

  double * worklin;
  double * workresp;
  double * workweight;


  // scaleparameter

  double sum = 0;

  worklin = (*linpred_current).getV();
  workresp = response.getV();
  workweight = weight.getV();

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {

    help = *workresp - *worklin;
    sum += *workweight*help*help;
    }

  if(uniformprior==true)
    {
    double help = 1000000;
    while (help > 200000)
      help = rand_invgamma(-0.5+0.5*nrobsmweightzero,0.5*sum);
    scale(0,0) = help;
    }
  else
    {
    if(shrinkage)
      {
      scale(0,0) = rand_invgamma(a_invgamma+0.5*nrobsmweightzero + 0.5*nrridge + 0.5*nrlasso,
                     b_invgamma+0.5*sum + 0.5*ridgesum + 0.5*lassosum);
      }
    else
      {
      scale(0,0) = rand_invgamma(a_invgamma+0.5*nrobsmweightzero,
                     b_invgamma+0.5*sum);
      }
    }

  // Prediction

  if (predictresponse==true)
    {
    worklin = (*linpred_current).getV();
    workresp = response.getV();
    workweight = weight.getV();
    double * workpredictind = predictindicator.getV();
    double sscale = sqrt(scale(0,0));

    for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++,workpredictind++)
      {
      if (*workpredictind == 0)
        {
        *workresp = *worklin + (sscale/(*workweight))*rand_normal();
        }

      }
    }


  DISTRIBUTION::update();

  }


double DISTRIBUTION_gaussian_re::loglikelihood(double * res,
                       double * lin,
                       double * w,
                       const int & i) const
  {
  double help = *res-*lin;
  return  - *w * ( help * help )/(2* scale(0,0));
  }


bool DISTRIBUTION_gaussian_re::posteriormode(void)
  {

  unsigned i;

  double * worklin = (*linpred_current).getV();
  double * workresp = response.getV();
  double * workweight = weight.getV();

  double sum = 0;
  double help;
  double sumweight=0;

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {
    help = *workresp - *worklin;
    sum += *workweight*help*help;
    sumweight += *workweight;
    }


  if (constscale==false)
    scale(0,0) = (1.0/sumweight)*sum;


  return true;

  }


// constructor without offset
DISTRIBUTION_gaussian_re::DISTRIBUTION_gaussian_re(const double & a,
                                             const double & b,
                                             MCMCoptions * o,
                                             const datamatrix & r,
                                             const ST::string & p,
                                             const ST::string & ps,
                                             const datamatrix & w)
  : DISTRIBUTION(o,r,w,p,ps)

  {
// auskommentiert um nichtinformative Prioris zu erlauben
//  assert (a > 0);
//  assert (b > 0);

  constscale=false;
  uniformprior=false;

  a_invgamma = a;
  b_invgamma = b;
  family = "Gaussian_re";

  acceptancescale=100;

  constant_iwlsweights=true;
  iwlsweights_notchanged_df = true;

  }


// constructor with offset
DISTRIBUTION_gaussian_re::DISTRIBUTION_gaussian_re(const datamatrix & offset,
                      const double & a, const double & b, MCMCoptions * o,
                      const datamatrix & r,const ST::string & p,
                      const ST::string & ps,const datamatrix & w)
  : DISTRIBUTION(offset,o,r,w,p,ps)

  {

// auskommentiert um nichtinformative Prioris zu erlauben
//  assert (a > 0);
//  assert (b > 0);

  constscale=false;
  uniformprior=false;

  a_invgamma = a;
  b_invgamma = b;
  family = "Gaussian_re";

  acceptancescale=100;

  constant_iwlsweights=true;
  iwlsweights_notchanged_df = true;

  }


void DISTRIBUTION_gaussian_re::set_distrpointer(DISTRIBUTION * dp)
  {
  distrp = dp;
  trmult = distrp->get_trmultmat();

  datamatrix tr(1,1,trmult(0,0)*trmult(0,0));
  Scalesave.set_transformmult(tr);

  }

const DISTRIBUTION_gaussian_re & DISTRIBUTION_gaussian_re::operator=(
                                      const DISTRIBUTION_gaussian_re & nd)
  {
  if (this==&nd)
    return *this;
  DISTRIBUTION::operator=(DISTRIBUTION(nd));
  a_invgamma = nd.a_invgamma;
  b_invgamma = nd.b_invgamma;
  constscale = nd.constscale;
  uniformprior = nd.uniformprior;
  distrp = nd.distrp;
  return *this;
  }


DISTRIBUTION_gaussian_re::DISTRIBUTION_gaussian_re(const DISTRIBUTION_gaussian_re & nd)
   : DISTRIBUTION(DISTRIBUTION(nd))
     {
     a_invgamma = nd.a_invgamma;
     b_invgamma = nd.b_invgamma;
     constscale = nd.constscale;
     uniformprior = nd.uniformprior;
     distrp = nd.distrp;
     }


void DISTRIBUTION_gaussian_re::compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           const datamatrix & scale,const int & i) const
  {
  if ((*weight != 0))
    {
    double s = scale(0,0)*pow(trmult(0,0),2);
    double r = *response*trmult(0,0)-*mu;
    *deviance =  (*weight/s)*r*r+log(2*M_PI*s/(*weight));
    *deviancesat = (*weight/s)*r*r;
    }
  else
    {
    *deviance = 0;
    *deviancesat = 0;
    }
  }


double DISTRIBUTION_gaussian_re::compute_weight(double * linpred, double * weight,
                        const int & i, const unsigned & col) const
  {
  return *weight;
  }


void DISTRIBUTION_gaussian_re::compute_mu(const double * linpred,double * mu) const
  {
  *mu = trmult(0,0)* *linpred;
  }


void DISTRIBUTION_gaussian_re::compute_mu_notransform(
const double * linpred,double * mu) const
  {
  *mu = *linpred;
  }


void  DISTRIBUTION_gaussian_re::tr_nonlinear(vector<double *> b,
                                          vector<double *> br,
                                          vector<FULLCOND*> & fcp,
                                          unsigned & nr,
                                          unsigned & it,ST::string & trtype)
  {
  if (trtype == "exp")
    DISTRIBUTION::tr_nonlinear(b,br,fcp,nr,it,trtype);
  else if(trtype == "lognormal")
    {
    datamatrix help(1,1);
    Scalesave.readsample2(help,it);
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = exp(interceptsample(it,0)+*b[i]+help(0,0)/2.0);
      }
    }
  else if (trtype == "elasticity")
    {
    if (b.size() == 2)
      {
      *br[1] = *b[1] * fcp[0]->get_data(nr,0) /(interceptsample(it,0)+ *b[0]);
      }
    }
  else if (trtype=="marginal")
    {
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = interceptsample(it,0)+ *b[i];
      }
    }
  else if (trtype=="marginalintercept")
    {
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = interceptsample(it,0);
      }
    }
  else if (trtype == "lognormalintercept")
    {
    datamatrix help(1,1);
    Scalesave.readsample2(help,it);
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = exp(interceptsample(it,0)+help(0,0)/2.0);
      }
    }

  }
*/


} // end: namespace MCMC



