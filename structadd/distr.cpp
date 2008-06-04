

#include "distr.h"

namespace MCMC
{


bool DISTR::check_workingweights_one(void)
  {
  unsigned i=0;
  double * work_workingweight = workingweight.getV();
  bool one = true;
  while (i<nrobs && one == true)
    {
    if (*work_workingweight != 1)
      one = false;
    work_workingweight++;
    i++;
    }
  return one;  
  }


DISTR::DISTR(GENERAL_OPTIONS * o, const datamatrix & r,
             const datamatrix & w)
  {

  sigma2=1;

  optionsp = o;

  family = "unknown";
  updateIWLS = false;

  response = r;
  workingresponse = r;
  responsename = "Y";

  nrobs = response.rows();

  if (w.rows() == 1)
    {
    weight = datamatrix(r.rows(),1,1);
    weights_one=true;
    }
  else
    {
    weight = w;
    }

  workingweight = weight;
  weights_one = check_workingweights_one();

  changingweight = false;

  weightname = "W";

  linearpred1 = datamatrix(nrobs,1,0);
  linearpred2 = datamatrix(nrobs,1,0);

  linpred_current = 1;

  trmult=1;

  }



DISTR::DISTR(const DISTR & d)
  {

  sigma2 = d.sigma2;

  optionsp = d.optionsp;
  nrobs = d.nrobs;

  response = d.response;
  workingresponse = d.workingresponse;
  responsename = d.responsename;

  weight = d.weight;
  weightname = d.weightname;

  workingweight = d.workingweight;
  changingweight = d.changingweight;
  weights_one = d.weights_one;

  linearpred1 = d.linearpred1;
  linearpred2 = d.linearpred2;
  linpred_current = d.linpred_current;

  updateIWLS = d.updateIWLS;
  family = d.family;

  trmult=d.trmult;
  }


const DISTR & DISTR::operator=(const DISTR & d)
  {
  if (this == &d)
    return *this;

  sigma2 = d.sigma2;

  optionsp = d.optionsp;
  nrobs = d.nrobs;

  response = d.response;
  workingresponse = d.workingresponse;
  responsename = d.responsename;

  weight = d.weight;
  weightname = d.weightname;

  workingweight = d.workingweight;
  changingweight = d.changingweight;
  weights_one = d.weights_one;

  linearpred1 = d.linearpred1;
  linearpred2 = d.linearpred2;
  linpred_current = d.linpred_current;

  updateIWLS = d.updateIWLS;  
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
    {
    if (linpred_current==1)
      worklin = linearpred1.getV();
    else
      worklin = linearpred2.getV();
    }
  else
    {
    if (linpred_current==1)
      worklin = linearpred2.getV();
    else
      worklin = linearpred1.getV();
    }

  for (i=0;i<nrobs;i++,workweight++,worklin++,workres++)
    help += loglikelihood(workres,worklin,workweight);

  return help;

  }


void DISTR::compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           double * scale) const
  {

  }



void DISTR::compute_mu(const double * linpred,double * mu)
  {

  }


void DISTR::swap_linearpred(void)
  {

  if (linpred_current==1)
    linpred_current=2;
  else
    linpred_current=1;
  }


void DISTR::update(void)
  {
  } // end: update


bool DISTR::posteriormode(void)
  {
  double h = compute_iwls(true,false);
  return true;
  }


double DISTR::compute_iwls(const bool & current, const bool & like)
  {

  register unsigned  i;

  double * workweight = weight.getV();
  double * workresponse = response.getV();

  double * worklin;
  if (current)
    {
    if (linpred_current == 1)
      worklin = linearpred1.getV();
    else
      worklin = linearpred2.getV();
    }
  else  // use porposed
    {
    if (linpred_current == 1)
      worklin = linearpred2.getV();
    else
      worklin = linearpred1.getV();
    }

  double * work_workingresponse=workingresponse.getV();
  double * work_workingweight = workingweight.getV();

  double likelihood = 0;

  for (i=0;i<nrobs;i++,workweight++,work_workingweight++,workresponse++,
          work_workingresponse++,worklin++)
    {

    likelihood += compute_iwls(workresponse,worklin,
                               workweight,work_workingweight,
                               work_workingresponse,like);
    }

  // TEST
  /*
  ofstream out("c:\\bayesx\\test\\results\\workresponse.res");
  workingresponse.prettyPrint(out);

  ofstream out2("c:\\bayesx\\test\\results\\workweight.res");
  workingweight.prettyPrint(out2);
  */
  // TEST

  return likelihood;
  }


void DISTR::outresults(ST::string pathresults)
  {
  optionsp->out("\n");
  }


void DISTR::reset(void)
  {
  linearpred1 = datamatrix(nrobs,1,0);
  linearpred2 = datamatrix(nrobs,1,0);
  linpred_current = 1;
  }


double DISTR::get_scale(void)
  {
  return sigma2;
  }


double DISTR::get_scalemean(void)
  {
  return sigma2;
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
//  trmult=1;

  FCsigma2 = FC(o,"Gaussian variance parameter",1,1,ps);
  FCsigma2.transform(0,0) = pow(trmult,2);

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
  return *this;
  }



DISTR_gaussian::DISTR_gaussian(const DISTR_gaussian & nd)
   : DISTR(DISTR(nd))
  {
  a_invgamma = nd.a_invgamma;
  b_invgamma = nd.b_invgamma;
  FCsigma2 = nd.FCsigma2;
  }



void DISTR_gaussian::standardise(void)
  {

  trmult = sqrt(response.var(0,weight));

  unsigned i;
  double * workresp = workingresponse.getV();
  double * worklin = linearpred1.getV();
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

  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  workresp = workingresponse.getV();
  workweight = weight.getV();

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {
    help = *workresp - *worklin;
    sum += *workweight*pow(help,2);
    }

  sigma2  = rand_invgamma(a_invgamma+0.5*nrobs,
                          b_invgamma+0.5*sum);

  FCsigma2.beta(0,0) = sigma2;
  FCsigma2.acceptance++;
  FCsigma2.update();

  DISTR::update();

  }


void DISTR_gaussian::compute_mu(const double * linpred,double * mu)
  {

    *mu = trmult * (*linpred);
  }


void DISTR_gaussian::compute_deviance(const double * response,
                                 const double * weight, const double * mu,
                                 double * deviance, double * deviancesat,
                                 double * scale) const
  {

  double s = *scale*pow(trmult,2);
  double r = *response*trmult-*mu;
  *deviance =  (*weight/s)*r*r+log(2*M_PI*s/(*weight));
  *deviancesat = (*weight/s)*r*r;
  }



double DISTR_gaussian::loglikelihood(double * res, double * lin,
                                     double * w) const
  {
  double help = *res-*lin;
  return  - *w * (pow(help,2))/(2* sigma2);
  }


double DISTR_gaussian::compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse, const bool & like)
  {
  *workingweight=*weight;
  *workingresponse = *response;
  }

bool DISTR_gaussian::posteriormode(void)
  {

  unsigned i;

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * workresp = workingresponse.getV();
  double * workweight = weight.getV();

  double sum = 0;
  double sumweight=0;
  double help;

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {
    help = *workresp - *worklin;
    sum += *workweight*pow(help,2);
    sumweight+=*workweight;
    }

  sigma2 = (1.0/sumweight)*sum;

  FCsigma2.beta(0,0) = sigma2;

  FCsigma2.posteriormode_betamean();

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

  if (optionsp->samplesize > 1)
    {
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
    }


  if (pathresults.isvalidfile() != 1)
    {

    optionsp->out("\n");

    optionsp->out("  Results for variance parameter are also stored in file\n");
    optionsp->out("  " +  pathresults + "\n");
    optionsp->out("\n");

    ofstream outscale(pathresults.strtochar());

    if (optionsp->samplesize > 1)
      {
      outscale << "pmean   pstddev   pqu" << nl1 << "   pqu" << nl2 <<
                "   pqu50   pqu" << nu1 << "   pqu" << nu2 << endl;
      }
    else
      {
      outscale << "pmean" << endl;
      }

    outscale << FCsigma2.betamean(0,0) << "  ";

    if (optionsp->samplesize > 1)
      {
      if (FCsigma2.betavar(0,0) < 0)
        help = 0;
      else
        help = sqrt(FCsigma2.betavar(0,0));
      outscale << help << "  ";

      outscale << FCsigma2.betaqu_l1_lower(0,0) << "  ";

      outscale << FCsigma2.betaqu_l2_lower(0,0) << "  ";

      outscale << FCsigma2.betaqu50(0,0) << "  ";

      outscale << FCsigma2.betaqu_l2_upper(0,0) << "  ";

      outscale << FCsigma2.betaqu_l1_upper(0,0) << "  ";
      }

    outscale << endl;
    }

  optionsp->out("\n");


  }


//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_gaussian_re -----------------------
//------------------------------------------------------------------------------


DISTR_gaussian_re::DISTR_gaussian_re(GENERAL_OPTIONS * o,const datamatrix & r,
                                     const datamatrix & w)
  : DISTR_gaussian(1,1,o,r,"",w)

  {

  family = "Gaussian_random_effect";

  }


const DISTR_gaussian_re & DISTR_gaussian_re::operator=(
                                      const DISTR_gaussian_re & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gaussian::operator=(DISTR_gaussian(nd));
  return *this;
  }



DISTR_gaussian_re::DISTR_gaussian_re(const DISTR_gaussian_re & nd)
   : DISTR_gaussian(DISTR_gaussian(nd))
  {
  }



void DISTR_gaussian_re::update(void)
  {

  // TEST
  //  ofstream out("c:\\bayesx\\test\\results\\response_RE.res");
  //  response.prettyPrint(out);
  // ENDE TEST

  }


bool DISTR_gaussian_re::posteriormode(void)
  {

  // TEST
  /*
    ofstream out("c:\\bayesx\\test\\results\\response_RE.res");
    response.prettyPrint(out);
  */  
  // ENDE TEST

  return true;
  }

void DISTR_gaussian_re::outoptions(void)
  {
  optionsp->out("RANDOM EFFECTS DISTRIBUTION:\n",true);
  optionsp->out("\n");
  optionsp->out("  Family: " + family + "\n");
  optionsp->out("  Number of clusters: " + ST::inttostring(nrobs) + "\n");
  optionsp->out("\n");  
  }


void DISTR_gaussian_re::outresults(ST::string pathresults)
  {
  }

} // end: namespace MCMC



