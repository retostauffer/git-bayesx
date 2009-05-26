

#include "distr.h"

namespace MCMC
{


bool DISTR::check_weightsone(void)
  {
  unsigned i=0;
  double * work_weight = weight.getV();
  bool one = true;
  while (i<nrobs && one == true)
    {
    if (*work_weight != 1)
      one = false;
    work_weight++;
    i++;
    }
  return one;
  }


unsigned DISTR::compute_nrzeroweights(void)
  {
  unsigned i=0;
  double * work_weight = weight.getV();
  unsigned nr=0;
  for (i=0;i<weight.rows();i++,work_weight++)
    {
    if (*work_weight == 0)
      nr++;
    }
  return nr;
  }
  

DISTR::DISTR(GENERAL_OPTIONS * o, const datamatrix & r,
             const datamatrix & w)
  {

  option1 = "";

  sigma2=1;

  optionsp = o;

  family = "unknown";
  updateIWLS = false;

  response = r;
  response_untransformed = r;
  workingresponse = r;
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

  weightsone = check_weightsone();
  nrzeroweights = compute_nrzeroweights();

  wtype = wweightschange_weightsneqone;

  weightname = "W";

  linearpred1 = datamatrix(nrobs,1,0);
  linearpred2 = datamatrix(nrobs,1,0);

  linpred_current = 1;

  trmult=1;

  meaneffect = 0;

  }



DISTR::DISTR(const DISTR & d)
  {

  option1 = d.option1;
  optionbool1 = d.optionbool1;

  sigma2 = d.sigma2;

  optionsp = d.optionsp;
  nrobs = d.nrobs;

  response = d.response;
  response_untransformed = d.response_untransformed;
  workingresponse = d.workingresponse;
  responsename = d.responsename;

  weight = d.weight;
  weightname = d.weightname;
  nrzeroweights = d.nrzeroweights;  

  workingweight = d.workingweight;
  weightsone = d.weightsone;
  wtype = d.wtype;

  linearpred1 = d.linearpred1;
  linearpred2 = d.linearpred2;
  linpred_current = d.linpred_current;

  updateIWLS = d.updateIWLS;
  family = d.family;

  trmult=d.trmult;

  meaneffect = d.meaneffect;
  }


const DISTR & DISTR::operator=(const DISTR & d)
  {
  if (this == &d)
    return *this;

  option1 = d.option1;
  optionbool1 = d.optionbool1;

  sigma2 = d.sigma2;

  optionsp = d.optionsp;
  nrobs = d.nrobs;

  response = d.response;
  response_untransformed = d.response_untransformed;
  workingresponse = d.workingresponse;
  responsename = d.responsename;

  weight = d.weight;
  weightname = d.weightname;
  nrzeroweights = d.nrzeroweights;

  workingweight = d.workingweight;
  weightsone = d.weightsone;
  wtype = d.wtype;

  linearpred1 = d.linearpred1;
  linearpred2 = d.linearpred2;
  linpred_current = d.linpred_current;

  updateIWLS = d.updateIWLS;
  family = d.family;

  trmult=d.trmult;

  meaneffect = d.meaneffect;

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

  if (weightsone==true)
    {
    for (i=0;i<nrobs;i++,worklin++,workres++)
      help += loglikelihood_weightsone(workres,worklin);
    }
  else
    {
    for (i=0;i<nrobs;i++,workweight++,worklin++,workres++)
      help += loglikelihood(workres,worklin,workweight);
    }

  return help;

  }


double DISTR::loglikelihood(int & begin,
int & end, statmatrix<double *> & responsep,
statmatrix<double *> & workingweightp, statmatrix<double *> & linpredp) const
  {
  double help=0;
  int i;

  double * * workresponsep = responsep.getV()+begin;
  double * * work_workingweightp = workingweightp.getV()+begin;
  double * * work_linpredp  = linpredp.getV()+begin;

  for (i=begin;i<=end;i++,work_workingweightp++,work_linpredp++,workresponsep++)
    help += loglikelihood(*workresponsep,*work_linpredp,*work_workingweightp);

  return help;


  }


double DISTR::compute_iwls_loglikelihood(int & begin,
int & end, statmatrix<double *> & responsep,
statmatrix<double *> & workingresponsep,
statmatrix<double *> & weightp,
statmatrix<double *> & workingweightp, statmatrix<double *> & linpredp
)
  {
  double help=0;
  int i;

  double * * workresponsep = responsep.getV()+begin;
  double * * work_workingresponsep = workingresponsep.getV()+begin;
  double * * work_weightp = weightp.getV()+begin;
  double * * work_workingweightp = workingweightp.getV()+begin;
  double * * work_linpredp  = linpredp.getV()+begin;

  for (i=begin;i<=end;i++,work_workingweightp++,work_linpredp++,
       workresponsep++,work_weightp++,work_workingresponsep++)
    {

    help += compute_iwls(*workresponsep,*work_linpredp,*work_weightp,
                         *work_workingweightp,*work_workingresponsep,true);
    }

  return help;


  }



double DISTR::compute_iwls_loglikelihood_sumworkingweight(int & begin,
int & end, statmatrix<double *> & responsep,
statmatrix<double *> & workingresponsep,
statmatrix<double *> & weightp,
statmatrix<double *> & workingweightp, statmatrix<double *> & linpredp,
datamatrix & intvar2,
double & sumworkingweight)
  {
  double help=0;
  int i;

  double * * workresponsep = responsep.getV()+begin;
  double * * work_workingresponsep = workingresponsep.getV()+begin;
  double * * work_weightp = weightp.getV()+begin;
  double * * work_workingweightp = workingweightp.getV()+begin;
  double * * work_linpredp  = linpredp.getV()+begin;

  sumworkingweight = 0;


  if (wtype==wweightschange_weightsneqone)
    {

    for (i=begin;i<=end;i++,work_workingweightp++,work_linpredp++,
         workresponsep++,work_weightp++,work_workingresponsep++)
      {

      help += compute_iwls(*workresponsep,*work_linpredp,*work_weightp,
                         *work_workingweightp,*work_workingresponsep,true);

      sumworkingweight+= *(*work_workingweightp);
      }

    }
  else if (wtype==wweightschange_weightsone)
    {

    for (i=begin;i<=end;i++,work_workingweightp++,work_linpredp++,
         workresponsep++,work_workingresponsep++)
      {

      compute_iwls_wweightschange_weightsone(*workresponsep,
                         *work_linpredp, *work_workingweightp,
                         *work_workingresponsep,help,true);

      sumworkingweight+= *(*work_workingweightp);
      }

    }
  else if (wtype==wweightsnochange_constant)
    {

    for (i=begin;i<=end;i++,work_workingweightp++,work_linpredp++,
         workresponsep++,work_workingresponsep++)
      {

      compute_iwls_wweightsnochange_constant(*workresponsep,
                         *work_linpredp, *work_workingweightp,
                         *work_workingresponsep,help,true);

      sumworkingweight+= *(*work_workingweightp);
      }

    }
  else if (wtype==wweightsnochange_one)
    {


    for (i=begin;i<=end;i++,work_linpredp++,
         workresponsep++,work_workingresponsep++)
      {

      compute_iwls_wweightsnochange_one(*workresponsep,
                         *work_linpredp, *work_workingresponsep,help,true);

      }

    sumworkingweight = begin-end+1;

    }

  return help;

  }


void DISTR::compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           double * scale) const
  {

  }



void DISTR::compute_mu(const double * linpred,double * mu,
                       bool notransform)
  {

  }


double DISTR::compute_MSE(const double * response, const double * weight,
                          const double * linpred)
  {
  return pow(*response-*linpred,2); 

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

  if (wtype==wweightschange_weightsneqone)
    {

    for (i=0;i<nrobs;i++,workweight++,work_workingweight++,workresponse++,
          work_workingresponse++,worklin++)
      {

      likelihood += compute_iwls(workresponse,worklin,
                                 workweight,work_workingweight,
                                 work_workingresponse,like);
      }

    }
  else if (wtype==wweightschange_weightsone)
    {

    for (i=0;i<nrobs;i++,work_workingweight++,workresponse++,
          work_workingresponse++,worklin++)
      {

      compute_iwls_wweightschange_weightsone(workresponse,worklin,
                                 work_workingweight,work_workingresponse,
                                 likelihood,like);
      }


    }
  else if (wtype==wweightsnochange_constant)
    {

    for (i=0;i<nrobs;i++,work_workingweight++,workresponse++,
          work_workingresponse++,worklin++)
      {

      compute_iwls_wweightsnochange_constant(workresponse,worklin,
                                 work_workingweight,work_workingresponse,
                                 likelihood,like);
      }

    }
  else if (wtype==wweightsnochange_one)
    {

    for (i=0;i<nrobs;i++,workresponse++,
          work_workingresponse++,worklin++)
      {

      compute_iwls_wweightsnochange_one(workresponse,worklin,
                                        work_workingresponse,
                                        likelihood,like);
      }

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


double DISTR::get_scale(bool transform)
  {
  if (!transform)
    return sigma2;
  else
    return sigma2*pow(trmult,2);
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


  if (check_weightsone())
    wtype = wweightsnochange_one;
  else
    wtype = wweightsnochange_constant;

  a_invgamma = a;
  b_invgamma = b;
  family = "Gaussian";

  standardise();

  FCsigma2 = FC(o,"",1,1,ps,false);
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
  double * resp_p = response.getV();
  double * worklin = linearpred1.getV();
  for (i=0;i<nrobs;i++,workresp++,worklin++,resp_p++)
   {
   *workresp = *workresp/trmult;
   *resp_p = (*resp_p)/trmult;
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

  sigma2  = rand_invgamma(a_invgamma+0.5*(nrobs-nrzeroweights),
                          b_invgamma+0.5*sum);

  FCsigma2.beta(0,0) = sigma2;
  FCsigma2.acceptance++;
  FCsigma2.update();

  DISTR::update();

  }


void DISTR_gaussian::compute_mu(const double * linpred,double * mu,
                                bool notransform)
  {
    if (!notransform)
      *mu = trmult * (*linpred);
    else
      *mu = (*linpred);
  }


void DISTR_gaussian::compute_deviance(const double * response,
                                 const double * weight, const double * mu,
                                 double * deviance, double * deviancesat,
                                 double * scale) const
  {
  if (*weight == 0)
    {
    *deviance = 0;
    *deviancesat = 0;
    }
  else
    {
    double r = *response-*mu;
    *deviance =  (*weight/(*scale))*r*r+log(2*M_PI*(*scale)/(*weight));
    *deviancesat = (*weight/(*scale))*r*r;
    }
  }



double DISTR_gaussian::loglikelihood(double * res, double * lin,
                                     double * w) const
  {
  if (*w==0)
    return 0;
  else
    {
    double help = *res-*lin;
    return  - *w * (pow(help,2))/(2* sigma2);
    }
  }


double DISTR_gaussian::loglikelihood_weightsone(double * res, double * lin) const
  {
  double help = *res-*lin;
  return  - (pow(help,2))/(2* sigma2);
  }


double DISTR_gaussian::compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse, const bool & like)
  {
  *workingweight=*weight;
  *workingresponse = *response;
  if (like && (*weight != 0))
    return  - *weight * (pow(*response-(*linpred),2))/(2* sigma2);
  else
    return 0;
  }


void DISTR_gaussian::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  *workingweight=1;
  *workingresponse = *response;
  if (compute_like)
    like -=  (pow(*response-(*linpred),2))/(2* sigma2);

  }


void DISTR_gaussian::compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  *workingresponse = *response;
  if (compute_like && *workingweight!=0)
    like  =- *workingweight * (pow(*response-(*linpred),2))/(2* sigma2);
  }


void DISTR_gaussian::compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  *workingresponse = *response;
  if (compute_like)
    like -=  (pow(*response-(*linpred),2))/(2* sigma2);
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

  ofstream out1;
  ofstream out2;

  FCsigma2.outresults(out1,out2,"");


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

  optionsp->out("  SCALE PARAMETER:\n",true);
  optionsp->out("\n");


  ST::string vstr;

  vstr = "    Mean:         ";
  optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(FCsigma2.betamean(0,0),6) + "\n");

  if (optionsp->samplesize > 1)
    {
    vstr = "    Std. dev.:    ";
    if (FCsigma2.betavar(0,0) < 0)
      help = 0;
    else
      help = sqrt(FCsigma2.betavar(0,0));
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(help,6) + "\n");

    vstr = "    " + l1 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(FCsigma2.betaqu_l1_lower(0,0),6) + "\n");

    vstr = "    " + l2 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(FCsigma2.betaqu_l2_lower(0,0),6) + "\n");

    vstr = "    50% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(FCsigma2.betaqu50(0,0),6) + "\n");

    vstr = "    " + u1 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(FCsigma2.betaqu_l2_upper(0,0),6) + "\n");

    vstr = "    " + u2 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(FCsigma2.betaqu_l1_upper(0,0),6) + "\n");
    }


  if (pathresults.isvalidfile() != 1)
    {

    optionsp->out("\n");

    optionsp->out("    Results for variance parameter are also stored in file\n");
    optionsp->out("    " +  pathresults + "\n");
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


double DISTR_gaussian::get_scalemean(void)
  {
  return FCsigma2.betamean(0,0);
  }


//------------------------------------------------------------------------------
//---------------------- CLASS DISTRIBUTION_loggaussian ------------------------
//------------------------------------------------------------------------------


DISTR_loggaussian::DISTR_loggaussian(const double & a,
                                             const double & b,
                                             GENERAL_OPTIONS * o,
                                             const datamatrix & r,
                                             const ST::string & ps,
                                             const datamatrix & w)
  : DISTR_gaussian(a,b,o,r,ps,w)

  {
  family = "log-Gaussian";
  }


const DISTR_loggaussian & DISTR_loggaussian::operator=(
                                      const DISTR_loggaussian & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gaussian::operator=(DISTR_gaussian(nd));
  return *this;
  }


DISTR_loggaussian::DISTR_loggaussian(const DISTR_loggaussian & nd)
   : DISTR_gaussian(DISTR_gaussian(nd))
  {
  }



void DISTR_loggaussian::compute_mu(const double * linpred,double * mu,
                                bool notransform)
  {
//  double scale = FCsigma2.beta(0,0);

  if (!notransform)
    *mu = exp(trmult * (*linpred) + sigma2*pow(trmult,2)/2.0);
  else
    *mu = exp((*linpred) + sigma2*pow(trmult,2)/2.0);
  }


double DISTR_loggaussian::compute_MSE(const double * response,
                                      const double * weight,
                                      const double * linpred)
  {
  double diff = exp(*response) - exp(*linpred +  FCsigma2.betamean(0,0)/2);
  return pow(diff,2);
  }


void DISTR_loggaussian::compute_deviance(const double * response,
                                 const double * weight, const double * mu,
                                 double * deviance, double * deviancesat,
                                 double * scale) const
  {
  double r = *response-*mu;

  if (*weight != 0)
    {
    *deviance =  (*weight/(*scale))*r*r+log(2*M_PI*(*scale)/(*weight));
    *deviancesat = (*weight/(*scale))*r*r;
    }

  else
    {
    *deviance = 0;
    *deviancesat = 0;
    }

  }

//------------------------------------------------------------------------------
//--------------------- CLASS DISTRIBUTION_gaussian_exp ------------------------
//------------------------------------------------------------------------------


DISTR_gaussian_exp::DISTR_gaussian_exp(const double & a,
                                             const double & b,
                                             GENERAL_OPTIONS * o,
                                             const datamatrix & r,
                                             const ST::string & ps,
                                             const datamatrix & w)
  : DISTR_gaussian(a,b,o,r,ps,w)

  {
  standardise();
//  changingworkingweights = true;
  updateIWLS = true;
  }



const DISTR_gaussian_exp & DISTR_gaussian_exp::operator=(
                                      const DISTR_gaussian_exp & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gaussian::operator=(DISTR_gaussian(nd));
  return *this;
  }


DISTR_gaussian_exp::DISTR_gaussian_exp(const DISTR_gaussian_exp & nd)
   : DISTR_gaussian(DISTR_gaussian(nd))
  {
  }



void DISTR_gaussian_exp::standardise(void)
  {

  trmult = 1;

  unsigned i;
  double * workresp = workingresponse.getV();
  double * worklin = linearpred1.getV();
  for (i=0;i<nrobs;i++,workresp++,worklin++)
    {
    *workresp = response(i,0);
    *worklin = 0;
    }

  FCsigma2.transform(0,0) = pow(trmult,2);    

  }



void DISTR_gaussian_exp::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: exponential\n");

  optionsp->out("  Hyperparameter a: " + ST::doubletostring(a_invgamma,6) + "\n");
  optionsp->out("  Hyperparameter b: " + ST::doubletostring(b_invgamma,6) + "\n");

  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussian_exp::update(void)
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

  workresp = response.getV();
  workweight = weight.getV();

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {
    help = *workresp - exp(*worklin);
    sum += *workweight*pow(help,2);
    }

  sigma2  = rand_invgamma(a_invgamma+0.5*nrobs,
                          b_invgamma+0.5*sum);

  FCsigma2.beta(0,0) = sigma2;
  FCsigma2.acceptance++;
  FCsigma2.update();

  DISTR::update();

  }


void DISTR_gaussian_exp::compute_mu(const double * linpred,double * mu,
                                    bool notransform)
  {
    if (!notransform)
      *mu = trmult * exp(*linpred);
    else
      *mu = exp(*linpred);
  }


double DISTR_gaussian_exp::loglikelihood(double * res, double * lin,
                                         double * w) const
  {
  double help = *res-exp(*lin);
  return  - *w * (pow(help,2))/(2* sigma2);
  }


double DISTR_gaussian_exp::compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse, const bool & like)
  {
  double mu = exp(*linpred);
  double mu2 = pow(mu,2);
  *workingweight=mu2 * (*weight)/sigma2;

  *workingresponse = (*response-mu)/mu + (*linpred);

  if (like==true)
    {
    double h = *response-mu;
    return  - (*weight) * pow(h,2)/(2* sigma2);
    }
  else
    return 0;
  }


bool DISTR_gaussian_exp::posteriormode(void)
  {

  unsigned i;

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * workresp = response.getV();
  double * workweight = weight.getV();

  double sum = 0;
  double sumweight=0;
  double help;

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {
    help = *workresp - exp(*worklin);
    sum += *workweight*pow(help,2);
    sumweight+=*workweight;
    }

  sigma2 = (1.0/sumweight)*sum;

  FCsigma2.beta(0,0) = sigma2;

  FCsigma2.posteriormode_betamean();

  return true;

  }


//------------------------------------------------------------------------------
//--------------------- CLASS DISTRIBUTION_gaussian_mult -----------------------
//------------------------------------------------------------------------------


DISTR_gaussian_mult::DISTR_gaussian_mult(const double & a,
                                             const double & b,
                                             GENERAL_OPTIONS * o,
                                             const datamatrix & r,
                                             const ST::string & ps,
                                             const datamatrix & w)
  : DISTR_gaussian_exp(a,b,o,r,ps,w)

  {
  standardise();
  optionbool1 = false;
//  changingworkingweights = false;
  updateIWLS = false;  
  }


void DISTR_gaussian_mult::set_mult(bool & m)
  {

  if (m==true)
    {
    optionbool1 = true;
//    changingworkingweights = true;
    updateIWLS = true;

    }
  else
    {
    optionbool1 = false;
//    changingworkingweights = false;
    updateIWLS = false;
    }

  }


const DISTR_gaussian_mult & DISTR_gaussian_mult::operator=(
                                      const DISTR_gaussian_mult & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gaussian_exp::operator=(DISTR_gaussian_exp(nd));
  return *this;
  }


DISTR_gaussian_mult::DISTR_gaussian_mult(const DISTR_gaussian_mult & nd)
   : DISTR_gaussian_exp(DISTR_gaussian_exp(nd))
  {
  }



void DISTR_gaussian_mult::standardise(void)
  {

  trmult = 1;

  unsigned i;
  double * workresp = workingresponse.getV();
  double * worklin = linearpred1.getV();
  for (i=0;i<nrobs;i++,workresp++,worklin++)
    {
    *workresp = response(i,0);
    *worklin = 0;
    }

  FCsigma2.transform(0,0) = pow(trmult,2);

  }



void DISTR_gaussian_mult::outoptions(void)
  {
  DISTR_gaussian::outoptions();
  }


void DISTR_gaussian_mult::update(void)
  {
  if (optionbool1==true)
    DISTR_gaussian_exp::update();
  else
    DISTR_gaussian::update();
  }


void DISTR_gaussian_mult::compute_mu(const double * linpred,double * mu,
                                    bool notransform)
  {
  if (!optionbool1)
    {
    DISTR_gaussian::compute_mu(linpred,mu,notransform);
    }
 else
   {
   DISTR_gaussian_exp::compute_mu(linpred,mu,notransform);
   }

  }


double DISTR_gaussian_mult::loglikelihood(double * res, double * lin,
                                         double * w) const
  {

  if (!optionbool1)
    {
    return DISTR_gaussian::loglikelihood(res,lin,w);
    }
 else
   {
   return DISTR_gaussian_exp::loglikelihood(res,lin,w);
   }

  }


double DISTR_gaussian_mult::compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse, const bool & like)
  {


  if (!optionbool1)
    {
    return DISTR_gaussian::compute_iwls(response,linpred,weight,workingweight,
                                         workingresponse,like);
    }
 else
   {
    return DISTR_gaussian_exp::compute_iwls(response,linpred,weight,workingweight,
                                         workingresponse,like);
   }

  }


bool DISTR_gaussian_mult::posteriormode(void)
  {

  if (!optionbool1)
    {
    return DISTR_gaussian::posteriormode();

    }
 else
   {
    return DISTR_gaussian_exp::posteriormode();
   }


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



