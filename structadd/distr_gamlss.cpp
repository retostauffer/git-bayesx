/* BayesX - Software for Bayesian Inference in
Structured Additive Regression Models.
Copyright (C) 2011  Christiane Belitz, Andreas Brezger,
Thomas Kneib, Stefan Lang, Nikolaus Umlauf

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */

#include "distr_gamlss.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_negbin_delta --------------------------
//------------------------------------------------------------------------------


DISTR_negbin_delta::DISTR_negbin_delta(GENERAL_OPTIONS * o,
                                       const datamatrix & r,
                                       double & stpsum, int & strmax,
                                       const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Negative_Binomial - delta";

  double responsemax = response.max(0);

  stopsum = stpsum;
  stoprmax = strmax;
  if (stoprmax < responsemax)
    stoprmax = responsemax;
  }


DISTR_negbin_delta::DISTR_negbin_delta(const DISTR_negbin_delta & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  E_dig_y_delta = nd.E_dig_y_delta;
  E_trig_y_delta = nd.E_trig_y_delta;
  delta = nd.delta;
  log_delta_div_delta_plus_mu = nd.log_delta_div_delta_plus_mu;
  lngamma_delta = nd.lngamma_delta;
  delta_plus_mu = nd.lngamma_delta;

  stopsum = nd.stopsum;
  stoprmax = nd.stoprmax;
  }


const DISTR_negbin_delta & DISTR_negbin_delta::operator=(
                            const DISTR_negbin_delta & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  E_dig_y_delta = nd.E_dig_y_delta;
  E_trig_y_delta = nd.E_trig_y_delta;
  delta = nd.delta;
  log_delta_div_delta_plus_mu = nd.log_delta_div_delta_plus_mu;
  lngamma_delta = nd.lngamma_delta;
  delta_plus_mu = nd.lngamma_delta;

  stopsum = nd.stopsum;
  stoprmax = nd.stoprmax;

  return *this;
  }



double DISTR_negbin_delta::get_intercept_start(void)
  {
  return 0;
  }


double DISTR_negbin_delta::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double delta;

  if (*linpred <= linpredlimit)
    delta = explinpredlimit;
  else
    delta = exp(*linpred);

  double resp_plus_delta = (*response) + delta;

  double log_mu_plus_delta = log((*worktransformlin[0]) + delta);

  modify_worklin();

  return  randnumbers::lngamma_exact(resp_plus_delta) -
          randnumbers::lngamma_exact(delta) -
          (resp_plus_delta)* log_mu_plus_delta  +
          delta*log(delta);

  }


void DISTR_negbin_delta::compute_expectation(void)
  {

  int k=1;
  double k_delta;
  double kplus1;
  double psum;

  double L = exp(delta*log_delta_div_delta_plus_mu);
  E_dig_y_delta = randnumbers::digamma_exact(delta)*L;
  E_trig_y_delta = randnumbers::trigamma_exact(delta)*L;

  psum = L;

  while ((psum < stopsum) && (k <=stoprmax))
    {
    k_delta = k + delta;
    kplus1 = k + 1;

    L = exp(randnumbers::lngamma_exact(k_delta) - randnumbers::lngamma_exact(kplus1)
            - lngamma_delta + delta*log_delta_div_delta_plus_mu +
            k* log(*worktransformlin[0]/delta_plus_mu) );

    psum += L;

    E_dig_y_delta += randnumbers::digamma_exact(k_delta)*L;

    E_trig_y_delta += randnumbers::trigamma_exact(k_delta)*L;

    k++;
    }


  }


void DISTR_negbin_delta::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  if (*linpred <= linpredlimit)
    delta = explinpredlimit;
  else
    delta = exp(*linpred);


  delta_plus_mu = delta + (*worktransformlin[0]);

  log_delta_div_delta_plus_mu = log(delta/delta_plus_mu);

  lngamma_delta = randnumbers::lngamma_exact(delta);

  double delta_plus_response = delta + (*response);

  double nu = delta*(randnumbers::digamma_exact(delta_plus_response) -
                    randnumbers::digamma_exact(delta) +
                     log_delta_div_delta_plus_mu +
                    (*worktransformlin[0]-(*response))/delta_plus_mu);

  compute_expectation();

  *workingweight = -delta*(log_delta_div_delta_plus_mu + (*worktransformlin[0])/delta_plus_mu)
                   -delta*(E_dig_y_delta - randnumbers::digamma_exact(delta))
                   -pow(delta,2)*(E_trig_y_delta - randnumbers::trigamma_exact(delta));


  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    double resp_plus_delta = (*response) + delta;

    like += randnumbers::lngamma_exact(resp_plus_delta) -
            lngamma_delta -
            resp_plus_delta*log(delta_plus_mu) +
            delta*log(delta);


    }

  modify_worklin();

  }


void DISTR_negbin_delta::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (delta): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_negbin_delta::update_end(void)
  {
  DISTR_gamlss::update_end();
  }



//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_negbin_mu -----------------------------
//------------------------------------------------------------------------------


DISTR_negbin_mu::DISTR_negbin_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Negative_Binomial - mu";
  }


DISTR_negbin_mu::DISTR_negbin_mu(const DISTR_negbin_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_negbin_mu & DISTR_negbin_mu::operator=(
                            const DISTR_negbin_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_negbin_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<double> scale) const
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_delta
   // *linpred[1] = eta_mu

   double delta = exp(*linpred[0]);
   double mu = exp(*linpred[1]);
   double resp_plus_one = (*response[1]) + 1;
   double delta_plus_mu = delta + mu;
   double delta_plus_response = delta+(*response[1]);

   double l = randnumbers::lngamma_exact(delta_plus_response) -
              randnumbers::lngamma_exact(resp_plus_one) -
              randnumbers::lngamma_exact(delta) +
              delta*log(delta/delta_plus_mu)+
              (*response[1])*log(mu/delta_plus_mu);


  *deviance = -2*l;


  }


double DISTR_negbin_mu::get_intercept_start(void)
  {
  return log(response.mean(0));
  }


double DISTR_negbin_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of delta equation
  // *worktransformlin[0] = exp(eta_delta);

  if (counter==0)
    {
    set_worklin();
    }

  double mu;

  if (*linpred <= linpredlimit)
    mu = explinpredlimit;
  else
    mu = exp(*linpred);

  modify_worklin();

  return - ((*worktransformlin[0]) + (*response))*
           log((*worktransformlin[0])+mu) +(*response)*log(mu);

  }


void DISTR_negbin_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of delta equation
  // *worktransformlin[0] = exp(eta_delta);

  if (counter==0)
    {
    set_worklin();
    }

  double mu;


  if (*linpred <= linpredlimit)
    mu = explinpredlimit;
  else
    mu = exp(*linpred);


  double delta_plus_mu = (*worktransformlin[0]) + mu;

  double nu = (*worktransformlin[0])*((*response)-mu)/delta_plus_mu;

  *workingweight = (*worktransformlin[0])*mu/delta_plus_mu;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    like += -((*worktransformlin[0])+(*response))*log(delta_plus_mu) +
            (*response)*log(mu);

    }

  modify_worklin();

  }


void DISTR_negbin_mu::compute_mu_mult(vector<double *> linpred,double * mu)
  {
  *mu = exp(*linpred[1]);
  }


void DISTR_negbin_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_negbin_mu::update_end(void)
  {
  DISTR_gamlss::update_end();
  }



//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_zip_cloglog_pi --------------------------
//------------------------------------------------------------------------------


DISTR_zip_cloglog_pi::DISTR_zip_cloglog_pi(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Zero_Inflated_Poisson - pi";
  helpmat1 = datamatrix(nrobs,1,1-exp(-exp(0)));
  }


DISTR_zip_cloglog_pi::DISTR_zip_cloglog_pi(const DISTR_zip_cloglog_pi & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_zip_cloglog_pi & DISTR_zip_cloglog_pi::operator=(
                            const DISTR_zip_cloglog_pi & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_zip_cloglog_pi::get_intercept_start(void)
  {
  return 0;
  }


double DISTR_zip_cloglog_pi::loglikelihood_weightsone(double * response,
                                                      double * linpred)
  {

  if (counter==0)
    set_worklin();

  double explinpi = exp(*linpred);
  double oneminuspi = 1 - exp(-explinpi);
  double pi = 1-oneminuspi;
  double expminuslambda = exp(-(*worktransformlin[0]));
  double denompart = pi+oneminuspi*expminuslambda;

  modify_worklin();

  if (*response == 0)
    return log(denompart);
  else
    return log(oneminuspi);

  }


void DISTR_zip_cloglog_pi::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu);

  if (counter==0)
    set_worklin();

//  double t = *worktransformlin[0];
//  double t2 = *worklin[0];

  double explinpi = exp(*linpred);
  double oneminuspi = 1 - exp(-explinpi);
  double pi = 1-oneminuspi;
  double expminuslambda = exp(-(*worktransformlin[0]));
  double denompart = pi+oneminuspi*expminuslambda;
  double denom =  denompart*oneminuspi;
  double explinpi_pi = explinpi*pi;

  double nu = explinpi_pi/oneminuspi;
  if (*response == 0)
    nu -= explinpi_pi/denom;

  *workingweight =  pow(explinpi,2)*pow(pi,2)*(1-expminuslambda)/denom;

  *workingresponse = *linpred + nu/(*workingweight);


  if (compute_like)
    {

    if (*response == 0)
      like += log(denompart);
    else
      like += log(oneminuspi);

    }


  modify_worklin();

  }


void DISTR_zip_cloglog_pi::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (pi): complementary log log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_zip_cloglog_pi::update_end(void)
  {

//  ofstream out("d:\\_sicher\\papzip\\results\\etapi.raw");

  // helpmat1 stores 1-pi

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();


  double * ppi = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,ppi++,worklin++)
    {
    *ppi = 1-exp(-exp(*worklin));
//    out << *worklin << "  " << *ppi << endl;
    }

  }


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_zip_cloglog_mu --------------------------
//------------------------------------------------------------------------------


DISTR_zip_cloglog_mu::DISTR_zip_cloglog_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Zero_Inflated_Poisson - lambda";
  }


DISTR_zip_cloglog_mu::DISTR_zip_cloglog_mu(const DISTR_zip_cloglog_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_zip_cloglog_mu & DISTR_zip_cloglog_mu::operator=(
                            const DISTR_zip_cloglog_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_zip_cloglog_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<double> scale) const
  {

  double l;
  double explinpi = exp(*linpred[0]);
  double pi = exp(-explinpi);
  double lambda = exp(*linpred[1]);

  if (*response[1]==0)
    {
    l=  log(pi+(1-pi)*exp(-lambda));
    }
  else // response > 0
    {
    double help1 = *response[1]+1;
    l= log(1-pi) + (*response[1])*(*linpred[1])- lambda
       - randnumbers::lngamma_exact(help1);
    }

  *deviance = -2*l;

  }


double DISTR_zip_cloglog_mu::get_intercept_start(void)
  {
  return log(response.mean(0));
  }


double DISTR_zip_cloglog_mu::loglikelihood_weightsone(double * response,
                                                      double * linpred)
  {

  // *worklin[0] = linear predictor of pi equation
  // *worktransformlin[0] = 1-pi = 1-exp(-exp(eta_pi));

  if (counter==0)
    {
    set_worklin();
    }

  double lambda;
  double expminuslambda;

  if (*linpred <= linpredlimit)
    {
    lambda  = explinpredlimit;
    expminuslambda = expminusexplinpredlimit;
    }
  else
    {
    lambda = exp(*linpred);
    expminuslambda = exp(-lambda);
    }

  double denom = 1-(*worktransformlin[0])+(*worktransformlin[0])*expminuslambda;

  modify_worklin();

  if (*response==0)
    return log(denom);
  else
    return (*response)*(*linpred)-lambda;

  }


void DISTR_zip_cloglog_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of pi equation
  // *worktransformlin[0] = 1-pi = 1-exp(-exp(eta_pi));

  if (counter==0)
    {
    set_worklin();
    }

  double lambda;
  double expminuslambda;

  if (*linpred <= linpredlimit)
    {
    lambda  = explinpredlimit;
    expminuslambda = expminusexplinpredlimit;
    }
  else
    {
    lambda = exp(*linpred);
    expminuslambda = exp(-lambda);
    }


  double pi = 1-(*worktransformlin[0]);
  double denom = pi+(*worktransformlin[0])*expminuslambda;

  double nu = (*response) - lambda;
  if (*response == 0)
    nu += pi*lambda/denom;

  *workingweight = (lambda* (*worktransformlin[0])*(denom-expminuslambda*lambda*pi))/denom;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    if (*response==0)
      {
      like += log(denom);
      }
    else // response > 0
      {
      like += (*response)*(*linpred)-lambda;
      }

    }

  modify_worklin();

  }


void DISTR_zip_cloglog_mu::compute_mu_mult(vector<double *> linpred,double * mu)
  {
  double el = exp(-exp(*linpred[0]));

  *mu = (1-el)*exp(*linpred[1]);
  }


void DISTR_zip_cloglog_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (lambda): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_zip_cloglog_mu::update_end(void)
  {
  DISTR_gamlss::update_end();
  }


//------------------------------------------------------------------------------
//----------------------------- CLASS DISTR_gamlss -----------------------------
//------------------------------------------------------------------------------


DISTR_gamlss::DISTR_gamlss(GENERAL_OPTIONS * o, const datamatrix & r,
                           unsigned nrdistr,
                           const datamatrix & w)
  : DISTR(o,r,w)

  {

  predict_mult = true;

  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  counter = 0;

  helpmat1 = datamatrix(nrobs,1,1);

  worklin = vector<double*>(nrdistr);
  worktransformlin = vector<double*>(nrdistr);

  updateIWLS = true;

  }


const DISTR_gamlss & DISTR_gamlss::operator=(
const DISTR_gamlss & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  counter = nd.counter;
  worklin = nd.worklin;
  worktransformlin = nd.worktransformlin;
  distrp = nd.distrp;
  return *this;
  }


DISTR_gamlss::DISTR_gamlss(const DISTR_gamlss & nd)
   : DISTR(DISTR(nd))
  {
  counter = nd.counter;
  worklin = nd.worklin;
  worktransformlin = nd.worktransformlin;
  distrp = nd.distrp;
  }


void DISTR_gamlss::outoptions(void)
  {
  DISTR::outoptions();
  }


void DISTR_gamlss::set_worklin(void)
  {

  unsigned i;
  for (i=0;i<worklin.size();i++)
    {

    if (distrp[i]->linpred_current==1)
      worklin[i] = distrp[i]->linearpred1.getV();
    else
      worklin[i] = distrp[i]->linearpred2.getV();

    worktransformlin[i] = distrp[i]->helpmat1.getV();
    }

  }


void DISTR_gamlss::modify_worklin(void)
  {

  if (counter<nrobs-1)
    {
    counter++;
    unsigned i;
    for (i=0;i<worklin.size();i++)
      {
      worklin[i]++;
      worktransformlin[i]++;
      }
    }
  else
    {
    counter=0;
    }

  }


double DISTR_gamlss::get_intercept_start(void)
  {
  return 0;
  }


double DISTR_gamlss::loglikelihood(double * response, double * linpred,
                                         double * weight)
  {
  return loglikelihood_weightsone(response,linpred);
  }


double DISTR_gamlss::loglikelihood_weightsone(double * response,
                                                    double * linpred)
  {
  return 0;
  }


void DISTR_gamlss::compute_mu_mult(vector<double *> linpred,double * mu)
  {
  }


void DISTR_gamlss::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<double> scale) const
  {
  //  *deviance = -2*l;
  }


void DISTR_gamlss::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  if (counter==0)
    set_worklin();

//  modify_worklin();

  }


void DISTR_gamlss::posteriormode_end(void)
  {
  update_end();
  }


void DISTR_gamlss::update_end(void)
  {

//  ofstream out("d:\\_sicher\\papzip\\results\\etamu.raw");

  // helpmat1 stores exp(linpred)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    if (*worklin <= linpredlimit)
      *pmu  = explinpredlimit;
    else
      *pmu = exp(*worklin);
//    out << *worklin << "  " << *pmu << endl;
    }

  }

//------------------------------------------------------------------------------
//------------------------- CLASS DISTR_negbinzip_mu ---------------------------
//------------------------------------------------------------------------------


DISTR_negbinzip_mu::DISTR_negbinzip_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                                 const datamatrix & w)
  : DISTR(o,r,w)

  {

  predict_mult = true;

  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  counter = 0;

  helpmat1 = datamatrix(nrobs,1,1);

  family = "Zero_Inflated_Negative_Binomial - mu";
  updateIWLS = true;
  }


const DISTR_negbinzip_mu & DISTR_negbinzip_mu::operator=(
const DISTR_negbinzip_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  worklinpi = nd.worklinpi;
  workexplinpi = nd.workexplinpi;
  workonempi = nd.workonempi;
  worklindelta = nd.worklindelta;
  workexplindelta = nd.workexplindelta;
  distrpi = nd.distrpi;
  distrdelta = nd.distrdelta;
  counter = nd.counter;
  return *this;
  }


DISTR_negbinzip_mu::DISTR_negbinzip_mu(const DISTR_negbinzip_mu & nd)
   : DISTR(DISTR(nd))
  {
  worklinpi = nd.worklinpi;
  workonempi = nd.workonempi;
  workexplinpi = nd.workexplinpi;
  worklindelta = nd.worklindelta;
  workexplindelta = nd.workexplindelta;
  distrpi = nd.distrpi;
  distrdelta = nd.distrdelta;
  counter = nd.counter;
  }


void DISTR_negbinzip_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): exponential\n");
  optionsp->out("  Response function (pi): logistic distribution function\n");
  optionsp->out("  Response function (delta): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_negbinzip_mu::set_worklinpidelta(void)
  {

  if (distrpi->linpred_current==1)
    worklinpi = distrpi->linearpred1.getV();
  else
    worklinpi = distrpi->linearpred2.getV();

  workonempi = distrpi->helpmat1.getV();
  workexplinpi = distrpi->helpmat2.getV();

  if (distrdelta->linpred_current==1)
    worklindelta = distrdelta->linearpred1.getV();
  else
    worklindelta = distrdelta->linearpred2.getV();

  workexplindelta = distrdelta->helpmat1.getV();

  }


void DISTR_negbinzip_mu::modify_worklinpidelta(void)
  {

  if (counter<nrobs-1)
    {
    counter++;
    worklinpi++;
    workonempi++;
    workexplinpi++;
    worklindelta++;
    workexplindelta++;
    }
  else
    {
    counter=0;
    }

  }


double DISTR_negbinzip_mu::get_intercept_start(void)
  {
  return log(response.mean(0));
  }


double DISTR_negbinzip_mu::loglikelihood(double * response, double * linpred,
                                         double * weight)
  {
  return loglikelihood_weightsone(response,linpred);
  }


double DISTR_negbinzip_mu::loglikelihood_weightsone(double * response,
                                                    double * linpred)
  {

  if (counter==0)
    set_worklinpidelta();

  double mu;

  if (*linpred <= linpredlimit)
    mu  = explinpredlimit;
  else
    mu = exp(*linpred);

  double deltaplusmu = (*workexplindelta)+mu;
  double pot = pow((*workexplindelta)/deltaplusmu,(*workexplindelta));

  double l;

  if (*response==0)
    {
    l = log((*workexplinpi)+pot);
    }
  else // response > 0
    {
    double deltaplusy = (*response)+(*workexplindelta);
    l = (*response)*(*linpred) - deltaplusy * log((*workexplindelta)+mu);
    }

  modify_worklinpidelta();

  return l;

  }


void DISTR_negbinzip_mu::compute_mu_mult(vector<double *> linpred,double * mu)
  {

  double lambda;
  if (*linpred[2] <= linpredlimit)
    lambda  = explinpredlimit;
  else
    lambda = exp(*linpred[2]);

  double explinpi;
  if (*linpred[1] > linpredlimit)
    explinpi= exp(*linpred[1]);
  else
    explinpi = explinpredlimit;

  *mu = 1/(1+explinpi)*lambda;

  }


void DISTR_negbinzip_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<double> scale) const
  {

  double mu;
  if (*linpred[2] <= linpredlimit)
    mu  = explinpredlimit;
  else
    mu = exp(*linpred[2]);

  double explinpi;
  if (*linpred[1] > linpredlimit)
    explinpi= exp(*linpred[1]);
  else
    explinpi = explinpredlimit;

  double explindelta;
  if (*linpred[0] > linpredlimit)
    explindelta= exp(*linpred[0]);
  else
    explindelta = explinpredlimit;

  double l= -log(1+explinpi);

  if (*response[2]==0)
    {
    double pot = pow(explindelta/(explindelta+mu),explindelta);
    l += log(explinpi+pot);
    }
  else // response > 0
    {
    double help1 = *response[2]+explindelta;
    double help2 = *response[2]+1;
    l +=    randnumbers::lngamma_exact(help1)
          - randnumbers::lngamma_exact(help2)
          - randnumbers::lngamma_exact(explindelta)
          + explindelta*(*linpred[0])
          + (*response[2])*(*linpred[2])
          - (explindelta+(*response[2]))*log(explindelta+mu);
    }

  *deviance = -2*l;
  }


void DISTR_negbinzip_mu::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  if (counter==0)
    set_worklinpidelta();

  double mu;

  if (*linpred <= linpredlimit)
    mu  = explinpredlimit;
  else
    mu = exp(*linpred);

  double pi = 1-(*workonempi);

  double deltaplusmu = (*workexplindelta)+mu;

  double pot = pow((*workexplindelta)/deltaplusmu,(*workexplindelta));
  double denom = pi+(*workonempi)*pot;
  double denomfull = denom*deltaplusmu;
  double denomfull2 = denomfull*deltaplusmu;

  double nu = ((*response)*(*workexplindelta) - (*workexplindelta)*mu)/
               deltaplusmu;
  if (*response == 0)
    nu += (pi*(*workexplindelta)*mu ) /denomfull;

  *workingweight = (*workexplindelta)*mu*(*workonempi)/deltaplusmu -
                    (pi* (*workonempi) * pow((*workexplindelta),2) * pow(mu,2) * pot)  /denomfull2;

  *workingresponse = *linpred + nu/(*workingweight);


  if (compute_like)
    {

    if (*response==0)
      {
      like += log((*workexplinpi)+pot);
      }
    else // response > 0
      {
      double deltaplusy = (*response)+(*workexplindelta);
      like += (*response)*(*linpred) - deltaplusy * log((*workexplindelta)+mu);
      }

    }

  modify_worklinpidelta();

  }


void DISTR_negbinzip_mu::posteriormode_end(void)
  {
  update_end();
  }


void DISTR_negbinzip_mu::update_end(void)
  {

  // helpmat1 stores mu, i.e. exp(linpred)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    if (*worklin <= linpredlimit)
      *pmu  = explinpredlimit;
    else
      *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS DISTR_negbinzip_pi ---------------------------
//------------------------------------------------------------------------------


DISTR_negbinzip_pi::DISTR_negbinzip_pi(GENERAL_OPTIONS * o, const datamatrix & r,
                                 const datamatrix & w)
  : DISTR(o,r,w)

  {

  predict_mult = true;

  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  counter = 0;
  helpmat1=datamatrix(nrobs,1,0.5);
  helpmat2=datamatrix(nrobs,1,1.0);

  family = "Zero_Inflated_Negative_Binomial - pi";
  updateIWLS = true;
  }


const DISTR_negbinzip_pi & DISTR_negbinzip_pi::operator=(
const DISTR_negbinzip_pi & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  worklinmu = nd.worklinmu;
  workexplinmu = nd.workexplinmu;
  worklindelta = nd.worklindelta;
  workexplindelta = nd.workexplindelta;
  distrmu = nd.distrmu;
  distrdelta = nd.distrdelta;
  counter = nd.counter;
  return *this;
  }


DISTR_negbinzip_pi::DISTR_negbinzip_pi(const DISTR_negbinzip_pi & nd)
   : DISTR(DISTR(nd))
  {
  worklinmu = nd.worklinmu;
  workexplinmu = nd.workexplinmu;
  worklindelta = nd.worklindelta;
  workexplindelta = nd.workexplindelta;
  distrmu = nd.distrmu;
  distrdelta = nd.distrdelta;
  counter = nd.counter;
  }


void DISTR_negbinzip_pi::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): exponential\n");
  optionsp->out("  Response function (pi): logistic distribution function\n");
  optionsp->out("  Response function (delta): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_negbinzip_pi::set_worklinmudelta(void)
  {
  if (distrmu->linpred_current==1)
    worklinmu = distrmu->linearpred1.getV();
  else
    worklinmu = distrmu->linearpred2.getV();

  workexplinmu = distrmu->helpmat1.getV();


  if (distrdelta->linpred_current==1)
    worklindelta = distrdelta->linearpred1.getV();
  else
    worklindelta = distrdelta->linearpred2.getV();

  workexplindelta = distrdelta->helpmat1.getV();
  }


void DISTR_negbinzip_pi::modify_worklinmudelta(void)
  {
  if (counter<nrobs-1)
    {
    counter++;
    worklinmu++;
    workexplinmu++;
    worklindelta++;
    workexplindelta++;
    }
  else
    {
    counter=0;
    }
  }


double DISTR_negbinzip_pi::get_intercept_start(void)
  {
  unsigned i;
  double * responsep = response.getV();
  double m = 0;
  for (i=0;i<nrobs;i++,responsep++)
    {
    if (*responsep==0)
      m += 1;
    }

  m /= nrobs;

  return log(m/(1-m));
  }


double DISTR_negbinzip_pi::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {
  return loglikelihood_weightsone(response,linpred);
  }


double DISTR_negbinzip_pi::loglikelihood_weightsone(
                                  double * response, double * linpred)
  {

  if (counter==0)
    set_worklinmudelta();

  double explinpredpi;
  if (*linpred <= linpredlimit)
    explinpredpi  = explinpredlimit;
  else
    explinpredpi = exp(*linpred);

  double l = -log(1+explinpredpi);

  if (*response==0)
    {
    double deltaplusmu = (*workexplindelta)+(*workexplinmu);
    double pot = pow((*workexplindelta)/deltaplusmu,(*workexplindelta));

    l += log(explinpredpi+pot);
    }

  modify_worklinmudelta();

  return l;

  }


void DISTR_negbinzip_pi::compute_iwls_wweightschange_weightsone(
                                         double * response,
                                         double * linpred,
                                         double * workingweight,
                                         double * workingresponse,
                                         double & like,
                                         const bool & compute_like)
  {

  if (counter==0)
    set_worklinmudelta();

  double explinpredpi;
  if (*linpred <= linpredlimit)
    explinpredpi  = explinpredlimit;
  else
    explinpredpi = exp(*linpred);

  double oneminuspi = 0.001+0.998/(1+explinpredpi);
  double pi = 1-oneminuspi;

  double deltaplusmu = (*workexplindelta)+(*workexplinmu);

  double pot = pow((*workexplindelta)/deltaplusmu,(*workexplindelta));
  double denom = pi+oneminuspi*pot;

  double nu = -pi;
  if (*response == 0)
    nu +=  pi/denom;

  *workingweight = (pow(pi,2)*oneminuspi*(1-pot)) /denom;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    like -= log(1+explinpredpi);

    if (*response==0)
      like += log(explinpredpi+pot);
    }

  modify_worklinmudelta();

  }


void DISTR_negbinzip_pi::posteriormode_end(void)
  {
  update_end();
  }


void DISTR_negbinzip_pi::update_end(void)
  {

  // helpmat1 stores 1-pi
  // helpmat2 stores exp(eta_pi)

//  ofstream out("d:\\_sicher\\papzip\\results\\etapi.raw");

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * ppi = helpmat1.getV();
  double * workexplin = helpmat2.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,ppi++,worklin++,workexplin++)
    {
//    out << *worklin << endl;
    if (*worklin > linpredlimit)
      *workexplin = exp(*worklin);
    else
      *workexplin = explinpredlimit;
    *ppi = 0.001+0.998/(1+(*workexplin));
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS DISTR_negbinzip_delta ------------------------
//------------------------------------------------------------------------------


DISTR_negbinzip_delta::DISTR_negbinzip_delta(GENERAL_OPTIONS * o,
                                    const datamatrix & r,
                                    double & stpsum, int & strmax,
                                    const datamatrix & w)
  : DISTR(o,r,w)

  {

  responsemax = response.max(0);

  predict_mult = true;

  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  counter = 0;
  helpmat1=datamatrix(nrobs,1,1);

  family = "Zero_Inflated_Negative_Binomial - delta";
  updateIWLS = true;

  stopsum = stpsum;
  stoprmax = strmax;
  if (stoprmax < responsemax)
    stoprmax = responsemax;

  }


const DISTR_negbinzip_delta & DISTR_negbinzip_delta::operator=(
const DISTR_negbinzip_delta & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));

  responsemax = nd.responsemax;

  worklinmu = nd.worklinmu;
  workexplinmu = nd.workexplinmu;
  worklinpi = nd.worklinpi;
  workonempi = nd.workonempi;
  workexplinpi = nd.workexplinpi;;

  distrmu = nd.distrmu;
  distrpi = nd.distrpi;
  counter = nd.counter;

  delta = nd.delta;
  delta2 = nd.delta2;
  deltay = nd.deltay;
  dig_deltay = nd.dig_deltay;
  dig_delta = nd.dig_delta;
  trig_deltay = nd.trig_deltay;
  trig_delta = nd.trig_deltay;

  deltamu = nd.deltamu;
  delta_div_deltamu = nd.delta_div_deltamu;
  log_delta_div_deltamu = nd.log_delta_div_deltamu;
  mu_m_y_div_delta_m_mu = nd.mu_m_y_div_delta_m_mu;
  pi = nd.pi;
  mu_div_deltamu = nd.mu_div_deltamu;
  pot = nd.pot;
  denom = nd.denom;
  sum = nd.sum;

  log_one_explinpi = nd.log_one_explinpi;
  log_explinpi_pot = nd.log_explinpi_pot;
  lng_delta = nd.lng_delta;
  delta_linpred = nd.delta_linpred;
  log_delta_mu = nd.log_delta_mu;

  E_dig_y_delta = nd.E_dig_y_delta;
  E_trig_y_delta = nd.E_trig_y_delta;

  stopsum = nd.stopsum;
  stoprmax = nd.stoprmax;

  return *this;
  }


DISTR_negbinzip_delta::DISTR_negbinzip_delta(const DISTR_negbinzip_delta & nd)
   : DISTR(DISTR(nd))
  {

  responsemax = nd.responsemax;

  worklinmu = nd.worklinmu;
  workexplinmu = nd.workexplinmu;
  worklinpi = nd.worklinpi;
  workonempi = nd.workonempi;
  workexplinpi = nd.workexplinpi;;

  distrmu = nd.distrmu;
  distrpi = nd.distrpi;
  counter = nd.counter;

  delta = nd.delta;
  delta2 = nd.delta2;
  deltay = nd.deltay;
  dig_deltay = nd.dig_deltay;
  dig_delta = nd.dig_delta;
  trig_deltay = nd.trig_deltay;
  trig_delta = nd.trig_deltay;

  deltamu = nd.deltamu;
  delta_div_deltamu = nd.delta_div_deltamu;
  log_delta_div_deltamu = nd.log_delta_div_deltamu;
  mu_m_y_div_delta_m_mu = nd.mu_m_y_div_delta_m_mu;
  pi = nd.pi;
  mu_div_deltamu = nd.mu_div_deltamu;
  pot = nd.pot;
  denom = nd.denom;
  sum = nd.sum;

  log_one_explinpi = nd.log_one_explinpi;
  log_explinpi_pot = nd.log_explinpi_pot;
  lng_delta = nd.lng_delta;
  delta_linpred = nd.delta_linpred;
  log_delta_mu = nd.log_delta_mu;

  E_dig_y_delta = nd.E_dig_y_delta;
  E_trig_y_delta = nd.E_trig_y_delta;

  stopsum = nd.stopsum;
  stoprmax = nd.stoprmax;

  }


void DISTR_negbinzip_delta::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): exponential\n");
  optionsp->out("  Response function (pi): logistic distribution function\n");
  optionsp->out("  Response function (delta): exponential\n");
  optionsp->out("  Stop criteria for approximating expected values\n");
  optionsp->out("  in working weights of delta equation:\n");
  optionsp->out("    Cumulative probabilities:"  + ST::doubletostring(stopsum) +  "\n");
  optionsp->out("    Maximum responses:"  + ST::inttostring(stoprmax) +  "\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_negbinzip_delta::set_worklinmupi(void)
  {

  if (distrmu->linpred_current==1)
    worklinmu = distrmu->linearpred1.getV();
  else
    worklinmu = distrmu->linearpred2.getV();

  workexplinmu = distrmu->helpmat1.getV();


  if (distrpi->linpred_current==1)
    worklinpi = distrpi->linearpred1.getV();
  else
    worklinpi = distrpi->linearpred2.getV();

  workonempi = distrpi->helpmat1.getV();
  workexplinpi = distrpi->helpmat2.getV();

  }


void DISTR_negbinzip_delta::modify_worklinmupi(void)
  {

  if (counter<nrobs-1)
    {
    counter++;
    worklinmu++;
    workexplinmu++;
    worklinpi++;
    workonempi++;
    workexplinpi++;
    }
  else
    {
    counter=0;
    }

  }

double DISTR_negbinzip_delta::get_intercept_start(void)
  {
  return log(response.mean(0));
  }


double DISTR_negbinzip_delta::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {
  return loglikelihood_weightsone(response,weight);
  }


double DISTR_negbinzip_delta::loglikelihood_weightsone(double * response,
                                                       double * linpred)
  {

  if (counter==0)
    set_worklinmupi();

  double delta;
  if (*linpred <= linpredlimit)
    delta  = explinpredlimit;
  else
    delta = exp(*linpred);

  double deltay = delta+(*response);

  double deltamu = delta + (*workexplinmu);
  double delta_div_deltamu = delta/deltamu;
  double pot = pow(delta_div_deltamu,delta);

  double l;
  if (*response==0)
    {
    l = log((*workexplinpi)+pot);
    }
  else // response > 0
    {

    l =    randnumbers::lngamma_exact(deltay)
         - randnumbers::lngamma_exact(delta)
         + delta*(*linpred)
         - deltay * log(delta+(*workexplinmu));
    }

  modify_worklinmupi();

  return l;

  }


void DISTR_negbinzip_delta::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  if (counter==0)
    set_worklinmupi();

  if (*linpred <= linpredlimit)
    delta  = explinpredlimit;
  else
    delta = exp(*linpred);

  delta2 = delta*delta;
  deltay = delta+(*response);
  dig_deltay = randnumbers::digamma_exact(deltay);
  dig_delta = randnumbers::digamma_exact(delta);
  trig_deltay = randnumbers::trigamma_exact(deltay);
  trig_delta = randnumbers::trigamma_exact(delta);

  deltamu = delta + (*workexplinmu);
  delta_div_deltamu = delta/deltamu;
  log_delta_div_deltamu = log(delta_div_deltamu);
  mu_m_y_div_delta_m_mu = ((*workexplinmu) - (*response))/deltamu;
  pi = 1-(*workonempi);
  mu_div_deltamu = (*workexplinmu)/deltamu;
  pot = pow(delta_div_deltamu,delta);
  denom = pi+(*workonempi)*pot;
  sum = log_delta_div_deltamu+mu_div_deltamu;

  log_one_explinpi = log(1+(*workexplinpi));
  log_explinpi_pot = log((*workexplinpi)+pot);

  lng_delta = randnumbers::lngamma_exact(delta);
  delta_linpred = delta*(*linpred);
  log_delta_mu = log(delta+(*workexplinmu));

  double nu = delta*
  (dig_deltay-dig_delta + log_delta_div_deltamu + mu_m_y_div_delta_m_mu);

  if (*response == 0)
    nu -= delta*pi*(log_delta_div_deltamu + mu_div_deltamu)/denom;

  compute_expectation();

  *workingweight = -delta*(*workonempi)*(log_delta_div_deltamu+mu_div_deltamu)
                   - ((*workonempi)*pi*delta2*pot*sum*sum)/denom
                   - delta*(E_dig_y_delta - dig_delta)
                   - delta2*(E_trig_y_delta - trig_delta);
  if (*workingweight<0)
    *workingweight = 0.001;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    if (*response==0)
      {
      like += log_explinpi_pot;
      }
    else // response > 0
      {

      like +=   randnumbers::lngamma_exact(deltay)
              - lng_delta
              + delta_linpred
              - deltay * log_delta_mu;
      }

    }

  modify_worklinmupi();

  }


void DISTR_negbinzip_delta::compute_expectation(void)
  {

  int k=1;
  double k_delta;
  double kplus1;
  double psum;

  double L = exp(-log_one_explinpi + log_explinpi_pot);
         E_dig_y_delta = dig_delta*L;
         E_trig_y_delta = trig_delta*L;
  psum = L;

  while ((psum < stopsum) && (k <=stoprmax))
    {
    k_delta = k + delta;
    kplus1 = k + 1;

    L = exp(-log_one_explinpi
        +randnumbers::lngamma_exact(k_delta)
        -randnumbers::lngamma_exact(kplus1)
        -lng_delta+delta_linpred+k*(*worklinmu)-(delta+k)*log_delta_mu);

    psum += L;

    E_dig_y_delta += randnumbers::digamma_exact(k_delta)*L;

    E_trig_y_delta += randnumbers::trigamma_exact(k_delta)*L;

    k++;
    }

  }


void DISTR_negbinzip_delta::posteriormode_end(void)
  {
  update_end();
  }


void DISTR_negbinzip_delta::update_end(void)
  {


//  ofstream out("d:\\_sicher\\papzip\\results\\etadelta.raw");
  // helpmat1 stores exp(linpred)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * l = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,worklin++,l++)
    {
//    out << *worklin << endl;
    if (*worklin <= linpredlimit)
      *l  = explinpredlimit;
    else
      *l = exp(*worklin);
    }

  }



//------------------------------------------------------------------------------
//---------------------- CLASS DISTRIBUTION_ziplambda --------------------------
//------------------------------------------------------------------------------


DISTR_ziplambda::DISTR_ziplambda(GENERAL_OPTIONS * o, const datamatrix & r,
                                 const datamatrix & w)
  : DISTR(o,r,w)

  {

  predict_mult = true;

  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  counter = 0;

  family = "ZIP";
  updateIWLS = true;
  }


const DISTR_ziplambda & DISTR_ziplambda::operator=(const DISTR_ziplambda & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  distrpi = nd.distrpi;
  counter = nd.counter;
  worklinpi = nd.worklinpi;
  workonempi = nd.workonempi;
  workexplinpi = nd.workexplinpi;

  return *this;
  }


DISTR_ziplambda::DISTR_ziplambda(const DISTR_ziplambda & nd)
   : DISTR(DISTR(nd))
  {
  distrpi = nd.distrpi;
  counter = nd.counter;
  worklinpi = nd.worklinpi;
  workonempi = nd.workonempi;
  workexplinpi = nd.workexplinpi;
  }


void DISTR_ziplambda::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (lambda): exponential\n");
  optionsp->out("  Response function (pi): logistic distribution function\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTR_ziplambda::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {

  if (counter==0)
    {
    set_worklinpi();
    }

  double lambda;
  double exptildeeta;
  double expminuslambda;


  if (*linpred <= linpredlimit)
    {
    lambda  = explinpredlimit;
    expminuslambda = expminusexplinpredlimit;
    }
  else
    {
    lambda = exp(*linpred);
    expminuslambda = exp(-lambda);
    }

  if (*worklinpi > linpredlimit)
    exptildeeta= exp(*worklinpi);
  else
    exptildeeta = explinpredlimit;

  double l;

  if (*response==0)
    {
    l = -log(1+exptildeeta) + log(exptildeeta+expminuslambda);
    }
  else // response > 0
    {
    l = -log(1+exptildeeta) + (*response)*(*linpred)-lambda;
    }

  modify_worklinpi();
  return l;


  }


void DISTR_ziplambda::set_worklinpi(void)
  {
  if (distrpi->linpred_current==1)
    worklinpi = distrpi->linearpred1.getV();
  else
    worklinpi = distrpi->linearpred2.getV();

  workexplinpi = distrpi->helpmat1.getV();
  workonempi = distrpi->helpmat2.getV();

  }


void DISTR_ziplambda::modify_worklinpi(void)
  {
  if (counter<nrobs-1)
    {
    counter++;
    worklinpi++;
    workexplinpi++;
    workonempi++;
    }
  else
    {
    counter=0;
    }
  }


double DISTR_ziplambda::get_intercept_start(void)
  {
  return log(response.mean(0));
  }


double DISTR_ziplambda::loglikelihood_weightsone(
                                  double * response, double * linpred)
  {

  if (counter==0)
    {
    set_worklinpi();
    }

  double lambda;
  double expminuslambda;

  if (*linpred <= linpredlimit)
    {
    lambda  = explinpredlimit;
    expminuslambda = expminusexplinpredlimit;
    }
  else
    {
    lambda = exp(*linpred);
    expminuslambda = exp(-lambda);
    }

  double l;

  if (*response==0)
    l = -log(1+(*workexplinpi)) + log((*workexplinpi)+expminuslambda);
  else // response > 0
    l = -log(1+(*workexplinpi)) + (*response)*(*linpred)-lambda;

  modify_worklinpi();

  return l;

  }


void DISTR_ziplambda::compute_mu_mult(vector<double *> linpred,double * mu)
  {

  double lambda;
  if (*linpred[1] <= linpredlimit)
    lambda  = explinpredlimit;
  else
    lambda = exp(*linpred[1]);

  double explinpi;
  if (*linpred[0] > linpredlimit)
    explinpi= exp(*linpred[0]);
  else
    explinpi = explinpredlimit;

  *mu = 1/(1+explinpi)*lambda;

  }


void DISTR_ziplambda::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<double> scale) const
  {

  double l;

  double lambda;
  double expminuslambda;
  double explinpi;

  if (*linpred[1] <= linpredlimit)
    {
    lambda  = explinpredlimit;
    expminuslambda = expminusexplinpredlimit;
    }
  else
    {
    lambda = exp(*linpred[1]);
    expminuslambda = exp(-lambda);
    }

  if (*linpred[0] > linpredlimit)
    explinpi= exp(*linpred[0]);
  else
    explinpi = explinpredlimit;


  if (*response[1]==0)
    {
    l= -log(1+ explinpi) + log(explinpi+ expminuslambda);
    }
  else // response > 0
    {
    double help1 = *response[1]+1;
    l= -log(1+ explinpi) + (*response[1])*(*linpred[1])- lambda
       - randnumbers::lngamma_exact(help1);
    }

  *deviance = -2*l;
  }


void DISTR_ziplambda::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  if (counter==0)
    {
    set_worklinpi();
    }

  double lambda;
  double expminuslambda;


  if (*linpred <= linpredlimit)
    {
    lambda  = explinpredlimit;
    expminuslambda = expminusexplinpredlimit;
    }
  else
    {
    lambda = exp(*linpred);
    expminuslambda = exp(-lambda);
    }

  double pi = 1-(*workonempi);
  double denom = pi+(*workonempi)*expminuslambda;

  double nu = (*response) - lambda;
  if (*response == 0)
    nu += pi*lambda/denom;

  *workingweight = (lambda* (*workonempi)*(denom-expminuslambda*lambda*pi))/denom;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    if (*response==0)
      {
      like += -log(1+ (*workexplinpi)) + log((*workexplinpi)+expminuslambda);
      }
    else // response > 0
      {
      like += -log(1+(*workexplinpi)) + (*response)*(*linpred)-lambda;
      }

    }

  modify_worklinpi();

  }


void DISTR_ziplambda::posteriormode_end(void)
  {
  update_end();
  }


void DISTR_ziplambda::update_end(void)
  {

  // helpmat1 stores exp(-lambda)
  // helpmat2 stores lambda

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  if (helpmat1.rows() == 1)
    {
    helpmat1 = datamatrix(nrobs,1,0);
    helpmat2 = datamatrix(nrobs,1,0);
    }

  double * ph = helpmat1.getV();
  double * l = helpmat2.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,ph++,worklin++,l++)
    {
    if (*worklin <= linpredlimit)
      {
      *l  = explinpredlimit;
      *ph = expminusexplinpredlimit;
      }
    else
      {
      *l = exp(*worklin);
      *ph = exp(-(*l));
      }
    }

  }


//------------------------------------------------------------------------------
//------------------------ CLASS DISTRIBUTION_zippi ----------------------------
//------------------------------------------------------------------------------


DISTR_zippi::DISTR_zippi(GENERAL_OPTIONS * o, const datamatrix & r,
                                 const datamatrix & w)
  : DISTR(o,r,w)

  {

  maindistribution=false;
  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  counter = 0;
  helpmat1=datamatrix(nrobs,1,1);
  helpmat2=datamatrix(nrobs,1,0.5);


  family = "ZIP_pi";
  updateIWLS = true;
  }


const DISTR_zippi & DISTR_zippi::operator=(const DISTR_zippi & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  distrlambda = nd.distrlambda;
  counter = nd.counter;
  worklinlambda = nd.worklinlambda;
  workexpmlambda = nd.workexpmlambda;
  return *this;
  }


DISTR_zippi::DISTR_zippi(const DISTR_zippi & nd)
   : DISTR(DISTR(nd))
  {
  distrlambda = nd.distrlambda;
  counter = nd.counter;
  worklinlambda = nd.worklinlambda;
  workexpmlambda = nd.workexpmlambda;
  }


void DISTR_zippi::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (pi): logistic distribution function\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTR_zippi::get_intercept_start(void)
  {
  unsigned i;
  double * responsep = response.getV();
  double m = 0;
  for (i=0;i<nrobs;i++,responsep++)
    {
    if (*responsep==0)
      m += 1;
    }

  m /= nrobs;

  return log(m/(1-m));
  }


double DISTR_zippi::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {
  return loglikelihood_weightsone(response,linpred);;
  }



void DISTR_zippi::set_worklinlambda(void)
  {
  if (distrlambda->linpred_current==1)
    worklinlambda = distrlambda->linearpred1.getV();
  else
    worklinlambda = distrlambda->linearpred2.getV();

  workexpmlambda = distrlambda->helpmat1.getV();
  worklambda = distrlambda->helpmat2.getV();
  }


void DISTR_zippi::modify_worklinlambda(void)
  {
  if (counter<nrobs-1)
    {
    counter++;
    worklinlambda++;
    worklambda++;
    workexpmlambda++;
    }
  else
    {
    counter=0;
    }
  }


double DISTR_zippi::loglikelihood_weightsone(
                                  double * response, double * linpred)
  {

  if (counter==0)
    set_worklinlambda();

  double exptildeeta;

  if (*linpred > linpredlimit)
    exptildeeta= exp(*linpred);
  else
    exptildeeta = explinpredlimit;

  double l;

  if (*response==0)
    {
    l = -log(1+exptildeeta) + log(exptildeeta+(*workexpmlambda));
    }
  else // response > 0
    {
    l = -log(1+exptildeeta) + (*response)*(*worklinlambda)-(*worklambda);
    }

  modify_worklinlambda();

  return l;

  }


void DISTR_zippi::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  if (counter==0)
    {
    set_worklinlambda();
    }

  double exptildeeta;

  if (*linpred > linpredlimit)
    exptildeeta= exp(*linpred);
  else
    exptildeeta = explinpredlimit;

  double oneminuspi = 0.001+0.998/(1+exptildeeta);
  double pi = 1-oneminuspi;
  double denom = pi+oneminuspi* (*workexpmlambda);

  double nu = - pi;
  if (*response == 0)
    nu += pi/denom;

  *workingweight = (pi*pi*(1-(*workexpmlambda))*oneminuspi)/denom;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {
    if (*response==0)
      like += -log(1+exptildeeta) + log(exptildeeta+ (*workexpmlambda));
    else // response > 0
      like += -log(1+exptildeeta) + (*response)*(*worklinlambda)-(*worklambda);
    }

  modify_worklinlambda();

  }


void DISTR_zippi::posteriormode_end(void)
  {
  update_end();
  }


void DISTR_zippi::update_end(void)
  {

  // helpmat1 stores exp(tildeeta)
  // helpmat2 stores 1-pi

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  if (helpmat1.rows() == 1)
    {
    helpmat1 = datamatrix(nrobs,1,0);
    helpmat2 = datamatrix(nrobs,1,0);
    }

  double * ete = helpmat1.getV();
  double * wpi = helpmat2.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,ete++,worklin++,wpi++)
    {

    if (*worklin > linpredlimit)
      *ete= exp(*worklin);
    else
      *ete = explinpredlimit;

    *wpi = 0.001+0.998/(1+(*ete));
    }

  }


} // end: namespace MCMC



