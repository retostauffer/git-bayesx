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
    worklindelta++;
    workexplindelta++;
    }
  else
    {
    counter=0;
    }

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

  if (*linpred <= -10)
    mu  = 0.0000454;
  else
    mu = exp(*linpred);

  double deltaplusmu = (*workexplindelta)+mu;
  double pot = pow((*workexplindelta)/deltaplusmu,(*workexplindelta));

  double l;

  if (*response==0)
    {
    l = log(exp(*worklinpi)+pot);
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

  double el = exp(*linpred[1]);

  *mu = 1/(1+el)*exp(*linpred[2]);

  }


void DISTR_negbinzip_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             double * deviancesat,
                             vector<double> scale) const
  {

  double explindelta = exp(*linpred[0]);
  double explinpi = exp(*linpred[1]);
  double lambda = exp(*linpred[2]);

  double l= -log(1+explinpi);

  if (*response[2]==0)
    {
    double pot = pow(explindelta/(explindelta+lambda),explindelta);
    l += log(explinpi+pot);
    }
  else // response > 0
    {
    double help1 = *response[2]+explindelta;
    double help2 = *response[2]+1.0;
    l +=    randnumbers::lngamma_exact(help1)
          - randnumbers::lngamma_exact(help2)
          - randnumbers::lngamma_exact(explindelta)
          + explindelta*(*linpred[0])
          + (*response[2])*(*linpred[2])
          - (explindelta+(*response[2]))*log(explindelta+lambda);
    }


  *deviance = -2*l;
  *deviancesat = *deviance;

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

  if (*linpred <= -10)
    mu  = 0.0000454;
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
      double explinpredpi = exp(*worklinpi);

      like += log(explinpredpi+pot);
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

  if (helpmat1.rows() == 1)
    {
    helpmat1 = datamatrix(nrobs,1,0);
    }

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    if (*worklin <= -10)
      {
      *pmu  = 0.0000454;
      }
    else
      {
      *pmu = exp(*worklin);
      }

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

  double explinpredpi = exp(*linpred);

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

  double explinpredpi = exp(*linpred);
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


  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  if (helpmat1.rows() == 1)
    {
    helpmat1 = datamatrix(nrobs,1,0);
    }

  double * ppi = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,ppi++,worklin++)
    {
    *ppi = 0.001+0.998/(1+exp(*worklin));
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

  explinpi = nd.explinpi;
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

  explinpi = nd.explinpi;
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
    }
  else
    {
    counter=0;
    }

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

  double explinpi = exp(*worklinpi);

  double delta = exp(*linpred);
  double deltay = delta+(*response);

  double deltamu = delta + (*workexplinmu);
  double delta_div_deltamu = delta/deltamu;
  double pot = pow(delta_div_deltamu,delta);

  double l;
  if (*response==0)
    {
    l = log(explinpi+pot);
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

  explinpi = exp(*worklinpi);
  log_one_explinpi = log(1+(explinpi));
  log_explinpi_pot = log(explinpi+pot);

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
        +randnumbers::lngamma(k_delta)
        -randnumbers::lngamma(kplus1)
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

  // helpmat1 stores exp(linpred)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  if (helpmat1.rows() == 1)
    {
    helpmat1 = datamatrix(nrobs,1,0);
    }

  double * l = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,worklin++,l++)
    {
    if (*worklin <= -10)
      {
      *l  = 0.0000454;
      }
    else
      {
      *l = exp(*worklin);
      }
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

  double lambda = exp(*linpred);
  double exptildeeta = exp(*worklinpi);

  double l;

  if (*response==0)
    {
    l = -log(1+exptildeeta) + log(exptildeeta+exp(-lambda));

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


double DISTR_ziplambda::loglikelihood_weightsone(
                                  double * response, double * linpred)
  {

  if (counter==0)
    {
    set_worklinpi();
    }

  double lambda = exp(*linpred);

  double l;

  if (*response==0)
    {
    l = -log(1+(*workexplinpi)) + log((*workexplinpi)+exp(-lambda));

    }
  else // response > 0
    {
    l = -log(1+(*workexplinpi)) + (*response)*(*linpred)-lambda;
    }

  modify_worklinpi();

  return l;

  }


void DISTR_ziplambda::compute_mu_mult(vector<double *> linpred,double * mu)
  {

  double el = exp(*linpred[0]);

  *mu = 1/(1+el)*exp(*linpred[1]);

  }


void DISTR_ziplambda::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             double * deviancesat,
                             vector<double> scale) const
  {

  double l;
  double explinpi = exp(*linpred[0]);
  double lambda = exp(*linpred[1]);

  if (*response[1]==0)
    {
    l= -log(1+ explinpi) + log(explinpi+ exp(-lambda));
    }
  else // response > 0
    {
    double help = *response[1]+1.0;
    l= -log(1+ explinpi) + (*response[1])*(*linpred[1])- lambda
       + randnumbers::lngamma_exact(help);
    }

  *deviance = -2*l;
  *deviancesat = *deviance;

  }


double DISTR_ziplambda::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {

  if (counter==0)
    {
    set_worklinpi();
    }


  double lambda;
  double expminuslambda;

  if (*linpred <= -10)
    {
    lambda  = 0.0000454;
    expminuslambda = 0.9999546;
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

  *workingweight = (lambda*(*workonempi)*(denom-expminuslambda*lambda*pi))/denom;

  *workingresponse = *linpred + nu/(*workingweight);

  double l=0;

  if (like)
    {

    if (*response==0)
      {
      l= -log(1+ (*workexplinpi)) + log((*workexplinpi)+expminuslambda);
      }
    else // response > 0
      {
      l= -log(1+(*workexplinpi)) + (*response)*(*linpred)-lambda;
      }

    }

  modify_worklinpi();

  return l;
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

  if (*linpred <= -10)
    {
    lambda  = 0.0000454;
    expminuslambda = 0.9999546;
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
    if (*worklin <= -10)
      {
      *l  = 0.0000454;
      *ph = 0.9999546;
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


double DISTR_zippi::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {

  if (counter==0)
    set_worklinlambda();

  double exptildeeta = exp(*linpred);

  double l;

  if (*response==0)
    {
    l = -log(1+exptildeeta) + log(exptildeeta+ (*workexpmlambda));
    }
  else // response > 0
    {
    l = -log(1+exptildeeta) + (*response)*(*worklinlambda)-(*worklambda);
    }

  modify_worklinlambda();

  return l;

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

  double exptildeeta = exp(*linpred);

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


void DISTR_zippi::compute_mu_mult(vector<double *> linpred,double * mu)
  {


  }


void DISTR_zippi::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             double * deviancesat,
                             vector<double> scale) const
  {


  }




double DISTR_zippi::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {

  if (counter==0)
    {
    set_worklinlambda();
    }

  double oneminuspi = 0.001+0.998/(1+exp(*linpred));
  double pi = 1-oneminuspi;
  double denom = pi+oneminuspi*(*workexpmlambda);
  double nu = - pi;

  if (*response == 0)
    nu += pi/denom;

  *workingweight = (pi*pi*(1-(*workexpmlambda))*oneminuspi)/denom;

  *workingresponse = *linpred + nu/(*workingweight);

  double l=0;
  if (like)
    {

    double exptildeeta = exp(*linpred);

    if (*response==0)
      {
      l = -log(1+exptildeeta) + log(exptildeeta+ (*workexpmlambda));

      }
    else // response > 0
      {
      l= -log(1+exptildeeta) + (*response)*(*worklinlambda)-(*worklambda);
      }

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

  double oneminuspi = 0.001+0.998/(1+exp(*linpred));
  double pi = 1-oneminuspi;
  double denom = pi+oneminuspi* (*workexpmlambda);

  double nu = - pi;
  if (*response == 0)
    nu += pi/denom;

  *workingweight = (pi*pi*(1-(*workexpmlambda))*oneminuspi)/denom;

//  if (*workingweight < 0.000001)
//    *workingweight = 0.000001;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    double exptildeeta = exp(*linpred);

    if (*response==0)
      {
      like += -log(1+exptildeeta) + log(exptildeeta+ (*workexpmlambda));
      }
    else // response > 0
      {
      like += -log(1+exptildeeta) + (*response)*(*worklinlambda)-(*worklambda);
      }

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
    *ete = exp(*worklin);
    *wpi = 0.001+0.998/(1+(*ete));
    }

  }


} // end: namespace MCMC



