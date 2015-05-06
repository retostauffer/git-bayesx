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

#include "distr_gamlss_copula.h"
//#include "gsl/gsl_randist.h"
//#include "gsl/gsl_cdf.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_gausscopula --------------------
//------------------------------------------------------------------------------
void DISTR_gausscopula::check_errors(void)
  {
  // Note: check marginal distributions & weights.
  if (errors==false)
    {
    errors=false;
    }
  }


DISTR_gausscopula::DISTR_gausscopula(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Gauss Copula - rho";

  outpredictor = true;
  outexpectation = true;
  predictor_name = "rho";
    linpredminlimit=-100;
  linpredmaxlimit=100;
  check_errors();
  }


DISTR_gausscopula::DISTR_gausscopula(const DISTR_gausscopula & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_gausscopula & DISTR_gausscopula::operator=(
                            const DISTR_gausscopula & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_gausscopula::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gausscopula::compute_param_mult(vector<double *>  linpred,double * param)
  {
// check
  *param = (*linpred[2]) / pow(1+ pow((*linpred[2]), 2), 0.5);
  }

void DISTR_gausscopula::set_worklin(void)
  {
  DISTR_gamlss::set_worklin();
  responsep = response.getV();
  response2p = response2.getV();
  }

void DISTR_gausscopula::modify_worklin(void)
  {
  DISTR_gamlss::modify_worklin();
  if (counter<nrobs)
    {
    responsep++;
    response2p++;
    }
  }

void DISTR_gausscopula::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {
   if (*weight[0] == 0)
     *deviance=0;
   else
     {
     double rho = (*linpred[0]) / pow(1 + pow((*linpred[0]), 2), 0.5);
     if (*linpred[0] <= -100)
        rho  = -0.99995;
     else if (*linpred[0] >= 100)
        rho  = 0.99995;

     double orho = 1 - pow(rho, 2);
     double phinvu = randnumbers::invPhi2(distrp[0]->cdf(*response[0]));
     double phinvv = randnumbers::invPhi2(distrp[1]->cdf(*response[1]));
     double l;

      l = - 0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho;


    *deviance = -2*l;
    }

  }

double DISTR_gausscopula::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {
  if (counter==0)
    {
    set_worklin();
    }
    double rho = (*linpred) / pow(1 + pow((*linpred), 2), 0.5);
    if (*linpred <= -100)
        rho  = -0.99995;
    else if (*linpred >= 100)
        rho  = 0.99995;

    double orho = 1 - pow(rho, 2);
    double phinvu = randnumbers::invPhi2(distrp[0]->cdf(*response));
    double phinvv = randnumbers::invPhi2(distrp[1]->cdf(*response2p));

    double l;

    l = - 0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho;

  modify_worklin();

  return l;

  }

void DISTR_gausscopula::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  if (counter==0)
    {
    set_worklin();
    }
    double hilfs = pow(1 + pow((*linpred), 2), 0.5);
    double rho = (*linpred) / hilfs;
    if (*linpred <= -100)
        rho  = -0.99995;
    else if (*linpred >= 100)
        rho  = 0.99995;

    double orho = 1 - pow(rho, 2);
    double phinvu = randnumbers::invPhi2(distrp[0]->cdf(*response));
    double phinvv = randnumbers::invPhi2(distrp[1]->cdf(*response2p));

    double nu = rho * pow(orho, 0.5) + (hilfs + rho * (*linpred)) * (phinvu * phinvv) - (*linpred) * (pow(phinvu, 2) + pow(phinvv, 2));

    *workingweight = 1 - pow(rho, 4);

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {
      like += - 0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho;
      }

  modify_worklin();

  }

void DISTR_gausscopula::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 2 * std::asin((*linpred[predstart_mumult]) / (pow(1 + pow((*linpred[predstart_mumult]), 2), 0.5))) / PI ;
  }


void DISTR_gausscopula::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (rho): \n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gausscopula::update_end(void)
  {

  // helpmat1 stores rho

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = (*worklin) / pow(1 + pow((*worklin), 2), 0.5);
    }

  }

vector<double> DISTR_gausscopula::derivative(double * linpred)
  {
  if(counter==0)
    set_worklin();
  vector<double> res;
  double F1 = distrp[0]->cdf(*responsep);
  double F2 = distrp[1]->cdf(*response2p);
// do something with F1, F2
  res.push_back(F1);
  res.push_back(F2);
  modify_worklin();
  return res;
  }

double DISTR_gausscopula::logc(double & F, int & copulapos)
  {
  double Fa;
  if(copulapos==0)
    {
    Fa = distrp[1]->cdf(*response2p);
    return logc(F, Fa);
    }
  else
    {
    Fa = distrp[0]->cdf(*responsep);
    return logc(Fa, F);
    }
  }

double DISTR_gausscopula::logc(double & F1, double & F2)
  {
  // do something;
  return 0.0;
  }

} // end: namespace MCMC
