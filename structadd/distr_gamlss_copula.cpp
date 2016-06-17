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

//------------------------------------------------------------------------
//------------------------- CLASS: DISTR_copula_basis --------------------
//------------------------------------------------------------------------
void DISTR_copula_basis::check_errors(void)
  {
  if (errors==false)
    {
    errors=false;
    }
  }


DISTR_copula_basis::DISTR_copula_basis(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  check_errors();
  }


DISTR_copula_basis::DISTR_copula_basis(const DISTR_copula_basis & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  response1 = nd.response1;
  response1p = nd.response1p;
  linpredp = nd.linpredp;
  }


const DISTR_copula_basis & DISTR_copula_basis::operator=(
                            const DISTR_copula_basis & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  response1 = nd.response1;
  response1p = nd.response1p;
  linpredp = nd.linpredp;
  return *this;
  }


double DISTR_copula_basis::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_copula_basis::compute_param_mult(vector<double *>  linpred,double * param)
  {
  }

void DISTR_copula_basis::set_worklin(void)
  {
  response1p = response1.getV();
  response2p = response2.getV();
  }

void DISTR_copula_basis::modify_worklin(void)
  {
  if (counter<nrobs-1)
    {
    counter++;
    }
  else
    {
    counter=0;
    }
  if (counter<nrobs)
    {
    response1p++;
    response2p++;
    }
  }

void DISTR_copula_basis::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {
  }

double DISTR_copula_basis::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {
  return 0.0;
  }

void DISTR_copula_basis::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  }

void DISTR_copula_basis::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  }


void DISTR_copula_basis::outoptions(void)
  {
  DISTR::outoptions();
  }


void DISTR_copula_basis::update_end(void)
  {
  }

vector<double> DISTR_copula_basis::derivative(double & F1, double & F2, double * linpred)
  {
  vector<double> res;
  return res;
  }

vector<double> DISTR_copula_basis::logc(double & F, int & copulapos, const bool & deriv)
  {
  vector<double> res;
  if (counter==0)
    {
    if (linpred_current==1)
      linpredp = linearpred1.getV();
    else
      linpredp = linearpred2.getV();

    response1p = response1.getV();
    response2p = response2.getV();
    }
  double Fa;
 // double eta = (*linpred);
  if(copulapos==0)
    {
    //implement Fa
    //cdf virtual in distr hat nur ein Argument!
    Fa = distrp[1]->cdf(*response1p,true);
    if(optionsp->rotation == 90)
      Fa = 1-Fa;
    else if(optionsp->rotation == 180)
      {
      Fa = 1-Fa;
      F = 1-F;
      }
    else if(optionsp->rotation == 270)
      F = 1-F;

    res.push_back(logc(Fa, F, linpredp));
    }
  else
    {
    // implement Fa
    Fa = distrp[0]->cdf(*response2p,true);
    if(optionsp->rotation == 90)
      F = 1-F;
    else if(optionsp->rotation == 180)
      {
      Fa = 1-Fa;
      F = 1-F;
      }
    else if(optionsp->rotation == 270)
      Fa = 1-Fa;

    res.push_back(logc(F, Fa, linpredp));
    }

  if(deriv)
    {
    vector<double> derivs = derivative(F, Fa, linpredp);

    if(copulapos==0)
      {
      if((optionsp->rotation == 180) || (optionsp->rotation == 270))
        derivs[0] = -derivs[0];
      }
    else
      {
      if((optionsp->rotation == 180) || (optionsp->rotation == 90))
        derivs[0] = -derivs[0];
      }
    res.push_back(derivs[0]);
    res.push_back(derivs[1]);
    }
  linpredp++;
  response1p++;
  response2p++;
  if (counter<nrobs-1)
    counter++;
  else
    counter=0;

  return res;
  }

double DISTR_copula_basis::logc(double & F1, double & F2, double * linpred)
  {
  return 0.0;
  }

double DISTR_copula_basis::condfc(double & x, double & linpred_F, double & y, int & copulapos)
  {
  if (counter==0)
    {
    if (linpred_current==1)
      linpredp = linearpred1.getV();
    else
      linpredp = linearpred2.getV();

    response1p = response1.getV();
    response2p = response2.getV();
    weightp = weight.getV();
    }

  double Fa;
  if(copulapos==0)
    Fa = distrp[1]->cdf(*response1p,true);
  else
    Fa = distrp[0]->cdf(*response2p,true);

  double res;

  if(*weightp==0)
    {
    if(y==0)
      {
      res = trunc_normal2(-20,0,linpred_F,1);
      }
    else
      {
      res = trunc_normal2(0,20,linpred_F,1);
      }
    }
  else
    {
    res = condfc(x, linpred_F, y, Fa, linpredp);
    }

  linpredp++;
  response1p++;
  response2p++;
  weightp++;
  if (counter<nrobs-1)
    counter++;
  else
    counter=0;

  return res;
  }

double DISTR_copula_basis::condfc(double & x, double & linpred_F, double & y, double & F2, double * linpred)
  {
  return 0.0;
  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_gausscopula ---------------------------
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
  : DISTR_copula_basis(o,r,w)
  {
  family = "Gauss Copula - rho";

  outpredictor = true;
  outexpectation = true;
  predictor_name = "rho";
  linpredminlimit=-15;
  linpredmaxlimit=15;
  check_errors();
  }


DISTR_gausscopula::DISTR_gausscopula(const DISTR_gausscopula & nd)
   : DISTR_copula_basis(DISTR_copula_basis(nd))
  {
  }


const DISTR_gausscopula & DISTR_gausscopula::operator=(
                            const DISTR_gausscopula & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_copula_basis::operator=(DISTR_copula_basis(nd));
  return *this;
  }


double DISTR_gausscopula::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gausscopula::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[(linpred.size()-1)]) / pow(1+ pow((*linpred[(linpred.size()-1)]), 2), 0.5);
  }

void DISTR_gausscopula::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {
   if ((*weight[0] == 0) || (*weight[weight.size()-2] == 0))
     *deviance=0;
   else
     {
     double rho = (*linpred[(linpred.size()-1)]) / pow(1 + pow((*linpred[(linpred.size()-1)]), 2), 0.5);
     if (*linpred[(linpred.size()-1)] <= -100)
        rho  = -0.99995;
     else if (*linpred[(linpred.size()-1)] >= 100)
        rho  = 0.99995;

     double orho = 1 - pow(rho, 2);

    int s1 = distrp[1]->distrp.size();
    int s2 = distrp[0]->distrp.size();

    vector<double*> linpredvec1;
    vector<double*> responsevec1;
    vector<double*> weightvec1;
    vector<double*> linpredvec2;
    vector<double*> responsevec2;
    vector<double*> weightvec2;


    int j;
    for (j=0;j<(s2+1);j++)
      {
      linpredvec2.push_back(linpred[j]);
      weightvec2.push_back(weight[j]);
      responsevec2.push_back(response[j]);
      }
    int k;
    for (k=0;k<(s1+1);k++)
      {
      linpredvec1.push_back(linpred[s2+1+k]);
      weightvec1.push_back(weight[s2+1+k]);
      responsevec1.push_back(response[s2+1+k]);
      }

    double d1;
    double d2;
    distrp[0]->compute_deviance_mult(responsevec2,weightvec2,linpredvec2,&d2,aux);
    distrp[1]->compute_deviance_mult(responsevec1,weightvec1,linpredvec1,&d1,aux);

    double phinvu = randnumbers::invPhi2(distrp[1]->cdf(*response[response.size()-2],linpredvec1));
    double phinvv = randnumbers::invPhi2(distrp[0]->cdf(*response[0],linpredvec2));

    double l;

    l =  -0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho;

    *deviance = -2*l+d1+d2;
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
    double phinvu = randnumbers::invPhi2(distrp[1]->cdf(*response1p,true));
    double phinvv = randnumbers::invPhi2(distrp[0]->cdf(*response2p,true));

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

//  if(counter>=997)
//    cout << "counter: " << counter << "\n";

  double hilfs = pow(1 + pow((*linpred), 2), 0.5);
  double rho = (*linpred) / hilfs;
  if (*linpred <= -100)
      rho  = -0.99995;
  else if (*linpred >= 100)
      rho  = 0.99995;

  double orho = 1 - pow(rho, 2);
  double phinvu = randnumbers::invPhi2(distrp[1]->cdf(*response1p,true));
  double phinvv = randnumbers::invPhi2(distrp[0]->cdf(*response2p,true));

  double nu = rho * pow(orho, 0.5) + (hilfs + rho * (*linpred)) * (phinvu * phinvv) - (*linpred) * (pow(phinvu, 2) + pow(phinvv, 2));

  *workingweight = 1 - pow(rho, 4);

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {
    like += - 0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho;
    }

  modify_worklin();

  }

void DISTR_gausscopula::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
	// Kendall tau
  *mu = 2 * std::asin((*linpred[predstart_mumult+(linpred.size()-1)]) / (pow(1 + pow((*linpred[predstart_mumult+(linpred.size()-1)]), 2), 0.5))) / PI ;
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

vector<double> DISTR_gausscopula::derivative(double & F1, double & F2, double * linpred)
  {

  vector<double> res;

  double rho = (*linpred)/sqrt(1+(*linpred)*(*linpred));
  double phiinvu = randnumbers::invPhi2(F1);
  double phiinvv = randnumbers::invPhi2(F2);

    //first and second derivative of Phi^-1
  double dphiinvu = sqrt(2*PI)/exp(-0.5*phiinvu*phiinvu);
  double ddphiinvu = 2*PI*phiinvu/pow(exp(-0.5*phiinvu*phiinvu),2);

    // first derivative
  double dlc = rho*dphiinvu*(phiinvv-rho*phiinvu)/(1-rho*rho);
    // second derivative
  double ddlc = rho*ddphiinvu*(phiinvv-rho*phiinvu)/(1-rho*rho) - rho*rho*dphiinvu*dphiinvu/(1-rho*rho);
    // return first and second derivative.
  res.push_back(dlc);
  res.push_back(ddlc);

  return res;
  }


double DISTR_gausscopula::logc(double & F1, double & F2, double * linpred)
  {
  double rho = (*linpred)/sqrt(1+(*linpred)*(*linpred));
  double phiinvu = randnumbers::invPhi2(F1);
  double phiinvv = randnumbers::invPhi2(F2);
  double lc = -0.5*log(1-rho*rho) + rho*phiinvu*phiinvv/(1-rho*rho)-0.5*rho*rho*(phiinvu*phiinvu+phiinvv*phiinvv)/(1-rho*rho);
  return lc;
  }

double DISTR_gausscopula::condfc(double & x, double & linpred_F, double & y, double & F2, double * linpred)
  {
  double res = randnumbers::invPhi2(x);
  double rho = (*linpred)/sqrt(1+(*linpred)*(*linpred));
  double help2 = sqrt(1-rho*rho);
  double help3 = randnumbers::invPhi2(F2);

  double help = randnumbers::Phi2( (-linpred_F - rho*help3) / (help2));

  double xstar = x;
  if(y>0)
    xstar = x*(1-help) + help;
  else
    xstar = x*help;
  res = randnumbers::invPhi2(xstar)*help2 + rho*help3 + linpred_F;
  return res;
  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_gausscopula2 --------------------------
//------------------------------------------------------------------------------
void DISTR_gausscopula2::check_errors(void)
  {
  // Note: check marginal distributions & weights.
  if (errors==false)
    {
    errors=false;
    }
  }


DISTR_gausscopula2::DISTR_gausscopula2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_copula_basis(o,r,w)
  {
  family = "Gauss Copula - rho=tanh(eta)";

  outpredictor = true;
  outexpectation = true;
  predictor_name = "rho";
  linpredminlimit=-100;
  linpredmaxlimit=100;
  check_errors();
  }


DISTR_gausscopula2::DISTR_gausscopula2(const DISTR_gausscopula2 & nd)
   : DISTR_copula_basis(DISTR_copula_basis(nd))
  {
  }


const DISTR_gausscopula2 & DISTR_gausscopula2::operator=(
                            const DISTR_gausscopula2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_copula_basis::operator=(DISTR_copula_basis(nd));
  return *this;
  }


double DISTR_gausscopula2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gausscopula2::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = tanh(*linpred[(linpred.size()-1)]);
  }

void DISTR_gausscopula2::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {
   if ((*weight[0] == 0) || (*weight[weight.size()-2] == 0))
     *deviance=0;
   else
     {
     double rho = tanh((*linpred[(linpred.size()-1)]));
     if (*linpred[(linpred.size()-1)] <= -100)
        rho  = -0.99995;
     else if (*linpred[(linpred.size()-1)] >= 100)
        rho  = 0.99995;

     double orho = 1 - pow(rho, 2);

    int s1 = distrp[1]->distrp.size();
    int s2 = distrp[0]->distrp.size();

    vector<double*> linpredvec1;
    vector<double*> responsevec1;
    vector<double*> weightvec1;
    vector<double*> linpredvec2;
    vector<double*> responsevec2;
    vector<double*> weightvec2;


    int j;
    for (j=0;j<(s2+1);j++)
      {
      linpredvec2.push_back(linpred[j]);
      weightvec2.push_back(weight[j]);
      responsevec2.push_back(response[j]);
      }
    int k;
    for (k=0;k<(s1+1);k++)
      {
      linpredvec1.push_back(linpred[s2+1+k]);
      weightvec1.push_back(weight[s2+1+k]);
      responsevec1.push_back(response[s2+1+k]);
      }

    double d1;
    double d2;
    distrp[0]->compute_deviance_mult(responsevec2,weightvec2,linpredvec2,&d2,aux);
    distrp[1]->compute_deviance_mult(responsevec1,weightvec1,linpredvec1,&d1,aux);

    double phinvu = randnumbers::invPhi2(distrp[1]->cdf(*response[response.size()-2],linpredvec1));
    double phinvv = randnumbers::invPhi2(distrp[0]->cdf(*response[0],linpredvec2));

    double l;

     l =  -0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho;

    *deviance = -2*l+d1+d2;
    }

  }

double DISTR_gausscopula2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {
   if (counter==0)
    {
    set_worklin();
    }
    double rho = tanh(*linpred);
    if (*linpred <= -100)
        rho  = -0.99995;
    else if (*linpred >= 100)
        rho  = 0.99995;

    double orho = 1 - pow(rho, 2);
    double phinvu = randnumbers::invPhi2(distrp[1]->cdf(*response1p,true));
    double phinvv = randnumbers::invPhi2(distrp[0]->cdf(*response2p,true));

    double l;

    l = - 0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho;

  modify_worklin();

  return l;

  }

void DISTR_gausscopula2::compute_iwls_wweightschange_weightsone(
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
  double rho = tanh(*linpred);
  if (*linpred <= -100)
      rho  = -0.99995;
  else if (*linpred >= 100)
      rho  = 0.99995;

  double orho = 1 - pow(rho, 2);
  double phinvu = randnumbers::invPhi2(distrp[1]->cdf(*response1p,true));
  double phinvv = randnumbers::invPhi2(distrp[0]->cdf(*response2p,true));

  double nu = rho + ((1+rho*rho)/(1-rho*rho)) * (phinvu * phinvv) -
              0.25*( exp(*linpred)*exp(*linpred) - exp(-(*linpred))*exp(-(*linpred)) ) * (pow(phinvu, 2) + pow(phinvv, 2));

  double drho = 4/((exp((*linpred))+exp(-(*linpred)))*(exp((*linpred))+exp(-(*linpred))));

  *workingweight = 2-drho;//-drho - ((4*rho*drho)/pow(1-rho*rho,2)) * (phinvu * phinvv) +
              //0.5*(exp((*linpred))*exp((*linpred))+exp(-(*linpred))*exp(-(*linpred))) * (pow(phinvu, 2) + pow(phinvv, 2));

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {
    like += - 0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho;
    }

  modify_worklin();

  }

void DISTR_gausscopula2::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
	//Kendall tau
  *mu = 2 * std::asin(tanh(*linpred[predstart_mumult+(linpred.size()-1)]) ) / PI ;
  }


void DISTR_gausscopula2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (rho): tanh(eta)\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gausscopula2::update_end(void)
  {

  }

vector<double> DISTR_gausscopula2::derivative(double & F1, double & F2, double * linpred)
  {

  vector<double> res;

  double rho = tanh(*linpred);
  double phiinvu = randnumbers::invPhi2(F1);
  double phiinvv = randnumbers::invPhi2(F2);

  //first and second derivative of Phi^-1
  double dphiinvu = sqrt(2*PI)/exp(-0.5*phiinvu*phiinvu);
  double ddphiinvu = 2*PI*phiinvu/pow(exp(-0.5*phiinvu*phiinvu),2);

  // first derivative
  double dlc = rho*dphiinvu*(phiinvv-rho*phiinvu)/(1-rho*rho);
  // second derivative
  double ddlc = rho*ddphiinvu*(phiinvv-rho*phiinvu)/(1-rho*rho) - rho*rho*dphiinvu*dphiinvu/(1-rho*rho);
  // return first and second derivative.
  res.push_back(dlc);
  res.push_back(ddlc);
  return res;
  }


double DISTR_gausscopula2::logc(double & F1, double & F2, double * linpred)
  {
  double rho = tanh(*linpred);
  double phiinvu = randnumbers::invPhi2(F1);
  double phiinvv = randnumbers::invPhi2(F2);
  double lc = -0.5*log(1-rho*rho) + rho*phiinvu*phiinvv/(1-rho*rho)-0.5*rho*rho*(phiinvu*phiinvu+phiinvv*phiinvv)/(1-rho*rho);
  return lc;
  }

double DISTR_gausscopula2::condfc(double & x, double & linpred_F, double & y, double & F2, double * linpred)
  {
  double rho = tanh(*linpred);
  double help2 = sqrt(1-rho*rho);
  double help3 = randnumbers::invPhi2(F2);
  double help = randnumbers::Phi2( (-linpred_F - rho*help3) / (help2));
  double xstar = x;
  if(y>0)
     xstar = x*(1-help) + help;
  else
     xstar = x*help;
  double res = randnumbers::invPhi2(xstar)*help2 + rho*help3 + linpred_F;
  return res;
  }



//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_clayton_copula ------------------------
//------------------------------------------------------------------------------
void DISTR_clayton_copula::check_errors(void)
  {
  if (errors==false)
    {
    errors=false;
    }
  }


DISTR_clayton_copula::DISTR_clayton_copula(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_copula_basis(o,r,w)
  {
  family = "Claytoncopula - rho";

  outpredictor = true;
  outexpectation = true;
  predictor_name = "rho";
  linpredminlimit=-10;
  linpredmaxlimit=15;
  check_errors();
  }


DISTR_clayton_copula::DISTR_clayton_copula(const DISTR_clayton_copula & nd)
   : DISTR_copula_basis(DISTR_copula_basis(nd))
  {
  }


const DISTR_clayton_copula & DISTR_clayton_copula::operator=(
                            const DISTR_clayton_copula & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_copula_basis::operator=(DISTR_copula_basis(nd));
  return *this;
  }


double DISTR_clayton_copula::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_clayton_copula::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp((*linpred[(linpred.size()-1)]));
  }

void DISTR_clayton_copula::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {
   if ((*weight[0] == 0) || (*weight[response.size()-2] == 0))
     *deviance=0;
   else
     {
     double rho = exp((*linpred[(linpred.size()-1)]));

     int s1 = dynamic_cast<DISTR_gamlss *>(distrp[1])->distrp.size();
     int s2 = dynamic_cast<DISTR_gamlss *>(distrp[0])->distrp.size();
     vector<double*> linpredvec1;
     vector<double*> responsevec1;
     vector<double*> weightvec1;
     vector<double*> linpredvec2;
     vector<double*> responsevec2;
     vector<double*> weightvec2;
     int j;
     for (j=0;j<(s2+1);j++)
       {
       linpredvec2.push_back(linpred[j]);
       weightvec2.push_back(weight[j]);
       responsevec2.push_back(response[j]);
       }
     int k;
     for (k=0;k<(s1+1);k++)
       {
       linpredvec1.push_back(linpred[s2+1+k]);
       weightvec1.push_back(weight[s2+1+k]);
       responsevec1.push_back(response[s2+1+k]);
       }

     double d1;
     double d2;
     distrp[0]->compute_deviance_mult(responsevec2,weightvec2,linpredvec2,&d2,aux);
     distrp[1]->compute_deviance_mult(responsevec1,weightvec1,linpredvec1,&d1,aux);

     double u = distrp[1]->cdf(*response[response.size()-2],linpredvec1);
     double v = distrp[0]->cdf(*response[0],linpredvec2);
     if(optionsp->rotation==90)
       {
       u = 1-u;
       //rho=-rho;
       }
     else if(optionsp->rotation==270)
       {
       v = 1-v;
       //rho=-rho;
       }
     else if(optionsp->rotation==180)
       {
       u = 1-u;
       v = 1-v;
       }
     double logu = log(u);
     double logv = log(v);
     double urho = pow(u, -rho);
     double vrho = pow(v, -rho);
     double arg = urho + vrho - 1;
     double l;

     l = log(rho + 1) - (1 + rho) * (logu + logv) - (2 + 1 / rho) * log(arg);
         //+distrp[0]->logpdf(*response[0])+distrp[1]->logpdf(*response[response.size()-2]);

    *deviance = -2*l+d1+d2;
    }

  }

double DISTR_clayton_copula::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {
  if (counter==0)
    {
    set_worklin();
    }
  double rho = exp((*linpred));
  double u = distrp[1]->cdf(*response1p,true);
  double v = distrp[0]->cdf(*response2p,true);
  if(optionsp->rotation==90)
    {
    u = 1-u;
    }
  else if(optionsp->rotation==270)
    {
    v = 1-v;
    }
  else if(optionsp->rotation==180)
    {
    u = 1-u;
    v = 1-v;
    }
  double logu = log(u);
  double logv = log(v);
  double urho = pow(u, -rho);
  double vrho = pow(v, -rho);
  double arg = urho + vrho - 1;
  double l;

  l = log(rho + 1) - (1 + rho) * (logu + logv) - (2 + 1 / rho) * log(arg);

  modify_worklin();

  return l;

  }

void DISTR_clayton_copula::compute_iwls_wweightschange_weightsone(
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
  double rho = exp((*linpred));
  double u = distrp[1]->cdf(*response1p,true);
  double v = distrp[0]->cdf(*response2p,true);
  if(optionsp->rotation==90)
    {
    u = 1-u;
    }
  else if(optionsp->rotation==270)
    {
    v = 1-v;
    }
  else if(optionsp->rotation==180)
    {
    u = 1-u;
    v = 1-v;
    }
  double logu = log(u);
  double logv = log(v);
  double urho = pow(u, -rho);
  double vrho = pow(v, -rho);
  double arg = urho + vrho - 1;

  double nu = rho / (rho + 1) - rho * (logu + logv) + log(arg) / rho + (2 * rho + 1) * (logu * urho + logv * vrho) / arg;
  *workingweight = -rho / pow(rho + 1, 2) + rho * (logu + logv) + log(arg) / rho + (1 - 2 * rho) * (logu * urho + logv * vrho) / arg
                        - ((pow(rho, 2) * (2 + 1 / rho)) / (arg)) * (pow((logu * urho + logv * vrho), 2) / arg - pow(logu, 2) * urho - pow(logv, 2) * vrho );

  if((*workingweight) <= 0)
    *workingweight = 0.0001;

  *workingresponse = *linpred + nu/(*workingweight);


  if (compute_like)
    {
    like += log(rho + 1) - (1 + rho) * (logu + logv) - (2 + 1 / rho) * log(arg);
    }

  modify_worklin();

  }

void DISTR_clayton_copula::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double arg = exp((*linpred[predstart_mumult+(linpred.size()-1)]));
  *mu = arg / (arg + 2);
  }


void DISTR_clayton_copula::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (rho): exp(eta)\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_clayton_copula::update_end(void)
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
    *pmu = exp(*worklin);
    }

  }

vector<double> DISTR_clayton_copula::derivative(double & F1, double & F2, double * linpred)
  {
  vector<double> res;
//////////////////////////////////////////////////////////

  double rho = exp(*linpred);

  double logu = log(F1);
  double logv = log(F2);

  double arg = pow(F1, -rho) + pow(F2, -rho) - 1;
  // first derivative
  double dlc = -(1+rho)/F1+(2+1/rho)*rho*pow(F1,(-rho-1))/arg;
  // second derivative
  double ddlc = (1+rho)/(F1*F1)+(2+1/rho)*pow(rho*pow(F1,(-rho-1))/arg,2)-(2+1/rho)*rho*(rho+1)*pow(F1,(-rho-2))/arg;


  // return first and second derivative.
  res.push_back(dlc);
  res.push_back(ddlc);
  return res;
  }

double DISTR_clayton_copula::logc(double & F1, double & F2, double * linpred)
  {
  double rho = exp((*linpred));

  double arg = pow(F1, -rho) + pow(F2, -rho) - 1;


  double lc = log(rho + 1) - (1 + rho) * (log(F1) + log(F2)) - (2 + 1 / rho) * log(arg);
  return lc;
  }

double DISTR_clayton_copula::condfc(double & x, double & linpred_F, double & y, double & F2, double * linpred)
  {
  // ToDo
  double rho = exp(*linpred);
  double help1 = 0;

  double xstar = x;
  double argPhi=0;

  if(optionsp->rotation==90)
    {
    help1 = 1-pow(F2,-rho-1)*pow((pow(1-randnumbers::Phi2(-linpred_F),-rho)+pow(F2,-rho)-1),-1/rho-1);

    if(y>0)
      xstar = x*(1-help1) + help1;
    else
      xstar = x*help1;

    argPhi = 1-pow(pow(((1-xstar)/pow(F2,-rho-1)),-rho/(rho+1)) + 1 - pow(F2,-rho),-1/rho);
    }
  else if(optionsp->rotation==270)
    {
    help1 = pow(1-F2,-rho-1)*pow((pow(randnumbers::Phi2(-linpred_F),-rho)+pow(1-F2,-rho)-1),-1/rho-1);
    if(y>0)
      xstar = x*(1-help1) + help1;
    else
      xstar = x*help1;

    argPhi = pow(pow((xstar/pow(1-F2,-rho-1)),-rho/(rho+1)) + 1 - pow(1-F2,-rho),-1/rho);
    }
  else if(optionsp->rotation==180)
    {
    help1 = 1-pow(1-F2,-rho-1)*pow((pow(1-randnumbers::Phi2(-linpred_F),-rho)+pow(1-F2,-rho)-1),-1/rho-1);
    if(y>0)
      xstar = x*(1-help1) + help1;
    else
      xstar = x*help1;

    argPhi = 1-pow(pow(((1-xstar)/pow(1-F2,-rho-1)),-rho/(rho+1)) + 1 - pow(1-F2,-rho),-1/rho);
    }
  else
    {
    help1 = pow(F2,-rho-1)*pow((pow(randnumbers::Phi2(-linpred_F),-rho)+pow(F2,-rho)-1),-1/rho-1);
    if(y>0)
      xstar = x*(1-help1) + help1;
    else
      xstar = x*help1;

    argPhi = pow( ( pow( (xstar/pow(F2,-rho-1)) ,-rho/(rho+1) ) + 1 - pow(F2,-rho)) ,-1/rho);
    }


  double res = randnumbers::invPhi2(argPhi)  + linpred_F;
  return res;
  }

} // end: namespace MCMC
