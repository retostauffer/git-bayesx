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
  response1 = nd.response1;
  response1p = nd.response1p;
  linpredp = nd.linpredp;
  }


const DISTR_gausscopula & DISTR_gausscopula::operator=(
                            const DISTR_gausscopula & nd)
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


double DISTR_gausscopula::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gausscopula::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[0]) / pow(1+ pow((*linpred[0]), 2), 0.5);
  }

void DISTR_gausscopula::set_worklin(void)
  {
  DISTR_gamlss::set_worklin();
  response1p = response1.getV();
  response2p = response2.getV();
  }

void DISTR_gausscopula::modify_worklin(void)
  {
  DISTR_gamlss::modify_worklin();
  if (counter<nrobs)
    {
    response1p++;
    response2p++;
    }
  }

void DISTR_gausscopula::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {
/*  cout << linpred.size() << endl;
  cout << response.size() << endl;
  cout << *weight[(linpred.size()-1)] << endl;
  cout << (*linpred[(linpred.size()-1)]) <<endl;*/
   if ((*weight[0] == 0) || (*weight[weight.size()-2] == 0))
     *deviance=0;
   else
     {
 /*    cout << (*response[(response.size()-1)]) <<endl;
     cout << (*response[(response.size()-2)]) <<endl;
     cout << (*response[(response.size()-3)]) <<endl;
     cout << (*response[(response.size()-4)]) <<endl;
     cout << (*response[(response.size()-5)]) <<endl;*/
     double rho = (*linpred[(linpred.size()-1)]) / pow(1 + pow((*linpred[(linpred.size()-1)]), 2), 0.5);
  //   cout << rho <<endl;
     if (*linpred[(linpred.size()-1)] <= -100)
        rho  = -0.99995;
     else if (*linpred[(linpred.size()-1)] >= 100)
        rho  = 0.99995;

     double orho = 1 - pow(rho, 2);
   //  double phinvu = randnumbers::invPhi2(distrp[1]->cdf(*response[response.size()-2],true));
   //  double phinvv = randnumbers::invPhi2(distrp[0]->cdf(*response[0],true));
 //    cout << "phiinvu: " << randnumbers::invPhi2(distrp[1]->cdf(*response[response.size()-2],true)) <<endl;
 //    cout << "phiinvv: " << randnumbers::invPhi2(distrp[0]->cdf(*response[0],true)) <<endl;
    /* cout << "linpred: " << (*linpred[(linpred.size()-1)]) <<endl;
     cout << "rho: " << rho <<endl;
     cout << "orho: " << orho <<endl;
     cout << "pdf(y1): " << distrp[1]->logpdf(*response[response.size()-2]) <<endl;
     cout << "pdf(y2): " << distrp[0]->logpdf(*response[0]) <<endl;*/
     int s1 = dynamic_cast<DISTR_gamlss *>(distrp[1])->distrp.size();
     int s2 = dynamic_cast<DISTR_gamlss *>(distrp[0])->distrp.size();
    // cout << "s1: " << s1 <<endl;
    // cout << "s2: " << s2 <<endl;
    vector<double*> linpredvec1;
    vector<double*> responsevec1;
    vector<double*> weightvec1;
    vector<double*> linpredvec2;
    vector<double*> responsevec2;
    vector<double*> weightvec2;

  //  vector<double> linpredvecval1;
  //  vector<double> linpredvecval2;

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
/*    for (j=0;j<(s2+1);j++)
      {
      linpredvecval2.push_back(*linpredvec2[j]);
      }
    for (k=0;k<(s1+1);k++)
      {
      linpredvecval1.push_back(*linpredvec1[s2+1+k]);
      }*/
    double d1;
    double d2;
    distrp[0]->compute_deviance_mult(responsevec2,weightvec2,linpredvec2,&d2,aux);
    distrp[1]->compute_deviance_mult(responsevec1,weightvec1,linpredvec1,&d1,aux);

    double phinvu = randnumbers::invPhi2(distrp[1]->cdf(*response[response.size()-2],linpredvec1));
    double phinvv = randnumbers::invPhi2(distrp[0]->cdf(*response[0],linpredvec2));

 //   cout << "phiinvu: " << phinvu <<endl;
 //    cout << "phiinvv: " << phinvv <<endl;

 /*   cout << "resp1 0 : " << *responsevec1[0] <<endl;
    cout << "resp1 1 : " << *responsevec1[1] <<endl;
    cout << "resp2 0 : " << *responsevec2[0] <<endl;
    cout << "resp2 1 : " << *responsevec2[1] <<endl;
    cout << "d1 : " << d1 <<endl;
    cout << "d2 : " << d2 <<endl;*/

     double l;

      l =  -0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho +d1+d2;
         //  +distrp[0]->logpdf(*response[0])+distrp[1]->logpdf(*response[response.size()-2]);

    // cout << "l: " << l <<endl;
    *deviance = -2*l;
    }

  }

double DISTR_gausscopula::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {
 /* if(counter<=3)
    {
    cout << "DISTR_gausscopula::loglikelihood_weightsone\n";
    cout << "counter " << counter << "\n";
    }*/
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

 /*  if(counter<=3)
    {
    cout << "DISTR_gausscopula::compute_iwls_wweightschange_weightsone\n";
    cout << "counter " << counter << "\n";
    }*/
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

/*vector<double> DISTR_gausscopula::derivative(double & F, int & copulapos)
  {

  double Fa;
  if(copulapos==0)
    {
    Fa = distrp[1]->cdf(*response1p,false);
    }
  else
    {
    Fa = distrp[0]->cdf(*response2p,false);
    }
    //only works when copula is symmetric in its arguments!!
    return derivative(F, Fa, linpredp);
  }
*/

vector<double> DISTR_gausscopula::derivative(double & F1, double & F2, double * linpred)
  {
 /*  if(counter<=3)
    {
    cout << "DISTR_gausscopula::derivative\n";
    cout << "counter " << counter << "\n";
    }*/

  vector<double> res;
//////////////////////////////////////////////////////////

  double rho = (*linpred)/sqrt(1+(*linpred)*(*linpred));
  double phiinvu = randnumbers::invPhi2(F1);
  double phiinvv = randnumbers::invPhi2(F2);

  //first and second derivative of Phi^-1
  double dphiinvu = sqrt(2*PI)/exp(-0.5*phiinvu*phiinvu);
  double ddphiinvu = 2*PI*phiinvu/pow(exp(-0.5*phiinvu*phiinvu),2);

  // first derivative
  double dlc = rho*dphiinvu*(phiinvv-rho*phiinvu)/(1-rho*rho);
  // second derivative
  double ddlc = rho*ddphiinvu*(phiinvv-rho*phiinvu)/(1-rho*rho) - rho*rho*dphiinvu*dphiinvu/(1-rho*rho);;
  // return first and second derivative.
  res.push_back(dlc);
  res.push_back(ddlc);
  return res;
  }

vector<double> DISTR_gausscopula::logc(double & F, int & copulapos, const bool & deriv)
  {
 /* if(counter<=3)
    {
    cout << "DISTR_gausscopula::logc(double & F, int & copulapos, const bool & deriv)\n";
    cout << "counter " << counter << "\n";
    }*/
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
  vector<double> res;
 // double eta = (*linpred);
  if(copulapos==0)
    {
    //implement Fa
    //cdf virtual in distr hat nur ein Argument!
    Fa = distrp[1]->cdf(*response1p,true);

    res.push_back(logc(Fa, F, linpredp));
    }
  else
    {
    // implement Fa
    Fa = distrp[0]->cdf(*response2p,true);

    res.push_back(logc(F, Fa, linpredp));
    }

  if(deriv)
    {
    vector<double> derivs = derivative(F, Fa, linpredp);
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

double DISTR_gausscopula::logc(double & F1, double & F2, double * linpred)
  {
 /* if(counter<=3)
    {
    cout << "DISTR_gausscopula::logc(double & F1, double & F2, double * linpred)\n";
    cout << "counter " << counter << "\n";
    }*/
  double rho = (*linpred)/sqrt(1+(*linpred)*(*linpred));
  double phiinvu = randnumbers::invPhi2(F1);
  double phiinvv = randnumbers::invPhi2(F2);
  double lc = -0.5*log(1-rho*rho) + rho*phiinvu*phiinvv/(1-rho*rho)-0.5*rho*rho*(phiinvu*phiinvu+phiinvv*phiinvv)/(1-rho*rho);
  return lc;
  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_clayton_copula ------------------------
//------------------------------------------------------------------------------
void DISTR_clayton_copula::check_errors(void)
  {
  // Note: check marginal distributions & weights.
  if (errors==false)
    {
    errors=false;
    }
  }


DISTR_clayton_copula::DISTR_clayton_copula(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
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
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  response1 = nd.response1;
  response1p = nd.response1p;
  linpredp = nd.linpredp;
  }


const DISTR_clayton_copula & DISTR_clayton_copula::operator=(
                            const DISTR_clayton_copula & nd)
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


double DISTR_clayton_copula::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_clayton_copula::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp((*linpred[0]));
  }

void DISTR_clayton_copula::set_worklin(void)
  {
  DISTR_gamlss::set_worklin();
  response1p = response1.getV();
  response2p = response2.getV();
  }

void DISTR_clayton_copula::modify_worklin(void)
  {
  DISTR_gamlss::modify_worklin();
  if (counter<nrobs)
    {
    response1p++;
    response2p++;
    }
  }

void DISTR_clayton_copula::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {
   if (*weight[0] == 0)
     *deviance=0;
   else
     {
     double rho = exp((*linpred[(linpred.size()-1)]));

     double u = distrp[1]->cdf(*response[response.size()-2],true);
     double v = distrp[0]->cdf(*response[0],true);
     double logu = log(u);
     double logv = log(v);
     double urho = pow(u, -rho);
     double vrho = pow(v, -rho);
     double arg = urho + vrho - 1;
     double l;

     l = log(rho + 1) - (1 + rho) * (logu + logv) - (2 + 1 / rho) * log(arg)
         +distrp[0]->logpdf(*response[0])+distrp[1]->logpdf(*response[response.size()-2]);

    *deviance = -2*l;
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
  double arg = exp((*linpred[predstart_mumult]));
  *mu = arg / (arg + 2);
  }


void DISTR_clayton_copula::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (rho): \n");
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

  double rho = (*linpred)/sqrt(1+(*linpred)*(*linpred));
  double phiinvu = randnumbers::invPhi2(F1);
  double phiinvv = randnumbers::invPhi2(F2);

  //first and second derivative of Phi^-1
  double dphiinvu = sqrt(2*PI)/exp(-0.5*phiinvu*phiinvu);
  double ddphiinvu = 2*PI*phiinvu/pow(exp(-0.5*phiinvu*phiinvu),2);

  // first derivative
  double dlc = rho*dphiinvu*(phiinvv-rho*phiinvu)/(1-rho*rho);
  // second derivative
  double ddlc = rho*ddphiinvu*(phiinvv-rho*phiinvu)/(1-rho*rho) - rho*rho*dphiinvu*dphiinvu/(1-rho*rho);;
  // return first and second derivative.
  res.push_back(dlc);
  res.push_back(ddlc);
  return res;
  }

vector<double> DISTR_clayton_copula::logc(double & F, int & copulapos, const bool & deriv)
  {
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
  vector<double> res;
 // double eta = (*linpred);
  if(copulapos==0)
    {
    //implement Fa
    //cdf virtual in distr hat nur ein Argument!
    Fa = distrp[1]->cdf(*response1p,true);

    res.push_back(logc(Fa, F, linpredp));
    }
  else
    {
    // implement Fa
    Fa = distrp[0]->cdf(*response2p,true);

    res.push_back(logc(F, Fa, linpredp));
    }

  if(deriv)
    {
    vector<double> derivs = derivative(F, Fa, linpredp);
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

double DISTR_clayton_copula::logc(double & F1, double & F2, double * linpred)
  {
  double rho = (*linpred)/sqrt(1+(*linpred)*(*linpred));
  double phiinvu = randnumbers::invPhi2(F1);
  double phiinvv = randnumbers::invPhi2(F2);
  double lc = -0.5*log(1-rho*rho) + rho*phiinvu*phiinvv/(1-rho*rho)-0.5*rho*rho*(phiinvu*phiinvu+phiinvv*phiinvv)/(1-rho*rho);
  return lc;
  }


} // end: namespace MCMC
