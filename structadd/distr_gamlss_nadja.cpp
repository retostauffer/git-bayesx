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

#include "distr_gamlss_nadja.h"


namespace MCMC
{

//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_invgaussian_sigma2 ----------------------
//------------------------------------------------------------------------------


DISTR_invgaussian_sigma2::DISTR_invgaussian_sigma2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Inverse Gaussian - sigma2";
  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma2";

  }


DISTR_invgaussian_sigma2::DISTR_invgaussian_sigma2(const DISTR_invgaussian_sigma2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_invgaussian_sigma2 & DISTR_invgaussian_sigma2::operator=(
                            const DISTR_invgaussian_sigma2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_invgaussian_sigma2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_invgaussian_sigma2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sigma_2 = exp((*linpred));

  double l;

     l = -0.5*log(sigma_2)-pow((((*response))-(*worklin[0])),2)/(2*(*response)*pow((*worktransformlin[0]),2)*sigma_2);


  modify_worklin();

  return l;

  }

void DISTR_invgaussian_sigma2::compute_iwls_wweightschange_weightsone(
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

    double sigma_2 = exp(*linpred);


    double nu = -0.5 + (pow(((*response)-(*worklin[0])),2))/(2*(*response)*(pow((*worktransformlin[0]),2))*sigma_2);



    *workingweight = 0.5;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(sigma_2)-pow((((*response))-(*worklin[0])),2)/(2*(*response)*pow((*worktransformlin[0]),2)*sigma_2);

      }

  modify_worklin();

  }


void DISTR_invgaussian_sigma2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_invgaussian_sigma2::update_end(void)
  {

  // helpmat1 stores sigma2

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


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_invgaussian_mu ------------------------
//------------------------------------------------------------------------------
void DISTR_invgaussian_mu::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight > 0)
        {

        if (*workresp <= 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative/zero response values encountered\n");
          }


        }
      else if (*workweight == 0)
        {
        }
      else
        {
        errors=true;
        errormessages.push_back("ERROR: negative weights encountered\n");
        }

      i++;
      workresp++;
      workweight++;

      }

    }

  }


DISTR_invgaussian_mu::DISTR_invgaussian_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Inverse Gaussian - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  }


DISTR_invgaussian_mu::DISTR_invgaussian_mu(const DISTR_invgaussian_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_invgaussian_mu & DISTR_invgaussian_mu::operator=(
                            const DISTR_invgaussian_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_invgaussian_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_sigma2
   // *linpred[1] = eta_mu

   if (*weight[1] == 0)
     *deviance=0;
   else
     {
     double sigma_2 = exp(*linpred[0]);
     double mu = exp(*linpred[1]);

     double l;

       l = -0.5*log(2*PI)-0.5*log(sigma_2)-1.5*log((*response[1]))-(1/(2*(*response[1])))*pow((((*response[0]))-mu),2)/(pow(mu,2)*sigma_2);


    *deviance = -2*l;
    }

  }


double DISTR_invgaussian_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_invgaussian_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = exp(*linpred);

  double l;

     l = -pow((((*response))-mu),2)/(2*(*response)*pow(mu,2)*(*worktransformlin[0]));

  modify_worklin();

  return l;

  }


void DISTR_invgaussian_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double mu = exp(*linpred);

    double nu = ((*response)-mu)/(pow(mu,2)*(*worktransformlin[0]));

    *workingweight = 1/(mu*(*worktransformlin[0]));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow((((*response))-mu),2)/(2*(*response)*pow(mu,2)*(*worktransformlin[0]));

      }


  modify_worklin();

  }


void DISTR_invgaussian_mu::compute_mu_mult(vector<double *> linpred,double * mu)
  {

  *mu = exp((*linpred[predstart_mumult+1]));

  }


void DISTR_invgaussian_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_invgaussian_mu::update_end(void)
  {


  // helpmat1 stores (eta_mu)

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
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_betainf_mu --------------------------
//------------------------------------------------------------------------------


DISTR_betainf_mu::DISTR_betainf_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Beta Inflated - mu";
  }


DISTR_betainf_mu::DISTR_betainf_mu(const DISTR_betainf_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_betainf_mu & DISTR_betainf_mu::operator=(
                            const DISTR_betainf_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_betainf_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_betainf_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = (1-sigma2)/sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double exp_lin = exp((*linpred));
  double mu = exp_lin/(1+exp_lin);
  double mu_worktrans = mu*(*worktransformlin[0]);
  double one_minus_mu_worktrans = (1-mu)*(*worktransformlin[0]);

  double l;

     l = mu*(*worktransformlin[0])*log(*response) +
		 (1-mu)*(*worktransformlin[0])*log(1-(*response))-
		 randnumbers::lngamma_exact(mu_worktrans)-
		 randnumbers::lngamma_exact(one_minus_mu_worktrans);

  modify_worklin();

  return l;

  }


void DISTR_betainf_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = (1-sigma2)/sigma2;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double exp_lin = exp(*linpred);

    double mu = exp_lin/(1+exp_lin);

    double mu_worktrans = mu*(*worktransformlin[0]);

    double one_minus_mu_worktrans = (1-mu)*(*worktransformlin[0]);

    double nu = mu*one_minus_mu_worktrans*log((*response))-mu*one_minus_mu_worktrans*log((1-(*response)))+
				mu*one_minus_mu_worktrans*(randnumbers::digamma_exact(one_minus_mu_worktrans)-randnumbers::digamma_exact(mu_worktrans));

    *workingweight = pow((*worktransformlin[0]),2)*pow(mu,2)*pow((1-mu),2)*(randnumbers::trigamma_exact(one_minus_mu_worktrans)+randnumbers::trigamma_exact(mu_worktrans));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += mu_worktrans*log((*response))+one_minus_mu_worktrans*log(1-(*response))-
		        randnumbers::lngamma_exact(mu_worktrans)-randnumbers::lngamma_exact((one_minus_mu_worktrans));

      }


  modify_worklin();

  }



void DISTR_betainf_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): logit\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_betainf_mu::update_end(void)
  {


  // helpmat1 stores exp(eta_mu)/(1+exp(eta_mu))

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double exp_lin;
  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin)/(1+exp(*worklin));
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_betainf_sigma2 ------------------------
//------------------------------------------------------------------------------


DISTR_betainf_sigma2::DISTR_betainf_sigma2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Beta Inflated - sigma2";
  }


DISTR_betainf_sigma2::DISTR_betainf_sigma2(const DISTR_betainf_sigma2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_betainf_sigma2 & DISTR_betainf_sigma2::operator=(
                            const DISTR_betainf_sigma2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_betainf_sigma2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_betainf_sigma2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu)/(1+exp(eta_mu));

  if (counter==0)
    {
    set_worklin();
    }

  double exp_lin = exp((*linpred));
  double sigma_2 = exp_lin/(1+exp_lin);
  double help = (1-sigma_2)/sigma_2;
  double mu_help = (*worktransformlin[0])*help;
  double one_minus_mu_help = (1-(*worktransformlin[0]))*help;

  double l;

     l = mu_help*log(*response) +
         one_minus_mu_help*log(1-(*response)) -
         randnumbers::lngamma_exact(mu_help) - randnumbers::lngamma_exact(one_minus_mu_help) +
         randnumbers::lngamma_exact(help);


  modify_worklin();

  return l;

  }

void DISTR_betainf_sigma2::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu)/(1+exp(eta_mu));

  if (counter==0)
    {
    set_worklin();
    }

    double exp_lin = exp(*linpred);

    double sigma_2 = exp_lin/(1+exp_lin);

    double help = (1-sigma_2)/sigma_2;
    double mu_help = (*worktransformlin[0])*help;
    double one_minus_mu_help = (1-(*worktransformlin[0]))*help;

    double nu = -help*( -(*worktransformlin[0])*randnumbers::digamma_exact(mu_help)
                       - (1-(*worktransformlin[0]))*randnumbers::digamma_exact(one_minus_mu_help) +
                       randnumbers::digamma_exact(help) + (*worktransformlin[0])*log(*response)
                       + (1-(*worktransformlin[0]))*log(1-(*response)) );



    *workingweight =  pow(help,2)*( pow((1-(*worktransformlin[0])),2)*randnumbers::trigamma_exact(one_minus_mu_help)+
                                   pow(((*worktransformlin[0])),2)*randnumbers::trigamma_exact(mu_help) -
                                   randnumbers::trigamma_exact(help) );

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += mu_help*log(*response) +
         one_minus_mu_help*log(1-(*response)) -
         randnumbers::lngamma_exact(mu_help) - randnumbers::lngamma_exact(one_minus_mu_help) +
         randnumbers::lngamma_exact(help);

      }

  modify_worklin();

  }


void DISTR_betainf_sigma2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): logit\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_betainf_sigma2::update_end(void)
  {

  // helpmat1 stores (1-sigma2)/sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double exp_lin;
  double sigma_2;
  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    exp_lin = exp(*worklin);
    sigma_2 = exp_lin/(1+exp_lin);
    *pmu = (1-sigma_2)/sigma_2;
    }

  }

//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_betainf_nu --------------------------
//------------------------------------------------------------------------------


DISTR_betainf_nu::DISTR_betainf_nu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Beta Inflated - nu";
  }


DISTR_betainf_nu::DISTR_betainf_nu(const DISTR_betainf_nu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_betainf_nu & DISTR_betainf_nu::operator=(
                            const DISTR_betainf_nu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_betainf_nu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_betainf_nu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of tau equation
  // *worktransformlin[0] = exp(eta_tau);

  if (counter==0)
    {
    set_worklin();
    }

  double nup = exp((*linpred));

  double l;

    if ((*response)==0)
     {
         l = log(nup) - log(1+(*worktransformlin[0])+nup);
     }
     else
     {
        l =  - log(1+(*worktransformlin[0])+nup);
     }

  modify_worklin();

  return l;

  }


void DISTR_betainf_nu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of tau equation
  // *worktransformlin[0] = exp(eta_tau);

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double nup = exp(*linpred);

    double hilfs = 1+ (*worktransformlin[0]) + nup;

    double nu = -nup/hilfs;

    if ((*response)==0)
    {
        nu += 1;
    }

    *workingweight = nup*(1+(*worktransformlin[0]))/(pow(hilfs,2));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        if ((*response)==0)
        {
            like += log(nup) - log(1+(*worktransformlin[0])+nup);
        }
        else
        {
            like +=  - log(1+(*worktransformlin[0])+nup);
        }
      }


  modify_worklin();

  }



void DISTR_betainf_nu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (nu): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_betainf_nu::update_end(void)
  {


  // helpmat1 stores exp(eta_tau)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double exp_lin;
  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_betainf_tau --------------------------
//------------------------------------------------------------------------------


DISTR_betainf_tau::DISTR_betainf_tau(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Beta Inflated - tau";
  }


DISTR_betainf_tau::DISTR_betainf_tau(const DISTR_betainf_tau & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_betainf_tau & DISTR_betainf_tau::operator=(
                            const DISTR_betainf_tau & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_betainf_tau::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_tau
   // *linpred[1] = eta_nu
   // *linpred[2] = eta_sigma2
   // *linpred[3] = eta_mu

   if (*weight[3] == 0)
     *deviance=0;
   else
     {
     double tau =  exp(*linpred[3]);
     double nup = exp(*linpred[2]);
     double sigma_2 = exp(*linpred[1])/(1+exp(*linpred[1]));
     double mu = exp(*linpred[0])/(1+exp(*linpred[0]));
     double help = (1-sigma_2)/sigma_2;
     double one_minus_mu_help = (1-mu)*help;
     double mu_help = mu*help;
     double one_nup_tau = (1+nup+tau);

     double l;

     if ((*response[3])==0)
     {
         l = log(nup) - log(one_nup_tau);
     }
     else if ((*response[3])==1)
     {
        l = log(tau) - log(one_nup_tau);
     }
     else
       l = (mu_help-1)*log(*response[3]) +
			(one_minus_mu_help-1)*log(1-(*response[3]))-
			randnumbers::lngamma_exact(mu_help)-
			randnumbers::lngamma_exact(one_minus_mu_help)+
			randnumbers::lngamma_exact(help)- log(one_nup_tau);


    *deviance = -2*l;
    }

  }


double DISTR_betainf_tau::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_betainf_tau::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of nu equation
  // *worktransformlin[0] = exp(eta_nu);

  if (counter==0)
    {
    set_worklin();
    }

  double tau = exp((*linpred));

  double l;

    if ((*response)==1)
     {
         l = log(tau) - log(1+(*worktransformlin[0])+tau);
     }
     else
     {
        l =  - log(1+(*worktransformlin[0])+tau);
     }

  modify_worklin();

  return l;

  }


void DISTR_betainf_tau::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of nu equation
  // *worktransformlin[0] = exp(eta_nu);

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double tau = exp(*linpred);

    double hilfs = 1+ (*worktransformlin[0]) + tau;

    double nu = -tau/hilfs;

    if ((*response)==1)
    {
        nu += 1;
    }

    *workingweight = tau*(1+(*worktransformlin[0]))/(pow(hilfs,2));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        if ((*response)==1)
        {
            like += log(tau) - log(1+(*worktransformlin[0])+tau);
        }
        else
        {
            like +=  - log(1+(*worktransformlin[0])+tau);
        }
      }


  modify_worklin();

  }


void DISTR_betainf_tau::compute_mu_mult(vector<double *> linpred,double * mu)
  {
  double exp_lin_tau = exp(*linpred[3]);
  double exp_lin_nup = exp(*linpred[2]);
  double exp_lin_mu = exp(*linpred[0]);
  double hilfs = (1+exp_lin_tau+exp_lin_nup);
  *mu = (1-(exp_lin_tau+exp_lin_nup)/hilfs)*(exp_lin_mu/(1+exp_lin_mu))+exp_lin_tau/(hilfs);
  }


void DISTR_betainf_tau::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (tau): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_betainf_tau::update_end(void)
  {


  // helpmat1 stores exp(eta_nu)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double exp_lin;
  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_pareto_p ---------------------------
//------------------------------------------------------------------------------


DISTR_pareto_p::DISTR_pareto_p(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Pareto - p";
  }


DISTR_pareto_p::DISTR_pareto_p(const DISTR_pareto_p & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_pareto_p & DISTR_pareto_p::operator=(
                            const DISTR_pareto_p & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_pareto_p::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_pareto_p::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of b equation
  // *worktransformlin[0] = exp(eta_b);

  if (counter==0)
    {
    set_worklin();
    }

  double p = exp((*linpred));

  double l;

     l = log(p) + p*log(*worktransformlin[0]) - (p)*log((*response)+(*worktransformlin[0]));


  modify_worklin();

  return l;

  }

void DISTR_pareto_p::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of b equation
  // *worktransformlin[0] = exp(eta_b);

  if (counter==0)
    {
    set_worklin();
    }

    double p = exp((*linpred));

    double nu = 1  + p*log(*worktransformlin[0]) - (p)*log((*response)+(*worktransformlin[0]));



    *workingweight = 1;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  log(p) + p*log(*worktransformlin[0]) - (p)*log((*response)+(*worktransformlin[0]));

      }

  modify_worklin();

  }


void DISTR_pareto_p::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (p): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_pareto_p::update_end(void)
  {

  // helpmat1 stores sigma

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


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_pareto_b ------------------------
//------------------------------------------------------------------------------

void DISTR_pareto_b::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight > 0)
        {

        if (*workresp <= 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative/zero response values encountered\n");
          }


        }
      else if (*workweight == 0)
        {
        }
      else
        {
        errors=true;
        errormessages.push_back("ERROR: negative weights encountered\n");
        }

      i++;
      workresp++;
      workweight++;

      }

    }

  }

DISTR_pareto_b::DISTR_pareto_b(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Pareto - b";
  }


DISTR_pareto_b::DISTR_pareto_b(const DISTR_pareto_b & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_pareto_b & DISTR_pareto_b::operator=(
                            const DISTR_pareto_b & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_pareto_b::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_p
   // *linpred[1] = eta_b

   if (*weight[1] == 0)
     *deviance=0;
   else
     {
     double p = exp(*linpred[0]);
     double b = exp(*linpred[1]);

     double l;

       l =   log(p) + p*log(b) - (p+1)*log((*response[1])+b);


    *deviance = -2*l;
    }

  }


double DISTR_pareto_b::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_pareto_b::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of p equation
  // *worktransformlin[0] = p;

  if (counter==0)
    {
    set_worklin();
    }

  double b = exp(*linpred);

  double l;

     l =  (*worktransformlin[0])*log(b) - ((*worktransformlin[0])+1)*log((*response)+b) ;

  modify_worklin();


  return l;

  }


void DISTR_pareto_b::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of p equation
  // *worktransformlin[0] = p;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double b = exp(*linpred);

    double nu = (*worktransformlin[0]) - (((*worktransformlin[0])+1)*b)/((*response)+b);

    *workingweight = ((*worktransformlin[0]))/((*worktransformlin[0])+2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += (*worktransformlin[0])*log(b) - ((*worktransformlin[0])+1)*log((*response)+b) ;

      }


  modify_worklin();

  }


void DISTR_pareto_b::compute_mu_mult(vector<double *> linpred,double * mu)
  {

   double p = exp((*linpred[0]));
   double b = exp((*linpred[1]));
  *mu = b/(p-1);

  }


void DISTR_pareto_b::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (b): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_pareto_b::update_end(void)
  {


  // helpmat1 stores (eta_mu)

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
//    double t = 0;
    }

  }



//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_dagum_p --------------------------
//------------------------------------------------------------------------------

DISTR_dagum_p::DISTR_dagum_p(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Dagum - p";
  }


DISTR_dagum_p::DISTR_dagum_p(const DISTR_dagum_p & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_dagum_p & DISTR_dagum_p::operator=(
                            const DISTR_dagum_p & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_dagum_p::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_dagum_p::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of b equation
  // *worktransformlin[0] = exp(eta_b);
  // *worklin[1] = linear predictor of a equation
  // *worktransformlin[1] = exp(eta_a);

  if (counter==0)
    {
    set_worklin();
    }

  double p = exp((*linpred));
  double hilfs = pow((*response)/(*worktransformlin[0]),(*worktransformlin[1]));

  double l;

     l = log(p) + (*worktransformlin[1])*p*log((*response)) - (*worktransformlin[1])*p*log((*worktransformlin[0]))
        -p*log(1+hilfs);


  modify_worklin();

  return l;

  }

void DISTR_dagum_p::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of b equation
  // *worktransformlin[0] = exp(eta_b);
  // *worklin[1] = linear predictor of a equation
  // *worktransformlin[1] = exp(eta_a);


  if (counter==0)
    {
    set_worklin();
    }

    double p = exp((*linpred));
    double hilfs = pow((*response)/(*worktransformlin[0]),(*worktransformlin[1]));

    double nu = 1 + (*worktransformlin[1])*p*log((*response)) - (*worktransformlin[1])*p*log((*worktransformlin[0]))
                -p*log(1+hilfs);

	double exp_linsigma_plus1 = ((*worktransformlin[0])+1);

    *workingweight = 1;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += log(p) + (*worktransformlin[1])*p*log((*response)) - (*worktransformlin[1])*p*log((*worktransformlin[0]))
        -p*log(1+hilfs);

      }

  modify_worklin();

  }


void DISTR_dagum_p::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (p): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_dagum_p::update_end(void)
  {

  // helpmat1 stores tau

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


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_dagum_b ---------------------------
//------------------------------------------------------------------------------


DISTR_dagum_b::DISTR_dagum_b(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Dagum - b";
  }


DISTR_dagum_b::DISTR_dagum_b(const DISTR_dagum_b & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_dagum_b & DISTR_dagum_b::operator=(
                            const DISTR_dagum_b & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_dagum_b::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_dagum_b::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of p equation
  // *worktransformlin[0] = exp(eta_p);
  // *worklin[1] = linear predictor of a equation
  // *worktransformlin[1] = exp(eta_a);

  if (counter==0)
    {
    set_worklin();
    }

  double b = exp((*linpred));
  double hilfs = pow((*response)/b,(*worktransformlin[1]));
  double l;

     l = - (*worktransformlin[1])*(*worktransformlin[0])*log(b) - ((*worktransformlin[0])+1)*log(1+hilfs);


  modify_worklin();

  return l;

  }

void DISTR_dagum_b::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of p equation
  // *worktransformlin[0] = exp(eta_p);
  // *worklin[1] = linear predictor of a equation
  // *worktransformlin[1] = exp(eta_a);

  if (counter==0)
    {
    set_worklin();
    }

    double b = exp((*linpred));
    double hilfs = pow((*response)/b,(*worktransformlin[1]));

    double nu = (*worktransformlin[1]) - (((*worktransformlin[0])+1)*(*worktransformlin[1]))/(1+hilfs) ;



    *workingweight = (((*worktransformlin[0])+1)*pow((*worktransformlin[1]),2)*hilfs)/pow((1+hilfs),2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  - (*worktransformlin[1])*(*worktransformlin[0])*log(b) - ((*worktransformlin[0])+1)*log(1+hilfs);

      }

  modify_worklin();

  }


void DISTR_dagum_b::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (b): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_dagum_b::update_end(void)
  {

  // helpmat1 stores b

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


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_dagum_a -----------------------------
//------------------------------------------------------------------------------


void DISTR_dagum_a::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight > 0)
        {

        if (*workresp <= 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative/zero response values encountered\n");
          }


        }
      else if (*workweight == 0)
        {
        }
      else
        {
        errors=true;
        errormessages.push_back("ERROR: negative weights encountered\n");
        }

      i++;
      workresp++;
      workweight++;

      }

    }

  }


DISTR_dagum_a::DISTR_dagum_a(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Dagum - a";
  }


DISTR_dagum_a::DISTR_dagum_a(const DISTR_dagum_a & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_dagum_a & DISTR_dagum_a::operator=(
                            const DISTR_dagum_a & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_dagum_a::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_p
   // *linpred[1] = eta_b
   // *linpred[2] = eta_a

   if (*weight[2] == 0)
     *deviance=0;
   else
     {
	 double p = exp(*linpred[0]);
     double b = exp(*linpred[1]);
     double a = exp(*linpred[2]);
     double hilfs = pow((*response[2])/b,a);

     double l;

       l = log(a) + log(p) +(a*p-1)*log((*response[2])) - a*p*log(b) - (p+1)*log(1+hilfs);


    *deviance = -2*l;
    }

  }


double DISTR_dagum_a::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_dagum_a::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of p equation
  // *worktransformlin[0] = p;
  // *worklin[1] = linear predictor of b equation
  // *worktransformlin[1] = b;

  if (counter==0)
    {
    set_worklin();
    }

  double a = exp(*linpred);
  double hilfs = pow((*response)/(*worktransformlin[1]),a);
  double l;

     l =  log(a) +(a*(*worktransformlin[0]))*log((*response)) - a*(*worktransformlin[0])*log((*worktransformlin[1]))
            - ((*worktransformlin[0])+1)*log(1+hilfs);
  modify_worklin();

  return l;

  }


void DISTR_dagum_a::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of p equation
  // *worktransformlin[0] = p;
  // *worklin[1] = linear predictor of b equation
  // *worktransformlin[1] = b;


  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double a = exp((*linpred));

    double hilfs = pow((*response)/(*worktransformlin[1]),a);

    double nu = 1 + a*(*worktransformlin[0])*log((*response)/(*worktransformlin[1]))
                - (((*worktransformlin[0])+1)*a*hilfs*log((*response)/(*worktransformlin[1])))/(1+hilfs);

    *workingweight = 1 + (((*worktransformlin[0])+1)*pow(a,2)*hilfs*pow(log((*response)/(*worktransformlin[1])),2))/pow((1+hilfs),2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  log(a) +(a*(*worktransformlin[0]))*log((*response)) - a*(*worktransformlin[0])*log((*worktransformlin[1]))
            - ((*worktransformlin[0])+1)*log(1+hilfs);

      }


  modify_worklin();

  }


void DISTR_dagum_a::compute_mu_mult(vector<double *> linpred,double * mu)
  {

  double exp_lin_a = exp((*linpred[2]));
  double exp_lin_b = exp((*linpred[1]));
  double exp_lin_p = exp((*linpred[0]));
  double help1 = -1/exp_lin_a;
  double help2 = -help1 + exp_lin_p;

  *mu = 0;

  if (exp_lin_a>1)
  {
      *mu = -(exp_lin_b/exp_lin_a)*(randnumbers::gamma_exact(help1)*randnumbers::gamma_exact(help2))/(randnumbers::gamma_exact(exp_lin_p));
  }

  }


void DISTR_dagum_a::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (a): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_dagum_a::update_end(void)
  {


  // helpmat1 stores exp(eta_a)

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
//    double t = 0;
    }

  }




//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_weibull_alpha ---------------------------
//------------------------------------------------------------------------------


DISTR_weibull_alpha::DISTR_weibull_alpha(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Weibull - sigma";
  }


DISTR_weibull_alpha::DISTR_weibull_alpha(const DISTR_weibull_alpha & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_weibull_alpha & DISTR_weibull_alpha::operator=(
                            const DISTR_weibull_alpha & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_weibull_alpha::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_weibull_alpha::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sig = exp((*linpred));

  double l;

     l = (sig)*log(*response) - pow((*response)/(*worktransformlin[0]),sig) -sig*log((*worktransformlin[0])) + log(sig);


  modify_worklin();

  return l;

  }

void DISTR_weibull_alpha::compute_iwls_wweightschange_weightsone(
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

    double sig = exp((*linpred));
    double hilfs1 = pow((*response)/(*worktransformlin[0]),sig);

    double nu = 1 + sig*log((*response)/(*worktransformlin[0]))*(1-hilfs1);



    *workingweight = sig*log((*response)/(*worktransformlin[0]))*( hilfs1 + sig*hilfs1*log((*response)/(*worktransformlin[0])) -1 ) ;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  (sig)*log(*response) - hilfs1 -sig*log((*worktransformlin[0])) + log(sig);

      }

  modify_worklin();

  }


void DISTR_weibull_alpha::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_weibull_alpha::update_end(void)
  {

  // helpmat1 stores sigma

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


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_weibull_lambda ------------------------
//------------------------------------------------------------------------------

void DISTR_weibull_lambda::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight > 0)
        {

        if (*workresp < 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative response values encountered\n");
          }


        }
      else if (*workweight == 0)
        {
        }
      else
        {
        errors=true;
        errormessages.push_back("ERROR: negative weights encountered\n");
        }

      i++;
      workresp++;
      workweight++;

      }

    }

  }


DISTR_weibull_lambda::DISTR_weibull_lambda(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Weibull - mu";
  }


DISTR_weibull_lambda::DISTR_weibull_lambda(const DISTR_weibull_lambda & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_weibull_lambda & DISTR_weibull_lambda::operator=(
                            const DISTR_weibull_lambda & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_weibull_lambda::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_sigma
   // *linpred[1] = eta_mu

   if (*weight[1] == 0)
     *deviance=0;
   else
     {
     double sig = exp(*linpred[0]);
     double mu = exp(*linpred[1]);

     double l;

       l =   (sig-1)*log(*response[1]) - pow((*response[1])/mu,sig) -sig*log(mu) + log(sig) ;


    *deviance = -2*l;
    }

  }


double DISTR_weibull_lambda::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_weibull_lambda::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = sigma;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = exp(*linpred);

  double l;

     l = - pow((*response)/mu,(*worktransformlin[0])) - (*worktransformlin[0])*log(mu) ;

  modify_worklin();


  return l;

  }


void DISTR_weibull_lambda::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = sigma;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double mu = exp(*linpred);
    double hilfs1 = pow((*response)/mu,(*worktransformlin[0]));

    double nu = (*worktransformlin[0])*( hilfs1-1 );

    *workingweight = pow((*worktransformlin[0]),2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += - hilfs1 - (*worktransformlin[0])*log(mu) ;

      }


  modify_worklin();

  }


void DISTR_weibull_lambda::compute_mu_mult(vector<double *> linpred,double * mu)
  {

   double hilfs = 1+1/exp((*linpred[0]));
  *mu = exp((*linpred[1]))*randnumbers::gamma_exact(hilfs);

  }


void DISTR_weibull_lambda::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_weibull_lambda::update_end(void)
  {


  // helpmat1 stores (eta_mu)

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
//    double t = 0;
    }

  }



//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_zinb2_delta --------------------------
//------------------------------------------------------------------------------

/*
DISTR_zinb2_delta::DISTR_zinb2_delta(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Zero-inflated negative binomial - delta";
  }


DISTR_zinb2_delta::DISTR_zinb2_delta(const DISTR_zinb2_delta & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_zinb2_delta & DISTR_zinb2_delta::operator=(
                            const DISTR_zinb2_delta & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_zinb2_delta::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_zinb2_delta::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of pi equation
  // *worktransformlin[0] = pi;
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }


  double delta;


     double delta = exp(*linpred);

     double pot = pow(delta/((*worktransformlin[1])+delta),delta);

     double l;

      if (*response==0)
      {
        l = log((*worktransformlin[0])+(1-(*worktransformlin[0]))*pot);
      }
      else
      {
        double help1 = (*response) + delta;
        double help2 = (*response) + 1;
        l = randnumbers::lngamma_exact(help1)
          - randnumbers::lngamma_exact(help2)
          - randnumbers::lngamma_exact(delta)
          + delta*(*linpred)
          - (delta+(*response))*log(delta+(*worktransformlin[1]));
      }

  modify_worklin();

  return l;

  }

void DISTR_zinb2_delta::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of pi equation
  // *worktransformlin[0] = pi;
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);


  if (counter==0)
    {
    set_worklin();
    }


  double delta;


    double delta = exp(*linpred);
    double pot = pow(delta/((*worktransformlin[1])+delta),delta);
    double hilfs1 = (*response) + delta;

    double nu = delta*(randnumbers::digamma_exact(hilfs1)-randnumbers::digamma_exact(delta)+log(delta/(delta+(*worktransformlin[1])))+((*worktransformlin[1])-(*response))/(delta+(*worktransformlin[1])));

	if (*response==0)
    {
        nu -= (delta*(*worktransformlin[0])*(log(delta/(delta+(*worktransformlin[1])))+(*worktransformlin[1])/(delta+(*worktransformlin[1]))))/((*worktransformlin[0])+(1-(*worktransformlin[0]))*pot);
    }

    *workingweight = pow(nu,2);

    *workingresponse = *linpred + 1/nu;

    if (compute_like)
      {

      if (*response==0)
      {
        like += log((*worktransformlin[0])+(1-(*worktransformlin[0]))*pot);
      }
      else
      {
        double help1 = (*response) + delta;
        double help2 = (*response) + 1;
        like += randnumbers::lngamma_exact(help1)
          - randnumbers::lngamma_exact(help2)
          - randnumbers::lngamma_exact(delta)
          + delta*(*linpred)
          - (delta+(*response))*log(delta+(*worktransformlin[1]));
      }

      }

  modify_worklin();

  }


void DISTR_zinb2_delta::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (delta): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_zinb2_delta::update_end(void)
  {

  // helpmat1 stores tau

  double * worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_zinb2_pi ---------------------------
//------------------------------------------------------------------------------


DISTR_zinb2_pi::DISTR_zinb2_pi(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Zero-inflated negative binomial - pi";
  }


DISTR_zinb2_pi::DISTR_zinb2_pi(const DISTR_zinb2_pi & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_zinb2_pi & DISTR_zinb2_pi::operator=(
                            const DISTR_zinb2_pi & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_zinb2_pi::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_zinb2_pi::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of delta equation
  // *worktransformlin[0] = exp(eta_delta);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }



  double pi = exp((*linpred))/(1+exp((*linpred)));


  double pot = pow((*worktransformlin[0])/((*worktransformlin[1])+(*worktransformlin[0])),(*worktransformlin[0]));

  double l;

      if (*response==0)
      {
        l = log(pi+(1-pi)*pot);
      }
      else
      {
        l = log(1-pi);
      }

  modify_worklin();

  return l;

  }

void DISTR_zinb2_pi::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of delta equation
  // *worktransformlin[0] = exp(eta_delta);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }


  double pi = exp(*linpred)/(1+exp(*linpred));
  double ci = 0.998*explinpredpi/pow((1+explinpredpi),2);
   double pot = pow((*worktransformlin[0])/((*worktransformlin[1])+(*worktransformlin[0])),(*worktransformlin[0]));

    double nu = -ci/(1-pi);
    if (*response==0)
    {
        nu += (ci)/((pi+(1-pi)*pot)*(1-pi));
    }


    *workingweight = pow(nu,2);

    *workingresponse = *linpred + 1/nu;

 //   *workingweight = ( pow(pi,2)*(1-pi)*(1-pot) )/( pi+(1-pi)*pot );

  //  *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

      if (*response==0)
      {
        like += log(pi+(1-pi)*pot);
      }
      else
      {
        like += log(1-pi);
      }

      }

  modify_worklin();

  }


void DISTR_zinb2_pi::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (pi): logit\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_zinb2_pi::update_end(void)
  {

  // helpmat1 stores sigma

  double * worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin)/(1+exp(*worklin));
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_zinb2_mu ------------------------
//------------------------------------------------------------------------------


DISTR_zinb2_mu::DISTR_zinb2_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Zero-inflated negative binomial - mu";
  }


DISTR_zinb2_mu::DISTR_zinb2_mu(const DISTR_zinb2_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_zinb2_mu & DISTR_zinb2_mu::operator=(
                            const DISTR_zinb2_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_zinb2_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_delta
   // *linpred[1] = eta_pi
   // *linpred[2] = eta_mu

  double mu = exp(*linpred[2]);

  double pi = exp(*linpred[1])/(1+exp(*linpred[1]));

  double  delta= exp(*linpred[0]);

   if (*weight[1] == 0)
     *deviance=0;
   else
     {

     double pot = pow(delta/(mu+delta),delta);

     double l;

      if ((*response[1])==0)
      {
        l = log(pi+(1-pi)*pot);
      }
      else
      {
        double help1 = (*response[2]) + delta;
        double help2 = (*response[2]) + 1;
        l = log(1-pi) + randnumbers::lngamma_exact(help1)
          - randnumbers::lngamma_exact(help2)
          - randnumbers::lngamma_exact(delta)
          + delta*(*linpred[0])
          + (*response[2])*(*linpred[2])
          - (delta+(*response[2]))*log(delta+mu);
      }



    *deviance = -2*l;
    }

  }


double DISTR_zinb2_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_zinb2_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of delta equation
  // *worktransformlin[0] = delta;
  // *worklin[1] = linear predictor of pi equation
  // *worktransformlin[1] = pi;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = exp(*linpred);
  double pot = pow((*worktransformlin[0])/(mu+(*worktransformlin[0])),(*worktransformlin[0]));

  double l;

      if (*response==0)
      {
        l = log((*worktransformlin[1])+(1-(*worktransformlin[1]))*pot);
      }
      else
      {
        l = (*response)*(*linpred)
          - ((*worktransformlin[0])+(*response))*log((*worktransformlin[0])+mu);
      }

  modify_worklin();

  return l;

  }




void DISTR_zinb2_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of delta equation
  // *worktransformlin[0] = delta;
  // *worklin[1] = linear predictor of pi equation
  // *worktransformlin[1] = pi;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = exp(*linpred);

    double pot = pow((*worktransformlin[0])/(mu+(*worktransformlin[0])),(*worktransformlin[0]));

    double nu = (*worktransformlin[0])*((*response)-mu)/((*worktransformlin[0])+mu);

    if (*response==0)
    {
        nu += ((*worktransformlin[0])*(*worktransformlin[1])*mu)/(((*worktransformlin[1])+(1-(*worktransformlin[1]))*pot)*(mu+(*worktransformlin[0])));
    }

    *workingweight = pow(nu,2);

    *workingresponse = *linpred + 1/nu;

//     *workingweight = ( (*worktransformlin[0])*mu*(1-(*worktransformlin[1])) )/( (*worktransformlin[0])+mu ) -
//            ( (*worktransformlin[1])*(1-(*worktransformlin[1]))*pow((*worktransformlin[0]),2)*pow(mu,2)*pot )/( ((*worktransformlin[1])+(1-(*worktransformlin[1]))*pot)*pow((mu+(*worktransformlin[0])),2) );

 //   *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

     if (*response==0)
      {
        like += log((*worktransformlin[1])+(1-(*worktransformlin[1]))*pot);
      }
      else
      {
        like += (*response)*(*linpred)
          - ((*worktransformlin[0])+(*response))*log((*worktransformlin[0])+mu);
      }

      }


  modify_worklin();

  }


void DISTR_zinb2_mu::compute_mu_mult(vector<double *> linpred,double * mu)
  {

  double exp_linmu = exp((*linpred[2]));
  double exp_linpi = exp((*linpred[1]))/(1+ exp((*linpred[1])));

  *mu = (1-exp_linpi)*exp_linmu;

  }


void DISTR_zinb2_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_zinb2_mu::update_end(void)
  {


  // helpmat1 stores (eta_mu)

  double * worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }

*/

//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_gengamma_tau --------------------------
//------------------------------------------------------------------------------


DISTR_gengamma_tau::DISTR_gengamma_tau(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Generalized gamma - tau";
  }


DISTR_gengamma_tau::DISTR_gengamma_tau(const DISTR_gengamma_tau & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_gengamma_tau & DISTR_gengamma_tau::operator=(
                            const DISTR_gengamma_tau & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_gengamma_tau::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_gengamma_tau::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = exp(eta_sigma);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double tau = exp((*linpred));

  double l;

     l = log(tau) + (tau*(*worktransformlin[0]))*log((*response)) -pow((((*worktransformlin[0])/(*worktransformlin[1]))*(*response)),tau) +
     	 tau*(*worktransformlin[0])*log((*worktransformlin[0])/(*worktransformlin[1]));


  modify_worklin();

  return l;

  }

void DISTR_gengamma_tau::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = exp(eta_sigma);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);


  if (counter==0)
    {
    set_worklin();
    }

    double tau = exp((*linpred));


    double nu = 1 + (*worktransformlin[0])*tau*log((*response)) -
				pow(((*worktransformlin[0])/(*worktransformlin[1]))*(*response),tau)*tau*log(((*worktransformlin[0])/(*worktransformlin[1]))*(*response)) +
				tau*(*worktransformlin[0])*log((*worktransformlin[0])/(*worktransformlin[1]));

	double exp_linsigma_plus1 = ((*worktransformlin[0])+1);

    *workingweight = (*worktransformlin[0])*((randnumbers::trigamma_exact(exp_linsigma_plus1))+pow((randnumbers::digamma_exact(exp_linsigma_plus1)),2));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += log(tau) + (tau*(*worktransformlin[0]))*log((*response)) -pow((((*worktransformlin[0])/(*worktransformlin[1]))*(*response)),tau) +
     	 tau*(*worktransformlin[0])*log((*worktransformlin[0])/(*worktransformlin[1]));

      }

  modify_worklin();

  }


void DISTR_gengamma_tau::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (tau): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gengamma_tau::update_end(void)
  {

  // helpmat1 stores tau

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


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_gengamma_sigma ---------------------------
//------------------------------------------------------------------------------


DISTR_gengamma_sigma::DISTR_gengamma_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Generalized gamma - sigma";
  }


DISTR_gengamma_sigma::DISTR_gengamma_sigma(const DISTR_gengamma_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_gengamma_sigma & DISTR_gengamma_sigma::operator=(
                            const DISTR_gengamma_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_gengamma_sigma::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_gengamma_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of tau equation
  // *worktransformlin[0] = exp(eta_tau);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sig = exp((*linpred));

  double l;

     l = (sig*(*worktransformlin[0])-1)*log((*response)) -pow(((sig/(*worktransformlin[1]))*(*response)),(*worktransformlin[0])) +
     	 sig*(*worktransformlin[0])*log(sig/(*worktransformlin[1])) -randnumbers::lngamma_exact(sig);


  modify_worklin();

  return l;

  }

void DISTR_gengamma_sigma::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of tau equation
  // *worktransformlin[0] = exp(eta_tau);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

    double sig = exp((*linpred));


    double nu = sig*(*worktransformlin[0])*log((*response)) - (*worktransformlin[0])*pow((sig/(*worktransformlin[1]))*(*response),(*worktransformlin[0])) +
				sig*(*worktransformlin[0])*log(sig/((*worktransformlin[1]))) + sig*(*worktransformlin[0]) - sig*(randnumbers::digamma_exact(sig));



    *workingweight = sig*(sig*randnumbers::trigamma_exact(sig) - 2*(*worktransformlin[0]) + pow((*worktransformlin[0]),2));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  (sig*(*worktransformlin[0]))*log((*response)) -pow(((sig/(*worktransformlin[1]))*(*response)),(*worktransformlin[0])) +
     	 sig*(*worktransformlin[0])*log(sig/(*worktransformlin[1])) -randnumbers::lngamma_exact(sig);

      }

  modify_worklin();

  }


void DISTR_gengamma_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gengamma_sigma::update_end(void)
  {

  // helpmat1 stores sigma

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


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_gengamma_mu ------------------------
//------------------------------------------------------------------------------

void DISTR_gengamma_mu::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight > 0)
        {

        if (*workresp <= 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative/zero response values encountered\n");
          }


        }
      else if (*workweight == 0)
        {
        }
      else
        {
        errors=true;
        errormessages.push_back("ERROR: negative weights encountered\n");
        }

      i++;
      workresp++;
      workweight++;

      }

    }

  }

DISTR_gengamma_mu::DISTR_gengamma_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Generalized gamma - mu";
  }


DISTR_gengamma_mu::DISTR_gengamma_mu(const DISTR_gengamma_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_gengamma_mu & DISTR_gengamma_mu::operator=(
                            const DISTR_gengamma_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_gengamma_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_tau
   // *linpred[1] = eta_sigma
   // *linpred[2] = eta_mu

   if (*weight[2] == 0)
     *deviance=0;
   else
     {
	 double tau = exp(*linpred[0]);
     double sig = exp(*linpred[1]);
     double mu = exp(*linpred[2]);

     double l;

       l = log(tau) + (sig*tau-1)*log((*response[2])) -pow(((sig/mu)*(*response[2])),tau) + sig*tau*log(sig/mu) -randnumbers::lngamma_exact(sig);


    *deviance = -2*l;
    }

  }


double DISTR_gengamma_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_gengamma_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of tau equation
  // *worktransformlin[0] = tau;
  // *worklin[1] = linear predictor of sigma equation
  // *worktransformlin[1] = sigma;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = exp(*linpred);

  double l;

     l =  -pow((((*worktransformlin[1])/mu)*(*response)),(*worktransformlin[0])) - (*worktransformlin[1])*(*worktransformlin[0])*log(mu) ;
  modify_worklin();

  return l;

  }


void DISTR_gengamma_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of tau equation
  // *worktransformlin[0] = tau;
  // *worklin[1] = linear predictor of sigma equation
  // *worktransformlin[1] = sigma;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double mu = exp((*linpred));

    double exponent = (*worktransformlin[0]);

    double nu = (*worktransformlin[0])*pow(((*worktransformlin[1])/mu)*(*response),exponent) - ((*worktransformlin[1]))*((*worktransformlin[0]));

    *workingweight = (*worktransformlin[1])*pow((*worktransformlin[0]),2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -pow((((*worktransformlin[1])/mu)*(*response)),(*worktransformlin[0])) - (*worktransformlin[1])*(*worktransformlin[0])*log(mu);

      }


  modify_worklin();

  }


void DISTR_gengamma_mu::compute_mu_mult(vector<double *> linpred,double * mu)
  {

  double exp_linmu = exp((*linpred[2]));
  double exp_linsigma = exp((*linpred[1]));
  double exp_lintau = exp((*linpred[0]));
  double help = exp_linsigma+1/exp_lintau;
  *mu = (randnumbers::lngamma_exact(help))*exp_linmu/(exp_linsigma*(randnumbers::lngamma_exact(exp_linsigma)));

  }


void DISTR_gengamma_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gengamma_mu::update_end(void)
  {


  // helpmat1 stores (eta_mu)

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
//    double t = 0;
    }

  }



//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_gamma_sigma ---------------------------
//------------------------------------------------------------------------------


DISTR_gamma_sigma::DISTR_gamma_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "gamma - sigma";
  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma";
  }


DISTR_gamma_sigma::DISTR_gamma_sigma(const DISTR_gamma_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_gamma_sigma & DISTR_gamma_sigma::operator=(
                            const DISTR_gamma_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_gamma_sigma::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_gamma_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sig = exp((*linpred));

  double l;

     l = sig*log(sig) - sig*log((*worktransformlin[0])) - randnumbers::lngamma_exact(sig) + (sig-1)*log(*response) - (sig/(*worktransformlin[0]))*(*response);


  modify_worklin();

  return l;

  }

void DISTR_gamma_sigma::compute_iwls_wweightschange_weightsone(
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

    double sig = exp((*linpred));


    double nu = sig*log(sig) +  sig - sig*log((*worktransformlin[0])) - sig*(randnumbers::digamma_exact(sig)) +
				sig*log((*response)) - (sig/(*worktransformlin[0]))*(*response);



    *workingweight = sig*(sig*randnumbers::trigamma_exact(sig) - 1);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += sig*log(sig) - sig*log((*worktransformlin[0])) - randnumbers::lngamma_exact(sig) + (sig-1)*log(*response) - (sig/(*worktransformlin[0]))*(*response);

      }

  modify_worklin();

  }


void DISTR_gamma_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gamma_sigma::update_end(void)
  {

  // helpmat1 stores sigma

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


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_gamma_mu ------------------------
//------------------------------------------------------------------------------

void DISTR_gamma_mu::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight > 0)
        {

        if (*workresp <= 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative/zero response values encountered\n");
          }


        }
      else if (*workweight == 0)
        {
        }
      else
        {
        errors=true;
        errormessages.push_back("ERROR: negative weights encountered\n");
        }

      i++;
      workresp++;
      workweight++;

      }

    }

  }

DISTR_gamma_mu::DISTR_gamma_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "gamma - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  }


DISTR_gamma_mu::DISTR_gamma_mu(const DISTR_gamma_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_gamma_mu & DISTR_gamma_mu::operator=(
                            const DISTR_gamma_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_gamma_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_sigma
   // *linpred[1] = eta_mu

   if (*weight[1] == 0)
     *deviance=0;
   else
     {
     double sig = exp(*linpred[0]);
     double mu = exp(*linpred[1]);

     double l;

       l = sig*log(sig) - sig*log(mu) - randnumbers::lngamma_exact(sig) + (sig-1)*log(*response[1]) - (sig/mu)*(*response[1]);


    *deviance = -2*l;
    }

  }


double DISTR_gamma_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_gamma_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = sigma;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = exp(*linpred);

  double l;

     l =  - (*worktransformlin[0])*log(mu) - ((*worktransformlin[0])/mu)*(*response);
  modify_worklin();

  return l;

  }


void DISTR_gamma_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = sigma;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double mu = exp(*linpred);

    double nu = -(*worktransformlin[0]) + ((*worktransformlin[0])/mu)*(*response);

    *workingweight = (*worktransformlin[0]);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  - (*worktransformlin[0])*log(mu) - ((*worktransformlin[0])/mu)*(*response);

      }


  modify_worklin();

  }


void DISTR_gamma_mu::compute_mu_mult(vector<double *> linpred,double * mu)
  {

  *mu = exp((*linpred[predstart_mumult+1]));

  }


void DISTR_gamma_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gamma_mu::update_end(void)
  {


  // helpmat1 stores (eta_mu)

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
//    double t = 0;
    }

  }



//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_lognormal_sigma2 ----------------------
//------------------------------------------------------------------------------


DISTR_lognormal_sigma2::DISTR_lognormal_sigma2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Lognormal - sigma2";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma2";

  }


DISTR_lognormal_sigma2::DISTR_lognormal_sigma2(const DISTR_lognormal_sigma2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_lognormal_sigma2 & DISTR_lognormal_sigma2::operator=(
                            const DISTR_lognormal_sigma2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_lognormal_sigma2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_lognormal_sigma2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sigma_2 = exp((*linpred));

  double l;

     l = -0.5*log(sigma_2)-pow((log((*response))-(*worklin[0])),2)/(2*sigma_2);


  modify_worklin();

  return l;

  }

void DISTR_lognormal_sigma2::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu)/(1+exp(eta_mu));

  if (counter==0)
    {
    set_worklin();
    }

    double sigma_2 = exp(*linpred);


    double nu = -0.5 + (pow((log(*response)-(*worklin[0])),2))/(2*sigma_2);



    *workingweight = 0.5;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(sigma_2)-pow((log((*response))-(*worklin[0])),2)/(2*sigma_2);

      }

  modify_worklin();

  }


void DISTR_lognormal_sigma2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_lognormal_sigma2::update_end(void)
  {

  // helpmat1 stores sigma2

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


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_lognormal_mu ------------------------
//------------------------------------------------------------------------------
void DISTR_lognormal_mu::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight > 0)
        {

        if (*workresp <= 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative/zero response values encountered\n");
          }


        }
      else if (*workweight == 0)
        {
        }
      else
        {
        errors=true;
        errormessages.push_back("ERROR: negative weights encountered\n");
        }

      i++;
      workresp++;
      workweight++;

      }

    }

  }


DISTR_lognormal_mu::DISTR_lognormal_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Lognormal - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  }


DISTR_lognormal_mu::DISTR_lognormal_mu(const DISTR_lognormal_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_lognormal_mu & DISTR_lognormal_mu::operator=(
                            const DISTR_lognormal_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_lognormal_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_sigma2
   // *linpred[1] = eta_mu

   if (*weight[1] == 0)
     *deviance=0;
   else
     {
     double sigma_2 = exp(*linpred[0]);
     double mu = (*linpred[1]);

     double l;

       l = -0.5*log(2*PI)-0.5*log(sigma_2)-log((*response[0]))-pow((log((*response[0]))-mu),2)/(2*sigma_2);


    *deviance = -2*l;
    }

  }


double DISTR_lognormal_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_lognormal_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double l;

     l = -pow((log((*response))-mu),2)/(2*(*worktransformlin[0]));

  modify_worklin();

  return l;

  }


void DISTR_lognormal_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double mu = (*linpred);

    double nu = (log(*response)-mu)/(*worktransformlin[0]);

    *workingweight = 1/(*worktransformlin[0]);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow((log((*response))-mu),2)/(2*(*worktransformlin[0]));

      }


  modify_worklin();

  }


void DISTR_lognormal_mu::compute_mu_mult(vector<double *> linpred,double * mu)
  {
  double s = exp(*linpred[predstart_mumult]);
  *mu = exp((*linpred[predstart_mumult+1])+s/2);
  }


void DISTR_lognormal_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_lognormal_mu::update_end(void)
  {


  // helpmat1 stores (eta_mu)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = (*worklin);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_beta_sigma2 ---------------------------
//------------------------------------------------------------------------------


DISTR_beta_sigma2::DISTR_beta_sigma2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Beta - sigma2";
  }


DISTR_beta_sigma2::DISTR_beta_sigma2(const DISTR_beta_sigma2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_beta_sigma2 & DISTR_beta_sigma2::operator=(
                            const DISTR_beta_sigma2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_beta_sigma2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_beta_sigma2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu)/(1+exp(eta_mu));

  if (counter==0)
    {
    set_worklin();
    }

  double exp_lin = exp((*linpred));
  double sigma_2 = exp_lin/(1+exp_lin);
  double help = (1-sigma_2)/sigma_2;
  double mu_help = (*worktransformlin[0])*help;
  double one_minus_mu_help = (1-(*worktransformlin[0]))*help;

  double l;

     l = mu_help*log(*response) +
         one_minus_mu_help*log(1-(*response)) -
         randnumbers::lngamma_exact(mu_help) - randnumbers::lngamma_exact(one_minus_mu_help) +
         randnumbers::lngamma_exact(help);


  modify_worklin();

  return l;

  }

void DISTR_beta_sigma2::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu)/(1+exp(eta_mu));

  if (counter==0)
    {
    set_worklin();
    }

    double exp_lin = exp(*linpred);

    double sigma_2 = exp_lin/(1+exp_lin);

    double help = (1-sigma_2)/sigma_2;
    double mu_help = (*worktransformlin[0])*help;
    double one_minus_mu_help = (1-(*worktransformlin[0]))*help;

    double nu = -help*( -(*worktransformlin[0])*randnumbers::digamma_exact(mu_help)
                       - (1-(*worktransformlin[0]))*randnumbers::digamma_exact(one_minus_mu_help) +
                       randnumbers::digamma_exact(help) + (*worktransformlin[0])*log(*response)
                       + (1-(*worktransformlin[0]))*log(1-(*response)) );



    *workingweight =  pow(help,2)*( pow((1-(*worktransformlin[0])),2)*randnumbers::trigamma_exact(one_minus_mu_help)+
                                   pow(((*worktransformlin[0])),2)*randnumbers::trigamma_exact(mu_help) -
                                   randnumbers::trigamma_exact(help) );

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += mu_help*log(*response) +
         one_minus_mu_help*log(1-(*response)) -
         randnumbers::lngamma_exact(mu_help) - randnumbers::lngamma_exact(one_minus_mu_help) +
         randnumbers::lngamma_exact(help);

      }

  modify_worklin();

  }


void DISTR_beta_sigma2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): logit\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_beta_sigma2::update_end(void)
  {

  // helpmat1 stores (1-sigma2)/sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double exp_lin;
  double sigma_2;
  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    exp_lin = exp(*worklin);
    sigma_2 = exp_lin/(1+exp_lin);
    *pmu = (1-sigma_2)/sigma_2;
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_beta_mu -----------------------------
//------------------------------------------------------------------------------

void DISTR_beta_mu::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight > 0)
        {

        if (*workresp <= 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative/zero response values encountered\n");
          }

        if (*workresp >= 1)
          {
          errors=true;
          errormessages.push_back("ERROR: response values greater or equal to one encountered\n");
          }



        }
      else if (*workweight == 0)
        {
        }
      else
        {
        errors=true;
        errormessages.push_back("ERROR: negative weights encountered\n");
        }

      i++;
      workresp++;
      workweight++;

      }

    }

  }


DISTR_beta_mu::DISTR_beta_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Beta - mu";
  }


DISTR_beta_mu::DISTR_beta_mu(const DISTR_beta_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_beta_mu & DISTR_beta_mu::operator=(
                            const DISTR_beta_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_beta_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_sigma2
   // *linpred[1] = eta_mu

   if (*weight[1] == 0)
     *deviance=0;
   else
     {
     double sigma_2 = exp(*linpred[0])/(1+exp(*linpred[0]));
     double mu = exp(*linpred[1])/(1+exp(*linpred[1]));
     double help = (1-sigma_2)/sigma_2;
     double one_minus_mu_help = (1-mu)*help;
     double mu_help = mu*help;

     double l;

       l = (mu_help-1)*log(*response[1]) +
			(one_minus_mu_help-1)*log(1-(*response[1]))-
			randnumbers::lngamma_exact(mu_help)-
			randnumbers::lngamma_exact(one_minus_mu_help)+
			randnumbers::lngamma_exact(help);


    *deviance = -2*l;
    }

  }


double DISTR_beta_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_beta_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = (1-sigma2)/sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double exp_lin = exp((*linpred));
  double mu = exp_lin/(1+exp_lin);
  double mu_worktrans = mu*(*worktransformlin[0]);
  double one_minus_mu_worktrans = (1-mu)*(*worktransformlin[0]);

  double l;

     l = mu*(*worktransformlin[0])*log(*response) +
		 (1-mu)*(*worktransformlin[0])*log(1-(*response))-
		 randnumbers::lngamma_exact(mu_worktrans)-
		 randnumbers::lngamma_exact(one_minus_mu_worktrans);

  modify_worklin();

  return l;

  }


void DISTR_beta_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = (1-sigma2)/sigma2;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double exp_lin = exp(*linpred);

    double mu = exp_lin/(1+exp_lin);

    double mu_worktrans = mu*(*worktransformlin[0]);

    double one_minus_mu_worktrans = (1-mu)*(*worktransformlin[0]);

    double nu = mu*one_minus_mu_worktrans*log((*response))-mu*one_minus_mu_worktrans*log((1-(*response)))+
				mu*one_minus_mu_worktrans*(randnumbers::digamma_exact(one_minus_mu_worktrans)-randnumbers::digamma_exact(mu_worktrans));

    *workingweight = pow((*worktransformlin[0]),2)*pow(mu,2)*pow((1-mu),2)*(randnumbers::trigamma_exact(one_minus_mu_worktrans)+randnumbers::trigamma_exact(mu_worktrans));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += mu_worktrans*log((*response))+one_minus_mu_worktrans*log(1-(*response))-
		        randnumbers::lngamma_exact(mu_worktrans)-randnumbers::lngamma_exact((one_minus_mu_worktrans));

      }


  modify_worklin();

  }


void DISTR_beta_mu::compute_mu_mult(vector<double *> linpred,double * mu)
  {
  double exp_lin = exp(*linpred[1]);
  *mu = exp_lin/(1+exp_lin);
  }


void DISTR_beta_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): logit\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_beta_mu::update_end(void)
  {


  // helpmat1 stores exp(eta_mu)/(1+exp(eta_mu))

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double exp_lin;
  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    exp_lin = exp(*worklin);
    *pmu = exp_lin/(1+exp_lin);
//    double t = 0;
    }

  }





} // end: namespace MCMC



