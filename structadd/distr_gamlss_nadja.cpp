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
  optionsp->out("  Response function (sigma2): logit\n");
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
  optionsp->out("  Response function (mu): logit\n");
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



