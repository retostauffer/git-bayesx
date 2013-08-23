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
//--------------------------- CLASS: DISTR_betainf0_nu --------------------------
//------------------------------------------------------------------------------


DISTR_betainf0_nu::DISTR_betainf0_nu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,0,w)
  {
  family = "Beta Zero Inflated - nu";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "tau";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_betainf0_nu::DISTR_betainf0_nu(const DISTR_betainf0_nu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_betainf0_nu & DISTR_betainf0_nu::operator=(
                            const DISTR_betainf0_nu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_betainf0_nu::compute_deviance_mult(vector<double *> response,
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

   if (*weight[2] == 0)
     *deviance=0;
   else
     {
     double nup = exp(*linpred[2]);
     double sigma_2 = exp(*linpred[1])/(1+exp(*linpred[1]));
     double mu = exp(*linpred[0])/(1+exp(*linpred[0]));
     double help = (1-sigma_2)/sigma_2;
     double one_minus_mu_help = (1-mu)*help;
     double mu_help = mu*help;
     double one_nup = (1+nup);

     double l;

     if ((*response[2])==0)
     {
         l = log(nup) - log(one_nup);
     }
      else
       l = (mu_help-1)*log(*response[1]) +
			(one_minus_mu_help-1)*log(1-(*response[1]))-
			randnumbers::lngamma_exact(mu_help)-
			randnumbers::lngamma_exact(one_minus_mu_help)+
			randnumbers::lngamma_exact(help)- log(one_nup);


    *deviance = -2*l;
    }

  }


double DISTR_betainf0_nu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

 double DISTR_betainf0_nu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_betainf0_nu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
//    double a =  (*param[2])*(*param[3]);
//    double b = (*param[2])*(1-(*param[3]));
 //   double frac = 1 + (*param[0]) + (*param[1]);
    return 0;
//
//    return ( ((*param[0])+(*param[1]))/frac+((1-(*param[0])-(*param[1]))/frac)*randnumbers::incomplete_beta(a,b,(*response[3])) );
    }


double DISTR_betainf0_nu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of nu equation
  // *worktransformlin[0] = exp(eta_nu);

  if (counter==0)
    {
    set_worklin();
    }

  double nup = exp((*linpred));

  double l;

    if ((*response)==0)
     {
         l = log(nup) - log(1+nup);
     }
     else
     {
        l =  - log(1+nup);
     }

  modify_worklin();

  return l;

  }


void DISTR_betainf0_nu::compute_iwls_wweightschange_weightsone(
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

    double nup = exp(*linpred);

    double hilfs = 1+ nup;

    double nu = -nup/hilfs;

    if ((*response)==0)
    {
        nu += 1;
    }

    *workingweight = nup/(pow(hilfs,2));

    *workingresponse = *linpred + nu/(*workingweight);

        if (compute_like)
      {

        if ((*response)==0)
        {
            like += log(nup) - log(1+nup);
        }
        else
        {
            like -=   log(1+nup);
        }
        }


  modify_worklin();

  }


void DISTR_betainf0_nu::compute_mu_mult(vector<double *> linpred,double * mu)
  {
  double exp_lin_nup = exp(*linpred[predstart_mumult+2]);
  double exp_lin_mu = exp(*linpred[predstart_mumult]);
  double hilfs = (1+exp_lin_nup);
  *mu = (1-(exp_lin_nup)/hilfs)*(exp_lin_mu/(1+exp_lin_mu));
  }


void DISTR_betainf0_nu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (nu): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_betainf0_nu::update_end(void)
  {


  // helpmat1 stores exp(eta_nu)

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
//------------------------- CLASS: DISTR_t_df ----------------------------------
//------------------------------------------------------------------------------


DISTR_t_df::DISTR_t_df(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "t-distribution - degrees of freedom";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "df";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_t_df::DISTR_t_df(const DISTR_t_df & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_t_df & DISTR_t_df::operator=(
                            const DISTR_t_df & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_t_df::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_t_df::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = exp(eta_sigma2);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

    double degf = exp(*linpred);
    double arg = (degf+1)/2;
    double arg1 = (degf)/2;
    double l;

     l = (randnumbers::lngamma_exact(arg))-(randnumbers::lngamma_exact(arg1))-
            0.5*log(degf) - (arg)*log(1+pow((*response)-(*worklin[1]),2)/((*worktransformlin[0])*degf));


  modify_worklin();

  return l;

  }

void DISTR_t_df::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = exp(eta_sigma2);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);


  if (counter==0)
    {
    set_worklin();
    }

    double degf = exp((*linpred));
    double denom = (*worktransformlin[0])*degf+pow((*response)-(*worklin[1]),2);
    double frac = pow((*response)-(*worktransformlin[1]),2)/((*worktransformlin[0])*degf);
    double arg = (degf+1)/2;
    double arg1 = degf/2;
    double nu = -0.5 + arg*pow(((*response)-(*worklin[1])),2)/denom + 0.5*degf*(randnumbers::digamma_exact(arg)-randnumbers::digamma_exact(arg1)-log(1+frac));

  //  *workingweight = pow(nu,2);
    *workingweight =  - 0.25*pow(degf,2)*(randnumbers::trigamma_exact(arg)-randnumbers::trigamma_exact(arg1)) - degf/((degf+1)) + degf/(2*(degf+3));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += (randnumbers::lngamma_exact(arg))-(randnumbers::lngamma_exact(arg1))-
            0.5*log(degf) - (arg)*log(1+frac);

      }

  modify_worklin();

  }


void DISTR_t_df::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (df): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_t_df::update_end(void)
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
//------------------------- CLASS: DISTR_t_sigma2 ---------------------------
//------------------------------------------------------------------------------


DISTR_t_sigma2::DISTR_t_sigma2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "t-distribution - sigma2";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma2";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_t_sigma2::DISTR_t_sigma2(const DISTR_t_sigma2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_t_sigma2 & DISTR_t_sigma2::operator=(
                            const DISTR_t_sigma2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_t_sigma2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_t_sigma2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of n equation
  // *worktransformlin[0] = exp(eta_n);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sigma_2 = exp(*linpred);

  double l;

     l =  - 0.5*log(sigma_2) - (((*worktransformlin[0])+1)/(2))*log(1+pow((*response)-(*worklin[1]),2)/(sigma_2*(*worktransformlin[0])));


  modify_worklin();

  return l;

  }

void DISTR_t_sigma2::compute_iwls_wweightschange_weightsone(
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

    double sig = exp(*linpred);
    double denom = sig*(*worktransformlin[0])+pow(((*response)-(*worklin[1])),2);

    double nu = -0.5 + (((*worktransformlin[0])+1)/(2)*pow(((*response)-(*worklin[1])),2))/denom;

  //  *workingweight = pow(nu,2);

 //  *workingweight = 0.5*sig*(*worktransformlin[0])*((*worktransformlin[0])+1)*pow((*response)-(*worktransformlin[1]),2)/pow(denom,2);

    *workingweight = (*worktransformlin[0])/(2*((*worktransformlin[0])+3));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  - 0.5*log(sig) - (((*worktransformlin[0])+1)/(2))*log(1+pow((*response)-(*worklin[1]),2)/(sig*(*worktransformlin[0])));

      }

  modify_worklin();

  }


void DISTR_t_sigma2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_t_sigma2::update_end(void)
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
//--------------------------- CLASS: DISTR_t_mu -----------------------------
//------------------------------------------------------------------------------


void DISTR_t_mu::check_errors(void)
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


DISTR_t_mu::DISTR_t_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "t - mu";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
   // linpredminlimit=-10;
  //linpredmaxlimit=15;
  }


DISTR_t_mu::DISTR_t_mu(const DISTR_t_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_t_mu & DISTR_t_mu::operator=(
                            const DISTR_t_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_t_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_df
   // *linpred[1] = eta_sigma2
   // *linpred[2] = eta_mu

   if (*weight[2] == 0)
     *deviance=0;
   else
     {
	 double degf = exp(*linpred[0]);
     double sigma_2 = exp(*linpred[1]);
     double mu = (*linpred[2]);
     double hilfs = (degf+1)/2;
     double hilfs2 = degf/2;

     double l;

       l = randnumbers::lngamma_exact((hilfs)) -log(sqrt(PI)) -randnumbers::lngamma_exact(hilfs2) -0.5*log(degf) - 0.5*log(sigma_2) -
            hilfs*log(1+(pow((*response[2])-mu,2))/(sigma_2*degf));


    *deviance = -2*l;
    }

  }


double DISTR_t_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

  void DISTR_t_mu::compute_param(const double * linpred,double * param)
  {
  *param = (*linpred);
  }

 double DISTR_t_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_t_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    double a = 0.5*(*param[0]);
    double b = 0.5;
    double x = (*param[0])/(pow((((*response[2])-(*param[2]))/pow((*param[1]),0.5)),2)+(*param[0]));
    return 0;
//    return ( 1- 0.5*randnumbers::incomplete_beta(a,b,x));
    }

double DISTR_t_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of degf equation
  // *worktransformlin[0] = degf;
  // *worklin[1] = linear predictor of sigma_2 equation
  // *worktransformlin[1] = sigma_2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);
  double hilfs =((*worktransformlin[0])+1)/2;
  double l;

     l =  - hilfs*log(1+(pow((*response)-mu,2))/((*worktransformlin[1])*(*worktransformlin[0])));

  modify_worklin();

  return l;

  }


void DISTR_t_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of degf equation
  // *worktransformlin[0] = degf;
  // *worklin[1] = linear predictor of sigma_2 equation
  // *worktransformlin[1] = sigma_2;


  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double mu = (*linpred);
    double hilfs = ((*worktransformlin[0])+1)/2;
    double denom1 = (*worktransformlin[1])*(*worktransformlin[0])-pow((*response)-mu,2);
    double denom2 = (*worktransformlin[1])*(*worktransformlin[0])+pow((*response)-mu,2);

    double nu = ((*worktransformlin[0])+1)*((*response)-mu)/denom2;

   // *workingweight = pow(nu,2);
    *workingweight = -1/(*worktransformlin[1]) + 2*((*worktransformlin[0])+2)/((*worktransformlin[1])*((*worktransformlin[0])+3));
    //*workingweight = (((*worktransformlin[0])+1)*denom1)/(denom2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=   - hilfs*log(1+(pow((*response)-mu,2))/((*worktransformlin[1])*(*worktransformlin[0])));

      }


  modify_worklin();

  }


void DISTR_t_mu::compute_mu_mult(vector<double *> linpred,double * mu)
  {

   if(exp((*linpred[predstart_mumult]))>1)
   {
       *mu = (*linpred[predstart_mumult+2]);
   } else
   {
       *mu = 0;
   }



  }


void DISTR_t_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_t_mu::update_end(void)
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
    *pmu = (*worklin);
//    double t = 0;
    }

  }



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
    linpredminlimit=-10;
  linpredmaxlimit=15;

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

     l = -0.5*log(sigma_2)-pow((((*response))-(*worktransformlin[0])),2)/(2*(*response)*pow((*worktransformlin[0]),2)*sigma_2);


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


    double nu = -0.5 + (pow(((*response)-(*worktransformlin[0])),2))/(2*(*response)*(pow((*worktransformlin[0]),2))*sigma_2);



    *workingweight = 0.5;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(sigma_2)-pow((((*response))-(*worktransformlin[0])),2)/(2*(*response)*pow((*worktransformlin[0]),2)*sigma_2);

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
    linpredminlimit=-10;
  linpredmaxlimit=15;
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

 double DISTR_invgaussian_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_invgaussian_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    double arg = pow(1/((*param[0])*(*response[1])),0.5);
    double arg0 = (*response[1])/((*param[1]));
    double arg1 = arg*(arg0-1);
    double arg2 = -arg*(arg0+1);
    double arg3 = 2/((*param[0])*(*param[0]));

    return (randnumbers::Phi2(arg1)+exp(arg3)*randnumbers::Phi2(arg2));
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
    outpredictor = true;
  outexpectation = false;
  predictor_name = "mu";
    linpredminlimit=-10;
  linpredmaxlimit=10;
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

void DISTR_betainf_mu::compute_param(const double * linpred,double * param)
  {
  double arg = exp(*linpred);
  *param = arg/(1+arg);
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
    outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma2";
    linpredminlimit=-10;
  linpredmaxlimit=10;
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

void DISTR_betainf_sigma2::compute_param(const double * linpred,double * param)
  {
   double arg = exp(*linpred);
  *param = arg/(1+arg);
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
    outpredictor = true;
  outexpectation = false;
  predictor_name = "nu";
    linpredminlimit=-10;
  linpredmaxlimit=15;
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
    outpredictor = true;
  outexpectation = true;
  predictor_name = "tau";
    linpredminlimit=-10;
  linpredmaxlimit=15;
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
       l = (mu_help-1)*log(*response[2]) +
			(one_minus_mu_help-1)*log(1-(*response[2]))-
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

 double DISTR_betainf_tau::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_betainf_tau::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    double a =  (*param[2])*(*param[3]);
    double b = (*param[2])*(1-(*param[3]));
    double frac = 1 + (*param[0]) + (*param[1]);
    return 0;
//
//    return ( ((*param[0])+(*param[1]))/frac+((1-(*param[0])-(*param[1]))/frac)*randnumbers::incomplete_beta(a,b,(*response[3])) );
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
  double exp_lin_tau = exp(*linpred[predstart_mumult+3]);
  double exp_lin_nup = exp(*linpred[predstart_mumult+2]);
  double exp_lin_mu = exp(*linpred[predstart_mumult]);
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
    outpredictor = true;
  outexpectation = false;
  predictor_name = "p";
    linpredminlimit=-10;
  linpredmaxlimit=15;
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
    outpredictor = true;
  outexpectation = true;
  predictor_name = "b";
    linpredminlimit=-10;
  linpredmaxlimit=15;
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

 double DISTR_pareto_b::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_pareto_b::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    double b = *param[1];
    double p = *param[0];
    double res = 1- pow(b,p)*pow((*response[1])+b,-p);
    return ( res );
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

   double p = exp((*linpred[predstart_mumult]));
   double b = exp((*linpred[predstart_mumult+1]));
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
    outpredictor = true;
  outexpectation = false;
  predictor_name = "p";
    linpredminlimit=-10;
  linpredmaxlimit=15;
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
    outpredictor = true;
  outexpectation = false;
  predictor_name = "b";
    linpredminlimit=-10;
  linpredmaxlimit=15;
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


   // *workingweight = (pow((*worktransformlin[1]),2)*pow((*worktransformlin[0]),2)*((*worktransformlin[0])+1))/((*worktransformlin[0])+2);
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
    outpredictor = true;
  outexpectation = true;
  predictor_name = "a";
    linpredminlimit=-10;
  linpredmaxlimit=15;
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

 double DISTR_dagum_a::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_dagum_a::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {

    return ( pow((1+pow((*response[1])/(*param[1]),-(*param[2]))),-(*param[0])) );
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

  double exp_lin_a = exp((*linpred[predstart_mumult+2]));
  double exp_lin_b = exp((*linpred[predstart_mumult+1]));
  double exp_lin_p = exp((*linpred[predstart_mumult]));
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
  family = "Weibull - alpha";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "alpha";
    linpredminlimit=-10;
  linpredmaxlimit=15;
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
//1 + a*log(y/b) - a*log(y/b)*(y/b)^a
//(1 + a^2*log(y/b)^2*(y/b)^a)

  //  *workingweight = 1 + pow(sig,2)*pow((log((*response)/(*worktransformlin[0]))),2)*hilfs1;
    *workingweight =  1.823681;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {
//log(a) + (a)*log(y/b) - (y/b)^a
        like +=  (sig)*log(*response) - hilfs1 -sig*log((*worktransformlin[0])) + log(sig);

      }

  modify_worklin();

  }


void DISTR_weibull_alpha::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (alpha): exponential\n");
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
  family = "Weibull - lambda";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "lambda";
    linpredminlimit=-10;
  linpredmaxlimit=15;
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

 double DISTR_weibull_lambda::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_weibull_lambda::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {

    return ( 1 - exp(-pow((*response[1])*(*param[1]),(*param[0]))) );
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
 //   double hilfs2 = (*worktransformlin[0])+1;

    double nu = (*worktransformlin[0])*( hilfs1-1 );

  //  *workingweight = (*worktransformlin[0])*randnumbers::gamma_exact(hilfs2)*((*worktransformlin[0])-1)+(*worktransformlin[0]);

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

   double hilfs = 1+1/exp((*linpred[predstart_mumult+1]));
  *mu = (exp((*linpred[predstart_mumult])))*randnumbers::gamma_exact(hilfs);

  }


void DISTR_weibull_lambda::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (lambda): exponential\n");
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
//------------------------- CLASS: DISTR_gengamma_tau --------------------------
//------------------------------------------------------------------------------


DISTR_gengamma_tau::DISTR_gengamma_tau(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Generalized gamma - tau";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "tau";
    linpredminlimit=-10;
  linpredmaxlimit=15;
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

    *workingweight = (*worktransformlin[0])*((randnumbers::trigamma_exact(exp_linsigma_plus1))+pow((randnumbers::digamma_exact(exp_linsigma_plus1)),2))+1;

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
    outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma";
    linpredminlimit=-10;
  linpredmaxlimit=15;
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
    outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
    linpredminlimit=-10;
  linpredmaxlimit=15;
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

 double DISTR_gengamma_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_gengamma_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {

    return ( 0 );
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

  double exp_linmu = exp((*linpred[predstart_mumult+2]));
  double exp_linsigma = exp((*linpred[predstart_mumult+1]));
  double exp_lintau = exp((*linpred[predstart_mumult]));
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
    linpredminlimit=-10;
  linpredmaxlimit=15;
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
    linpredminlimit=-10;
  linpredmaxlimit=15;
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

 double DISTR_gamma_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_gamma_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {

    return ( 0 );
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
//------------------------- CLASS: DISTR_lognormal2_sigma ----------------------
//------------------------------------------------------------------------------


DISTR_lognormal2_sigma::DISTR_lognormal2_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Lognormal - sigma";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma2";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_lognormal2_sigma::DISTR_lognormal2_sigma(const DISTR_lognormal2_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_lognormal2_sigma & DISTR_lognormal2_sigma::operator=(
                            const DISTR_lognormal2_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_lognormal2_sigma::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_lognormal2_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sigma_2 = pow(exp((*linpred)),2);

  double l;

     l = -0.5*log(sigma_2)-pow((log((*response))-(*worklin[0])),2)/(2*sigma_2);


  modify_worklin();

  return l;

  }

void DISTR_lognormal2_sigma::compute_iwls_wweightschange_weightsone(
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

    double sigma_2 = pow(exp(*linpred),2);


    double nu = -1 + (pow((log(*response)-(*worklin[0])),2))/(sigma_2);



    *workingweight = 2;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(sigma_2)-pow((log((*response))-(*worklin[0])),2)/(2*sigma_2);

      }

  modify_worklin();

  }


void DISTR_lognormal2_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_lognormal2_sigma::update_end(void)
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
    *pmu = pow(exp(*worklin),2);
    }

  }



//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_lognormal2_mu ------------------------
//------------------------------------------------------------------------------
void DISTR_lognormal2_mu::check_errors(void)
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


DISTR_lognormal2_mu::DISTR_lognormal2_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Lognormal - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
 //   linpredminlimit=-10;
 // linpredmaxlimit=15;
  }


DISTR_lognormal2_mu::DISTR_lognormal2_mu(const DISTR_lognormal2_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_lognormal2_mu & DISTR_lognormal2_mu::operator=(
                            const DISTR_lognormal2_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_lognormal2_mu::compute_deviance_mult(vector<double *> response,
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
     double sigma_2 = pow(exp(*linpred[0]),2);
     double mu = (*linpred[1]);

     double l;

       l = -0.5*log(2*PI)-0.5*log(sigma_2)-log((*response[0]))-pow((log((*response[0]))-mu),2)/(2*sigma_2);


    *deviance = -2*l;
    }

  }


double DISTR_lognormal2_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_lognormal2_mu::compute_param(const double * linpred,double * param)
  {
  *param = (*linpred);
  }

 double DISTR_lognormal2_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_lognormal2_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    double arg = (log(*response[1])-(*param[1]))/(*param[0]) ;

    return (randnumbers::Phi2(arg));
    }

double DISTR_lognormal2_mu::loglikelihood_weightsone(double * response,
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


void DISTR_lognormal2_mu::compute_iwls_wweightschange_weightsone(
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


void DISTR_lognormal2_mu::compute_mu_mult(vector<double *> linpred,double * mu)
  {
  double s2 = pow(exp(*linpred[predstart_mumult]),2);
  *mu = exp((*linpred[predstart_mumult+1])+s2/2);
  }


void DISTR_lognormal2_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_lognormal2_mu::update_end(void)
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
    linpredminlimit=-10;
  linpredmaxlimit=15;

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
   // linpredminlimit=-10;
  //linpredmaxlimit=15;
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

void DISTR_lognormal_mu::compute_param(const double * linpred,double * param)
  {
  *param = (*linpred);
  }

 double DISTR_lognormal_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_lognormal_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    double arg = (log(*response[1])-(*param[1]))/pow((*param[0]),0.5) ;

    return (randnumbers::Phi2(arg));
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
//------------------------- CLASS: DISTR_normal2_sigma ----------------------
//------------------------------------------------------------------------------


DISTR_normal2_sigma::DISTR_normal2_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "normal2 - sigma";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_normal2_sigma::DISTR_normal2_sigma(const DISTR_normal2_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_normal2_sigma & DISTR_normal2_sigma::operator=(
                            const DISTR_normal2_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_normal2_sigma::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_normal2_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sigma_2 = pow(exp((*linpred)),2);

  double l;

     l = -0.5*log(sigma_2)-pow((((*response))-(*worklin[0])),2)/(2*sigma_2);


  modify_worklin();

  return l;

  }

void DISTR_normal2_sigma::compute_iwls_wweightschange_weightsone(
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

    double sigma_2 = pow(exp(*linpred),2);


    double nu = -1 + (pow(((*response)-(*worklin[0])),2))/(sigma_2);



    *workingweight = 2;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(sigma_2)-pow((((*response))-(*worklin[0])),2)/(2*sigma_2);

      }

  modify_worklin();

  }


void DISTR_normal2_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_normal2_sigma::update_end(void)
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
    *pmu = pow(exp(*worklin),2);
    }

  }

  //------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_normal2_mu ------------------------
//------------------------------------------------------------------------------
void DISTR_normal2_mu::check_errors(void)
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


DISTR_normal2_mu::DISTR_normal2_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "normal2 - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
//    linpredminlimit=-10;
 // linpredmaxlimit=15;
  }


DISTR_normal2_mu::DISTR_normal2_mu(const DISTR_normal2_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_normal2_mu & DISTR_normal2_mu::operator=(
                            const DISTR_normal2_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_normal2_mu::compute_deviance_mult(vector<double *> response,
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
     double sigma_2 = pow(exp(*linpred[0]),2);
     double mu = (*linpred[1]);

     double l;

       l = -0.5*log(2*PI)-0.5*log(sigma_2)-pow((((*response[0]))-mu),2)/(2*sigma_2);


    *deviance = -2*l;
    }

  }


double DISTR_normal2_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_normal2_mu::compute_param(const double * linpred,double * param)
  {
  *param = (*linpred);
  }

 double DISTR_normal2_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_normal2_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    double arg = ((*response[1])-(*param[1]))/(*param[0]) ;

    return (randnumbers::Phi2(arg));
    }

double DISTR_normal2_mu::loglikelihood_weightsone(double * response,
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

     l = -pow((((*response))-mu),2)/(2*(*worktransformlin[0]));

  modify_worklin();

  return l;

  }


void DISTR_normal2_mu::compute_iwls_wweightschange_weightsone(
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

    double nu = ((*response)-mu)/(*worktransformlin[0]);

    *workingweight = 1/(*worktransformlin[0]);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow((((*response))-mu),2)/(2*(*worktransformlin[0]));

      }


  modify_worklin();

  }


void DISTR_normal2_mu::compute_mu_mult(vector<double *> linpred,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+1]));
  }


void DISTR_normal2_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_normal2_mu::update_end(void)
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
//------------------------- CLASS: DISTR_normal_sigma2 ---------------------------
//--------------------------------------------------------------------------------


DISTR_normal_sigma2::DISTR_normal_sigma2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "normal - sigma2";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma2";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_normal_sigma2::DISTR_normal_sigma2(const DISTR_normal_sigma2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_normal_sigma2 & DISTR_normal_sigma2::operator=(
                            const DISTR_normal_sigma2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_normal_sigma2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


double DISTR_normal_sigma2::loglikelihood_weightsone(double * response,
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

     l = -0.5*log(sigma_2)-pow((((*response))-(*worklin[0])),2)/(2*sigma_2);


  modify_worklin();

  return l;

  }

void DISTR_normal_sigma2::compute_iwls_wweightschange_weightsone(
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


    double nu = -0.5 + (pow(((*response)-(*worklin[0])),2))/(2*sigma_2);



    *workingweight = 0.5;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(sigma_2)-pow((((*response))-(*worklin[0])),2)/(2*sigma_2);

      }

  modify_worklin();

  }


void DISTR_normal_sigma2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_normal_sigma2::update_end(void)
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
//--------------------------- CLASS: DISTR_normal_mu ---------------------------
//------------------------------------------------------------------------------
void DISTR_normal_mu::check_errors(void)
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


DISTR_normal_mu::DISTR_normal_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "normal - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
 //   linpredminlimit=-10;
 // linpredmaxlimit=15;
  }


DISTR_normal_mu::DISTR_normal_mu(const DISTR_normal_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_normal_mu & DISTR_normal_mu::operator=(
                            const DISTR_normal_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_normal_mu::compute_deviance_mult(vector<double *> response,
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

       l = -0.5*log(2*PI)-0.5*log(sigma_2)-pow((((*response[0]))-mu),2)/(2*sigma_2);


    *deviance = -2*l;
    }

  }


double DISTR_normal_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_normal_mu::compute_param(const double * linpred,double * param)
  {
  *param = (*linpred);
  }

 double DISTR_normal_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_normal_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    double arg = ((*response[1])-(*param[1]))/pow((*param[0]),0.5) ;

    return (randnumbers::Phi2(arg));
    }

//double DISTR_normal_mu::compute_quantile_residual_mult(vector<double *> response,
//                                             vector<double *> param,
//                                             vector<double *> weight,
//                                             vector<datamatrix *> aux)
//  {
//  double u_est;
//  u_est = cdf_mult(response,param,weight,aux);
//  double res_est = randnumbers::invPhi2(u_est);
//  return res_est;
//  }

double DISTR_normal_mu::loglikelihood_weightsone(double * response,
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

     l = -pow((((*response))-mu),2)/(2*(*worktransformlin[0]));

  modify_worklin();

  return l;

  }


void DISTR_normal_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;


  if (counter==0)
    {
    set_worklin();
    }

    double mu = (*linpred);

    double nu = ((*response)-mu)/(*worktransformlin[0]);

    *workingweight = 1/(*worktransformlin[0]);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow((((*response))-mu),2)/(2*(*worktransformlin[0]));

      }


  modify_worklin();

  }


void DISTR_normal_mu::compute_mu_mult(vector<double *> linpred,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+1]));
  }


void DISTR_normal_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_normal_mu::update_end(void)
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
    outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma2";
    linpredminlimit=-10;
  linpredmaxlimit=10;
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

void DISTR_beta_sigma2::compute_param(const double * linpred,double * param)
  {
   double arg = exp((*linpred));
  *param = arg/(1+arg);
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
    outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  linpredminlimit=-10;
  linpredmaxlimit=10;
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

 /*   double a = sigma_2*mu;
    double b = sigma_2*(1-mu);
    std::ofstream out;
  // helpmat1.prettyPrint(out);
    for (int i=0; i<(a+b); i++) {
    out.open ("C:\\Urs\\nchoosek.raw", std::ofstream::out | std::ofstream::app);
    out << a ;
    out << " " ;
    out << b ;
    out << " " ;
    out << i ;
    out << " " ;
    out << randnumbers::n_choose_k(a+b-1,i) ;
    out << " " ;
    out << randnumbers::incomplete_beta(a,b,(*response[1])) << endl;
    out.close();
    }*/

    }

  }


double DISTR_beta_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

 double DISTR_beta_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

void DISTR_beta_mu::compute_param(const double * linpred,double * param)
  {
   double arg = exp((*linpred));
  *param = arg/(1+arg);
  }

double DISTR_beta_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    double a =  (*param[0])*(*param[1]);
    double b = (*param[0])*(1-(*param[1]));

     return 0;
//   return ( randnumbers::incomplete_beta(a,b,(*response[1])) );
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
  double exp_lin = exp(*linpred[predstart_mumult+1]);
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


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_cloglog -----------------------------
//------------------------------------------------------------------------------

void DISTR_cloglog::check_errors(void)
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

        if ((*workresp!= 0) | (*workresp!= 1) )
          {
          errors=true;
          errormessages.push_back("ERROR: response has to be equal to zero or one\n");
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


DISTR_cloglog::DISTR_cloglog(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,0,w)
  {
  family = "binomial-cloglog";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  linpredminlimit=-7.5;
  linpredmaxlimit=2.2;
  }


DISTR_cloglog::DISTR_cloglog(const DISTR_cloglog & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_cloglog & DISTR_cloglog::operator=(
                            const DISTR_cloglog & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_cloglog::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response


   if (*weight[0] == 0)
     *deviance=0;
   else
     {
     double oneminuspi = exp(-exp((*linpred[0])));
     double pi = 1 - oneminuspi;

     double l;

     if ((*response[0])==0) {
        l = log(oneminuspi);
     } else {
        l = log(pi);
     }


    *deviance = -2*l;
    }

  }


double DISTR_cloglog::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_cloglog::compute_param(const double * linpred,double * param)
  {
   double arg = exp((*linpred));
  *param = 1-exp(-arg);
  }

double DISTR_cloglog::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
        double Fy = 0;
        if (*response[0]>0) {
            Fy = 1;
            if((*response[0]<1)&(*response[0]>0)) {
                Fy = 1-(*param[0]);
            }
        }
    return Fy;
    }

double DISTR_cloglog::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_cloglog::compute_quantile_residual_mult(vector<double *> response,
                                         vector<double *> param,
                                         vector<double *> weight,
                                          vector<datamatrix *> aux)
    {
    double u_est;
    if(*response[0]==0) {
        u_est = randnumbers::uniform_ab(0,(*param[0]));
    }
    else {
        u_est = randnumbers::uniform();
    }
    double res_est = randnumbers::invPhi2(u_est);
    return res_est;
    }


double DISTR_cloglog::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = (1-sigma2)/sigma2;

  if (counter==0)
    {
    set_worklin();
    }

     double oneminuspi = exp(-exp((*linpred)));
     double pi = 1 - oneminuspi;

  double l;

   if ((*response)==0) {
        l = log(oneminuspi);
     } else {
        l = log(pi);
     }

  modify_worklin();

  return l;

  }


void DISTR_cloglog::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

     double oneminuspi = exp(-exp((*linpred)));
     double pi = 1 - oneminuspi;


    double nu = exp((*linpred))*oneminuspi/pi;

    if ((*response)==0) {
        nu -= exp((*linpred))/pi;
    }

    *workingweight = pow(exp((*linpred)),2)*oneminuspi/pi;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

   if (*response==0)
      like += log(oneminuspi);
    else // response == 1
      like += log(pi);

      }


  modify_worklin();

  }


void DISTR_cloglog::compute_mu_mult(vector<double *> linpred,double * mu)
  {
  double exp_lin = exp(*linpred[predstart_mumult]);
  *mu = 1-exp(-exp_lin);
  }


void DISTR_cloglog::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu):cloglog\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_cloglog::update_end(void)
  {


  // helpmat1 stores 1-exp(-exp_lin)

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
    *pmu = 1-exp(-exp_lin);
//    double t = 0;
    }

  }

//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_dirichlet -----------------------------
//------------------------------------------------------------------------------

void DISTR_dirichlet::check_errors(void)
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

        if ((*workresp>1) | (*workresp<0) )
          {
          errors=true;
          errormessages.push_back("ERROR: response has to be between zero and one\n");
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


DISTR_dirichlet::DISTR_dirichlet(GENERAL_OPTIONS * o,
                                           const datamatrix & r, int & nrc, unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,nrc-1,w)
  {
  family = "dirichlet";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "alpha";
  linpredminlimit=-10;
  linpredmaxlimit=15;
  nrcat = nrc;
  pos = p;
  }


DISTR_dirichlet::DISTR_dirichlet(const DISTR_dirichlet & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
    nrcat = nd.nrcat;
    pos = nd.pos;
  }


const DISTR_dirichlet & DISTR_dirichlet::operator=(
                            const DISTR_dirichlet & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  nrcat = nd.nrcat;
  pos = nd.pos;
  return *this;
  }


void DISTR_dirichlet::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response


   if (*weight[0] == 0)
     *deviance=0;
   else
     {
     double sum_alpha = 0;
     double sum_log_gamma = 0;
     double sum_rest = 0;

     unsigned i;
     for(i=0;i<nrcat;i++) {
        double hilfs = exp(*linpred[i]);
        sum_log_gamma += randnumbers::lngamma_exact(hilfs);
        sum_alpha += hilfs;
        sum_rest += (hilfs-1)*log(*response[i]);
     }

     double l = -sum_log_gamma + randnumbers::lngamma_exact(sum_alpha) + sum_rest;


    *deviance = -2*l;
    }

  }


 double DISTR_dirichlet::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

//  void DISTR_dirichlet::compute_param(const double * linpred,double * param)
//  {
//  double el = exp(*linpred);
//  *param = el/(1+el);
//  }

 double DISTR_dirichlet::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

 double DISTR_dirichlet::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_dirichlet::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of alpha_i equation

  if (counter==0)
    {
    set_worklin();
    }
     double alpha_current = exp(*linpred);
     double sum_alpha = alpha_current;

     unsigned i;
     for(i=0;i<(nrcat-1);i++) {
        sum_alpha += (*worktransformlin[i]);
     }

     double l = -randnumbers::lngamma_exact(alpha_current) + randnumbers::lngamma_exact(sum_alpha) + alpha_current*log(*response);

  modify_worklin();

  return l;

  }


void DISTR_dirichlet::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

     double alpha_current = exp(*linpred);
     double sum_alpha = alpha_current;

     unsigned i;
     for(i=0;i<(nrcat-1);i++) {
        sum_alpha += (*worktransformlin[i]);
     }

    double nu = alpha_current*( -randnumbers::digamma_exact(alpha_current) +randnumbers::digamma_exact(sum_alpha) + log(*response) );

    *workingweight = pow(alpha_current,2)*( randnumbers::trigamma_exact(alpha_current) - randnumbers::trigamma_exact(sum_alpha) );

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

      like += -randnumbers::lngamma_exact(alpha_current) + randnumbers::lngamma_exact(sum_alpha) + alpha_current*log(*response);

      }


  modify_worklin();

  }


void DISTR_dirichlet::compute_mu_mult(vector<double *> linpred,double * mu)
  {
     double alpha_current = exp(*linpred[predstart_mumult+pos]);
     double sum_alpha = 0;

     unsigned i;
     for(i=0;i<nrcat;i++) {
        double hilfs = exp(*linpred[predstart_mumult+i]);
        sum_alpha += hilfs;
        }

     *mu = alpha_current/sum_alpha;
  }


void DISTR_dirichlet::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (alpha): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_dirichlet::update_end(void)
  {


  // helpmat1 stores 1-exp(-exp_lin)

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
//------------------------- CLASS: DISTR_bivt_df ----------------------------
//------------------------------------------------------------------------------


DISTR_bivt_df::DISTR_bivt_df(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,5,w)
  {
  family = "bivariate t - degrees of freedom";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "df";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_bivt_df::DISTR_bivt_df(const DISTR_bivt_df & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivt_df & DISTR_bivt_df::operator=(
                            const DISTR_bivt_df & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_bivt_df::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


void DISTR_bivt_df::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_bivt_df::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs-1)
    {
    response2p++;
    }

  }


double DISTR_bivt_df::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of sigma_2 equation
  // *worktransformlin[1] = sigma_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of mu_2 equation
  // *worktransformlin[3] = mu_2;
  // *worklin[4] = linear predictor of mu_1 equation
  // *worktransformlin[4] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }

  double degf = exp((*linpred));

  double oneminusrho2 = 1- pow((*worktransformlin[0]),2);
  double np2d2 = (degf+2)/2;
  double nd2 = degf/2;
  double X_1 = ((*response)-(*worklin[4]))/(*worktransformlin[2]);
  double X_2 = ((*response2p)-(*worklin[3]))/(*worktransformlin[1]);
  double l;


     l = randnumbers::lngamma_exact(np2d2)-randnumbers::lngamma_exact(nd2)-log(degf)
     -np2d2*log( 1 + (1/(degf*oneminusrho2))*( pow(X_1,2)-2*(*worktransformlin[0])*X_1*X_2+pow(X_2,2) ) );


  modify_worklin();

  return l;

  }

void DISTR_bivt_df::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of sigma_2 equation
  // *worktransformlin[1] = sigma_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of mu_2 equation
  // *worktransformlin[3] = mu_2;
  // *worklin[4] = linear predictor of mu_1 equation
  // *worktransformlin[4] = mu_1;

  if (counter==0)
    {
    set_worklin();
    }

  double degf = exp((*linpred));

  double oneminusrho2 = 1- pow((*worktransformlin[0]),2);
  double np2d2 = (degf+2)/2;
  double nd2 = degf/2;
  double X_1 = ((*response)-(*worklin[4]))/(*worktransformlin[2]);
  double X_2 = ((*response2p)-(*worklin[3]))/(*worktransformlin[1]);
  double nenner_C = 1 + (1/(degf*oneminusrho2))*( pow(X_1,2)-2*(*worktransformlin[0])*X_1*X_2+pow(X_2,2) );

    double nu = nd2*( randnumbers::digamma_exact(np2d2)-randnumbers::digamma_exact(nd2)-log(nenner_C) ) - 1 + np2d2*(nenner_C-1)/nenner_C;

    *workingweight =  - pow(nd2,2)*( randnumbers::trigamma_exact(np2d2) - randnumbers::trigamma_exact(nd2) ) - degf/(degf+2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  randnumbers::lngamma_exact(np2d2)-randnumbers::lngamma_exact(nd2)-log(degf)
     -np2d2*log( nenner_C );

      }

  modify_worklin();

  }


void DISTR_bivt_df::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (df): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivt_df::update_end(void)
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
//------------------------- CLASS: DISTR_bivt_rho ------------------------------
//------------------------------------------------------------------------------


DISTR_bivt_rho::DISTR_bivt_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,5,w)
  {
  family = "bivariate t - rho";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "rho";
    linpredminlimit=-100;
  linpredmaxlimit=100;

  }


DISTR_bivt_rho::DISTR_bivt_rho(const DISTR_bivt_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivt_rho & DISTR_bivt_rho::operator=(
                            const DISTR_bivt_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_bivt_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivt_rho::compute_param(const double * linpred,double * param)
  {
   double arg = (*linpred);
  *param = arg/pow(1+pow(arg,2),0.5);
  }

void DISTR_bivt_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_bivt_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs-1)
    {
    response2p++;
    }

  }



double DISTR_bivt_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

 // *worklin[0] = linear predictor of df equation
  // *worktransformlin[0] = df;
  // *worklin[1] = linear predictor of sigma_2 equation
  // *worktransformlin[1] = sigma_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of mu_2 equation
  // *worktransformlin[3] = mu_2;
  // *worklin[4] = linear predictor of mu_1 equation
  // *worktransformlin[4] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }


  double rho;

  if (*linpred <= -100)
    rho  = -0.99995;
  else if (*linpred >= 100)
    rho  = 0.99995;
  else
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);

  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;
  double X_1 = ((*response)-(*worklin[4]))/(*worktransformlin[2]);
  double X_2 = ((*response2p)-(*worklin[3]))/(*worktransformlin[1]);
  double nenner_C = 1 + (1/((*worktransformlin[0])*oneminusrho2))*( pow(X_1,2)-2*rho*X_1*X_2+pow(X_2,2) );
  double l;


     l = -0.5*log(oneminusrho2)-(((*worktransformlin[0])+2)/2)*log( nenner_C );



  modify_worklin();

  return l;

  }

void DISTR_bivt_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

 // *worklin[0] = linear predictor of df equation
  // *worktransformlin[0] = df;
  // *worklin[1] = linear predictor of sigma_2 equation
  // *worktransformlin[1] = sigma_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of mu_2 equation
  // *worktransformlin[3] = mu_2;
  // *worklin[4] = linear predictor of mu_1 equation
  // *worktransformlin[4] = mu_1;

  if (counter==0)
    {
    set_worklin();
    }

  double rho;

  if (*linpred <= -100)
    rho  = -0.99995;
  else if (*linpred >= 100)
    rho  = 0.99995;
  else
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);

  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;
  double X_1 = ((*response)-(*worklin[4]))/(*worktransformlin[2]);
  double X_2 = ((*response2p)-(*worklin[3]))/(*worktransformlin[1]);
  double nenner_C = 1 + (1/((*worktransformlin[0])*oneminusrho2))*( pow(X_1,2)-2*rho*X_1*X_2+pow(X_2,2) );

  double nu = oneminusrho2*(*linpred) - (((*worktransformlin[0])+2)/((*worktransformlin[0])*nenner_C))*( (*linpred)*(pow(X_1,2)+pow(X_2,2)) -
                                                                                                      (pow(1/oneminusrho2,0.5)+rho*(*linpred))*X_1*X_2 );



    *workingweight = 1-pow(rho2,2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(oneminusrho2)-(((*worktransformlin[0])+2)/2)*log( nenner_C );

      }

  modify_worklin();

  }


void DISTR_bivt_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (rho):\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivt_rho::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = (*worklin)/pow((1+pow((*worklin),2)),0.5);
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_bivt_sigma ----------------------------
//------------------------------------------------------------------------------


DISTR_bivt_sigma::DISTR_bivt_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,5,w)
  {
  family = "bivariate t - sigma";

  pos = p;
  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_bivt_sigma::DISTR_bivt_sigma(const DISTR_bivt_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivt_sigma & DISTR_bivt_sigma::operator=(
                            const DISTR_bivt_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_bivt_sigma::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


void DISTR_bivt_sigma::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_bivt_sigma::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs-1)
    {
    response2p++;
    }

  }


double DISTR_bivt_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of df equation
  // *worktransformlin[0] = df;
  // *worklin[1] = linear predictor of rho equation
  // *worktransformlin[1] = rho;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of mu_2 equation
  // *worktransformlin[3] = mu_2;
  // *worklin[4] = linear predictor of mu_1 equation
  // *worktransformlin[4] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }

  double sig_current = exp((*linpred));

  double oneminusrho2 = 1- pow((*worktransformlin[1]),2);
  double X_1 = ((*response)-(*worklin[3]))/(sig_current);
  double X_2 = ((*response2p)-(*worklin[4]))/(*worktransformlin[2]);
  double nenner_C = 1 + (1/((*worktransformlin[0])*oneminusrho2))*( pow(X_1,2)-2*(*worktransformlin[1])*X_1*X_2+pow(X_2,2) );
  double l;


     l = -log(sig_current)
     -(((*worktransformlin[0])+2)/2)*log( nenner_C );


  modify_worklin();

  return l;

  }

void DISTR_bivt_sigma::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of df equation
  // *worktransformlin[0] = df;
  // *worklin[1] = linear predictor of rho equation
  // *worktransformlin[1] = rho;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of mu_2 equation
  // *worktransformlin[3] = mu_2;
  // *worklin[4] = linear predictor of mu_1 equation
  // *worktransformlin[4] = mu_1;

  if (counter==0)
    {
    set_worklin();
    }

  double sig_current = exp((*linpred));

  double oneminusrho2 = 1- pow((*worktransformlin[1]),2);
  double X_1 = ((*response)-(*worklin[3]))/(sig_current);
  double X_2 = ((*response2p)-(*worklin[4]))/(*worktransformlin[2]);
  double nenner_C = 1 + (1/((*worktransformlin[0])*oneminusrho2))*( pow(X_1,2)-2*(*worktransformlin[1])*X_1*X_2+pow(X_2,2) );

    double nu = -1 - (((*worktransformlin[0])+2)/(nenner_C*(*worktransformlin[0])*oneminusrho2))*((*worktransformlin[1])*X_1*X_2-pow(X_1,2));



    *workingweight = 1+1/oneminusrho2 ;
//- 2/(*worktransformlin[0]) - 2*pow((*worktransformlin[1]),2)/((*worktransformlin[0])*oneminusrho2)
    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -log(sig_current)-(((*worktransformlin[0])+2)/2)*log( nenner_C );

      }

  modify_worklin();

  }


void DISTR_bivt_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivt_sigma::update_end(void)
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
//--------------------------- CLASS: DISTR_bivt_mu -----------------------------
//------------------------------------------------------------------------------
void DISTR_bivt_mu::check_errors(void)
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


DISTR_bivt_mu::DISTR_bivt_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,5,w)
  {
  pos =p;
  family = "bivariate t - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
//    linpredminlimit=-10;
//  linpredmaxlimit=15;
  }


DISTR_bivt_mu::DISTR_bivt_mu(const DISTR_bivt_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivt_mu & DISTR_bivt_mu::operator=(
                            const DISTR_bivt_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


void DISTR_bivt_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[2] = *response[4] = first component of two dimensional reponse
   // *linpred[0] = eta_df
   // *linpred[1] = eta_rho
   // *linpred[2] = eta_sigma_2
   // *linpred[3] = eta_mu_2
   // *linpred[4] = eta_sigma_1
   // *linpred[5] = eta_mu_1

   if (*weight[5] == 0)
     *deviance=0;
   else
     {
     double rho = (*linpred[1])/pow(1+pow((*linpred[1]),2),0.5);
     double degf = exp(*linpred[0]);
     double np2d2 = (degf+2)/2;
     double nd2 = degf/2;
     double sigma_2 = exp(*linpred[2]);
     double mu_2 = (*linpred[4]);
     double sigma_1 = exp(*linpred[3]);
     double mu_1 = (*linpred[5]);
     double hilfs1 = 1-pow(rho,2);
     double X_1 = ((*response[5])-mu_1)/sigma_1;
     double X_2 = ((*response[4])-mu_2)/sigma_2;
     double nenner_C = 1+(1/(degf*hilfs1))*(pow(X_1,2)-2*rho*X_1*X_2+pow(X_2,2));
     double l;

       l = randnumbers::lngamma_exact(np2d2)-randnumbers::lngamma_exact(nd2)-log(degf)-log(PI)-log(sigma_1)-log(sigma_2)-0.5*log(hilfs1)-
           np2d2*log(nenner_C);


    *deviance = -2*l;
    }

  }


double DISTR_bivt_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivt_mu::compute_param(const double * linpred,double * param)
  {
  *param = (*linpred);
  }

void DISTR_bivt_mu::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_bivt_mu::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs-1)
    {
    response2p++;
    }

  }



double DISTR_bivt_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of df equation
  // *worktransformlin[0] = df;
  // *worklin[1] = linear predictor of rho equation
  // *worktransformlin[1] = rho;
  // *worklin[2] = linear predictor of mu_1 equation
  // *worktransformlin[2] = mu_1;
  // *worklin[3] = linear predictor of sigma_2 equation
  // *worktransformlin[3] = sigma_2;
  // *worklin[4] = linear predictor of sigma_1 equation
  // *worktransformlin[4] = sigma_1;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double oneminusrho2 = 1- pow((*worktransformlin[1]),2);
  double X_1 = ((*response)-mu)/((*worktransformlin[3]));
  double X_2 = ((*response2p)-(*worklin[2]))/(*worktransformlin[4]);
  double nenner_C = 1 + (1/((*worktransformlin[0])*oneminusrho2))*( pow(X_1,2)-2*(*worktransformlin[1])*X_1*X_2+pow(X_2,2) );
  double l;


     l = -(((*worktransformlin[0])+2)/2)*log( nenner_C );

  modify_worklin();

  return l;

  }


void DISTR_bivt_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of df equation
  // *worktransformlin[0] = df;
  // *worklin[1] = linear predictor of rho equation
  // *worktransformlin[1] = rho;
  // *worklin[2] = linear predictor of mu_1 equation
  // *worktransformlin[2] = mu_1;
  // *worklin[3] = linear predictor of sigma_2 equation
  // *worktransformlin[3] = sigma_2;
  // *worklin[4] = linear predictor of sigma_1 equation
  // *worktransformlin[4] = sigma_1;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double oneminusrho2 = 1- pow((*worktransformlin[1]),2);
  double X_1 = ((*response)-mu)/((*worktransformlin[3]));
  double X_2 = ((*response2p)-(*worklin[2]))/(*worktransformlin[4]);
  double nenner_C = 1 + (1/((*worktransformlin[0])*oneminusrho2))*( pow(X_1,2)-2*(*worktransformlin[1])*X_1*X_2+pow(X_2,2) );

    double nu =(((*worktransformlin[0])+2)/((*worktransformlin[0])*oneminusrho2*(*worktransformlin[3])*nenner_C))*(X_1-(*worktransformlin[1])*X_2);

    *workingweight = (1)/(oneminusrho2*pow((*worktransformlin[3]),2)) ;
    //- ((*worktransformlin[0])+2)/(pow((*worktransformlin[0]),2)*oneminusrho2*pow((*worktransformlin[3]),2));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -(((*worktransformlin[0])+2)/2)*log( nenner_C );

      }


  modify_worklin();

  }


void DISTR_bivt_mu::compute_mu_mult(vector<double *> linpred,double * mu)
  {
    if(exp((*linpred[predstart_mumult]))>2)
    {
        *mu = ((*linpred[predstart_mumult+4+pos]));
    } else
    {
        *mu = 0;
    }

  }


void DISTR_bivt_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivt_mu::update_end(void)
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
//------------------------- CLASS: DISTR_bivnormal_rho -------------------------
//------------------------------------------------------------------------------


DISTR_bivnormal_rho::DISTR_bivnormal_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "bivariate normal - rho";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "rho";
    linpredminlimit=-100;
  linpredmaxlimit=100;

  }


DISTR_bivnormal_rho::DISTR_bivnormal_rho(const DISTR_bivnormal_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivnormal_rho & DISTR_bivnormal_rho::operator=(
                            const DISTR_bivnormal_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_bivnormal_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivnormal_rho::compute_param(const double * linpred,double * param)
  {
   double arg = (*linpred);
  *param = arg/pow(1+pow(arg,2),0.5);
  }

void DISTR_bivnormal_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_bivnormal_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs-1)
    {
    response2p++;
    }

  }



double DISTR_bivnormal_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma_2 equation
  // *worktransformlin[0] = sigma_2;
  // *worklin[1] = linear predictor of sigma_1 equation
  // *worktransformlin[1] = sigma_1;
  // *worklin[2] = linear predictor of mu_2 equation
  // *worktransformlin[2] = mu_2;
  // *worklin[3] = linear predictor of mu_1 equation
  // *worktransformlin[3] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }
  double rho;

  if (*linpred <= -100)
    rho  = -0.99995;
  else if (*linpred >= 100)
    rho  = 0.99995;
  else
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);

  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;
  double l;


     l = -0.5*log(oneminusrho2) -(1/(2*oneminusrho2))*( pow((((*response))-(*worklin[3])),2)/pow((*worktransformlin[1]),2) -
                                 2*rho*(((*response)-(*worklin[3]))/((*worktransformlin[1])))*(((*response2p)-(*worklin[2]))/((*worktransformlin[0])))
                                +  pow((((*response2p))-(*worklin[2])),2)/pow((*worktransformlin[0]),2) );


  modify_worklin();

  return l;

  }

void DISTR_bivnormal_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma_2 equation
  // *worktransformlin[0] = sigma_2;
  // *worklin[1] = linear predictor of sigma_1 equation
  // *worktransformlin[1] = sigma_1;
  // *worklin[2] = linear predictor of mu_2 equation
  // *worktransformlin[2] = mu_2;
  // *worklin[3] = linear predictor of mu_1 equation
  // *worktransformlin[3] = mu_1;

  if (counter==0)
    {
    set_worklin();
    }

  double rho;
  double hilfs;

  if (*linpred <= -100) {
    rho  = -0.99995;
    hilfs = 100.05;
  }
  else if (*linpred >= 100) {
    rho  = 0.99995;
    hilfs = 100.05;
  }
  else {
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);
    hilfs = pow((1+pow((*linpred),2)),0.5);
  }


  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;


    double nu = oneminusrho2*(*linpred) - (*linpred)*( pow((((*response))-(*worktransformlin[3])),2)/pow((*worktransformlin[1]),2)
                                                      +  pow((((*response2p))-(*worktransformlin[2])),2)/pow((*worktransformlin[0]),2) )
                +(hilfs+rho*(*linpred))*( (((*response)-(*worktransformlin[3]))/((*worktransformlin[1])))*(((*response2p)-(*worktransformlin[2]))/((*worktransformlin[0]))) );



    *workingweight = 1-pow(rho2,2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(oneminusrho2) -(1/(2*oneminusrho2))*( pow((((*response))-(*worklin[3])),2)/pow((*worktransformlin[1]),2) -
                                 2*rho*(((*response)-(*worklin[3]))/((*worktransformlin[1])))*(((*response2p)-(*worklin[2]))/((*worktransformlin[0])))
                                +  pow((((*response2p))-(*worklin[2])),2)/pow((*worktransformlin[0]),2) );

      }

  modify_worklin();

  }


void DISTR_bivnormal_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (rho):\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivnormal_rho::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = (*worklin)/pow((1+pow((*worklin),2)),0.5);
    }

  }

//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_bivnormal_sigma ----------------------
//------------------------------------------------------------------------------


DISTR_bivnormal_sigma::DISTR_bivnormal_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "bivariate normal - sigma";

  pos = p;
  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_bivnormal_sigma::DISTR_bivnormal_sigma(const DISTR_bivnormal_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivnormal_sigma & DISTR_bivnormal_sigma::operator=(
                            const DISTR_bivnormal_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_bivnormal_sigma::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }


void DISTR_bivnormal_sigma::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_bivnormal_sigma::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs-1)
    {
    response2p++;
    }

  }


double DISTR_bivnormal_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of sigma_2 equation
  // *worktransformlin[3] = sigma_2;


  if (counter==0)
    {
    set_worklin();
    }

  double sigma = exp((*linpred));

  double rho2 = pow((*worktransformlin[0]),2);
  double oneminusrho2 = 1- pow((*worktransformlin[0]),2);
  double l;

//hier ist jetzt das problem, dass ich die aktuelle response gleichung brauche und die von dem mu2 oder sigma2

     l = -log(sigma) -(1/(2*oneminusrho2))*( pow((((*response))-(*worktransformlin[2])),2)/pow(sigma,2) -
                                 2*(*worktransformlin[0])*(((*response)-(*worktransformlin[2]))/(sigma))*(((*response2p)-(*worktransformlin[3]))/((*worktransformlin[1]))) );


  modify_worklin();

  return l;

  }

void DISTR_bivnormal_sigma::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of sigma_2 equation
  // *worktransformlin[3] = sigma_2;

  if (counter==0)
    {
    set_worklin();
    }

    double sigma = exp((*linpred));

    double rho2 = pow((*worktransformlin[0]),2);
    double oneminusrho2 = 1- rho2;


    double nu = -1 + (1/oneminusrho2)*(pow(((*response)-(*worklin[2])),2))/pow(sigma,2)
                - ((*worktransformlin[0])/oneminusrho2)*(((*response)-(*worktransformlin[2]))/(sigma))*(((*response2p)-(*worklin[3]))/((*worktransformlin[1])));



    *workingweight = 1+1/oneminusrho2;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -log(sigma) -(1/(2*oneminusrho2))*( pow((((*response))-(*worktransformlin[2])),2)/pow(sigma,2) -
                                 2*(*worktransformlin[0])*(((*response)-(*worktransformlin[2]))/(sigma))*(((*response2p)-(*worklin[3]))/((*worktransformlin[1]))) );

      }

  modify_worklin();

  }


void DISTR_bivnormal_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivnormal_sigma::update_end(void)
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
//--------------------------- CLASS: DISTR_bivnormal_mu ------------------------
//------------------------------------------------------------------------------
void DISTR_bivnormal_mu::check_errors(void)
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


DISTR_bivnormal_mu::DISTR_bivnormal_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  pos =p;
  family = "bivariate normal - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
//    linpredminlimit=-10;
//  linpredmaxlimit=15;
  }


DISTR_bivnormal_mu::DISTR_bivnormal_mu(const DISTR_bivnormal_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivnormal_mu & DISTR_bivnormal_mu::operator=(
                            const DISTR_bivnormal_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


void DISTR_bivnormal_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[2] = *response[4] = first component of two dimensional reponse
   // *linpred[0] = eta_rho
   // *linpred[1] = eta_sigma_2
   // *linpred[2] = eta_mu_2
   // *linpred[3] = eta_sigma_1
   // *linpred[4] = eta_mu_1

   if (*weight[4] == 0)
     *deviance=0;
   else
     {
     double rho = (*linpred[0])/pow(1+pow((*linpred[0]),2),0.5);
     double sigma_2 = exp(*linpred[1]);
     double mu_2 = (*linpred[3]);
     double sigma_1 = exp(*linpred[2]);
     double mu_1 = (*linpred[4]);
     double hilfs1 = 1-pow(rho,2);
     double l;

       l = -log(2*PI)-log(sigma_1)-log(sigma_2)-0.5*log(hilfs1)-(1/(2*hilfs1))*( pow((((*response[4]))-mu_1),2)/pow(sigma_1,2) -
                                                                                2*rho*(((*response[4])-mu_1)/(sigma_1))*(((*response[3])-mu_2)/(sigma_2))
                                                                                + pow((((*response[3]))-mu_2),2)/pow(sigma_2,2) );


    *deviance = -2*l;
    }

  }


double DISTR_bivnormal_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivnormal_mu::compute_param(const double * linpred,double * param)
  {
  *param = (*linpred);
  }

void DISTR_bivnormal_mu::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_bivnormal_mu::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs-1)
    {
    response2p++;
    }

  }



double DISTR_bivnormal_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of sigma_2 equation
  // *worktransformlin[3] = sigma_2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);
  double rho2 = pow((*worktransformlin[0]),2);
  double oneminusrho2 = 1- rho2;
  double l;


     l = -(1/(2*oneminusrho2))*( pow((((*response))-mu),2)/pow((*worktransformlin[2]),2) -
                                 2*(*worktransformlin[0])*(((*response)-mu)/((*worktransformlin[2])))*(((*response2p)-(*worktransformlin[1]))/((*worktransformlin[3]))) );

  modify_worklin();

  return l;

  }


void DISTR_bivnormal_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of sigma_2 equation
  // *worktransformlin[3] = sigma_2;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double mu = (*linpred);
    double rho2 = pow((*worktransformlin[0]),2);
   double oneminusrho2 = 1- rho2;


    double nu = (1/(oneminusrho2))*( (((*response))-mu)/pow((*worktransformlin[2]),2) -
                                 ((*worktransformlin[0])/(*worktransformlin[2]))*(((*response2p)-(*worktransformlin[1]))/((*worktransformlin[3]))) );

    *workingweight = 1/(oneminusrho2*pow((*worktransformlin[2]),2));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -(1/(2*oneminusrho2))*( pow((((*response))-mu),2)/pow((*worktransformlin[2]),2) -
                                 2*(*worktransformlin[0])*(((*response)-mu)/((*worktransformlin[2])))*(((*response2p)-(*worktransformlin[1]))/((*worktransformlin[3]))) );

      }


  modify_worklin();

  }


void DISTR_bivnormal_mu::compute_mu_mult(vector<double *> linpred,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+3+pos]));
  }


void DISTR_bivnormal_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivnormal_mu::update_end(void)
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
//------------------------- CLASS: DISTR_bivprobit_rho ----------------------
//------------------------------------------------------------------------------


DISTR_bivprobit_rho::DISTR_bivprobit_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "bivariate probit - rho";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "rho";
    linpredminlimit=-100;
  linpredmaxlimit=100;
 // responseorig = response;

  }


DISTR_bivprobit_rho::DISTR_bivprobit_rho(const DISTR_bivprobit_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
 // responseorig = nd.responseorig;
//  response2 = nd.response2;
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  response1p = nd.response1p;
  workingresponse1p = nd.workingresponse1p;
  }


const DISTR_bivprobit_rho & DISTR_bivprobit_rho::operator=(
                            const DISTR_bivprobit_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  responseorig = nd.responseorig;
 // response2 = nd.response2;
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
    response1p = nd.response1p;
  workingresponse1p = nd.workingresponse1p;
  return *this;
  }


double DISTR_bivprobit_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivprobit_rho::compute_param(const double * linpred,double * param)
  {
   double arg = (*linpred);
  *param = arg/pow(1+pow(arg,2),0.5);
  }

void DISTR_bivprobit_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

   response2p = workingresponse2p->getV();
   response1p = workingresponse1p->getV();
  }



void DISTR_bivprobit_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs-1)
    {
    response2p++;
    response1p++;
    }

  }



double DISTR_bivprobit_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu_2 equation
  // *worktransformlin[0] = mu_2;
  // *worklin[1] = linear predictor of mu_1 equation
  // *worktransformlin[1] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }
  double rho;

  if (*linpred <= -100)
    rho  = -0.99995;
  else if (*linpred >= 100)
    rho  = 0.99995;
  else
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);

  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;
  double l;


     l = -0.5*log(oneminusrho2) -(1/(2*oneminusrho2))*( pow((((*response1p))-(*worktransformlin[1])),2) -
                                 2*rho*(((*response1p)-(*worktransformlin[1])))*(((*response2p)-(*worktransformlin[0])))
                                +  pow((((*response2p))-(*worktransformlin[0])),2) );


  modify_worklin();

  return l;

  }

void DISTR_bivprobit_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of mu_2 equation
  // *worktransformlin[0] = mu_2;
  // *worklin[1] = linear predictor of mu_1 equation
  // *worktransformlin[1] = mu_1;

  if (counter==0)
    {
    set_worklin();
    }

  double rho;
  double hilfs;

  if (*linpred <= -100) {
    rho  = -0.99995;
    hilfs = 100.05;
  }
  else if (*linpred >= 100) {
    rho  = 0.99995;
    hilfs = 100.05;
  }
  else {
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);
    hilfs = pow((1+pow((*linpred),2)),0.5);
  }

//          std::ofstream out;
//  // helpmat1.prettyPrint(out);
//    out.open ("C:\\tmp\\bivprobit2.raw", std::ofstream::out | std::ofstream::app);
//    out << *workingresponse ;
//    out << " " ;
//    out << *response1p ;
//    out << " " ;
//    out << *response  ;
//    out << " " ;
//    out << *response2p << endl;
//    out.close();

  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;


    double nu = oneminusrho2*(*linpred) - (*linpred)*( pow((((*response1p))-(*worktransformlin[1])),2)
                                                      +  pow((((*response2p))-(*worktransformlin[0])),2) )
                +(hilfs+rho*(*linpred))*( (((*response1p)-(*worktransformlin[1])))*(((*response2p)-(*worktransformlin[0]))) );



    *workingweight = 1-pow(rho2,2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(oneminusrho2) -(1/(2*oneminusrho2))*( pow((((*response1p))-(*worktransformlin[1])),2) -
                                 2*rho*(((*response1p)-(*worktransformlin[1])))*(((*response2p)-(*worktransformlin[0])))
                                +  pow((((*response2p))-(*worktransformlin[0])),2) );

      }


  modify_worklin();



  }


void DISTR_bivprobit_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (rho):\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivprobit_rho::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = (*worklin)/pow((1+pow((*worklin),2)),0.5);
    }

  }

//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_bivprobit_mu ------------------------
//------------------------------------------------------------------------------
void DISTR_bivprobit_mu::check_errors(void)
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
            if((*workresp)!=0 | (*workresp)!=1) {
                errors=true;
                errormessages.push_back("ERROR: response has to be zero or one\n");
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


DISTR_bivprobit_mu::DISTR_bivprobit_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  pos =p;
  family = "bivariate probit - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  responseorig = response;
//    linpredminlimit=-10;
//  linpredmaxlimit=15;
  }


DISTR_bivprobit_mu::DISTR_bivprobit_mu(const DISTR_bivprobit_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  responseorig = nd.responseorig;
//  response2 = nd.response2;
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  }


const DISTR_bivprobit_mu & DISTR_bivprobit_mu::operator=(
                            const DISTR_bivprobit_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  responseorig = nd.responseorig;
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  return *this;
  }


void DISTR_bivprobit_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[2] = *response[4] = first component of two dimensional reponse
   // *linpred[0] = eta_rho
   // *linpred[1] = eta_mu_2
   // *linpred[2] = eta_mu_1

   if (*weight[2] == 0)
     *deviance=0;
   else
     {
     double rho = (*linpred[0])/pow(1+pow((*linpred[0]),2),0.5);
     double mu_2 = (*linpred[1]);
     double mu_1 = (*linpred[2]);
     double hilfs1 = 1-pow(rho,2);
     double l;

       l = -log(2*PI)-0.5*log(hilfs1)-(1/(2*hilfs1))*( pow((((*response[4]))-mu_1),2) -
                                                                                2*rho*(((*response[4])-mu_1))*(((*response[3])-mu_2))
                                                                                + pow((((*response[3]))-mu_2),2) );


    *deviance = -2*l;
    }

  }


double DISTR_bivprobit_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivprobit_mu::compute_param(const double * linpred,double * param)
  {
  *param = (*linpred);
  }

void DISTR_bivprobit_mu::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

    response2p = workingresponse2p->getV();

  }



void DISTR_bivprobit_mu::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs-1)
    response2p++;



  }



void DISTR_bivprobit_mu::update(void)
  {


  double * workresp = response.getV();
  double * workresporig =responseorig.getV();
  double * weightwork = weight.getV();

  double * worklin_current;
  if (linpred_current==1)
    worklin_current = linearpred1.getV();
  else
    worklin_current = linearpred2.getV();

  set_worklin();

  unsigned i;
  for(i=0;i<nrobs;i++,worklin_current++,workresp++,weightwork++,
           response2p++,workresporig++,worktransformlin[0]++,worklin[1]++)
    {

    if (*weightwork != 0)
      {
      if (*workresporig > 0)
        *workresp = trunc_normal2(0,20,*worklin_current+(*worktransformlin[0])*((*response2p)-(*worklin[1])),pow(1-pow(*worktransformlin[0],2),0.5));
      else
        *workresp = trunc_normal2(-20,0,*worklin_current+(*worktransformlin[0])*((*response2p)-(*worklin[1])),pow(1-pow(*worktransformlin[0],2),0.5));
      }

//          std::ofstream out;
//  // helpmat1.prettyPrint(out);
//    out.open ("C:\\tmp\\bivprobit.raw", std::ofstream::out | std::ofstream::app);
//    out << *workresp ;
//    out << " " ;
//    out << *workresporig ;
//    out << " " ;
//    out << *worklin[1] ;
//    out << " " ;
//    out << *worktransformlin[0] ;
//    out << " " ;
//    out << *worktransformlin[1] ;
//    out << " " ;
//    out << *response2p << endl;

    }




  }



double DISTR_bivprobit_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;


  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);
  double rho2 = pow((*worktransformlin[0]),2);
  double oneminusrho2 = 1- rho2;
  double l;


     l = -(1/(2*oneminusrho2))*( pow((((*response))-mu),2) -
                                 2*(*worktransformlin[0])*(((*response)-mu))*(((*response2p)-(*worktransformlin[1]))) );

  modify_worklin();

  return l;

  }


void DISTR_bivprobit_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;


  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double mu = (*linpred);
    double rho2 = pow((*worktransformlin[0]),2);
   double oneminusrho2 = 1- rho2;


    double nu = (1/(oneminusrho2))*( (((*response))-mu) -
                                 ((*worktransformlin[0]))*(((*response2p)-(*worktransformlin[1]))) );

    *workingweight = 1/(oneminusrho2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -(1/(2*oneminusrho2))*( pow((((*response))-mu),2) -
                                 2*(*worktransformlin[0])*(((*response)-mu))*(((*response2p)-(*worktransformlin[1]))) );

      }

  modify_worklin();

  }


void DISTR_bivprobit_mu::compute_mu_mult(vector<double *> linpred,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+1+pos]));
  }


void DISTR_bivprobit_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivprobit_mu::update_end(void)
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




} // end: namespace MCMC



