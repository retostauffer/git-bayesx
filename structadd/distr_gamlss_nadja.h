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



#if !defined (DISTRgamlssnadja_INCLUDED)
#define DISTRgamlssnadja_INCLUDED

#include"../export_type.h"
#include"statmat.h"
#include"Random.h"
#include"GENERAL_OPTIONS.h"
#include"FC.h"
#include"distr_gamlss.h"


namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_betainf1_tau -------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_betainf1_tau : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_betainf1_tau(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_betainf1_tau(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_betainf1_tau(const DISTR_betainf1_tau & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_betainf1_tau & operator=(const DISTR_betainf1_tau & nd);

   // DESTRUCTOR

  ~DISTR_betainf1_tau() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };




//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_betainf0_nu -------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_betainf0_nu : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_betainf0_nu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_betainf0_nu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_betainf0_nu(const DISTR_betainf0_nu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_betainf0_nu & operator=(const DISTR_betainf0_nu & nd);

   // DESTRUCTOR

  ~DISTR_betainf0_nu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };



//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_t_df ------------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_t_df : public DISTR_gamlss
  {

  protected:


  public:


   // DEFAULT CONSTRUCTOR

  DISTR_t_df(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_t_df(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_t_df(const DISTR_t_df & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_t_df & operator=(const DISTR_t_df & nd);

   // DESTRUCTOR

  ~DISTR_t_df() {}

  double get_intercept_start(void);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_t_sigma2 --------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_t_sigma2 : public DISTR_gamlss
  {

  protected:


  public:


   // DEFAULT CONSTRUCTOR

  DISTR_t_sigma2(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_t_sigma2(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_t_sigma2(const DISTR_t_sigma2 & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_t_sigma2 & operator=(const DISTR_t_sigma2 & nd);

   // DESTRUCTOR

  ~DISTR_t_sigma2() {}

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_t_mu ------------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_t_mu : public DISTR_gamlss
  {

  protected:


  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_t_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_t_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_t_mu(const DISTR_t_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_t_mu & operator=(const DISTR_t_mu & nd);

   // DESTRUCTOR

  ~DISTR_t_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };



//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_invgaussian_sigma2 ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_invgaussian_sigma2 : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_invgaussian_sigma2(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_invgaussian_sigma2(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_invgaussian_sigma2(const DISTR_invgaussian_sigma2 & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_invgaussian_sigma2 & operator=(const DISTR_invgaussian_sigma2 & nd);

   // DESTRUCTOR

  ~DISTR_invgaussian_sigma2() {}

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_invgaussian_mu --------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_invgaussian_mu : public DISTR_gamlss
  {

  protected:


  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_invgaussian_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_invgaussian_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_invgaussian_mu(const DISTR_invgaussian_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_invgaussian_mu & operator=(const DISTR_invgaussian_mu & nd);

   // DESTRUCTOR

  ~DISTR_invgaussian_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_pareto_p --------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_pareto_p : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_pareto_p(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_pareto_p(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_pareto_p(const DISTR_pareto_p & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_pareto_p & operator=(const DISTR_pareto_p & nd);

   // DESTRUCTOR

  ~DISTR_pareto_p() {}

  double get_intercept_start(void);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_pareto_b ------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_pareto_b : public DISTR_gamlss
  {

  protected:


  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_pareto_b(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_pareto_b(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_pareto_b(const DISTR_pareto_b & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_pareto_b & operator=(const DISTR_pareto_b & nd);

   // DESTRUCTOR

  ~DISTR_pareto_b() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_betainf_mu -------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_betainf_mu : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_betainf_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_betainf_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_betainf_mu(const DISTR_betainf_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_betainf_mu & operator=(const DISTR_betainf_mu & nd);

   // DESTRUCTOR

  ~DISTR_betainf_mu() {}

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };



//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_betainf_nu -------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_betainf_nu : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_betainf_nu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_betainf_nu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_betainf_nu(const DISTR_betainf_nu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_betainf_nu & operator=(const DISTR_betainf_nu & nd);

   // DESTRUCTOR

  ~DISTR_betainf_nu() {}

  double get_intercept_start(void);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };



//------------------------------------------------------------------------------
//------------------ CLASS: DISTR_betainf_sigma2 -------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_betainf_sigma2 : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_betainf_sigma2(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_betainf_sigma2(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_betainf_sigma2(const DISTR_betainf_sigma2 & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_betainf_sigma2 & operator=(const DISTR_betainf_sigma2 & nd);

   // DESTRUCTOR

  ~DISTR_betainf_sigma2() {}

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };



//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_betainf_tau -----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_betainf_tau : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_betainf_tau(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_betainf_tau(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_betainf_tau(const DISTR_betainf_tau & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_betainf_tau & operator=(const DISTR_betainf_tau & nd);

   // DESTRUCTOR

  ~DISTR_betainf_tau() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };





//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_dagum_p ---------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_dagum_p : public DISTR_gamlss
  {

  protected:


  public:


   // DEFAULT CONSTRUCTOR

  DISTR_dagum_p(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_dagum_p(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_dagum_p(const DISTR_dagum_p & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_dagum_p & operator=(const DISTR_dagum_p & nd);

   // DESTRUCTOR

  ~DISTR_dagum_p() {}

  double get_intercept_start(void);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_dagum_b ---------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_dagum_b : public DISTR_gamlss
  {

  protected:


  public:


   // DEFAULT CONSTRUCTOR

  DISTR_dagum_b(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_dagum_b(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_dagum_b(const DISTR_dagum_b & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_dagum_b & operator=(const DISTR_dagum_b & nd);

   // DESTRUCTOR

  ~DISTR_dagum_b() {}

  double get_intercept_start(void);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_dagum_a ---------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_dagum_a : public DISTR_gamlss
  {

  protected:


  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_dagum_a(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_dagum_a(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_dagum_a(const DISTR_dagum_a & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_dagum_a & operator=(const DISTR_dagum_a & nd);

   // DESTRUCTOR

  ~DISTR_dagum_a() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };




//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_weibull_alpha ---------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_weibull_alpha : public DISTR_gamlss
  {

  protected:


  public:


   // DEFAULT CONSTRUCTOR

  DISTR_weibull_alpha(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_weibull_alpha(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_weibull_alpha(const DISTR_weibull_alpha & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_weibull_alpha & operator=(const DISTR_weibull_alpha & nd);

   // DESTRUCTOR

  ~DISTR_weibull_alpha() {}

  double get_intercept_start(void);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_weibull_lambda --------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_weibull_lambda : public DISTR_gamlss
  {

  protected:


  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_weibull_lambda(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_weibull_lambda(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_weibull_lambda(const DISTR_weibull_lambda & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_weibull_lambda & operator=(const DISTR_weibull_lambda & nd);

   // DESTRUCTOR

  ~DISTR_weibull_lambda() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };




//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_gengamma_tau ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gengamma_tau : public DISTR_gamlss
  {

  protected:


  public:


   // DEFAULT CONSTRUCTOR

  DISTR_gengamma_tau(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_gengamma_tau(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_gengamma_tau(const DISTR_gengamma_tau & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_gengamma_tau & operator=(const DISTR_gengamma_tau & nd);

   // DESTRUCTOR

  ~DISTR_gengamma_tau() {}

  double get_intercept_start(void);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_gengamma_sigma --------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gengamma_sigma : public DISTR_gamlss
  {

  protected:


  public:


   // DEFAULT CONSTRUCTOR

  DISTR_gengamma_sigma(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_gengamma_sigma(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_gengamma_sigma(const DISTR_gengamma_sigma & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_gengamma_sigma & operator=(const DISTR_gengamma_sigma & nd);

   // DESTRUCTOR

  ~DISTR_gengamma_sigma() {}

  double get_intercept_start(void);


  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_gengamma_mu -----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gengamma_mu : public DISTR_gamlss
  {

  protected:


  public:

   void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_gengamma_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_gengamma_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_gengamma_mu(const DISTR_gengamma_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_gengamma_mu & operator=(const DISTR_gengamma_mu & nd);

   // DESTRUCTOR

  ~DISTR_gengamma_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };





//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_gamma_sigma -----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gamma_sigma : public DISTR_gamlss
  {

  protected:


  public:


   // DEFAULT CONSTRUCTOR

  DISTR_gamma_sigma(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_gamma_sigma(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_gamma_sigma(const DISTR_gamma_sigma & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_gamma_sigma & operator=(const DISTR_gamma_sigma & nd);

   // DESTRUCTOR

  ~DISTR_gamma_sigma() {}

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_gamma_mu --------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gamma_mu : public DISTR_gamlss
  {

  protected:


  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_gamma_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_gamma_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_gamma_mu(const DISTR_gamma_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_gamma_mu & operator=(const DISTR_gamma_mu & nd);

   // DESTRUCTOR

  ~DISTR_gamma_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_lognormal2_sigma ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_lognormal2_sigma : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_lognormal2_sigma(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_lognormal2_sigma(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_lognormal2_sigma(const DISTR_lognormal2_sigma & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_lognormal2_sigma & operator=(const DISTR_lognormal2_sigma & nd);

   // DESTRUCTOR

  ~DISTR_lognormal2_sigma() {}

  double get_intercept_start(void);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_lognormal2_mu ---------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_lognormal2_mu : public DISTR_gamlss
  {

  protected:


  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_lognormal2_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_lognormal2_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_lognormal2_mu(const DISTR_lognormal2_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_lognormal2_mu & operator=(const DISTR_lognormal2_mu & nd);

   // DESTRUCTOR

  ~DISTR_lognormal2_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  void update(void);

  };




//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_lognormal_sigma2 ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_lognormal_sigma2 : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_lognormal_sigma2(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_lognormal_sigma2(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_lognormal_sigma2(const DISTR_lognormal_sigma2 & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_lognormal_sigma2 & operator=(const DISTR_lognormal_sigma2 & nd);

   // DESTRUCTOR

  ~DISTR_lognormal_sigma2() {}

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_lognormal_mu ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_lognormal_mu : public DISTR_gamlss
  {

  protected:


  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_lognormal_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_lognormal_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_lognormal_mu(const DISTR_lognormal_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_lognormal_mu & operator=(const DISTR_lognormal_mu & nd);

   // DESTRUCTOR

  ~DISTR_lognormal_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  void update(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_truncnormal2_sigma ----------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_truncnormal2_sigma : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_truncnormal2_sigma(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_truncnormal2_sigma(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_truncnormal2_sigma(const DISTR_truncnormal2_sigma & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_truncnormal2_sigma & operator=(const DISTR_truncnormal2_sigma & nd);

   // DESTRUCTOR

  ~DISTR_truncnormal2_sigma() {}

  double get_intercept_start(void);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };

//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_truncnormal2_mu -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_truncnormal2_mu : public DISTR_gamlss
  {

  protected:


  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_truncnormal2_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_truncnormal2_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_truncnormal2_mu(const DISTR_truncnormal2_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_truncnormal2_mu & operator=(const DISTR_truncnormal2_mu & nd);

   // DESTRUCTOR

  ~DISTR_truncnormal2_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_normal2_sigma ---------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_normal2_sigma : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_normal2_sigma(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_normal2_sigma(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_normal2_sigma(const DISTR_normal2_sigma & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_normal2_sigma & operator=(const DISTR_normal2_sigma & nd);

   // DESTRUCTOR

  ~DISTR_normal2_sigma() {}

  double get_intercept_start(void);


  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };

  //-----------------------------------------------------------------------------
//----------------------- CLASS: DISTR_normal2_mu -------------------------------
//-------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_normal2_mu : public DISTR_gamlss
  {

  protected:


  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_normal2_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_normal2_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_normal2_mu(const DISTR_normal2_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_normal2_mu & operator=(const DISTR_normal2_mu & nd);

   // DESTRUCTOR

  ~DISTR_normal2_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  void update(void);

  };




//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_normal_sigma2 ---------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_normal_sigma2 : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_normal_sigma2(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_normal_sigma2(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_normal_sigma2(const DISTR_normal_sigma2 & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_normal_sigma2 & operator=(const DISTR_normal_sigma2 & nd);

   // DESTRUCTOR

  ~DISTR_normal_sigma2() {}

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_normal_mu -------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_normal_mu : public DISTR_gamlss
  {

  protected:


  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_normal_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_normal_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_normal_mu(const DISTR_normal_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_normal_mu & operator=(const DISTR_normal_mu & nd);

   // DESTRUCTOR

  ~DISTR_normal_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  void update(void);

  };





//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_beta_sigma2 -----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_beta_sigma2 : public DISTR_gamlss
  {

  protected:


  public:


   // DEFAULT CONSTRUCTOR

  DISTR_beta_sigma2(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_beta_sigma2(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_beta_sigma2(const DISTR_beta_sigma2 & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_beta_sigma2 & operator=(const DISTR_beta_sigma2 & nd);

   // DESTRUCTOR

  ~DISTR_beta_sigma2() {}

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };



//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_beta_mu ---------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_beta_mu : public DISTR_gamlss
  {

  protected:


  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_beta_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_beta_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_beta_mu(const DISTR_beta_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_beta_mu & operator=(const DISTR_beta_mu & nd);

   // DESTRUCTOR

  ~DISTR_beta_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_cloglog ---------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_cloglog : public DISTR_gamlss
  {

  protected:


  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_cloglog(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_cloglog(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_cloglog(const DISTR_cloglog & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_cloglog & operator=(const DISTR_cloglog & nd);

   // DESTRUCTOR

  ~DISTR_cloglog() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);


  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double compute_quantile_residual_mult(vector<double *> response,
                                         vector<double *> param,
                                         vector<double *> weight,
                                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };

//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_gumbelcopula_rho ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gumbelcopula_rho : public DISTR_gamlss
  {

  protected:

  void set_worklin(void);
  void modify_worklin(void);

  public:

  datamatrix response2;
  double * response2p;

   // DEFAULT CONSTRUCTOR

  DISTR_gumbelcopula_rho(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_gumbelcopula_rho(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_gumbelcopula_rho(const DISTR_gumbelcopula_rho & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_gumbelcopula_rho & operator=(const DISTR_gumbelcopula_rho & nd);

   // DESTRUCTOR

  ~DISTR_gumbelcopula_rho() {}

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_gumbelcopula ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gumbelcopula : public DISTR_gamlss
  {

  protected:

  unsigned pos;

  void set_worklin(void);
  void modify_worklin(void);


  public:

  datamatrix response2;
  double * response2p;

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_gumbelcopula(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_gumbelcopula(GENERAL_OPTIONS * o, const datamatrix & r, unsigned & p,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_gumbelcopula(const DISTR_gumbelcopula & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_gumbelcopula & operator=(const DISTR_gumbelcopula & nd);

   // DESTRUCTOR

  ~DISTR_gumbelcopula() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_gaussiancopula_rho ----------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gaussiancopula_rho : public DISTR_gamlss
  {

  protected:

  void set_worklin(void);
  void modify_worklin(void);

  public:

  datamatrix response2;
  double * response2p;

   // DEFAULT CONSTRUCTOR

  DISTR_gaussiancopula_rho(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_gaussiancopula_rho(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_gaussiancopula_rho(const DISTR_gaussiancopula_rho & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_gaussiancopula_rho & operator=(const DISTR_gaussiancopula_rho & nd);

   // DESTRUCTOR

  ~DISTR_gaussiancopula_rho() {}

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_gaussiancopula --------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gaussiancopula : public DISTR_gamlss
  {

  protected:

  unsigned pos;

  void set_worklin(void);
  void modify_worklin(void);


  public:

  datamatrix response2;
  double * response2p;

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_gaussiancopula(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_gaussiancopula(GENERAL_OPTIONS * o, const datamatrix & r, unsigned & p,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_gaussiancopula(const DISTR_gaussiancopula & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_gaussiancopula & operator=(const DISTR_gaussiancopula & nd);

   // DESTRUCTOR

  ~DISTR_gaussiancopula() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_dirichlet -------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_dirichlet : public DISTR_gamlss
  {

  protected:

  int nrcat;

  unsigned pos;

  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_dirichlet(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_dirichlet(GENERAL_OPTIONS * o, const datamatrix & r, int & nrc, unsigned & p,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_dirichlet(const DISTR_dirichlet & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_dirichlet & operator=(const DISTR_dirichlet & nd);

   // DESTRUCTOR

  ~DISTR_dirichlet() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

 // void compute_param(const double * linpred,double * param);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };




//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_bivt_df ---------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_bivt_df : public DISTR_gamlss
  {

  protected:

  void set_worklin(void);
  void modify_worklin(void);

  public:

  datamatrix response2;
  double * response2p;

   // DEFAULT CONSTRUCTOR

  DISTR_bivt_df(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_bivt_df(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_bivt_df(const DISTR_bivt_df & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_bivt_df & operator=(const DISTR_bivt_df & nd);

   // DESTRUCTOR

  ~DISTR_bivt_df() {}

  double get_intercept_start(void);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };

//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_bivt_rho --------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_bivt_rho : public DISTR_gamlss
  {

  protected:

  void set_worklin(void);
  void modify_worklin(void);

  public:

  datamatrix response2;
  double * response2p;

   // DEFAULT CONSTRUCTOR

  DISTR_bivt_rho(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_bivt_rho(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_bivt_rho(const DISTR_bivt_rho & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_bivt_rho & operator=(const DISTR_bivt_rho & nd);

   // DESTRUCTOR

  ~DISTR_bivt_rho() {}

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_bivt_sigma ------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_bivt_sigma : public DISTR_gamlss
  {

  protected:

  unsigned pos;

  void set_worklin(void);
  void modify_worklin(void);

  public:

  datamatrix response2;
  double * response2p;

   // DEFAULT CONSTRUCTOR

  DISTR_bivt_sigma(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_bivt_sigma(GENERAL_OPTIONS * o, const datamatrix & r,unsigned & p,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_bivt_sigma(const DISTR_bivt_sigma & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_bivt_sigma & operator=(const DISTR_bivt_sigma & nd);

   // DESTRUCTOR

  ~DISTR_bivt_sigma() {}

  double get_intercept_start(void);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_bivt_mu ---------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_bivt_mu : public DISTR_gamlss
  {

  protected:

  unsigned pos;

  void set_worklin(void);
  void modify_worklin(void);

  public:

  datamatrix response2;
  double * response2p;

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_bivt_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_bivt_mu(GENERAL_OPTIONS * o, const datamatrix & r,unsigned & p,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_bivt_mu(const DISTR_bivt_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_bivt_mu & operator=(const DISTR_bivt_mu & nd);

   // DESTRUCTOR

  ~DISTR_bivt_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);


  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_bivnormal_rho ---------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_bivnormal_rho : public DISTR_gamlss
  {

  protected:

  void set_worklin(void);
  void modify_worklin(void);

  public:

  datamatrix response2;
  double * response2p;

   // DEFAULT CONSTRUCTOR

  DISTR_bivnormal_rho(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_bivnormal_rho(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_bivnormal_rho(const DISTR_bivnormal_rho & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_bivnormal_rho & operator=(const DISTR_bivnormal_rho & nd);

   // DESTRUCTOR

  ~DISTR_bivnormal_rho() {}

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_bivnormal_sigma ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_bivnormal_sigma : public DISTR_gamlss
  {

  protected:

  unsigned pos;

  void set_worklin(void);
  void modify_worklin(void);

  public:

  datamatrix response2;
  double * response2p;


   // DEFAULT CONSTRUCTOR

  DISTR_bivnormal_sigma(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_bivnormal_sigma(GENERAL_OPTIONS * o, const datamatrix & r,unsigned & p,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_bivnormal_sigma(const DISTR_bivnormal_sigma & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_bivnormal_sigma & operator=(const DISTR_bivnormal_sigma & nd);

   // DESTRUCTOR

  ~DISTR_bivnormal_sigma() {}

  double get_intercept_start(void);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_bivnormal_mu -------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_bivnormal_mu : public DISTR_gamlss
  {

  protected:

  unsigned pos;

  void set_worklin(void);
  void modify_worklin(void);


  public:

  datamatrix response2;
  double * response2p;

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_bivnormal_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_bivnormal_mu(GENERAL_OPTIONS * o, const datamatrix & r, unsigned & p,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_bivnormal_mu(const DISTR_bivnormal_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_bivnormal_mu & operator=(const DISTR_bivnormal_mu & nd);

   // DESTRUCTOR

  ~DISTR_bivnormal_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  void update(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_bivprobit_rho ---------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_bivprobit_rho : public DISTR_gamlss
  {

  protected:

  void set_worklin(void);
  void modify_worklin(void);

  public:
  datamatrix responseorig;
  //datamatrix response2;
  double * response2p;

  double * response1p;

  datamatrix * workingresponse2p;

  datamatrix * workingresponse1p;

   // DEFAULT CONSTRUCTOR

  DISTR_bivprobit_rho(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_bivprobit_rho(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_bivprobit_rho(const DISTR_bivprobit_rho & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_bivprobit_rho & operator=(const DISTR_bivprobit_rho & nd);

   // DESTRUCTOR

  ~DISTR_bivprobit_rho() {}

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_bivprobit_mu ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_bivprobit_mu : public DISTR_gamlss
  {

  protected:

  unsigned pos;

  void set_worklin(void);
  void modify_worklin(void);


  public:

  datamatrix responseorig;
  //datamatrix response2;
  double * response2p;

  datamatrix * workingresponse2p;

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_bivprobit_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_bivprobit_mu(GENERAL_OPTIONS * o, const datamatrix & r, unsigned & p,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_bivprobit_mu(const DISTR_bivprobit_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_bivprobit_mu & operator=(const DISTR_bivprobit_mu & nd);

   // DESTRUCTOR

  ~DISTR_bivprobit_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update(void);

  void update_end(void);


  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_bivlogit_or -----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_bivlogit_or : public DISTR_gamlss
  {

  protected:

  void set_worklin(void);
  void modify_worklin(void);

  public:

  datamatrix response2;
  double * response2p;

   // DEFAULT CONSTRUCTOR

  DISTR_bivlogit_or(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_bivlogit_or(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_bivlogit_or(const DISTR_bivlogit_or & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_bivlogit_or & operator=(const DISTR_bivlogit_or & nd);

   // DESTRUCTOR

  ~DISTR_bivlogit_or() {}

  double get_intercept_start(void);


  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };



//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_bivlogit_mu ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_bivlogit_mu : public DISTR_gamlss
  {

  protected:

  unsigned pos;

  void set_worklin(void);
  void modify_worklin(void);


  public:

  datamatrix response2;
  double * response2p;

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_bivlogit_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_bivlogit_mu(GENERAL_OPTIONS * o, const datamatrix & r, unsigned & p,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_bivlogit_mu(const DISTR_bivlogit_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_bivlogit_mu & operator=(const DISTR_bivlogit_mu & nd);

   // DESTRUCTOR

  ~DISTR_bivlogit_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_BCCG_nu ---------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_BCCG_nu : public DISTR_gamlss
  {

  protected:


  public:


   // DEFAULT CONSTRUCTOR

  DISTR_BCCG_nu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_BCCG_nu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_BCCG_nu(const DISTR_BCCG_nu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_BCCG_nu & operator=(const DISTR_BCCG_nu & nd);

   // DESTRUCTOR

  ~DISTR_BCCG_nu() {}

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_BCCG_sigma ------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_BCCG_sigma : public DISTR_gamlss
  {

  protected:


  public:


   // DEFAULT CONSTRUCTOR

  DISTR_BCCG_sigma(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_BCCG_sigma(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_BCCG_sigma(const DISTR_BCCG_sigma & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_BCCG_sigma & operator=(const DISTR_BCCG_sigma & nd);

   // DESTRUCTOR

  ~DISTR_BCCG_sigma() {}

  double get_intercept_start(void);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_BCCG_mu ---------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_BCCG_mu : public DISTR_gamlss
  {

  protected:


  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_BCCG_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_BCCG_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_BCCG_mu(const DISTR_BCCG_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_BCCG_mu & operator=(const DISTR_BCCG_mu & nd);

   // DESTRUCTOR

  ~DISTR_BCCG_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);



  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_sfa_alpha -------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_sfa_alpha : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_sfa_alpha(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_sfa_alpha(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_sfa_alpha(const DISTR_sfa_alpha & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_sfa_alpha & operator=(const DISTR_sfa_alpha & nd);

   // DESTRUCTOR

  ~DISTR_sfa_alpha() {}

  double get_intercept_start(void);

//  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };

//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_sfa_sigma_v -----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_sfa_sigma_v : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_sfa_sigma_v(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_sfa_sigma_v(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_sfa_sigma_v(const DISTR_sfa_sigma_v & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_sfa_sigma_v & operator=(const DISTR_sfa_sigma_v & nd);

   // DESTRUCTOR

  ~DISTR_sfa_sigma_v() {}

  double get_intercept_start(void);

//  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };

//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_sfa_sigma_u -----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_sfa_sigma_u : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_sfa_sigma_u(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_sfa_sigma_u(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_sfa_sigma_u(const DISTR_sfa_sigma_u & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_sfa_sigma_u & operator=(const DISTR_sfa_sigma_u & nd);

   // DESTRUCTOR

  ~DISTR_sfa_sigma_u() {}

  double get_intercept_start(void);

//  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };

//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_sfa_mu_u --------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_sfa_mu_u : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_sfa_mu_u(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_sfa_mu_u(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_sfa_mu_u(const DISTR_sfa_mu_u & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_sfa_mu_u & operator=(const DISTR_sfa_mu_u & nd);

   // DESTRUCTOR

  ~DISTR_sfa_mu_u() {}

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };

//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_sfa_mu_y --------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_sfa_mu_y : public DISTR_gamlss
  {

  protected:


  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_sfa_mu_y(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_sfa_mu_y(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_sfa_mu_y(const DISTR_sfa_mu_y & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_sfa_mu_y & operator=(const DISTR_sfa_mu_y & nd);

   // DESTRUCTOR

  ~DISTR_sfa_mu_y() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_sfa2_sigma_v -----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_sfa2_sigma_v : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_sfa2_sigma_v(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_sfa2_sigma_v(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_sfa2_sigma_v(const DISTR_sfa2_sigma_v & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_sfa2_sigma_v & operator=(const DISTR_sfa2_sigma_v & nd);

   // DESTRUCTOR

  ~DISTR_sfa2_sigma_v() {}

  double get_intercept_start(void);

//  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };

//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_sfa2_sigma_u -----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_sfa2_sigma_u : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_sfa2_sigma_u(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_sfa2_sigma_u(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_sfa2_sigma_u(const DISTR_sfa2_sigma_u & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_sfa2_sigma_u & operator=(const DISTR_sfa2_sigma_u & nd);

   // DESTRUCTOR

  ~DISTR_sfa2_sigma_u() {}

  double get_intercept_start(void);

//  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };

//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_sfa2_mu_u --------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_sfa2_mu_u : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_sfa2_mu_u(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_sfa2_mu_u(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_sfa2_mu_u(const DISTR_sfa2_mu_u & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_sfa2_mu_u & operator=(const DISTR_sfa2_mu_u & nd);

   // DESTRUCTOR

  ~DISTR_sfa2_mu_u() {}

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void update_end(void);

  };

//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_sfa2_mu_y --------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_sfa2_mu_y : public DISTR_gamlss
  {

  protected:


  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_sfa2_mu_y(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_sfa2_mu_y(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_sfa2_mu_y(const DISTR_sfa2_mu_y & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_sfa2_mu_y & operator=(const DISTR_sfa2_mu_y & nd);

   // DESTRUCTOR

  ~DISTR_sfa2_mu_y() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  void compute_param(const double * linpred,double * param);

  double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);


  void outoptions(void);

  void update_end(void);

  };


} // end: namespace MCMC


#endif
