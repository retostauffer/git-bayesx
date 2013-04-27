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



#if !defined (DISTRgamlss_INCLUDED)
#define DISTRgamlss_INCLUDED

#include"../export_type.h"
#include"statmat.h"
#include"Random.h"
#include"GENERAL_OPTIONS.h"
#include"FC.h"
#include"distr.h"


namespace MCMC
{


//------------------------------------------------------------------------------
//-------------------------- CLASS: DISTR_gamlss -------------------------------
//------------------------------------------------------------------------------

/*
const double linpredlimit = -10.0;
const double linpredlimitmax = 15.0;
const double explinpredlimit = exp(linpredlimit);
const double explinpredlimitmax = exp(linpredlimitmax);
const double expminusexplinpredlimit = 0.9999546;     // for exp(-exp(eta))
*/

class __EXPORT_TYPE DISTR_gamlss : public DISTR
  {

  protected:

  unsigned counter;

  vector<double *> worklin;
  vector<double *> worktransformlin;

  virtual void set_worklin(void);
  virtual void modify_worklin(void);


  public:

  vector<DISTR*> distrp;

   // DEFAULT CONSTRUCTOR

  DISTR_gamlss(void) : DISTR()
    {
    }

   // CONSTRUCTOR

  DISTR_gamlss(GENERAL_OPTIONS * o, const datamatrix & r,unsigned nrdistr,
                   const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_gamlss(const DISTR_gamlss & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_gamlss & operator=(const DISTR_gamlss & nd);

   // DESTRUCTOR

  ~DISTR_gamlss() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,double * mu);

  void outoptions(void);

  void posteriormode_end(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_negbin_delta ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_negbin_delta : public DISTR_gamlss
  {

  protected:

  double E_dig_y_delta;
  double E_trig_y_delta;
  double delta;
  double log_delta_div_delta_plus_mu;
  double lngamma_delta;
  double delta_plus_mu;

  datamatrix E_dig_y_delta_m;
  datamatrix E_trig_y_delta_m;  
  double * Ep;
  double * Ep_trig;  


  double stopsum;
  int stoprmax;
  int nrbetween;

  bool slow;

  void compute_expectation(void);

  public:

   // DEFAULT CONSTRUCTOR

  DISTR_negbin_delta(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_negbin_delta(GENERAL_OPTIONS * o, const datamatrix & r,
                     double & ss, int & strmax, int & sts, bool & sl,
                     const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_negbin_delta(const DISTR_negbin_delta & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_negbin_delta & operator=(const DISTR_negbin_delta & nd);

   // DESTRUCTOR

  ~DISTR_negbin_delta() {}

  double get_intercept_start(void);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

/*
  void compute_iwls_wweightschange_weightsone_slow(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);
*/

  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_negbin_mu -------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_negbin_mu : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_negbin_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_negbin_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_negbin_mu(const DISTR_negbin_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_negbin_mu & operator=(const DISTR_negbin_mu & nd);

   // DESTRUCTOR

  ~DISTR_negbin_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,double * mu);

  void outoptions(void);

  void update_end(void);

  };



//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_zip_cloglog_mu --------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_zip_cloglog_mu : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_zip_cloglog_mu(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_zip_cloglog_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_zip_cloglog_mu(const DISTR_zip_cloglog_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_zip_cloglog_mu & operator=(const DISTR_zip_cloglog_mu & nd);

   // DESTRUCTOR

  ~DISTR_zip_cloglog_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,double * mu);

  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_zip_cloglog_pi --------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_zip_cloglog_pi : public DISTR_gamlss
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

  DISTR_zip_cloglog_pi(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_zip_cloglog_pi(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_zip_cloglog_pi(const DISTR_zip_cloglog_pi & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_zip_cloglog_pi & operator=(const DISTR_zip_cloglog_pi & nd);

   // DESTRUCTOR

  ~DISTR_zip_cloglog_pi() {}

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
//----------------------- CLASS: DISTR_negbinzip_mu ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_negbinzip_mu : public DISTR
  {

  protected:

  unsigned counter;

  double * worklinpi;
  double * workexplinpi;
  double * workonempi;

  double * worklindelta;
  double * workexplindelta;

  void set_worklinpidelta(void);
  void modify_worklinpidelta(void);


  public:

  DISTR*  distrpi;
  DISTR*  distrdelta;

   // DEFAULT CONSTRUCTOR

   DISTR_negbinzip_mu(void) : DISTR()
     {
     }

   // CONSTRUCTOR

   DISTR_negbinzip_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                   const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_negbinzip_mu(const DISTR_negbinzip_mu & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_negbinzip_mu & operator=(const DISTR_negbinzip_mu & nd);

   // DESTRUCTOR

   ~DISTR_negbinzip_mu() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix *> aux);

  double get_intercept_start(void);

  double loglikelihood(double * response, double * linpred,
                                         double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);


  void compute_mu_mult(vector<double *> linpred,double * mu);

  void outoptions(void);

  void posteriormode_end(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_negbinzip_pi ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_negbinzip_pi : public DISTR
  {

  protected:

  unsigned counter;

  double * worklinmu;
  double * workexplinmu;

  double * worklindelta;
  double * workexplindelta;

  void set_worklinmudelta(void);
  void modify_worklinmudelta(void);


  public:

  DISTR*  distrmu;
  DISTR*  distrdelta;

   // DEFAULT CONSTRUCTOR

   DISTR_negbinzip_pi(void) : DISTR()
     {
     }

   // CONSTRUCTOR

   DISTR_negbinzip_pi(GENERAL_OPTIONS * o, const datamatrix & r,
                   const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_negbinzip_pi(const DISTR_negbinzip_pi & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_negbinzip_pi & operator=(const DISTR_negbinzip_pi & nd);

   // DESTRUCTOR

   ~DISTR_negbinzip_pi() {}

  double get_intercept_start(void);

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);



  void outoptions(void);

  void posteriormode_end(void);

  void update_end(void);


  };



//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_negbinzip_delta -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_negbinzip_delta : public DISTR
  {

  protected:

  double stopsum;
  double fraclimit;
  int stoprmax;
  bool slow;
  int nrbetween;

  //----------------------------------------------------------------------------
/*
  double delta;
  double delta2;
  double deltay;
  double dig_deltay;
  double dig_delta;
  double trig_deltay;
  double trig_delta;
  double deltamu;
  double delta_div_deltamu;
  double log_delta_div_deltamu;
  double mu_m_y_div_delta_m_mu;
  double pi;
  double mu_div_deltamu;
  double pot;
  double denom;
  double sum;
  double log_one_explinpi;
  double log_explinpi_pot;
  double lng_delta;
  double delta_linpred;
  double log_delta_mu;
  double E_dig_y_delta;
  double E_trig_y_delta;
*/
  datamatrix E_dig_y_delta_m;
  datamatrix E_trig_y_delta_m;  
  double * Ep;
  double * Ep_trig;

  //----------------------------------------------------------------------------

  unsigned counter;

  double responsemax;

  double * worklinmu;
  double * workexplinmu;

  double * worklinpi;
  double * workonempi;
  double * workexplinpi;

  void set_worklinmupi(void);
  void modify_worklinmupi(void);



  public:

  DISTR*  distrmu;
  DISTR*  distrpi;

   // DEFAULT CONSTRUCTOR

   DISTR_negbinzip_delta(void) : DISTR()
     {
     }

   // CONSTRUCTOR

   DISTR_negbinzip_delta(GENERAL_OPTIONS * o, const datamatrix & r,
                                          double & ss, int & strmax,
                                          int & sts,bool & sl,
                                          const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_negbinzip_delta(const DISTR_negbinzip_delta & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_negbinzip_delta & operator=(const DISTR_negbinzip_delta & nd);

   // DESTRUCTOR

   ~DISTR_negbinzip_delta() {}

  double get_intercept_start(void);

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

/*
  void compute_iwls_wweightschange_weightsone_slow(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);
*/

  double compute_iwls(double * response, double * linpred,
                                     double * weight, double * workingweight,
                                     double * workingresponse,
                                     const bool & like);


  void outoptions(void);

  void posteriormode_end(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//---------------------- CLASS: DISTRIBUTION_ziplambda -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_ziplambda : public DISTR
  {

  protected:

  double * worklinpi;
  double * workonempi;
  double * workexplinpi;

  unsigned counter;

  void set_worklinpi(void);
  void modify_worklinpi(void);

  public:

  DISTR*  distrpi;

   // DEFAULT CONSTRUCTOR

   DISTR_ziplambda(void) : DISTR()
     {
     }

   // CONSTRUCTOR

   DISTR_ziplambda(GENERAL_OPTIONS * o, const datamatrix & r,
                  const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_ziplambda(const DISTR_ziplambda & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_ziplambda & operator=(const DISTR_ziplambda & nd);

   // DESTRUCTOR

   ~DISTR_ziplambda() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  double get_intercept_start(void);

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like);

  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);


  void compute_mu_mult(vector<double *> linpred,double * mu);

  void outoptions(void);


  void posteriormode_end(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//------------------------ CLASS: DISTRIBUTION_zipi ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_zippi : public DISTR
  {

  protected:

  double * worklinlambda;
  double * worklambda;
  double * workexpmlambda;

  unsigned counter;

  void set_worklinlambda(void);
  void modify_worklinlambda(void);


  public:

  DISTR*  distrlambda;

   // DEFAULT CONSTRUCTOR

   DISTR_zippi(void) : DISTR()
     {
     }

   // CONSTRUCTOR

   DISTR_zippi(GENERAL_OPTIONS * o, const datamatrix & r,
                  const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_zippi(const DISTR_zippi & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_zippi & operator=(const DISTR_zippi & nd);

   // DESTRUCTOR

   ~DISTR_zippi() {}

  double get_intercept_start(void);

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);


  void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like);


  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);


  void outoptions(void);

  void posteriormode_end(void);

  void update_end(void);

  };



} // end: namespace MCMC


#endif
