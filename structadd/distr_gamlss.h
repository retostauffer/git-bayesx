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
//----------------------- CLASS: DISTR_negbinzip_mu ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_negbinzip_mu : public DISTR
  {

  protected:

  unsigned counter;

  double * worklinpi;
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
                             double * deviancesat,
                             vector<double> scale) const;

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

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

//  void sample_responses(unsigned i,datamatrix & sr);

//  void sample_responses_cv(unsigned i,datamatrix & linpred,datamatrix & sr);

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

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);


  void outoptions(void);

  void posteriormode_end(void);

  void update_end(void);

//  void sample_responses(unsigned i,datamatrix & sr);

//  void sample_responses_cv(unsigned i,datamatrix & linpred,datamatrix & sr);

  };


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_negbinzip_delta -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_negbinzip_delta : public DISTR
  {

  protected:

  double stopsum;
  int stoprmax;

  // auxiliary variables for iwls

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

  double explinpi;
  double log_one_explinpi;
  double log_explinpi_pot;

  double lng_delta;
  double delta_linpred;
  double log_delta_mu;

  double E_dig_y_delta;
  double E_trig_y_delta;


  // end auxiliary variables for iwls

  unsigned counter;

  double responsemax;

  double * worklinmu;
  double * workexplinmu;

  double * worklinpi;
  double * workonempi;

  void set_worklinmupi(void);
  void modify_worklinmupi(void);

  void compute_expectation(void);

  public:

  DISTR*  distrmu;
  DISTR*  distrpi;

   // DEFAULT CONSTRUCTOR

   DISTR_negbinzip_delta(void) : DISTR()
     {
     }

   // CONSTRUCTOR

   DISTR_negbinzip_delta(GENERAL_OPTIONS * o, const datamatrix & r,
                         double & stpsum, int & strmax,
                         const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_negbinzip_delta(const DISTR_negbinzip_delta & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_negbinzip_delta & operator=(const DISTR_negbinzip_delta & nd);

   // DESTRUCTOR

   ~DISTR_negbinzip_delta() {}

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void outoptions(void);

  void posteriormode_end(void);

  void update_end(void);

//  void sample_responses(unsigned i,datamatrix & sr);

//  void sample_responses_cv(unsigned i,datamatrix & linpred,datamatrix & sr);

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
                             double * deviancesat,
                             vector<double> scale) const;

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);

  void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,double * mu);

  void outoptions(void);


  void posteriormode_end(void);

  void update_end(void);


//  void sample_responses(unsigned i,datamatrix & sr);

//  void sample_responses_cv(unsigned i,datamatrix & linpred,datamatrix & sr);

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

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             double * deviancesat,
                             vector<double> scale) const;


  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);

  void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like);

  void compute_mu_mult(vector<double *> linpred,double * mu);

  void outoptions(void);

  void posteriormode_end(void);

  void update_end(void);



//  void sample_responses(unsigned i,datamatrix & sr);

//  void sample_responses_cv(unsigned i,datamatrix & linpred,
//                                   datamatrix & sr);

  };



} // end: namespace MCMC


#endif
