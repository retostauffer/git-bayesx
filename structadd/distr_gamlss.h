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

  void compute_mu(const double * linpred,double * mu);

  void compute_deviance(const double * response, const double * weight,
                        const double * mu,double * deviance,
                        double * deviancesat, double * scale) const;

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

  void compute_mu(const double * linpred,double * mu);

  void compute_deviance(const double * response, const double * weight,
                        const double * mu,double * deviance,
                        double * deviancesat, double * scale) const;

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

  void outoptions(void);

  void posteriormode_end(void);

  void update_end(void);



//  void sample_responses(unsigned i,datamatrix & sr);

//  void sample_responses_cv(unsigned i,datamatrix & linpred,
//                                   datamatrix & sr);

  };



} // end: namespace MCMC


#endif