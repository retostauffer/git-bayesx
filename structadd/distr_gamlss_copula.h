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



#if !defined (DISTR_gamlss_copula_INCLUDED)
#define DISTR_gamlss_copula_INCLUDED

#include"../export_type.h"
#include"statmat.h"
#include"Random.h"
#include"GENERAL_OPTIONS.h"
#include"FC.h"
#include"distr_gamlss.h"
#include"distr_gamlss_nadja.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//-------------------------- CLASS: DISTR_gausscopula --------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gausscopula : public DISTR_gamlss
  {

  protected:

  void set_worklin(void);
  void modify_worklin(void);

  public:

  datamatrix response2;
  double * responsep;
  double * response2p;

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_gausscopula(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_gausscopula(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_gausscopula(const DISTR_gausscopula & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_gausscopula & operator=(const DISTR_gausscopula & nd);

   // DESTRUCTOR

  ~DISTR_gausscopula() {}

  double get_intercept_start(void);

  void compute_param_mult(vector<double *>  linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  vector<double> derivative(double * linpred);

  double logc(double & F, int & copulapos, double * linpred);

  double logc(double & F1, double & F2, double & eta);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);
  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };

} // end: namespace MCMC


#endif
