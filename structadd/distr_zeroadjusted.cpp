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

#include "distr_zeroadjusted.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS DISTR_zeroadjusted -------------------------
//------------------------------------------------------------------------------


DISTR_zeroadjusted::DISTR_zeroadjusted(GENERAL_OPTIONS * o,DISTR* dpi,
                                       DISTR* dmu)


  {
  optionsp = o;
  distrp_pi = dpi;
  distrp_mu = dmu;

  response = distrp_pi->response;
  workingresponse = response;

  maindistribution=true;
  predict_mult=false;

  nrobs = response.rows();

/*
  if (w.rows() == 1)
    {
    weight = datamatrix(r.rows(),1,1);
    }
  else
    {
    weight = w;
    }

  workingweight = weight;

  weightsone = check_weightsone();

  nrzeroweights = compute_nrzeroweights();

  wtype = wweightschange_weightsneqone;

  weightname = "W";
*/
  linearpred1 = datamatrix(nrobs,1,0);
  linearpred2 = datamatrix(nrobs,1,0);

  linpred_current = 1;

  }


const DISTR_zeroadjusted & DISTR_zeroadjusted::operator=(
const DISTR_zeroadjusted & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  distrp_pi = nd.distrp_pi;
  distrp_mu = nd.distrp_mu;
  return *this;
  }


DISTR_zeroadjusted::DISTR_zeroadjusted(const DISTR_zeroadjusted & nd)
   : DISTR(DISTR(nd))
  {
  distrp_pi = nd.distrp_pi;
  distrp_mu = nd.distrp_mu;
  }


void DISTR_zeroadjusted::outoptions(void)
  {

  }


void DISTR_zeroadjusted::compute_mu_mult(vector<double *> linpred,double * mu)
  {
  }


void DISTR_zeroadjusted::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {
  //  *deviance = -2*l;
  }





} // end: namespace MCMC



