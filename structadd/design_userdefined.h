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



#if !defined (DESIGNuserdefinedINCLUDED)

#define DESIGNuserdefinedINCLUDED

#include"../export_type.h"
#include<deque>
#include<design.h>

using std::deque;

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: DESIGN_userdefined ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN_userdefined : public DESIGN
  {

  protected:

  vector< vector<double> > Zout2;            // Nonzero Elements of Z
  vector< vector<int> > index_Zout2;         // Columns of nonzero elements of Z

  datamatrix mK;                             // matrix to adjust for nonzero prior mean

  public:

  // FUNCTION: compute_f
  // TASK: compute Zout*beta, i.e. the estimated/current function evaluated at
  //       the different observations in data

  void compute_f(datamatrix & beta,datamatrix & betalin,
                       datamatrix & f, datamatrix & ftot);

  void compute_Zout(datamatrix & Z);
  void compute_Zout_transposed(datamatrix & Z);

  void init_data(datamatrix & dm, datamatrix & iv);

  void compute_basisNull(void);

  double round;
  double binning;

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN_userdefined(void);

  // CONSTRUCTOR

  DESIGN_userdefined(datamatrix & dm, datamatrix & iv, datamatrix & designmat,
             datamatrix & penmat, datamatrix & priormean,
             GENERAL_OPTIONS * o, DISTR * dp, FC_linear * fcl,
             vector<ST::string> & op,
             vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  DESIGN_userdefined(const DESIGN_userdefined & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN_userdefined & operator=(const DESIGN_userdefined & m);

  void compute_precision(double l);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  void outoptions(GENERAL_OPTIONS * op);

  void outbasis_R(ofstream & out);

// DESTRUCTOR

  ~DESIGN_userdefined() {}

  // FUNCTION: computes XWres
  // TASK: computes XWres, res is the partial residual
  //       l is the inverse smoothing variance (1/tau2)

  void compute_XtransposedWres(datamatrix & partres, double l, double t2);

  void compute_Zout_transposed(void);
  };


} // end: namespace MCMC

#endif


