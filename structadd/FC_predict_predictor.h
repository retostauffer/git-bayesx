/* BayesX - Software for Bayesian Inference in
Structured Additive Regression Models.
Copyright (C) 2019 Christiane Belitz, Andreas Brezger,
Nadja Klein, Thomas Kneib, Stefan Lang, Nikolaus Umlauf

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



#if !defined (FCpredictpredictorINCLUDED)

#define FCpredictpredictorINCLUDED

#include"../export_type.h"
#include"statmat.h"
#include"sparsemat.h"

#include"Random.h"
#include"../values.h"
#include<fstream>
#include<vector>
#include<bitset>
#include"GENERAL_OPTIONS.h"
#include"distr.h"
#include"clstring.h"
#include<cmath>

namespace MCMC
{

using std::vector;
using std::bitset;

//------------------------------------------------------------------------------
//---------------------- CLASS: FC_predict_predictor ---------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_predict_predictor   : public FC
  {

  protected:

  DISTR * likep;
  datamatrix designmatrix;
  vector<ST::string> varnames;

  void get_predictor(void);

  public:

  // DEFAULT CONSTRUCTOR

  FC_predict_predictor(void);

  // CONSTRUCTOR

  FC_predict_predictor(GENERAL_OPTIONS * o,DISTR * lp,const ST::string & t,
     const ST::string & fp,const ST::string & fpd, datamatrix & dm,
     vector<ST::string> & dn);

  // COPY CONSTRUCTOR

  FC_predict_predictor(const FC_predict_predictor & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_predict_predictor & operator=(const FC_predict_predictor & m);

  // DESTRUCTOR

  ~FC_predict_predictor()
    {
    }


  void update(void);

  bool posteriormode(void);

  void outresults(ofstream & out_stata, ofstream & out_R, ofstream & out_R2BayesX,
                  const ST::string & pathresults);

  void compute_autocorr_all(const ST::string & path, unsigned lag,
                                    ofstream & outg) const;

  void reset(void);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  };


} // end: namespace MCMC

#endif


