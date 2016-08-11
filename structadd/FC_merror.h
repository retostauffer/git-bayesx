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



#if !defined (FCmerrorINCLUDED)

#define FCmerrorINCLUDED

#include"../export_type.h"
#include"../values.h"
#include<fstream>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include"FC_nonp.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_merror ---------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_merror  : public FC
  {

  protected:

  datamatrix xobs;
  datamatrix xvar;
  datamatrix xmean;

  double minbin;
  double maxbin;
  double deltabin;
  double binning;

  datamatrix mevar;

  FC_nonp * FCp;

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_merror(void);

  // CONSTRUCTOR

  FC_merror(GENERAL_OPTIONS * o, const ST::string & t,
            const ST::string & fp, vector<ST::string> & op,
            vector<ST::string> & vn, datamatrix & xo, datamatrix & mv,
            datamatrix & xd,
            FC_nonp * fcn);

  // COPY CONSTRUCTOR

  FC_merror(const FC_merror & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_merror & operator=(const FC_merror & m);

  // DESTRUCTOR

  ~FC_merror()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

  void outresults(ofstream & out_stata,ofstream & out_R, ofstream & out_R2BayesX,
                 const ST::string & pathresults);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  };

} // end: namespace MCMC

#endif


