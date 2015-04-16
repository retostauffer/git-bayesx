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



#include "design_userdefined.h"

using std::deque;

namespace MCMC
{


//------------------------------------------------------------------------------
//--------- CLASS: DESIGN_userdefined implementation of member functions -------
//------------------------------------------------------------------------------

void DESIGN_userdefined::compute_f(datamatrix & beta,datamatrix & betalin,
                       datamatrix & f, datamatrix & ftot)
  {

  }

void DESIGN_userdefined::read_options(vector<ST::string> & op,
                                 vector<ST::string> & vn)
  {

  /*
  1       degree
  2       numberknots
  3       difforder
  4       lambda
  5       a
  6       b
  7       center
  8       map
  9       lambda_re
  10      a_re
  11      b_re
  12      internal_mult
  13      samplemult
  14      constraints
  15      round
  16      centermethod
  17      internal_mult
  18      pvalue
  19      meaneffect
  20      binning
  21      update
  22      nu
  23      maxdist
  24      ccovariate
  */


  int f;

  if (op[7] == "false")   //nocenter==false, i.e. center
    center = true;
  else
    center = false;

  f = op[15].strtodouble(round);

  if (op[16]=="meancoeff")
    centermethod = meancoeff;
  else if (op[16] == "meansimple")
    centermethod = meansimple;
  else if (op[16] == "integralsimple")
    centermethod = integralsimple;
  else if (op[16] == "nullspace")
    centermethod = nullspace;
  else if (op[16] == "meaninvvar")
    centermethod = cmeaninvvar;
  else if (op[16] == "meanintegral")
    centermethod = cmeanintegral;
  else if (op[16] == "meanf")
    centermethod = meanf;
  else if (op[16] == "meanfd")
    centermethod = meanfd;
  else if (op[16] == "meansum2")
    centermethod = meansum2;

  f = op[20].strtodouble(binning);

  f = op[59].strtodouble(rankK);

  datanames = vn;
  }


DESIGN_userdefined::DESIGN_userdefined(void) : DESIGN()
  {

  }

void DESIGN_userdefined::init_data(datamatrix & dm, datamatrix & iv)
  {

  unsigned j;

  // 1. Indexsort of data
  if (index_data.rows() <= 1)
    {
    index_data = statmatrix<int>(dm.rows(),1);
    index_data.indexinit();
    dm.indexsort(index_data,0,dm.rows()-1,0,0);
    }

  double dm_mean = dm.mean(0);

  //2. data = sorted observations, init intvar and intvar2
  data = datamatrix(dm.rows(),1);
  double * workdata = data.getV();
  int * workindex = index_data.getV();
  for (j=0;j<dm.rows();j++,workdata++,workindex++)
    {
    *workdata = dm(*workindex,0);
    }

  if (iv.rows() == dm.rows())
    {
    intvar = iv;
    intvar2 = datamatrix(iv.rows(),1);
    double * workintvar2 = intvar2.getV();
    double * workintvar = intvar.getV();

    for (j=0;j<iv.rows();j++,workintvar++,workintvar2++)
      {
      *workintvar2 = pow(*workintvar,2);
      }
    }

  // 3. Creates posbeg, posend

  posbeg.erase(posbeg.begin(),posbeg.end());
  posend.erase(posend.begin(),posend.end());
  posbeg.push_back(0);
  workdata = data.getV()+1;
  double help = data(0,0);
  for(j=1;j<data.rows();j++,workdata++)
    {
    if (  *workdata != help)
      {
      posend.push_back(j-1);
      if (j < data.rows())
        posbeg.push_back(j);
      }

    help = *workdata;

    }

  if (posend.size() < posbeg.size())
    posend.push_back(data.rows()-1);

  // 4. initializes ind
  int k;
  workindex = index_data.getV();
  ind = statmatrix<unsigned>(dm.rows(),1);
  for (j=0;j<posend.size();j++)
    {
    for (k=posbeg[j];k<=posend[j];k++,workindex++)
      ind(*workindex,0) = j;

    }

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\ind.res");
  // ind.prettyPrint(out);
  // TEST


  // 5. Compute meaneffectnr, mclosest, effectvalues
  effectvalues.erase(effectvalues.begin(),effectvalues.end());
  double d;
  meaneffectnr = 0;
  double mclosest = data(posbeg[0],0);
  for(j=0;j<posbeg.size();j++)
    {
    d = data(posbeg[j],0);
    if ( fabs(d-dm_mean) < fabs(mclosest-dm_mean) )
      {
      meaneffectnr = j;
      mclosest = d;
      }

    effectvalues.push_back(ST::doubletostring(d,0));
    }

  compute_meaneffectintvar();
  }

void DESIGN_userdefined::compute_precision(double l)
  {
  if (precisiondeclared==false)
    {
    precision = envmatdouble(K.computeMaxXenv(XWX),1.0,nrpar);
    precisiondeclared = true;
    }

  precision.addto(XWX,K,1.0,l);
  }

void DESIGN_userdefined::compute_Zout(datamatrix & Z)
  {
  unsigned i,j;

  for(i=0; i<posbeg.size() ; i++)
    {
    Zout2.push_back(vector<double>());
    for(j=0; j<Z.cols(); j++)
      {
      if(Z(index_data(posbeg[i],0),j)!=0)
        {
        Zout2[i].push_back(Z(index_data(posbeg[i],0),j));
        }
      }
    }
  }

void DESIGN_userdefined::compute_Zout_transposed(datamatrix & Z)
  {

  vector<double> h;
  ZoutT = vector<vector<double> >(nrpar,h);

  vector<int> h2;
  index_ZoutT = vector<vector<int> >(nrpar,h2);


  unsigned i,j;

  for(i=0;i<Z.cols();i++)
    for(j=0;j<Z.rows();j++)
      {
      if(Z(j,i)!=0)
        {
        ZoutT[i].push_back(Z(j,i));
        index_ZoutT[i].push_back(j);
        }
      }
  }

  // CONSTRUCTOR

DESIGN_userdefined::DESIGN_userdefined(datamatrix & dm,datamatrix & iv,
                       datamatrix & designmat, datamatrix & penmat,
                       GENERAL_OPTIONS * o,DISTR * dp,FC_linear * fcl,
                       vector<ST::string> & op,
                       vector<ST::string> & vn)
                      : DESIGN(o,dp,fcl)
  {

  read_options(op,vn);

  // Adjust to multiple columns in dm
  init_data(dm,iv);

  nrpar = designmat.cols();

  // to be implemented
  compute_Zout(designmat);
  compute_Zout_transposed(designmat);

  K = envmatdouble(penmat, 0.0);
  if(rankK==-1)
    {
    // compute rankK;
    }

  Wsum = datamatrix(posbeg.size(),1,1);
  compute_XtransposedWX();
  XWres = datamatrix(nrpar,1);

  compute_precision(1.0);
  }


  // COPY CONSTRUCTOR

DESIGN_userdefined::DESIGN_userdefined(const DESIGN_userdefined & m)
    : DESIGN(DESIGN(m))
  {
  round = m.round;
  binning = m.binning;
  Zout2 = m.Zout2;
  index_Zout2 = m.index_Zout2;
  }


  // OVERLOADED ASSIGNMENT OPERATOR

const DESIGN_userdefined & DESIGN_userdefined::operator=(const DESIGN_userdefined & m)
  {
  if (this == &m)
    return *this;
  DESIGN::operator=(DESIGN(m));
  round = m.round;
  binning = m.binning;
  Zout2 = m.Zout2;
  index_Zout2 = m.index_Zout2;
  return *this;
  }

void DESIGN_userdefined::outbasis_R(ofstream & out)
  {
  }

void DESIGN_userdefined::outoptions(GENERAL_OPTIONS * op)
  {
  ST::string centerm;

  if (center==false)
    centerm = "uncentered sampling";
  else
    centerm = "centered sampling";

  op->out("  Prior: user defined\n");
  op->out("  Rank of penalty matrix:" + ST::inttostring(rankK) + "\n");
  op->out("  " + centerm + "\n" );

  op->out("\n");

  }

} // end: namespace MCMC


