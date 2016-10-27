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



#include "FC_merror.h"


//------------------------------------------------------------------------------
//-------------- CLASS: FC_hrandom implementation of member functions ----------
//------------------------------------------------------------------------------


namespace MCMC
{


FC_merror::FC_merror(void)
  {
  }

void FC_merror::read_options(vector<ST::string> & op,vector<ST::string> & vn)
  {
  int f;
  f = op[20].strtodouble(binning);
  }

FC_merror::FC_merror(GENERAL_OPTIONS * o, const ST::string & t,
            const ST::string & fp, vector<ST::string> & op,
            vector<ST::string> & vn, datamatrix & xo, datamatrix & mv,
            datamatrix & xd, FC_nonp * fcn)
     : FC(o,t,xo.rows(),1,fp)
  {
  read_options(op,vn);

  xobs = xo;
  merror = (double)(xo.cols());
  FCp = fcn;
  mevar = mv;
  mesd = mv;
  for(unsigned i=0; i< mesd.rows(); i++)
    mesd(i,0) = sqrt(mesd(i,0));
  xmean = xd;

  setbeta(xd);

  FCp->designp->changingdesign=true;

  minbin = xobs.min();
  maxbin = xobs.max();
  deltabin = (maxbin-minbin)/binning;

  a_tau2_x = 0.001;
  b_tau2_x = 0.001;
  m_mu_x = 0.0;
  s_mu_x = 1000.0*1000.0;
  }

FC_merror::FC_merror(const FC_merror & m)
  : FC(FC(m))
  {
  xobs = m.xobs;
  xmean = m.xmean;
  FCp = m.FCp;
  mevar = m.mevar;
  binning = m.binning;
  minbin = m.minbin;
  maxbin=m.maxbin;
  deltabin = m.deltabin;
  merror = m.merror;
  a_tau2_x = m.a_tau2_x;
  b_tau2_x = m.b_tau2_x;
  m_mu_x = m.m_mu_x;
  s_mu_x = m.s_mu_x;
  mu_x = m.mu_x;
  tau2_x = m.tau2_x;
  }

const FC_merror & FC_merror::operator=(const FC_merror & m)
  {
  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  xobs = m.xobs;
  xmean = m.xmean;
  FCp = m.FCp;
  mevar = m.mevar;
  binning = m.binning;
  minbin = m.minbin;
  maxbin=m.maxbin;
  deltabin = m.deltabin;
  merror = m.merror;
  a_tau2_x = m.a_tau2_x;
  b_tau2_x = m.b_tau2_x;
  m_mu_x = m.m_mu_x;
  s_mu_x = m.s_mu_x;
  mu_x = m.mu_x;
  tau2_x = m.tau2_x;
  return *this;
  }

void FC_merror::update(void)
  {
  unsigned i;
  double u1 = minbin+deltabin/2;
  double h;
  double prop;

  double sqrtM = sqrt(merror);
  double lognew, logold;

  for(i=0; i<xmean.rows(); i++)
    {
    // generate proposal
    prop = beta(i,0) + 2*mesd(i,0)*randnumbers::rand_normal()/sqrtM;
    // start binning of proposal
    if(prop < minbin)
      prop = minbin;
    if(prop > maxbin)
      prop = maxbin;
    h = floor((prop - minbin)/deltabin);
    if (h >= binning)
      h -= 1.0;
    prop = u1+h*deltabin;
    // end binning of proposal


    }




  statmatrix<int> countmat((unsigned)binning,1,0);

  double dmrw;

  for(i=0; i<xmean.rows(); i++)
    {
    dmrw = xmean(i,0) + 0.01*randnumbers::rand_normal();
    if(dmrw < minbin)
      dmrw = minbin;
    if(dmrw > maxbin)
      dmrw = maxbin;

    h = floor((dmrw - minbin)/deltabin);
    if (h >= binning)
      h -= 1.0;

    countmat((unsigned)h,0)++;
    xmean(i,0) = u1+h*deltabin;
    }


  //

  unsigned j;

  // 1. Indexsort of data
  FCp->designp->index_data.indexinit();
  xmean.indexsort(FCp->designp->index_data,0,xmean.rows()-1,0,0);

  //2. data = sorted observations
  double * workdata = FCp->designp->data.getV();
  int * workindex = FCp->designp->index_data.getV();
  for (j=0;j<xmean.rows();j++,workdata++,workindex++)
    {
    *workdata = xmean(*workindex,0);
    }

  // 3. Creates posbeg, posend
  int countsum=0;
  for(j=0; j<FCp->designp->posbeg.size(); j++)
    {
    if(countmat(j,0)!=0)
      {
      FCp->designp->posbeg[j] = countsum;
      FCp->designp->posend[j] = countsum + countmat(j,0)-1;
      countsum += countmat(j,0);
      }
    else
      {
      FCp->designp->posbeg[j] = -1;
      FCp->designp->posend[j] = -1;
      }
    }

  // 4. initializes ind
  int k;
  workindex = FCp->designp->index_data.getV();
  for (j=0;j<FCp->designp->posend.size();j++)
    {
    if(FCp->designp->posbeg[j]!=-1)
      {
      for (k=FCp->designp->posbeg[j];k<=FCp->designp->posend[j];k++,workindex++)
        FCp->designp->ind(*workindex,0) = j;
      }
    }

//  FCp->designp->update_covs_merror(xmean);
  }

bool FC_merror::posteriormode(void)
  {
  }

void FC_merror::outresults(ofstream & out_stata,ofstream & out_R, ofstream & out_R2BayesX,
                            const ST::string & pathresults)
  {
  }

} // end: namespace MCMC




