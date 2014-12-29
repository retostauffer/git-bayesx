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



#include "FC_nonp_variance.h"
//#include "gsl_randist.h"
//#include "gsl_cdf.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//---------- CLASS: FC_nonp_variance implementation of member functions ---------
//------------------------------------------------------------------------------

void FC_nonp_variance::read_options(vector<ST::string> & op,
                                   vector<ST::string> & vn)
  {

  /*
  1       degree
  2       numberknots
  3       difforder
  4       lambda
  5       a
  6       b
  */

  int f;

  f = op[4].strtodouble(lambdastart);
  f = op[5].strtodouble(a_invgamma);
  f = op[6].strtodouble(b_invgamma_orig);

  if (op[38] == "true")
    {
    lambdaconst = true;
    nosamples =true;
    }


  f = op[47].strtodouble(tildea);
  f = op[48].strtodouble(tildeb);
  f = op[51].strtodouble(scaletau2);
  if (op[49] == "true")
    {
    cauchy = true;
    }
  else
    cauchy = false;

  if (op[50] == "true")
    {
    wei = true;
    }
  else
    wei = false;

  if (op[58] == "iwls_tau")
    {
    proposal = 2;
    }
  else if(op[58] == "gamma")
    {
    proposal = 1;
    }
  else if(op[58] == "iwls_logtau2")
    {
    proposal = 3;
    }
  else
    {
    proposal = 0;
    }

  }





FC_nonp_variance::FC_nonp_variance(void)
  {

  }



FC_nonp_variance::FC_nonp_variance(MASTER_OBJ * mp, unsigned & enr, GENERAL_OPTIONS * o,
                 DISTR * lp, const ST::string & t,const ST::string & fp,
                 DESIGN * Dp,FC_nonp * FCn,vector<ST::string> & op,
                 vector<ST::string> & vn)
     : FC(o,t,1,2,fp)
  {
  FCnonpp = FCn;
  likep = lp;
  designp = Dp;
  masterp = mp;
  equationnr = enr,
  lambdaconst = false;

  read_options(op,vn);

  datamatrix betanew(1,2);
  betanew(0,0) = likep->get_scale()/lambdastart;
  betanew(0,1) = lambdastart;
  setbeta(betanew);

  FCnonpp->tau2 = beta(0,0);
  FCnonpp->lambda = beta(0,1);

  }


FC_nonp_variance::FC_nonp_variance(const FC_nonp_variance & m)
    : FC(FC(m))
  {
  FCnonpp = m.FCnonpp;
  likep = m.likep;
  designp = m.designp;
  masterp = m.masterp;
  equationnr = m.equationnr;
  a_invgamma = m.a_invgamma;
  b_invgamma_orig = m.b_invgamma_orig;
  b_invgamma = m.b_invgamma;
  lambdastart = m.lambdastart;
  lambdaconst = m.lambdaconst;
  tildea = m.tildea;
  tildeb = m.tildeb;
  cauchy = m.cauchy;
  wei = m.wei;
  scaletau2 = m.scaletau2;
  proposal = m.proposal;
  }


const FC_nonp_variance & FC_nonp_variance::operator=(const FC_nonp_variance & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  FCnonpp = m.FCnonpp;
  likep = m.likep;
  designp = m.designp;
  masterp = m.masterp;
  equationnr = m.equationnr;
  a_invgamma = m.a_invgamma;
  b_invgamma_orig = m.b_invgamma_orig;
  b_invgamma = m.b_invgamma;
  lambdastart = m.lambdastart;
  lambdaconst = m.lambdaconst;
  tildea = m.tildea;
  tildeb = m.tildeb;
  cauchy = m.cauchy;
  wei = m.wei;
  scaletau2 = m.scaletau2;
  proposal = m.proposal;
  return *this;
  }


void FC_nonp_variance::update(void)
  {

  // TEST
  //  ofstream out("c:\\bayesx\\test\\results\\param.res");
  //  (FCnonpp->param).prettyPrint(out);
  // END: TEST

   /* double clk = (double)CLK_TCK;
    clock_t clkbegin2 = clock();

    double counti = 0;
    double nrs = 1;*/
//    double tmp = 0;
//    std::ofstream out2;
//    out2.open ("C:\\tmp\\GIG2.raw", std::ofstream::out | std::ofstream::app);
    //out2 << "GIG2"  << endl;
   // for (counti=0;counti<nrs;counti++)
  //  {
//    tmp = randnumbers::GIG2(0, 0.5, 1.5);
//    out2 << tmp << endl;
  //  }
  //  out2 << endl;
  //  clock_t clkend2 = clock();
  //  double sec2 = (clkend2-clkbegin2)/clk;

   // cout << "time GIG2: " << sec2 << endl;

   // clock_t clkbegin = clock();
   // counti = 0;
 //   std::ofstream out;
//    out.open ("C:\\tmp\\GIG.raw", std::ofstream::out | std::ofstream::app);
    //out << "GIG"  << endl;
    //for (counti=0;counti<nrs;counti++)
    //{
//    tmp = randnumbers::GIG(0, 0.5, 1.5);
 //   out << tmp << endl;
    //}
    //out << endl;

    /*clock_t clkend = clock();
    double sec = (clkend-clkbegin)/clk;

    cout << "time GIG: " << sec << endl;*/

  b_invgamma = masterp->level1_likep[equationnr]->trmult*b_invgamma_orig;

  if (lambdaconst == false)
    {

    if (cauchy == true)
      {
      double quadf = designp->penalty_compute_quadform(FCnonpp->param);
      double gamma = rand_invgamma(designp->rankK/2 +tildea ,0.5*quadf+tildeb);

      double u = log(uniform());

      double fcold = -(0.5*designp->rankK+0.5)*log(beta(0,0))-1/(2*beta(0,0))*quadf-log(1+beta(0,0)) ;
      double fcnew = -(0.5*designp->rankK+0.5)*log(gamma)-1/(2*gamma)*quadf-log(1+gamma);
      double proposalold = -(tildea+designp->rankK/2+1)*log(beta(0,0)) - (0.5*quadf+tildeb) / beta(0,0);
      double proposalnew = -(tildea+designp->rankK/2+1)*log(gamma) - (0.5*quadf+tildeb) / gamma;
      if (u <= (fcnew - fcold - proposalnew + proposalold))
        {

        beta(0,0) = gamma;
        acceptance++;
        }

      }
    else if (wei == true)
      {
      double quadf = designp->penalty_compute_quadform(FCnonpp->param);
      double u = log(uniform());

      if(proposal == 0) //simply draw proposal as if tau^2 would be IG distributed with tildea and tildeb
        {
        double gamma = rand_invgamma(designp->rankK/2 +tildea ,0.5*quadf+tildeb);

        double fcold = -(0.5*designp->rankK+0.5)*log(beta(0,0))-pow(beta(0,0)/scaletau2, 0.5)-1/(2*beta(0,0))*quadf;
        double fcnew = -(0.5*designp->rankK+0.5)*log(gamma)-pow(gamma/scaletau2, 0.5)-1/(2*gamma)*quadf;
        double proposalold = -(designp->rankK/2 +tildea+1)*log(beta(0,0)) - (0.5*quadf+tildeb) / beta(0,0);
        double proposalnew = -(designp->rankK/2 +tildea+1)*log(gamma) - (0.5*quadf+tildeb) / gamma;

        if (u <= (fcnew - fcold - proposalnew + proposalold))
          {

          beta(0,0) = gamma;
          acceptance++;
          }
        }
      else if(proposal == 1) // gamma approximation of weibull distribution
        {
        double tmp1 = 3.0;
        double tmp2 = 5.0;
        double ew = randnumbers::gamma_exact(tmp1) *scaletau2;
        double varw = scaletau2 * scaletau2 * (randnumbers::gamma_exact(tmp2) - pow(randnumbers::gamma_exact(tmp1), 2));
        double bgamma = ew / varw;
        double agamma = bgamma * ew;
        double p = -designp->rankK/2 + agamma;
        double a = 2 * bgamma;
        double b = quadf;
        double gamma = randnumbers::GIG2(p, a, b);

        double fcold = -(0.5*designp->rankK+0.5)*log(beta(0,0))-pow(beta(0,0)/scaletau2, 0.5)-1/(2*beta(0,0))*quadf;
        double fcnew = -(0.5*designp->rankK+0.5)*log(gamma)-pow(gamma/scaletau2, 0.5)-1/(2*gamma)*quadf;
        double proposalold = (p-1)*log(beta(0,0)) - 0.5*(a*beta(0,0)+b/beta(0,0));
        double proposalnew = (p-1)*log(gamma) - 0.5*(a*gamma+b/gamma);

        if (u <= (fcnew - fcold - proposalnew + proposalold))
          {

          beta(0,0) = gamma;
          acceptance++;
          }
        }
      else if(proposal == 2)// iwls proposal for tau
        {
        double vartau = 1/ ((3 * quadf / beta(0,0) - designp->rankK) / beta(0,0));
        double mutau = sqrt(beta(0,0)) + vartau * (-designp->rankK / sqrt(beta(0,0)) + quadf / pow(beta(0, 0), 1.5) - 1/sqrt(scaletau2));

        double gamma = mutau + rand_normal() * sqrt(vartau);
        double fcold = -(0.5*designp->rankK)*log(beta(0,0))-pow(beta(0,0)/scaletau2, 0.5)-1/(2*beta(0,0))*quadf;
        double fcnew = -(designp->rankK)*log(gamma)-pow(gamma*gamma/scaletau2, 0.5)-1/(2*gamma*gamma)*quadf;

        double vartauold = 1/ ((3 * quadf / (gamma*gamma) - designp->rankK) / (gamma*gamma));
        double mutauold = gamma + vartau * (-designp->rankK / gamma + quadf / pow(gamma, 3) - 1/sqrt(scaletau2));
        double proposalold = -0.5*log(vartauold)-0.5*pow((sqrt(beta(0,0))-mutauold), 2)/vartauold;
        double proposalnew = -0.5*log(vartau)-0.5*pow((gamma-mutau), 2)/vartau;

        if (u <= (fcnew - fcold - proposalnew + proposalold))
          {
          gamma = gamma*gamma;
          beta(0,0) = gamma;
          acceptance++;
          }
        }
      else // iwls proposal for log(tau^2)
        {
        double vartau = 1/(0.5*quadf/beta(0,0) + 0.25*sqrt(beta(0,0))/sqrt(scaletau2));
        double mutau = log(beta(0,0)) + vartau * (1 - 0.5*(designp->rankK+1) + 0.5*quadf/beta(0,0) - 0.5*sqrt(beta(0,0))/sqrt(scaletau2));

        double gamma = mutau + rand_normal() * sqrt(vartau);
        double fcold = log(beta(0,0)) - 0.5*(designp->rankK+1)*log(beta(0,0)) - 0.5*pow(beta(0,0)/scaletau2, 0.5) - 1/(2*beta(0,0))*quadf;
        double fcnew = gamma- 0.5*(designp->rankK+1)*gamma - 0.5*pow(exp(gamma)/scaletau2, 0.5) - 1/(2*exp(gamma))*quadf;

        double vartauold = 1/(0.5*quadf/exp(gamma) + 0.25*sqrt(exp(gamma))/sqrt(scaletau2));
        double mutauold = gamma + vartau * (1 - 0.5*(designp->rankK+1) + 0.5*quadf/exp(gamma) - 0.5*sqrt(exp(gamma))/sqrt(scaletau2));
        double proposalold = -0.5*log(vartauold)-0.5*pow((log(beta(0,0))-mutauold), 2)/vartauold;
        double proposalnew = -0.5*log(vartau)-0.5*pow((gamma-mutau), 2)/vartau;

        if (u <= (fcnew - fcold - proposalnew + proposalold))
          {
          gamma = exp(gamma);
          beta(0,0) = gamma;
          acceptance++;
          }
        }
      }
    else
      {
      beta(0,0) = rand_invgamma(a_invgamma+0.5*designp->rankK,
                  b_invgamma+0.5*designp->penalty_compute_quadform(FCnonpp->param));
      acceptance++;
      }

    beta(0,1) = likep->get_scale()/beta(0,0);

    FCnonpp->tau2 = beta(0,0);

    FC::update();

    }

  }


bool FC_nonp_variance::posteriormode(void)
  {

  b_invgamma = masterp->level1_likep[equationnr]->trmult*b_invgamma_orig;

  beta(0,0) = likep->get_scale()/beta(0,1);

  FCnonpp->tau2 = beta(0,0);

  posteriormode_betamean();

  return true;
  }



void FC_nonp_variance::outresults(ofstream & out_stata,ofstream & out_R,
                                  const ST::string & pathresults)
  {

  FC::outresults(out_stata,out_R,"");

  ST::string l1 = ST::doubletostring(optionsp->lower1,4);
  ST::string l2 = ST::doubletostring(optionsp->lower2,4);
  ST::string u1 = ST::doubletostring(optionsp->upper1,4);
  ST::string u2 = ST::doubletostring(optionsp->upper2,4);

  ST::string nl1 = ST::doubletostring(optionsp->lower1,4);
  ST::string nl2 = ST::doubletostring(optionsp->lower2,4);
  ST::string nu1 = ST::doubletostring(optionsp->upper1,4);
  ST::string nu2 = ST::doubletostring(optionsp->upper2,4);
  nl1 = nl1.replaceallsigns('.','p');
  nl2 = nl2.replaceallsigns('.','p');
  nu1 = nu1.replaceallsigns('.','p');
  nu2 = nu2.replaceallsigns('.','p');

  ST::string vstr;


  if (optionsp->samplesize > 1)
    {

    FC::outresults_acceptance();

    vstr = "    Mean:         ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betamean(0,0),6) + "\n");

    vstr = "    Std. dev.:    ";

    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(sqrt(betavar(0,0)),6) + "\n");

    vstr = "    " + l1 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l1_lower(0,0),6) + "\n");

    vstr = "    " + l2 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l2_lower(0,0),6) + "\n");

    vstr = "    50% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu50(0,0),6) + "\n");

    vstr = "    " + u1 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l2_upper(0,0),6) + "\n");

    vstr = "    " + u2 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l1_upper(0,0),6) + "\n");

    optionsp->out("\n");

    }
  else
    {
    optionsp->out("    Smoothing parameter: " +
    ST::doubletostring(betamean(0,1),6) + "\n");

    optionsp->out("\n");
    }

//  out_R << "term=" << title <<  ";" << endl;

  if (pathresults.isvalidfile() != 1)
    {

    optionsp->out("    Results for the variance component are also stored in file\n");
    optionsp->out("    " +  pathresults + "\n");
    optionsp->out("\n");

    ST::string paths = pathresults.substr(0,pathresults.length()-4) +
                                 "_sample.raw";

    out_R << "pathvarsample=" << paths << endl;
//    out_R << "filetype=param; path=" << pathresults << ";" <<  endl;

    ofstream ou(pathresults.strtochar());

    if (optionsp->samplesize > 1)
      {
      ou << "pmean  pstd  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
      nu1 << "   pqu" << nu2 << endl;
      }
    else
      {
      ou << "pmean" << endl;
      }

    ou << betamean(0,0) << "  ";
    if (optionsp->samplesize > 1)
      {
      if (betavar(0,0) < 0.0000000000001)
        ou << 0 << "  ";
      else
        ou << sqrt(betavar(0,0)) << "  ";
      ou << betaqu_l1_lower(0,0) << "  ";
      ou << betaqu_l2_lower(0,0) << "  ";
      ou << betaqu50(0,0) << "  ";
      ou << betaqu_l2_upper(0,0) << "  ";
      ou << betaqu_l1_upper(0,0) << "  " << endl;
      }

    optionsp->out("\n");
    }

  }


void FC_nonp_variance::outoptions(void)
  {
  if (cauchy)
    {
    optionsp->out("  Cauchy prior\n");

    optionsp->out("  Hyperparameter tildea for proposal density: " +
                ST::doubletostring(tildea) + "\n" );

    optionsp->out("  Hyperparameter tildeb for proposal density: " +
                ST::doubletostring(tildeb) + "\n" );

    }
  else if (wei)
    {
    optionsp->out("  Weibull prior\n");

    optionsp->out("  Scale parameter: " +
                ST::doubletostring(scaletau2) + "\n" );

    if(proposal == 0)
      {
      optionsp->out("  Inverse gamma proposal density \n" );
      optionsp->out("  Hyperparameter tildea for proposal density: " +
                ST::doubletostring(tildea) + "\n" );

      optionsp->out("  Hyperparameter tildeb for proposal density: " +
                ST::doubletostring(tildeb) + "\n" );
      }
    else if(proposal == 1)
      {
      optionsp->out("  Generalised inverse Gaussian proposal density \n" );
      }
    else if(proposal == 2)
      {
      optionsp->out("  IWLS proposal density for tau \n" );
      }
    else
      {
      optionsp->out("  IWLS proposal density for log(tau^2) \n" );
      }
    }
  else
    {

    optionsp->out("  Inverse gamma prior\n");

    b_invgamma = masterp->level1_likep[equationnr]->trmult*b_invgamma_orig;

    optionsp->out("  Hyperprior a for variance parameter: " +
                ST::doubletostring(a_invgamma) + "\n" );
    optionsp->out("  Hyperprior b for variance parameter: " +
                ST::doubletostring(b_invgamma) + "\n" );
    optionsp->out("\n");
    }
  }


void FC_nonp_variance::reset(void)
  {

  datamatrix betanew(1,2);
  betanew(0,0) = likep->get_scale()/lambdastart;
  betanew(0,1) = lambdastart;
  setbeta(betanew);

  FCnonpp->tau2 = beta(0,0);
  FCnonpp->lambda = beta(0,1);
//  transform(0,0) = 1;
//  transform(1,0) = 1;

  }



//------------------------------------------------------------------------------
//--- CLASS: FC_nonp_variance_varselection implementation of member functions --
//------------------------------------------------------------------------------

void FC_nonp_variance_varselection::read_options(vector<ST::string> & op,
                                   vector<ST::string> & vn)
  {
  FC_nonp_variance::read_options(op,vn);

  int f;

  f = op[39].strtodouble(a_omega);
  f = op[40].strtodouble(b_omega);

  f = op[41].strtodouble(r);
  f = op[52].strtodouble(r2);

  f = op[53].strtodouble(v1);
  f = op[54].strtodouble(v2);
  f = op[55].strtodouble(tildev1);
  f = op[56].strtodouble(tildev2);
  if (op[50] == "true")
    wei = true;
  else
    wei = false;

  if (op[57] == "true")
    gig = true;
  else
    gig = false;
  }


FC_nonp_variance_varselection::FC_nonp_variance_varselection(void)
  {

  }



FC_nonp_variance_varselection::FC_nonp_variance_varselection(MASTER_OBJ * mp,
                 unsigned & enr, GENERAL_OPTIONS * o,
                 DISTR * lp, bool so, const ST::string & t,const ST::string & fp,
                 DESIGN * Dp,FC_nonp * FCn,vector<ST::string> & op,
                 vector<ST::string> & vn)
     : FC_nonp_variance(mp,enr,o,lp,t,fp,Dp,FCn,op,vn)
  {

  read_options(op,vn);

  singleomega = so;

//  FC_delta = FC(o,"",1,1,"");
//  FC_delta.setbeta(1,1,0);
  FC_delta = FC(o,"",1,2,"");
  FC_delta.setbeta(1,1,0);
  FC_delta.setbeta(1,2,0.5);

  FC_psi2 = FC(o,"",1,1,"");
  FC_psi2.setbeta(1,1,0.5);

  //tauold = rand_normal()*sqrt(0.5);

  if(!singleomega)
    {
    FC_omega = FC(o,"",1,1,"");
    FC_omega.setbeta(1,1,0.5);
    }
  else
    {
    omega=0.5;
    }
  }


FC_nonp_variance_varselection::FC_nonp_variance_varselection(const FC_nonp_variance_varselection & m)
    : FC_nonp_variance(FC_nonp_variance(m))
  {
  singleomega = m.singleomega;
  FC_delta = m.FC_delta;
  FC_psi2 = m.FC_psi2;
  FC_omega = m.FC_omega;
  a_omega = m.a_omega;
  b_omega = m.b_omega;
  omega = m.omega;
  tauold = m.tauold;
  wei = m.wei;
  v1 = m.v1;
  v2 = m.v2;
  tildev1 = m.tildev1;
  tildev2 = m.tildev2;
  r = m.r;
  r2 = m.r2;
  X = m.X;
  gig = m.gig;
  }


const FC_nonp_variance_varselection & FC_nonp_variance_varselection::operator=(const FC_nonp_variance_varselection & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  singleomega = m.singleomega;
  FC_delta = m.FC_delta;
  FC_psi2 = m.FC_psi2;
  FC_omega = m.FC_omega;
  a_omega = m.a_omega;
  b_omega = m.b_omega;
  omega = m.omega;
  tauold = m.tauold;
  wei = m.wei;
  v1 = m.v1;
  v2 = m.v2;
  tildev1 = m.tildev1;
  tildev2 = m.tildev2;
  r = m.r;
  r2 = m.r2;
  X = m.X;
  gig = m.gig;
  return *this;
  }


void FC_nonp_variance_varselection::update(void)
  {

  unsigned i;

  // updating psi2

  //if more than one omega, set omegas
  if(singleomega == false)
    {
    omega = FC_omega.beta(0,0) ;
    }

  double r_delta;

  if(wei==true)
    {
    if (FC_delta.beta(0,0) < 1)
      r_delta = r2;
    else
      r_delta = 1;

    double gamma = rand_invgamma(tildev1+0.5,tildev2+0.5*beta(0,0)/r_delta);
    //double gamma = rand_invgamma(v+0.5,Q/r_delta);
    double logu = log(uniform());

    double fcold = -log(FC_psi2.beta(0,0))-pow(FC_psi2.beta(0,0)/scaletau2, 0.5)-beta(0,0)/(2*r_delta*FC_psi2.beta(0,0));
    double fcnew = -log(gamma)-pow(gamma/scaletau2, 0.5)-beta(0,0)/(2*r_delta*gamma);
    double proposalold = -(tildev1+0.5+1)*log(FC_psi2.beta(0,0)) - (tildev2+0.5*FC_psi2.beta(0,0)/r_delta) / FC_psi2.beta(0,0);
    double proposalnew = -(tildev1+0.5+1)*log(gamma) - (tildev2+0.5*beta(0,0)/r_delta) / gamma;

    if (logu <= (fcnew - fcold + proposalold - proposalnew))
      {

      FC_psi2.beta(0,0) = gamma;
      FC_psi2.acceptance++;
      }

    FC_psi2.update();

    // end: updating psi2

    // updating delta

    double u = uniform();
    double L = 1/sqrt(r)*exp(- beta(0,0)/(2*FC_psi2.beta(0,0))*(1/r2-1));
    double pr1 = 1/(1+ ((1-omega)/omega)*L);

    if (u <=pr1)
      {
      FC_delta.beta(0,0) = 1;
      r_delta = 1;
      }
    else
      {
      FC_delta.beta(0,0) = 0;
      r_delta = r2;
      }
    FC_delta.beta(0,1) = pr1;

    FC_delta.update();

    //end: updating delta
    }
  else
   {
    if (FC_delta.beta(0,0) < 1)
      r_delta = r;
    else
      r_delta = 1;

    FC_psi2.beta(0,0) = rand_invgamma(v1+0.5,v2+0.5*beta(0,0)/r_delta);

    FC_psi2.update();


    // end: updating psi2

    // updating delta

    double u = uniform();
    double L = 1/sqrt(r)*exp(- beta(0,0)/(2*FC_psi2.beta(0,0))*(1/r-1));
    double pr1 = 1/(1+ ((1-omega)/omega)*L);

    if (u <=pr1)
      {
      FC_delta.beta(0,0) = 1;
      r_delta = 1;
      }
    else
      {
      FC_delta.beta(0,0) = 0;
      r_delta = r;
      }

    FC_delta.beta(0,1) = pr1;
    FC_delta.update();

  // end: updating delta
   }


  // updating w

  if(singleomega == false) {

    FC_omega.beta(0,0) = randnumbers::rand_beta(a_omega+FC_delta.beta(0,0),
                                          b_omega+1-FC_delta.beta(0,0));
//    omega = FC_omega.beta(0,0);

    FC_omega.update();
  }

   // end: updating w

   if(gig)
     {

   // updating tau2

    //double w = designp->penalty_compute_quadform(FCnonpp->param)/(r_delta*FC_psi2.beta(0,0));
     double tau2 = randnumbers::GIG2(0, 1/(r_delta*FC_psi2.beta(0,0)), designp->penalty_compute_quadform(FCnonpp->param));

     beta(0,0) = tau2;
     beta(0,1) = likep->get_scale()/beta(0,0);
     FCnonpp->tau2 = beta(0,0);

     acceptance++;
     FC::update();
     }
    else
    {

  //variante stefan, für GAUSS fall
    double Sigmatau;
    double mutau;

    FCnonpp->designp->compute_effect(X,FCnonpp->beta);

    double * worklin;
    if (likep->linpred_current==1)
      worklin = likep->linearpred1.getV();
    else
      worklin = likep->linearpred2.getV();

    double * Xp = X.getV();
    double * responsep = likep->workingresponse.getV();
    double varinv = 1/(likep->get_scale()*beta(0,0));
    double xtx=0;
    for (i=0;i<X.rows();i++,Xp++,responsep++,worklin++)
      {
      xtx += pow(*Xp,2);
      mutau += (*Xp) * ((*responsep) - (*worklin)+(*Xp));
      }
    Sigmatau = 1/(varinv*xtx + 1/(r_delta*FC_psi2.beta(0,0)));
    mutau *= Sigmatau/(likep->get_scale()*sqrt(beta(0,0)));

    double tau = mutau + sqrt(Sigmatau) * rand_normal();

    double tau2 = tau*tau;
    if (tau2 < 0.000000001)
      tau2 = 0.000000001;

    beta(0,0) = tau2;

    beta(0,1) = likep->get_scale()/beta(0,0);

    FCnonpp->tau2 = beta(0,0);

    // end: updating tau2

    acceptance++;
    FC::update();

//---------------------------------------------------------------------------------------
/*
//  extension to GAMLSS case

  Sigmatau = 1/(FCnonpp->designp->XWX_p->compute_quadform(FCnonpp->beta, 0) +
                      1/(r_delta*FC_psi2.beta(0,0)));
  double * p1 = FCnonpp->beta.getV();
  double * p2 = FCnonpp->designp->XWres_p->getV();
  for(i=0; i<FCnonpp->beta.rows(); i++, p1++, p2++)
    {
    mutau += *p1 * *p2;
    }

  mutau *= Sigmatau;

  double stand_new = rand_normal();
  double tau = mutau + sqrt(Sigmatau) * stand_new;
  double tau2 = tau*tau;
  if (tau2 < 0.000000001)
    tau2 = 0.000000001;

  double tau2old = tauold*tauold;

  u = log(uniform());

  double fcold = likep->compute_iwls(true, true) - 0.5*pow(tauold,2)/(r_delta*FC_psi2.beta(0,0)) +
                 log(randnumbers::phi(stand_new));

  // ToDO: calculate fcnew!

  double fcnew = 0;

  //calculate proposal evaluated at new tau2 with parameter  based on tauold2
  double qold = get_log_proposal();

  //calculate proposal evaluated at tauold with parameters based on tau2
  double qnew = get_log_proposal();

  if (u <= (fcnew - fcold + qold - qnew))
    {
    acceptance++;
    }
  else
    {
    tau = tauold;
    tau2 = tau2old;
    }

  beta(0,0) = tau;
  beta(0,1) = likep->get_scale()/beta(0,0);
  FCnonpp->tau2 = beta(0,0);

  FC::update();*/

  // end: updating tau2
    }
  }


bool FC_nonp_variance_varselection::posteriormode(void)
  {
  bool t = FC_nonp_variance::posteriormode();

  FC_psi2.beta(0,0) = beta(0,0);

  return true;
  }



void FC_nonp_variance_varselection::outresults(ofstream & out_stata,ofstream & out_R,
                                  const ST::string & pathresults)
  {

  if (pathresults.isvalidfile() != 1)
    {

    ST::string pathresults_delta = pathresults.substr(0,pathresults.length()-4) + "_delta.res";
    ST::string pathresults_omega = pathresults.substr(0,pathresults.length()-4) + "_omega.res";
    ST::string pathresults_psi2 = pathresults.substr(0,pathresults.length()-4) + "_psi2.res";

    if(singleomega == false)
    {
        FC_omega.outresults(out_stata,out_R,pathresults_omega);
    }


    FC_nonp_variance::outresults(out_stata,out_R,pathresults);


    optionsp->out("    Inclusion probability: " + ST::doubletostring(FC_delta.betamean(0,0),6)  + "\n");
    optionsp->out("\n");
	optionsp->out("    Rao-Blackwellised inclusion probability: " + ST::doubletostring(FC_delta.betamean(0,1),6)  + "\n");
    optionsp->out("\n");
    optionsp->out("    Results for the inclusion probabilities are also stored in file\n");
    optionsp->out("    " +  pathresults_delta + "\n");
    optionsp->out("\n");
    optionsp->out("\n");

    if(singleomega == false)
    {
        optionsp->out("    Inclusion probability parameter omega:\n");
        optionsp->out("\n");
        FC_omega.outresults_singleparam(out_stata,out_R,"");
        optionsp->out("    Results for the inclusion probability parameter omega are also stored in file\n");
        optionsp->out("    " +  pathresults_omega + "\n");
        optionsp->out("\n");
        optionsp->out("\n");
    }

    // deltas
    ofstream ou(pathresults_delta.strtochar());

    ou << "pmean_delta " << "pmean_prob" << endl;
    ou << FC_delta.betamean(0,0) << " " << FC_delta.betamean(0,1) << endl;
   /* ou<< "pmean_delta ";
    ou << "pmean_prob" << endl;
    ou << FC_delta.betamean(0,0) ;
    ou << " ";
    ou << FC_delta.betamean(0,1) << endl;*/

    FC_psi2.outresults(out_stata,out_R,pathresults_psi2);


    optionsp->out("    Psi2: " + ST::doubletostring(FC_psi2.betamean(0,0),6)  + "\n");
    optionsp->out("\n");

    if(wei==true)
      {
      FC_psi2.outresults_acceptance();
      optionsp->out("\n");
      optionsp->out("    Results for psi2 are also stored in file\n");
      optionsp->out("    " +  pathresults_psi2 + "\n");
      optionsp->out("\n");
      optionsp->out("\n");
      double rate;
      if (nrtrials == 0)
        {
        rate = (double(FC_psi2.acceptance)/double(optionsp->nriter))*100;
        }
      else
        {
        rate = (double(FC_psi2.acceptance)/double(nrtrials))*100;
        }
      ST::string pathresults_psi2_acceptance = pathresults.substr(0,pathresults.length()-4) + "_psi2_acceptance.res";
      ofstream ou2(pathresults_psi2_acceptance.strtochar());

      ou2 << "acceptance ";
      ou2 << "r2" << endl;
      ou2 << rate ;
      ou2 << " ";
      ou2 << r2 << endl;
      }
    else
      {
      optionsp->out("\n");
      optionsp->out("    Results for psi2 are also stored in file\n");
      optionsp->out("    " +  pathresults_psi2 + "\n");
      optionsp->out("\n");
      optionsp->out("\n");
      }

    }
  }


void FC_nonp_variance_varselection::get_samples(
   const ST::string & filename,ofstream & outg) const
  {
  FC_nonp_variance::get_samples(filename,outg);

  ST::string filename_delta = filename.substr(0,filename.length()-11) + "_delta_sample.raw";
  FC_delta.get_samples(filename_delta,outg);

  if(singleomega == false)
  {
    ST::string filename_omega = filename.substr(0,filename.length()-11) + "_omega_sample.raw";
    FC_omega.get_samples(filename_omega,outg);

  }
  ST::string filename_psi2 = filename.substr(0,filename.length()-11) + "_psi2_sample.raw";
  FC_psi2.get_samples(filename_psi2,outg);
  }

void FC_nonp_variance_varselection::outoptions(void)
  {
  if (wei)
    {
    optionsp->out("  Weibull prior\n");

    optionsp->out("  Scale parameter: " +
                ST::doubletostring(scaletau2) + "\n" );

    optionsp->out("  Hyperparameter tildev1 for prior: " +
                ST::doubletostring(tildev1) + "\n" );

    optionsp->out("  Hyperparameter tildev2 for prior: " +
                ST::doubletostring(tildev2) + "\n" );
    }
  else
    {
    FC_nonp_variance::outoptions();
    }
  if(singleomega==false)
    {
    optionsp->out("\n");
    optionsp->out("  Inclusion probability omega\n");
    optionsp->out("  Parameter a_omega: " +
                ST::doubletostring(a_omega) + "\n" );

    optionsp->out("  Parameter b_omega: " +
                ST::doubletostring(b_omega) + "\n" );

    }
  }


void FC_nonp_variance_varselection::reset(void)
  {

  FC_nonp_variance::reset();

  }

//-------------------------------------------------------------------------------------------------------------
//------------------------------------  FC_varselection_omega -------------------------------------------------
//-------------------------------------------------------------------------------------------------------------



FC_varselection_omega::FC_varselection_omega(void)
  {

  }


FC_varselection_omega::FC_varselection_omega(MASTER_OBJ * mp,unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,
           const ST::string & t)
     : FC(o,t,1,1,"")

  {
  likep = lp;
  masterp = mp;

  setbeta(1,1,0.5);
  a_omega = 1;
  b_omega = 1;
  }


FC_varselection_omega::FC_varselection_omega(const FC_varselection_omega & m)
    : FC(FC(m))
   {
   FC_tau2s = m.FC_tau2s;
   a_omega = m.a_omega;
   b_omega = m.b_omega;
   likep = m.likep;
   masterp = m.masterp;
   }


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_varselection_omega & FC_varselection_omega::operator=(const FC_varselection_omega & m)
   {

   if (this==&m)
   	 return *this;
   FC::operator=(FC(m));
   FC_tau2s = m.FC_tau2s;
   a_omega = m.a_omega;
   b_omega = m.b_omega;
   likep = m.likep;
   masterp = m.masterp;
   return *this;
   }


  void FC_varselection_omega::update(void)
    {
        //update omega when having one single omega for each equation.
        double sumdelta=0;
        unsigned nfunc = FC_tau2s.size();
        unsigned j;
        for (j=0;j<nfunc;j++)
        {
            sumdelta +=  FC_tau2s[j]->FC_delta.beta(0,0);
        }

        //draw new omega
        beta(0,0) = randnumbers::rand_beta(a_omega+sumdelta,
                                          b_omega+nfunc-sumdelta);


        for (j=0;j<nfunc;j++)
        {
        FC_tau2s[j]->omega=beta(0,0);
        }


        FC::update();
    }

//bool FC_varselection_omega::posteriormode(void)
//  {
//
//  }
//
//
  void FC_varselection_omega::outoptions(void)
    {
    optionsp->out("  Inclusion probability omega\n");
    optionsp->out("  Parameter a_omega: " +
                ST::doubletostring(a_omega) + "\n" );

    optionsp->out("  Parameter b_omega: " +
                ST::doubletostring(b_omega) + "\n" );
    }


  void FC_varselection_omega::outresults(ofstream & out_stata,ofstream & out_R,
                  const ST::string & pathresults)
   {
        ST::string pathresults_omega = pathresults.substr(0,pathresults.length()-4) + "_omega.res";
        optionsp->out("    Inclusion probability parameter omega: " + ST::doubletostring(FC::betamean(0,0),6)  + "\n");
        optionsp->out("\n");
        FC::outresults(out_stata,out_R,pathresults_omega);

        optionsp->out("\n");
        optionsp->out("    Results for the inclusion probability parameter omega are also stored in file\n");
        optionsp->out("    " +  pathresults_omega + "\n");
        optionsp->out("\n");
        optionsp->out("\n");
   }

  void FC_varselection_omega::reset(void)
    {
        FC::reset();
    }


//  void FC_varselection_omega::read_options(vector<ST::string> & op,vector<ST::string> & vn)
//    {
//
//    }
//

  void FC_varselection_omega::get_samples(const ST::string & filename,ofstream & outg) const
   {
        ST::string filename_omega = filename.substr(0,filename.length()-11) + "_omega_sample.raw";
        FC::get_samples(filename_omega,outg);

   }



} // end: namespace MCMC



