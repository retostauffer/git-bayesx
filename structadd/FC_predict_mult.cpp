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



#include "FC_predict_mult.h"
#include "clstring.h"

//------------------------------------------------------------------------------
//--------- CLASS: FC_predict_mult implementation of member functions ----------
//------------------------------------------------------------------------------


namespace MCMC
{



FC_predict_mult::FC_predict_mult(void)
  {
  }


FC_predict_mult::FC_predict_mult(GENERAL_OPTIONS * o,vector<DISTR *> lp,
                                 const ST::string & t, const ST::string & fp,
                                 const ST::string & fpd, datamatrix & dm,
                                 vector<ST::string> & dn)
  : FC(o,t,1,1,fp)
  {
  likep = lp;
  designmatrix= dm;
  varnames = dn;

  setbeta((likep[0])->nrobs,likep.size()*2,0);
//  setbeta(1000,likep.size()*2,0);

  FC_deviance = FC(o,"",2,1,fpd);

  }


FC_predict_mult::FC_predict_mult(const FC_predict_mult & m)
  : FC(FC(m))
  {
  likep = m.likep;
  designmatrix = m.designmatrix;
  varnames = m.varnames;
  FC_deviance = m.FC_deviance;
  deviance = m.deviance;
  deviancesat = m.deviancesat;
  }


const FC_predict_mult & FC_predict_mult::operator=(const FC_predict_mult & m)
  {
  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  likep = m.likep;
  designmatrix = m.designmatrix;
  varnames = m.varnames;
  FC_deviance = m.FC_deviance;
  deviance = m.deviance;
  deviancesat = m.deviancesat;
  return *this;
  }


void  FC_predict_mult::update(void)
  {

  get_predictor();

  acceptance++;

  FC::update();

  FC_deviance.beta(0,0) = deviance;
  FC_deviance.beta(1,0) = deviancesat;
  FC_deviance.acceptance++;
  FC_deviance.update();

  }



void FC_predict_mult::get_predictor(void)
  {

    unsigned i,j;

    double * betap = beta.getV();

    vector<double *> worklinp;
    vector<double *> workresponse;
    vector<double *> workweight;
    vector<double>   scalehelp;

    for (j=0;j<likep.size();j++)
      {
      if (likep[j]->linpred_current==1)
        worklinp.push_back(likep[j]->linearpred1.getV());
      else
        worklinp.push_back(likep[j]->linearpred2.getV());

      workresponse.push_back(likep[j]->response_untransformed.getV());

      workweight.push_back(likep[j]->weight.getV());

      scalehelp.push_back(likep[j]->get_scale());
      }

    deviance=0;
    deviancesat=0;
    double deviancehelp;
    double deviancesathelp;


    for(i=0;i<likep[0]->nrobs;i++)
      {

      for (j=0;j<likep.size();j++,betap++)
        {
        *betap = *worklinp[j];
        }


      for (j=0;j<likep.size();j++,betap++)
        {
        likep[j]->compute_mu_mult(worklinp,betap);

        worklinp[j]++;
        }



      }

   // TEST
   // ofstream out("c:\\bayesx\\testh\\results\\beta.raw");
   // beta.prettyPrint(out);
   // out.close();
   // END TEST

  }


bool FC_predict_mult::posteriormode(void)
  {

  get_predictor();

  posteriormode_betamean();

  return true;
  }


void FC_predict_mult::outoptions(void)
  {

  }



void FC_predict_mult::outresults_DIC(const ST::string & pathresults)
  {
/*
  ST::string pathresultsdic = pathresults.substr(0,pathresults.length()-4) + "_DIC.res";
  ofstream out(pathresultsdic.strtochar());

  double deviance2=0;
  double deviance2_sat=0;


  double scalehelp = likep->get_scalemean();


  double devhelp;
  double devhelp_sat;

  double mu_meanlin;
  double * workmeanlin = betamean.getV();
  double * workresponse = likep->response_untransformed.getV();
  double * workweight = likep->weight.getV();

  unsigned i;
  for (i=0;i<likep->nrobs;i++,workmeanlin+=beta.cols(),workresponse++,workweight++)
    {

    likep->compute_mu(workmeanlin,&mu_meanlin);

    likep->compute_deviance(workresponse,workweight,&mu_meanlin,
    &devhelp,&devhelp_sat,&scalehelp);

    deviance2 += devhelp;
    deviance2_sat += devhelp_sat;
    }


  double devhelpm = FC_deviance.betamean(0,0);

  unsigned d;
  if (devhelpm > 1000000000)
    d = 14;
  else if (devhelpm > 1000000)
    d = 11;
  else
    d = 8;

  out << "deviance   pd dic" << endl;

  optionsp->out("  ESTIMATION RESULTS FOR THE DIC: \n",true);
  optionsp->out("\n");

  optionsp->out("    DIC based on the unstandardized deviance\n");
  optionsp->out("\n");

  optionsp->out("    Deviance(bar_mu):           " +
  ST::doubletostring(deviance2,d) + "\n");
  out << deviance2 << "   ";

  optionsp->out("    pD:                         " +
  ST::doubletostring(devhelpm-deviance2,d) + "\n");
  out << (devhelpm-deviance2) << "   ";

  optionsp->out("    DIC:                        " +
  ST::doubletostring(2*devhelpm-deviance2,d) + "\n");
  optionsp->out("\n");
  out << (2*devhelpm-deviance2) << "   " << endl;

  optionsp->out("    DIC based on the saturated deviance\n");
  optionsp->out("\n");


  double devhelpm_sat = FC_deviance.betamean(1,0);

  optionsp->out("    Deviance(bar_mu):           " +
  ST::doubletostring(deviance2_sat,d) + "\n");

  optionsp->out("    pD:                         " +
  ST::doubletostring(devhelpm_sat-deviance2_sat,d) + "\n");

  optionsp->out("    DIC:                        " +
  ST::doubletostring(2*devhelpm_sat-deviance2_sat,d) + "\n");

  optionsp->out("\n");
  */
  }


void FC_predict_mult::outresults_deviance(void)
    {
    /*
    ST::string l1 = ST::doubletostring(optionsp->lower1,4);
    ST::string l2 = ST::doubletostring(optionsp->lower2,4);
    ST::string u1 = ST::doubletostring(optionsp->upper1,4);
    ST::string u2 = ST::doubletostring(optionsp->upper2,4);
    l1 = l1.replaceallsigns('.','p');
    l2 = l2.replaceallsigns('.','p');
    u1 = u1.replaceallsigns('.','p');
    u2 = u2.replaceallsigns('.','p');


    ST::string meanstr = "    Mean:          ";
    unsigned l_meanstr = meanstr.length();

    ST::string stdstr =  "    Std. Dev:      ";
    unsigned l_stdstr = stdstr.length();

    ST::string l1str = "    " + l1 + "% Quantile: ";
    unsigned l_l1str = l1str.length();

    ST::string l2str = "    " + l2 + "% Quantile: ";
    unsigned l_l2str = l2str.length();

    ST::string medianstr = "    50% Quantile: ";
    unsigned l_medianstr = medianstr.length();

    ST::string u1str = "    " + u1 + "% Quantile: ";
    unsigned l_u1str = u1str.length();

    ST::string u2str = "    " + u2 + "% Quantile: ";
    unsigned l_u2str = u2str.length();


    optionsp->out("  ESTIMATION RESULT FOR THE DEVIANCE: \n",true);
    optionsp->out("\n");

    optionsp->out("    Unstandardized Deviance (-2*Loglikelihood(y|mu))\n");
    optionsp->out("\n");

    double devhelpm = FC_deviance.betamean(0,0);
    double devhelp;

    unsigned d;
    if (devhelpm > 1000000000)
      d = 14;
    else if (devhelpm > 1000000)
      d = 11;
    else
      d = 8;

    optionsp->out(meanstr + ST::string(' ',20-l_meanstr) +
    ST::doubletostring(devhelpm,d) + "\n");


    devhelp = sqrt(FC_deviance.betavar(0,0));
    optionsp->out(stdstr + ST::string(' ',20-l_stdstr) +
    ST::doubletostring(devhelp,d) +  "\n");

    devhelp = FC_deviance.betaqu_l1_lower(0,0);
    optionsp->out(l1str +  ST::string(' ',20-l_l1str) +
    ST::doubletostring(devhelp,d) +  "\n");

    devhelp = FC_deviance.betaqu_l2_lower(0,0);
    optionsp->out(l2str +  ST::string(' ',20-l_l2str) +
    ST::doubletostring(devhelp,d) +  "\n");

    devhelp = FC_deviance.betaqu50(0,0);
    optionsp->out(medianstr +  ST::string(' ',20-l_medianstr) +
    ST::doubletostring(devhelp,d) +  "\n");


    devhelp = FC_deviance.betaqu_l2_upper(0,0);
    optionsp->out(u1str +  ST::string(' ',20-l_u1str) +
    ST::doubletostring(devhelp,d) +  "\n");


    devhelp = FC_deviance.betaqu_l1_upper(0,0);
    optionsp->out(u2str +  ST::string(' ',20-l_u2str) +
    ST::doubletostring(devhelp,d) +  "\n");

    optionsp->out("\n");


    optionsp->out("  Saturated Deviance (-2*Loglikelihood(y|mu) + 2*Loglikelihood(y|mu=y))\n");
    optionsp->out("\n");

    devhelpm = FC_deviance.betamean(1,0);

    optionsp->out(meanstr + ST::string(' ',20-l_meanstr) +
    ST::doubletostring(devhelpm,d) + "\n");

    devhelp = sqrt(FC_deviance.betavar(1,0));
    optionsp->out(stdstr + ST::string(' ',20-l_stdstr) +
    ST::doubletostring(devhelp,d) +  "\n");

    devhelp = FC_deviance.betaqu_l1_lower(1,0);;
    optionsp->out(l1str +  ST::string(' ',20-l_l1str) +
    ST::doubletostring(devhelp,d) +  "\n");

    devhelp = FC_deviance.betaqu_l2_lower(1,0);
    optionsp->out(l2str +  ST::string(' ',20-l_l2str) +
    ST::doubletostring(devhelp,d) +  "\n");

    devhelp = FC_deviance.betaqu50(1,0);
    optionsp->out(medianstr +  ST::string(' ',20-l_medianstr) +
    ST::doubletostring(devhelp,d) +  "\n");

    devhelp = FC_deviance.betaqu_l2_upper(1,0);
    optionsp->out(u1str +  ST::string(' ',20-l_u1str) +
    ST::doubletostring(devhelp,d) +  "\n");

    devhelp = FC_deviance.betaqu_l1_upper(1,0);
    optionsp->out(u2str +  ST::string(' ',20-l_u2str) +
    ST::doubletostring(devhelp,d) +  "\n");

    optionsp->out("\n");

    optionsp->out("\n");
     */
    }



void FC_predict_mult::outresults(ofstream & out_stata, ofstream & out_R,
                            const ST::string & pathresults)
  {

  if (pathresults.isvalidfile() != 1)
    {

    FC::outresults(out_stata,out_R,"");

    /*
    if (likep->maindistribution == true)
      {
      FC_deviance.outresults(out_stata,out_R,"");
      }
    */

    optionsp->out("  PREDICTED VALUES: \n",true);
    optionsp->out("\n");

    optionsp->out("    Results for the predictor, mean are stored in file\n");
    optionsp->out("    " +  pathresults + "\n");
    optionsp->out("\n");

    ofstream outres(pathresults.strtochar());

    optionsp->out("\n");

    unsigned i,j;

    ST::string l1 = ST::doubletostring(optionsp->lower1,4);
    ST::string l2 = ST::doubletostring(optionsp->lower2,4);
    ST::string u1 = ST::doubletostring(optionsp->upper1,4);
    ST::string u2 = ST::doubletostring(optionsp->upper2,4);
    l1 = l1.replaceallsigns('.','p');
    l2 = l2.replaceallsigns('.','p');
    u1 = u1.replaceallsigns('.','p');
    u2 = u2.replaceallsigns('.','p');

    outres << "intnr" << "   ";

    for (i=0;i<varnames.size();i++)
      outres << varnames[i] << "   ";

    for (i=0;i<likep.size();i++)
      {
      outres << "pmean_pred_" << i << "   ";

      if (optionsp->samplesize > 1)
        {
        outres << "pqu"  << l1  << "_pred_" << i << "   ";
        outres << "pqu"  << l2  << "_pred_" << i << "   ";
        outres << "pmed_pred_" << i << "   ";
        outres << "pqu"  << u1  << "_pred_" << i << "   ";
        outres << "pqu"  << u2  << "_pred_" << i << "   ";
        }
      }

    for (i=0;i<likep.size();i++)
      {
      outres << "pmean_mu_" << i << "   ";

      if (optionsp->samplesize > 1)
        {
        outres << "pqu"  << l1  << "_mu_" << i << "   ";
        outres << "pqu"  << l2  << "_mu_" << i << "   ";
        outres << "pmed_mu_" << i << "   ";
        outres << "pqu"  << u1  << "_mu_" << i << "   ";
        outres << "pqu"  << u2  << "_mu_" << i << "   ";
        }
      }

    outres << endl;

    double * workmean = betamean.getV();
    double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
    double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
    double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
    double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
    double * workbetaqu50 = betaqu50.getV();


    for(i=0;i<designmatrix.rows();i++)
      {

      outres << (i+1) << "   ";

      for (j=0;j<designmatrix.cols();j++)
        outres << designmatrix(i,j) << "   ";

      for (j=0;j<likep.size();j++,
          workmean++,
          workbetaqu_l1_lower_p++,
          workbetaqu_l2_lower_p++,
          workbetaqu50++,
          workbetaqu_l1_upper_p++,
          workbetaqu_l2_upper_p++)
        {
        outres << *workmean << "   ";

        if (optionsp->samplesize > 1)
          {
          outres << *workbetaqu_l1_lower_p << "   ";
          outres << *workbetaqu_l2_lower_p << "   ";
          outres << *workbetaqu50 << "   ";
          outres << *workbetaqu_l2_upper_p << "   ";
          outres << *workbetaqu_l1_upper_p << "   ";
          }

      // mu
        workmean++;
        workbetaqu_l1_lower_p++;
        workbetaqu_l2_lower_p++;
        workbetaqu50++;
        workbetaqu_l1_upper_p++;
        workbetaqu_l2_upper_p++;

        outres << *workmean << "   ";

        if (optionsp->samplesize > 1)
          {
          outres << *workbetaqu_l1_lower_p << "   ";
          outres << *workbetaqu_l2_lower_p << "   ";
          outres << *workbetaqu50 << "   ";
          outres << *workbetaqu_l2_upper_p << "   ";
          outres << *workbetaqu_l1_upper_p << "   ";
          }
      // end mu
      }

      outres << endl;

     }

//      outresults_deviance();
//      outresults_DIC(pathresults);


    }   // end if (pathresults.isvalidfile() != 1)

  }


void FC_predict_mult::reset(void)
  {

  }



} // end: namespace MCMC


