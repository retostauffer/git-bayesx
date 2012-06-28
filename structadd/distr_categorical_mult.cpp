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




#include "distr_categorical_mult.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_multinomprobit --------------------
//------------------------------------------------------------------------------


void DISTR_multinomprobit::assign_othercat(DISTR* o)
  {
  othercat.push_back(o);
  nrcat++;
  nrothercat = othercat.size();
  }



void DISTR_multinomprobit::create_responsecat(void)
  {

  responsecat = datamatrix(nrobs,1,nrothercat);

  unsigned i,j;

  bool found;
  for (i=0;i<nrobs;i++)
    {
    found = false;
    j=0;
    while ((found==false) && (j < nrothercat) )
      {
      if (othercat[j]->response(i,0) == 1)
        {
        responsecat(i,0) = j;
        found=true;
        }

      j++;
      }

    if ((found==false) && response(i,0)==1)    // master
      {
//      responsecat(i,0) = nrothercat;
      found=true;
      }

    if (found==false)   // reference
      responsecat(i,0) = -1;

    }

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\responsecat.raw");
  // responsecat.prettyPrint(out);
  // ENDE: TEST

  }


DISTR_multinomprobit::DISTR_multinomprobit(GENERAL_OPTIONS * o,
                                           bool mast,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR(o,r,w)

  {

  master = mast;

  if (master==true)
    {

    nrcat=2;
    nrothercat = 0;

    if (check_weightsone() == true)
      wtype = wweightschange_weightsone;
    else
      wtype = wweightschange_weightsneqone;

    }
  else
    wtype = wweightschange_weightsone;

  family = "Multinomial probit";

  updateIWLS = false;
  }




const DISTR_multinomprobit & DISTR_multinomprobit::operator=(
                                      const DISTR_multinomprobit & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  responsecat = nd.responsecat;
  master=nd.master;
  othercat = nd.othercat;
  nrcat = nd.nrcat;
  nrothercat = nd.nrothercat;
  return *this;
  }


DISTR_multinomprobit::DISTR_multinomprobit(const DISTR_multinomprobit & nd)
   : DISTR(DISTR(nd))
  {
  responsecat = nd.responsecat;
  master=nd.master;
  othercat = nd.othercat;
  nrcat = nd.nrcat;
  nrothercat = nd.nrothercat;
  }


void DISTR_multinomprobit::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: standard normal (probit link)\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTR_multinomprobit::loglikelihood(double * response, double * linpred,
                                     double * weight) const
  {
/*
  if (*weight!=0)
    {
    double mu = randnumbers::Phi2(*linpred);
    if (*response > 0)
      return log(mu);
    else
      return log(1-mu);
    }
  else
*/    return 0;

  }



double DISTR_multinomprobit::loglikelihood_weightsone(
                                  double * response, double * linpred) const
  {
/*
  double mu = randnumbers::Phi2(*linpred);
  if (*response > 0)
    return log(mu);
  else
    return log(1-mu);
  */
  }


void DISTR_multinomprobit::compute_mu(const double * linpred,double * mu)
  {
  // *mu = randnumbers::Phi2(*linpred);
  }


void DISTR_multinomprobit::compute_deviance(const double * response,
                   const double * weight,const double * mu,double * deviance,
                   double * deviancesat, double * scale) const
  {
/*
  if (*weight !=  0)
    {

    if (*response<=0)
      {
      *deviance = -2*log(1-*mu);
      *deviancesat = *deviance;
      }
    else if (*response > 0)
      {
      *deviance = -2*log(*mu);
      *deviancesat = *deviance;
      }

    }
  else
    {
    *deviance = 0;
    *deviancesat = 0;
    }
  */
  }


double DISTR_multinomprobit::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {

  double  mu = randnumbers::Phi2(*linpred);

  double h = 0.39894228*exp(-0.5 * *linpred * *linpred);
  double g = 1/pow(h,2);

  *workingweight = *weight / (mu*(1-mu) * g);


  *workingresponse = *linpred + (*response - mu)/h;

  if (like)
    {

    if (*response > 0)
      return log(mu);
    else
      return log(1-mu);
    }
  else
    {
    return 0;
    }

  }



void DISTR_multinomprobit::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  double  mu = randnumbers::Phi2(*linpred);
  double h = 0.39894228*exp(-0.5 * *linpred * *linpred);
  double g = 1/pow(h,2);

  *workingweight = 1.0 / (mu*(1-mu) * g);

  *workingresponse = *linpred + (*response - mu)/h;

  if (compute_like)
    {

    if (*response > 0)
      like+= log(mu);
    else
      like+= log(1-mu);
    }

  }


void DISTR_multinomprobit::compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  }

void DISTR_multinomprobit::compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  }





double DISTR_multinomprobit::maxutility(vector<datamatrix*> responsep,
const unsigned & i, const unsigned & cat)
  {
  unsigned j;
  double max = 0;
  double help;


  for (j=0;j<=nrothercat;j++)
    {
    help = (*responsep[j])(i,0);
    if ( (j != cat) && (help > max) )
      max = help;
    }

  return max;
  }



void DISTR_multinomprobit::update(void)
  {

  if (optionsp->nriter==1)
    {
    workingweight = weight;
    }

  if (optionsp->nriter==2)
    {

    if (check_weightsone() == true)
      wtype = wweightsnochange_one;
    else
      wtype = wweightsnochange_constant;

    }

  if (master==true)
    {
    unsigned i,j;


    vector<datamatrix *> worklin;
    for (j=0;j<nrothercat;j++)
      {
      if (othercat[j]->linpred_current==1)
        worklin.push_back(&othercat[j]->linearpred1);
      else
        worklin.push_back(&othercat[j]->linearpred2);
      }

    if (linpred_current==1)
      worklin.push_back(&linearpred1);
    else
      worklin.push_back(&linearpred2);


    vector<datamatrix *> responsep;
    for (j=0;j<nrothercat;j++)
      {
      responsep.push_back(&othercat[j]->workingresponse);
      }
    responsep.push_back(&workingresponse);


    double lin;

    for (i=0;i<nrobs;i++)
      {
      if (responsecat(i,0) == -1)   // reference category
        {

        for (j=0;j<=nrothercat;j++)
          {
          lin = (*worklin[j])(i,0);
          (*responsep[j])(i,0) = lin+truncnormal(-20-lin,-lin);
          }

        }
      else
        {
        lin = (*worklin[responsecat(i,0)])(i,0);
        (*responsep[responsecat(i,0)])(i,0) = lin + truncnormal(maxutility(responsep,i,responsecat(i,0)) - lin,20-lin);

        for (j=0;j<=nrothercat;j++)
          {
          if (j != responsecat(i,0))
            {
            lin = (*worklin[j])(i,0);
            (*responsep[j])(i,0) = lin + truncnormal(-20-lin,(*responsep[responsecat(i,0)])(i,0) - lin);
            }
          }

        }

      }

    // TEST
    /*
    ofstream out("c:\\bayesx\\testh\\results\\utility.raw");
    for (i=0;i<nrobs;i++)
      {
      for (j=0;j<nrcat-1;j++)
        out << (*responsep[j])(i,0) << "   ";
      out << endl;
      }
    */
    // TEST

    } // end: if (master==true)


  DISTR::update();

  }



//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_multinomlogit- --------------------
//------------------------------------------------------------------------------
// Logit bisher: ähnlich wie frühwirth eindimensional

DISTR_multinomlogit::DISTR_multinomlogit(GENERAL_OPTIONS * o,
                                         bool mast,
                                         const datamatrix & r,
                                         const datamatrix & w)
 	: DISTR_multinomprobit(o, mast, r, w)
	{

  family = "Multinomial logit";
  H = 6; // vorläufig fix vorgegeben
// Tabellenwerte aus Frühwirth-Schnatter
  SQ = datamatrix(6,5,0); // entspricht 1/s_r^2
  SQ(0,1) = 1/1.2131;
  SQ(1,1) = 1/2.9955;
  SQ(2,1) = 1/7.5458;

  SQ(0,4) = 1/0.68159;
  SQ(1,4) = 1/1.2419;
  SQ(2,4) = 1/2.2388;
  SQ(3,4) = 1/4.0724;
  SQ(4,4) = 1/7.4371;
  SQ(5,4) = 1/13.772;

  //müsste um H = 5 und H= 4  und H = 2 erweitern
  SQ(0,3) = 1/0.79334;
  SQ(1,3) = 1/1.5474;
  SQ(2,3) = 1/3.012;
  SQ(3,3) = 1/5.9224;
  SQ(4,3) = 1/11.77;

  SQ(0,2) = 1/0.95529;
  SQ(1,2) = 1/2.048;
  SQ(2,2) = 1/4.4298;
  SQ(3,2) = 1/9.701;

  SQ(0,0) = 1/1.6927;
  SQ(1,0) = 1/5.2785;


  weights_mixed = datamatrix(6,5,0);
  weights_mixed(0,1) = 0.2522;  // enstpricht w_r
  weights_mixed(1,1) = 0.58523;
  weights_mixed(2,1) = 0.16257;


  weights_mixed(0,4) = 0.018446;
  weights_mixed(1,4) = 0.17268;
  weights_mixed(2,4) = 0.37393;
  weights_mixed(3,4) = 0.31697;
  weights_mixed(4,4) = 0.1089;
  weights_mixed(5,4) = 0.0090745;

  // ebenso hier die erweitern für H= 5 und H = 4 und H = 2

  weights_mixed(0,3) = 0.044333;
  weights_mixed(1,3) = 0.29497;
  weights_mixed(2,3) = 0.42981;
  weights_mixed(3,3) = 0.20759;
  weights_mixed(4,3) = 0.023291;

  weights_mixed(0,2) = 0.1065;
  weights_mixed(1,2) = 0.45836;
  weights_mixed(2,2) = 0.37419;
  weights_mixed(3,2) = 0.060951;

  weights_mixed(0,0) = 0.56442;
  weights_mixed(1,0) = 0.43558;

	}

DISTR_multionomlogit::DISTR_multinomlogit(const DISTR_multinomlogit & nd)
	: DISTR_multinomprobit(DISTR_multinomprobit(nd))
	{
    H = nd.H;
	SQ = nd.SQ;
	weights_mixed = nd.weights_mixed;
	}

const DISTR_multinomlogit & DISTR_multinomlogit::operator=(const DISTR_multinomlogit & nd)
	{
	if (this==&nd)
  	return *this;
  DISTR_multinomprobit::operator=(DISTR_multinomprobit(nd));
  H = nd.H;
  SQ = nd.SQ;
  weights_mixed = nd.weights_mixed;
  return *this;
	}

void DISTR_multinomlogit::outoptions()   //output for logfile
  {
  DISTR::outoptions();
  optionsp->out("  Response function: logistic distribution function\n");
  optionsp->out("  Number of mixture components: " + ST::inttostring(H) + "\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


// Update nach Vorbild fruehwirthschnatter
void DISTR_multinomlogit::update(void)
  {
/*
  double * workresp;  // yi
  double * workwresp; // wki datamatrix workwresp(Anzahlkat,nrobs)
  double * weightwork;
  double * wweightwork; // Gewicht sj mit dem dann beta gesampelt wird

  workresp = response.getV();  // yi
  workwresp = workingresponse.getV(); // wki
  weightwork = weight.getV();
  wweightwork = workingweight.getV(); //

  double * worklin;  // sollte xi*betak sein
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double lambda;
  double lambda_k;
  double U;
  bool kategoriek;
  datamatrix weights_aux(H,1);

    for(int i=0;i<nrobs;i++,worklin++,workresp++,weightwork++,workwresp++,wweightwork++)  // workwresp nicht raufzählen (w(k,i)) weightwork wozu???
     {
      lambda = exp(*worklin);
      lambda_k = sum(exp(*worklin)); // summe worklin außer für k te Kategorie - wie machen?
      U = uniform();
      if (workresp = k)
        kategoriek = 1;
      else
        kategoriek = 0;

      *workwresp(k,i) = log(lambda*U/lambda_k+kategoriek)-log(1-U+lambda/lambda_k*(1+ kategoriek));

      //weights_mixed  // H- 2 da H eigentlich von 2 - 6 läuft und nicht von 0 - 4
      for(int j=0; j < H; j++)
    	  {
        weights_aux(j,0) = weights_mixed(j,H-2)*sqrt(SQ(j,H-2)) * exp(-1/2*  pow((*workwresp - *worklin), 2)*SQ(j,H-2) );
        }

      //distribution function
      for(int j=1; j <H; j++)  // alle zusammenzählen
    	  {
         weights_aux(j,0) = weights_aux(j-1,0) + weights_aux(j,0);
        }

      U = uniform();
      U = U*weights_aux(H-1,0);	//scale to [0, max]


      int iaux = 0; // ziehe rki damit man dann das Gewicht bekommt
      while (U > weights_aux(iaux,0))
        {
        iaux++;
        }

      *wweightwork =  SQ(iaux,H-2); // herauslesen von sj aus der Tabelle von Frühwirth und Schnatter; von k und i abhängig!

      }
*/
    }




} // end: namespace MCMC



