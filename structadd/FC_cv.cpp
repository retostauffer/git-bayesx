
#include "FC_cv.h"
#include "clstring.h"

//------------------------------------------------------------------------------
//----------------- CLASS: FC_cv implementation of member functions ------------
//------------------------------------------------------------------------------


namespace MCMC
{



void FC_cv::read_options(vector<ST::string> & op,
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
  */

  }



void FC_cv::get_ind(void)
  {

  vector<FC_hrandom>::iterator it = hrandoms->begin();
  ind = (*it).designp->ind;
  nrcat = (*it).beta.rows();
  effectvalues = (*it).designp->effectvalues;

  }



FC_cv::FC_cv(void)
  {
  }


FC_cv::FC_cv(GENERAL_OPTIONS * o,DISTR * lp,
                 const ST::string & t, const ST::string & fp,
                 vector<FC_hrandom> * FChs)
  : FC(o,t,1,1,fp)
  {

  likep = lp;

  hrandoms = FChs;
  sampled_etas = datamatrix(likep->nrobs,o->compute_samplesize(),0);
  sampled_responses = datamatrix(likep->nrobs,o->compute_samplesize(),0);
  effect = datamatrix(likep->nrobs,1,0);
  linpred = datamatrix(likep->nrobs,1,0);
  size =   hrandoms->size();
  get_ind();
  }


FC_cv::FC_cv(const FC_cv & m)
  : FC(FC(m))
  {
  ind = m.ind;
  nrcat = m.nrcat;
  sampled_etas = m.sampled_etas;
  sampled_responses = m.sampled_responses;
  likep = m.likep;
  hrandoms = m.hrandoms;
  effect = m.effect;
  linpred = m.linpred;
  size = m.size;
  e_score = m.e_score;
  log_score = m.log_score;
  effectvalues = m.effectvalues;
  }


const FC_cv & FC_cv::operator=(
const FC_cv & m)
  {
  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  ind = m.ind;
  nrcat = m.nrcat;
  sampled_etas = m.sampled_etas;
  sampled_responses = m.sampled_responses;
  likep = m.likep;
  hrandoms = m.hrandoms;
  effect = m.effect;
  linpred = m.linpred;
  size = m.size;
  e_score = m.e_score;
  log_score = m.log_score;
  effectvalues = m.effectvalues;
  return *this;
  }


void  FC_cv::update(void)
  {

  if(
     (optionsp->nriter > optionsp->burnin)
     &&
     ((optionsp->nriter-optionsp->burnin-1) % (optionsp->step) == 0)
    )
    {

    unsigned samplesize = optionsp->samplesize;

    unsigned i;

    vector<FC_hrandom>::iterator it = hrandoms->begin();

    if (likep->linpred_current==1)
      linpred.assign(likep->linearpred1);
    else
      linpred.assign(likep->linearpred2);


    for(i=0;i<size;i++,++it)
      {
      (*it).get_effect(effect);

      linpred.minus(linpred,effect);

      (*it).sample_for_cv(effect);

      linpred.plus(linpred,effect);
      }

    double * workse = sampled_etas.getV()+samplesize-1;
    double * worklin = linpred.getV();
    for (i=0;i<sampled_etas.rows();i++,workse+=sampled_etas.cols())
      *workse = likep->trmult*(*worklin);

    likep->sample_responses(samplesize-1,sampled_responses);
    }


  }


bool FC_cv::posteriormode(void)
  {

  return true;
  }


void FC_cv::outoptions(void)
  {

  }


void FC_cv::outresults(ofstream & out_stata, ofstream & out_R,
                            const ST::string & pathresults)
  {


  if (pathresults.isvalidfile() != 1)
    {

    optionsp->out("  CV: \n",true);
    optionsp->out("\n");

    optionsp->out("    Samples for CV are stored in\n");
    optionsp->out("    " +  pathresults + "\n");
    optionsp->out("\n");

    ofstream outres(pathresults.strtochar());

    unsigned nrobs = sampled_etas.rows();
    unsigned i,j;

    for(j=0;j<sampled_etas.cols();j++)
      outres << "s_eta_" << (j+1) << "  ";

    for(j=0;j<sampled_etas.cols();j++)
      outres << "s_resp_" << (j+1) << "  ";

    outres << endl;

    for (i=0;i<nrobs;i++)
      {

      for(j=0;j<sampled_etas.cols();j++)
        outres << sampled_etas(i,j) << "  ";

      for(j=0;j<sampled_responses.cols();j++)
        outres << sampled_responses(i,j) << "  ";

      outres << endl;
      }


    // Energy score

    double es = compute_energyscore();

    ST::string pathresults_e = pathresults.substr(0,pathresults.length()-4)+
                               "_energy.res";

    ofstream out2(pathresults_e.strtochar());

    for (i=0;i<e_score.rows();i++)
      out2 << effectvalues[i] << "  " << e_score(i,0) << endl;

    optionsp->out("    Energy score: " + ST::doubletostring(es,8) + "\n");


    }   // end if (pathresults.isvalidfile() != 1)


  }


double FC_cv::compute_energyscore(void)
  {


  unsigned s,i,j;

  unsigned S = sampled_responses.cols();
  unsigned I = sampled_responses.rows();

  double * srp = sampled_responses.getV();
  double * esp1;
  double * esp2;


  datamatrix es1 = datamatrix(nrcat,S,0);
  datamatrix es2 = datamatrix(nrcat,S,0);

/*
  for (i=0;i<I;i++)
    {
    esp1 = es1.getV()+ind(i,0)*S;
    esp2 = es2.getV()+ind(i,0)*S;

    for (s=0;s<S-1;s++,srp++,esp1++,esp2++)
      {
      *esp1 += pow(sampled_responses(i,s)-likep->response_untransformed(i,0),2);
      *esp2 += pow(sampled_responses(i,s+1)-sampled_responses(i,s),2);
      }
    esp1++;
    *esp1 += pow(sampled_responses(i,s)-likep->response_untransformed(i,0),2);

    }
*/




  unsigned in;

  for (i=0;i<I;i++)
    {

    in = ind(i,0);

    for (s=0;s<S-1;s++)
      {
      es1(in,s) += pow(sampled_responses(i,s)-likep->response_untransformed(i,0),2);
      es2(in,s) += pow(sampled_responses(i,s+1)-sampled_responses(i,s),2);
      }

    es1(in,S-1) += pow(sampled_responses(i,s)-likep->response_untransformed(i,0),2);

    }

  ofstream out("c:\\bayesx\\testh\\results\\es1.res");
  es1.prettyPrint(out);



  if (e_score.rows() != nrcat)
    e_score = datamatrix(nrcat,1,0);

  double h;
  for (j=0;j<nrcat;j++)
    {
    for (s=0;s<S;s++)
      e_score(j,0) += sqrt(es1(j,s));

    e_score(j,0) /= S;

    h=0;
    for(s=0;s<S-1;s++)
      h += sqrt(es2(j,s));

    e_score(j,0) -= h/(2*(S-1));
    }


  return e_score.mean(0);
  }


double FC_cv::compute_logscore(void)
  {

  return 0;
  }


void FC_cv::reset(void)
  {

  }



} // end: namespace MCMC



