
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
  }


FC_cv::FC_cv(const FC_cv & m)
  : FC(FC(m))
  {
  sampled_etas = m.sampled_etas;
  sampled_responses = m.sampled_responses;
  likep = m.likep;
  hrandoms = m.hrandoms;
  effect = m.effect;
  linpred = m.linpred;
  size = m.size;
  }


const FC_cv & FC_cv::operator=(
const FC_cv & m)
  {
  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  sampled_etas = m.sampled_etas;
  sampled_responses = m.sampled_responses;
  likep = m.likep;
  hrandoms = m.hrandoms;
  effect = m.effect;
  linpred = m.linpred;
  size = m.size;
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
      *workse = *worklin;

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

    }   // end if (pathresults.isvalidfile() != 1)

  }


void FC_cv::reset(void)
  {

  }



} // end: namespace MCMC



