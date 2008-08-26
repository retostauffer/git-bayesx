
#include "FC_hrandom.h"


//------------------------------------------------------------------------------
//-------------- CLASS: FC_hrandom implementation of member functions ----------
//------------------------------------------------------------------------------


namespace MCMC
{


FC_hrandom::FC_hrandom(void)
  {
  }


void FC_hrandom::read_options(vector<ST::string> & op,vector<ST::string> & vn)
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

  if (op[14] == "increasing")
    stype = increasing;
  else if (op[14] == "decreasing")
    stype = decreasing;
  else
    stype = unconstrained;

  if (op[12] == "true")
    mult = true;
  else
    mult = false;  

  }

FC_hrandom::FC_hrandom(GENERAL_OPTIONS * o,DISTR * lp,DISTR * lp_RE,
                 const ST::string & t,const ST::string & fp,
                 const ST::string & fp2, DESIGN * Dp,
                 vector<ST::string> & op, vector<ST::string> & vn)
     : FC_nonp(o,lp,t,fp,Dp,op,vn)
  {
  read_options(op,vn);
  likep_RE = lp_RE;
  likep_RE->trmult = likep->trmult;
  FCrcoeff = FC(o,t + "_random coefficients",beta.rows(),beta.cols(),fp2);
  }


FC_hrandom::FC_hrandom(const FC_hrandom & m)
  : FC_nonp(FC_nonp(m))
  {
  mult = m.mult;
  likep_RE = m.likep_RE;
  FCrcoeff = m.FCrcoeff;
  }


const FC_hrandom & FC_hrandom::operator=(const FC_hrandom & m)
  {

  if (this==&m)
	 return *this;
  FC_nonp::operator=(FC_nonp(m));
  mult = m.mult;
  likep_RE = m.likep_RE;
  FCrcoeff = m.FCrcoeff;
  logold = m.logold;
  lognew = m.lognew;
  return *this;
  }


void FC_hrandom::set_rcoeff(void)
  {
  unsigned i;
  double * betap = beta.getV();
  double * betarcoeffp = FCrcoeff.beta.getV();


  double * linpredREp;
  if (likep_RE->linpred_current==1)
    linpredREp = likep_RE->linearpred1.getV();
  else
    linpredREp = likep_RE->linearpred2.getV();

  for (i=0;i<beta.rows();i++,betap++,betarcoeffp++,linpredREp++)
    *betarcoeffp = *betap - *linpredREp;

  FCrcoeff.transform(0,0) = transform(0,0);
  }


void FC_hrandom::transform_beta(void)
  {
  if (mult)
    transform(0,0) = 1.0;
  else
    FC_nonp::transform_beta();
  }



void FC_hrandom::update_linpred(int & begin, int & end, double  & value)
  {

  unsigned i;

  double * * linpredp;

  if (likep->linpred_current==1)
    linpredp = designp->linpredp1.getV()+begin;
  else
    linpredp = designp->linpredp2.getV()+begin;

  for (i=begin;i<=end;i++,linpredp++)
    *(*linpredp) += value;

  }




void FC_hrandom::update_IWLS(void)
  {
  unsigned i;
  double logold;

  vector<int>::iterator itbeg = designp->posbeg.begin();
  vector<int>::iterator itend = designp->posend.begin();

  statmatrix<double *> * linpredp;

  if (likep->linpred_current==1)
    linpredp = &(designp->linpredp1);
  else
    linpredp = &(designp->linpredp2);

  if (optionsp->nriter == 1)
    {
    paramold.assign(param);
    betaold.assign(beta);
    lognew = datamatrix(beta.rows(),1);
    logold = datamatrix(beta.rows(),1);    
    }

  double h = likep->compute_iwls(true,false);

  for (i=0;i<beta.rows();i++,++itbeg,++itend)
    {
    logold = likep->loglikelihood(*itbeg,*itend,designp->responsep,
                                  designp->weightp,*linpredp);


    update_linpred( );
    }

  }


void FC_hrandom::update(void)
  {

  if (IWLS)
    {
    update_IWLS();
    }
  else
    {
    FC_nonp::update();
    }

  set_rcoeff();

  FCrcoeff.acceptance++;
  FCrcoeff.update();

  likep_RE->workingresponse.assign(beta);
  likep_RE->response.assign(beta);
  likep_RE->trmult = likep->trmult;
  likep_RE->sigma2 = likep->get_scale();
  }


bool FC_hrandom::posteriormode(void)
  {

  bool conv;
  conv= FC_nonp::posteriormode();

  set_rcoeff();

  bool conv2 = FCrcoeff.posteriormode();

  likep_RE->workingresponse.assign(beta);
  likep_RE->response.assign(beta);
  likep_RE->trmult = likep->trmult;

  // TEST
  /*
  ofstream out5("c:\\bayesx\\test\\results\\fhrandom.res");
  beta.prettyPrint(out5);
  */
  // TEST

  return conv;
  }


void FC_hrandom::outresults(const ST::string & pathresults)
  {

  if (pathresults.isvalidfile() != 1)
    {

    FC::outresults(pathresults);
    FCrcoeff.outresults("");

    optionsp->out("  Results are stored in file\n");
    optionsp->out("  " +  pathresults + "\n");
    optionsp->out("\n");

    ofstream outres(pathresults.strtochar());

    optionsp->out("\n");

    unsigned i;

    ST::string l1 = ST::doubletostring(optionsp->lower1,4);
    ST::string l2 = ST::doubletostring(optionsp->lower2,4);
    ST::string u1 = ST::doubletostring(optionsp->upper1,4);
    ST::string u2 = ST::doubletostring(optionsp->upper2,4);
    l1 = l1.replaceallsigns('.','p');
    l2 = l2.replaceallsigns('.','p');
    u1 = u1.replaceallsigns('.','p');
    u2 = u2.replaceallsigns('.','p');

    outres << "intnr" << "   ";
    outres << designp->datanames[designp->datanames.size()-1] << "   ";
    outres << "pmean_tot   ";

    if (optionsp->samplesize > 1)
      {
      outres << "pqu"  << l1  << "_tot   ";
      outres << "pqu"  << l2  << "_tot   ";
      outres << "pmed_tot   ";
      outres << "pqu"  << u1  << "_tot   ";
      outres << "pqu"  << u2  << "_tot   ";
      outres << "pcat" << optionsp->level1 << "_tot   ";
      outres << "pcat" << optionsp->level2 << "_tot   ";
      }


    outres << "pmean   ";

    if (optionsp->samplesize > 1)
      {
      outres << "pqu"  << l1  << "   ";
      outres << "pqu"  << l2  << "   ";
      outres << "pmed   ";
      outres << "pqu"  << u1  << "   ";
      outres << "pqu"  << u2  << "   ";
      outres << "pcat" << optionsp->level1 << "   ";
      outres << "pcat" << optionsp->level2 << "   ";
      }


    outres << endl;

    double * workmean = betamean.getV();
    double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
    double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
    double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
    double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
    double * workbetaqu50 = betaqu50.getV();

    double * workmean_rcoeff = FCrcoeff.betamean.getV();
    double * workbetaqu_l1_lower_p_rcoeff = FCrcoeff.betaqu_l1_lower.getV();
    double * workbetaqu_l2_lower_p_rcoeff = FCrcoeff.betaqu_l2_lower.getV();
    double * workbetaqu_l1_upper_p_rcoeff = FCrcoeff.betaqu_l1_upper.getV();
    double * workbetaqu_l2_upper_p_rcoeff = FCrcoeff.betaqu_l2_upper.getV();
    double * workbetaqu50_rcoeff = FCrcoeff.betaqu50.getV();

    unsigned nrpar = beta.rows();
    for(i=0;i<nrpar;i++,workmean++,workbetaqu_l1_lower_p++,
                              workbetaqu_l2_lower_p++,workbetaqu50++,
                              workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++,
                              workmean_rcoeff++,workbetaqu_l1_lower_p_rcoeff++,
                              workbetaqu_l2_lower_p_rcoeff++,
                              workbetaqu_l1_upper_p_rcoeff++,
                              workbetaqu_l2_upper_p_rcoeff++,
                              workbetaqu50_rcoeff++)
      {
      outres << (i+1) << "   ";
      outres << designp->effectvalues[i] << "   ";
      outres << *workmean << "   ";

      if (optionsp->samplesize > 1)
        {
        outres << *workbetaqu_l1_lower_p << "   ";
        outres << *workbetaqu_l2_lower_p << "   ";
        outres << *workbetaqu50 << "   ";
        outres << *workbetaqu_l2_upper_p << "   ";
        outres << *workbetaqu_l1_upper_p << "   ";

        if (*workbetaqu_l1_lower_p > 0)
          outres << 1 << "   ";
        else if (*workbetaqu_l1_upper_p < 0)
          outres << -1 << "   ";
        else
          outres << 0 << "   ";

        if (*workbetaqu_l2_lower_p > 0)
          outres << 1 << "   ";
        else if (*workbetaqu_l2_upper_p < 0)
          outres << -1 << "   ";
        else
          outres << 0 << "   ";
        }

      outres << *workmean_rcoeff << "   ";

      if (optionsp->samplesize > 1)
        {
        outres << *workbetaqu_l1_lower_p_rcoeff << "   ";
        outres << *workbetaqu_l2_lower_p_rcoeff << "   ";
        outres << *workbetaqu50_rcoeff << "   ";
        outres << *workbetaqu_l2_upper_p_rcoeff << "   ";
        outres << *workbetaqu_l1_upper_p_rcoeff << "   ";

        if (*workbetaqu_l1_lower_p_rcoeff > 0)
          outres << 1 << "   ";
        else if (*workbetaqu_l1_upper_p_rcoeff < 0)
          outres << -1 << "   ";
        else
          outres << 0 << "   ";

        if (*workbetaqu_l2_lower_p_rcoeff > 0)
          outres << 1 << "   ";
        else if (*workbetaqu_l2_upper_p_rcoeff < 0)
          outres << -1 << "   ";
        else
          outres << 0 << "   ";
        }

      outres << endl;
      }

    }



  }

} // end: namespace MCMC



