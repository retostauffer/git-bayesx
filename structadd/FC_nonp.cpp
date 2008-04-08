
#include "FC_nonp.h"


//------------------------------------------------------------------------------
//----------------- CLASS: FC implementation of member functions ---------------
//------------------------------------------------------------------------------


namespace MCMC
{


FC_nonp::FC_nonp(void)
  {
  }


FC_nonp::FC_nonp(GENERAL_OPTIONS * o,DISTR * lp,
                 const ST::string & t,const ST::string & fp,
                 DESIGN * Dp)
     : FC(o,t,Dp->Zout.rows(),1,fp)
  {
  likep = lp;
  designp = Dp;
  param = datamatrix(designp->nrpar,1,0);
  betaold = beta;
  betadiff = beta;
  partres = datamatrix(designp->posbeg.size(),1,0);
  lambda=1;
  tau2 = likep->get_scale()/lambda;
  betahelp = beta;
  }


FC_nonp::FC_nonp(const FC_nonp & m)
  : FC(FC(m))
  {

  likep = m.likep;
  designp = m.designp;
  param = m.param;
  betaold = m.betaold;
  betadiff = m.betadiff;
  partres = m.partres;
  lambda=m.lambda;
  tau2 = m.tau2;
  betahelp = m.betahelp;
  }


const FC_nonp & FC_nonp::operator=(const FC_nonp & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  likep = m.likep;
  designp = m.designp;
  param = m.param;
  betaold = m.betaold;
  betadiff = m.betadiff;
  partres = m.partres;
  lambda=m.lambda;
  tau2 = m.tau2;
  betahelp = m.betahelp;
  return *this;
  }



void FC_nonp::update(void)
  {

  bool lambdaconst = false;

  betaold.assign(beta);


  designp->compute_partres(partres,beta);


  if ((likep->changingweight) || (changingdesign))
    designp->compute_XtransposedWX_XtransposedWres(partres,lambda);
  else
    designp->compute_XtransposedWres(partres, lambda);

  if ((likep->changingweight) || (changingdesign) || (!lambdaconst))
    designp->compute_precision(lambda);


  double sigmaresp = sqrt(likep->get_scale());

  double * work = betahelp.getV();
  unsigned i;
  unsigned nrpar = beta.rows();
  for(i=0;i<nrpar;i++,work++)
    *work = sigmaresp*rand_normal();

  designp->precision.solveU(betahelp);

  designp->precision.solve(designp->XWres,betahelp,param);



//  if(center)
//    centerparam();

  designp->compute_f(param,beta);

  betadiff.minus(beta,betaold);


  designp->update_linpred(betadiff,true);

  acceptance++;

  transform(0,0) = likep->trmult;

  FC::update();
  }


bool FC_nonp::posteriormode(void)
  {

  betaold.assign(beta);

  bool lambdaconst = false;

  designp->compute_partres(partres,beta);

  /*
  // TEST
  ofstream out("c:\\bayesx\\test\\results\\partres.res");
  partres.prettyPrint(out);
  // TEST
  */

  if ((likep->changingweight) || (changingdesign))
    designp->compute_XtransposedWX_XtransposedWres(partres, lambda);
  else
    designp->compute_XtransposedWres(partres, lambda);

  if ((likep->changingweight) || (changingdesign) || (!lambdaconst))
    designp->compute_precision(lambda);

  // TEST
  // ofstream out1("c:\\bayesx\\test\\results\\precision.res");
  // designp->precision.print1(out1);
  // TEST


  designp->precision.solve(designp->XWres,param);

  // TEST
  // ofstream out2("c:\\bayesx\\test\\results\\param.res");
  // param.prettyPrint(out2);
  // TEST

//  if(center)
//    centerparam();

  designp->compute_f(param,beta);

  betadiff.minus(beta,betaold);

  designp->update_linpred(betadiff,true);

  transform(0,0) = likep->trmult;
  return FC::posteriormode();

  }



void FC_nonp::outresults(const ST::string & pathresults)
  {

  if (pathresults.isvalidfile() != 1)
    {

    FC::outresults(pathresults);

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
    for (i=0;i<designp->datanames.size();i++)
      outres << designp->datanames[i] << "   ";
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


//    unsigned j;
    unsigned nrpar = beta.rows();
    for(i=0;i<nrpar;i++,workmean++,workbetaqu_l1_lower_p++,
                              workbetaqu_l2_lower_p++,workbetaqu50++,
                              workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++)
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

      outres << endl;
      }

    }

  }


void FC_nonp::centerparam(void)
  {

  unsigned i;
  double sum=0;
  double * workparam = param.getV();
  unsigned nrparam = param.rows();

  for (i=0;i<nrparam;i++,workparam++)
    {
    sum+= *workparam;
    }

  workparam = param.getV();

  sum /= double(nrparam);

  for (i=0;i<nrparam;i++,workparam++)
    *workparam-= sum;

  }

  
void FC_nonp::reset(void)
  {

  }


} // end: namespace MCMC



