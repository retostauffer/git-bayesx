
#include "FC_nonp.h"


//------------------------------------------------------------------------------
//----------------- CLASS: FC implementation of member functions ---------------
//------------------------------------------------------------------------------


namespace MCMC
{


void FC_nonp::read_options(vector<ST::string> & op,vector<ST::string> & vn)
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
  */

  if (op[14] == "increasing")
    stype = increasing;
  else if (op[14] == "decreasing")
    stype = decreasing;
  else
    stype = unconstrained;

  if (op[18] == "true")
    pvalue = true;
  else
    pvalue = false;
  }


FC_nonp::FC_nonp(void)
  {
  }


FC_nonp::FC_nonp(GENERAL_OPTIONS * o,DISTR * lp,
                 const ST::string & t,const ST::string & fp,
                 DESIGN * Dp,vector<ST::string> & op,
                 vector<ST::string> & vn,bool sstore)
     : FC(o,t,Dp->Zout.rows(),1,fp,sstore)
  {
  read_options(op,vn);
  likep = lp;
  designp = Dp;
  param = datamatrix(designp->nrpar,1,0);
  paramold = param;
  parammode = param;
  paramhelp = param;
  betaold = beta;
  betadiff = beta;
  partres = datamatrix(designp->posbeg.size(),1,0);
  lambda=1;
  tau2 = likep->get_scale()/lambda;
  IWLS = likep->updateIWLS;

  if (Dp->position_lin != -1)
    {
    fsample = FC(o,"",beta.rows(),beta.cols(),fp + ".lin",sstore);
    paramlin = datamatrix(Dp->designlinear.cols(),1,0);
    }

  if (pvalue==true)
    {
    pvalue_sample = FC(o,"",param.rows()*2+6,1,fp + ".pvalue",sstore);
    }

   paramsample = FC(o,"",param.rows(),1,fp + ".param",sstore);

  }


void FC_nonp::update_pvalue(void)
  {

  if( (optionsp->nriter > optionsp->burnin) &&
      ((optionsp->nriter-optionsp->burnin-1) % (optionsp->step) == 0) )
    {
    unsigned i,j;
    unsigned nrpar = param.rows();

    double * contourp = pvalue_sample.beta.getV();
    double * parammodep = parammode.getV();

    for(i=0;i<nrpar;i++,contourp++,parammodep++)
      *contourp = *parammodep;

    *contourp = 1/likep->get_scale();                                // contour(,nrpar) : 1/sigma^2(nu)

    // TEST
    // ofstream out("c:\\bayesx\\testh\\results\\beta_pvalue.raw");
    // pvalue_sample.beta.prettyPrint(out);
    // END: TEST

    contourp++;
    *contourp = lambda;                                              // contour(,nrpar+1) : lambda(nu)

    contourp++;
    *contourp =                                                      // contour(,nrpar+2) : beta(t)' X'WX beta(t)
    designp->XWX.compute_quadform(param,0);
    contourp++;
    *contourp = designp->K.compute_quadform(param,0);                // contour(,nrpar+3) : beta(t)' K beta(t)
    contourp++;
    *contourp = designp->precision.compute_quadform(parammode,0);    // contour(,nrpar+4) : m(nu)' P(nu) m(nu)
    contourp++;
    *contourp = designp->precision.getLogDet();                      // contour(,nrpar+5)  : log(det(P(v))


    datamatrix mPhelp = datamatrix(nrpar,1,0);
    for(i=0;i<nrpar;i++)
      {
      parammodep = parammode.getV();
      for(j=0;j<nrpar;j++,parammodep++)
        mPhelp(i,0) += (*parammodep)*(designp->precision)(j,i);
      }

    contourp++;
    for(i=0;i<nrpar;i++,contourp++)
      *contourp = mPhelp(i,0);

  // TEST
  //   ofstream out("c:\\bayesx\\testh\\results\\beta_pvalue.res");
  //   pvalue_sample.beta.prettyPrint(out);
  // END: TEST

    }

  pvalue_sample.update();

  }


void FC_nonp::compute_pvalue(ST::string & pathresults)
  {

  unsigned nrpar = param.rows();

  unsigned k,t,nu;

  datamatrix mPhelp(nrpar,1,0);

  datamatrix * contour;
  datamatrix contourm;

  if (pvalue_sample.samplestore==true)
    {
    contourm = datamatrix(optionsp->samplesize,pvalue_sample.beta.rows(),0);
    pvalue_sample.readsample3(contourm);
    contour = &contourm;
    }
  else
    contour = &pvalue_sample.sampled_beta;

  // TEST
  // ofstream out2("c:\\bayesx\\testh\\results\\contour.res");
  // contour.prettyPrint(out2);
  // END: TEST

  double exponent,mPbeta;

  double pbeta_0;
  datamatrix pbeta_t(optionsp->samplesize,1,0);

  datamatrix RB(optionsp->samplesize,1,0);

  datamatrix * parameter;
  datamatrix parameterm;

  if (samplestore==true)
    {
    parameterm=datamatrix(optionsp->samplesize,nrpar,0);  // samplesize x nrpar
    paramsample.readsample3(parameterm);
    parameter = &parameterm;
    }
  else
    {
    parameter = &paramsample.sampled_beta;
    }


  // TEST
  // ofstream out3("c:\\bayesx\\testh\\results\\parameter.res");
  // parameter.prettyPrint(out3);
  // TEST


// 0-nrpar-1 : pmode
// nrpar     : 1/scale
// nrpar+1   : lambda
// nrpar+2   : param' XWX param
// nrpar+3   : param'K param
// nrpar+4   : parammode' precision parammode
// nrpar+5   : logdetprecison
// rest      : m'P

  double * paramp;
  double * pbeta_tp = pbeta_t.getV();
  for(t=0;t<optionsp->samplesize;t++,pbeta_tp++)
    {
    for(nu=0;nu<optionsp->samplesize;nu++)
      {
      mPbeta = 0.0;
      paramp = (*parameter).getV()+t*nrpar;
      for(k=0;k<nrpar;k++,paramp++)
        mPbeta += (*contour)(nu,nrpar+6+k)* (*paramp);
//        mPbeta += (*contour)(nu,nrpar+6+k)* (*parameter)(t,k);

      exponent =  (*contour)(t,nrpar+2)
               + (*contour)(nu,nrpar+1) * (*contour)(t,nrpar+3)
               + (*contour)(nu,nrpar+4) - 2*mPbeta/transform(0,0);
      RB(nu,0) = 0.5* (*contour)(nu,nrpar+5) - (*contour)(nu,nrpar)   * 0.5*(exponent);
      }

    *pbeta_tp = RB.quantile(50,0);
    }

// compute p(beta_0)
  for(nu=0;nu<optionsp->samplesize;nu++)
    {

    exponent = (*contour)(nu,nrpar+4);
    RB(nu,0) = 0.5*(*contour)(nu,nrpar+5) - (*contour)(nu,nrpar)*0.5*(exponent);
    }

  pbeta_0 = RB.quantile(50,0);

//------------------------------------------------------------------------------

  double contourprob=0;

  pbeta_tp = pbeta_t.getV();
  for(t=0;t<optionsp->samplesize;t++,pbeta_tp++)
    {
    if((*pbeta_tp) < pbeta_0)
      contourprob++;
    }
  contourprob = contourprob/optionsp->samplesize;


  optionsp->out("    Bayesian p-value: " +
  ST::doubletostring(contourprob,2) + "\n");

  optionsp->out("\n");

  ST::string path = pathresults.substr(0,pathresults.length()-4)+"_contour.res";

  ofstream out(path.strtochar());
  out << "contourprob" << endl;
  out << ST::doubletostring(contourprob) << endl;

  out.close();

  }



FC_nonp::FC_nonp(const FC_nonp & m)
  : FC(FC(m))
  {

  fsample = m.fsample;

  paramsample = m.paramsample;

  pvalue_sample = m.pvalue_sample;
  pvalue = m.pvalue;

  stype = m.stype;
  likep = m.likep;
  designp = m.designp;
  param = m.param;
  paramlin = m.paramlin;
  parammode = m.parammode;
  paramold = m.paramold;
  paramhelp = m.paramhelp;
  paramKparam = paramKparam;
  betaold = m.betaold;
  betadiff = m.betadiff;
  partres = m.partres;
  lambda=m.lambda;
  tau2 = m.tau2;
  IWLS = m.IWLS;

  Vcenter = m.Vcenter;
  Wcenter = m.Wcenter;
  Ucenter = m.Ucenter;
  Utc = m.Utc;
  ccenter = m.ccenter;

  }


const FC_nonp & FC_nonp::operator=(const FC_nonp & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));

  fsample = m.fsample;

  paramsample = m.paramsample;

  pvalue_sample = m.pvalue_sample;
  pvalue = m.pvalue;

  stype = m.stype;
  likep = m.likep;
  designp = m.designp;
  param = m.param;
  paramlin = m.paramlin;
  paramKparam = paramKparam;
  parammode = m.parammode;
  paramold = m.paramold;
  paramhelp = m.paramhelp;
  betaold = m.betaold;
  betadiff = m.betadiff;
  partres = m.partres;
  lambda=m.lambda;
  tau2 = m.tau2;
  IWLS = m.IWLS;

  Vcenter = m.Vcenter;
  Wcenter = m.Wcenter;
  Ucenter = m.Ucenter;
  Utc = m.Utc;
  ccenter = m.ccenter;

  return *this;
  }


void FC_nonp::get_linparam(void)
  {
  int pos = designp->position_lin;
  int i;
  for (i=0;i<paramlin.rows();i++)
    paramlin(i,0) = designp->FClinearp->beta(pos+i,0);
  }

void FC_nonp::update_IWLS(void)
  {

  // TEST
//  ofstream out("c:\\bayesx\\testh\\results\\beta.res");
//  beta.prettyPrint(out);
  // TEST


  unsigned i;
  double * workparam;

//  lambda = likep->get_scale()/tau2;
  lambda = 1/tau2;

  if (optionsp->nriter == 1)
    {
    paramold.assign(param);
    betaold.assign(beta);
    paramKparam=designp->K.compute_quadform(param,0);
    }


  // Compute log-likelihood with old param, computes workingweight and
  // workingresponse
  double logold = likep->compute_iwls(true,true);
  logold -= 0.5*paramKparam*lambda;

  designp->compute_partres(partres,beta);
  designp->compute_XtransposedWX_XtransposedWres(partres,lambda);

  designp->compute_precision(lambda);

  designp->precision.solve(designp->XWres,paramhelp);

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\paramhelp_v.res");
  // paramhelp.prettyPrint(out);
  // TEST

  workparam = param.getV();
  unsigned nrpar = param.rows();
  for(i=0;i<nrpar;i++,workparam++)
    *workparam = rand_normal();

  designp->precision.solveU(param,paramhelp); // param contains now the proposed
                                              // new parametervector


  if(designp->center)
    centerparam_sample();
//    centerparam();


  paramhelp.minus(param,paramhelp);

  double qold = 0.5*designp->precision.getLogDet()-
                0.5*designp->precision.compute_quadform(paramhelp,0);

  if (designp->position_lin!=-1)
    {
    get_linparam();
    }


  designp->compute_f(param,paramlin,beta,fsample.beta);

  betadiff.minus(beta,betaold);

  designp->update_linpred(betadiff);

  // Compute new log-likelihood

  double lognew = likep->compute_iwls(true,true);
  lognew  -= 0.5*designp->K.compute_quadform(param,0)*lambda;

  designp->compute_partres(partres,beta);
  designp->compute_XtransposedWX_XtransposedWres(partres,lambda);

  designp->compute_precision(lambda);

  designp->precision.solve(designp->XWres,paramhelp);

  // TEST
  // ofstream out2("c:\\bayesx\\testh\\results\\paramhelp_n.res");
  // paramhelp.prettyPrint(out2);
  // TEST


  paramhelp.minus(paramold,paramhelp);
  double qnew = 0.5*designp->precision.getLogDet() -
                0.5*designp->precision.compute_quadform(paramhelp,0);


  double u = log(uniform());
  if (u <= (lognew - logold  + qnew - qold) )
    {
    acceptance++;

    paramKparam=designp->K.compute_quadform(param,0);

    betaold.assign(beta);
    paramold.assign(param);
    }
  else
    {

    betadiff.minus(betaold,beta);
    designp->update_linpred(betadiff);


    param.assign(paramold);
    beta.assign(betaold);
    }


  // TEST
  /*
  ofstream out("c:\\bayesx\\test\\results\\param.res");
  param.prettyPrint(out);

  ofstream out2("c:\\bayesx\\test\\results\\beta.res");
  beta.prettyPrint(out2);
  */
  // TEST

  transform_beta();

  if (designp->position_lin!=-1)
    {
    fsample.update();
    }

 /*
   // w�rde nur f�r pvalue gebraucht!!
  paramsample.beta.assign(param);
  paramsample.transform(0,0) = likep->trmult;
  paramsample.update();
*/

  FC::update();

  }


void FC_nonp::update(void)
  {
  if (IWLS)
    {
    update_IWLS();
    }
  else
    {
    if ((stype == increasing) || (stype==decreasing))
      {
      update_isotonic();
      }
    else
      {
      update_gaussian();
      }
    }
  }


void FC_nonp::update_gaussian(void)
  {

  // TEST
   //ofstream out0("c:\\bayesx\\testh\\results\\beta_re.res");
   //beta.prettyPrint(out0);

   //ofstream out00("c:\\bayesx\\testh\\results\\intvar_re.res");
   //designp->intvar.prettyPrint(out00);
  // TEST


  bool lambdaconst = false;

  betaold.assign(beta);

  double sigmaresp = sqrt(likep->get_scale());
  lambda = likep->get_scale()/tau2;

  // TEST

//  ofstream out("c:\\bayesx\\testh\\results\\responseRE.res");
//  likep->response.prettyPrint(out);

  // TEST

  designp->compute_partres(partres,beta);


  if ((likep->changingweight) || (designp->changingdesign))
    designp->compute_XtransposedWX_XtransposedWres(partres,lambda);
  else
    designp->compute_XtransposedWres(partres, lambda);

  if ((likep->changingweight) || (designp->changingdesign) || (!lambdaconst))
    designp->compute_precision(lambda);


  double * work = paramhelp.getV();
  unsigned i;
  unsigned nrpar = param.rows();
  for(i=0;i<nrpar;i++,work++)
    *work = sigmaresp*rand_normal();

  designp->precision.solveU(paramhelp);

  if (pvalue)
    {
    designp->precision.solve(designp->XWres,parammode);
    param.plus(parammode,paramhelp);
    }
  else
    designp->precision.solve(designp->XWres,paramhelp,param);

  if(designp->center)
//    centerparam();
    centerparam_sample();

  if (designp->position_lin!=-1)
    {
    get_linparam();
    }

  designp->compute_f(param,paramlin,beta,fsample.beta);

  betadiff.minus(beta,betaold);

  designp->update_linpred(betadiff);

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\beta_re.res");
  // beta.prettyPrint(out);
  // TEST


  acceptance++;

  transform_beta();

  if (designp->position_lin!=-1)
    {
    fsample.update();
    }


  if (pvalue)
    update_pvalue();

  paramsample.beta.assign(param);
  paramsample.transform(0,0) = likep->trmult;
  paramsample.update();

  FC::update();

  }


void FC_nonp::update_isotonic(void)
  {

//  if (optionsp->nriter==1)
  {

  // TEST
   // ofstream out0("c:\\bayesx\\testh\\results\\beta_f.res");
   // beta.prettyPrint(out0);

   // ofstream out00("c:\\bayesx\\testh\\results\\intvar_f.res");
   // designp->intvar.prettyPrint(out00);
  // TEST

  unsigned i,j;

  bool lambdaconst = false;

  double sigma2resp = likep->get_scale();
  lambda = likep->get_scale()/tau2;

  betaold.assign(beta);

  designp->compute_partres(partres,beta);


  if ((likep->changingweight) || (designp->changingdesign))
    designp->compute_XtransposedWX_XtransposedWres(partres,lambda);
  else
    designp->compute_XtransposedWres(partres, lambda);

  if ((likep->changingweight) || (designp->changingdesign) || (!lambdaconst))
    designp->compute_precision(lambda);

  int count = 0;
  int maxit = 20;
  double mu;
  double s;


  while(count < maxit)
    {

    for (i=0;i<param.rows();i++)
      {

      mu = 0;
      for (j=0;j<i;j++)    // links
        {
        mu+= param(j,0)*designp->precision(i,j);
        }

      for (j=i+1;j<param.rows();j++)  // rechts
        {
        mu+= param(j,0)*designp->precision(i,j);
        }

      mu = (designp->XWres(i,0) -mu)/designp->precision(i,i);

      s = sqrt(sigma2resp/designp->precision(i,i));

      if(i==0)
        {
        if(stype==increasing)
          param(i,0) = trunc_normal2(-20,param(1,0),mu,s);
        else
          param(i,0) = trunc_normal2(param(1,0),20,mu,s);
        }
      else if(i==param.rows()-1)
        {
        if(stype==increasing)
          param(i,0) = trunc_normal2(param(param.rows()-2,0),20,mu,s);
        else
          param(i,0) = trunc_normal2(-20,param(param.rows()-2,0),mu,s);
        }
      else
        {
        if(stype==increasing)
          {
          param(i,0) = trunc_normal2(param(i-1,0),param(i+1,0),mu,s);
          }
        else
          param(i,0) = trunc_normal2(param(i+1,0),param(i-1,0),mu,s);
        }

      }

    count++;
    }

  /*
  TEST
  ofstream out("c:\\bayesx\\test\\results\\paramhelp.res");
  paramhelp.prettyPrint(out);

  ofstream out2("c:\\bayesx\\test\\results\\param.res");
  param.prettyPrint(out2);
  TEST
  */

  if(designp->center)
    centerparam();
//    centerparam_sample();

  if (designp->position_lin!=-1)
    {
    get_linparam();
    }

  designp->compute_f(param,paramlin,beta,fsample.beta);


  betadiff.minus(beta,betaold);

  designp->update_linpred(betadiff);

  acceptance++;

  transform_beta();

  if (designp->position_lin!=-1)
    {
    fsample.update();
    }

  } //samplesize

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\beta_f.res");
  // beta.prettyPrint(out);
  // TEST


  FC::update();

  }


void FC_nonp::transform_beta(void)
  {
  transform(0,0) = likep->trmult;
  if (designp->position_lin != -1)
  fsample.transform(0,0) = likep->trmult;
  }


bool FC_nonp::posteriormode(void)
  {

  // TEST

  // ofstream out("c:\\bayesx\\test\\results\\data.res");
  // designp->data.prettyPrint(out);

//  ofstream out2("c:\\bayesx\\testh\\results\\intvar.res");
//  designp->intvar.prettyPrint(out2);

  /*
  ofstream out3("c:\\bayesx\\test\\results\\index.res");
  designp->index_data.prettyPrint(out3);
  */
  // TEST

  double h = likep->compute_iwls(true,false);

  betaold.assign(beta);

  lambda = likep->get_scale()/tau2;

  bool lambdaconst = false;

  designp->compute_partres(partres,beta);


  // TEST

  // ofstream out4("c:\\bayesx\\test\\results\\partres.res");
  // partres.prettyPrint(out4);

  // TEST


  if ((likep->changingweight) || (designp->changingdesign))
    designp->compute_XtransposedWX_XtransposedWres(partres, lambda);
  else
    designp->compute_XtransposedWres(partres, lambda);

  if ((likep->changingweight) || (designp->changingdesign) || (!lambdaconst))
    designp->compute_precision(lambda);

  // TEST
//   ofstream out1("c:\\bayesx\\test\\results\\precision.res");
//   designp->precision.print1(out1);
  // TEST


  designp->precision.solve(designp->XWres,param);

  // TEST
  // ofstream out8("c:\\bayesx\\test\\results\\param.res");
  // param.prettyPrint(out8);
  // TEST

  if(designp->center)
    centerparam();
//    centerparam_sample();

  if (designp->position_lin!=-1)
    {
    get_linparam();
    }

  designp->compute_f(param,paramlin,beta,fsample.beta);

  // TEST
//    ofstream out5("c:\\bayesx\\testh\\results\\f.res");
//    beta.prettyPrint(out5);
  // TEST

  betadiff.minus(beta,betaold);

  designp->update_linpred(betadiff);

  transform_beta();

  if (designp->position_lin!=-1)
    {
    fsample.posteriormode_betamean();
    }

  return FC::posteriormode();

  }

void FC_nonp::outoptions(void)
  {
  optionsp->out("  " + title + "\n",true);
  optionsp->out("\n");
  designp->outoptions(optionsp);
  }

void FC_nonp::outresults(const ST::string & pathresults)
  {

  if (pathresults.isvalidfile() != 1)
    {

    FC::outresults(pathresults);
    if (designp->position_lin != -1)
      fsample.outresults(pathresults);

    outresults_acceptance();

    optionsp->out("    Results are stored in file\n");
    optionsp->out("    " +  pathresults + "\n");
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
//    for (i=0;i<designp->datanames.size();i++)
      outres << designp->datanames[designp->datanames.size()-1] << "   ";
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


    if (designp->position_lin!=-1)
      {

      outres << "pmean_d   ";

      if (optionsp->samplesize > 1)
        {
        outres << "pqu"  << l1  << "_d   ";
        outres << "pqu"  << l2  << "_d   ";
        outres << "pmed_d   ";
        outres << "pqu"  << u1  << "_d   ";
        outres << "pqu"  << u2  << "_d   ";
        outres << "pcat" << optionsp->level1 << "_d   ";
        outres << "pcat" << optionsp->level2 << "_d   ";
        }

      }

    outres << endl;


    double * workmean;
    double * workbetaqu_l1_lower_p;
    double * workbetaqu_l2_lower_p;
    double * workbetaqu_l1_upper_p;
    double * workbetaqu_l2_upper_p;
    double * workbetaqu50;

    double * dworkmean;
    double * dworkbetaqu_l1_lower_p;
    double * dworkbetaqu_l2_lower_p;
    double * dworkbetaqu_l1_upper_p;
    double * dworkbetaqu_l2_upper_p;
    double * dworkbetaqu50;


    if (designp->position_lin!=-1)
      {
      workmean = fsample.betamean.getV();
      workbetaqu_l1_lower_p = fsample.betaqu_l1_lower.getV();
      workbetaqu_l2_lower_p = fsample.betaqu_l2_lower.getV();
      workbetaqu_l1_upper_p = fsample.betaqu_l1_upper.getV();
      workbetaqu_l2_upper_p = fsample.betaqu_l2_upper.getV();
      workbetaqu50 = fsample.betaqu50.getV();

      dworkmean = betamean.getV();
      dworkbetaqu_l1_lower_p = betaqu_l1_lower.getV();
      dworkbetaqu_l2_lower_p = betaqu_l2_lower.getV();
      dworkbetaqu_l1_upper_p = betaqu_l1_upper.getV();
      dworkbetaqu_l2_upper_p = betaqu_l2_upper.getV();
      dworkbetaqu50 = betaqu50.getV();

      }
    else
      {
      workmean = betamean.getV();
      workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
      workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
      workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
      workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
      workbetaqu50 = betaqu50.getV();
      }

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


      if (designp->position_lin!=-1)
        {

        outres << *dworkmean << "   ";

        if (optionsp->samplesize > 1)
          {
          outres << *dworkbetaqu_l1_lower_p << "   ";
          outres << *dworkbetaqu_l2_lower_p << "   ";
          outres << *dworkbetaqu50 << "   ";
          outres << *dworkbetaqu_l2_upper_p << "   ";
          outres << *dworkbetaqu_l1_upper_p << "   ";

          if (*dworkbetaqu_l1_lower_p > 0)
            outres << 1 << "   ";
          else if (*dworkbetaqu_l1_upper_p < 0)
            outres << -1 << "   ";
          else
            outres << 0 << "   ";

          if (*dworkbetaqu_l2_lower_p > 0)
            outres << 1 << "   ";
          else if (*dworkbetaqu_l2_upper_p < 0)
            outres << -1 << "   ";
          else
            outres << 0 << "   ";

          }

        if (i <nrpar-1)
          {
          dworkmean++;
          dworkbetaqu_l1_lower_p++;
          dworkbetaqu_l2_lower_p++;
          dworkbetaqu50++;
          dworkbetaqu_l1_upper_p++;
          dworkbetaqu_l2_upper_p++;
          }

        }

      outres << endl;
      }

    if (pvalue)
      compute_pvalue(pathresults);

    }

  }


void FC_nonp::initialize_center(void)
  {
  int nrrest = designp->basisNull.rows();
  int nrpar = param.rows();
  Vcenter = datamatrix(nrpar,nrrest);
  Wcenter = datamatrix(nrrest,nrrest);
  Ucenter = datamatrix(nrrest,nrpar);
  ccenter = datamatrix(nrrest,1);
  Utc = datamatrix(nrpar,1);
  }


void FC_nonp::centerparam_sample(void)
  {

  int nrrest = designp->basisNull.rows();
  int nrpar = param.rows();

  if ((int(Vcenter.rows()) != nrpar) ||
      (Vcenter.cols() != designp->basisNull.rows()))
    initialize_center();

  datamatrix help(nrpar,1);

  int i,j;
//  double sigma2 = likep->get_scale());
  for (i=0;i<nrrest;i++)
    {
    designp->precision.solve(designp->basisNullt[i],help);
    for (j=0;j<int(help.rows());j++)
      Vcenter(j,i) = help(j,0);
    }

  Wcenter.mult(designp->basisNull,Vcenter);
  Ucenter = Wcenter.inverse()*Vcenter.transposed();
  ccenter.mult(designp->basisNull,param);
  Utc = Ucenter.transposed()*ccenter;

//  TEST
//  ofstream out("c:\\bayesx\\test\\results\\Utc.res");
//  Utc.prettyPrint(out);
//  TEST

  param.minus(param,Utc);

//  TEST
//  ofstream out2("c:\\bayesx\\test\\results\\param.res");
//  param.prettyPrint(out2);
//  TEST
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



