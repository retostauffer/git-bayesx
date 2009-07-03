
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
  19      meaneffect
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

  if (op[19] == "true")
    computemeaneffect = true;
  else
    computemeaneffect = false;


  }


FC_nonp::FC_nonp(void)
  {
  }


FC_nonp::FC_nonp(MASTER_OBJ * mp,GENERAL_OPTIONS * o,DISTR * lp,
                 const ST::string & t,const ST::string & fp,
                 DESIGN * Dp,vector<ST::string> & op,
                 vector<ST::string> & vn,bool sstore)
     : FC(o,t,Dp->Zout.rows(),1,fp,sstore)
  {
  read_options(op,vn);
  masterp = mp;
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

  paramsample = FC(o,"",param.rows(),1,fp + ".param",sstore);

  helpcenter = datamatrix(designp->nrpar,1);

  if (pvalue==true)
    {
    pvalue_sample = FC(o,"",param.rows()*2+6,1,fp + ".pvalue",sstore);
    mPhelp = datamatrix(param.rows(),1,0);
    }

  if (computemeaneffect==true)
    {
    meaneffect_sample = FC(o,"",beta.rows(),1,fp+".meaneffect",sstore);
    }

  }


FC_nonp::FC_nonp(const FC_nonp & m)
  : FC(FC(m))
  {

  masterp = m.masterp;

  fsample = m.fsample;

  paramsample = m.paramsample;

  pvalue_sample = m.pvalue_sample;
  pvalue = m.pvalue;
  mPhelp = m.mPhelp;

  computemeaneffect = m.computemeaneffect;
  meaneffect_sample = m.meaneffect_sample;

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
  Vcentert = m.Vcentert;
  Wcenter = m.Wcenter;
  Ucenter = m.Ucenter;
  Utc = m.Utc;
  ccenter = m.ccenter;
  helpcenter = m.helpcenter;

  }


const FC_nonp & FC_nonp::operator=(const FC_nonp & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));

  masterp = m.masterp;

  fsample = m.fsample;

  paramsample = m.paramsample;

  pvalue_sample = m.pvalue_sample;
  pvalue = m.pvalue;
  mPhelp = m.mPhelp;

  computemeaneffect = m.computemeaneffect;
  meaneffect_sample = m.meaneffect_sample;

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
  Vcentert = m.Vcentert;
  Wcenter = m.Wcenter;
  Ucenter = m.Ucenter;
  Utc = m.Utc;
  ccenter = m.ccenter;
  helpcenter = m.helpcenter;

  return *this;
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


    double * mPhelp_p=mPhelp.getV();
    for(i=0;i<nrpar;i++,mPhelp_p++)
      {
      *mPhelp_p = 0;
      parammodep = parammode.getV();
      for(j=0;j<nrpar;j++,parammodep++)
        *mPhelp_p += (*parammodep)*(designp->precision)(j,i);
      }

    contourp++;
    mPhelp_p=mPhelp.getV();
    for(i=0;i<nrpar;i++,contourp++,mPhelp_p++)
      *contourp = *mPhelp_p;

  // TEST
  //   ofstream out("c:\\bayesx\\testh\\results\\beta_pvalue.res");
  //   pvalue_sample.beta.prettyPrint(out);
  // END: TEST

    }

  pvalue_sample.update();

  }


void FC_nonp::compute_pvalue(const ST::string & pathresults)
  {

  if (optionsp->nriter == optionsp->iterations)
    {
    unsigned nrpar = param.rows();

    unsigned k,t,nu;

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

   //---------------------------------------------------------------------------

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

  designp->precision.solve(*(designp->XWres_p),paramhelp);

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
    {
    if (designp->centermethod==meansimple)
      centerparam();
    else
      centerparam_sample();
    }

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

  designp->precision.solve(*(designp->XWres_p),paramhelp);

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


  paramsample.beta.assign(param);
  paramsample.transform(0,0) = likep->trmult;
  paramsample.update();

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


  designp->compute_meaneffect(masterp->level1_likep,meaneffect,beta,
                             meaneffect_sample.beta,computemeaneffect);

  if (computemeaneffect == true)
    {
    meaneffect_sample.update();
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

  if ( (likep->wtype==wweightschange_weightsneqone)  ||
       (likep->wtype==wweightschange_weightsone) ||
       (designp->changingdesign)
     )
    designp->compute_XtransposedWX_XtransposedWres(partres,lambda);
  else
    designp->compute_XtransposedWres(partres, lambda);

  if ((likep->wtype==wweightschange_weightsneqone) ||
      (likep->wtype==wweightschange_weightsone) ||
      (designp->changingdesign) ||
      (!lambdaconst)
      )
    designp->compute_precision(lambda);

  double * work = paramhelp.getV();
  unsigned i;
  unsigned nrpar = param.rows();
  for(i=0;i<nrpar;i++,work++)
    *work = sigmaresp*rand_normal();

  designp->precision.solveU(paramhelp);

  if (pvalue)
    {
    designp->precision.solve(*(designp->XWres_p),parammode);
    param.plus(parammode,paramhelp);
    }
  else
    designp->precision.solve(*(designp->XWres_p),paramhelp,param);

  if(designp->center)
    {
    if (designp->centermethod==meansimple)
      centerparam();
    else
      centerparam_sample();
    }

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
    {
    update_pvalue();
    }

  paramsample.beta.assign(param);
  paramsample.transform(0,0) = likep->trmult;
  paramsample.update();

  FC::update();

  }


void FC_nonp::update_isotonic(void)
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

  if ( (likep->wtype==wweightschange_weightsneqone)  ||
       (likep->wtype==wweightschange_weightsone) ||
       (designp->changingdesign)
     )
    designp->compute_XtransposedWX_XtransposedWres(partres,lambda);
  else
    designp->compute_XtransposedWres(partres, lambda);

  if ((likep->wtype==wweightschange_weightsneqone) ||
      (likep->wtype==wweightschange_weightsone) ||
      (designp->changingdesign) ||
      (!lambdaconst)
      )
    designp->compute_precision(lambda);

  int count = 0;
  int maxit = 20;
  double mu;
  double s;
  double * paramp;
  double * parampi;
  double * XWresp;


  while(count < maxit)
    {

    XWresp = (*(designp->XWres_p)).getV();
    parampi = param.getV();
    for (i=0;i<param.rows();i++,XWresp++,parampi++)
      {

      mu = 0;
      paramp = param.getV();
      for (j=0;j<i;j++,paramp++)    // links
        {
//        mu+= param(j,0)*designp->precision(i,j);
        mu+= (*paramp) * designp->precision(i,j);
        }



      paramp = param.getV()+i+1;
      for (j=i+1;j<param.rows();j++,paramp++)  // rechts
        {
//        mu+= param(j,0)*designp->precision(i,j);
        mu+= (*paramp)*designp->precision(i,j);
        }

//      mu = ((*(designp->XWres_p))(i,0) -mu)/designp->precision(i,i);
      mu = (*XWresp -mu)/designp->precision(i,i);

      s = sqrt(sigma2resp/designp->precision(i,i));

      if(i==0)
        {
        if(stype==increasing)
          *parampi = trunc_normal2(-20,param(1,0),mu,s);
        else
          *parampi = trunc_normal2(param(1,0),20,mu,s);
        }
      else if(i==param.rows()-1)
        {
        if(stype==increasing)
          *parampi = trunc_normal2(param(param.rows()-2,0),20,mu,s);
        else
          *parampi = trunc_normal2(-20,param(param.rows()-2,0),mu,s);
        }
      else
        {
        if(stype==increasing)
          {
          *parampi = trunc_normal2(param(i-1,0),param(i+1,0),mu,s);
          }
        else
          *parampi = trunc_normal2(param(i+1,0),param(i-1,0),mu,s);
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
    {
    if (designp->centermethod==meansimple)
      centerparam();
    else
      centerparam_sample();
    }

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

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\beta_f.res");
  // beta.prettyPrint(out);
  // TEST

  paramsample.beta.assign(param);
  paramsample.transform(0,0) = likep->trmult;
  paramsample.update();

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

  if ( (likep->wtype==wweightschange_weightsneqone)  ||
       (likep->wtype==wweightschange_weightsone) ||
       (designp->changingdesign)
     )
    designp->compute_XtransposedWX_XtransposedWres(partres, lambda);
  else
    designp->compute_XtransposedWres(partres, lambda);

  if ((likep->wtype==wweightschange_weightsneqone) ||
      (likep->wtype==wweightschange_weightsone) ||
      (designp->changingdesign) ||
      (!lambdaconst)
      )
    designp->compute_precision(lambda);

  // TEST
//   ofstream out1("c:\\bayesx\\test\\results\\precision.res");
//   designp->precision.print1(out1);
  // TEST


  designp->precision.solve(*(designp->XWres_p),param);

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
  // ofstream out5("c:\\bayesx\\testh\\results\\f.res");
  // beta.prettyPrint(out5);
  // TEST

  betadiff.minus(beta,betaold);

  designp->update_linpred(betadiff);

  transform_beta();

  if (designp->position_lin!=-1)
    {
    fsample.posteriormode_betamean();
    }

  ST::string n = designp->datanames[0];

  designp->compute_meaneffect(masterp->level1_likep,meaneffect,beta,
                             meaneffect_sample.beta,computemeaneffect);

  return FC::posteriormode();

  }

void FC_nonp::outoptions(void)
  {
  optionsp->out("  " + title + "\n",true);
  optionsp->out("\n");
  designp->outoptions(optionsp);
  }


void FC_nonp::outgraphs(ofstream & out_stata, ofstream & out_R,const ST::string & path)
  {

  ST::string pathps = path.substr(0,path.length()-4) + "_statagraph";

  double u = optionsp->level1;
  double o = optionsp->level2;
  double u1 = optionsp->lower1;
  double u2 = optionsp->upper2;
  double o1 = optionsp->lower2;
  double o2 = optionsp->upper1;
  ST::string u_str = ST::doubletostring(u,0);
  ST::string o_str = ST::doubletostring(o,0);
  ST::string u1_str = ST::doubletostring(u1,5);
  ST::string u2_str = ST::doubletostring(u2,5);
  ST::string o1_str = ST::doubletostring(o1,5);
  ST::string o2_str = ST::doubletostring(o2,5);

  ST::string pu1_str = u1_str.replaceallsigns('.','p');
  ST::string pu2_str = u2_str.replaceallsigns('.','p');
  ST::string po1_str = o1_str.replaceallsigns('.','p');
  ST::string po2_str = o2_str.replaceallsigns('.','p');
  ST::string pu_str = u_str.replaceallsigns('.','p');
  ST::string po_str = o_str.replaceallsigns('.','p');

  ST::string xvar = designp->datanames[designp->datanames.size()-1];


  out_stata << "clear" << endl
            << "infile intnr " << xvar
            << " pmean pqu" << pu1_str
            << " pqu" << po1_str << " pmed pqu" << po2_str << " pqu" << pu2_str
            << " pcat" << pu_str << " pcat" << po_str;


  if (designp->position_lin!=-1)
    {
    out_stata << " pmean_d pqu"
              << pu1_str << "_d"
              << " pqu" << po1_str << "_d"
              << " pmed_d pqu" << po2_str << "_d"
              << " pqu" << pu2_str << "_d"
              << " pcat" << pu_str << "_d"
              << " pcat" << po_str << "_d";
    }


  if (computemeaneffect==true)
    {
    out_stata << " pmean_mu pqu"
              << pu1_str << "_mu"
              << " pqu" << po1_str << "_mu"
              << " pmed_d pqu" << po2_str << "_mu"
              << " pqu" << pu2_str << "_mu";
    }


  out_stata << " using "
            << path << endl
            << "drop in 1" << endl
            << "graph twoway rarea pqu" << pu1_str << " pqu" << pu2_str
            << " " << xvar << ", bcolor(gs13) || rarea pqu" << po1_str
            << " pqu" << po2_str << " " << xvar << " , bcolor(gs10) || /*"
            << endl << " */ scatter pmean "
            << xvar << ", c(l) m(i) clpattern(l) clcolor(gs0) /* "
            << endl << " */ ytitle(\"Effect of "
            << xvar << "\") xtitle(\"" << xvar
            << "\") xlab(,grid) ylab(,grid) legend(off)"
            << endl << "graph export " << pathps << ".eps, replace"
                    << endl << endl;


  if (designp->position_lin!=-1)
    {

    out_stata << "graph twoway rarea pqu" << pu1_str << "_d"
              << " pqu" << pu2_str << "_d"
              << " " << xvar << ", bcolor(gs13) || rarea pqu" << po1_str << "_d"
              << " pqu" << po2_str << "_d"
              << " " << xvar << " , bcolor(gs10) || /*"
              << endl << " */ scatter pmean_d "
              << xvar << ", c(l) m(i) clpattern(l) clcolor(gs0) /* "
              << endl << " */ ytitle(\"Effect of "
              << xvar << "\") xtitle(\"" << xvar
              << "\") xlab(,grid) ylab(,grid) legend(off)"
              << endl << "graph export " << pathps << "_linexluded.eps, replace"
                    << endl << endl;

    }


  if (computemeaneffect==true)
    {

    out_stata << "graph twoway rarea pqu" << pu1_str << "_mu"
              << " pqu" << pu2_str << "_mu"
              << " " << xvar << ", bcolor(gs13) || rarea pqu" << po1_str << "_mu"
              << " pqu" << po2_str << "_mu"
              << " " << xvar << " , bcolor(gs10) || /*"
              << endl << " */ scatter pmean_mu "
              << xvar << ", c(l) m(i) clpattern(l) clcolor(gs0) /* "
              << endl << " */ ytitle(\"Effect of "
              << xvar << "\") xtitle(\"" << xvar
              << "\") xlab(,grid) ylab(,grid) legend(off)"
              << endl << "graph export " << pathps << "_mu.eps, replace"
                    << endl << endl;

    }


  }


void FC_nonp::outresults(ofstream & out_stata, ofstream & out_R,
                        const ST::string & pathresults)
  {

  if (pathresults.isvalidfile() != 1)
    {

    outgraphs(out_stata,out_R,pathresults);

    FC::outresults(out_stata,out_R,pathresults);
    if (designp->position_lin != -1)
      fsample.outresults(out_stata,out_R,pathresults);

   if (computemeaneffect==true)
      meaneffect_sample.outresults(out_stata,out_R,pathresults);

    outresults_acceptance();

    optionsp->out("    Results are stored in file\n");
    optionsp->out("    " +  pathresults + "\n");
    optionsp->out("\n");

    optionsp->out("    Mean effects evaluated at " +
                  designp->datanames[designp->datanames.size()-1] + "=" +
                  designp->effectvalues[designp->meaneffectnr]);

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


    if (computemeaneffect==true)
      {

      outres << "pmean_mu   ";

      if (optionsp->samplesize > 1)
        {
        outres << "pqu"  << l1  << "_mu   ";
        outres << "pqu"  << l2  << "_mu   ";
        outres << "pmed_mu   ";
        outres << "pqu"  << u1  << "_mu   ";
        outres << "pqu"  << u2  << "_mu   ";
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


    double * mu_workmean;
    double * mu_workbetaqu_l1_lower_p;
    double * mu_workbetaqu_l2_lower_p;
    double * mu_workbetaqu_l1_upper_p;
    double * mu_workbetaqu_l2_upper_p;
    double * mu_workbetaqu50;


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


    if (computemeaneffect==true)
      {
      mu_workmean = meaneffect_sample.betamean.getV();
      mu_workbetaqu_l1_lower_p = meaneffect_sample.betaqu_l1_lower.getV();
      mu_workbetaqu_l2_lower_p = meaneffect_sample.betaqu_l2_lower.getV();
      mu_workbetaqu_l1_upper_p = meaneffect_sample.betaqu_l1_upper.getV();
      mu_workbetaqu_l2_upper_p = meaneffect_sample.betaqu_l2_upper.getV();
      mu_workbetaqu50 = meaneffect_sample.betaqu50.getV();
      }


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


      if (computemeaneffect==true)
        {

        outres << *mu_workmean << "   ";

        if (optionsp->samplesize > 1)
          {
          outres << *mu_workbetaqu_l1_lower_p << "   ";
          outres << *mu_workbetaqu_l2_lower_p << "   ";
          outres << *mu_workbetaqu50 << "   ";
          outres << *mu_workbetaqu_l2_upper_p << "   ";
          outres << *mu_workbetaqu_l1_upper_p << "   ";
          }

        if (i <nrpar-1)
          {
          mu_workmean++;
          mu_workbetaqu_l1_lower_p++;
          mu_workbetaqu_l2_lower_p++;
          mu_workbetaqu50++;
          mu_workbetaqu_l1_upper_p++;
          mu_workbetaqu_l2_upper_p++;
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
  Vcentert = datamatrix(nrrest,nrpar);
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

  int i,j;
  double * helpcenterp = helpcenter.getV();

  double * Vcentertp = Vcentert.getV();
  double * Vcenterp;
  for (i=0;i<nrrest;i++)
    {
    Vcenterp = Vcenter.getV()+i;

    designp->precision.solve(designp->basisNullt[i],helpcenter);

    helpcenterp = helpcenter.getV();

    for (j=0;j<nrpar;j++,helpcenterp++,Vcentertp++,Vcenterp+=nrrest)
      {
      *Vcenterp = *helpcenterp;
      *Vcentertp = *helpcenterp;
      }
    }

  // TEST
   /*
   ofstream out0("c:\\bayesx\\testh\\results\\basisnull.res");
   (designp->basisNullt[0]).prettyPrint(out0);

   ofstream out("c:\\bayesx\\testh\\results\\Vcenter.res");
   Vcenter.prettyPrint(out);

   ofstream out2("c:\\bayesx\\testh\\results\\praecision.res");
   designp->precision.print4(out2);
   */
  // TEST

  Wcenter.mult(designp->basisNull,Vcenter);
  Ucenter = Wcenter.inverse()*Vcentert;
  ccenter.mult(designp->basisNull,param);
  Utc = Ucenter.transposed()*ccenter;

  //  TEST

  // ofstream out4("c:\\bayesx\\testh\\results\\param.res");
  // param.prettyPrint(out4);

  //  TEST

  param.minus(param,Utc);

  //  TEST
  // ofstream out5("c:\\bayesx\\testh\\results\\paramneu.res");
  // param.prettyPrint(out5);

  // ofstream out6("c:\\bayesx\\testh\\results\\Utc.res");
  // Utc.prettyPrint(out6);
  //  TEST


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


void FC_nonp::compute_autocorr_all(const ST::string & path,
                                      unsigned lag, ofstream & outg) const
  {
  paramsample.compute_autocorr_all(path,lag,outg);
  }


void FC_nonp::get_samples(const ST::string & filename,ofstream & outg) const
  {
  paramsample.get_samples(filename,outg);
  }



void FC_nonp::reset(void)
  {

  }


} // end: namespace MCMC



