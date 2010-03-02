
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

  17     internal_multexp
  */

/*  if (op[14] == "increasing")
    stype = increasing;
  else if (op[14] == "decreasing")
    stype = decreasing;
  else
*/
    stype = unconstrained;

  rtype = additive;
  if (op[12] == "true")
    rtype = mult;

  if (op[17] == "true")
    rtype = multexp;

  if (op[18] == "true")
    pvalue = true;
  else
    pvalue = false;

  }

FC_hrandom::FC_hrandom(MASTER_OBJ * mp,GENERAL_OPTIONS * o,DISTR * lp,DISTR * lp_RE,
                 const ST::string & t,const ST::string & fp,
                 const ST::string & fp2, DESIGN * Dp,
                 vector<ST::string> & op, vector<ST::string> & vn)
     : FC_nonp(mp,o,lp,t,fp,Dp,op,vn)
  {
  read_options(op,vn);
  likep_RE = lp_RE;
  likep_RE->trmult = likep->trmult;
  FCrcoeff = FC(o,"",beta.rows(),beta.cols(),fp2);
  }


FC_hrandom::FC_hrandom(const FC_hrandom & m)
  : FC_nonp(FC_nonp(m))
  {
  rtype = m.rtype;
  likep_RE = m.likep_RE;
  FCrcoeff = m.FCrcoeff;
  response_o = m.response_o;
  linpred_o = m.linpred_o;
  likelihoodc = m.likelihoodc;
  likelihoodn = m.likelihoodn;
  }


const FC_hrandom & FC_hrandom::operator=(const FC_hrandom & m)
  {

  if (this==&m)
	 return *this;
  FC_nonp::operator=(FC_nonp(m));
  rtype = m.rtype;
  likep_RE = m.likep_RE;
  FCrcoeff = m.FCrcoeff;
  response_o = m.response_o;
  linpred_o = m.linpred_o;
  likelihoodc = m.likelihoodc;
  likelihoodn = m.likelihoodn;
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
  if (rtype == mult)
    transform(0,0) = 1.0;
  else if (rtype== multexp)
    transform(0,0) = 1.0;
  else
    FC_nonp::transform_beta();
  }


void FC_hrandom::update_IWLS(void)
  {
  unsigned i;

  lambda = likep->get_scale()/tau2;

  if (optionsp->nriter == 1)
    {
    betaold.assign(beta);
    }

  double * betap = beta.getV();
  double * betaoldp = betaold.getV();

  double * linpredREp;
  if (likep_RE->linpred_current==1)
    linpredREp = likep_RE->linearpred1.getV();
  else
    linpredREp = likep_RE->linearpred2.getV();

  if (likelihoodc.rows() <=1)
    {
    likelihoodc = datamatrix(beta.rows(),1,0);
    likelihoodn = datamatrix(beta.rows(),1,0);
    }

  double postmode;
  double diff;
  double var;
  double u;
  double xwres;


  likep->compute_iwls(true,likelihoodc,designp->ind);

  designp->compute_partres(partres,beta);

  double * workpartres = partres.getV();
  double * worklikelihoodc = likelihoodc.getV();
  double * workWsum = designp->Wsum.getV();

  for (i=0;i<beta.rows();i++,betap++,linpredREp++,
       workpartres++,worklikelihoodc++,workWsum++)

    {

    *worklikelihoodc  -= 0.5*pow((*betap)-(*linpredREp),2)/tau2;

    xwres =  lambda*(*linpredREp)+ (*workpartres);

    var = 1/(*workWsum+lambda);
    postmode =  var * xwres;
    *betap = postmode + sqrt(var)*rand_normal();
    diff = *betap - postmode;
    *worklikelihoodc += -1.0/(2*var)* pow(diff,2)-0.5*log(var);
    }


  betadiff.minus(beta,betaold);
  designp->update_linpred(betadiff);

  likep->compute_iwls(true,likelihoodn,designp->ind);
  designp->compute_partres(partres,beta);

  workpartres = partres.getV();
  double * worklikelihoodn = likelihoodn.getV();
  worklikelihoodc = likelihoodc.getV();
  workWsum = designp->Wsum.getV();

  betap = beta.getV();
  betaoldp = betaold.getV();

  if (likep_RE->linpred_current==1)
    linpredREp = likep_RE->linearpred1.getV();
  else
    linpredREp = likep_RE->linearpred2.getV();

  double * betadiffp = betadiff.getV();

  for (i=0;i<beta.rows();i++,betap++,linpredREp++,betadiffp++,
       betaoldp++,workpartres++,worklikelihoodn++,workWsum++,worklikelihoodc++)

    {

    *worklikelihoodn -= 0.5*pow((*betap)-(*linpredREp),2)/tau2;

    xwres =  lambda*(*linpredREp)+ (*workpartres);


    var = 1/(*workWsum+lambda);
    diff = *betaoldp - var * xwres;

    *worklikelihoodn += -1.0/(2*var)* pow(diff,2)-0.5*log(var);


    nrtrials++;
    u = log(uniform());
    if (u <= (*worklikelihoodn) - (*worklikelihoodc))
      {
      acceptance++;
      *betaoldp = *betap;
      *betadiffp = 0;
      }
    else
      {
      *betadiffp = *betaoldp - *betap;
      *betap = *betaoldp;
      }

    }

  designp->update_linpred(betadiff);

  transform_beta();

  FC::update();

  }



void FC_hrandom::update(void)
  {

  if (IWLS)
    {
    update_IWLS();
    }
  else
    FC_nonp::update();

  set_rcoeff();

  FCrcoeff.acceptance++;
  FCrcoeff.update();

  likep_RE->workingresponse.assign(beta);
  likep_RE->response.assign(beta);
  likep_RE->trmult = likep->trmult;

//  ofstream out("c:\\bayesx\\testh\\results\\response_h.res");
//  beta.prettyPrint(out);

  }


void FC_hrandom::update_response_multexp(void)
  {
/*
  unsigned i,j;

  int size = designp->posbeg.size();

  vector<int>::iterator itbeg = designp->posbeg.begin();
  vector<int>::iterator itend = designp->posend.begin();

  double * * linpredp;

  if (likep->linpred_current==1)
    {
    linpredp = designp->linpredp1.getV();
    linpred_o.assign(likep->linearpred1);
    }
  else
    {
    linpredp = designp->linpredp2.getV();
    linpred_o.assign(likep->linearpred2);
    }


  double ** responsepp = designp->responsep.getV();

  double * workintvar = designp->intvar.getV();

  double * betap = beta.getV();

  for (j=0;j<size;j++,++itbeg,++itend,betap++)
    {
    for (i=*itbeg;i<=*itend;i++,linpredp++,responsepp++,workintvar++)
      {
      *(*responsepp) = *(*responsepp) - *(*linpredp) + exp((*betap)  * (*workintvar));
      *(*linpredp) = *betap * (*workintvar);
      }
    }

  */
  }


void FC_hrandom::update_linpred_multexp(void)
  {

  /*
  unsigned i,j;

  int size = designp->posbeg.size();

  vector<int>::iterator itbeg = designp->posbeg.begin();
  vector<int>::iterator itend = designp->posend.begin();

  double * * linpredp;

  if (likep->linpred_current==1)
    {
    linpredp = designp->linpredp1.getV();
    }
  else
    {
    linpredp = designp->linpredp2.getV();
    }

  double * workintvar = designp->intvar.getV();

  double * betap = beta.getV();
  double * betaoldp = betaold.getV();

  for (j=0;j<size;j++,++itbeg,++itend,betap++,betaoldp++)
    {
    for (i=*itbeg;i<=*itend;i++,linpredp++,workintvar++)
      {
      *(*linpredp) =  linpred_o(designp->index_data(i,0),0) + exp((*betap)  * (*workintvar))
                      - exp(*(betaoldp)  * (*workintvar));

      }
    }
   */
  }


bool FC_hrandom::posteriormode_multexp(void)
  {
/*
  if (response_o.rows()==1)
    {
    response_o = likep->response;
    linpred_o = datamatrix(response_o.rows(),1);
    }
  // intvar = log(f)
  // linpred = etarest + exp(random_effect)* intvar

  likep->optionbool1 = true;
//  likep->changingworkingweights = true;
  likep->updateIWLS = true;

  update_response_multexp();     // linpred = random_effect*intvar

  bool h = posteriormode_additive();

  update_linpred_multexp();

  likep->optionbool1 = false;
//  likep->changingworkingweights = false;
  likep->updateIWLS = false;

  likep->response.assign(response_o);
 */
 
  return true;

  }


bool FC_hrandom::posteriormode_additive(void)
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


bool FC_hrandom::posteriormode(void)
  {
  if (rtype==multexp)
    {
    return posteriormode_multexp();
    }
  else
    {
    return posteriormode_additive();
    }
  }


void FC_hrandom::compute_autocorr_all(const ST::string & path,
                                      unsigned lag, ofstream & outg) const
  {
  FC::compute_autocorr_all(path,lag,outg);
  ST::string path2 = path.substr(0,path.length()-4) + "2.raw";

  FCrcoeff.compute_autocorr_all(path2,lag,outg);
  }


void FC_hrandom::get_samples(const ST::string & filename,ofstream & outg) const
  {
  FC::get_samples(filename,outg);
  ST::string path2 = filename.substr(0,filename.length()-4) + "2.raw";
  FCrcoeff.get_samples(path2,outg);
  }


void FC_hrandom::outgraphs(ofstream & out_stata, ofstream & out_R,
                          const ST::string & path)
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
            << " pmean_tot pqu" << pu1_str << "_tot "
            << " pqu" << po1_str << "_tot" << " pmed_tot pqu" << po2_str
            << "_tot" << " pqu" << pu2_str << "_tot"
            << " pcat" << pu_str << "_tot" << " pcat" << po_str << "_tot"
            << " pmean pqu" << pu1_str
            << " pqu" << po1_str << " pmed pqu" << po2_str << " pqu" <<
            pu2_str << " pcat" << pu_str  << " pcat" << po_str;


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
            << "drop in 1" << endl;

  out_stata << "kdensity pmean_tot" << endl
            << "graph export " << pathps << "_tot.eps, replace"
            << endl << endl;

  out_stata << "kdensity pmean" << endl
            << "graph export " << pathps << ".eps, replace"
            << endl << endl;

  if (computemeaneffect==true)
    {
    out_stata << "kdensity pmean_mu" << endl
              << "graph export " << pathps << "_mu.eps, replace"
              << endl << endl;
    }

  }


void FC_hrandom::outresults(ofstream & out_stata,ofstream & out_R,
                            const ST::string & pathresults)
  {

  if (pathresults.isvalidfile() != 1)
    {

    outgraphs(out_stata,out_R,pathresults);

    FC::outresults(out_stata,out_R,"");
    FCrcoeff.outresults(out_stata,out_R,"");

   if (computemeaneffect==true)
      meaneffect_sample.outresults(out_stata,out_R,"");

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

    double * mu_workmean;
    double * mu_workbetaqu_l1_lower_p;
    double * mu_workbetaqu_l2_lower_p;
    double * mu_workbetaqu_l1_upper_p;
    double * mu_workbetaqu_l2_upper_p;
    double * mu_workbetaqu50;

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
      FC_nonp::compute_pvalue(pathresults);


    }



  }

} // end: namespace MCMC



