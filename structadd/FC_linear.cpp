
#include "FC_linear.h"


//------------------------------------------------------------------------------
//----------------- CLASS: FC implementation of member functions ---------------
//------------------------------------------------------------------------------


namespace MCMC
{


void FC_linear::read_options(vector<ST::string> & op,vector<ST::string> & vn)
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


FC_linear::FC_linear(void)
  {
  }


int FC_linear::add_variable(const datamatrix & d,ST::string & name)
  {
  datanames.push_back(name);
  designhelp.push_back(d);
  return designhelp.size()-1;
  }


FC_linear::FC_linear(MASTER_OBJ * mp,GENERAL_OPTIONS * o,DISTR * lp,
                    datamatrix & d,
                 vector<ST::string> & vn, const ST::string & t,
                 const ST::string & fp,bool sstore,bool cent)
     : FC(o,t,1,1,fp,sstore)
  {

  masterp = mp;
  likep = lp;
  int i;
  datanames = vn;
  if (datanames.size() > 0)
    {
    for (i=0;i<d.cols();i++)
      designhelp.push_back(d.getCol(i));
    }
  initialize = false;
  IWLS = likep->updateIWLS;

  center = cent;

  }


FC_linear::FC_linear(const FC_linear & m)
  : FC(FC(m))
  {
  constposition = m.constposition;
  masterp = m.masterp;
  IWLS = m.IWLS;
  likep = m.likep;
  design = m.design;
  designhelp = m.designhelp;
  meaneffectdesign = m.meaneffectdesign;
  XWX = m.XWX;
  XWXold = m.XWXold;
  XWXroot = m.XWXroot;
  Xt = m.Xt;
  initialize = m.initialize;
  residual = m.residual;
  Xtresidual = m.Xtresidual;
  betaold = m.betaold;
  betadiff = m.betadiff;
  betam = m.betam;
  mode = m.mode;
  help = m.help;
  linold = m.linold;
  linnew = m.linnew;
  linmode = m.linmode;
  proposal=m.proposal;
  diff = m.diff;
  linoldp = m.linoldp;
  linnewp = m.linnewp;
  datanames = m.datanames;
  mean_designcols = m.mean_designcols;
  center = m.center;
  }


const FC_linear & FC_linear::operator=(const FC_linear & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  constposition = m.constposition;
  masterp = m.masterp;
  IWLS = m.IWLS;
  likep = m.likep;
  design = m.design;
  designhelp = m.designhelp;
  meaneffectdesign = m.meaneffectdesign;
  XWX = m.XWX;
  XWXold = m.XWXold;
  XWXroot = m.XWXroot;
  Xt = m.Xt;
  initialize = m.initialize;
  residual = m.residual;
  Xtresidual = m.Xtresidual;
  betaold = m.betaold;
  betadiff = m.betadiff;
  betam = m.betam;
  mode = m.mode;
  proposal=m.proposal;
  help = m.help;
  linold = m.linold;
  linnew = m.linnew;
  linmode = m.linmode;
  diff = m.diff;
  linoldp = m.linoldp;
  linnewp = m.linnewp;
  mean_designcols = m.mean_designcols;
  datanames = m.datanames;
  center = m.center;
  return *this;
  }


void FC_linear::add_linpred(datamatrix & l)
  {

  if (likep->linpred_current==1)
    likep->linearpred1.plus(l);
  else
    likep->linearpred2.plus(l);
  }




void FC_linear::update_IWLS(void)
  {

  double qoldbeta;
  double qnewbeta;

  if (!initialize)
    create_matrices();


  if (design.cols() == 1)
    {
    double moold;

    if (optionsp->nriter == 1)
      {
      mode(0,0) = beta(0,0);
      }

    double logold = likep->loglikelihood(true);

    diff.mult_scalar(design,mode(0,0)-beta(0,0));
    add_linpred(diff);

    double h = likep->compute_iwls(true,false);

    compute_XWX(XWXold);

    Xtresidual(0,0) = compute_XtWpartres(mode(0,0));

    moold = mode(0,0);
    mode(0,0) = Xtresidual(0,0)/XWXold(0,0);

    qoldbeta = -0.5*pow(beta(0,0)-mode(0,0),2)*XWXold(0,0);

    proposal(0,0) = mode(0,0) + sqrt(1/XWXold(0,0)) *rand_normal();

    qnewbeta = -0.5* pow(proposal(0,0)-mode(0,0),2)*XWXold(0,0);

    diff.mult_scalar(design,proposal(0,0)-moold);
    add_linpred(diff);                           // (mit proposed)

    double logprop = likep->loglikelihood();     // mit proposed

    double u = log(uniform());

    if (u <= (logprop + qoldbeta - logold - qnewbeta) )
      {
      beta.assign(proposal);
      acceptance++;
      }
    else
      {
      diff.mult_scalar(design,beta(0,0)-proposal(0,0));
      add_linpred(diff);
      }

    FC::update();

    }
  else

    {

    if (optionsp->nriter == 1)
      {
      linold.mult(design,beta);
      mode.assign(beta);
      }

    double logold = likep->loglikelihood(true);

    linmode.mult(design,mode);
    diff.minus(linmode,*linoldp);
    add_linpred(diff);

    double h = likep->compute_iwls(true,false);

    compute_XWX(XWXold);

    XWXroot = XWXold.root();

    compute_Wpartres(linmode);
    Xtresidual.mult(Xt,residual);

    XWXroot.solveroot(Xtresidual,help,mode);

    help.minus(beta,mode);
    qoldbeta = -0.5*XWXold.compute_quadform(help);


    unsigned i;
    double * workh = help.getV();
    for(i=0;i<help.rows();i++,workh++)
      *workh = rand_normal();

    XWXroot.solveroot_t(help,proposal);

    proposal.plus(mode);


    help.minus(proposal,mode);

    qnewbeta = -0.5*XWXold.compute_quadform(help);


    linnewp->mult(design,proposal);

    diff.minus(*linnewp,linmode);

    add_linpred(diff);                           // (mit proposed)

    double logprop = likep->loglikelihood();     // mit proposed


    double u = log(uniform());

    if (u <= (logprop + qoldbeta - logold - qnewbeta) )
      {
      datamatrix * mp = linoldp;
      linoldp = linnewp;
      linnewp = mp;

      beta.assign(proposal);

      acceptance++;
      }
    else
      {
      diff.minus(*linoldp,*linnewp);
      add_linpred(diff);
      }

    FC::update();
    }

  }

void FC_linear::update(void)
  {
  if (datanames.size() > 0)
    {
    if (IWLS)
      update_IWLS();
    else
      update_gaussian();

    masterp->level1_likep->meaneffect -= meaneffect;
    meaneffect = (meaneffectdesign*beta)(0,0);
    masterp->level1_likep->meaneffect += meaneffect;

    }
  else
    nosamples = true;
  }


void FC_linear::update_gaussian(void)
  {
  if (datanames.size() > 0)
    {
    if (!initialize)
      create_matrices();

    if ((likep->wtype==wweightschange_weightsneqone) ||
        (likep->wtype==wweightschange_weightsone) ||
        (optionsp->nriter==1))
      {
      compute_XWX(XWX);
      XWXroot = XWX.root();
      }

    linold.mult(design,beta);
    compute_Wpartres(linold);
    Xtresidual.mult(Xt,residual);

    XWXroot.solveroot(Xtresidual,help,betam);

    double sigmaresp = sqrt(likep->get_scale());
    unsigned i;
    double * workh = help.getV();
    for(i=0;i<help.rows();i++,workh++)
      *workh = sigmaresp*rand_normal();

    XWXroot.solveroot_t(help,beta);
    beta.plus(betam);

    betadiff.minus(beta,betaold);

    if (likep->linpred_current==1)
      likep->linearpred1.addmult(design,betadiff);
    else
      likep->linearpred2.addmult(design,betadiff);

    betaold.assign(beta);

    transform(0,0) = likep->trmult;
    acceptance++;

    FC::update();
    }
  }




void FC_linear::compute_XWX(datamatrix & r)
  {
  unsigned i,j,k;
  unsigned nrconst = beta.rows();
  unsigned nrobs = Xt.cols();
  double * Xt_ip;
  double * Xt_jp;
  double * workingweightp;
  double help;

  if (likep->wtype==wweightsnochange_one)
    {
    for (i=0;i<nrconst;i++)
      for (j=i;j<nrconst;j++)
        {
        help = 0;
        Xt_ip = Xt.getV()+i*nrobs;
        Xt_jp = Xt.getV()+j*nrobs;

        for (k=0;k<nrobs;k++,Xt_ip++,Xt_jp++)
          help += (*Xt_ip)*(*Xt_jp);

        r(i,j) = help;
        if (i!=j)
          r(j,i) = help;
        }
    }
  else
    {
    for (i=0;i<nrconst;i++)
      for (j=i;j<nrconst;j++)
        {
        help = 0;
        Xt_ip = Xt.getV()+i*nrobs;
        Xt_jp = Xt.getV()+j*nrobs;
        workingweightp = likep->workingweight.getV();

        for (k=0;k<nrobs;k++,Xt_ip++,Xt_jp++,workingweightp++)
          help += (*workingweightp) * (*Xt_ip)*(*Xt_jp);

        r(i,j) = help;
        if (i!=j)
          r(j,i) = help;
        }
    }

  }


void FC_linear::compute_meaneffect_design(void)
  {
  unsigned i,j;

  meaneffectdesign = datamatrix(1,design.cols(),0);

  double  mhelp;

  double bestdiff;
  double currentdiff;

  for (j=0;j<design.cols();j++)
    {
    mhelp = design.mean(j);
    bestdiff = fabs(design(0,j) - mhelp);
    meaneffectdesign(0,j) = design(0,j);
    for (i=1;i<design.rows();i++)
      {
      currentdiff = fabs(design(i,j) - mhelp);
      if (currentdiff < bestdiff)
        {
        bestdiff = currentdiff;
        meaneffectdesign(0,j) = design(i,j);
        }

      }

    }

  }


void FC_linear::find_const(datamatrix & design)
  {
  constposition = -1;
  bool constfound = false;
  int i=0;
  int j;
  while (constfound==false && i < design.cols())
    {
    j = 0;
    bool allone = true;
    while (allone == true && j < design.rows())
      {
      if (design(j,i) != 1)
        allone = false;
      j++;
      }
    if (allone == true)
      {
      constfound = true;
      constposition = i;
      }
    i++;
    }

  }


void FC_linear::create_matrices(void)
  {

  int i,j;
  design = datamatrix(designhelp[0].rows(),designhelp.size());
  for(i=0;i<designhelp.size();i++)
    design.putCol(i,designhelp[i]);

  find_const(design);

  if (center == true)
    {
    double m;
    mean_designcols = datamatrix(1,design.cols(),1);
    for (i=0;i<design.cols();i++)
      {
      if (i!= constposition)
        {
        m = design.mean(i);
        for (j=0;j<design.rows();j++)
          design(j,i) -= m;
        mean_designcols(0,i) = m;
        }
      }

    }

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\design.res");
  // design.prettyPrint(out);
  // TEST

  compute_meaneffect_design();

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\meandesign.res");
  // meaneffectdesign.prettyPrint(out);
  // TEST

  Xt = design.transposed();
  XWX = datamatrix(design.cols(),design.cols(),0);


  residual = datamatrix(design.rows(),1,0);
  Xtresidual = datamatrix(design.cols(),1,0);

  setbeta(design.cols(),1,0);
  betaold=datamatrix(beta.rows(),1,0);
  betadiff = betaold;
  betam = beta;
  help = beta;
  linold = datamatrix(design.rows(),1,0);
  initialize=true;

  // For IWLS
  linnew = datamatrix(design.rows(),1,0);
  linmode = datamatrix(design.rows(),1,0);
  diff = datamatrix(design.rows(),1,0);
  linnewp = &linnew;
  linoldp = &linold;
  mode = beta;
  proposal = beta;
  XWXold = datamatrix(design.cols(),design.cols(),0);
  }


void FC_linear::compute_Wpartres(datamatrix & linpred)
  {
  unsigned i;
  double * workingweightp = likep->workingweight.getV();
  double * workingresponsep = likep->workingresponse.getV();
  double * residualp = residual.getV();
  double * linpredp = linpred.getV();

  double * predictorp;
  if (likep->linpred_current==1)
    predictorp = likep->linearpred1.getV();
  else
    predictorp = likep->linearpred2.getV();

  if (likep->wtype == wweightsnochange_one)
    {
    for (i=0;i<likep->nrobs;i++,workingresponsep++,
                                residualp++,linpredp++,predictorp++)
      *residualp = ((*workingresponsep)  - (*predictorp)
                   + (*linpredp));
    }
  else
    {
    for (i=0;i<likep->nrobs;i++,workingweightp++,workingresponsep++,
                                residualp++,linpredp++,predictorp++)
      *residualp = *workingweightp * ((*workingresponsep)  - (*predictorp)
                   + (*linpredp));
    }
  }


double FC_linear::compute_XtWpartres(double & mo)
  {

  unsigned i;
  double * workingweightp = likep->workingweight.getV();
  double * workingresponsep = likep->workingresponse.getV();
  double res=0;


  double * predictorp;
  if (likep->linpred_current==1)
    predictorp = likep->linearpred1.getV();
  else
    predictorp = likep->linearpred2.getV();

  double * workdesign = design.getV();

  for (i=0;i<likep->nrobs;i++,workingweightp++,workingresponsep++,
                              predictorp++,workdesign++)
    res += *workdesign * (*workingweightp) *
           ((*workingresponsep)  - (*predictorp) + *workdesign*mo);

  return res;
  }



bool FC_linear::posteriormode(void)
  {

  if (datanames.size() > 0)
    {
    if (!initialize)
      create_matrices();

    double h = likep->compute_iwls(true,false);

    compute_XWX(XWX);

    linold.mult(design,beta);
    compute_Wpartres(linold);

    Xtresidual.mult(Xt,residual);

    beta = XWX.solve(Xtresidual);

    betadiff.minus(beta,betaold);

    if (likep->linpred_current==1)
      likep->linearpred1.addmult(design,betadiff);
    else
      likep->linearpred2.addmult(design,betadiff);

    betaold.assign(beta);

    transform(0,0) = likep->trmult;

    masterp->level1_likep->meaneffect -= meaneffect;
    meaneffect = (meaneffectdesign*beta)(0,0);
    masterp->level1_likep->meaneffect += meaneffect;

    return FC::posteriormode();
    }

  return true;
  }


void FC_linear::compute_autocorr_all(const ST::string & path,
                              unsigned lag, ofstream & outg) const
  {
  if (datanames.size() > 0)
    {
    FC::compute_autocorr_all(path,lag,outg);
    }
  }



void FC_linear::outoptions(void)
  {
//  optionsp->out("  OPTIONS FOR TERM: " + title + "\n",true);
//  optionsp->out("\n");
  }

void FC_linear::outresults(ofstream & out_stata,ofstream & out_R,
                           const ST::string & pathresults)
  {
  if (datanames.size() > 0)
    {
    unsigned i;
    unsigned nrconst = design.cols();

    FC::outresults(out_stata,out_R,pathresults);

    ST::string l1 = ST::doubletostring(optionsp->lower1,4);
    ST::string l2 = ST::doubletostring(optionsp->lower2,4);
    ST::string u1 = ST::doubletostring(optionsp->upper1,4);
    ST::string u2 = ST::doubletostring(optionsp->upper2,4);
    l1 = l1.replaceallsigns('.','p');
    l2 = l2.replaceallsigns('.','p');
    u1 = u1.replaceallsigns('.','p');
    u2 = u2.replaceallsigns('.','p');

    ofstream outp(pathresults.strtochar());

    if (pathresults.isvalidfile() != 1)
      outp << "paramnr varname pmean pstd pqu" << l1 << " pqu" << l2 <<
              " pmed pqu" << u1 << " pqu" << u2 << " pcat" << optionsp->level1
               << " pcat" << optionsp->level2 << endl;

    optionsp->out("\n");

    ST::string l;
    int maxvarnamelength = 0;
    int len;

    for(i=0;i<nrconst;i++)
      {
      len = datanames[i].length();
      if (len > maxvarnamelength)
        maxvarnamelength = len;
      }

    if (maxvarnamelength>10)
      l = ST::string(' ',maxvarnamelength-4);
    else
      l = "  ";

      ST::string help =  ST::doubletostring(optionsp->lower1,4) + "% quant.";
      ST::string levell = help + ST::string(' ',15-help.length());
      help = ST::doubletostring(optionsp->upper2,4) + "% quant.";
      ST::string levelu = help + ST::string(' ',15-help.length());

      optionsp->out("    Variable" + l +
                    "mean           " +
                    "Std. Dev.      " +
                    levell +
                    "median         " +
                    levelu + "\n");

      ST::string mean;
      ST::string std;
      ST::string qu10;
      ST::string qu50;
      ST::string qu90;

      double m,stddouble;

      unsigned nsp;

      for (i=0;i<nrconst;i++)
        {

        if (maxvarnamelength  > 10)
          nsp = 4+maxvarnamelength-datanames[i].length();
        else
          nsp = 10-datanames[i].length();

        m= betamean(i,0);

        if (betavar(i,0) == 0)
            stddouble = 0;
        else
          stddouble = sqrt(betavar(i,0));

        if (pathresults.isvalidfile() != 1)
          {
          outp << (i+1) << "   ";
          outp << datanames[i] << "   ";
          outp << m << "   ";
          outp << stddouble << "   ";



          outp << betaqu_l1_lower(i,0) << "   ";
          outp << betaqu_l2_lower(i,0) << "   ";
          outp << betaqu50(i,0) << "   ";
          outp << betaqu_l2_upper(i,0) << "   ";
          outp << betaqu_l1_upper(i,0) << "   ";
          if (betaqu_l1_lower(i,0) > 0)
            outp << "1   ";
          else if (betaqu_l1_upper(i,0) < 0)
            outp << "-1   ";
          else
            outp << "0   ";

          if (betaqu_l2_lower(i,0) > 0)
            outp << "1   ";
          else if (betaqu_l2_upper(i,0) < 0)
            outp << "-1   ";
          else
            outp << "0   ";

          outp << endl;


          optionsp->out(ST::outresults(nsp,datanames[i],m,
                        stddouble,betaqu_l1_lower(i,0),
                        betaqu50(i,0),betaqu_l1_upper(i,0)) + "\n");

          }

        }

      optionsp->out("\n");

      optionsp->out("    Results for fixed effects are also stored in file\n");
      optionsp->out("    " + pathresults + "\n");

      if (center==true)
        {
        optionsp->out("\n");
        optionsp->out("    Note: Covariates with linear effects are centered around zero before estimation\n");
        optionsp->out("          Centering of covariates may improve the mixing of the MCMC sampler while\n");
        optionsp->out("          the regression coefficents are unchanged\n");
        optionsp->out("          However the intercept is changed due to the centering of covariates.\n");
        optionsp->out("          The means of the covariates are:\n");
        unsigned k;
        for (k=0;k<mean_designcols.cols();k++)
          {
          if (k != constposition)
            {
            optionsp->out("          " + datanames[k] + ": " + ST::doubletostring(mean_designcols(0,k),6) + "\n");  

            }

          }

        }

      optionsp->out("\n");
    }

  }


void FC_linear::reset(void)
  {

  }


} // end: namespace MCMC



