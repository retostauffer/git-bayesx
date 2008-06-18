
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


FC_linear::FC_linear(GENERAL_OPTIONS * o,DISTR * lp,datamatrix & d,
                 vector<ST::string> & vn, const ST::string & t,
                 const ST::string & fp)
     : FC(o,t,1,1,fp)
  {

  likep = lp;
  design = d;
  datanames = vn;
  initialize = false;
  IWLS = likep->updateIWLS;
  }


FC_linear::FC_linear(const FC_linear & m)
  : FC(FC(m))
  {
  IWLS = m.IWLS;
  likep = m.likep;
  design = m.design;
  XWX = m.XWX;
  XWXroot = m.XWXroot;
  Xt = m.Xt;
  initialize = m.initialize;
  residual = m.residual;
  Xtresidual = m.Xtresidual;
  betaold = m.betaold;
  betadiff = m.betadiff;
  linold = m.linold;
  datanames = m.datanames;
  }


const FC_linear & FC_linear::operator=(const FC_linear & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  IWLS = m.IWLS;
  likep = m.likep;
  design = m.design;
  XWX = m.XWX;
  XWXroot = m.XWXroot;  
  Xt = m.Xt;
  initialize = m.initialize;
  residual = m.residual;
  Xtresidual = m.Xtresidual;
  betaold = m.betaold;
  betadiff = m.betadiff;
  linold = m.linold;
  datanames = m.datanames;  
  return *this;
  }


void FC_linear::update_IWLS(void)
  {
  FC::update();
  }

void FC_linear::update(void)
  {
  if (IWLS)
    update_IWLS();
  else
    update_gaussian();
  }


void FC_linear::update_gaussian(void)
  {

  if (!initialize)
    create_matrices();

  if (likep->changingweight || optionsp->nriter==1)
    {
    compute_XWX();
    XWXroot = XWX.root();
    ofstream out("c:\\bayesx\\test\\results\\XWXW.root.res");
    XWXroot.prettyPrint(out);

    }

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

  FC::update();
  }


void FC_linear::compute_XWX(void)
  {
  unsigned i,j,k;
  unsigned nrconst = beta.rows();
  unsigned nrobs = Xt.cols();
  double * Xt_ip;
  double * Xt_jp;
  double * workingweightp;
  double help;
  for (i=0;i<nrconst;i++)
    for (j=i;j<nrconst;j++)
      {
      help = 0;
      Xt_ip = Xt.getV()+i*nrobs;
      Xt_jp = Xt.getV()+j*nrobs;
      workingweightp = likep->workingweight.getV();

      for (k=0;k<nrobs;k++,Xt_ip++,Xt_jp++,workingweightp++)
        help += (*workingweightp) * (*Xt_ip)*(*Xt_jp);

      XWX(i,j) = help;
      if (i!=j)
        XWX(j,i) = help;
      }

  }


void FC_linear::create_matrices(void)
  {
  Xt = design.transposed();
  XWX = datamatrix(design.cols(),design.cols(),0);
  

  residual = datamatrix(design.rows(),1,0);
  Xtresidual = datamatrix(design.cols(),1,0);

  setbeta(design.cols(),1,0);
  betaold=datamatrix(beta.rows(),1,0);
  betadiff = betaold;
  linold = datamatrix(design.rows(),1,0);
  initialize=true;
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

  for (i=0;i<likep->nrobs;i++,workingweightp++,workingresponsep++,
                              residualp++,linpredp++,predictorp++)
    *residualp = *workingweightp * ((*workingresponsep)  - (*predictorp)
                 + (*linpredp));
  }


bool FC_linear::posteriormode(void)
  {

  if (!initialize)
    create_matrices();

  compute_XWX();

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

  return FC::posteriormode();

  }


void FC_linear::outoptions(void)
  {
//  optionsp->out("  OPTIONS FOR TERM: " + title + "\n",true);
//  optionsp->out("\n");
  }

void FC_linear::outresults(const ST::string & pathresults)
  {

  unsigned i;
  unsigned nrconst = design.cols();

  FC::outresults(pathresults);

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
    l = ST::string(' ',maxvarnamelength-6);
  else
    l = "  ";

    ST::string help =  ST::doubletostring(optionsp->lower1,4) + "% quant.";
    ST::string levell = help + ST::string(' ',15-help.length());
    help = ST::doubletostring(optionsp->upper2,4) + "% quant.";
    ST::string levelu = help + ST::string(' ',15-help.length());

    optionsp->out("  Variable" + l +
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
        nsp = 2+maxvarnamelength-datanames[i].length();
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
        }

      optionsp->out(ST::outresults(nsp,datanames[i],m,
                      stddouble,betaqu_l1_lower(i,0),
                      betaqu50(i,0),betaqu_l1_upper(i,0)) + "\n");

      }

    optionsp->out("\n");

    optionsp->out("  Results for fixed effects are also stored in file\n");
    optionsp->out("  " + pathresults + "\n");

    optionsp->out("\n");

  }


void FC_linear::reset(void)
  {

  }


} // end: namespace MCMC



