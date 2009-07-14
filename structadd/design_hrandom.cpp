
#include "design_hrandom.h"
#include "clstring.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//------------ CLASS: DESIGN_hrandom implementation of member functions --------
//------------------------------------------------------------------------------


void DESIGN_hrandom::read_options(vector<ST::string> & op,
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
  */


  datanames = vn;

  }


DESIGN_hrandom::DESIGN_hrandom(void)
  {

  }


DESIGN_hrandom::DESIGN_hrandom(const datamatrix & dm, const datamatrix & iv,
                               DISTR * dp, FC_linear * fcl, DISTR * dp_RE,
                               vector<ST::string> & op,vector<ST::string> & vn)
      : DESIGN(dp,fcl)
  {

  read_options(op,vn);

  likep_RE = dp_RE;

  type = Hrandom;

  init_data(dm,iv);

  compute_penalty();

  datamatrix  help(Zout.rows(),1,1);

  compute_XtransposedWX();
  compute_XtransposedWres(help,1);
  Wsum = datamatrix(posbeg.size(),1,1);
  compute_precision(1.0);

  center = false;

  }


DESIGN_hrandom::DESIGN_hrandom(const DESIGN_hrandom & m)
  : DESIGN(DESIGN(m))
  {
  likep_RE = m.likep_RE;
  }


const DESIGN_hrandom & DESIGN_hrandom::operator=(const DESIGN_hrandom & m)
  {
  if (this == &m)
    return *this;
  DESIGN::operator=(DESIGN(m));
  likep_RE = m.likep_RE;
  return *this;

  }


void DESIGN_hrandom::init_data(const datamatrix & dm, const datamatrix & iv)
  {

  // TASK: sorts the data such that the precision has minimum envelope
  //       computes index_data
  //       computes Zout, posbeg, posend
  //       computes nrpar
  //       computes effectvalues
  //       initializes datanames

  // TASK of make_index: sorts the data,
  //                     creates sorted intvar, data2
  //                     initializes index_data,
  //                     posbeg, posend, effectvalues

  make_index(dm,iv);

  meaneffectnr = compute_modecategorie();
  compute_meaneffectintvar();

  nrpar = posbeg.size();

  Zout = datamatrix(posbeg.size(),1,1);
  index_Zout = statmatrix<int>(Zout.rows(),1);
  index_Zout.indexinit();

  consecutive = 1;

  compute_Zout_transposed();



  }



void DESIGN_hrandom::compute_penalty(void)
  {
  K =   envmatrix<double>(1,nrpar);
  rankK = nrpar;
  }



void DESIGN_hrandom::compute_XtransposedWX(void)
  {

  if (XWXdeclared == false)
    {
    XWX = envmatdouble(0,nrpar);
    XWXdeclared = true;
    }


  unsigned i;

  double * Wsump = Wsum.getV();
  vector<double>::iterator d = XWX.getDiagIterator();

  for (i=0;i<Wsum.rows();i++,++d,Wsump++)
    *Wsump = *d;

 }



void DESIGN_hrandom::compute_XtransposedWres(datamatrix & partres, double l)
  {

  if (XWresdeclared == false)
    {
    XWres = datamatrix(nrpar,1);
    XWresdeclared = true;
    }


  double * workXWres = XWres.getV();

  double * linpredREp;
  if (likep_RE->linpred_current==1)
    linpredREp = likep_RE->linearpred1.getV();
  else
    linpredREp = likep_RE->linearpred2.getV();

  double * partresp = partres.getV();

  unsigned i;

  for(i=0;i<nrpar;i++,workXWres++,linpredREp++,partresp++)
    *workXWres =  l*(*linpredREp)+(*partresp);

  XWres_p = &XWres;

  }


void DESIGN_hrandom::compute_precision(double l)
  {

  if (precisiondeclared==false)
    {
    precision = envmatdouble(K.getXenv(),0,nrpar);
    precisiondeclared = true;
    }

  precision.addtodiag(XWX,K,1.0,l);

  /*
  // TEST
  ofstream out2("c:\\bayesx\\test\\results\\XWX.res");
  XWX.print2(out2);

  ofstream out3("c:\\bayesx\\test\\results\\K.res");
  K.print2(out3);


  ofstream out("c:\\bayesx\\test\\results\\precision.res");
  precision.print2(out);
  // TEST
  */

  }



void DESIGN_hrandom::compute_meaneffect(DISTR * level1_likep,double & meaneffect,
                                datamatrix & beta,datamatrix & meaneffectbeta,
                                bool computemeaneffect)

  {

  level1_likep->meaneffect -= meaneffect;

  double * linpredREp;
  double linm;
  if (likep_RE->linpred_current==1)
    {
    linpredREp = likep_RE->linearpred1.getV();
    linm = likep_RE->linearpred1(meaneffectnr,0);
    }
  else
    {
    linpredREp = likep_RE->linearpred2.getV();
    linm = likep_RE->linearpred2(meaneffectnr,0);
    }

  meaneffect = meaneffectintvar*(beta(meaneffectnr,0) - linm);

  if (computemeaneffect==true)
    {
    unsigned i;
    double * betap = beta.getV();
    double * meffectp = meaneffectbeta.getV();
    double l;
    for(i=0;i<beta.rows();i++,meffectp++,betap++,linpredREp++)
      {
      l=level1_likep->meaneffect+meaneffectintvar*((*betap)- (*linpredREp));
      level1_likep->compute_mu(&l,meffectp);
      }
    }

  level1_likep->meaneffect += meaneffect;

  }


void DESIGN_hrandom::outoptions(GENERAL_OPTIONS * op)
  {

  }


} // end: namespace MCMC



