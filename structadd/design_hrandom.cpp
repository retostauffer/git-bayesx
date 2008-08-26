
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

  int f;

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
  compute_XtransposedWX_XtransposedWres(help,1);

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

  if (workingresponsep.rows() != data.rows())
    {
    make_pointerindex();
    }

  unsigned i;
  int j;

  double * * workingweightpp = workingweightp.getV();
  vector<double>::iterator d = XWX.getDiagIterator();


  if (intvar.rows() != data.rows())   // additive
    {
    for(i=0;i<nrpar;i++,++d)
      {
      *d=0;
      if (posbeg[i] != -1)
        {
        for (j=posbeg[i];j<=posend[i];j++,workingweightpp++)
          {
          *d += *(*workingweightpp);
          }

        }
      }

    }
  else                    // varying coefficients
    {

    double * workdata2 = intvar2.getV();
    for(i=0;i<nrpar;i++,++d)
      {
      *d=0;
      if (posbeg[i] != -1)
        {
        for (j=posbeg[i];j<=posend[i];j++,workdata2++,workingweightpp++)
          {
          *d += *(*workingweightpp) * (*workdata2);
          }
        }
      }

    }

  }


void DESIGN_hrandom::compute_XtransposedWX_XtransposedWres(
                                                         datamatrix & partres,
                                                         double l)
  {

  if (XWXdeclared == false)
    {
    XWX = envmatdouble(0,nrpar);
    XWXdeclared = true;
    }

  if (XWresdeclared == false)
    {
    XWres = datamatrix(nrpar,1);
    XWresdeclared = true;
    }

  if (workingresponsep.rows() != data.rows())
    {
    make_pointerindex();
    }

  unsigned i;
  int j;

  double * * workingweightpp = workingweightp.getV();
  vector<double>::iterator d = XWX.getDiagIterator();

  double * workXWres = XWres.getV();

  double * linpredREp;
  if (likep_RE->linpred_current==1)
    linpredREp = likep_RE->linearpred1.getV();
  else
    linpredREp = likep_RE->linearpred2.getV();

  double * partresp = partres.getV();

  if (intvar.rows() != data.rows())   // additive
    {
    for(i=0;i<nrpar;i++,++d,workXWres++,linpredREp++,partresp++)
      {
      *d=0;
      *workXWres =  l*(*linpredREp)+(*partresp);
      if (posbeg[i] != -1)
        {
        for (j=posbeg[i];j<=posend[i];j++,workingweightpp++)
          {
          *d += *(*workingweightpp);
          }

        }
      }

    }
  else                    // varying coefficients
    {

    double * workdata2 = intvar2.getV();
    for(i=0;i<nrpar;i++,++d,workXWres++,linpredREp++,partresp++)
      {
      *d=0;
      *workXWres =  l*(*linpredREp)+(*partresp);
      if (posbeg[i] != -1)
        {
        for (j=posbeg[i];j<=posend[i];j++,workdata2++,workingweightpp++)
          {
          *d += *(*workingweightpp) * (*workdata2);
          }
        }
      }

    }

// TEST
  /*
  ofstream out("c:\\bayesx\\test\\results\\XWX.res");
  XWX.print2(out);
  */
// TEST

  }


void DESIGN_hrandom::compute_XtransposedWres(datamatrix & partres, double l)
  {

  if (XWresdeclared == false)
    {
    XWres = datamatrix(nrpar,1);
    XWresdeclared = true;
    }

  if (workingresponsep.rows() != data.rows())
    {
    make_pointerindex();
    }

  double * workXWres = XWres.getV();


  double * linpredREp;
  if (likep_RE->linpred_current==1)
    linpredREp = likep_RE->linearpred1.getV();
  else
    linpredREp = likep_RE->linearpred2.getV();

  double * partresp = partres.getV();

  unsigned i,j;

  for(i=0;i<nrpar;i++,workXWres++,linpredREp++,partresp++)
    *workXWres =  l*(*linpredREp)+(*partresp);

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




void DESIGN_hrandom::outoptions(GENERAL_OPTIONS * op)
  {

  }


} // end: namespace MCMC



