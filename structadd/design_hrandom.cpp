
#include "design_hrandom.h"
#include "clstring.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//------------ CLASS: DESIGN_hrandom implementation of member functions --------
//------------------------------------------------------------------------------


DESIGN_hrandom::DESIGN_hrandom(void)
  {

  }


DESIGN_hrandom::DESIGN_hrandom(const datamatrix & dm, const datamatrix & iv,
                               DISTR * dp, DISTR * dp_RE)
      : DESIGN(dp)
  {


  likep_RE = dp_RE;

  type = hrandom;

  init_data(dm,iv);

  compute_penalty();

  compute_XtransposedWX_XtransposedWres(1);

  compute_precision(1.0);

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

  make_index(dm);

  nrpar = posbeg.size();


  datanames.push_back("X_1");

  Zout = datamatrix(posbeg.size(),1,1);
  index_Zout = statmatrix<int>(Zout.rows(),1);
  index_Zout.indexinit();

  }




void DESIGN_hrandom::compute_penalty(void)
  {
  K =   envmatrix<double>(1,nrpar);
  rankK = nrpar;
  }


void DESIGN_hrandom::compute_XtransposedWX_XtransposedWres(double l)
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


  unsigned i;
  int j;
  int *  workindex = index_data.getV();
  double * workXWres = XWres.getV();
  double w;

  vector<double>::iterator d = XWX.getDiagIterator();

  datamatrix * linpredRE = likep_RE->linpred_current;
  double * linpredREp = (*linpredRE).getV();

  if (intvar.rows() != data.rows())   // additive
    {
    for(i=0;i<posbeg.size();i++,++d,workXWres++,linpredREp++)
      {
      *d=0;
      *workXWres =  l*(*linpredREp);
      if (posbeg[i] != -1)
        {
        for (j=posbeg[i];j<=posend[i];j++,workindex++)
          {
          w = likep->workingweight(*workindex,0);
          *d += w;
          *workXWres+= w*likep->partres(*workindex,0);
          }

        }
      }
    }
  else                    // varying coefficients
    {
    /*
    double * workdata2 = data2.getV();
    double * workintvar = intvar.getV();
    for(i=0;i<posbeg.size();i++,++d,workXWres++)
      {
      *d=0;
      *workXWres = 0;
      if (posbeg[i] != -1)
        {
        for (j=posbeg[i];j<=posend[i];j++,workindex++,workdata2++,
             workintvar++)
          {
          w = likep->workingweight(*workindex,0);
          *d += w * (*workdata2);
          *workXWres+= w*(*workintvar)*res(*workindex,0);
          }
        }
      }
    */
    }

  }


void DESIGN_hrandom::compute_XtransposedWres(datamatrix & partres, double l)
  {

  if (XWresdeclared == false)
    {
    XWres = datamatrix(nrpar,1);
    XWresdeclared = true;
    }


  int * workindex = index_data.getV();
  double * workXWres = XWres.getV();

  datamatrix * linpredRE = likep_RE->linpred_current;
  double * linpredREp = (*linpredRE).getV();

  unsigned i,j;

  if (intvar.rows()!= data.rows())   // additive
    {
    for(i=0;i<nrpar;i++,workXWres++,linpredREp++)
      {
      *workXWres =  l*(*linpredREp);
      if (posbeg[i] != -1)
        {
        for(j=posbeg[i];j<=posend[i];j++,workindex++)
          {
          *workXWres+= likep->workingweight(*workindex,0)*likep->partres(*workindex,0);
          }
        }
      }
    }
  else                              // varying coefficient
    {
/*
    double * workintvar = intvar.getV();
    for(i=0;i<nrpar;i++,workXWres++)
      {
      *workXWres = 0;
      if (posbeg[i] != -1)
        for(j=posbeg[i];j<=posend[i];j++,workindex++,workintvar++)
          *workXWres+= likep->workingweight(*workindex,0)*(*workintvar)*res(*workindex,0);
      }
*/
    }

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



} // end: namespace MCMC



