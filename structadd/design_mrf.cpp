
#include "design_mrf.h"
#include "clstring.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//-------------- CLASS: DESIGN_mrf implementation of member functions ----------
//------------------------------------------------------------------------------


DESIGN_mrf::DESIGN_mrf(void) : DESIGN()
  {

  }

  // CONSTRUCTOR 1
  // Spatial covariates

DESIGN_mrf::DESIGN_mrf(const datamatrix & dm,const datamatrix & iv,
                       DISTR * dp, const MAP::map & m)
                      : DESIGN(dp)
  {
  ma = m;
  type = mrf;

  init_data(dm,iv);

  compute_penalty();

  compute_XtransposedWX_XtransposedWres(1);

  compute_precision(1.0);

  }



  // COPY CONSTRUCTOR

DESIGN_mrf::DESIGN_mrf(const DESIGN_mrf & m)
    : DESIGN(DESIGN(m))
  {
  ma = m.ma;
  data2 = m.data2;
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const DESIGN_mrf & DESIGN_mrf::operator=(const DESIGN_mrf & m)
  {
  if (this == &m)
    return *this;
  DESIGN::operator=(DESIGN(m));
  ma = m.ma;
  data2 = m.data2;
  return *this;
  }


void DESIGN_mrf::init_data(const datamatrix & dm, const datamatrix & iv)
  {

  datanames.push_back("X_1");


  Zout = datamatrix(ma.get_nrregions(),1,1);
  index_Zout = statmatrix<int>(Zout.rows(),1);
  index_Zout.indexinit();

  if (ma.get_bandsize() > 40)
    ma.reorderopt();

  ma.compute_reg(dm,posbeg,posend,effectvalues,index_data);

  unsigned i;
  data = datamatrix(dm.rows(),dm.cols());
  intvar = datamatrix(iv.rows(),iv.cols());
  int * workindex = index_data.getV();
  double * workdata = data.getV();
  for (i=0;i<data.rows();i++,workindex++,workdata++)
    *workdata = dm(*workindex,0);

  if (intvar.rows() == data.rows())
    {
    int * workindex = index_data.getV();
    double * workintvar = intvar.getV();
    data2 = datamatrix(dm.rows(),1);
    double * workdata2 = data2.getV();
    double h;
    for (i=0;i<intvar.rows();i++,workindex++,workdata2++,workintvar++)
      {
      h = iv(*workindex,0);
      *workintvar = h;
      *workdata2 = pow(h,2);
      }
    }

  make_partresindex();

  if (ma.get_errormessages().size() > 0)
    {
//  FEHLT!!
    }

  nrpar = ma.get_nrregions();

  }




void DESIGN_mrf::compute_penalty(void)
  {
  if (type==mrf)
    K = Kmrfenv(ma);
  rankK = ma.get_nrregions()-1;
  }


void DESIGN_mrf::compute_XtransposedWX_XtransposedWres(double l)
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

  if (intvar.rows() != data.rows())   // additive
    {
    for(i=0;i<posbeg.size();i++,++d,workXWres++)
      {
      *d=0;
      *workXWres = 0;
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
          *workXWres+= w*(*workintvar)*likep->partres(*workindex,0);
          }
        }
      }
    }

  }


void DESIGN_mrf::compute_XtransposedWres(datamatrix & partres, double l)
  {
  if (XWresdeclared == false)
    {
    XWres = datamatrix(nrpar,1);
    XWresdeclared = true;
    }

  int * workindex = index_data.getV();
  double * workXWres = XWres.getV();

  unsigned i,j;

  if (intvar.rows()!= data.rows())   // additive
    {
    for(i=0;i<nrpar;i++,workXWres++)
      {
      *workXWres = 0;
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
    double * workintvar = intvar.getV();
    for(i=0;i<nrpar;i++,workXWres++)
      {
      *workXWres = 0;
      if (posbeg[i] != -1)
        for(j=posbeg[i];j<=posend[i];j++,workindex++,workintvar++)
          *workXWres+= likep->workingweight(*workindex,0)*(*workintvar)*likep->partres(*workindex,0);
      }
    }

  }


void DESIGN_mrf::compute_precision(double l)
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



