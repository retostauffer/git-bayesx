
#include "design.h"
#include "clstring.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//---------------- CLASS: DESIGN implementation of member functions ------------
//------------------------------------------------------------------------------


// DEFAULT CONSTRUCTOR

DESIGN::DESIGN(void)
  {
  data = datamatrix(1,1,0);
  }

// CONSTRUCTOR

DESIGN::DESIGN(DISTR * lp)
  {

  likep = lp;

  XWXdeclared = false;
  XWresdeclared = false;
  precisiondeclared=false;
  }


// COPY CONSTRUCTOR

DESIGN::DESIGN(const DESIGN & m)
  {
  likep = m.likep;

  data = m.data;
  intvar = m.intvar;
  index_data = m.index_data;
  datanames = m.datanames;
  effectvalues = m.effectvalues;

  Zout = m.Zout;
  index_Zout=m.index_Zout;
  posbeg = m.posbeg;
  posend = m.posend;

  nrpar = m.nrpar;

  K = m.K;
  rankK = m.rankK;

  XWX = m.XWX;
  XWXdeclared = m.XWXdeclared;
  precision = m.precision;
  precisiondeclared = m.precisiondeclared;


  XWres = m.XWres;
  XWresdeclared = m.XWresdeclared;

  type=m.type;

  }


// OVERLOADED ASSIGNMENT OPERATOR

const DESIGN & DESIGN::operator=(const DESIGN & m)
  {
  if (this == &m)
    return *this;
  likep = m.likep;

  data = m.data;
  intvar = m.intvar;
  index_data = m.index_data;
  datanames = m.datanames;
  effectvalues = m.effectvalues;

  Zout = m.Zout;
  index_Zout=m.index_Zout;
  posbeg = m.posbeg;
  posend = m.posend;

  nrpar = m.nrpar;

  K = m.K;
  rankK = m.rankK;

  XWX = m.XWX;
  XWXdeclared = m.XWXdeclared;
  precision = m.precision;
  precisiondeclared = m.precisiondeclared;


  XWres = m.XWres;
  XWresdeclared = m.XWresdeclared;

  type=m.type;

  return *this;
  }


void DESIGN::init_data(datamatrix & dm, datamatrix & iv)
  {

  }




void DESIGN::compute_penalty(void)
  {

  }



void DESIGN::compute_XtransposedWX_XtransposedWres(const datamatrix & res)
  {

  }


void DESIGN::compute_XtransposedWres(const datamatrix & res)
  {


  }


void DESIGN::compute_f(datamatrix & beta,datamatrix & f)
  {
  double * workf = f.getV();

  int * workindex;
  double * workZ;
  unsigned rows;
  unsigned cols;

  workindex = index_Zout.getV();
  workZ = Zout.getV();
  rows = Zout.rows();
  cols = Zout.cols();

  unsigned i,j;
  for(i=0;i<rows;i++,workf++)
    {
    *workf = 0;
    for(j=0;j<cols;j++,workindex++,workZ++)
      {
      *workf += (*workZ) * beta(*workindex,0);
      }
    }

  }

void DESIGN::compute_precision(double l)
  {
  if (precisiondeclared==false)
    {
    precision = envmatdouble(K.getXenv(),0,nrpar);
    precisiondeclared = true;
    }

  precision.addtodiag(XWX,K,1.0,l);

  }


void DESIGN::update_linpred(datamatrix & f,bool add)
  {
  unsigned i,j;

  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();
  int * workindex = index_data.getV();
  double * workf = f.getV();
  double * workintvar = intvar.getV();

  datamatrix * linpredp = likep->linpred_current;

  if (add==true)
    {
    if (intvar.rows()==data.rows())   // varying coefficient
      {
      for (i=0;i<posbeg.size();i++,++itbeg,++itend,workf++)
        {
        if (*itbeg != -1)
          {
          for (j=*itbeg;j<=*itend;j++,workindex++,workintvar++)
//            if (*workintvar !=0)
              (*linpredp)(*workindex,0) += (*workintvar) * (*workf);
          }
        }
      }
    else                              // additive
      {
      for (i=0;i<posbeg.size();i++,++itbeg,++itend,workf++)
        {
        if (*itbeg != -1)
          {
          for (j=*itbeg;j<=*itend;j++,workindex++)
            (*linpredp)(*workindex,0) += *workf;
          }
        }
      }
    }
  else
    {

    if (intvar.rows()==data.rows())   // varying coefficient
      {
      for (i=0;i<posbeg.size();i++,++itbeg,++itend,workf++)
        {
        if (*itbeg != -1)
          {
          for (j=*itbeg;j<=*itend;j++,workindex++,workintvar++)
            (*linpredp)(*workindex,0) -= (*workintvar) * (*workf);
          }
        }
      }
    else                              // additive
      {
      for (i=0;i<posbeg.size();i++,++itbeg,++itend,workf++)
        {
        if (*itbeg != -1)
          {
          for (j=*itbeg;j<=*itend;j++,workindex++)
            (*linpredp)(*workindex,0) -= *workf;
          }
        }
      }

    }

  }



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

  datamatrix help(dm.rows(),1,1);
  compute_XtransposedWX_XtransposedWres(help);

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


void DESIGN_mrf::init_data(datamatrix & dm, datamatrix & iv)
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


void DESIGN_mrf::compute_XtransposedWX_XtransposedWres(const datamatrix & res)
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
          *workXWres+= w*res(*workindex,0);
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
          *workXWres+= w*(*workintvar)*res(*workindex,0);
          }
        }
      }
    }

  }


void DESIGN_mrf::compute_XtransposedWres(const datamatrix & res)
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
          *workXWres+= likep->workingweight(*workindex,0)*res(*workindex,0);
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
          *workXWres+= likep->workingweight(*workindex,0)*(*workintvar)*res(*workindex,0);
      }
    }

  }


} // end: namespace MCMC



