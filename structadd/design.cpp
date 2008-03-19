
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

DESIGN::DESIGN(const datamatrix & dm, DISTR * lp)
  {
  data = dm;
  likep = lp;

  unsigned i;
  for (i=0;i<dm.cols();i++)
    datanames.push_back("X_"+ST::inttostring(i));

  index_data = statmatrix<int>(data.rows(),1);
  index_data.indexinit();
  dm.indexsort(index_data,0,dm.rows()-1,0,0);

  XWXdeclared = false;
  XWresdeclared = false;
  precisiondeclared=false;
  }

// COPY CONSTRUCTOR

DESIGN::DESIGN(const DESIGN & m)
  {
  data = m.data;
  data2 = m.data2;
  index_data = m.index_data;
  Z=m.Z;
  index_Z = m.index_Z;
  Zout = m.Zout;
  index_Zout=m.index_Zout;
  type=m.type;
  K = m.K;
  XWX = m.XWX;
  XWXdeclared = m.XWXdeclared;
  precision = m.precision;
  precisiondeclared = m.precisiondeclared;
  likep = m.likep;
  nrpar = m.nrpar;
  XWres = m.XWres;
  XWresdeclared = m.XWresdeclared;
  posbeg = m.posbeg;
  posend = m.posend;
  datanames = m.datanames;
  effectvalues = m.effectvalues;
  }

// OVERLOADED ASSIGNMENT OPERATOR

const DESIGN & DESIGN::operator=(const DESIGN & m)
  {
  if (this == &m)
    return *this;
  data = m.data;
  data2 = m.data2;
  index_data = m.index_data;
  Z=m.Z;
  index_Z = m.index_Z;
  Zout = m.Zout;
  index_Zout=m.index_Zout;
  type=m.type;
  K = m.K;
  XWX = m.XWX;
  XWXdeclared = m.XWXdeclared;
  precision = m.precision;
  precisiondeclared = m.precisiondeclared;
  likep = m.likep;
  nrpar = m.nrpar;
  XWres = m.XWres;
  XWresdeclared = m.XWresdeclared;
  posbeg = m.posbeg;
  posend = m.posend;
  datanames = m.datanames;
  effectvalues = m.effectvalues;
  return *this;
  }


void DESIGN::compute_design(void)
  {

  }


void DESIGN::compute_penalty(void)
  {

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


void DESIGN::compute_XtransposedWX_XtransposedWres(const datamatrix & res)
  {

  }


void DESIGN::compute_XtransposedWres(const datamatrix & res)
  {


  }


void DESIGN::compute_f(datamatrix & beta,datamatrix & f)
  {
  double * workf = f.getV();
  int * workindex = index_Zout.getV();

  double * workZout = Zout.getV();

  unsigned i,j;
  for(i=0;i<Zout.rows();i++,workf++)
    {
    *workf = 0;
    for(j=0;j<Zout.cols();j++,workindex++,workZout++)
      {
      *workf += (*workZout) * beta(*workindex,0);
      }
    }

  }

void DESIGN::update_linpred(datamatrix & f,bool add)
  {
  unsigned i,j;

  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();
  int * workindex = index_data.getV();
  double * workf = f.getV();

  datamatrix * linpredp = likep->linpred_current;

  if (add==true)
    {
    for (i=0;i<posbeg.size();i++,++itbeg,++itend,workf++)
      {
      for (j=*itbeg;j<=*itend;j++,workindex++)
        if (*itbeg != -1)
          (*linpredp)(*workindex,0) += *workf;
      }
    }
  else
    {
    for (i=0;i<posbeg.size();i++,++itbeg,++itend,workf++)
      {
      for (j=*itbeg;j<=*itend;j++,workindex++)
        if (*itbeg != -1)
          (*linpredp)(*workindex,0) -= *workf;
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

DESIGN_mrf::DESIGN_mrf(const datamatrix & dm,DISTR * dp, const MAP::map & m)
                      : DESIGN(dm,dp)
  {
  ma = m;
  type = mrf;

  if (data.cols() == 2)                    // interaction variable in second col
    {
    data2 = datamatrix(data.rows(),1);
    unsigned i;
    for(i=0;i<data.rows();i++)
      data2(i,0) = pow(data(i,1),2);
    }

  compute_design();

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
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const DESIGN_mrf & DESIGN_mrf::operator=(const DESIGN_mrf & m)
  {
  if (this == &m)
    return *this;
  DESIGN::operator=(DESIGN(m));
  ma = m.ma;
  return *this;
  }


void DESIGN_mrf::compute_design(void)
  {

  Zout = datamatrix(ma.get_nrregions(),1,1);
  index_Zout = statmatrix<int>(Zout.rows(),1);
  index_Zout.indexinit();

  if (ma.get_bandsize() > 40)
    ma.reorderopt();

  ma.compute_reg(data,posbeg,posend,effectvalues,index_data);

  if (ma.get_errormessages().size() > 0)
    {
//  FEHLT!!
    }

  nrpar = ma.get_nrregions();

  compute_Z();

  }


void DESIGN_mrf::compute_Z(void)
  {

  if (data.cols() > 1)
    {
    unsigned i,j;
    double vorher;
    double help;

    for (i=0;i<posbeg.size();i++)
      {

      data.indexsort(index_data,posbeg[i],posend[i],1,0);

      vorher = data(index_data(posbeg[i],0),1);
      j= posbeg[i]+1;
      posbeg_Z.push_back(posbeg[i]);
      while (j<=posend[i])
        {
        help =  data(index_data(j,0),1);
        if (help != vorher)
          {
          posend_Z.push_back(j);
          if (j !=posend[i])
            posbeg_Z.push_back(j+1);
          }

        if ((j==posend[i]) && (help == vorher))
          posend_Z.push_back(j);

        vorher = help;
        j++;

        }

      }

    Z = datamatrix(posbeg_Z.size(),1);
    index_Z = statmatrix<int>(Z.rows(),1);

    ofstream out0("c:\\bayesx\\test\\results\\posbeg_Z.res");

    int t1 = posbeg_Z.size();
    int t2 = posend_Z.size();

    for (i=0;i<posbeg_Z.size();i++)
      {
      out0 << posbeg_Z[i] << "   " << posend_Z[i] << endl;
      Z(i,0) = data(index_data(posbeg_Z[i],0),1);
      }



    int posc = 0;
    for (i=0;i<posbeg_Z.size();i++)
      {
      for (j=posbeg_Z[i];j<=posend_Z[i];j++)
        {
        if ((posc < posbeg.size()) && (j==posbeg[posc+1]))
          posc++;
        index_Z(i,0) = posc;
        }
      }

    ofstream out("c:\\bayesx\\test\\results\\index_Z.res");
    index_Z.prettyPrint(out);

    }
  else
    {
    Z = Zout;
    index_Z = index_Zout;
    posbeg_Z = posbeg;
    posend_Z = posend;
    }

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

  if (data.cols() == 1)   // additive
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
    double * workdata = data.getV()+1;
    for(i=0;i<posbeg.size();i++,++d,workXWres++)
      {
      *d=0;
      *workXWres = 0;
      if (posbeg[i] != -1)
        {
        for (j=posbeg[i];j<=posend[i];j++,workindex++,workdata2++,
             workdata+=2)
          {
          w = likep->workingweight(*workindex,0);
          *d += w*(*workdata2);
          *workXWres+= w*(*workdata)*res(*workindex,0);
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

  if (data.cols()==1)
    {
    for(i=0;i<nrpar;i++,workXWres++)
      {
      *workXWres = 0;
      if (posbeg[i] != -1)
        for(j=posbeg[i];j<=posend[i];j++,workindex++)
          {
          int h = *workindex;
          *workXWres+= likep->workingweight(*workindex,0)*res(*workindex,0);

          }
      }
    }
  else
    {
    double * workdata = data.getV()+1;
    for(i=0;i<nrpar;i++,workXWres++)
      {
      *workXWres = 0;
      if (posbeg[i] != -1)
        for(j=posbeg[i];j<=posend[i];j++,workindex++,workdata+=2)
          *workXWres+= likep->workingweight(*workindex,0)*(*workdata)*res(*workindex,0);
      }
    }

  }


} // end: namespace MCMC



