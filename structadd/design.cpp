
#include "design.h"
#include "clstring.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//---------------- CLASS: DESIGN implementation of member functions ------------
//------------------------------------------------------------------------------

void DESIGN::make_index(const datamatrix & dm)
  {

  unsigned j;

  index_data = statmatrix<int>(dm.rows(),1);
  index_data.indexinit();
  dm.indexsort(index_data,0,dm.rows()-1,0,0);

  data = datamatrix(dm.rows(),1);
  double * workdata = data.getV();
  int * workindex = index_data.getV();
  for (j=0;j<dm.rows();j++,workdata++,workindex++)
    {
    *workdata = dm(*workindex,0);
    }

  posbeg.push_back(0);
  workdata = data.getV()+1;
  double help = data(0,0);
  for(j=1;j<data.rows();j++,workdata++)
    {
    if (  *workdata != help)
      {
      posend.push_back(j-1);
      if (j < data.rows())
        posbeg.push_back(j);
      }

    help = *workdata;

    }


  if (posend.size() < posbeg.size())
    posend.push_back(data.rows()-1);

  for(j=0;j<posbeg.size();j++)
    effectvalues.push_back(ST::doubletostring(data(posbeg[j],0)));

  int ev = effectvalues.size();

  /*
  // TEST
  ofstream out("c:\\bayesx\\test\\results\\posbeg.res");
  for(j=0;j<posbeg.size();j++)
    out << posbeg[j] << "  " << posend[j] << endl;
  // TEST
  */

  }


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

  ZoutT = m.ZoutT;
  index_ZoutT = m.index_ZoutT;

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

  ZoutT = m.ZoutT;
  index_ZoutT = m.index_ZoutT;

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


void DESIGN::init_data(const datamatrix & dm, const datamatrix & iv)
  {

  }




void DESIGN::compute_penalty(void)
  {

  }


void DESIGN::compute_Zout_transposed(void)
  {

  vector<double> h;
  ZoutT = vector<vector<double> >(nrpar,h);

  vector<int> h2;
  index_ZoutT = vector<vector<int> >(nrpar,h2);


  unsigned i,j;

  for (i=0;i<Zout.rows();i++)
    for(j=0;j<Zout.cols();j++)
      {
      ZoutT[index_Zout(i,j)].push_back(Zout(i,j));
      index_ZoutT[index_Zout(i,j)].push_back(i);
      }

   /*
  // TEST
  ofstream out("c:\\bayesx\\test\\results\\ZoutT.res");
  for (i=0;i<ZoutT.size();i++)
    {
    for(j=0;j<ZoutT[i].size();j++)
      out <<  ZoutT[i][j] << "  ";
    out << endl;
    }

  ofstream out2("c:\\bayesx\\test\\results\\ZoutT_index.res");
  for (i=0;i<index_ZoutT.size();i++)
    {
    for(j=0;j<index_ZoutT[i].size();j++)
      out2 <<  index_ZoutT[i][j] << "  ";
    out2 << endl;
    }
  // TEST
  */

  }



void DESIGN::compute_XtransposedWX_XtransposedWres(const datamatrix & res,double l)
  {

  unsigned i,j,k;

  vector<double>::iterator diag = XWX.getDiagIterator();
  double help;
  int ip;
  double wsum;
  double t;
  for (i=0;i<nrpar;i++,++diag)
    {
    *diag=0;

    for (j=0;j<ZoutT[i].size();j++)
      {
      help=pow(ZoutT[i][j],2);
      ip = index_ZoutT[i][j];
      wsum=0;
      for (k=posbeg[ip];k<=posend[ip];k++)
        wsum+= likep->workingweight(index_data(k,0),0);

      *diag += help*wsum;
      }

    }

  vector<double>::iterator env = XWX.getEnvIterator();
  vector<unsigned>::iterator xenv = XWX.getXenvIterator();
  unsigned start = *xenv;
  unsigned nrnnull;
  xenv++;

  unsigned envs = XWX.getXenv(nrpar);
  for(i=0;i<nrpar;i++,++xenv)
    {
    nrnnull = *xenv-start;
    if (nrnnull > 0)
      {
      for (j=i-nrnnull;j<i;j++,++env)
        {
        *env = compute_ZtZ(i,j);
        }

      }
    start = *xenv;
    }


  XWX.setDecomposed(false);

  /*
  // TEST
  vector<double> env = XWX.getEnv();

  ofstream out2("c:\\bayesx\\test\\results\\XWXenv.res");
  for (j=0;j<env.size();j++)
    out2 << env[j] << endl;


  vector<unsigned> Xenv = XWX.getXenv();

  ofstream out3("c:\\bayesx\\test\\results\\XWX_Xenv.res");
  for (j=0;j<Xenv.size();j++)
    out3 << Xenv[j] << endl;

  ofstream out("c:\\bayesx\\test\\results\\XWX.res");
  XWX.print2(out);
  // TEST
  */

  }


double DESIGN::compute_ZtZ(unsigned & i, unsigned & j)
  {
  unsigned l;
  unsigned k_i=0;
  unsigned k_j=0;
  int beg_i;
  int end_i;
  int beg_j;
  int end_j;
  int beg;
  int end;
  double result = 0;
  double sumw;

  while (k_i<ZoutT[i].size() && k_j < ZoutT[j].size())
    {

    beg_i = posbeg[index_ZoutT[i][k_i]];
    end_i = posend[index_ZoutT[i][k_i]];
    beg_j = posbeg[index_ZoutT[j][k_j]];
    end_j = posend[index_ZoutT[j][k_j]];

    if (beg_j > end_i)
      {
      k_i++;
      }
    else if (end_j < beg_i)
      {
      k_j++;
      }
    else  // overlapp
      {
      sumw=0;
      if (beg_j >= beg_i)
        {
        beg = beg_j;
        }
      else
        {
        beg = beg_i;
        }

      if (end_j >= end_i)
        {
        end = end_i;
        }
      else
        {
        end = end_j;
        }


      for (l=beg;l<=end;l++)
        sumw += likep->workingweight(index_data(l,0),0);

      result += sumw * ZoutT[i][k_i]* ZoutT[j][k_j];
      if (end_i <= end_j)
        k_i++;
      else
        k_j++;
      }


    }

  return result;
  }

void DESIGN::compute_XtransposedWres(const datamatrix & res,double l)
  {

  unsigned i,j,k;

  if (XWresdeclared == false)
    {
    XWres = datamatrix(nrpar,1);
    XWresdeclared = true;
    }

   /*
  // TEST
  datamatrix Zoutm(data.rows(),nrpar,0);
  for (i=0;i<posbeg.size();i++)
    {
    for (j=posbeg[i];j<=posend[i];j++)
      {
      for(k=0;k<Zout.cols();k++)
        Zoutm(j,index_Zout(i,k)) = Zout(i,k);
      }

    }

  ofstream out4("c:\\bayesx\\test\\results\\Zoutm.res");
  Zoutm.prettyPrint(out4);

  datamatrix ressort(res.rows(),1);
  for (i=0;i<res.rows();i++)
    ressort(i,0) = res(index_data(i,0),0);

  datamatrix Ztres =  Zoutm.transposed()*ressort;

  ofstream out3("c:\\bayesx\\test\\results\\XWres_alt.res");
  Ztres.prettyPrint(out3);
  // END TEST
  */

  if (ZoutT.size() != nrpar)
    compute_Zout_transposed();

  int * workindex = index_data.getV();
  double * workXWres = XWres.getV();


  int beg;
  int end;
  for(i=0;i<nrpar;i++,workXWres++)
    {
    *workXWres=0;
    for (j=0;j<ZoutT[i].size();j++)
      {
      beg = posbeg[index_ZoutT[i][j]];
      end = posend[index_ZoutT[i][j]];
      for (k=beg;k<=end;k++)
        {
        *workXWres += ZoutT[i][j] *
        likep->workingweight(index_data(k,0),0) * res(index_data(k,0),0);

        }
      }

    }

  /*
  // TEST
  ofstream out("c:\\bayesx\\test\\results\\XWres.res");
  XWres.prettyPrint(out);
  // TEST
  */

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

  /*
  // TEST
  ofstream out("c:\\bayesx\\test\\results\\f.res");
  f.prettyPrint(out);

  datamatrix Zoutm(data.rows(),nrpar,0);
  unsigned k;
  for (i=0;i<posbeg.size();i++)
    {
    for (j=posbeg[i];j<=posend[i];j++)
      {
      for(k=0;k<Zout.cols();k++)
        Zoutm(j,index_Zout(i,k)) = Zout(i,k);
      }

    }

  ofstream out2("c:\\bayesx\\test\\results\\falt.res");

  datamatrix falt = Zoutm*beta;

  falt.prettyPrint(out2);
  // TEST
  */

  }

void DESIGN::compute_precision(double l)
  {

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


  /*
  // TEST
  ofstream out("c:\\bayesx\\test\\results\\lin.res");
  (*linpredp).prettyPrint(out);
  // TEST
  */

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
  compute_XtransposedWX_XtransposedWres(help,1);

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


void DESIGN_mrf::compute_XtransposedWX_XtransposedWres(const datamatrix & res,
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


void DESIGN_mrf::compute_XtransposedWres(const datamatrix & res,double l)
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

  datamatrix help(dm.rows(),1,1);
  compute_XtransposedWX_XtransposedWres(help,1);

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


void DESIGN_hrandom::compute_XtransposedWX_XtransposedWres(const datamatrix & res,double l)
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
          *workXWres+= w*res(*workindex,0);
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


void DESIGN_hrandom::compute_XtransposedWres(const datamatrix & res,double l)
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
          *workXWres+= likep->workingweight(*workindex,0)*res(*workindex,0);
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



