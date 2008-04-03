
#include "design.h"
#include "clstring.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//---------------- CLASS: DESIGN implementation of member functions ------------
//------------------------------------------------------------------------------

void DESIGN::make_partresindex(void)
  {
  unsigned i;
  partres_pindex = statmatrix<double*> (data.rows(),1);
//  double work_pi =
  for(i=0;i<data.rows();i++)
    partres_pindex(i,0) = &(likep->partres(index_data(i,0),0));
  }

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

  make_partresindex();

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

  partres_pindex = m.partres_pindex;

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

  partres_pindex = m.partres_pindex;

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



void DESIGN::compute_XtransposedWX_XtransposedWres(datamatrix & partres, double l)
  {

  unsigned i,j,k;

  vector<double>::iterator diag = XWX.getDiagIterator();
  double help;
  int ip;
  double wsum;

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

  compute_XtransposedWres(partres, l);

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


void DESIGN::compute_XtransposedWres(datamatrix & partres, double l)
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


  double * workXWres = XWres.getV();
  double * workpartres = partres.getV();

  int beg;
  int end;
  double help;
  for(i=0;i<nrpar;i++,workXWres++,workpartres++)
    {
    *workXWres=0;

    for (j=0;j<ZoutT[i].size();j++)
      {
      *workXWres += ZoutT[i][j]* (*workpartres);
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




void DESIGN::compute_partres(datamatrix & res, datamatrix & f)
  {

  unsigned i,j;
  unsigned size = posbeg.size();
  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();
  int * workindex = index_data.getV();
  double * workf = f.getV();
  double * workintvar = intvar.getV();

  datamatrix * linpredp = likep->linpred_current;
  double * workres = res.getV();

  if (intvar.rows()==data.rows())   // varying coefficient
    {
    for (i=0;i<size;i++,++itbeg,++itend,workf++)
      {
      if (*itbeg != -1)
        {
        for (j=*itbeg;j<=*itend;j++,workindex++,workintvar++)
          *workres += likep->workingweight(*workindex,0)* (*workintvar)*
          (likep->response(*workindex,0)- (*linpredp)(*workindex,0) + *workf);
        }
      }
    }
  else                              // additive
    {
    for (i=0;i<size;i++,++itbeg,++itend,workf++,workres++)
      {
      *workres = 0;
      if (*itbeg != -1)
        {
        for (j=*itbeg;j<=*itend;j++,workindex++)
          {
          *workres += likep->workingweight(*workindex,0)*
          (likep->response(*workindex,0)- (*linpredp)(*workindex,0) + *workf);
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




} // end: namespace MCMC



