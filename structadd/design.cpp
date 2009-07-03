
#include "design.h"
#include "clstring.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//---------------- CLASS: DESIGN implementation of member functions ------------
//------------------------------------------------------------------------------


void DESIGN::make_pointerindex(void)
  {
  unsigned i;
  workingresponsep = statmatrix<double*> (data.rows(),1);
  for(i=0;i<data.rows();i++)
    workingresponsep(i,0) = &(likep->workingresponse(index_data(i,0),0));

  workingweightp = statmatrix<double*> (data.rows(),1);
  for(i=0;i<data.rows();i++)
    workingweightp(i,0) = &(likep->workingweight(index_data(i,0),0));

  responsep = statmatrix<double*> (data.rows(),1);
  for(i=0;i<data.rows();i++)
    responsep(i,0) = &(likep->response(index_data(i,0),0));

  weightp = statmatrix<double*> (data.rows(),1);
  for(i=0;i<data.rows();i++)
    weightp(i,0) = &(likep->weight(index_data(i,0),0));

  linpredp1 = statmatrix<double*> (data.rows(),1);
  linpredp2 = statmatrix<double*> (data.rows(),1);

  for(i=0;i<data.rows();i++)
    linpredp1(i,0) = &(likep->linearpred1(index_data(i,0),0));

  for(i=0;i<data.rows();i++)
    linpredp2(i,0) = &(likep->linearpred2(index_data(i,0),0));

  }


bool DESIGN::check_Zout_consecutive(void)
  {
  bool cons = true;

  unsigned rows;
  unsigned cols;
  rows = Zout.rows();
  cols = Zout.cols();

  unsigned i,j;
  for(i=0;i<rows;i++)
    {
    for(j=1;j<cols;j++)
      {
      if (index_Zout(i,j) - index_Zout(i,j-1) > 1)
        cons=false;
      }
    }

  return cons;

  }


bool DESIGN::check_ZoutT_consecutive(void)
  {
  bool cons = true;
  int size;
  int i,j;

  for(i=0;i<int(nrpar);i++)
    {
    size = ZoutT[i].size();
    for (j=1;j<size;j++)
      {
      if (index_ZoutT[i][j] - index_ZoutT[i][j-1] > 1)
        cons = false;
      }
    }

  return cons;

  }


void DESIGN::make_data(const datamatrix & dm,const datamatrix & iv)
  {
  unsigned j;
  data = datamatrix(dm.rows(),1);
  double * workdata = data.getV();
  int * workindex = index_data.getV();
  for (j=0;j<dm.rows();j++,workdata++,workindex++)
    {
    *workdata = dm(*workindex,0);
    }

  if (iv.rows() == dm.rows())
    {
    intvar = datamatrix(iv.rows(),1);
    intvar2 = datamatrix(iv.rows(),1);
    double * workintvar = intvar.getV();
    double * workintvar2 = intvar2.getV();
    workindex = index_data.getV();
    for (j=0;j<iv.rows();j++,workintvar++,workindex++,workintvar2++)
      {
      *workintvar = iv(*workindex,0);
      *workintvar2 = pow(*workintvar,2);
      }
    }

  }

  
double DESIGN::compute_sumBk(unsigned & k)
  {

  double sum=0;
  unsigned c;

  unsigned j;
    for (j=0;j<ZoutT[k].size();j++)
      {
      c = index_ZoutT[k][j];
      sum += ZoutT[k][j]*(posend[c]-posbeg[c]+1);
      }

  return sum;
  }

void DESIGN::make_index(const datamatrix & dm,const datamatrix & iv)
  {

  unsigned j;

  index_data = statmatrix<int>(dm.rows(),1);
  index_data.indexinit();
  dm.indexsort(index_data,0,dm.rows()-1,0,0);

  double dm_mean = dm.mean(0);
//  double iv_mean;
//  if (iv.rows() == dm.rows())
//    iv_mean = iv.mean(0);

  make_data(dm,iv);

  double * workdata;
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


  double d;
  meaneffectnr = 0;
  double mclosest = data(posbeg[0],0);
  for(j=0;j<posbeg.size();j++)
    {
    d = data(posbeg[j],0);
    if ( fabs(d-dm_mean) < fabs(mclosest-dm_mean) )
      {
      meaneffectnr = j;
      mclosest = d;
      }

    effectvalues.push_back(ST::doubletostring(d,0));

    }

/*
  meaneffectintvar=1;
  if (iv.rows() == dm.rows())
    {
    double * intvarp = intvar.getV();
    meaneffectnr_intvar = 0;
    mclosest = *intvarp;
    for (j=0;j<intvar.rows();j++,intvarp++)
      {

      if ( fabs(*intvarp-iv_mean) < fabs(mclosest-iv_mean) )
        {
        meaneffectnr_intvar = j;
        mclosest = *intvarp;
        meaneffectintvar= *intvarp;
        }

      }

    }
*/

  compute_meaneffectintvar();


  /*
  // TEST
  int ev = effectvalues.size();

  ofstream out("c:\\bayesx\\test\\results\\posbeg.res");
  for(j=0;j<posbeg.size();j++)
    out << posbeg[j] << "  " << posend[j] << endl;
  // TEST
  */

  }


void DESIGN::compute_meaneffectintvar(void)
  {
  meaneffectintvar=1;
  if (data.rows() == intvar.rows())
    {

    double iv_mean = intvar.mean(0);

    double * intvarp = intvar.getV();
    meaneffectnr_intvar = 0;
    double mclosest = *intvarp;
    unsigned j;
    for (j=0;j<intvar.rows();j++,intvarp++)
      {

      if ( fabs(*intvarp-iv_mean) < fabs(mclosest-iv_mean) )
        {
        meaneffectnr_intvar = j;
        mclosest = *intvarp;
        meaneffectintvar= *intvarp;
        }

      }

    }
  }

  
unsigned DESIGN::compute_modecategorie(void)
  {

  unsigned j;

  double d;
  unsigned modenr = 0;
  double mnr = 0;
  for(j=0;j<posbeg.size();j++)
    {
    d = (posend[j] - posbeg[j]+1);
    if ( d > mnr)
      {
      modenr = j;
      mnr = d;
      }

    }

  return modenr;
  }



// DEFAULT CONSTRUCTOR

DESIGN::DESIGN(void)
  {
  data = datamatrix(1,1,0);
  }

// CONSTRUCTOR

DESIGN::DESIGN(DISTR * lp,FC_linear * fcp)
  {

  changingdesign = false;

  likep = lp;
  FClinearp = fcp;

  XWXdeclared = false;
  XWresdeclared = false;
  precisiondeclared=false;
  consecutive = -1;
  consecutive_ZoutT = -1;
  identity = false;

  position_lin = -1;
  }


// COPY CONSTRUCTOR

DESIGN::DESIGN(const DESIGN & m)
  {
  changingdesign = m.changingdesign;
  likep = m.likep;

  data = m.data;
  intvar2=m.intvar2;
  intvar = m.intvar;
  index_data = m.index_data;
  datanames = m.datanames;
  effectvalues = m.effectvalues;
  meaneffectnr = m.meaneffectnr;
  meaneffectnr_intvar = m.meaneffectnr_intvar;
  meaneffectintvar = m.meaneffectintvar;  

  Zout = m.Zout;
  index_Zout=m.index_Zout;
  posbeg = m.posbeg;
  posend = m.posend;
  consecutive = m.consecutive;
  consecutive_ZoutT = m.consecutive_ZoutT;
  identity = m.identity;

  ZoutT = m.ZoutT;
  index_ZoutT = m.index_ZoutT;

  ZoutTZout = m.ZoutTZout;
  beg_ZoutTZout = m.beg_ZoutTZout;
  Wsump = m.Wsump;

  ZoutTZout_d = m.ZoutTZout_d;
  Wsump_d = m.Wsump_d;

  responsep = m.responsep;
  weightp = m.weightp;
  workingresponsep = m.workingresponsep;
  workingweightp = m.workingweightp;
  linpredp1 = m.linpredp1;
  linpredp2 = m.linpredp2;

  nrpar = m.nrpar;

  center = m.center;
  centermethod=m.centermethod;
  basisNull = m.basisNull;
  basisNullt = m.basisNullt;
  FClinearp = m.FClinearp;
  position_lin = m.position_lin;
  designlinear = m.designlinear;

  K = m.K;
  rankK = m.rankK;

  XWX = m.XWX;
  XWXdeclared = m.XWXdeclared;
  precision = m.precision;
  precisiondeclared = m.precisiondeclared;
  Wsum = m.Wsum;

  XWres = m.XWres;
  XWresdeclared = m.XWresdeclared;
  XWres_p = m.XWres_p;

  type=m.type;

  }


// OVERLOADED ASSIGNMENT OPERATOR

const DESIGN & DESIGN::operator=(const DESIGN & m)
  {
  if (this == &m)
    return *this;

  changingdesign = m.changingdesign;

  likep = m.likep;

  data = m.data;
  intvar2=m.intvar2;
  intvar = m.intvar;
  index_data = m.index_data;
  datanames = m.datanames;
  effectvalues = m.effectvalues;
  meaneffectnr = m.meaneffectnr;
  meaneffectnr_intvar = m.meaneffectnr_intvar;

  Zout = m.Zout;
  index_Zout=m.index_Zout;
  posbeg = m.posbeg;
  posend = m.posend;
  consecutive = m.consecutive;
  consecutive_ZoutT = m.consecutive_ZoutT;
  identity = m.identity;

  ZoutT = m.ZoutT;
  index_ZoutT = m.index_ZoutT;

  ZoutTZout = m.ZoutTZout;
  beg_ZoutTZout = m.beg_ZoutTZout;
  Wsump = m.Wsump;

  ZoutTZout_d = m.ZoutTZout_d;
  Wsump_d = m.Wsump_d;

  responsep = m.responsep;
  weightp = m.weightp;
  workingresponsep = m.workingresponsep;
  workingweightp = m.workingweightp;
  linpredp1 = m.linpredp1;
  linpredp2 = m.linpredp2;

  nrpar = m.nrpar;

  center = m.center;
  centermethod=m.centermethod;
  basisNull = m.basisNull;
  basisNullt = m.basisNullt;
  FClinearp = m.FClinearp;
  position_lin = m.position_lin;
  designlinear = m.designlinear;

  K = m.K;
  rankK = m.rankK;

  XWX = m.XWX;
  XWXdeclared = m.XWXdeclared;
  precision = m.precision;
  precisiondeclared = m.precisiondeclared;
  Wsum = m.Wsum;

  XWres = m.XWres;
  XWresdeclared = m.XWresdeclared;
  XWres_p = m.XWres_p;

  type=m.type;

  return *this;
  }



void DESIGN::init_data(const datamatrix & dm, const datamatrix & iv)
  {
  unsigned j;

  double * workintvar = intvar.getV();
  double * workintvar2 = intvar2.getV();
  int * workindex = index_data.getV();
  for (j=0;j<iv.rows();j++,workintvar++,workindex++,workintvar2++)
    {
    *workintvar = iv(*workindex,0);
    *workintvar2 = pow(*workintvar,2);
    }

  }


void DESIGN::compute_penalty(void)
  {

  }


void DESIGN::compute_basisNull(void)
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


  // TEST
  /*
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

  ofstream out4("c:\\bayesx\\testh\\results\\Z.res");
  datamatrix Z(Zout.rows(),ZoutT.size(),0);
  for (i=0;i<ZoutT.size();i++)
    {
    for(j=0;j<ZoutT[i].size();j++)
      Z(index_ZoutT[i][j],i) = ZoutT[i][j];
    }

  Z.prettyPrint(out4);
  */
  // TEST

  }


void DESIGN::compute_XtransposedWX_XtransposedWres(datamatrix & partres, double l)
  {

  compute_XtransposedWX();
  compute_XtransposedWres(partres, l);

  }


void DESIGN::compute_XtransposedWX(void)
  {

  // TEST
//  ofstream out("c:\\bayesx\\test\\results\\workingweight.res");
//  likep->workingweight.prettyPrint(out);
  // TEST

  if (workingresponsep.rows() != data.rows())
    {
    make_pointerindex();
    }

  unsigned int i,j;

  int size = posbeg.size();

  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();
  double * work_Wsum = Wsum.getV();

  double * * work_workingweightp = workingweightp.getV();


  if (intvar.rows() == data.rows())
    {
    double * work_intvar2=intvar2.getV();

    if (likep->wtype==wweightsnochange_one)
      {
      for (i=0;i<size;i++,work_Wsum++,itbeg++,itend++)
        {
        *work_Wsum=0;
        if (*itbeg != -1)
          {
          for (j=*itbeg;j<=*itend;j++,work_intvar2++)
            *work_Wsum +=  *work_intvar2 ;
          }
        }
      }
    else
      {

      for (i=0;i<size;i++,work_Wsum++,itbeg++,itend++)
        {
        *work_Wsum=0;
        if (*itbeg != -1)
          {
          for (j=*itbeg;j<=*itend;j++,work_workingweightp++,work_intvar2++)
            *work_Wsum +=  (*work_intvar2) * (**work_workingweightp);
          }
        }

      }

    }
  else
    {

    if (likep->wtype==wweightsnochange_one)
      {

      for (i=0;i<size;i++,work_Wsum++,itbeg++,itend++)
        {
        if (*itbeg != -1)
          {
          *work_Wsum = *itend - *itbeg + 1;
          }
        else
          *work_Wsum=0;
        }

      }
    else
      {

      for (i=0;i<size;i++,work_Wsum++,itbeg++,itend++)
        {
        *work_Wsum=0;
        if (*itbeg != -1)
          {
          for (j=*itbeg;j<=*itend;j++,work_workingweightp++)
            *work_Wsum += *(*work_workingweightp);
          }
        }

      }

    }

  // TEST
//  ofstream out2("c:\\bayesx\\test\\results\\Wsum.res");
//  Wsum.prettyPrint(out2);
  // TEST



// VARIANTE 1

  if (ZoutTZout_d.size() <= 1)
    {

    for (i=0;i<int(nrpar);i++)
      {
      for (j=0;j<ZoutT[i].size();j++)
        {
        ZoutTZout_d.push_back(pow(ZoutT[i][j],2));
        Wsump_d.push_back(index_ZoutT[i][j]);
        }

      }

    }

  vector<double>::iterator diag = XWX.getDiagIterator();
  vector<double>::iterator ZoutTZout_d_p = ZoutTZout_d.begin();
  vector<int>::iterator Wsump_d_p = Wsump_d.begin();
  int s;

  for (i=0;i<int(nrpar);i++,++diag)
    {
    *diag=0;
    s=ZoutT[i].size();
    for (j=0;j<s;j++,++ZoutTZout_d_p,++Wsump_d_p)
      {
      *diag += *ZoutTZout_d_p  * Wsum(*Wsump_d_p,0);
      }

    }


  vector<double>::iterator env = XWX.getEnvIterator();
  vector<unsigned>::iterator xenv = XWX.getXenvIterator();
  unsigned start = *xenv;
  unsigned nrnnull;
  xenv++;


  if (ZoutTZout.size() <= 1)
    compute_ZoutTZout();

  vector<double>::iterator ZoutTZoutp = ZoutTZout.begin();
  vector<int>::iterator beg_ZoutTZoutp = beg_ZoutTZout.begin();
  vector<int>::iterator Wsumpp = Wsump.begin();
  int nr=0;

  unsigned k;
  int beg, end;

  for(i=0;i<int(nrpar);i++,++xenv)
    {
    nrnnull = *xenv-start;
    if (nrnnull > 0)
      {
      for (j=i-nrnnull;j<i;j++,++env)
        {
        beg = *beg_ZoutTZoutp;
        if (nr < beg_ZoutTZout.size()-1)
          {
          beg_ZoutTZoutp++;
          end = *beg_ZoutTZoutp-1;
          nr++;
          }
        else
          {
          end =ZoutTZout.size()-1;

          }

        *env = 0;
        for (k=beg;k<=end;k++,++ZoutTZoutp,++Wsumpp)
          {
          *env += *ZoutTZoutp * Wsum(*Wsumpp,0);
          }

        }

      }
    start = *xenv;
    }


//  ofstream out("c:\\bayesx\\testh\\results\\XWX.res");
//  XWX.print4(out);





// VARIANTE 2 (alt)
/*
  vector<double>::iterator diag = XWX.getDiagIterator();
  double help;
  int ip;

  for (i=0;i<int(nrpar);i++,++diag)
    {
    *diag=0;

    for (j=0;j<ZoutT[i].size();j++)
      {
      help=pow(ZoutT[i][j],2);
      ip = index_ZoutT[i][j];
      *diag += help*Wsum(ip,0);
      }

    }


  vector<double>::iterator env = XWX.getEnvIterator();
  vector<unsigned>::iterator xenv = XWX.getXenvIterator();
  unsigned start = *xenv;
  unsigned nrnnull;
  xenv++;

  for(i=0;i<int(nrpar);i++,++xenv)
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
*/
//  ofstream out("c:\\bayesx\\testh\\results\\XWXalt.res");
//  XWX.print4(out);


  XWX.setDecomposed(false);

  // TEST
  /*
  vector<double> env = XWX.getEnv();

  ofstream out2("c:\\bayesx\\test\\results\\XWXenv.res");
  for (j=0;j<env.size();j++)
    out2 << env[j] << endl;


  vector<unsigned> Xenv = XWX.getXenv();

  ofstream out3("c:\\bayesx\\test\\results\\XWX_Xenv.res");
  for (j=0;j<Xenv.size();j++)
    out3 << Xenv[j] << endl;



  ofstream out3("c:\\bayesx\\test\\results\\XWX.res");
  XWX.print2(out3);

  unsigned k;
  datamatrix Zoutm(data.rows(),nrpar,0);
  for (i=0;i<posbeg.size();i++)
    {
    for (j=posbeg[i];j<=posend[i];j++)
      {
      for(k=0;k<Zout.cols();k++)
        Zoutm(j,index_Zout(i,k)) =
        sqrt(likep->workingweight(index_data(j,0),0))*Zout(i,k);
      }

    }

  datamatrix ZtZ = Zoutm.transposed()*Zoutm;
  ofstream out5("c:\\bayesx\\test\\results\\XWXmat.res");
  ZtZ.prettyPrint(out5);


  Zoutm = datamatrix(data.rows(),nrpar,0);
  for (i=0;i<posbeg.size();i++)
    {
    for (j=posbeg[i];j<=posend[i];j++)
      {
      for(k=0;k<Zout.cols();k++)
        Zoutm(j,index_Zout(i,k)) =
        Zout(i,k);
      }

    }


//  ofstream out7("c:\\bayesx\\test\\results\\Zoutm.res");
//  Zoutm.prettyPrint(out7);
  */

  // TEST

  }


void DESIGN::compute_ZoutTZout(void)
  {
  vector<unsigned>::iterator xenv = XWX.getXenvIterator();
  unsigned start = *xenv;
  unsigned nrnnull;
  xenv++;

  unsigned i,j;
  for(i=0;i<int(nrpar);i++,++xenv)
    {
    nrnnull = *xenv-start;
    if (nrnnull > 0)
      {
      for (j=i-nrnnull;j<i;j++)
        {
        compute_ZoutTZout(i,j);
        }

      }
    start = *xenv;
    }

/*
  ofstream out("c:\\bayesx\\testh\\results\\beg_ZoutTZout.res");
  for(i=0;i<beg_ZoutTZout.size();i++)
    out << beg_ZoutTZout[i] << endl;

  out.close();

  ofstream out2("c:\\bayesx\\testh\\results\\ZoutTZout.res");
  for(i=0;i<ZoutTZout.size();i++)
    out2 << ZoutTZout[i] << endl;

  out2.close();
*/

  }


void DESIGN::compute_ZoutTZout(unsigned & i, unsigned & j)
  {

  beg_ZoutTZout.push_back(ZoutTZout.size());

  unsigned k_i=0;
  unsigned k_j=0;
  int pos_i;
  int pos_j;

  while (k_i<ZoutT[i].size() && k_j < ZoutT[j].size())
    {

    pos_i = index_ZoutT[i][k_i];
    pos_j = index_ZoutT[j][k_j];

    if (pos_j > pos_i)
      {
      k_i++;
      }
    else if (pos_j < pos_i)
      {
      k_j++;
      }
    else  // equal
      {

      ZoutTZout.push_back(ZoutT[i][k_i]* ZoutT[j][k_j]);
      Wsump.push_back(pos_i);

      k_i++;
      k_j++;

      }

    }


  }


double DESIGN::compute_ZtZ(unsigned & i, unsigned & j)
  {

  unsigned k_i=0;
  unsigned k_j=0;
  int pos_i;
  int pos_j;

  double result = 0;

  while (k_i<ZoutT[i].size() && k_j < ZoutT[j].size())
    {

    pos_i = index_ZoutT[i][k_i];
    pos_j = index_ZoutT[j][k_j];

    if (pos_j > pos_i)
      {
      k_i++;
      }
    else if (pos_j < pos_i)
      {
      k_j++;
      }
    else  // equal
      {

      result += Wsum(pos_i,0) * ZoutT[i][k_i]* ZoutT[j][k_j];

      k_i++;
      k_j++;

      }

    }

  return result;
  }


void DESIGN::compute_XtransposedWres(datamatrix & partres, double l)
  {

  unsigned i,j;

  if (XWresdeclared == false)
    {
    XWres = datamatrix(nrpar,1);
    XWresdeclared = true;
    }


  if (ZoutT.size() != nrpar)
    compute_Zout_transposed();

  if (consecutive_ZoutT == -1)
    {
    bool c = check_ZoutT_consecutive();
    consecutive_ZoutT = c;
    }


  // TEST
  /*
  unsigned k;
  datamatrix Zoutm(Zout.rows(),nrpar,0);
  for (i=0;i<posbeg.size();i++)
    {
    for(k=0;k<Zout.cols();k++)
      Zoutm(i,index_Zout(i,k)) = Zout(i,k);
    }

  ofstream out4("c:\\bayesx\\test\\results\\Zoutm.res");
  Zoutm.prettyPrint(out4);


  datamatrix Ztres =  Zoutm.transposed()*partres;

  ofstream out3("c:\\bayesx\\test\\results\\XWres_alt.res");
  Ztres.prettyPrint(out3);
  */
  // END TEST


  double * workXWres = XWres.getV();
  unsigned size;
    vector<double>::iterator wZoutT;

  if (consecutive_ZoutT == 0)
    {

    vector<int>::iterator wZoutT_index;

    for(i=0;i<nrpar;i++,workXWres++)
      {
      *workXWres=0;
      wZoutT = ZoutT[i].begin();
      wZoutT_index = index_ZoutT[i].begin();
      size = ZoutT[i].size();
      for (j=0;j<size;j++,++wZoutT,++wZoutT_index)
        {
//      *workXWres += ZoutT[i][j]* partres(index_ZoutT[i][j],0);
        *workXWres+= (*wZoutT)* partres(*wZoutT_index,0);
        }
      }
    }
  else
    {
    double * wpartres;

    for(i=0;i<nrpar;i++,workXWres++)
      {
      *workXWres=0;
      wZoutT = ZoutT[i].begin();
      size = ZoutT[i].size();
      wpartres = partres.getV()+index_ZoutT[i][0];
      for (j=0;j<size;j++,++wZoutT,wpartres++)
        {
//      *workXWres += ZoutT[i][j]* partres(index_ZoutT[i][j],0);
        *workXWres+= (*wZoutT)* (*wpartres);
        }
      }

    }

  XWres_p = &XWres;  

  // TEST
//  ofstream out("c:\\bayesx\\test\\results\\XWres.res");
//  XWres.prettyPrint(out);
  // TEST


  }


void DESIGN::compute_effect(datamatrix & effect,datamatrix & f,
                      effecttype2 et)
  {

  // TEST
  /*
  ofstream out0("c:\\bayesx\\test\\results\\fueb.res");
  f.prettyPrint(out0);
  */
  // TEST

  int i,j;

  if (effect.rows() != data.rows())
    effect = datamatrix(data.rows(),1,0);

  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();

  double * workf = f.getV();
  double * workintvar = intvar.getV();

  int * workindex = index_data.getV();

  int size = posbeg.size();

  if (et==Function)
    {
    for (i=0;i<size;i++,++itbeg,++itend,workf++)
      {
      if (*itbeg != -1)
        {
        for (j=*itbeg;j<=*itend;j++,workindex++)
          effect(*workindex,0) = *workf;
        }
      }
    }
  else if (et==Varcoefftotal)
    {
    for (i=0;i<size;i++,++itbeg,++itend,workf++)
      {
      if (*itbeg != -1)
        {
        for (j=*itbeg;j<=*itend;j++,workindex++,workintvar++)
          effect(*workindex,0) = *workintvar * (*workf);
        }
      }
    }

  // TEST
  /*
  ofstream out("c:\\bayesx\\test\\results\\effect.res");
  for (i=0;i<effect.rows();i++)
    out << effect(i,0) << "  " << endl;

  ofstream out2("c:\\bayesx\\test\\results\\data.res");
  data.prettyPrint(out2);

  ofstream out3("c:\\bayesx\\test\\results\\index.res");
  index_data.prettyPrint(out3);
  */
  // TEST


  }


void DESIGN::set_intvar(datamatrix & iv,double add)
  {

  unsigned j;

  double * workintvar = intvar.getV();
  double * workintvar2 = intvar2.getV();
  int * workindex = index_data.getV();
  for (j=0;j<iv.rows();j++,workintvar++,workindex++,workintvar2++)
    {
    *workintvar = iv(*workindex,0)+add;
    *workintvar2 = pow(*workintvar,2);
    }

  }


void DESIGN::compute_f(datamatrix & beta,datamatrix & betalin,
                       datamatrix & f, datamatrix & ftot)
  {

  if (identity)
    {
    f.assign(beta);
    }
  else
    {

    if (consecutive == -1)
      {
      bool c = check_Zout_consecutive();
      consecutive = c;
      }

    unsigned i,j;

    unsigned rows;
    unsigned cols;
    rows = Zout.rows();
    cols = Zout.cols();

    double * workZ = Zout.getV();

    double * workf = f.getV();

    if (consecutive == 0)
      {

      int * workindex = index_Zout.getV();

      for(i=0;i<rows;i++,workf++)
        {
        *workf = 0;
        for(j=0;j<cols;j++,workindex++,workZ++)
          {
          *workf += (*workZ) * beta(*workindex,0);
          }
        }
      }
    else
      {

      double * workbeta;
      for(i=0;i<rows;i++,workf++)
        {
        *workf = 0;
        workbeta = beta.getV()+index_Zout(i,0);
        for(j=0;j<cols;j++,workZ++,workbeta++)
          {
          *workf += (*workZ) * (*workbeta);
          }
        }

      }

    }

  // TEST
  
  /*
  ofstream out("c:\\bayesx\\testh\\results\\f.res");
  f.prettyPrint(out);

  ofstream out2("c:\\bayesx\\testh\\results\\beta.res");
  beta.prettyPrint(out2);
  */

  /*
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

  if (position_lin!=-1)
    {
    ftot.mult(designlinear,betalin);
    ftot.plus(f);
    }


  }

void DESIGN::compute_precision(double l)
  {

  }


void DESIGN::compute_partres(datamatrix & res, datamatrix & f)
  {

  int i,j;
  int size = posbeg.size();
  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();
//  int * workindex = index_data.getV();
  double * workf = f.getV();

  double * workres = res.getV();

  if (workingresponsep.rows() != data.rows())
    {
    make_pointerindex();
    }

  double * * work_responsep = workingresponsep.getV();
  double * * work_workingweightp = workingweightp.getV();

  double * * worklinp;
  if (likep->linpred_current==1)
    worklinp = linpredp1.getV();
  else
    worklinp = linpredp2.getV();

  if (intvar.rows()==data.rows())   // varying coefficient
    {

    double * workintvar = intvar.getV();

    if (likep->wtype==wweightsnochange_one)
      {
      for (i=0;i<size;i++,++itbeg,++itend,workf++,workres++)
        {
        *workres = 0;
        if (*itbeg != -1)
          {
          for (j=*itbeg;j<=*itend;j++,work_responsep++,worklinp++,workintvar++)
            {
            *workres += (*workintvar) * (*(*work_responsep) - (*(*worklinp)) +
            (*workintvar) * (*workf));
            }
          }
        }
      }
    else
      {
      for (i=0;i<size;i++,++itbeg,++itend,workf++,workres++)
        {
        *workres = 0;
        if (*itbeg != -1)
          {
          for (j=*itbeg;j<=*itend;j++,work_responsep++,
               work_workingweightp++,worklinp++,workintvar++)
            {
            *workres += *(*work_workingweightp) * (*workintvar) *
            (*(*work_responsep) - (*(*worklinp)) + (*workintvar) * (*workf));
            }
          }
        }
      }

    }
  else                              // additive
    {

    if (likep->wtype==wweightsnochange_one)
      {

      for (i=0;i<size;i++,++itbeg,++itend,workf++,workres++)
        {
        *workres = 0;
        if (*itbeg != -1)
          {
          for (j=*itbeg;j<=*itend;j++,work_responsep++,worklinp++)
            {
            *workres += *(*work_responsep) - (*(*worklinp)) + *workf;

            }
          }
        }

      }
    else
      {
      for (i=0;i<size;i++,++itbeg,++itend,workf++,workres++)
        {
        *workres = 0;
        if (*itbeg != -1)
          {
          for (j=*itbeg;j<=*itend;j++,work_responsep++,
               work_workingweightp++,worklinp++)
            {
            *workres += *(*work_workingweightp) *
            (*(*work_responsep) - (*(*worklinp)) + *workf);
            }
          }
        }
      }

    }

  // TEST
  // ofstream out("c:\\bayesx\\test\\results\\tildey.res");
  // (likep->workingresponse).prettyPrint(out);
  // TEST

  }


void DESIGN::compute_partres(int begin,int end,double & res, double & f)
  {

  int j;

  if (workingresponsep.rows() != data.rows())
    {
    make_pointerindex();
    }

  double * * work_responsep = workingresponsep.getV()+begin;
  double * * work_workingweightp = workingweightp.getV()+begin;

  double * * worklinp;
  if (likep->linpred_current==1)
    worklinp = linpredp1.getV()+begin;
  else
    worklinp = linpredp2.getV()+begin;

  if (intvar.rows()==data.rows())   // varying coefficient
    {

    double * workintvar = intvar.getV()+begin;


    if (likep->wtype==wweightsnochange_one)
      {
      res = 0;
      if (begin != -1)
        {
        for (j=begin;j<=end;j++,work_responsep++,worklinp++,workintvar++)
          {
          res += (*workintvar) * (*(*work_responsep) - (*(*worklinp)) +
          (*workintvar) * f);
          }
        }
      }
    else
      {
      res = 0;
      if (begin != -1)
        {
        for (j=begin;j<=end;j++,work_responsep++,
             work_workingweightp++,worklinp++,workintvar++)
          {
          res += *(*work_workingweightp) * (*workintvar) *
          (*(*work_responsep) - (*(*worklinp)) + (*workintvar) * f);
          }
        }
      }
    }
  else                              // additive
    {

    if (likep->wtype==wweightsnochange_one)
      {
      res = 0;
      if (begin != -1)
        {
        for (j=begin;j<=end;j++,work_responsep++,worklinp++)
          {
          res += *(*work_responsep) - (*(*worklinp)) + f;
          }
        }
      }
    else
      {
      res = 0;
      if (begin != -1)
        {
        for (j=begin;j<=end;j++,work_responsep++,
               work_workingweightp++,worklinp++)
          {
          res += *(*work_workingweightp) *
          (*(*work_responsep) - (*(*worklinp)) + f);
          }
        }
      }

   }

  // TEST
  //    ofstream out("c:\\bayesx\\test\\results\\tildey.res");
  //    (likep->workingresponse).prettyPrint(out);
  // TEST

  }



void DESIGN::compute_meaneffect(DISTR * level1_likep,double & meaneffect,
                                datamatrix & beta,datamatrix & meaneffectbeta,
                                bool computemeaneffect)

  {

  level1_likep->meaneffect -= meaneffect;

  meaneffect = meaneffectintvar* beta(meaneffectnr,0);

  if (computemeaneffect==true)
    {

    unsigned i;
    double * betap = beta.getV();
    double * meffectp = meaneffectbeta.getV();
    double l;
    if (intvar.rows() == data.rows())
      {
      for(i=0;i<beta.rows();i++,meffectp++,betap++)
        {
        l=level1_likep->meaneffect+meaneffectintvar*(*betap);
        level1_likep->compute_mu(&l,meffectp);
        }
      }
    else
      {
      for(i=0;i<beta.rows();i++,meffectp++,betap++)
        {
        l=level1_likep->meaneffect+(*betap);
        level1_likep->compute_mu(&l,meffectp);
        }
      }

    }

  level1_likep->meaneffect += meaneffect;

  }



void DESIGN::update_linpred(datamatrix & f)
  {
  int i,j;

  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();

  double * workf = f.getV();
  double * workintvar = intvar.getV();

  double * * linpredp;
  if (likep->linpred_current==1)
    linpredp = linpredp1.getV();
  else
    linpredp = linpredp2.getV();

  int size = posbeg.size();


  if (intvar.rows()==data.rows())   // varying coefficient
    {

    for (i=0;i<size;i++,++itbeg,++itend,workf++)
      {
      if (*itbeg != -1)
        {
        for (j=*itbeg;j<=*itend;j++,workintvar++,linpredp++)
          *(*linpredp) += (*workintvar) * (*workf);
        }
      }

    }
  else                              // additive
    {
    for (i=0;i<size;i++,++itbeg,++itend,workf++)
      {
      if (*itbeg != -1)
        {
        for (j=*itbeg;j<=*itend;j++,linpredp++)
          *(*linpredp) += *workf;
        }
      }
    }


  // TEST
  /*
  ofstream out3("c:\\bayesx\\test\\results\\lin.res");
  likep->linearpred1.prettyPrint(out3);
  */
  // TEST


  }



void DESIGN::read_options(vector<ST::string> & op,vector<ST::string> & vn)
  {

  if (op[7] == "false")
    center = true;
  else
    center = false;

  datanames = vn;

  }


void DESIGN::outoptions(GENERAL_OPTIONS * op)
  {

  }


} // end: namespace MCMC




