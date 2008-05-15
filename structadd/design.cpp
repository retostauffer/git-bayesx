
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
  responsep = statmatrix<double*> (data.rows(),1);
  for(i=0;i<data.rows();i++)
    responsep(i,0) = &(likep->response(index_data(i,0),0));

  workingweightp = statmatrix<double*> (data.rows(),1);
  for(i=0;i<data.rows();i++)
    workingweightp(i,0) = &(likep->workingweight(index_data(i,0),0));

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

  for(i=0;i<nrpar;i++)
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

void DESIGN::make_index(const datamatrix & dm,const datamatrix & iv)
  {

  unsigned j;

  index_data = statmatrix<int>(dm.rows(),1);
  index_data.indexinit();
  dm.indexsort(index_data,0,dm.rows()-1,0,0);

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

  for(j=0;j<posbeg.size();j++)
    effectvalues.push_back(ST::doubletostring(data(posbeg[j],0)));



  /*
  // TEST
  int ev = effectvalues.size();

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

  changingdesign = false;

  likep = lp;

  XWXdeclared = false;
  XWresdeclared = false;
  precisiondeclared=false;
  consecutive = -1;
  consecutive_ZoutT = -1;
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

  Zout = m.Zout;
  index_Zout=m.index_Zout;
  posbeg = m.posbeg;
  posend = m.posend;
  consecutive = m.consecutive;
  consecutive_ZoutT = m.consecutive_ZoutT;

  ZoutT = m.ZoutT;
  index_ZoutT = m.index_ZoutT;

  responsep = m.responsep;
  workingweightp = m.workingweightp;
  linpredp1 = m.linpredp1;
  linpredp2 = m.linpredp2;

  nrpar = m.nrpar;
  center = m.center;

  K = m.K;
  rankK = m.rankK;

  XWX = m.XWX;
  XWXdeclared = m.XWXdeclared;
  precision = m.precision;
  precisiondeclared = m.precisiondeclared;
  Wsum = m.Wsum;

  XWres = m.XWres;
  XWresdeclared = m.XWresdeclared;

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

  Zout = m.Zout;
  index_Zout=m.index_Zout;
  posbeg = m.posbeg;
  posend = m.posend;
  consecutive = m.consecutive;
  consecutive_ZoutT = m.consecutive_ZoutT;

  ZoutT = m.ZoutT;
  index_ZoutT = m.index_ZoutT;

  responsep = m.responsep;
  workingweightp = m.workingweightp;
  linpredp1 = m.linpredp1;
  linpredp2 = m.linpredp2;

  nrpar = m.nrpar;
  center = m.center;

  K = m.K;
  rankK = m.rankK;

  XWX = m.XWX;
  XWXdeclared = m.XWXdeclared;
  precision = m.precision;
  precisiondeclared = m.precisiondeclared;
  Wsum = m.Wsum;

  XWres = m.XWres;
  XWresdeclared = m.XWresdeclared;

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

  compute_XtransposedWX();
  compute_XtransposedWres(partres, l);

  }


void DESIGN::compute_XtransposedWX(void)
  {

  if (responsep.rows() != data.rows())
    {
    make_pointerindex();
    }

  unsigned i,j;

  int size = posbeg.size();

  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();
  double * work_Wsum = Wsum.getV();

  double * * work_workingweightp = workingweightp.getV();

  if (intvar.rows() == data.rows())
    {
    double * work_intvar2=intvar2.getV();
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


  vector<double>::iterator diag = XWX.getDiagIterator();
  double help;
  int ip;

  for (i=0;i<nrpar;i++,++diag)
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

//  unsigned envs = XWX.getXenv(nrpar);
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

  ofstream out("c:\\bayesx\\test\\results\\XWX.res");
  XWX.print2(out);
 */
  // TEST


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

  /*
  // TEST
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
  // END TEST
  */

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

  /*
  // TEST
  ofstream out("c:\\bayesx\\test\\results\\XWres.res");
  XWres.prettyPrint(out);
  // TEST
  */

  }


void DESIGN::compute_effect(datamatrix & effect,datamatrix & f,
                      effecttype et)
  {

  // TEST
  /*
  ofstream out0("c:\\bayesx\\test\\results\\fueb.res");
  f.prettyPrint(out0);
  */
  // TEST

  unsigned i,j;

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


void DESIGN::compute_f(datamatrix & beta,datamatrix & f)
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

  int i,j;
  int size = posbeg.size();
  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();
//  int * workindex = index_data.getV();
  double * workf = f.getV();

  double * workres = res.getV();

  if (responsep.rows() != data.rows())
    {
    make_pointerindex();
    }

  double * * work_responsep = responsep.getV();
  double * * work_workingweightp = workingweightp.getV();

  double * * worklinp;
  if (likep->linpred_current==1)
    worklinp = linpredp1.getV();
  else
    worklinp = linpredp2.getV();

  if (intvar.rows()==data.rows())   // varying coefficient
    {

    double * workintvar = intvar.getV();

    if ((likep->changingweight==true) ||
    ((likep->changingweight==false) && (likep->weights_one==false)))
      {
      for (i=0;i<size;i++,++itbeg,++itend,workf++,workres++)
        {
        *workres = 0;
        if (*itbeg != -1)
          {
          for (j=*itbeg;j<=*itend;j++,work_responsep++,
               work_workingweightp,worklinp++,workintvar++)
            {
            *workres += *(*work_workingweightp) * (*workintvar) *
            (*(*work_responsep) - (*(*worklinp)) + (*workintvar) * (*workf));
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
          for (j=*itbeg;j<=*itend;j++,work_responsep++,worklinp++,workintvar++)
            {
            *workres += (*workintvar) * (*(*work_responsep) - (*(*worklinp)) +
            (*workintvar) * (*workf));
            }
          }
        }
      }
    }
  else                              // additive
    {

    if ((likep->changingweight==true) ||
    ((likep->changingweight==false) && (likep->weights_one==false)))
      {
      for (i=0;i<size;i++,++itbeg,++itend,workf++,workres++)
        {
        *workres = 0;
        if (*itbeg != -1)
          {
          for (j=*itbeg;j<=*itend;j++,work_responsep++,
               work_workingweightp,worklinp++)
            {
            *workres += *(*work_workingweightp) *
            (*(*work_responsep) - (*(*worklinp)) + *workf);
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
          for (j=*itbeg;j<=*itend;j++,work_responsep++,worklinp++)
            {
            *workres += *(*work_responsep) - (*(*worklinp)) + *workf;

            }
          }
        }
      }
    }



  // TEST
  //  ofstream out("c:\\bayesx\\test\\results\\lin.res");
  //  (likep->linearpred1).prettyPrint(out);
  // TEST


  }


void DESIGN::update_linpred(datamatrix & f,bool add)
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

  if (add==true)
    {
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
    }
  else
    {

    if (intvar.rows()==data.rows())   // varying coefficient
      {

      for (i=0;i<size;i++,++itbeg,++itend,workf++)
        {
        if (*itbeg != -1)
          {
          for (j=*itbeg;j<=*itend;j++,workintvar++,linpredp++)
            *(*linpredp) -= (*workintvar) * (*workf);
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
            *(*linpredp) -= *workf;
          }
        }
      }

    }


  // TEST
/*
  ofstream out("c:\\bayesx\\test\\results\\lin.res");
  (*linpredp).prettyPrint(out);
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


} // end: namespace MCMC



