
#include "design_pspline.h"


namespace MCMC
{


//------------------------------------------------------------------------------
//-------------- CLASS: DESIGN_mrf implementation of member functions ----------
//------------------------------------------------------------------------------


DESIGN_pspline::DESIGN_pspline(void) : DESIGN()
  {

  }

  // CONSTRUCTOR

DESIGN_pspline::DESIGN_pspline(const datamatrix & dm,const datamatrix & iv,
                       DISTR * dp)
                      : DESIGN(dp)
  {

  degree = 3;
  nrknots=20;
  type=RW2;

  init_data(dm,iv);

  weightK = vector<double>(nrpar,1);
  compute_penalty();

  datamatrix  help(Zout.rows(),1,1);
  compute_XtransposedWX_XtransposedWres(help, 1);

  compute_precision(1.0);

  }


  // COPY CONSTRUCTOR

DESIGN_pspline::DESIGN_pspline(const DESIGN_pspline & m)
    : DESIGN(DESIGN(m))
  {
  knot = m.knot;
  nrknots = m.nrknots;
  degree = m.degree;
  weightK = m.weightK;
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const DESIGN_pspline & DESIGN_pspline::operator=(const DESIGN_pspline & m)
  {
  if (this == &m)
    return *this;
  DESIGN::operator=(DESIGN(m));
  knot = m.knot;
  nrknots = m.nrknots;
  degree = m.degree;
  weightK = m.weightK;
  return *this;
  }


void DESIGN_pspline::init_data(const datamatrix & dm,const datamatrix & iv)
  {

  // TASK: sorts the data such that the precision has minimum envelope
  //       computes index_data, posbeg, posend
  //       computes effectvalues
  //       computes Zout, index_Zout
  //       computes nrpar
  //       initializes datanames

  make_index(dm,iv);

  nrpar = nrknots-1+degree;

  datanames.push_back("X_1");

  make_Bspline();

  compute_Zout_transposed();


  }



void DESIGN_pspline::make_Bspline(void)
  {

  unsigned i,j,k;
  double value;

  datamatrix help;

// berechne x_min, x_max
  double min;
  min = data(0,0);

  double max;
  max = data(data.rows()-1,0);

  double dist = max-min;

  min -= 0.01*dist;
  max += 0.01*dist;


// Knoten berechnen

  dist = (max - min)/(nrknots-1);
  knot.push_back(min - degree*dist);
  for(i=1;i<nrknots+2*degree;i++)
    knot.push_back(knot[i-1] + dist);


  Zout = datamatrix(posbeg.size(),degree+1,0.0);
  index_Zout = statmatrix<int>(posbeg.size(),degree+1,0.0);
  double * work = Zout.getV();
  int * work_index = index_Zout.getV();

  help = datamatrix(nrpar,1,0.0);



  for (i=0;i<posbeg.size();i++)
    {
    value = data(posbeg[i],0);
    j=0;
    while(knot[degree+j+1] <= value)
      j++;
    help.assign(bspline(value));
    for (k=0;k<Zout.cols();k++,work++,work_index++)
      {
      *work = help(k+j,0);
      *work_index = j+k;
      }

    }

  consecutive = 1;


  // TEST
  /*
  bool t = check_Zout_consecutive();

  ofstream out("c:\\bayesx\\test\\results\\Zout.res");
  Zout.prettyPrint(out);

  ofstream out2("c:\\bayesx\\test\\results\\index_Zout.res");
  index_Zout.prettyPrint(out2);

  datamatrix Zoutm(data.rows(),nrpar,0);
  for (i=0;i<posbeg.size();i++)
    {
    for (j=posbeg[i];j<=posend[i];j++)
      {
      for(k=0;k<Zout.cols();k++)
        Zoutm(j,index_Zout(i,k)) = sqrt(likep->workingweight(index_data(j,0),0))*Zout(i,k)*intvar(j,0);
      }

    }

  datamatrix h(data.rows(),1,1);
  datamatrix Zteins =  Zoutm.transposed()*h;

  datamatrix ZtZ = Zoutm.transposed()*Zoutm;

  ofstream out3("c:\\bayesx\\test\\results\\Zteins.res");
  Zteins.prettyPrint(out3);

  ofstream out4("c:\\bayesx\\test\\results\\ZtZ.res");
  ZtZ.prettyPrint(out4);
  */
  // TEST


  }



datamatrix DESIGN_pspline::bspline(const double & x)
  {
// nach Hämmerlin/Hoffmann
  datamatrix b(nrpar,1,0.0);
  datamatrix help(nrpar+degree,1,0.0);
  unsigned j;
  double * bwork;
  double * helpwork;

// Grad 0

  for(j=0;j<nrpar;j++)
    if( knot[j]<=x && x<knot[j+1])
      b(j,0) = 1.0;

  for(unsigned l=1;l<=degree;l++)
    {
    bwork = b.getV();
    helpwork = help.getV();
    for(j=0;j<nrpar;j++,++helpwork,++bwork)
//      help(j,0) = b(j,0);
      *helpwork = *bwork;
    bwork = b.getV();
    helpwork = help.getV();
    for(j=0;j<nrpar;j++,++helpwork,++bwork)
      {
//      b(j,0) = (x-knot[j])*help(j,0)/(knot[j+l]-knot[j])
//                  + (knot[j+l+1]-x)*help(j+1,0)/(knot[j+l+1]-knot[j+1]);
      *bwork = (x-knot[j])**helpwork/(knot[j+l]-knot[j])
                  + (knot[j+l+1]-x)**(helpwork+1)/(knot[j+l+1]-knot[j+1]);

      }
    }

  return b;

  }


void DESIGN_pspline::compute_penalty(void)
  {
  if (type==RW1)
    {
    K = Krw1env(weightK);
    rankK = nrpar-1;
    }
  else if (type==RW2)
    {
    K = Krw2env(nrpar);
//    ofstream out("c:\\bayesx\\test\\results\\K.res");
//    K.print2(out);
    rankK = nrpar-2;
    }
  else if (type==RW3)
    {
    K = Krw3env(nrpar);
    rankK = nrpar-3;
    }
  }


void DESIGN_pspline::compute_XtransposedWX_XtransposedWres(
datamatrix & partres, double l)
  {

  if (XWXdeclared==false)
    {
    XWX = envmatdouble(bandmatdouble(nrpar,degree,0));
    Wsum = datamatrix(posbeg.size(),1,0);
    XWXdeclared  = true;
    }

  DESIGN::compute_XtransposedWX_XtransposedWres(partres, l);


  }


void DESIGN_pspline::compute_XtransposedWX(datamatrix & partres)
  {
  DESIGN::compute_XtransposedWX(partres);
  }


void DESIGN_pspline::compute_XtransposedWres(datamatrix & partres, double l)
  {
  DESIGN::compute_XtransposedWres(partres, l);
  }


void DESIGN_pspline::compute_precision(double l)
  {
  if (precisiondeclared==false)
    {
    precision = envmatdouble(0.0,nrpar,degree>2?degree:2);
    precisiondeclared = true;
    }

  precision.addto(XWX,K,1.0,l);

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



