
#include "design_mrf.h"
#include "clstring.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//-------------- CLASS: DESIGN_mrf implementation of member functions ----------
//------------------------------------------------------------------------------


void DESIGN_mrf::read_options(vector<ST::string> & op,vector<ST::string> & vn)
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
  13      samplemult
  14      constraints
  15      round
  16      centermethod
  */


  if (op[7] == "false")   //nocenter==false, i.e. center
    center = true;
  else
    center = false;

  if (op[16]=="mean" || op[16] == "nullspace")
    centermethod = cmean;
  else if (op[16] == "meansimple")
    centermethod = meansimple;
  else if (op[16] == "meaninvvar")
    centermethod = cmeaninvvar;
  else if (op[16] == "meanintegral")
    centermethod = cmeanintegral;
  else if (op[16] == "meanf")
    centermethod = cmeanf;

  datanames = vn;

  }


DESIGN_mrf::DESIGN_mrf(void) : DESIGN()
  {

  }

  // CONSTRUCTOR 1
  // Spatial covariates

DESIGN_mrf::DESIGN_mrf(const datamatrix & dm,const datamatrix & iv,
                       DISTR * dp,FC_linear * fcl,
                       const MAP::map & m,vector<ST::string> & op,
                       vector<ST::string> & vn)
                      : DESIGN(dp,fcl)
  {

  read_options(op,vn);

  ma = m;
  type = Mrf;

  init_data(dm,iv);

  compute_penalty();

  datamatrix  help(Zout.rows(),1,1);
  compute_XtransposedWX_XtransposedWres(help,1);

  compute_precision(1.0);

  compute_basisNull();

  identity=true;
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


void DESIGN_mrf::init_data(const datamatrix & dm, const datamatrix & iv)
  {

  // FUNCTION: init_data
  // TASK: sorts the data such that the precision has minimum envelope
  //       computes index_data
  //       computes Zout, posbeg, posend
  //       computes nrpar
  //       computes effectvalues
  //       computes consecutive
  //       computes ZoutT, index_ZoutT

  nrpar = ma.get_nrregions();
  consecutive=true;

  Zout = datamatrix(nrpar,1,1);
  index_Zout = statmatrix<int>(Zout.rows(),1);
  index_Zout.indexinit();

  if (ma.get_bandsize() > 40)
    ma.reorderopt();

  ma.compute_reg(dm,posbeg,posend,effectvalues,index_data);

  make_data(dm,iv);

  meaneffectnr = compute_modecategorie();
  compute_meaneffectintvar();

  if (ma.get_errormessages().size() > 0)
    {
//  FEHLT!!
    }

  }



void DESIGN_mrf::compute_penalty(void)
  {
  if (type==Mrf)
    K = Kmrfenv(ma);
  rankK = ma.get_nrregions()-1;
  }


void DESIGN_mrf::compute_basisNull(void)
  {
  int i,j;

  basisNull = datamatrix(1,nrpar,1);

  if (centermethod==cmeanf || centermethod==cmeanintegral)
    {

    unsigned k;
    for (k=0;k<nrpar;k++)
      {
      if (posbeg[k] != -1)
        basisNull(0,k) = posend[k]-posbeg[k]+1;
      else
        basisNull(0,k) = 0;
      }

    }


  position_lin = -1;


  for(i=0;i<basisNull.rows();i++)
    {
    basisNullt.push_back(datamatrix(basisNull.cols(),1));
    for(j=0;j<basisNull.cols();j++)
      basisNullt[i](j,0) = basisNull(i,j);
    }


  // TEST
  /*
    ofstream out("c:\\bayesx\\test\\results\\data.res");
    data.prettyPrint(out);

    ofstream out2("c:\\bayesx\\test\\results\\designlin.res");
    designlinear.prettyPrint(out2);
   */
  // TEST

  }


void DESIGN_mrf::compute_XtransposedWX_XtransposedWres(datamatrix & partres,
                                                        double l)
  {

  compute_XtransposedWX();

  compute_XtransposedWres(partres,l);

  }


void DESIGN_mrf::compute_XtransposedWres(datamatrix & partres, double l)
  {
  if (XWresdeclared == false)
    {
    XWresdeclared = true;
    }

  XWres_p = &partres;

  }


void DESIGN_mrf::compute_XtransposedWX(void)
  {

  if (XWXdeclared == false)
    {
    XWX = envmatdouble(0,nrpar);
    XWXdeclared = true;
    }

  unsigned i;

  double * Wsump = Wsum.getV();
  vector<double>::iterator d = XWX.getDiagIterator();

  for (i=0;i<Wsum.rows();i++,++d,Wsump++)
    *Wsump = *d;

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


void DESIGN_mrf::outoptions(GENERAL_OPTIONS * op)
  {

  }

} // end: namespace MCMC



