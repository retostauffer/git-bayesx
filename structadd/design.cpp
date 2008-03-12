
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

DESIGN::DESIGN(const datamatrix & dm)
  {
  data = dm;

  index_data = statmatrix<int>(data.rows(),1);
  index_data.indexinit();
  dm.indexsort(index_data,0,dm.rows()-1,0,0);


  }

// COPY CONSTRUCTOR

DESIGN::DESIGN(const DESIGN & m)
  {
  data = m.data;
  index_data = m.index_data;
  Z=m.Z;
  Zout = m.Zout;
  index_Zout=m.index_Zout;
  type=m.type;
  K = m.K;
  XWX = m.XWX;  
  }

// OVERLOADED ASSIGNMENT OPERATOR

const DESIGN & DESIGN::operator=(const DESIGN & m)
  {
  if (this == &m)
    return *this;
  data = m.data;
  index_data = m.index_data;
  Z=m.Z;
  Zout = m.Zout;  
  index_Zout=m.index_Zout;
  type=m.type;
  K = m.K;
  XWX = m.XWX;    
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

  }


void DESIGN::compute_XtransposedW(void)
  {


  }


void DESIGN::compute_XtransposedWX(void)
  {

  }

//------------------------------------------------------------------------------
//-------------- CLASS: DESIGN_mrf implementation of member functions ----------
//------------------------------------------------------------------------------


DESIGN_mrf::DESIGN_mrf(void) : DESIGN()
  {

  }

  // CONSTRUCTOR 1
  // Spatial covariates

DESIGN_mrf::DESIGN_mrf(const datamatrix & dm,const MAP::map & m)
                      : DESIGN(dm)
  {
  ma = m;
  type = mrf;

  compute_design();  

  compute_penalty();

  }

  // COPY CONSTRUCTOR

DESIGN_mrf::DESIGN_mrf(const DESIGN_mrf & m)
    : DESIGN(DESIGN(m))
  {
  ma = m.ma;
  posbeg = m.posbeg;
  posend = m.posend;
  effectvalues = m.effectvalues;
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const DESIGN_mrf & DESIGN_mrf::operator=(const DESIGN_mrf & m)
  {
  if (this == &m)
    return *this;
  DESIGN::operator=(DESIGN(m));
  ma = m.ma;
  posbeg = m.posbeg;
  posend = m.posend;
  effectvalues = m.effectvalues;
  return *this;
  }


void DESIGN_mrf::compute_design(void)
  {

  Zout = datamatrix(ma.get_nrregions(),1,1);
  index_Zout = statmatrix<int>(Zout.rows(),1);

  ma.compute_reg(data,posbeg,posend,effectvalues,index_data);

  }


void DESIGN_mrf::compute_penalty(void)
  {
  if (type==mrf)
    K = Kmrfenv(ma);

  }


void DESIGN::compute_XtransposedWX(void)
  {

  }


} // end: namespace MCMC



