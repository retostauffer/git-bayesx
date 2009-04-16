#if defined(BORLAND_OUTPUT_WINDOW)
#include "StatResults.h"
#include "statwinframe.h"

#endif


#include "MASTER_obj.h"



namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------------- CLASS: MASTER_OBJ ------------------------------
//------------------------------------------------------------------------------


MASTER_OBJ::MASTER_OBJ(void)
  {
  }


MASTER_OBJ::MASTER_OBJ(const MASTER_OBJ & o)
  {
  level1_likep = o.level1_likep;
  const_pointer = o.const_pointer;
  constposition = o.constposition;

  }


const MASTER_OBJ & MASTER_OBJ::operator=(const MASTER_OBJ & o)
  {
  if (this == &o)
    return *this;
  level1_likep = o.level1_likep;
  const_pointer = o.const_pointer;
  constposition = o.constposition;
  return *this;
  }



} // end: namespace MCMC



