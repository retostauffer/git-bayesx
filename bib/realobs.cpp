  
#include "first.h"

#include "realobs.h"

namespace realob
{

//------------------------------------------------------------------------------
//------------ CLASS realobs: implementation of member functions ---------------
//------------------------------------------------------------------------------


realobs sqrt(realobs & o)
  {
  if ( (o.value < 0) || (o.value == NA) )
	 return NA;
  else
	 return std::sqrt(o.value);
  }


realobs abs(realobs & o)
  {
  if (o.value == NA)
	 return NA;
  else
	 return fabs(o.value);
  }


realobs exp(realobs & o)
  {
  if (o.value == NA)
	 return NA;
  else
	 return std::exp(o.value);
  }


realobs cos(realobs & o)
  {
  if (o.value == NA)
	 return NA;
  else
	 return std::cos(o.value);
  }


realobs sin(realobs & o)
  {
  if (o.value == NA)
	 return NA;
  else
	 return std::sin(o.value);
  }


realobs log(realobs & o)
  {
  if  ( (o.value <= 0) || (o.value == NA) )
	 return NA;
  else
	 return std::log(o.value);
  }


realobs log10(realobs & o)
  {
  if  ( (o.value <= 0) || (o.value == NA) )
	 return NA;
  else
	 return std::log10(o.value);
  }


realobs pow(const realobs & o,const realobs & p)
  {
  if ((o.value == NA) || (p.value == NA))
	 return NA;
  else
	 return std::pow(o.value,p.value);
  }


realobs pow(realobs & o, double & p)
  {
  if ((o.value == NA) || (p == NA))
	 return NA;
  else
	 return std::pow(o.value,p);
  }


realobs pow(double o,realobs & p)
  {
  if ((o == NA) || (p.value == NA))
	 return NA;
  else
	 return std::pow(o,p.value);
  }

realobs floor(realobs & o)
  {
  if (o.value == NA)
	 return NA;
  else
	 return std::floor(o.value);
  }

} // end: namespace realobs

