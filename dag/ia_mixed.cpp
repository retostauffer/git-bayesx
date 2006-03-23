
#include "first.h"

#include "ia_mixed.h"
#include <algorithm>
#include <iterator>



namespace MCMC
{

	// DEFAULT CONSTRUCTOR:
	IA_MIXED::IA_MIXED(void): IA()
	{
	}
	


	// CONSTRUCTOR_1
	IA_MIXED::IA_MIXED(const datamatrix & d):IA(d)
	{
		
	}




	// CONSTRUCTOR_2
	// for interactions of order>2 (some day in future....) 
	IA_MIXED::IA_MIXED(unsigned order, const datamatrix & d):IA(order, d)
	{
	}


	// COPY CONSTRUCTOR
	IA_MIXED::IA_MIXED(const IA_MIXED & a) : IA (IA(a))
	{	
	}


	 // OVERLOADED ASSIGNMENT OPERATOR
	const IA_MIXED & IA_MIXED::operator=(const IA_MIXED & a)
	{
		if (this==&a)
		  return *this;

		// x = a.x;

		return *this;
	}





	


} //namespace MCMC