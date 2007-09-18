// functions that do not belong to one class
// but are used by several ones



#include "functions.h"

namespace func
{

	
// FUNCTION: accept
// TASK: returns true with probability ratio

bool accept (double ratio)
{
	double u =  randnumbers::uniform();
		
	if( log(u) > ratio )
		return false;
	else 
		return true ;
}
	

}  
