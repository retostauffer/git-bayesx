// Date: 4.12.99

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif


#if !defined (MCMCnonpgaussian_INCLUDED)

#define MCMCnonpgaussian_INCLUDED

#include<mcmc_nonpbasis.h>
#include<statmat_penalty.h>

namespace MCMC
{
} // end: namespace MCMC

#endif
