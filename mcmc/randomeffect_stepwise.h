

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (MCMCrandom_INCLUDED)

#define MCMCrandom_INCLUDED


#include<mcmc.h>
#include<fullcond.h>
#include<distribution.h>
#include<mcmc_nonpbasis.h>
#include<mcmc_nonp.h>


namespace MCMC
{
}   // end: namespace MCMC


#endif

