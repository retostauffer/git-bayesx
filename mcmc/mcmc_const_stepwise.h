// DATE 29.01.98

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif


#if !defined (MCMCconststepwise_INCLUDED)

#define MCMCconststepwise_INCLUDED

#include<mcmc.h>
#include<fullcond.h>
#include<mcmc_const.h>
#include<distribution.h>
#include<nbinomial.h>
#include<zip.h>


namespace MCMC
{
} // end: namespace MCMC

#endif
