
#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif

#if !defined(STATMATRIXPENALTY_INCLUDED)

#define STATMATRIXPENALTY_INCLUDED

#include<statmat.h>
#include<map.h>

namespace STATMAT_PENALTY
{

statmatrix<double> __EXPORT_TYPE Kmrf(const MAP::map & m);

statmatrix<double> __EXPORT_TYPE K2dim_pspline(const unsigned & nknots);

}

#endif
