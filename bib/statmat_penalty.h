
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

statmatrix<double> __EXPORT_TYPE K2dim_pspline_rw2(const unsigned & nknots, const unsigned & ox, const unsigned & oy);

statmatrix<double> __EXPORT_TYPE K2dim_pspline_biharmonic(const unsigned & nknots);
}

#endif
