
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE
#else
#define __EXPORT_TYPE __import
#endif

#if !defined(BANDMATRIX_PENALTY)

#define BANDMATRIX_PENALTY

#include "bandmat.h"


  //----------------------------------------------------------------------------
  //----------------- Functions for computing penalty matrices -----------------
  //----------------------------------------------------------------------------

  // FUNCTION: Krw1
  // TASK: returns the penalty matrix for a first order random walk

  bandmatdouble Krw1band(const vector<double> & weight);

  // FUNCTION: Krw2
  // TASK: returns the penalty matrix for a second order random walk

  bandmatdouble Krw2band(const vector<double> & weight);

  // FUNCTION: Kseason
  // TASK: returns the penalty matrix for a sesonal component with period 'per'

  bandmatdouble Kseasonband(const unsigned & per,const unsigned & s);

  // FUNCTION Kmrfband
  // TASK: returns the penalty matrix for MRF with characteristics stored in map

  bandmatdouble Kmrfband(const MAP::map & m);

  // FUNCTION: Kmrflinear
  // TASK: returns the penalty matrix for a columnwise first order random
  //       walk

  bandmatdouble Kmrflinearband(const unsigned & nr1,const unsigned & nr2);

#endif
