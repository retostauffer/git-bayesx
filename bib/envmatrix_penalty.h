#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined(ENVMATRIX_PENALTY)

#define ENVMATRIX_PENALTY

#include "envmatrix.h"

  //----------------------------------------------------------------------------
  //----------------- Functions for computing penalty matrices -----------------
  //----------------------------------------------------------------------------

  // FUNCTION Kmrfenv
  // TASK: returns the penalty matrix for MRF with characteristics stored in map

  envmatrix<double> Kmrfenv(const MAP::map & m);

  // FUNCTION: Krw1env
  // TASK: returns the penalty matrix for a first order random walk

  envmatrix<double> Krw1env(const vector<double> & weight);

  // FUNCTION: Krw2env
  // TASK: returns the penalty matrix for a second order random walk

  envmatrix<double> Krw2env(const vector<double> & weight);

  // FUNCTION: Kseasonenv
  // TASK: returns the penalty matrix for a sesonal component with period 'per'

  envmatrix<double> Kseasonenv(const unsigned & per,const unsigned & s);

  // FUNCTION: Krw0env
  // TASK: returns the identity matrix (penalty for a random effect)

  envmatrix<double> Krw0env(const unsigned & nrpar);


#endif
