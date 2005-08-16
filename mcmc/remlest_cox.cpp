#if defined(BORLAND_OUTPUT_WINDOW)
#include "StatResults.h"
#include "statwinframe.h"
#endif

#include<remlest.h>
#include<remlest_cox.h>

//------------------------------------------------------------------------------
//----------------- AFT models with smooth error distribution ------------------
//------------------------------------------------------------------------------

  bool remlest::estimate_aft(datamatrix resp, const datamatrix & offset,
                    const datamatrix & weight)
  {
  unsigned i;

  outoptions();
  out("\n");

  for(i=0;i<fullcond.size();i++)
    fullcond[i]->outoptionsreml();

  out("\n");
  out("REML ESTIMATION STARTED\n",true);
  out("\n");

  bool stop = check_pause();
  if (stop)
    return true;

  return false;
  }


