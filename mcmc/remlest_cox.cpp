
#include "first.h"

#if defined(BORLAND_OUTPUT_WINDOW)
#include "StatResults.h"
#include "statwinframe.h"
#endif

#include"remlest.h"
#include"remlest_cox.h"

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


  // implement AFTs here



  out("ESTIMATION RESULTS:\n",true);
  out("\n");

  for(i=1;i<fullcond.size();i++)
    {
//    beta(0,0) += fullcond[i]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[i],zcut[i-1],i-1,false,xcut[i],X.cols()+zcut[i-1],0,false,i);
    }
//  beta(0,0) += fullcond[0]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[0],0,0,false,xcut[0],0,0,false,0);

  return false;
  }

//------------------------------------------------------------------------------
//----------------- AFT models with smooth error distribution ------------------
//---------------------------- fixed effects only ------------------------------
//------------------------------------------------------------------------------


  bool remlest::estimate_aft_glm(datamatrix resp, const datamatrix & offset,
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


  // implement AFTs here



  out("ESTIMATION RESULTS:\n",true);
  out("\n");

//  beta(0,0) += fullcond[0]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[0],0,0,false,xcut[0],0,0,false,0);

  return false;
  }

