#include "first.h"

#include "fullcond_mult.h"

namespace MCMC
{



FULLCOND_mult::FULLCOND_mult(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_random * rp,
                         FULLCOND_nonp_basis * ba,
                         const ST::string & ti,
                         const ST::string & fp, const ST::string & pres,
                         const unsigned & c)
  {


  }


  // COPY CONSTRUCTOR

FULLCOND_mult::FULLCOND_mult(const FULLCOND_mult & fc)
  {
  basis1p = fc.basis1p;
  basis2p = fc.basis2p;
  reffectp = fc.reffectp;
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const FULLCOND_mult & FULLCOND_mult::operator=(const FULLCOND_mult & fc)
  {
  if (this == &fc)
    return *this;
  basis1p = fc.basis1p;
  basis2p = fc.basis2p;
  reffectp = fc.reffectp;
  return *this;
  }

void FULLCOND_mult::update(void)
  {


  }

//  void update_linpred(const bool & add)
//    {
//    }


bool FULLCOND_mult::posteriormode(void)
  {


  }

bool FULLCOND_mult::posteriormode_converged(const unsigned & itnr)
  {


  }


void FULLCOND_mult::outresults(void)
  {

  }

void FULLCOND_mult::get_effectmatrix(datamatrix & e,vector<ST::string> & enames,
                        unsigned be, unsigned en,effecttype t)
  {

  }


unsigned FULLCOND_mult::get_nreffects(effecttype t)
  {


  }


void FULLCOND_mult::outoptions(void)
  {


  }

ST::string FULLCOND_mult::getinfo(void)
  {

  }


void FULLCOND_mult::init_name(const ST::string & na)
  {

  }

void FULLCOND_mult::init_names(const vector<ST::string> & na)
  {


  }

void FULLCOND_mult::init_priorassumptions(const ST::string & na)
  {


  }

  // FUNCTION: reset
  // TASK: resets all parameters

void FULLCOND_mult::reset(void)
  {

  }


  } // end: namespace MCMC
