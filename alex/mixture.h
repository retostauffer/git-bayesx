#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (mixture_INCLUDED)

#define mixture_INCLUDED


#include<mcmc.h>
#include<fullcond.h>
#include<distribution.h>
#include<mcmc_nonpbasis.h>
#include<mcmc_nonp.h>


namespace MCMC
{


class __EXPORT_TYPE FULLCOND_mixture : public FULLCOND
  {

  protected:

  int nrcomp;  // Number of mixture components
  datamatrix compweight;  // Weights of mixture components
  datamatrix cwprior;  // Prior parameter weights of mixture components
  statmatrix<unsigned> csize;  // Sizes of mixture components
  statmatrix<unsigned> compind;  // Indicator for mixture components

  datamatrix compmean;
  datamatrix compvar;
  double cmpriorm,cmpriorv;
  double cvpriorsh,cvpriorsc;

  bool checkorder;
  datamatrix temp; //


  FULLCOND_const * fcconst;
  DISTRIBUTION * likep;

  statmatrix<int> index;
  statmatrix<int> index2;
  vector<unsigned>  posbeg;
  vector<unsigned>  posend;
  datamatrix effvalues;

  double centerbeta(void);
  const update_weights(void);


  
  public:

  // DEFAULT CONSTRUCTOR:

  FULLCOND_mixture(void) : FULLCOND()
    {
    }

  // CONSTRUCTOR1
  // random intercept

  FULLCOND_mixture(MCMCoptions * o,DISTRIBUTION * dp,
                  FULLCOND_const * fcc,
                  const datamatrix & d, const ST::string & t,
                  const ST::string & fp,const ST::string & pr, const int & nrc,
                  const unsigned & c=0);


  void init_name(const ST::string & na);

  void init_names(const vector<ST::string> & na);

  // COPY CONSTRUCTOR

  FULLCOND_mixture(const FULLCOND_mixture & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_mixture & operator=(
                        const FULLCOND_mixture & fc);

  // DESTRUCTOR

  ~FULLCOND_mixture() {}


  // FUNCTION update
  // TASK: updates parameters (i.e. matrix beta)

  void update(void);

  bool posteriormode(void);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file (after estimation)

  void outresults(void);

  // FUNCTION: outoptions

  void outoptions(void);

  // FUNCTION: reset
  // TASK: resets all parameters for a new simulation

  void reset(void)
    {
    FULLCOND::reset();
    }

  };     // end: class FULLCOND_mixture


//------------------------------------------------------------------------------
//--------------------- class: FULLCOND_mixture_gaussian ------------------------
//------------------------------------------------------------------------------

/*
class __EXPORT_TYPE FULLCOND_mixture_gaussian : public FULLCOND_mixture
  {

  protected:

//  datamatrix compmean; // Means of normal mixture components
//  datamatrix compvar;  // Variances of normal mixture components

//  datamatrix mu;

  public:

  // DEFAULT CONSTRUCTOR:

  FULLCOND_mixture_gaussian(void) : FULLCOND_mixture()
    {
    }

  // CONSTRUCTOR1
  // random intercept

  FULLCOND_mixture_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                           FULLCOND_const * fcc,
                           const datamatrix & d, const ST::string & t,
                           const ST::string & fp,const ST::string & pr,
                           const int & nrc, const double & la,
                           const unsigned & c = 0)
                           : FULLCOND_mixture(o,dp,fcc,d,t,fp,pr,nrc,la,c)
    {
//    mu = datamatrix(index.rows(),1);

//    compmean = datamatrix(nrcomp,1,0);
//    compvar = datamatrix(nrcomp,1,1);
    }

  // COPY CONSTRUCTOR

  FULLCOND_mixture_gaussian(const FULLCOND_mixture_gaussian & fc)
  : FULLCOND_mixture(FULLCOND_mixture(fc))
    {
//    mu = fc.mu;
//    muy = fc.muy;

//    compmean=fc.compmean;
//    compvar=fc.compvar;
    }

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_mixture_gaussian & operator=(
                        const FULLCOND_mixture_gaussian & fc)
    {
    if (this==&fc)
      return *this;
    FULLCOND_mixture::operator=(FULLCOND_mixture(fc));
//    mu = fc.mu;
//    muy = fc.muy;

//    compmean=fc.compmean;
//    compvar=fc.compvar;
    return *this;
    }

  // DESTRUCTOR

  ~FULLCOND_mixture_gaussian() {}

  // FUNCTION: update
  // TASK: updates parameters (i.e. matrix beta)

  void update(void);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file (after estimation)

  void outresults(void)
    {
    FULLCOND_mixture::outresults();
    }

  // FUNCTION: outoptions

  void outoptions(void)
    {
    FULLCOND_mixture::outoptions();
    }

  // FUNCTION: reset
  // TASK: resets all parameters for a new simulation

  void reset(void)
    {
    FULLCOND_mixture::reset();
    }

  };     // end: class FULLCOND_mixture_gaussian
*/

}   // end: namespace MCMC


#endif

