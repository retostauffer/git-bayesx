

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
  statmatrix<int> csize;  // Sizes of mixture components
  statmatrix<int> compind;  // Indicator for mixture components

  datamatrix muy;


  FULLCOND_const * fcconst;

  DISTRIBUTION * likep;

  statmatrix<int> index;
  statmatrix<int> index2;

  vector<unsigned>     posbeg;
  vector<unsigned>     posend;

  datamatrix XX;
  datamatrix effvalues;
  double sigma2;                            // prior variance parameter
  double lambda;
  double lambdaold1;
  double lambdaold2;
  double df_lambdaold1;
  double df_lambdaold2;
  bool lambdaconst;

  bool randomslope;
  bool includefixed;

  datamatrix data2;

  bool spatialtotal;
  statmatrix<int> indextotal;
  ST::string pathsample_total;


  FULLCOND ftotal;


  bool changingweight;

  double centerbeta(void);

  void update_linpred(const bool & add);



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
                  const double & la, const unsigned & c=0);


  void init_name(const ST::string & na);

  void init_names(const vector<ST::string> & na);

  // COPY CONSTRUCTOR

  FULLCOND_mixture(const FULLCOND_mixture & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_mixture & operator=(
                        const FULLCOND_mixture & fc);

  // DESTRUCTOR

  ~FULLCOND_mixture() {}

  void compute_XWX(const datamatrix & weightmat,const unsigned & col);

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
    sigma2 = 10;
    }


  void update_sigma2(const double & s)
    {
    sigma2 = s;
    }

  double getlambda(void)
    {
    return lambda;
    }

  double get_sigma2(void)
    {
    return sigma2;
    }

  void set_lambdaconst(double la);

  };     // end: class FULLCOND_mixture


//------------------------------------------------------------------------------
//--------------------- class: FULLCOND_mixture_gaussian ------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_mixture_gaussian : public FULLCOND_mixture
  {

  protected:

//  int nrcomp;  // Number of mixture components
//  datamatrix compweight;  // Weights of mixture components
  statmatrix<int> compind;  // Indicator for mixture components
  datamatrix compmean; // Means of normal mixture components
  datamatrix compvar;  // Variances of normal mixture components

  datamatrix mu;

  FULLCOND_nonp_basis * fbasisp;

/*
//  Update functions for mixture component parameters

  void update_compmean(datamatrix betaakt, statmatrix<int> compindakt,datamatrix compvarakt, datamatrix priormean, datamatrix priorvar)
  {
     unsigned i;
     for(i=0;i<compmean.rows();i++)
     {
     compmean(i,0) = compmean(i,0)+rand_normal();
     }
  }

  void update_compvar(datamatrix betaakt, statmatrix<int> compindakt, datamatrix compmeanakt, datamatrix priorshape,datamatrix priorscale)
  {
     unsigned i;
     for(i=0;i<compvar.rows();i++)
     {
     compvar(i,0) = compvar(i,0)+rand_normal();
     }
  }
*/


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
    mu = datamatrix(index.rows(),1);

    compind = statmatrix<int>(d.rows(),1,1);
    compmean = datamatrix(nrcomp,1,0);
    compvar = datamatrix(nrcomp,1,1);
    }

  // COPY CONSTRUCTOR

  FULLCOND_mixture_gaussian(const FULLCOND_mixture_gaussian & fc)
  : FULLCOND_mixture(FULLCOND_mixture(fc))
    {
    mu = fc.mu;
    muy = fc.muy;
    fbasisp = fc.fbasisp;

    compind=fc.compind;
    compmean=fc.compmean;
    compvar=fc.compvar;
    }

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_mixture_gaussian & operator=(
                        const FULLCOND_mixture_gaussian & fc)
    {
    if (this==&fc)
      return *this;
    FULLCOND_mixture::operator=(FULLCOND_mixture(fc));
    mu = fc.mu;
    muy = fc.muy;
    fbasisp = fc.fbasisp;

    compind=fc.compind;
    compmean=fc.compmean;
    compvar=fc.compvar;
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


}   // end: namespace MCMC


#endif

