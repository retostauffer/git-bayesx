

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

  int nrcomp;

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

  void init_spatialtotal(vector<ST::string> & ev, const ST::string & pnt,
                         const ST::string & prt);


  bool changingweight;

  double centerbeta(void);

  void update_linpred(const bool & add);

  void update_linpred_diff(datamatrix & b1,datamatrix & b2);


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

  // CONSTRUCTOR2
  // random slope

  FULLCOND_mixture(MCMCoptions * o,DISTRIBUTION * dp,
                  FULLCOND_const * fcc,
                  const datamatrix & intvar,const datamatrix & effmod,
                  const ST::string & t,
                  const ST::string & fp,const ST::string & pr,
                  const ST::string & prf,
                  const double & la, const bool & inclfix,
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

  double compute_quadform(void);

  void update_sigma2(const double & s)
    {
    sigma2 = s;
    }

  unsigned get_rankK(void);

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

  datamatrix mu;

  FULLCOND_nonp_basis * fbasisp;

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
    }

  // CONSTRUCTOR2
  // random slope

  FULLCOND_mixture_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                           FULLCOND_const * fcc,
                           const datamatrix & intvar,
                           const datamatrix & effmod,
                           const ST::string & t,
                           const ST::string & fp,const ST::string & pr,
                           const ST::string & prf,
                           const double & la,
                           const bool & inclfixed,
                           const unsigned & c = 0)
                           : FULLCOND_mixture(o,dp,fcc,intvar,effmod,t,
                                             fp,pr,prf,la,inclfixed,c)
    {
    mu = datamatrix(index.rows(),1);
    }


  // COPY CONSTRUCTOR

  FULLCOND_mixture_gaussian(const FULLCOND_mixture_gaussian & fc)
  : FULLCOND_mixture(FULLCOND_mixture(fc))
    {
    mu = fc.mu;
    muy = fc.muy;
    fbasisp = fc.fbasisp;
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

  void init_spatialtotal(FULLCOND_nonp_basis * sp,const ST::string & pnt,
                         const ST::string & prt);

  // FUNCTION: reset
  // TASK: resets all parameters for a new simulation

  void reset(void)
    {
    FULLCOND_mixture::reset();
    }

  };     // end: class FULLCOND_mixture_gaussian


}   // end: namespace MCMC


#endif

