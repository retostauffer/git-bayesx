
#if !defined (mixture_INCLUDED)

#define mixture_INCLUDED

#include"../export_type.h"
#include"mcmc.h"
#include"fullcond.h"
#include"distribution.h"
#include"mcmc_nonpbasis.h"
#include"mcmc_nonp.h"


namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------- class: FULLCOND_mixture --------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE FULLCOND_mixture : public FULLCOND
  {

  protected:

  unsigned nrcomp;              // Number k of mixture components
  datamatrix compweight;        // Weights w_k of mixture components
  statmatrix<unsigned> csize;   // Sizes n_k of mixture components
  datamatrix compind;           // Indicators, class probabilities for mixture components
  datamatrix compmean;          // Means mu_k of mixture components
  datamatrix compvar;           // Variances sigma_k^2 of mixture components

  datamatrix cwprior;           // Prior parameter weights of mixture components
  double cmpriorm,cmpriorv;     // Prior parameters means of mixture components
  double cvpriora,cvpriorb;     // Prior parameters variances of mixture components
  bool cvpriorbunif;            // cvpriorbunif=true  => sigma_k^2~IG(a,b)
                                //                               b~U(0,cvpriorb)
  bool cvpriorbgamma;           // cvpriorbgamma=true => sigma_k^2~IG(a,b)
                                //                               b~G(cvpriorb,(100*cvpriorb)/(cvpriora*cmpriorv))

  bool nosamples;               // samples of mixture parameters
  unsigned aclag;               // Lag for autocorrelations of mixture component parameters
                                // aclag=0 => no autocorrelations written in file
  ST::string ordertype;         // type of labelling restriction
  datamatrix temp;

  FULLCOND cpar_fc;             // FULLCOND object for component parameters
  FULLCOND cind_fc;             // FULLCOND object for component indicators
  FULLCOND_const * fcconst;
  DISTRIBUTION * likep;

  statmatrix<int> index;
  statmatrix<int> index2;
  vector<unsigned> posbeg;
  vector<unsigned> posend;
  datamatrix effvalues;

  double centerbeta(void);
  bool checkorder;
  void update_weights(void);

  public:

  // DEFAULT CONSTRUCTOR:

  FULLCOND_mixture(void) : FULLCOND()
    {
    }

  // CONSTRUCTOR1
  // mixture for random intercept

  FULLCOND_mixture(MCMCoptions * o,DISTRIBUTION * dp,
                  FULLCOND_const * fcc,
                  const datamatrix & d, const ST::string & t,
                  const ST::string & fp,const ST::string & pr,
                  const unsigned & nrc, const double & pw,
                  const double & pmm,const double & pmv,
                  const double & pva,const double & pvb,
                  const bool & s, const unsigned & acl,
                  const ST::string & ot,
                  const bool & pvbu,const bool & pvbg,
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


}   // end: namespace MCMC


#endif

