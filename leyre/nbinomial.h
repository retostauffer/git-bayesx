// DATE: 22.06.99


#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif


#if !defined (DISTRIBUTION_nbinomial_INCLUDED)

#define DISTRIBUTION_nbinomial_INCLUDED

#include <distribution.h>
#include <random.h>




namespace MCMC
{


// full conditional classes using conditional prior proposals work properly
// provided that the following virtual functions are implemented:

//  double loglikelihood(const double * response,const double * linpred,
//                       const double * weight,unsigned & i) const


// full conditional class for fixed effects (non Gaussian case) work properly
// provided that the following virtual functions are implemented:

//  double compute_weight(double * linpred,double * weight,
//                        const unsigned & col=0) const

//  double compute_gmu(double * linpred) const


// predicted linearpred, deviance, mu etc. will be computed correctly if the
// following functions are implemented:

//  void compute_mu(const double * linpred,double * mu) const

//  void compute_devresidual(const double * response,const double * weight,
//                           const double * mu,
//                           double * residual,const datamatrix & scale, const int & i) const

enum vertopt {nb,poga,poig};
enum propscale {unif,gam};


class __EXPORT_TYPE DISTRIBUTION_nbinomial : public DISTRIBUTION
  {

   protected:

   // include here additional private variables needed for the new
   // distribution

   bool oversize;
   datamatrix accept;
   datamatrix nu;
   FULLCOND nusave;
   FULLCOND nusavekfz;   
   datamatrix hierint;
   FULLCOND hierintsave;
   datamatrix pvar;
   double a_pri;
   datamatrix b_pri;
   FULLCOND b_pri_save;
   double prop_var;
   vertopt ver;
   propscale pscale;
   datamatrix sum_nu;
   datamatrix sum2_nu;

   bool hierarchical;

   public:

   //  the following public functions must be implemented


   // DEFAULT CONSTRUCTOR

    DISTRIBUTION_nbinomial(void);


   // CONSTRUCTOR WITHOUT OFFSET

   DISTRIBUTION_nbinomial(const double & a,const double & b,
                          const double & pv,const vertopt & vo,
                          const propscale & psc,bool hie,
                          MCMCoptions * o,
                          const datamatrix & r,
                          const ST::string & p,const ST::string & ps,
                          const datamatrix & w=datamatrix());


   // CONSTRUCTOR WITH OFFSET
   DISTRIBUTION_nbinomial(const double & a,const double & b,
                          const double & pv,const vertopt & vo,
                          const propscale & psc,bool hie,
                          const datamatrix & offset, MCMCoptions * o,
                          const datamatrix & r,
                          const ST::string & p,const ST::string & ps,
                          const datamatrix & w=datamatrix());


void DISTRIBUTION_nbinomial::create(MCMCoptions * o, const double & a,
                                    const double & b, const double & pv,
                                    const vertopt & vo, const propscale & psc,
                                    bool hie, const ST::string & ps);




   // COPY CONSTRUCTOR

   DISTRIBUTION_nbinomial(const DISTRIBUTION_nbinomial & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_nbinomial & operator=(const DISTRIBUTION_nbinomial & nd);

   // DESTRUCTOR

   ~DISTRIBUTION_nbinomial() {}

   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * response,double * linpred,
                       double * weight,const int & i) const;

  // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred' and stores
  //       the result in 'mu'

  void compute_mu(const double * linpred,double * mu) const;

  // FUNCTION: compute_deviance
  // TASK: computes the individual deviance

  void compute_deviance(const double * response,const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           const datamatrix & scale, const int & i) const;


   // FUNCTION: compute_weight
   // TASK: computes the weights for iteratively weighted least squares
   //       w_i = weight_i/scale*[b''(theta_i) g'^2(mu_i)]^{-1}

  double compute_weight(double * linpred,double * weight,
                        const int & i,const unsigned & col=0) const;

  // FUNCTION: compute_gmu
  // TASK: compute g'(eta_i) = 1/h'(eta_i)

  double compute_gmu(double * linpred,const unsigned & col=0) const;

  void compute_mu_notransform(const double * linpred,double * mu) const;

  double compute_IWLS(double * response,double * linpred,double * weight,
                      const int & i,double * weightiwls,double * tildey,
                      bool weightyes, const unsigned & col=0);

  void compute_IWLS_weight_tildey(double * response,double * linpred,
                              double * weight,const int & i,double * weightiwls,
                              double * tildey,const unsigned & col=0);


  // FUNCTION: outoptions
  // TASK: writing options of the distribution

  void outoptions(void);

  void outresults(void);

  // FUNCTION: update
  // TASK: updates the scale parameter  and the intercept

  void update(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode for the scale parameter and the
  //       intercept

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  double update_nu(void);

  double update_hierint(void) const;

  double update_scale(void) const;

  double update_b_pri(void) const;

  double proposal_scale(void) const;

  double pwork_tunin(unsigned i) const;

    //Logarithm of the density f(x)~G(a, b)

  double log_gamma(double & x, double & a, double & b) const;

  double log_gamma_likelihood(double &s, double &s_neu) const;

  double log_gamma_likelihood_hier(double &s, double &s_neu) const;

  double log_gamma_quot(double & x, double & x_prop, const double & a,
                        const double & b) const;

  double log_gamma_prop(double & x, double & x_prop) const;

    // loggamma-function

  double lgamma(const double & xx) const;


  double log_nbin(const double & s_neu, const double & s) const;

  double log_gamma_function_quot(const double & g, const double & g_prop,
                                 const unsigned & y) const;

  const double & get_sum_nu(void) const
    {
    return sum_nu(0,0);
    }

  const double & get_sum2_nu(void) const
    {
    return sum2_nu(0,0);
    }

  const double & get_hierint(void) const
    {
    return hierint(0,0);
    }

  const int get_distopt(void) const
    {
    if(ver==poga) return 1;
    else return 2;
    }

  const double & get_pvar(void) const
    {
    return pvar(nrobs+1, 0); //oder nrobs+2???????
    }

  void initialize_hierint(double & inter)
    {
    hierint(0,0) = inter;
    }

  void exchange_hierint(double & inter)
    {
    hierint(0,0) += inter;
    }

  void exchange_accept(void)
    {
        double *acceptwork = accept.getV();
        acceptwork += nrobs +1;
        *acceptwork += 1;

    }

  void add_nu(double m) const;    


  };   // end: DISTRIBUTION_nbinomial


}  // end: namespace MCMC

#endif
