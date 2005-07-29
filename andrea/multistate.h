//---------------------------------------------------------------------------
#ifndef multistateH
#define multistateH

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif

#include <distribution.h>
#include <mcmc_pspline.h>
#include <baseline.h>

namespace MCMC
{

//------------------------------------------------------------------------------
//------------------------ CLASS: DISTRIBUTION_multistatemodel ------------------------
//------------------------------------------------------------------------------



class __EXPORT_TYPE DISTRIBUTION_multistatemodel : public DISTRIBUTION
  {

  protected:

  unsigned nrcat;
  unsigned nrtransition;

  datamatrix int_ti;
  datamatrix mean_int_ti;

  datamatrix ti;
  datamatrix state_i;
  datamatrix transition;
  datamatrix transition_help;
//  datamatrix relrisk;
//  bool offsetexisting;


  public:


  double*  get_integral_ti (void )
  {
  return int_ti.getV();
  }



   // DEFAULT CONSTRUCTOR


   DISTRIBUTION_multistatemodel(void) : DISTRIBUTION()
     {
     }

   // CONSTRUCTOR

   DISTRIBUTION_multistatemodel(MCMCoptions * o, const datamatrix & r,
   const datamatrix & t,  const datamatrix & dbeg,
   const datamatrix & state , const datamatrix & w=datamatrix());


   // COPY CONSTRUCTOR

   DISTRIBUTION_multistatemodel(const DISTRIBUTION_multistatemodel & nd)
   : DISTRIBUTION(DISTRIBUTION(nd))
     {

     nrcat = nd.nrcat;
     nrtransition = nd.nrtransition;
     ti=nd.ti;
     state_i = nd.state_i;
     transition_help = nd.transition_help;
     transition = nd.transition;

     int_ti = nd.int_ti;
     mean_int_ti = nd.mean_int_ti;
//     offsetexisting = nd.offsetexisting;
//     relrisk = nd.relrisk;

     }

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_multistatemodel & operator=(const DISTRIBUTION_multistatemodel & nd)
     {

     if (this==&nd)
	   return *this;

     DISTRIBUTION::operator=(DISTRIBUTION(nd));

     nrcat = nd.nrcat;
     nrtransition = nd.nrtransition;
     ti=nd.ti;
     state_i = nd.state_i;
     transition_help = nd.transition_help;
     transition = nd.transition;

     int_ti= nd.int_ti;
     mean_int_ti = nd.mean_int_ti;
//     offsetexisting = nd.offsetexisting;
//     relrisk = nd.relrisk;
     return *this;
     }

   // DESTRUCTOR

   ~DISTRIBUTION_multistatemodel() {}

   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation



  double loglikelihood(double * response,double * linpred,
                       double * weight,const int & i) const;




 // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred'

  void compute_mu(const double * linpred,double * mu) const;



  void compute_mu_notransform(const double * linpred,double * mu) const
    {
    for(unsigned i=0;i<linearpred.cols();i++,linpred++,mu++)
      *mu = exp(*linpred);
    }


  // FUNCTION: compute_devresidual
  // TASK: computes the deviance residual

  void compute_deviance(const double * response,const double * weight,
                           const double * mu,double * deviance,double * deviancesat,
                           const datamatrix & scale,const int & i) const;


  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

//  void update_predict(void);

  bool posteriormode(void);

  void outresults(void);

  double compute_weight(double * linpred, double * weight, const int & i,
                                                   const unsigned & col) const;

  void tilde_y(datamatrix & tildey,datamatrix & m, const unsigned & col,
                const bool & current, const datamatrix & w);

  void compute_iwls(void);

  void compute_IWLS_weight_tildey(double * response,double * linpred,
                              double * weight,const int & i,
                              double * weightiwls,double * tildey,
                              const unsigned & col=0);

  double compute_IWLS(double * response,double * linpred,double * weight,const int & i,
                              double * weightiwls, double * tildey,
                              bool weightyes,const unsigned & col=0);

  void assign_int_ti(const datamatrix & m)
    {
    int_ti = m;
    }

  void assign_mean_int_ti(void)
    {
    mean_int_ti = int_ti;
    }

  datamatrix & get_int_ti(void)
    {
    return int_ti;
    }

  

  unsigned get_nrcat(void)
    {
    return nrcat;
    }

  unsigned get_nrtransition(void)
    {
    return nrtransition;
    }


  };  // end: class DISTRIBUTION_multistatemodel

} // END: namespace MCMC


//---------------------------------------------------------------------------
#endif
