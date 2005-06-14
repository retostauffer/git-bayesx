#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (TVARIANCE2DIM_INCLUDED)
#define TVARIANCE2DIM_INCLUDED

#include<fullcond_pspline_surf_gaussian.h>
#include<fullcond_nonp_gaussian.h>

namespace MCMC
{



class __EXPORT_TYPE FULLCOND_tvariance2dim : public FULLCOND
  {


  protected:

  FULLCOND_pspline_surf_gaussian * Kp;  // pointer to psplines full conditional
  FULLCOND_nonp_gaussian * Kp_spat;  // pointer to psplines full conditional
  bool spatial;

  statmatrix<int> indexmat; 

  double * Kmatdiag;                    //
  double * Kmatupper;

  ST::string pathresults;               // path for results

  unsigned nu;                          // hyperparameter nu

  unsigned m;                           // number of parameters per row
                                        // (or column) in the psplines fc

  bool rowwise;
  datamatrix u;

  envmatdouble K11;
  double detalt;
  double detneu;

  unsigned nrrows;

  vector<double> deltapropvec;
  vector<double> rowvec;
  vector<double> colvec;
  vector<double> betakvec;


  public:


  // DEFAULT CONSTRUCTOR

  FULLCOND_tvariance2dim(void) : FULLCOND()
    {
    }


  // CONSTRUCTOR

  FULLCOND_tvariance2dim(MCMCoptions * o,FULLCOND_pspline_surf_gaussian * p,
                     unsigned & v,const ST::string & ti, const ST::string & fp,
                     const ST::string & pres,const unsigned & bs,const bool & rw = false);

  FULLCOND_tvariance2dim(MCMCoptions * o,FULLCOND_nonp_gaussian * p,
                     unsigned & v,const ST::string & ti, const ST::string & fp,
                     const ST::string & pres,const unsigned & bs,const bool & rw = false);


  // COPY CONSTRUCTOR

  FULLCOND_tvariance2dim(const FULLCOND_tvariance2dim & t);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_tvariance2dim & operator=(const FULLCOND_tvariance2dim & t);

  void update(void);

  void update_2dim(void); 

  void update_spat(void);

  bool posteriormode(void)
    {
    return true;
    }

  void outresults(void);

  void outoptions(void);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    FULLCOND::reset();
    setbeta(nrpar,1,1);
    }

  // FUNCTION: predict (virtual)
  // TASK: predicts the mean for a new observation Xnew

  void predict(const datamatrix & newX, datamatrix & linpred)
    {

    }

  // DESTRUCTOR

  ~FULLCOND_tvariance2dim() {}

  }; // end: class FULLCOND_tvariance



} // end: namespace MCMC

#endif




