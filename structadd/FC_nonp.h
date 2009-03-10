
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FCNONPINCLUDED)

#define FCNONPINCLUDED

#include"../values.h"
#include<fstream>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include"FC.h"
#include"design.h"
#include"MASTER_obj.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_nonp -----------------------------------
//------------------------------------------------------------------------------

enum sampletype {unconstrained,increasing,decreasing};

class __EXPORT_TYPE FC_nonp  : public FC
  {

  protected:

  MASTER_OBJ * masterp;

  FC fsample;

  FC paramsample;

  bool IWLS;

  sampletype stype;

  DISTR * likep;                             // Pointer to DISTR obejct
  DESIGN * designp;                          // Pointer to design object


  datamatrix betaold;
  datamatrix betadiff;

  void centerparam(void);
  void centerparam_sample(void);


  FC pvalue_sample;                          // stores required quantities for
                                             // computing p-values
  bool pvalue;                               // compute pvalues, yes, no
                                             // default no
  datamatrix mPhelp;                         // help matrix for updating p-value
                                             // information

  void compute_pvalue(ST::string & pathresults);
  void update_pvalue(void);


  bool computemeaneffect;
  FC meaneffect_sample;
  void compute_meaneffect(void);

  public:

  datamatrix param;                          // Parameters, beta stores hatf
  datamatrix paramlin;
  datamatrix paramold;
  datamatrix paramhelp;
  datamatrix parammode;
  double paramKparam;

  datamatrix partres;                        // sum of partial residuals

  double lambda;
  double tau2;

  //---------------------------- centering -------------------------------------

  datamatrix Vcenter;
  datamatrix Vcentert;
  datamatrix Wcenter;
  datamatrix Ucenter;
  datamatrix Utc;
  datamatrix ccenter;
  datamatrix helpcenter;

  void get_linparam(void);

  void initialize_center(void);

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_nonp(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_nonp(MASTER_OBJ * mp,GENERAL_OPTIONS * o,DISTR * lp, const ST::string & t,
           const ST::string & fp,DESIGN * dp,vector<ST::string> & op,
             vector<ST::string> & vn,bool sstore);

  // COPY CONSTRUCTOR

  FC_nonp(const FC_nonp & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_nonp & operator=(const FC_nonp & m);

  // DESTRUCTOR

  ~FC_nonp()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  void update_gaussian(void);
  void update_IWLS(void);
  void update_isotonic(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void outoptions(void);

  // FUNCTION: outgraphs
  // TASK: writes batch files for STATA and R for visualizing results

  void outgraphs(ofstream & out_stata, ofstream & out_R,ST::string & path);


  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(ofstream & out_stata, ofstream & out_R,
                  const ST::string & pathresults);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void);

  virtual void transform_beta(void);

  };


} // end: namespace MCMC

#endif


