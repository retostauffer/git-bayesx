
#if !defined (FCNONPVARIANCEVECINCLUDED)

#define FCNONPVARIANCEVECINCLUDED

#include"../export_type.h"
#include"../values.h"
#include<fstream>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include"FC_nonp.h"
#include"design.h"
#include"FC_nonp_variance.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_nonp_variance_vec ----------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_nonp_variance_vec  : public FC_nonp_variance
  {

  protected:


  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_nonp_variance_vec(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_nonp_variance_vec(MASTER_OBJ * mp,GENERAL_OPTIONS * o,DISTR * lp,
           const ST::string & t,
           const ST::string & fp,DESIGN * dp,FC_nonp * FCn,
           vector<ST::string> & op,vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  FC_nonp_variance_vec(const FC_nonp_variance_vec & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_nonp_variance_vec & operator=(const FC_nonp_variance_vec & m);

  // DESTRUCTOR

  ~FC_nonp_variance_vec()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void outoptions(void);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(ofstream & out_stata,ofstream & out_R,
                  const ST::string & pathresults);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  // virtual void transform_beta(void);

  };


} // end: namespace MCMC

#endif


