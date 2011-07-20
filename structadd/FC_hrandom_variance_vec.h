
#if !defined (FChrandomVARIANCEVECINCLUDED)

#define FChrandomVARIANCEVECINCLUDED

#include"../export_type.h"
#include"../values.h"
#include<fstream>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include"FC_nonp_variance_vec.h"
#include"design.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- CLASS: FC_hrandom_variance_vec -----------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_hrandom_variance_vec  : public FC_nonp_variance_vec
  {

  protected:

  DISTR * likepRE;

  bool mult;

  double hyperLambda;

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_hrandom_variance_vec(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_hrandom_variance_vec(MASTER_OBJ * mp,GENERAL_OPTIONS * o,DISTR * lp,
                          DISTR * lpRE, const ST::string & t,
                          const ST::string & fp,DESIGN * dp,
                          FC_nonp * FCn,vector<ST::string> & op,
                          vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  FC_hrandom_variance_vec(const FC_hrandom_variance_vec & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_hrandom_variance_vec & operator=(const FC_hrandom_variance_vec & m);

  // DESTRUCTOR

  ~FC_hrandom_variance_vec()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

  // void transform_beta(void);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  };


} // end: namespace MCMC

#endif


