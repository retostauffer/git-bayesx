
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FCpredictivecheckINCLUDED)

#define FCpredictivecheckINCLUDED

#include"statmat.h"
#include"sparsemat.h"

#include"Random.h"
#include"../values.h"
#include<fstream>
#include<vector>
#include<bitset>
#include"GENERAL_OPTIONS.h"
#include"distr.h"
#include"clstring.h"
#include<cmath>

namespace MCMC
{

using std::vector;
using std::bitset;

//------------------------------------------------------------------------------
//------------------------- CLASS: FC_predictive_check -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE FC_predictive_check   : public FC
  {

  protected:


  datamatrix sampled_responses;

  DISTR * likep;

  datamatrix designmatrix;
  vector<ST::string> varnames;

  void get_predictor(void);


  public:

  // DEFAULT CONSTRUCTOR

  FC_predictive_check(void);

  // CONSTRUCTOR

  FC_predictive_check(GENERAL_OPTIONS * o,DISTR * lp,const ST::string & t,
     const ST::string & fp,datamatrix & dm, vector<ST::string> & dn);

  // COPY CONSTRUCTOR

  FC_predictive_check(const FC_predictive_check & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_predictive_check & operator=(const FC_predictive_check & m);

  // DESTRUCTOR

  ~FC_predictive_check()
    {
    }


  void update(void);

  bool posteriormode(void);

  void outoptions(void);

  void outresults(ofstream & out_stata, ofstream & out_R,
                  const ST::string & pathresults);

  void reset(void);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  };


} // end: namespace MCMC

#endif


