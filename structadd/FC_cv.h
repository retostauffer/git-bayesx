
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FCcvINCLUDED)

#define FCcvINCLUDED

#include"statmat.h"
#include"sparsemat.h"

#include"Random.h"
#include"../values.h"
#include<fstream>
#include<vector>
#include<bitset>
#include"GENERAL_OPTIONS.h"
#include"distr.h"
#include"FC_hrandom.h"
#include"clstring.h"
#include<cmath>

namespace MCMC
{

using std::vector;
using std::bitset;

//------------------------------------------------------------------------------
//------------------------------ CLASS: FC_cv ----------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE FC_cv   : public FC
  {

  protected:


  datamatrix sampled_etas;
  datamatrix sampled_responses;

  DISTR * likep;

  vector<FC_hrandom> * hrandoms;
  unsigned size;

  datamatrix effect;
  datamatrix linpred;


  public:

  // DEFAULT CONSTRUCTOR

  FC_cv(void);

  // CONSTRUCTOR

  FC_cv(GENERAL_OPTIONS * o,DISTR * lp,const ST::string & t,
     const ST::string & fp,vector<FC_hrandom> * FChs);

  // COPY CONSTRUCTOR

  FC_cv(const FC_cv & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_cv & operator=(const FC_cv & m);

  // DESTRUCTOR

  ~FC_cv()
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


