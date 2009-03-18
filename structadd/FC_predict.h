
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FCpredictINCLUDED)

#define FCpredictINCLUDED

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
//--------------------------- CLASS: FC_predict --------------------------------
//------------------------------------------------------------------------------

enum msetype{noMSE,yesMSE};


class __EXPORT_TYPE FC_predict   : public FC
  {

  protected:

  FC FC_deviance;

  DISTR * likep;
  datamatrix designmatrix;
  vector<ST::string> varnames;


  double deviance;
  double deviancesat;


  void get_predictor(void);

  void compute_MSE(ST::string & pathresults);

  public:

  msetype MSE;

  // DEFAULT CONSTRUCTOR

  FC_predict(void);

  // CONSTRUCTOR

  FC_predict(GENERAL_OPTIONS * o,DISTR * lp,const ST::string & t,
     const ST::string & fp,const ST::string & fpd, datamatrix & dm,
     vector<ST::string> & dn);

  // COPY CONSTRUCTOR

  FC_predict(const FC_predict & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_predict & operator=(const FC_predict & m);

  // DESTRUCTOR

  ~FC_predict()
    {
    }


  void update(void);

  bool posteriormode(void);

  void outoptions(void);

  void outresults_deviance(void);
  void outresults_DIC(void);
  void outresults(ofstream & out_stata, ofstream & out_R,
                  const ST::string & pathresults);

  void reset(void);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  };


} // end: namespace MCMC

#endif


