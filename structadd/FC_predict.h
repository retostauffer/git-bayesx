
#if !defined (FCpredictINCLUDED)

#define FCpredictINCLUDED

#include"../export_type.h"
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

  void compute_MSE(const ST::string & pathresults);

  public:

  msetype MSE;
  double MSEparam;

  ST::string getloss(void);

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

  void compute_autocorr_all(const ST::string & path, unsigned lag,
                                    ofstream & outg) const;

  void reset(void);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  };


} // end: namespace MCMC

#endif


