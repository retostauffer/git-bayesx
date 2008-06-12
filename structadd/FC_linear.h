
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FClinearINCLUDED)

#define FClinearINCLUDED

#include"../values.h"
#include<fstream.h>
#include"GENERAL_OPTIONS.h"
#include"distr.h"
#include"clstring.h"
#include"FC.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_linear ---------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_linear  : public FC
  {

  protected:

  bool IWLS;

  DISTR * likep;                             // Pointer to DISTR obejct
  datamatrix design;

  datamatrix Xt;                             // transposed designmatrix
  datamatrix XWX;

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_linear(void);

  // CONSTRUCTOR

  FC_linear(GENERAL_OPTIONS * o,DISTR * lp, datamatrix & d,const ST::string & t,
           const ST::string & fp,vector<ST::string> & op,
             vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  FC_linear(const FC_linear & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_linear & operator=(const FC_linear & m);

  // DESTRUCTOR

  ~FC_linear()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  void update_gaussian(void);
  void update_IWLS(void);

  // FUNCTION: posteriormode

  bool posteriormode(void);

  // FUNCTION: outoptions

  void outoptions(void);

  // FUNCTION: outresults

  void outresults(const ST::string & pathresults);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  // FUNCTION: reset

  void reset(void);

  };


} // end: namespace MCMC

#endif


