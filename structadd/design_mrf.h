
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (DESIGNmrfINCLUDED)

#define DESIGNmrfINCLUDED

#include"statmat.h"
//#include"sparsemat.h"

#include"design.h"
#include"random.h"
#include"envmatrix_penalty.h"
#include"../values.h"
#include<fstream.h>
#include<vector>
#include<bitset>
#include"clstring.h"
#include<cmath>


namespace MCMC
{


//------------------------------------------------------------------------------
//--------------------------- CLASS: DESIGN_mrf --------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN_mrf : public DESIGN
  {

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  MAP::map ma;

  protected:

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN_mrf(void);

  // CONSTRUCTOR 1
  // Spatial covariates

  DESIGN_mrf(const datamatrix & dm, const datamatrix & iv,
             DISTR * dp, const MAP::map & m,vector<ST::string> & op,
             vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  DESIGN_mrf(const DESIGN_mrf & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN_mrf & operator=(const DESIGN_mrf & m);

  // VIRTUAL FUNCTIONS

  void init_data(const datamatrix & dm,const datamatrix & iv);

  void compute_penalty(void);

  void compute_XtransposedWX_XtransposedWres(datamatrix & partres, double l);

  void compute_XtransposedWres(datamatrix & partres, double l);

  void compute_XtransposedWX(void);

  void compute_precision(double l);

  // DESTRUCTOR

  ~DESIGN_mrf() {}

  };


} // end: namespace MCMC

#endif


