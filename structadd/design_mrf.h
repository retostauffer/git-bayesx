
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (DESIGNINCLUDED)

#define DESIGNINCLUDED

#include"statmat.h"
#include"sparsemat.h"

#include"random.h"
#include"envmatrix_penalty.h"
#include"../values.h"
#include<fstream.h>
#include<vector>
#include<bitset>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include"distr.h"
#include<cmath>


namespace MCMC
{


enum ttype {
                RE,
                RW1,
                RW2,
                RW3,
                RW1RW2,
                RW1RW2RW3,
                seasonal,
                mrf,
                mrfI,
                mrfkronecker,
                mrflinear,
                mrflinearband,
                mrfquadratic8,
                mrfquadratic12,
                mrfkr1,
                mrfkr2,
                npspline,
                smoothspline,
                kriging,
                hrandom
                };


//------------------------------------------------------------------------------
//--------------------------- CLASS: DESIGN_mrf --------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN_mrf : public DESIGN
  {

  MAP::map ma;
  datamatrix data2;                          // vor varying coefficients

  protected:



  public:


//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN_mrf(void);

  // CONSTRUCTOR 1
  // Spatial covariates

  DESIGN_mrf(const datamatrix & dm, const datamatrix & iv,
             DISTR * dp, const MAP::map & m);

  // COPY CONSTRUCTOR

  DESIGN_mrf(const DESIGN_mrf & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN_mrf & operator=(const DESIGN_mrf & m);

  // virtual functions

  void init_data(const datamatrix & dm,const datamatrix & iv);

  void compute_penalty(void);

  void compute_XtransposedWX_XtransposedWres(datamatrix & partres, double l);

  void compute_XtransposedWres(datamatrix & partres, double l);

  void compute_precision(double l);

  // DESTRUCTOR

  ~DESIGN_mrf() {}

  };


} // end: namespace MCMC

#endif


