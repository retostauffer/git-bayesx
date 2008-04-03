
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
//--------------------------- CLASS: DESIGN_hrandom ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN_hrandom : public DESIGN
  {


  protected:

  DISTR * likep_RE;

  public:


//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN_hrandom(void);

  // CONSTRUCTOR 1
  // Spatial covariates

  DESIGN_hrandom(const datamatrix & dm, const datamatrix & iv,
             DISTR * dp,DISTR * dp_RE);

  // COPY CONSTRUCTOR

  DESIGN_hrandom(const DESIGN_hrandom & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN_hrandom & operator=(const DESIGN_hrandom & m);

  // virtual functions

  void init_data(const datamatrix & dm,const datamatrix & iv);

  void compute_penalty(void);

  void compute_XtransposedWX_XtransposedWres(datamatrix & partres, double l);

  void compute_XtransposedWres(datamatrix & partres, double l);

  void compute_precision(double l);

  // DESTRUCTOR

  ~DESIGN_hrandom() {}

  };



} // end: namespace MCMC

#endif


