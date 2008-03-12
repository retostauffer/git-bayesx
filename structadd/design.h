
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
                kriging
                };


//------------------------------------------------------------------------------
//--------------------------- CLASS: DESIGN ------------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN
  {

  protected:


  public:

  datamatrix data;
  statmatrix<int> index_data;

  datamatrix Z;
  
  datamatrix Zout;
  statmatrix<int> index_Zout;

  // Penalty Matrix

  envmatdouble K;
  double rankK;

  // X'WX

  envmatdouble XWX;


  ttype type;



//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN(void);

  // CONSTRUCTOR

  DESIGN(const datamatrix & dm);

  // COPY CONSTRUCTOR

  DESIGN(const DESIGN & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN & operator=(const DESIGN & m);

  virtual void compute_design(void);

  virtual void compute_penalty(void);

  virtual void compute_XtransposedWX(void);

  void compute_precision(double l);

  void compute_XtransposedW(void);

  // DESTRUCTOR

  ~DESIGN() {}


  };



//------------------------------------------------------------------------------
//--------------------------- CLASS: DESIGN_mrf --------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN_mrf : public DESIGN
  {

  MAP::map ma;
  vector<int> posbeg;
  vector<int> posend;
  vector<ST::string> effectvalues;

  protected:



  public:


//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN_mrf(void);

  // CONSTRUCTOR 1
  // Spatial covariates

  DESIGN_mrf(const datamatrix & dm,const MAP::map & m);

  // COPY CONSTRUCTOR

  DESIGN_mrf(const DESIGN_mrf & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN_mrf & operator=(const DESIGN_mrf & m);

  void compute_design(void);

  void compute_penalty(void);

  // DESTRUCTOR

  ~DESIGN_mrf() {}

  };



} // end: namespace MCMC

#endif


