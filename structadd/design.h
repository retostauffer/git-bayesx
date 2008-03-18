
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
                kriging
                };


//------------------------------------------------------------------------------
//--------------------------- CLASS: DESIGN ------------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN
  {

  protected:


  DISTR * likep;                             // Pointer to DISTR obejct

  public:

  datamatrix data;                           // data matrix
  datamatrix data2;                          // vor varying coefficients
                                             //
  statmatrix<int> index_data;                // index for index sort of data
  vector<ST::string> datanames;              // names of covariates

  vector<ST::string> effectvalues;           // values of the different
                                             // covariates

  datamatrix Z;

  datamatrix Zout;                           // Design matrix (only non null
                                             // elements for output of results
  statmatrix<int> index_Zout;                // stores the columns of the
                                             // non null elements of Zout

  vector<int> posbeg;                        // begin and end of equal covariate 
                                             // values in data
  vector<int> posend;                        



  unsigned nrpar;                            // number of parameters                                    



  envmatdouble K;                            // Penalty Matrix
  double rankK;

  envmatdouble XWX;                          // X'WX
  bool XWXdeclared;                          // is true if X'WX is already
                                             // defined (i.e. memory allocated
                                             // etc.)


  envmatdouble precision;                    // precision matrix
  bool precisiondeclared;                    // true if precision is already
                                             // defined                         

  datamatrix XWres;                          // X'W(y-eta)
  bool XWresdeclared;

  ttype type;



//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN(void);

  // CONSTRUCTOR

  DESIGN(const datamatrix & dm, DISTR * dp);

  // COPY CONSTRUCTOR

  DESIGN(const DESIGN & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN & operator=(const DESIGN & m);


  // FUNCTION: compute_design
  // TASK:
  //       nrpar = number of parameters is defined

  virtual void compute_design(void);

  // FUNCTION: compute_penalty
  // TASK: computes the penalty matrix and determines rankK

  virtual void compute_penalty(void);

  virtual void compute_XtransposedWX_XtransposedWres(const datamatrix & res);

  virtual void compute_XtransposedWres(const datamatrix & res);

  void compute_precision(double l);  

  void compute_f(datamatrix & beta,datamatrix & f);

  void update_linpred(datamatrix & f,bool add);

  double compute_quadform(const datamatrix & beta);

  // DESTRUCTOR

  ~DESIGN() {}


  };



//------------------------------------------------------------------------------
//--------------------------- CLASS: DESIGN_mrf --------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN_mrf : public DESIGN
  {

  MAP::map ma;

  protected:



  public:


//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN_mrf(void);

  // CONSTRUCTOR 1
  // Spatial covariates

  DESIGN_mrf(const datamatrix & dm, DISTR * dp, const MAP::map & m);

  // COPY CONSTRUCTOR

  DESIGN_mrf(const DESIGN_mrf & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN_mrf & operator=(const DESIGN_mrf & m);

  void compute_design(void);

  void compute_penalty(void);

  void compute_XtransposedWX_XtransposedWres(const datamatrix & res);

  void compute_XtransposedWres(const datamatrix & res);  

  // DESTRUCTOR

  ~DESIGN_mrf() {}

  };



} // end: namespace MCMC

#endif


