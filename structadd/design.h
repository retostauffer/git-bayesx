
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
//--------------------------- CLASS: DESIGN ------------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN
  {

  protected:

  // FUNCTION: make_index
  // TASK: sorts the data,
  //       initializes index_data,
  //       posbeg, posend, effectvalues

  void make_index(const datamatrix & dm, const datamatrix & iv);



  DISTR * likep;                             // Pointer to DISTR obejct

  //----------------------------------------------------------------------------

  vector< vector<double> > ZoutT;            // Nonzero Elements of Z'
  vector< vector<int> > index_ZoutT;         // Columns of nonzero elements of
                                             // Z'

  void compute_Zout_transposed(void);        // Computes Z', i.e. ZoutT from
                                             // Zout

  int consecutive_ZoutT;                     // -1 = not tested
                                             //  0 = not consecutive
                                             //  1 = consecutive

  bool check_ZoutT_consecutive(void);

  //----------------------------------------------------------------------------

  public:

  // Variables determined by function init_data

  datamatrix data;                           // data matrix
  datamatrix intvar;                         // interaction variable for
                                             // varying coefficients


  statmatrix<int> index_data;                // index for sort of data
  vector<ST::string> datanames;              // names of covariates


  vector<ST::string> effectvalues;           // values of the different
                                             // covariates

//------------------------------------------------------------------------------
  datamatrix Zout;                           // Design matrix (only non null
                                             // elements for output of results
  statmatrix<int> index_Zout;                // stores the columns of the
                                             // non null elements of Zout

  vector<int> posbeg;                        // begin and end of equal covariate
                                             // values in data
  vector<int> posend;

  int consecutive;                           // -1 = not tested
                                             //  0 = not consecutive
                                             //  1 = consecutive

  bool check_Zout_consecutive(void);

//------------------------------------------------------------------------------


  statmatrix<double *> responsep;            // matrix of pointers to response
                                             // observations
  statmatrix<double *> workingweightp;

  statmatrix<double *> linpredp1;
  statmatrix<double *> linpredp2;    

  unsigned nrpar;                            // number of parameters

  // ---------------------------------------------------------------------------

  // Variables determined by function compute_penalty

  envmatdouble K;                            // Penalty Matrix
  double rankK;

  // ---------------------------------------------------------------------------

  // Variables determined by function compute_precision

  envmatdouble precision;                    // precision matrix
  bool precisiondeclared;                    // true if precision is already
                                             // defined

  // ---------------------------------------------------------------------------

  // Variables determined by function  compute_XtransposedWX_XtransposedWres
  // and compute_XtransposedWres

  datamatrix Wsum;

  envmatdouble XWX;                          // X'WX
  bool XWXdeclared;                          // is true if X'WX is already
                                             // defined (i.e. memory allocated
                                             // etc.)

  datamatrix XWres;                          // X'W(y-eta)
  bool XWresdeclared;

  // ---------------------------------------------------------------------------

  ttype type;


  //----------------------- CONSTRUCTORS, DESTRUCTOR ---------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN(void);

  // CONSTRUCTOR

  DESIGN(DISTR * dp);

  // COPY CONSTRUCTOR

  DESIGN(const DESIGN & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN & operator=(const DESIGN & m);

  //----------------------------------------------------------------------------


  void compute_f(datamatrix & beta,datamatrix & f);

  void update_linpred(datamatrix & f,bool add);

  void compute_partres(datamatrix & res,datamatrix & f);


  void make_pointerindex(void);

  // ------------------------- VIRTUAL FUNCTIONS -------------------------------

  // FUNCTION: init_data
  // TASK: sorts the data such that the precision has minimum envelope
  //       computes index_data
  //       computes partres_pindex
  //       computes Zout, posbeg, posend
  //       computes nrpar
  //       computes effectvalues
  //       initializes datanames

  virtual void init_data(const datamatrix & dm, const datamatrix & iv);

  // FUNCTION: compute_penalty
  // TASK: computes the penalty matrix and determines rankK

  virtual void compute_penalty(void);

  // FUNCTION: compute_XtransposedWX_XtransposedWres
  // TASK: computes XWX and XWres, res is the partial residual

  virtual void compute_XtransposedWX_XtransposedWres(datamatrix & partres, double l);

  double compute_ZtZ(unsigned & i, unsigned & j);

  // FUNCTION: computes XWres
  // TASK: computes XWres, res is the partial residual

  virtual void compute_XtransposedWres(datamatrix & partres, double l);

  virtual void compute_precision(double l);

  // --------------------- END: VIRTUAL FUNCTIONS ------------------------------

  // DESTRUCTOR

  ~DESIGN() {}


  };






} // end: namespace MCMC

#endif


