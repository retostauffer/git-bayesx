
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE  __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (DESIGNINCLUDED)

#define DESIGNINCLUDED

#include"statmat.h"
#include"sparsemat.h"

#include"Random.h"
#include"envmatrix_penalty.h"
#include"../values.h"
#include<fstream>
#include<vector>
#include<bitset>
#include"GENERAL_OPTIONS.h"
#include"FC_linear.h"
#include"clstring.h"
#include"distr.h"
#include<cmath>


namespace MCMC
{


enum ttype2 {   Rw1,
                Rw2,
                Rw3,
                Mrf,
                Hrandom
                };


enum effecttype2 {
                Function,
                Varcoefftotal
   };

enum centerm {cmean,nullspace};

//------------------------------------------------------------------------------
//--------------------------- CLASS: DESIGN ------------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN
  {

  protected:

  // FUNCTION: make_index
  // TASK: takes index_data (computed elsewhere)
  //       creates sorted data ntvar, data2

  void make_data(const datamatrix & dm,const datamatrix & iv);

  // FUNCTION: make_index
  // TASK: sorts the data,
  //       creates sorted intvar, data2
  //       initializes index_data,
  //       posbeg, posend, effectvalues

  void make_index(const datamatrix & dm, const datamatrix & iv);


  unsigned compute_modecategorie(void);  


  DISTR * likep;                             // Pointer to DISTR obejct

  //----------------------------------------------------------------------------

  vector< vector<double> > ZoutT;            // Nonzero Elements of Z'
  vector< vector<int> > index_ZoutT;         // Columns of nonzero elements of
                                             // Z'
  int consecutive_ZoutT;                     // -1 = not tested
                                             //  0 = not consecutive
                                             //  1 = consecutive

  void compute_Zout_transposed(void);        // Computes Z', i.e. ZoutT from
                                             // Zout

  bool check_ZoutT_consecutive(void);        // checks if non zero elements are
                                             // consecutive

  //----------------------------------------------------------------------------

  public:

  bool changingdesign;

  // Variables determined by function init_data

  datamatrix data;                           // data matrix
  datamatrix intvar;                         // interaction variable for
                                             // varying coefficients
  datamatrix intvar2;                        // intvar^2 for varying coefficients

  statmatrix<int> index_data;                // index for sort of data
  vector<ST::string> datanames;              // names of covariates


  vector<ST::string> effectvalues;           // values of the different
                                             // covariates

  unsigned meaneffectnr;
  unsigned meaneffectnr_intvar;
  double meaneffectintvar;

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

  bool identity;                             // true if identity matrix

  bool check_Zout_consecutive(void);

  //----------------------------------------------------------------------------

  statmatrix<double *> responsep;            // matrix of pointers to
                                             // response observations
  statmatrix<double *> weightp;              // matrix of pointers to
                                             // weights

  statmatrix<double *> workingresponsep;     // matrix of pointers to working
                                             // response observations
  statmatrix<double *> workingweightp;       // matrix of pointers to working
                                             // weights

  statmatrix<double *> linpredp1;            // matrix of pointers to linpred1
  statmatrix<double *> linpredp2;            // matrix of pointers to linpred2

  // FUNCTION: make_pointerindex
  // TASK: computes pointer matrices responsep, workingweightp,linpredp1,
  //       linpredp2

  void make_pointerindex(void);

  //----------------------------------------------------------------------------

  unsigned nrpar;                            // number of parameters

  // --------------------------- for center ------------------------------------

  bool center;
  centerm centermethod;


  // for nullspace centering
  datamatrix basisNull;                     // contains a basis of the null
                                            // space of the penalty K
  vector<datamatrix> basisNullt;            // contains the transposed of
                                            // basisNull

  FC_linear * FClinearp;                    // Pointer to linear effects
  int position_lin;                         // position in the designmatrix
                                            // of linear effects
  datamatrix designlinear;                  // designmatrix linear effects
  // end for nullpsace centering


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
  bool XWresdeclared;                        // is true if X'Wres is already
                                             // defined


  // ---------------------------------------------------------------------------

  ttype2 type;                                // Term type


  //----------------------- CONSTRUCTORS, DESTRUCTOR ---------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN(void);

  // CONSTRUCTOR

  DESIGN(DISTR * dp,FC_linear * fcl);

  // COPY CONSTRUCTOR

  DESIGN(const DESIGN & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN & operator=(const DESIGN & m);

  //----------------------------------------------------------------------------

  // FUNCTION: compute_f
  // TASK: compute Zout*beta, i.e. the estimated/current function evaluated at
  //       the different observations in data

  void compute_f(datamatrix & beta,datamatrix & betalin,
                       datamatrix & f, datamatrix & ftot);

  // FUNCTION: compute_effect
  // TASK: computes the effect vector

  void compute_effect(datamatrix & effect,datamatrix & f,
                      effecttype2 et = Function);

  void set_intvar(datamatrix & iv, double add=0);

  // FUNCTION: update_linpred
  // TASK: updates the predictor based on the current function f

  void update_linpred(datamatrix & f);

  // FUNCTION: compute_partres
  // TASK: computes

  void compute_partres(datamatrix & res,datamatrix & f);

  void compute_partres(int begin,int end,double & res, double & f);

//  void compute_partres_nopred(datamatrix & res, datamatrix & f);

  double compute_ZtZ(unsigned & i, unsigned & j);

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

  // FUNCTION: compute_basisNull
  // TASK: computes the basis of the null space of the penalty matrix

  virtual void compute_basisNull(void);

  // FUNCTION: computes XWres
  // TASK: computes XWres, res is the partial residual

  virtual void compute_XtransposedWres(datamatrix & partres, double l);

  // FUNCTION: compute_XtransposedWX_XtransposedWres
  // TASK: computes XWX and XWres, res is the partial residual

  virtual void compute_XtransposedWX(void);

  // FUNCTION: compute_XtransposedWX_XtransposedWres
  // TASK: computes XWX and XWres, res is the partial residual

  virtual void compute_XtransposedWX_XtransposedWres(datamatrix & partres, double l);


  virtual void compute_meaneffect(DISTR * level1_likep,double & meaneffect,
                                datamatrix & beta,datamatrix & meaneffectbeta,
                                bool computemeaneffect);


  virtual void compute_precision(double l);

  // FUNCTION: read_options
  // TASK: reads options and initializes varnames stored in datanames

  virtual void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  virtual void outoptions(GENERAL_OPTIONS * op);

  // --------------------- END: VIRTUAL FUNCTIONS ------------------------------

  // DESTRUCTOR

  ~DESIGN() {}


  };






} // end: namespace MCMC

#endif


