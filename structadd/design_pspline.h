
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (DESIGNpsplineINCLUDED)

#define DESIGNpsplineINCLUDED

#include<deque>
#include<design.h>

namespace MCMC
{


//------------------------------------------------------------------------------
//--------------------------- CLASS: DESIGN_pspline ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN_pspline : public DESIGN
  {

  protected:

  deque<double> knot;          // Vektor der Knoten (sichtbare und unsichtbare)

  vector<double> weightK;      // weights to compute penalty matrix K

  // FUNCTION: bspline
  // TASK: computes B-splines at position x

  datamatrix bspline(const double & x);

  // FUNCTION: make_Bspline
  // TASK: computes knot, Zout and index_Zout

  void make_Bspline(void);


  public:

  long nrknots;                     // Anzahl der (sichtbaren) Knoten
  long degree;                      // Grad des Splines
  long difforder;                   // Differenzenordnung (1,2,3)


//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN_pspline(void);

  // CONSTRUCTOR

  DESIGN_pspline(const datamatrix & dm, const datamatrix & iv,
             DISTR * dp,vector<ST::string> & op);

  // COPY CONSTRUCTOR

  DESIGN_pspline(const DESIGN_pspline & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN_pspline & operator=(const DESIGN_pspline & m);

  // virtual functions

  void init_data(const datamatrix & dm, const datamatrix & iv);

  void compute_penalty(void);

  void compute_XtransposedWX(datamatrix & partres);

  void compute_XtransposedWres(datamatrix & partres, double l);

  void compute_XtransposedWX_XtransposedWres(datamatrix & partres, double l);

  void compute_precision(double l);

  void read_options(vector<ST::string> & op);

  // DESTRUCTOR

  ~DESIGN_pspline() {}

  };


} // end: namespace MCMC

#endif

