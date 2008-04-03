
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

  unsigned nrknots;                     // Anzahl der (sichtbaren) Knoten
  unsigned degree;                      // Grad des Splines



//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN_pspline(void);

  // CONSTRUCTOR

  DESIGN_pspline(const datamatrix & dm, const datamatrix & iv,
             DISTR * dp);

  // COPY CONSTRUCTOR

  DESIGN_pspline(const DESIGN_pspline & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN_pspline & operator=(const DESIGN_pspline & m);

  // virtual functions

  void init_data(const datamatrix & dm, const datamatrix & iv);

  void compute_penalty(void);

  void compute_XtransposedWX_XtransposedWres(double l);

  void compute_XtransposedWres(double l);

  void compute_precision(double l);

  // DESTRUCTOR

  ~DESIGN_pspline() {}

  };


} // end: namespace MCMC

#endif


