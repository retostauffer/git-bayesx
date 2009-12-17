
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (DESIGNpsplineINCLUDED)

#define DESIGNpsplineINCLUDED

#include<deque>
#include<design.h>

using std::deque;

namespace MCMC
{


//------------------------------------------------------------------------------
//--------------------------- CLASS: DESIGN_pspline ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN_pspline : public DESIGN
  {

  protected:

  bool ccov;                    // true if covariate should be centered

  bool multeffect;

  deque<double> knot;          // Vektor der Knoten (sichtbare und unsichtbare)

  vector<double> weightK;      // weights to compute penalty matrix K

  // FUNCTION: bspline
  // TASK: computes B-splines at position x

  datamatrix bspline(const double & x);

  // FUNCTION: make_Bspline
  // TASK: computes knot, Zout and index_Zout

  void make_Bspline(void);

  void compute_betaweight(datamatrix & betaweight);

  public:

  long nrknots;                     // Anzahl der (sichtbaren) Knoten
  long degree;                      // Grad des Splines
  long difforder;                   // Differenzenordnung (1,2,3)
  double round;
  double binning;


//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN_pspline(void);

  // CONSTRUCTOR

  DESIGN_pspline(datamatrix & dm, datamatrix & iv,
             DISTR * dp,FC_linear * fcl, vector<ST::string> & op,
             vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  DESIGN_pspline(const DESIGN_pspline & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN_pspline & operator=(const DESIGN_pspline & m);

  // virtual functions

  void compute_penalty(void);

  void compute_basisNull(void);

  void compute_precision(double l);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  void outoptions(GENERAL_OPTIONS * op);

  // DESTRUCTOR

  ~DESIGN_pspline() {}

  };


} // end: namespace MCMC

#endif


