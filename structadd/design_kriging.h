
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPingE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (DESIGNkrigingINCLUDED)

#define DESIGNkrigingINCLUDED

#include"statmat.h"
#include"design.h"
#include"Random.h"
#include"envmatrix_penalty.h"
#include"../values.h"
#include<fstream>
#include<vector>
#include<bitset>
#include"clstring.h"
#include<cmath>


namespace MCMC
{


//------------------------------------------------------------------------------
//--------------------------- CLASS: DESIGN_kriging ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN_kriging : public DESIGN
  {

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  double compute_matern(double & nu,double & r);

  void compute_tildeZ(void);  

  protected:

  double rho;
  double maxdist;
  double nu;

  vector<double> xknots;              // x-und y-Koordinaten der Knoten
  vector<double> yknots;

  vector<double> xvalues;             // unterschiedliche Werte der Kovariablen
  vector<double> yvalues;

  void compute_knots(const vector<double> & xvals,
                     const vector<double> & yvals,
                     unsigned nrknots,double p,double q,
                     vector<double> & xknots,
                     vector<double> & yknots);

  long nrknots;

  datamatrix tildeZ_t;

  datamatrix XWXfull;
  datamatrix WsumtildeZ;

  datamatrix Kfull;

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN_kriging(void);

  // CONSTRUCTOR 1
  // Spatial covariates

  DESIGN_kriging(const datamatrix & dm, const datamatrix & iv,
                 GENERAL_OPTIONS * o,DISTR * dp,FC_linear * fcl,
                 vector<ST::string> & op,
                 vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  DESIGN_kriging(const DESIGN_kriging & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN_kriging & operator=(const DESIGN_kriging & m);

  // VIRTUAL FUNCTIONS

  void init_data(const datamatrix & dm,const datamatrix & iv);

  void compute_penalty(void);

  void compute_XtransposedWres(datamatrix & partres, double l);

  void compute_XtransposedWX(void);

  void compute_precision(double l);

  void outoptions(GENERAL_OPTIONS * op);

  void compute_orthogonaldecomp(void);

  double penalty_compute_quadform(datamatrix & beta);

  // DESTRUCTOR

  ~DESIGN_kriging() {}

  };


} // end: namespace MCMC

#endif


