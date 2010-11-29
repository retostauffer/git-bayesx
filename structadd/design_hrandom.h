
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (DESIGNhrandomINCLUDED)

#define DESIGNhrandomINCLUDED

#include"statmat.h"
#include"sparsemat.h"

#include"Random.h"
#include"envmatrix_penalty.h"
#include"../values.h"
#include<fstream>
#include<vector>
#include<bitset>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include"distr.h"
#include"design.h"
#include<cmath>


namespace MCMC
{


//------------------------------------------------------------------------------
//--------------------------- CLASS: DESIGN_hrandom ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN_hrandom : public DESIGN
  {

  // beta contains linpredRE+random effect


  protected:

  DISTR * likep_RE;

  public:


//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN_hrandom(void);

  // CONSTRUCTOR 1
  // Spatial covariates

  DESIGN_hrandom(const datamatrix & dm, const datamatrix & iv,
             GENERAL_OPTIONS * o,DISTR * dp,FC_linear * fcl, DISTR * dp_RE,
             vector<ST::string> & op,
             vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  DESIGN_hrandom(const DESIGN_hrandom & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN_hrandom & operator=(const DESIGN_hrandom & m);

  // virtual functions

  void init_data(const datamatrix & dm,const datamatrix & iv);

  void compute_penalty(void);

  void compute_XtransposedWres(datamatrix & partres, double l);

  void compute_precision(double l);

  void compute_meaneffect(DISTR * level1_likep,double & meaneffect,
                          datamatrix & beta,datamatrix & meaneffectbeta,
                          bool computemeaneffect, double meaneffectconstant);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  void outoptions(GENERAL_OPTIONS * op);

  // DESTRUCTOR

  ~DESIGN_hrandom() {}

  };



} // end: namespace MCMC

#endif


