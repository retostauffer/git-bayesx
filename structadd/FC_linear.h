
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FClinearINCLUDED)

#define FClinearINCLUDED

#include"../values.h"
#include<fstream>
#include"GENERAL_OPTIONS.h"
#include"distr.h"
#include"clstring.h"
#include"FC.h"
#include"MASTER_obj.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_linear ---------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_linear  : public FC
  {

  protected:

  int constposition;
  bool center;

  MASTER_OBJ * masterp;

  void compute_meaneffect_design(void);

  bool initialize;
  bool IWLS;

  DISTR * likep;                             // Pointer to DISTR obejct
  datamatrix design;                         // Designmatrix
  datamatrix mean_designcols;
  vector<datamatrix> designhelp;             // help vector for constructing the
                                             // designmatrix
  datamatrix meaneffectdesign;
  vector<ST::string> datanames;              // names of covariates

  datamatrix Xt;                             // transposed designmatrix
  datamatrix XWX;
  datamatrix XWXold;
  datamatrix XWXroot;
  datamatrix residual;
  datamatrix Xtresidual;

  datamatrix betam;
  datamatrix help;
  datamatrix betaold;
  datamatrix betadiff;
  datamatrix mode;
  datamatrix proposal;

  datamatrix linold;
  datamatrix linnew;
  datamatrix linmode;
  datamatrix diff;
  datamatrix * linoldp;
  datamatrix * linnewp;

  void find_const(datamatrix & design);

  void create_matrices(void);

  // FUNCTION: compute_XWX
  // TASK: computes XWX on the basis of the current working weight and stores
  //       the result in r

  void compute_XWX(datamatrix & r);
  void compute_Wpartres(datamatrix & linpred);
  double compute_XtWpartres(double & mo);

  void add_linpred(datamatrix & l);

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_linear(void);

  // CONSTRUCTOR

  FC_linear(MASTER_OBJ * mp, GENERAL_OPTIONS * o,DISTR * lp, datamatrix & d,
            vector<ST::string> & vn, const ST::string & t,
           const ST::string & fp,bool cent);

  // COPY CONSTRUCTOR

  FC_linear(const FC_linear & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_linear & operator=(const FC_linear & m);

  // DESTRUCTOR

  ~FC_linear()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  void update_gaussian(void);
  void update_IWLS(void);

  // FUNCTION: posteriormode

  bool posteriormode(void);

  // FUNCTION: outoptions

  void outoptions(void);

  // FUNCTION: outresults

  void outresults(ofstream & out_stata,ofstream & out_R,
                  const ST::string & pathresults);

  void compute_autocorr_all(const ST::string & path,
                              unsigned lag, ofstream & outg) const;


  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  // FUNCTION: reset

  void reset(void);

  // FUNCTION: add_variable

  int add_variable(const datamatrix & d,ST::string & name);

  };


//------------------------------------------------------------------------------
//------------------------ CLASS: FC_linear_pen --------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_linear_pen  : public FC_linear
  {

  protected:


  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_linear_pen(void);

  // CONSTRUCTOR

  FC_linear_pen(MASTER_OBJ * mp, GENERAL_OPTIONS * o,DISTR * lp, datamatrix & d,
            vector<ST::string> & vn, const ST::string & t,
           const ST::string & fp,bool cent);

  // COPY CONSTRUCTOR

  FC_linear_pen(const FC_linear_pen & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_linear_pen & operator=(const FC_linear_pen & m);

  // DESTRUCTOR

  ~FC_linear_pen()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  // FUNCTION: posteriormode

  bool posteriormode(void);

  // FUNCTION: outoptions

  void outoptions(void);

  // FUNCTION: outresults

  void outresults(ofstream & out_stata,ofstream & out_R,
                  const ST::string & pathresults);

  // FUNCTION: reset

  void reset(void);

  };



} // end: namespace MCMC

#endif


