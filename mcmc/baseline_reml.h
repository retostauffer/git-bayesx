//---------------------------------------------------------------------------
#ifndef baseline_remlH
#define baseline_remlH

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif

#include<mcmc_pspline.h>
#include<spline_basis.h>
#include<vector>

namespace MCMC
{

//------------------------------------------------------------------------------
//---------------------------- class: baseline_reml ----------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE baseline_reml : public spline_basis
  {

  protected:

  double tstep;
  unsigned tgrid;
  datamatrix tsteps;

  knotpos gridpos;
  unsigned nrquant;
  unsigned nrbetween;

  vector<unsigned>tstart;
  vector<unsigned>tend;

  datamatrix t_X;
  datamatrix t_Z;

  datamatrix interact_var;

  datamatrix tvalues;

  public:

  // DEFAULT CONSTRUCTOR

  baseline_reml(void) : spline_basis()
    {
    }

  // CONSTRUCTOR 1

  baseline_reml(MCMCoptions * o, const datamatrix & d, const datamatrix & lo,
               const unsigned & nrk, const unsigned & degr, const unsigned & tgr,
               const unsigned & nrq, const unsigned & nrb, const knotpos & kp,
               const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const knotpos & gp);

  // CONSTRUCTOR 2 (VCM)

  baseline_reml(MCMCoptions * o,const datamatrix & d1,
                      const datamatrix & d2, const unsigned & nrk,
                      const unsigned & degr, const unsigned & tgr,
                      const knotpos & kp, const fieldtype & ft,
                      const ST::string & ti, const ST::string & fp,
                      const ST::string & pres, const double & l,
                      const double & sl);

  // COPY CONSTRUCTOR

  baseline_reml(const baseline_reml & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const baseline_reml & operator=(const baseline_reml & fc);

  // DESTRUCTOR

  ~baseline_reml() {}

  void createreml(datamatrix & X,datamatrix & Z,const unsigned & Xpos,
                  const unsigned & Zpos);

  void multDG(datamatrix & res, const datamatrix & b);

  void initialize_baseline(unsigned j, datamatrix & tx, datamatrix & tz,
               vector<unsigned> & ts, vector<unsigned> & te, datamatrix & iv,
               statmatrix<double> & steps, statmatrix<int> & ind);

  void outoptionsreml();

  void init_name(const ST::string & na);

  unsigned & get_tgrid(void)
    {
    return tgrid;
    }

  };

}   // end: namespace MCMC

//---------------------------------------------------------------------------
#endif
