//---------------------------------------------------------------------------
#ifndef baselineH
#define baselineH

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif


#include<cox.h>
#include<mcmc_pspline.h>
#include<spline_basis.h>
#include<vector>

namespace MCMC
{


//------------------------------------------------------------------------------
//---------------------------- class: pspline_baseline -------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE pspline_baseline : public FULLCOND_pspline
  {

  protected:

   bool begin0;
   datamatrix int_knots;
   datamatrix int_D;
   MCMC::bsplinemat testmat;
   vector<MCMC::bsplinemat> gaussmat;
   vector<pspline_baseline*> baselinep;
   datamatrix zi;
   unsigned gauss_n;
   datamatrix coeff;
   datamatrix z_vc;
   datamatrix zi_ges;
   datamatrix beg_i;
   statmatrix<int> zi_index;
   statmatrix<int> ges_index;
   datamatrix spline_ges;
   datamatrix spline_ges2;
   datamatrix spline_zi;
   datamatrix gaussspline;
   datamatrix int_ti_help;
   bool vc_dummy1;

  public:


  // DEFAULT CONSTRUCTOR

  pspline_baseline(void) : FULLCOND_pspline()
    {
    }

  // CONSTRUCTOR 1

  pspline_baseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & d,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l,const unsigned & minb,const unsigned & maxb,
                    const fieldtype & ft,const ST::string & ti,
                    const ST::string & fp, const ST::string & pres,
                    const int & gs,const unsigned & c,const datamatrix & anfang);


// CONSTRUCTOR 2 (für zeitabhängige Effekte)

  pspline_baseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & time, const datamatrix & z,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l,const unsigned & minb,const unsigned & maxb,
                    const fieldtype & ft,const ST::string & ti,
                    const ST::string & fp, const ST::string & pres,
                    const int & gs,const unsigned & c,const datamatrix & anfang);


  // COPY CONSTRUCTOR

  pspline_baseline(const pspline_baseline & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const pspline_baseline & operator=(const pspline_baseline & fc);


  void update(void);


  void outoptions(void);
//  void outresults(void);
  void compute_int_ti(const datamatrix & b);
  void compute_int_ti_linear(const double & b);
  void compute_int_ti(unsigned beg);
  void compute_int_ti_vc_di0(const vector<double *>,const vector<double *>);
  void compute_int_ti_vc_di(const int,const vector<double *>,const vector<double *>);
  void compute_int_gauss(void);
  void compute_int_gauss_DIC(void);
  void update_baseline(void);
  void compute_int_ti_mean(void);
  void set_baselinep(vector<pspline_baseline*> bp)
    {
    baselinep=bp;
    }

  double * get_int_D(void)
    {
    return int_D.getV();
    }

  double * get_z_vc(void)
    {
    return z_vc.getV();
    }

  datamatrix get_z_vc_np(void)
    {
    return z_vc;
    }

  double * get_spline_zi(void)
    {
    multBS(spline_zi,beta);
    return spline_zi.getV();
    }

  double * get_spline_ges(void)
    {
    testmat.mult_index(spline_ges,beta);
    return spline_ges.getV();
    }

  double * get_spline_ges2(void)
    {
    testmat.mult_index(spline_ges2,beta);
    return spline_ges2.getV();
    }

  double * get_spline_ges_mean(void)
    {
    testmat.mult_index(spline_ges,betamean);
    return spline_ges.getV();
    }

  double * get_gaussspline(void);


  double * get_gaussspline_mean(void);


  double * get_betamean(void)
    {
    return betamean.getV();
    }

  double * get_beta(void)
    {
    return beta.getV();
    }

  double * get_spline_zi_mean(void)
    {
    multBS(spline_zi,betamean);
    return spline_zi.getV();
    }

  void set_fcconst(FULLCOND_const * fcc)
    {
    fcconst = fcc;
    }




  // DESTRUCTOR

  ~pspline_baseline() {}

  };



}   // end: namespace MCMC


//---------------------------------------------------------------------------
#endif
