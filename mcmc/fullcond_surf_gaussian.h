// Date: 4.12.99

#if !defined (MCMCsurfgaussian_INCLUDED)

#define MCMCsurfgaussian_INCLUDED

#include<mcmc_nonpbasis.h>
#include<fullcond_nonp_gaussian.h>

namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_surf_gaussian ------------------------
//------------------------------------------------------------------------------

class FULLCOND_surf_gaussian : public FULLCOND_nonp_basis
  {


  protected:

  FULLCOND fchelp;
  ST::string fchelprespath;

  bool singleblock;

  bool centertotal;

  vector<SparseMatrix> Kleft;
  vector<SparseMatrix> Kright;
  vector<bandmatdouble> Kab;
  unsigned sizeK1;
  unsigned sizeK2;
  unsigned sizeK;

  vector<bandmatdouble> XX;

  bandmatdouble prec;

  datamatrix standnormal;

  datamatrix mu;
  datamatrix muy;

  datamatrix betahelp;
  datamatrix betahelp2;

  double lambda;
  double lambdaold;

  datamatrix sumx1;
  datamatrix sumx2;

  FULLCOND_nonp_gaussian * mainp1;
  FULLCOND_nonp_gaussian * mainp2;

  SparseMatrix Krw1(const vector<double> & weight);
  SparseMatrix Krw2(const vector<double> & weight);

  SparseMatrix Kmrflinear(const unsigned & nr1,
                          const unsigned & nr2);


  void FULLCOND_surf_gaussian::make_moddata(const datamatrix & moddata1,
                                            const datamatrix & moddata2);


  // FUNCTION: make_categories
  // TASK: devides the data in moddata into categories, maximal number of
  //       categories is 'maxnrint'
  //       Initialices 'index', 'posbeg', 'posend', 'weight' and
  //                   'effectvalues'
  // ADDITIONAL NOTE: Implementation is independent of the type of MRF

  datamatrix make_categories(const datamatrix & moddata,
                             vector<ST::string> & effvalues);


  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_surf_gaussian(void) : FULLCOND_nonp_basis()
    {
    }

  // CONSTRUCTOR 1
  // o    : pointer to MCMCoptions object
  // dp   : pointer to distribution object
  // fp   : file where sampled parameters are stored
  // pres : file where results are stored
  // order in beta:
  // first value d1, first d2
  // first value d1, second d2
  // ...

  FULLCOND_surf_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                         const datamatrix & d1,const datamatrix & d2,
                         const double & a, const double & b,
                         const fieldtype & ft, const ST::string & ti,
                         const ST::string & fp,
                         const ST::string & pres,const unsigned & c,
                         bool sb = true);

  void init_maineffects(FULLCOND_nonp_gaussian * mp1,
                        FULLCOND_nonp_gaussian * mp2,
                        const ST::string & pnt,const ST::string & prt);

  // COPY CONSTRUCTOR

  FULLCOND_surf_gaussian(const FULLCOND_surf_gaussian & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_surf_gaussian & operator=(const FULLCOND_surf_gaussian & fc);

  void update(void);

  bool posteriormode(void);

  void outresults(void);

  void outoptions(void);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    FULLCOND_nonp_basis::reset();
    }

  // FUNCTION: predict (virtual)
  // TASK: predicts the mean for a new observation Xnew

  void predict(const datamatrix & newX, datamatrix & linpred)
    {
    }


  double compute_quadform(void)
    {
    return Ksp.compute_quadform(beta,0);
    }

  // DESTRUCTOR

  ~FULLCOND_surf_gaussian() {}

  };


} // end: namespace MCMC

#endif
