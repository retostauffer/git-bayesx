
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (VARIANCEPEN_VECTOR_INCLUDED)
#define VARIANCEPEN_VECTOR_INCLUDED

#include"FC_linear.h"

namespace MCMC
{


class __EXPORT_TYPE FC_variance_pen_vector : public FC
  {

  protected:

  vector<double> tau;               // tau^2
  vector<double> lambda;            // Inverse Varianceparameter: lambda=1/tau^2

  bool update_sigma2;

  FC_linear_pen * Cp;

  DISTR * distrp;

  FC fc_shrinkage;
  vector<bool> shrinkagefix;          //  Shrinkageparameter fix
  vector<double> a_shrinkagegamma;    //  Hyperparameter for Shrinkageparameter
  vector<double> b_shrinkagegamma;    //  Hyperparameter for Shrinkageparameter
  vector<double> shrinkagestart;


  double lassosum;                  //  sum(beta^2/tau^2)
  double ridgesum;                  //  sum(beta^2/tau^2)

  bool is_ridge;          //  The Components indicates if "true" the L2-penalty
                         //  and if "false" the L1-penalty is used

  void outresults_shrinkage(void);  //  Function to write results to output window and files

  public:

  //____________________________________________________________________________
  //
  // DEFAULT CONSTRUCTOR
  //____________________________________________________________________________

  FC_variance_pen_vector(void) : FC()
    {
    }


  //____________________________________________________________________________
  //
  // CONSTRUCTOR
  //____________________________________________________________________________

  FC_variance_pen_vector(MASTER_OBJ * mp,GENERAL_OPTIONS * o, FC_linear_pen * p,
                         DISTR * d,const ST::string & ti,
                         const ST::string & fp, bool isr,
                         vector<ST::string> & op,
                         vector<ST::string> & vn);

  //____________________________________________________________________________
  //
  // COPY CONSTRUCTOR
  //____________________________________________________________________________

  FC_variance_pen_vector(const FC_variance_pen_vector & t);


  //____________________________________________________________________________
  //
  // OVERLOADED ASSIGNMENT OPERATOR
  //____________________________________________________________________________

  const FC_variance_pen_vector & operator=(const FC_variance_pen_vector & t);


  void add_variable(datamatrix & x, double la, double shrink,
                    bool sfix, double a, double b);

  //____________________________________________________________________________
  //
  // OVERLOADED ASSIGNMENT OPERATOR
  //____________________________________________________________________________


  void read_options(vector<ST::string> & op,vector<ST::string> & vn);


  // Pointer auf das shrinkage-Parameter Fullcond-Objekt
  FC * get_shrinkagepointer();

  void get_samples(const ST::string & filename,ofstream & outg) const;
  //  void get_samples(const ST::string & filename, const unsigned & step = 1) const;

  //____________________________________________________________________________
  //
  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //____________________________________________________________________________

  void update(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);
  

  //____________________________________________________________________________
  //
  // FUNCTION: outresults
  // TASK: - write results to output window and files
  //____________________________________________________________________________

  void outresults(ofstream & out_stata, ofstream & out_R,
                  const ST::string & pathresults);

  //____________________________________________________________________________
  //
  // FUNCTION: outoptions
  // TASK: - write options to output window
  //____________________________________________________________________________

  void outoptions(void);


  //____________________________________________________________________________
  //
  // FUNCTION: reset
  // TASK: resets all parameters
  //____________________________________________________________________________

  void reset(void)
    {
    /*
    FC::reset();
    setbeta(beta.rows(),1,0.1);
    */
    }


  //____________________________________________________________________________
  //
  // DESTRUCTOR
  //____________________________________________________________________________

  ~FC_variance_pen_vector() {}

  }; // end: class FC_variance_pen_vector


} // end: namespace MCMC

#endif

