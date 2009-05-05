
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FCINCLUDED)

#define FCINCLUDED

#include"statmat.h"
#include"sparsemat.h"

#include"Random.h"
#include"../values.h"
#include<fstream>
#include<vector>
#include<bitset>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include<cmath>

namespace MCMC
{

using std::vector;
using std::bitset;
using std::ofstream;

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC ----------------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC
  {

  protected:

  ST::string samplepath;         // filename for storing sampled parameters
  ofstream samplestream;         // stream object for storing sampled parameters

  vector<ST::string> priorassumptions;


  public:

  bool samplestore;

  GENERAL_OPTIONS * optionsp;    // Pointer to general MCMC options


  bool nosamples;

  ST::string title;              // Title/name of the full conditional

  datamatrix beta;               // Matrix of current parameters
  datamatrix beta_mode;          // for posteriormode, the parameters
                                 // of the last iteration of iwls
                                 // used to check if convergence has been
                                 // already achieved
  datamatrix betamean;           // Sampling mean of parameters
  datamatrix betas2;             // Sampling sum of squares of parameters
  datamatrix betavar;            // Sampling variance of parameters
  datamatrix betamin;            // Sampling minimum
  datamatrix betamax;            // Sampling maximum


  datamatrix betaqu_l1_lower;    // (100-level1)/2 percent quantile
  datamatrix betaqu_l2_lower;    // (100-level2)/2 percent quantile
  datamatrix betaqu50;           //  50 percent quantile
  datamatrix betaqu_l1_upper;    // corresponding upper quantile (level1)
  datamatrix betaqu_l2_upper;    // corresponding upper quantile (level2)

  datamatrix betameanold;        // Will be initialized after the burnin period
  datamatrix betavarold;
  datamatrix betaminold;
  datamatrix betamaxold;

  datamatrix sampled_beta;       // sampled beta's, the i-th row contains the
                                 // i-th sample, i.e. the number of rows
                                 // corresponds to the number of stored samples
                                 // the number of columns correspond to the
                                 // number of parameters

  datamatrix transform;          // The factor with which all beta's will be
                                 // multiplied before storing them and computing
                                 // means, std, etc.
                                 // DEFAULT: transform = 1
  double addon;                  // An additive constant that will be added
                                 // on each component of beta before storing
                                 // DEFAULT: addon = 0;



  unsigned long acceptance;      // number of accepted iterations
  unsigned long nrtrials;        // number of trials


  unsigned column;               // the response category the fc belongs to


  double meaneffect;            // for results in original scale          

  // FUNCTION: readsample
  // TASK: reads sample of parameter 'nr' and stores the sample in datamatrix
  //       sample (sample must be optionsp->samplesize x 1 matrix)
  //       nr >= 0 and nr < nrpar required

  void readsample(datamatrix & sample,const unsigned & nr,
                  const unsigned & col=0) const;


  // FUNCTION readsample
  // TASK: reads the nr-th stored parameter matrix of the sample
  //       'b' must have proper size, i.e. b.rows() == beta.rows() and
  //       b.cols() == b.cols()

  void readsample2(datamatrix & b,const unsigned & nr) const;


  // FUNCTION: readsample
  // TASK: reads ALL stored paramters of the sample
  //       'b' must have proper size, i.e. b.rows() = samplesize
  //       b.cols() = nrpar

  void readsample3(datamatrix & b) const;

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // rows : number of rows of the beta matrix (i.e. number of parameters)
  // cols : number of columns of the beta matrix
  //        (i.e. number of categories of the response variable)
  // fp   : file path for storing sampled parameters

  FC(GENERAL_OPTIONS * o,const ST::string & t,
           const unsigned & rows, const unsigned & cols,
           const ST::string & fp,bool sstore);

  // COPY CONSTRUCTOR

  FC(const FC & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC & operator=(const FC & m);

  // DESTRUCTOR

  ~FC()
    {
    if (nosamples==false)
      {
      remove(samplepath.strtochar());
      }
    }


  //------------- ACCESS TO PARAMETERS AND OTHER CHARACTERISTICS ---------------

   // FUNCTION: setbeta
   // TASK: initializes beta matrices (i.e. beta, betaold, etc)

   void setbeta(const unsigned & rows,const unsigned & cols,const double & v);

   void setbeta(const datamatrix & betanew);

   void setbetavalue(const unsigned & row,const unsigned & col,const double & v);

  // FUNCTION: compute_autocorr
  // TASK: computes autocorrelation function for lags 1 - 'lag' for parameter
  //       beta(row,col). returns the result as  a column vector of
  //       autocorrelations

  datamatrix compute_autocorr(const unsigned & lag,const unsigned & row,
                              const unsigned & col) const;

  double compute_autocorr_single(const unsigned & lag,const unsigned & row,
                                      const unsigned & col) const;


  // FUNCTION: compute_autocorr
  // TASK: computes autocorrelation function for lags 1 - 'lag' for parameter
  //       beta. writes the autocorrelations in a file 'path'

  void compute_autocorr(const ST::string & path, unsigned lag) const;

  virtual void compute_autocorr_all(const ST::string & path, unsigned lag,
                                    ofstream & outg) const;

  // FUNCTION: get_samples
  // TASK: stores the sampled parameters in ASCII format

  virtual void get_samples(const ST::string & filename,ofstream & outg) const;

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  virtual void update(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  virtual bool posteriormode(void);

  // FUNCTION: psoteriormode_betamean
  // TASK: computes the retransformed current betamean

  void posteriormode_betamean(void);

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  virtual void outoptions(void)
    {
    }

  // FUNCTION: outgraphs
  // TASK: writes batch files for STATA and R for visualizing results

  virtual void outgraphs(ofstream & out_stata, ofstream & out_R,
                         const ST::string & path)
    {
    }


  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  virtual void outresults(ofstream & out_stata, ofstream & out_R,
               const ST::string & pathresults);

  void outresults_acceptance(void);

  // FUNCTION: reset
  // TASK: resets all parameters

  virtual void reset(void);

  virtual void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  };


} // end: namespace MCMC

#endif


