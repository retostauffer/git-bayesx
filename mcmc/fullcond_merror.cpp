
#include "first.h"

#include "fullcond_merror.h"

//------------------------------------------------------------------------------
//------------ CLASS: FULLCOND_merror implementation of member functions -------
//------------------------------------------------------------------------------


namespace MCMC
{

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // CONSTRUCTOR
  // o    : pointer to MCMCoptions object
  // t    : title of the full conditional (for example "fixed effects")
  //        (i.e. number of categories of the response variable)
  // fp   : file path for storing sampled parameters

  fullcond_merror::fullcond_merror(MCMCoptions * o, FULLCOND_nonp_gaussian * p,
           const datamatrix & d, const ST::string & t, const ST::string & fp)
                         : FULLCOND(o,d,t,d.rows(),d.cols(),fp)
    {
    designp = p;

    unsigned i, j;
    datamatrix meandata = datamatrix(d.rows(),1,0);
    unsigned dcols = d.cols();
    for(i=0; i<d.rows(); i++)
      {
      for(j=0; j<dcols; j++)
        {
        meandata(i,0) = d(i,j);
        }
      meandata(i,0) /= dcols;
      }
    setbeta(meandata);
    }

  // COPY CONSTRUCTOR

  fullcond_merror::fullcond_merror(const fullcond_merror & m)
    : FULLCOND(FULLCOND(m))
    {
    designp = m.designp;
    }

  // OVERLOADED ASSIGNMENT OPERATOR

  const fullcond_merror & fullcond_merror::operator=(const fullcond_merror & m)
    {
    if (this == &m)
      return *this;
    FULLCOND::operator=(FULLCOND(m));

    designp = m.designp;

    return *this;
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void fullcond_merror::update(void)
    {


    designp->init_data_varcoeff(data);
    }


  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool fullcond_merror::posteriormode(void)
    {
    }

  bool fullcond_merror::posteriormode_converged(const unsigned & itnr)
    {
    }

  void fullcond_merror::posteriormode_set_beta_mode(void)
    {
    }

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void fullcond_merror::outoptions(void)
    {
    }

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void fullcond_merror::outresults(void)
    {
    }

  // FUNCTION: reset
  // TASK: resets all parameters

  void fullcond_merror::reset(void)
    {
    }

  vector<ST::string> & fullcond_merror::get_results_latex(void)
    {
    }

} // end: namespace MCMC

