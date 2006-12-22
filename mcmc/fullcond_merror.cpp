
#include "first.h"

#include "fullcond_merror.h"

//------------------------------------------------------------------------------
//------------ CLASS: FULLCOND_merror implementation of member functions -------
//------------------------------------------------------------------------------


namespace MCMC
{

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // CONSTRUCTOR (Susi)
  // o    : pointer to MCMCoptions object
  // t    : title of the full conditional (for example "fixed effects")
  //        (i.e. number of categories of the response variable)
  // fp   : file path for storing sampled parameters

  fullcond_merror::fullcond_merror(MCMCoptions * o, FULLCOND_nonp_gaussian * p,
           DISTRIBUTION * dp, const datamatrix & d, const ST::string & t, const ST::string & fp)
                         : FULLCOND(o,d,t,d.rows(),d.cols(),fp)
    {
    designp = p;
    likep = dp;

    varcoeff=true;
    unsigned i, j;
    meandata = datamatrix(d.rows(),1,0);
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

// BEGIN: merror
  // CONSTRUCTOR (Thomas)
  fullcond_merror::fullcond_merror(MCMCoptions * o, spline_basis * p,
           DISTRIBUTION * dp, const datamatrix & d, const ST::string & t,
           const ST::string & fp, const ST::string & pres)
                         : FULLCOND(o,d,t,d.rows(),1,fp)
  {
  splinep = p;
  likep= dp;
  varcoeff=false;
  unsigned i, j;
  meandata = datamatrix(d.rows(),1,0);
  unsigned dcols = d.cols();
  for(i=0; i<d.rows(); i++)
    {
    for(j=0; j<dcols; j++)
      {
      meandata(i,0) += d(i,j);
      }
    meandata(i,0) /= dcols;
    }
  minx = meandata.min(0);
  maxx = meandata.max(0);
  setbeta(meandata);

  index = statmatrix<int>(beta.rows(),1,0);
  for(i=0; i<beta.rows(); i++)
    {
    index(i,0)=i;
    }
  }
// END: merror

  // COPY CONSTRUCTOR

  fullcond_merror::fullcond_merror(const fullcond_merror & m)
    : FULLCOND(FULLCOND(m))
    {
    designp = m.designp;
    likep = m.likep;
// BEGIN: merror
    splinep = m.splinep;
    minx = m.minx;
    maxx = m.maxx;
    meandata = m.meandata;
    old = m.old;
    index = m.index;
// END: merror
    }

  // OVERLOADED ASSIGNMENT OPERATOR

  const fullcond_merror & fullcond_merror::operator=(const fullcond_merror & m)
    {
    if (this == &m)
      return *this;
    FULLCOND::operator=(FULLCOND(m));

    designp = m.designp;
    likep = m.likep;
// BEGIN: merror
    splinep = m.splinep;
    minx = m.minx;
    maxx = m.maxx;
    meandata = m.meandata;
    old = m.old;
    index = m.index;
// END: merror

    return *this;
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void fullcond_merror::update(void)
    {
    if(varcoeff)
      {
      designp->init_data_varcoeff(beta);
      }
    else
      {
      // generate proposed values

      old = beta;              // old stores results from the previous iteration
      // Fake generation
      unsigned i;
      double u;
      double acchelp=0;
      double * work = beta.getV();
      double * workmean = meandata.getV();
      for(i=0;i<beta.rows();i++,work++,workmean++)
        {
        *work = *workmean + 0.02*rand_normal();
        while(*work<minx || *work>maxx)
          *work = *workmean + 0.02*rand_normal();
        }

      // extract current f(x) from spline_basis

      datamatrix currentspline = splinep->get_spline();

      // call update_merror and compute new values for f(x).
      splinep->update_merror(beta);
      datamatrix diffspline = splinep->get_spline()-currentspline;

      // Compute acceptance probabilities, update acceptance and overwrite non-accepted entries in beta with values in old
      datamatrix logold(beta.rows(),1,0);
      datamatrix propold(beta.rows(),1,0);
      datamatrix lognew(beta.rows(),1,0);
      datamatrix propnew(beta.rows(),1,0);
      // nichtinformative Priori für die wahren Werte
      for(i=0; i<beta.rows(); i++)
        {
        logold(i,0) = likep->loglikelihood(i,i,index);
        propold(i,0) = -0.5*(old(i,0)-meandata(i,0))*(old(i,0)-meandata(i,0)) / (0.02*0.02);
        }

      likep->addtocurrent(diffspline);
      for(i=0; i<beta.rows(); i++)
        {
        lognew(i,0) = likep->loglikelihood(i,i,index,false);
        propnew(i,0) = -0.5*(beta(i,0)-meandata(i,0))*(beta(i,0)-meandata(i,0)) / (0.02*0.02);
        }

      for(i=0; i<beta.rows(); i++)
        {
        u = log(uniform());
        double test1 = lognew(i,0) - logold(i,0);
        double test2 = propold(i,0) - propnew(i,0);
        if(u <= (test1+test2))
          {
          acchelp++;
          }
        else
          {
          beta(i,0) = old(i,0);
          }
        }
      acceptance += acchelp/beta.rows();

      // update spline_basis with the final information*/
      splinep->update_merror(beta);

      // Update the linear predictor
      diffspline = splinep->get_spline()-currentspline;

/*      ofstream out1("c:\\temp\\currentspline.raw");
      currentspline.prettyPrint(out1);
      out1.close();
      ofstream out2("c:\\temp\\propspline.raw");
      propspline.prettyPrint(out2);
      out2.close();*/

      likep->addtocurrent(diffspline);
      likep->swap_linearpred();

      FULLCOND::update();
      }
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
    FULLCOND::outresults();
//    optionsp->out("  No results for merror!\n");
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

