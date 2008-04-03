
#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include<StatwinFrame.h>

#endif

#include"mcmcsim.h"
#include<time.h>
#include"clstring.h"
#include <stdlib.h>
#include<math.h>


namespace MCMC
{

//------------------------------------------------------------------------------
//---------------------------- class equation  ---------------------------------
//------------------------------------------------------------------------------


equation::equation(const ST::string & h, DISTR * dp, const vector<FC*> fcp,
                   const ST::string & pd, const vector<ST::string> & paths)
  {
  header = h;
  distrp = dp;
  FCpointer = fcp;
  nrfc = FCpointer.size();
  FCpaths = paths;
  pathd = pd;
  }


equation::equation(const equation & s)
  {
  header = s.header;
  distrp = s.distrp;
  FCpointer = s.FCpointer;
  nrfc = s.nrfc;
  FCpaths = s.FCpaths;
  pathd = s.pathd;
  }


const equation & equation::operator=(const equation & s)
  {
  if (this == &s)
    return *this;
  header = s.header;
  distrp = s.distrp;
  FCpointer = s.FCpointer;
  nrfc = s.nrfc;
  FCpaths = s.FCpaths;
  pathd = s.pathd;  
  return *this;
  }


//------------------------------------------------------------------------------
//------------------------- end: class equation  -------------------------------
//------------------------------------------------------------------------------


MCMCsim::MCMCsim(GENERAL_OPTIONS * go,vector<equation> & equ)
  {
  genoptions = go;
  equations = equ;
  maxiterations = 1000;
  }


MCMCsim::MCMCsim(const MCMCsim & s)
  {
  genoptions = s.genoptions;
  equations = s.equations;
  maxiterations = s.maxiterations;
  }

  // OVERLOADED ASSIGNMENT CONSTRUCTOR

const MCMCsim & MCMCsim::operator=(const MCMCsim & s)
  {

  if (this == &s)
    return *this;
  genoptions = s.genoptions;
  equations = s.equations;
  maxiterations = s.maxiterations;
  return *this;
  }

  // FUNCTION: simulate
  // TASK: runs a MCMC simulation
  //       returns true, if simulation error or user break occured

bool MCMCsim::simulate(const int & seed, const bool & computemode)
  {
  unsigned i,j;

  unsigned nrmodels = equations.size();

  bool errors=false;

  if (errors==false)
  {

  unsigned it;
  unsigned iterations = genoptions->iterations;


  for (i=0;i<nrmodels;i++)
    {

    genoptions->out("\n");
    genoptions->out("\n");
    genoptions->out(equations[nrmodels-1-i].header +
    "\n",true,false,16);
    genoptions->out("\n");
    if (i==0)
      {
      genoptions->outoptions();
      genoptions->out("\n");
      }

    equations[nrmodels-1-i].distrp->outoptions();

    genoptions->out("OPTIONS FOR ESTIMATION:\n",true);
    genoptions->out("\n");
    for(j=0;j<equations[nrmodels-1-i].FCpointer.size();j++)
      equations[nrmodels-1-i].FCpointer[j]->outoptions();
    genoptions->out("\n");
    genoptions->out("MCMC SIMULATION STARTED\n",true);
    genoptions->out("\n");

    }

  //--------------- Compute posterior mode as starting value -------------------


  if (computemode)
    {

    genoptions->out("Computing starting values (may take some time)");

    bool c = posteriormode(true);
    }

  //-------------- end: Compute posterior mode as starting value ---------------


  #if defined(MICROSOFT_VISUAL)
    {

    }
  #elif!defined(__BUILDING_GNU)
    {
    randomize();
    }
  #else
    {
    srand(1);
    }
  #endif

  if(seed >= 0)
    srand(seed);

  clock_t beginsim = clock();
  clock_t it1per;
  clock_t endsim;
  bool runtime=false;

  #if defined(__BUILDING_GNU)
    double clk = (double)CLOCKS_PER_SEC;
  #else
    double clk = (double)CLK_TCK;
  #endif

  for (it=1;it<=iterations;it++)
    {

    if ( (runtime ==false) && (iterations/it == 100) )
      {
      runtime = true;
      it1per = clock();
      double sec = (it1per-beginsim)/clk;
      long timeleft = long(double(iterations-it)*(double(sec)/double(it)));
      long min = timeleft/60;
      sec = timeleft-min*60;
      if (min == 0)
        {
        genoptions->out("\n");
        genoptions->out
         ("  APPROXIMATE RUN TIME: " + ST::inttostring(long(sec)) +
                        " seconds\n");
        genoptions->out("\n");
        }
      else if (min == 1)
        {
        genoptions->out("\n");
        genoptions->out
                 ("  APPROXIMATE RUN TIME: " + ST::inttostring(min) +
                        " minute " + ST::inttostring(long(sec)) + " seconds\n");
        genoptions->out("\n");
        }
      else
        {
        genoptions->out("\n");
        genoptions->out
                     ("  APPROXIMATE RUN TIME: " + ST::inttostring(min)
                          + " minutes "
                          + ST::inttostring(long(sec)) + " seconds\n");
        genoptions->out("\n");
        }
      }  // end: if ( (runtime ==false) && (iterations/it == 100) )

    genoptions->update();

    for(i=0;i<nrmodels;i++)
      {

      equations[nrmodels-1-i].distrp->update();

      for(j=0;j<equations[nrmodels-1-i].FCpointer.size();j++)
        equations[nrmodels-1-i].FCpointer[j]->update();

      }



#if defined(BORLAND_OUTPUT_WINDOW)
    Application->ProcessMessages();

    if (Frame->stop)
      {
//      genoptions->out("USER BREAK\n");
      break;
      }

    if (Frame->pause)
      {
      genoptions->out("\n");
      genoptions->out("SIMULATION PAUSED\n");
      genoptions->out("Click CONTINUE to proceed\n");
      genoptions->out("\n");

      while (Frame->pause)
        {
        Application->ProcessMessages();
        }

      genoptions->out("SIMULATION CONTINUED\n");
      genoptions->out("\n");
      }
#elif defined(JAVA_OUTPUT_WINDOW)
      bool stop = genoptions->adminb_p->breakcommand();
      if(stop)
        break;
#endif

    } // end: for (i=1;i<=genoptions->iterations;i++)


#if defined(BORLAND_OUTPUT_WINDOW)
    if (!Frame->stop)
#elif defined(JAVA_OUTPUT_WINDOW)
    if (!genoptions->adminb_p->get_stop())
#endif
      {
      genoptions->out("\n");
      genoptions->out("SIMULATION TERMINATED\n",true);
      genoptions->out("\n");
      endsim = clock();
      long sec = (endsim-beginsim)/clk;
      long min = sec/60;
      sec = sec-min*60;
      if (min == 0)
        {
        genoptions->out("SIMULATION RUN TIME: "
                                + ST::inttostring(sec) +
                                " seconds\n");
        genoptions->out("\n");
        }
      else if (min == 1)
        {
        genoptions->out("SIMULATION RUN TIME: "
                                + ST::inttostring(min) + " minute "
                                +  ST::inttostring(sec) + " seconds\n");
        genoptions->out("\n");
        }
      else
        {
        genoptions->out("SIMULATION RUN TIME: " +
                                ST::inttostring(min) +  " minutes "
        + ST::inttostring(sec) + " seconds\n");
        genoptions->out("\n");
        }

      genoptions->out("\n");
      genoptions->out("ESTIMATION RESULTS:\n",true);
      genoptions->out("\n");

      for (i=0;i<nrmodels;i++)
        {

        if (nrmodels > 1)
          {
          genoptions->out("\n");
          genoptions->out(equations[nrmodels-1-i].header +
          "\n",true,false,16);
          genoptions->out("\n");
          }

        equations[nrmodels-1-i].distrp->outresults();

        for(j=0;j<equations[nrmodels-1-i].nrfc;j++)
          equations[nrmodels-1-i].FCpointer[j]->outresults(equations[nrmodels-1-i].FCpaths[j]);

        }

      return false;

      } // end: if Frame->stop

#if defined(BORLAND_OUTPUT_WINDOW)
    else
      {

      genoptions->out("\n");
      genoptions->out("SIMULATION TERMINATED BY USER BREAK\n");
      genoptions->out("\n");
      genoptions->out("Estimation results: none\n");
      genoptions->out("\n");

      for(i=0;i<nrmodels;i++)
        {

        equations[nrmodels-1-i].distrp->reset();

        for(j=0;j<equations[nrmodels-1-i].nrfc;j++)
          {
          equations[nrmodels-1-i].FCpointer[j]->reset())
          } // end: for(j=0;j<equations[nrmodels-1-i].nrfc;j++)

        }

      return true;
      }
#elif defined(JAVA_OUTPUT_WINDOW)
    else
      {

      genoptions->out("\n");
      genoptions->out("Estimation results: none\n");
      genoptions->out("\n");

      for(i=0;i<nrmodels;i++)
        {

        equations[nrmodels-1-i].distrp->reset();

        for(j=0;j<equations[nrmodels-1-i].nrfc;j++)
          {
          equations[nrmodels-1-i].FCpointer[j]->reset())
          } // end: for(j=0;j<equations[nrmodels-1-i].nrfc;j++)

        }

      return true;
      }
#endif

  } // end: no errors

  return true;

  }


bool MCMCsim::posteriormode(const bool & presim)
  {

  unsigned i,j;

  unsigned nrmodels = equations.size();

  bool errors=false;


  for (i=0;i<nrmodels;i++)
    {

    if (equations[nrmodels-1-i].header !="")
      {
      genoptions->out("\n");
      genoptions->out("\n");
      genoptions->out(equations[nrmodels-1-i].header + "\n",true,false,16);

      genoptions->out("\n");
      }


    if (!presim)
      {

      genoptions->out("RESPONSE DISTRIBUTION:\n",true);
      genoptions->out("\n");
      genoptions->out("  " + equations[nrmodels-1-i].distrp->family + "\n");
      genoptions->out("  Number of observations: " +
      ST::inttostring(equations[nrmodels-1-i].distrp->response.rows()) + "\n");
      genoptions->out("\n");

      } // end: if (!presim)

    }


    bool converged;
    bool allconverged;

    converged=false;

    unsigned it=1;

    while ((!converged) && (it <= maxiterations))
      {

      allconverged = true;

      for (i=0;i<nrmodels;i++)
        {

//        likep_mult[nrmodels-1-i]->compute_iwls();

        if (equations[nrmodels-1-i].distrp->posteriormode() == false)
          allconverged = false;

        for(j=0;j<equations[nrmodels-1-i].nrfc;j++)
          {
          if (equations[nrmodels-1-i].FCpointer[j]->posteriormode() == false)
              allconverged = false;
          } // end: for(j=0;j<equations[nrmodels-1-i].nrfc;j++)
        }


      if (allconverged)
        converged = true;

      #if defined(BORLAND_OUTPUT_WINDOW)

      Application->ProcessMessages();

      if (Frame->stop)
        {
        break;
        }

      if (Frame->pause)
        {
        genoptions->out("\n");
        genoptions->out("SIMULATION PAUSED\n");
        genoptions->out("Click CONTINUE to proceed\n");
        genoptions->out("\n");

      while (Frame->pause)
        {
        Application->ProcessMessages();
        }

      genoptions->out("SIMULATION CONTINUED\n");
      genoptions->out("\n");
      }
      #elif defined(JAVA_OUTPUT_WINDOW)
      bool stop = genoptions->adminb_p->breakcommand();
      if(stop)
        break;
      #endif

      it++;

      } // end: while ((!converged) && (it <= maxiterations))


    if (!presim)
      {
      #if defined(BORLAND_OUTPUT_WINDOW)
      if (!Frame->stop)
      #elif defined(JAVA_OUTPUT_WINDOW)
      if (!genoptions->adminb_p->get_stop())
      #endif
        {
        genoptions->out("\n");
        genoptions->out("ESTIMATION RESULTS:\n",true);
        genoptions->out("\n");

        genoptions->out("Number of Iterations: " + ST::inttostring(it) + "\n");
        if (!converged)
          genoptions->out("ALGORITHM DID NOT CONVERGE\n",true,true,12,255,0,0);
        genoptions->out("\n");


        for(i=0;i<nrmodels;i++)
          {

          equations[nrmodels-1-i].distrp->outresults(equations[nrmodels-1-i].pathd);

          for(j=0;j<equations[nrmodels-1-i].nrfc;j++)
            equations[nrmodels-1-i].FCpointer[j]->outresults(equations[nrmodels-1-i].FCpaths[j]);

          }

        } // end: if Frame->stop
      #if defined(BORLAND_OUTPUT_WINDOW)
      else
        {

        genoptions->out("\n");
        genoptions->out(
        "ESTIMATION TERMINATED BY USER BREAK\n");
        genoptions->out("\n");
        genoptions->out("Estimation results: none\n");
        genoptions->out("\n");


        for(i=0;i<nrmodels;i++)
          {

          equations[nrmodels-1-i].distrp->reset();

          for(j=0;j<equations[nrmodels-1-i].nrfc;j++)
            {
            equations[nrmodels-1-i].FCpointer[j]->reset())
            } // end: for(j=0;j<equations[nrmodels-1-i].nrfc;j++)

          }

        }
      #elif defined(JAVA_OUTPUT_WINDOW)
      else
        {

        genoptions->out("\n");
        genoptions->out("Estimation results: none\n");
        genoptions->out("\n");

        }
      #endif

      } // end: if (!presim)


  return converged;


  }


void MCMCsim::out_effects(const vector<ST::string> & paths)
  {


  }

  // FUNCTION: autocorr
  // TASK: computes autocorrelations for all samples parameters
  //      (i.e. for all beta matrices) and stores the result in datamatrix
  //      'cmat'

void MCMCsim::autocorr(const unsigned & lag,datamatrix & cmat)
  {

  }

  // FUNCTION: autocorr
  // TASK: computes autocorrelations for all samples parameters
  //      (i.e. for all beta matrices) and stores the result in file 'path'

void MCMCsim::autocorr(const unsigned & lag,const ST::string & path)
  {

  }

  // FUNCTION: get_samples
  // TASK: stores sampled parameters of all full conditionals in ASCII format
  //       for each full conditional one file will be created with filename
  //       'path' + title of the full conditional + "_sample.raw"

void MCMCsim::get_samples(
  #if defined(JAVA_OUTPUT_WINDOW)
  vector<ST::string> & newc,
  #endif
  const ST::string & path,const unsigned & step)
  {

  }


//------------------------------------------------------------------------------
//--------------------------- end: class MCMCsim -------------------------------
//------------------------------------------------------------------------------



} // end: namespace MCMC

#if defined(BORLAND_OUTPUT_WINDOW)
//---------------------------------------------------------------------------
#pragma package(smart_init)
#endif






