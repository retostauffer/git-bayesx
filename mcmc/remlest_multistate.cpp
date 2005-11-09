#include <remlest_multistate.h>

#if defined(BORLAND_OUTPUT_WINDOW)
#include "StatResults.h"
#include "statwinframe.h"
#endif

//------------------------------------------------------------------------------
//----------------------------- Constructor ------------------------------------
//------------------------------------------------------------------------------

remlest_multistate::remlest_multistate(
#if defined(JAVA_OUTPUT_WINDOW)
administrator_basic * adb,
#endif
vector<MCMC::FULLCOND*> & fc,datamatrix & re,
                const ST::string & family, const ST::string & ofile,
                const int & maxiter, const double & lowerlimit,
                const double & epsi, const double & maxch,
                const datamatrix & categories,
                const datamatrix & weight, ostream * lo)
  {

    #if defined(JAVA_OUTPUT_WINDOW)
    adminb_p = adb;
    #endif

    logout = lo;
    respfamily=family;
    outfile=ofile;

    maxit=maxiter;
    lowerlim=lowerlimit;
    eps=epsi;
    maxchange=maxch;

    fullcond = fc;
    unsigned i;

    }

//------------------------------------------------------------------------------
//----------------------------- REML estimation --------------------------------
//------------------------------------------------------------------------------

  // Function: estimate
  // Task: Perform REML-estimation for multi state models
  //       returns true if an error or user break occured

  bool remlest_multistate::estimate(const datamatrix resp,
               const datamatrix & offset, const datamatrix & weight)
    {
    return false;
    }

  // Function: estimate_glm
  // Task: Perform REML-estimation of multi state models if only fixed effects
  //       are specified

  bool remlest_multistate::estimate_glm(const datamatrix resp,
               const datamatrix & offset, const datamatrix & weight)
     {
     return false;
     }

//------------------------------------------------------------------------------
//----------------------------- Object description -----------------------------
//------------------------------------------------------------------------------

  void remlest_multistate::outoptions()
    {
    }

//------------------------------------------------------------------------------
//----------------------------- Writing results --------------------------------
//------------------------------------------------------------------------------

void remlest_multistate::make_plots(ofstream & outtex,ST::string path_batch,
                         ST::string path_splus)
  {
  }

void remlest_multistate::make_model(ofstream & outtex, const ST::string & rname)
  {
  }

void remlest_multistate::make_predictor(ofstream & outtex)
  {
  }

void remlest_multistate::make_prior(ofstream & outtex)
  {
  }

void remlest_multistate::make_options(ofstream & outtex)
  {
  }

void remlest_multistate::make_fixed_table(ofstream & outtex)
  {
  }

void remlest_multistate::make_graphics(const ST::string & title,
                     const ST::string & path_batch,
                     const ST::string & path_tex,
                     const ST::string & path_splus,
                     const ST::string & rname,
                     const bool & dispers)
  {
  }

bool remlest_multistate::check_pause()
  {
#if defined(BORLAND_OUTPUT_WINDOW)
  Application->ProcessMessages();
  if (Frame->stop)
    {
    return true;
    }

  if (Frame->pause)
    {
    out("\n");
    out("ESTIMATION PAUSED\n");
    out("Click CONTINUE to proceed\n");
    out("\n");

    while (Frame->pause)
      {
      Application->ProcessMessages();
      }

    out("ESTIMATION CONTINUED\n");
    out("\n");
    }
  return false;
#elif defined(JAVA_OUTPUT_WINDOW)
  return adminb_p->breakcommand();
#endif
  }

void remlest_multistate::out(const ST::string & s,bool thick,bool italic,
                      unsigned size,int r,int g, int b)
  {
#if defined(BORLAND_OUTPUT_WINDOW)
  ST::string sh = s;
  sh = sh.replaceallsigns('\n',' ');
  if (!Frame->suppoutput)
    Results->ResultsRichEdit->Lines->Append(sh.strtochar());
 if (!(logout->fail()))
    (*logout) << s << flush;
#elif defined(JAVA_OUTPUT_WINDOW)
  ST::string sh = s;
  sh = sh.replaceallsigns('\n',' ');
  sh = sh+"\n";
  if (!adminb_p->get_suppressoutput())
    adminb_p->Java->CallVoidMethod(adminb_p->BayesX_obj, adminb_p->javaoutput,
    adminb_p->Java->NewStringUTF(sh.strtochar()),thick,italic,size,r,g,b);
  if (!(logout->fail()))
    (*logout) << s << flush;
#else
  (*logout) << s << flush;
#endif
  }


void remlest_multistate::outerror(const ST::string & s)
  {
  out(s,true,true,12,255,0,0);
  }







