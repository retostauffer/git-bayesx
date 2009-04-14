
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

equation::equation(void)
  {
  nrfc=0;
  header="";
  paths ="";
  FCpointer.reserve(30);
  FCpaths.reserve(30);
  equationnr = 1;
  equationtype = "mean";
  hlevel = 1;
  }


equation::equation(int enr, int hl,ST::string t)
  {
  nrfc=0;
  header="";
  paths ="";
  FCpointer.reserve(30);
  FCpaths.reserve(30);
  equationnr = enr;
  equationtype = t;
  hlevel = hl;
  }


equation::equation(const ST::string & h, DISTR * dp, const vector<FC*> fcp,
                   const ST::string & pd, const vector<ST::string> & ps)
  {
  equationnr = 1;
  equationtype = "mean";
  hlevel = 1;
  header = h;
  paths="";
  distrp = dp;
  FCpointer = fcp;
  nrfc = FCpointer.size();
  FCpaths = ps;
  pathd = pd;

  }


equation::equation(const equation & s)
  {
  equationnr = s.equationnr;
  equationtype = s.equationtype;
  hlevel = s.hlevel;
  header = s.header;
  paths = s.paths;
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
  paths = s.paths;
  distrp = s.distrp;
  FCpointer = s.FCpointer;
  nrfc = s.nrfc;
  FCpaths = s.FCpaths;
  pathd = s.pathd;
  return *this;
  }


void equation::add_FC(FC * FCp,const ST::string & p)
  {
  FCpointer.push_back(FCp);
  FCpaths.push_back(p);
  nrfc++;
  }

//------------------------------------------------------------------------------
//------------------------- end: class equation  -------------------------------
//------------------------------------------------------------------------------


MCMCsim::MCMCsim(GENERAL_OPTIONS * go,vector<equation> & equ)
  {
  genoptions = go;
  equations = equ;
  maxiterations = 1000;
//  maxiterations = 10;
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

bool MCMCsim::simulate(ST::string & pathgraphs, const int & seed, const bool & computemode)
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
    genoptions->out(equations[nrmodels-1-i].header +
    "\n",true,false,16);
    genoptions->out("\n");
    if (i==0)
      {
      genoptions->outoptions();
      genoptions->out("\n");
      }

    equations[nrmodels-1-i].distrp->outoptions();

    if (equations[nrmodels-1-i].FCpointer.size() > 0)
      {
      genoptions->out("OPTIONS FOR ESTIMATION:\n",true);
      genoptions->out("\n");

      for(j=0;j<equations[nrmodels-1-i].FCpointer.size();j++)
        {
        equations[nrmodels-1-i].FCpointer[j]->outoptions();
        }
      genoptions->out("\n");
      }
    }

  genoptions->out("MCMC SIMULATION STARTED\n",true);
  genoptions->out("\n");

  //--------------- Compute posterior mode as starting value -------------------

  if (computemode)
    {

    genoptions->out("  COMPUTING STARTING VALUES (MAY TAKE SOME TIME)");
    genoptions->out("\n");
    ST::string h = "";
    bool c = posteriormode(h,true);
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
      genoptions->out("USER BREAK\n");
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

      ST::string pathstata = pathgraphs + "_stata.do";
      ST::string pathR = pathgraphs + "_R.r";

      ofstream out_stata(pathstata.strtochar());
      ofstream out_R(pathR.strtochar());

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
          equations[nrmodels-1-i].FCpointer[j]->outresults(out_stata,out_R,equations[nrmodels-1-i].FCpaths[j]);

        }

      genoptions->out("  FILES FOR VISULAZING RESULTS:\n",true,true,12,255,0,0);
      genoptions->out("\n");
      genoptions->out("    STATA DO-FILE\n");
      genoptions->out("\n");
      genoptions->out("    " + pathstata);
      genoptions->out("\n");


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

/*
      for(i=0;i<nrmodels;i++)
        {

        equations[nrmodels-1-i].distrp->reset();

        for(j=0;j<equations[nrmodels-1-i].nrfc;j++)
          {
          equations[nrmodels-1-i].FCpointer[j]->reset();
          } // end: for(j=0;j<equations[nrmodels-1-i].nrfc;j++)

        }
*/
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


bool MCMCsim::posteriormode(ST::string & pathgraphs, const bool & presim)
  {

  unsigned i,j;

  unsigned nrmodels = equations.size();

  bool errors=false;


  for (i=0;i<nrmodels;i++)
    {

    if (!presim)
      {

      if (equations[nrmodels-1-i].header !="")
        {
        genoptions->out("\n");
        genoptions->out("\n");
        genoptions->out(equations[nrmodels-1-i].header + "\n",true,false,16);

        genoptions->out("\n");
        }

      genoptions->out("RESPONSE DISTRIBUTION:\n",true);
      genoptions->out("\n");
      genoptions->out("  " + equations[nrmodels-1-i].distrp->family + "\n");
      genoptions->out("  Number of observations: " +
      ST::inttostring(equations[nrmodels-1-i].distrp->response.rows()) + "\n");
      genoptions->out("  Number of observations with positive Weights: " +
      ST::inttostring(equations[nrmodels-1-i].distrp->response.rows()-
      equations[nrmodels-1-i].distrp->nrzeroweights
      ) + "\n");

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
      else
        it++;

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


        ST::string pathstata = pathgraphs + "_stata.do";
        ST::string pathR = pathgraphs + "_R.r";

        ofstream out_stata(pathstata.strtochar());
        ofstream out_R(pathR.strtochar());


        for(i=0;i<nrmodels;i++)
          {

          equations[nrmodels-1-i].distrp->outresults(equations[nrmodels-1-i].pathd);

          for(j=0;j<equations[nrmodels-1-i].nrfc;j++)
            equations[nrmodels-1-i].FCpointer[j]->outresults(out_stata,out_R,equations[nrmodels-1-i].FCpaths[j]);

          }

        genoptions->out("  FILES FOR VISULAZING RESULTS:\n",true,true,12,255,0,0);
        genoptions->out("  STATA DO-FILE\n");
        genoptions->out("\n");
        genoptions->out(pathstata);
        genoptions->out("\n");



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
            equations[nrmodels-1-i].FCpointer[j]->reset();
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


unsigned MCMCsim::compute_nrpar(void)
  {
  unsigned j,l;
  unsigned nrpar=0;
  unsigned nrmodels = equations.size();

  for (l=0;l<nrmodels;l++)
    {
    for(j=0;j<equations[l].FCpointer.size();j++)
      {
      if (equations[l].FCpointer[j]->nosamples == false)
        nrpar += equations[l].FCpointer[j]->beta.rows()*
                 equations[l].FCpointer[j]->beta.cols();
      }
    }

  return nrpar;
  }


  // FUNCTION: autocorr
  // TASK: computes autocorrelations for all samples parameters
  //      (i.e. for all beta matrices) and stores the result in datamatrix
  //      'cmat'

void MCMCsim::autocorr(const unsigned & lag,datamatrix & cmat)
  {

  unsigned p = compute_nrpar();

  cmat = datamatrix(lag,p);

  unsigned j,i,k,l;
  unsigned col =0;

  unsigned nrmodels = equations.size();


  for (l=0;l<nrmodels;l++)
    {

    for(j=0;j<equations[l].FCpointer.size();j++)
      {
      if (equations[l].FCpointer[j]->nosamples==false)
        {
        for (k=0;k<equations[l].FCpointer[j]->beta.cols();k++)
          for(i=0;i<equations[l].FCpointer[j]->beta.rows();i++)
            {
            cmat.putCol(col,equations[l].FCpointer[j]->compute_autocorr(lag,i,k));
            col++;

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

            }

        #if defined(BORLAND_OUTPUT_WINDOW)
        Application->ProcessMessages();
        if (Frame->stop)
          {
          cmat = datamatrix(1,1);
          break;
          }
        #elif defined(JAVA_OUTPUT_WINDOW)

        if (genoptions->adminb_p->get_stop())
          {
          cmat = datamatrix(1,1);
          break;
          }
        #endif
        }
      } // end:  for(j=0;j<fullcondp.size();j++)

    } // for (l=0;l<nrmodels;l++)

  }

  // FUNCTION: autocorr
  // TASK: computes autocorrelations for all samples parameters
  //      (i.e. for all beta matrices) and stores the result in file 'path'

/*
void MCMCsim::autocorr(const unsigned & lag,const ST::string & path)
  {

  unsigned nrmodels = equations.size();

  ofstream out(path.strtochar());
  assert(!out.fail());
  assert(!out.bad());

  ST::string name;
  datamatrix cmat;
  unsigned i,j,k,l,nr,s;
  double min,mean,max,n;
  bool miss,misstot;
  misstot = false;
  genoptions->out("Computing autocorrelation functions...\n");
  autocorr(lag,cmat);

  #if defined(BORLAND_OUTPUT_WINDOW)
  if (!Frame->stop)
  #elif defined(JAVA_OUTPUT_WINDOW)
  if (!genoptions->adminb_p->get_stop())
  #endif
    {
    out << "lag ";

    for (s=0;s<nrmodels;s++)
      {

      for(j=0;j<equations[s].FCpointer.size();j++)
        {
        if (equations[s].FCpointer[j]->nosamples == false)
          {
          name = equations[s].FCpointer[j]->title;

          for (k=0;k<equations[s].FCpointer[j]->beta.cols();k++)
            for(i=0;i<equations[s].FCpointer[j]->beta.rows();i++)
              {
              if (equations[s].FCpointer[j]->beta.cols() == 1)
                out << name << "_" << (i+1) << " ";
              else
                out << name << (i+1) << "_" << (k+1) << " ";
              }

          out << name << "_min " << name << "_mean " << name << "_max ";
          }

        }  // end: for(j=0;j<fullcondp.size();j++)
      }

      out << endl;

    for(l=0;l<lag;l++)
      {
      nr = 0;
      out << (l+1) << " ";

      for (s=0;s<equations.size();s++)
        {

        for(j=0;j<equations[s].FCpointer.size();j++)
          {
          if (equations[s].FCpointer[j]->nosamples == false)
            {
            min = 1;
            max = -1;
            mean = 0;
            miss = true;

            for (k=0;k<equations[s].FCpointer[j]->beta.cols();k++)
              for(i=0;i<equations[s].FCpointer[j]->beta.rows();i++)
                {
                if (cmat(l,nr) <= 1)
                  {
                  miss = false;
                  if (cmat(l,nr) > max)
                    max = cmat(l,nr);
                  if (cmat(l,nr) < min)
                    min = cmat(l,nr);
                  mean+=cmat(l,nr);
                  out << cmat(l,nr) << " ";
                  }
                else
                  {
                  misstot = true;
                  out << "NA ";
                  }

                nr++;
                }

            if (miss)
              out << "NA NA NA ";
            else
              {
              n = equations[s].FCpointer[j]->beta.cols() * equations[s].FCpointer[j]->beta.rows();
              out << min << " " << (mean/n) << " " << max << " ";
              }
            } // end: if (equations[s].FCpointer[j]->nosamples == false)
          }  // end: for(j=0;j<equations[s].FCpointer.size();j++)
        } // for (s=0;s<equations.size();s++)
        out << endl;
      } // end: for(l=0;l<lag;l++)

      genoptions->out("Autocorrelation functions computed and stored in file\n");
      genoptions->out(path + ".\n");
      genoptions->out("\n");
      if (misstot)
        {
        genoptions->out("WARNING: There were undefined autocorrelations\n",true,true);
        genoptions->out("\n");
        }

    #if !defined(JAVA_OUTPUT_WINDOW)
      genoptions->out("They may be visualized using the R function 'plotautocor'.\n");
      genoptions->out("\n");
    #endif
      } // end: if (!Frame->stop)
    #if defined(BORLAND_OUTPUT_WINDOW)
    else
      {
      genoptions->out("USER BREAK\n");
      genoptions->out("No autocorrelation functions computed\n");
      genoptions->out("\n");
      out.close();
      remove(path.strtochar());
      }
    #elif defined(JAVA_OUTPUT_WINDOW)
    else
      {
//      genoptions->out("SIMULATION TERMINATED BY USER BREAK\n");
      genoptions->out("No autocorrelation functions computed\n");
      genoptions->out("\n");
      out.close();
      remove(path.strtochar());
      }
    #endif

  }   // end: autocorr
*/


void MCMCsim::autocorr(const unsigned & lag)
  {

  unsigned i,j,k,l,s,r,c,nrpar;
  double min,mean,max;
  unsigned nrmodels = equations.size();
  double autoc;
  ST::string path;
  bool misstot;

  genoptions->out("Computing autocorrelation functions...\n");
  genoptions->out("Autocorrelations are stored in file(s):\n");
  genoptions->out("\n");


  for (s=0;s<nrmodels;s++)
    {

    for(j=0;j<equations[s].FCpointer.size();j++)
      {
      if (equations[s].FCpointer[j]->nosamples == false)
        {
        path = equations[s].FCpaths[j].substr(0,equations[s].FCpaths[j].length()-4) + "_autocor.raw";
        ofstream out(path.strtochar());

        genoptions->out(path);
        genoptions->out("\n");

        out << "lag  ";

        for (k=0;k<equations[s].FCpointer[j]->beta.cols();k++)
          for(i=0;i<equations[s].FCpointer[j]->beta.rows();i++)
            {
            if (equations[s].FCpointer[j]->beta.cols() == 1)
              out << "b_" << (i+1) << " ";
            else
              out << "b_" << (i+1) << "_" << (k+1) << " ";
            }

        out  << "b_min " << "b_mean " << "b_max " << endl;

        misstot = false;
        for(l=1;l<=lag;l++)
          {

          nrpar = 0;

          out << l << "  ";

          min = 1;
          max = -1;
          mean = 0;

          for(c=0;c<equations[s].FCpointer[j]->beta.cols();c++)
            for (r=0;r<equations[s].FCpointer[j]->beta.rows();r++)
              {
              autoc = equations[s].FCpointer[j]->compute_autocorr_single(l,r,c);
              if ( (autoc <= 1) && (autoc >= -1) )
                {
                nrpar++;
                if (autoc < min)
                  min = autoc;
                if (autoc > max)
                  max = autoc;
                mean += autoc;
                out << autoc << "  ";
                }
              else
                {
                out << "NA  " << endl;
                misstot = true;
                }
              }

          out << min << "  ";
          out << max << "  ";
          out << mean/nrpar << "  ";
          out << endl;
          }  // end: for(l=0;l<lag;l++)

        if (misstot==true)
          {
          genoptions->out("WARNING: There were undefined autocorrelations\n",true,true);
          genoptions->out("\n");
          }

        }
      }  // end: for(j=0;j<fullcondp.size();j++)
    } // end: for (s=0;s<nrmodels;s++)

  }   // end: autocorr



  // FUNCTION: get_samples
  // TASK: stores sampled parameters of all full conditionals in ASCII format
  //       for each full conditional one file will be created with filename
  //       'path' + title of the full conditional + "_sample.raw"

void MCMCsim::get_samples(
  #if defined(JAVA_OUTPUT_WINDOW)
  vector<ST::string> & newc
  #endif
  )
  {

  unsigned i,j;
  ST::string filename;
  ST::string help;
  ST::string psname;

  genoptions->out("Storing sampled parameters...\n");
  genoptions->out("Sampled parameters are stored in file(s):\n");
  genoptions->out("\n");

  for(j=0;j<equations.size();j++)
    {
    for(i=0;i<equations[j].FCpointer.size();i++)
      {
      if (equations[j].FCpointer[i]->nosamples == false)
        {
        filename =  equations[j].FCpaths[i].substr(0,equations[j].FCpaths[i].length()-4) + "_sample.raw";
        equations[j].FCpointer[i]->get_samples(filename);
        genoptions->out(filename + "\n");
        #if defined(JAVA_OUTPUT_WINDOW)

        psname = equations[j].FCpaths[i].substr(0,equations[j].FCpaths[i].length()-4) + "_sample.ps";
        newc.push_back("dataset _dat");
        newc.push_back("_dat.infile , nonote using " + filename);
        newc.push_back("graph _g");
        newc.push_back("_g.plotsample , replace outfile=" +
                      psname  + " using _dat");
        genoptions->out(psname + " (graphs)\n");
        newc.push_back("drop _dat _g");

        #endif
        genoptions->out("\n");
        }
      }
    }

  /*
  if (likepexisting)
    {
    for(i=0;i<likep_mult.size();i++)
      {

      if (likep_mult[i]->get_scaleexisting())
        {
        genoptions->out("\n");
        filename = likep_mult[i]->get_scale_sample();
        genoptions->out(filename+"\n");
        genoptions->out("\n");
        #if defined(JAVA_OUTPUT_WINDOW)

        psname = filename.substr(0,filename.length()-4) +   + ".ps";
        newc.push_back("dataset _dat");
        newc.push_back("_dat.infile , nonote using " + filename);
        newc.push_back("graph _g");
        newc.push_back("_g.plotsample , replace outfile=" +
                      psname  + " using _dat");
        genoptions->out(psname + " (graphs)\n");
        newc.push_back("drop _dat _g");

        #endif

        }
      }
    }
   */

  genoptions->out("\n");
  genoptions->out("Storing completed\n");
  genoptions->out("\n");
  #if defined(BORLAND_OUTPUT_WINDOW)
  genoptions->out(
  "Sampled parameters may be visualized using the R\n");
  genoptions->out("function 'plotsample'.\n");
  genoptions->out("\n");
  #endif

  }


//------------------------------------------------------------------------------
//--------------------------- end: class MCMCsim -------------------------------
//------------------------------------------------------------------------------



} // end: namespace MCMC

#if defined(BORLAND_OUTPUT_WINDOW)
//---------------------------------------------------------------------------
#pragma package(smart_init)
#endif






