// DATE: 01.02.99

//---------------------------------------------------------------------------
#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include<StatwinFrame.h>
#endif

#include<mcmcsimul2.h>
#include<time.h>
#include<clstring.h>
#include <stdlib.h>
#include<math.h>


namespace MCMC
{


STEPWISErun::STEPWISErun(MCMCoptions * go,DISTRIBUTION * dp,
vector<FULLCOND*> & fc)
: MCMCsimulate(go,dp,fc)
  {
  }


STEPWISErun::STEPWISErun(const STEPWISErun & st)
  : MCMCsimulate(MCMCsimulate(st))
  {
  }


const STEPWISErun & STEPWISErun::operator=(
const STEPWISErun & st)
  {
  if (this==&st)
    return *this;
  MCMCsimulate::operator=(MCMCsimulate(st));
  return *this;
  }


double STEPWISErun::compute_criterion(void)
  {
  double df = 0;
  unsigned i;
  for(i=0;i<fullcondp.size();i++)
    {
    df = df + fullcondp[i]->compute_df();
    }
  double kriterium;
  double nrobs = likep_mult[0]->get_nrobs();
  //double nrobs = likep_mult[0]->get_nrobs_wpw();

  // Versuch: mit Devianz
  double deviance = 0;
  double deviancesat = 0;
  likep_mult[0]->compute_deviance(deviance,deviancesat);

  if(likep_mult[0]->get_family()=="Gaussian")
    deviance = deviance - nrobs*(1+log(2*M_PI));

  if (criterion=="GCV")
    {
    if(likep_mult[0]->get_family()=="Gaussian")
      kriterium = likep_mult[0]->compute_rss();
    else
      kriterium = deviancesat;
    kriterium = kriterium / (nrobs*(1-df/nrobs)*(1-df/nrobs));
    }
    // kriterium = likep_mult[0]->compute_gcv(df);
  else if(criterion=="AIC")
    kriterium = deviance + 2*df;
    // kriterium = likep_mult[0]->compute_aic(df);
  else if(criterion=="BIC")
    kriterium = deviance + log(nrobs)*df;
    // kriterium = likep_mult[0]->compute_bic(df);
  else  //if(criterion=="AIC_imp")
    kriterium = deviance + 2*df + 2*df*(df+1)/(nrobs-df-1);
    // kriterium = likep_mult[0]->compute_improvedaic(df);

  return kriterium;
  }

bool STEPWISErun::stepwise(const ST::string & crit, const int & stp,
                   const ST::string & trac, const int & number,
                   const ST::string & stam, const int & inc,
                   const bool & finet, const datamatrix & D,
                   const vector<ST::string> & modelv, const ST::string & name,
                   vector<FULLCOND*> & fullcond_z, ST::string & path)
  {

  criterion = crit;
  increment = inc;
  steps = stp;
  startmodel = stam;
  fine_tuning = finet;
  trace = trac;

  vector<double> modell_alt;
  double kriterium_alt;
  ST::string text_alt;
  vector<double> modell_final;
  double kriterium_final;
  ST::string text_final;

  ST::string tr_akt = "trace_on";

  vector<vector<double> > lambdavec;

  vector<unsigned> anfang;
  vector<unsigned> ende;
  vector<ST::string> names_fixed;
  vector<vector<ST::string> > names_nonparametric;

  initialise_lambdas(names_nonparametric,names_fixed,lambdavec,anfang,ende,number);

  unsigned i;
  unsigned j;

  /*         Kontrolle
  for(i=0;i<lambdavec.size();i++)
     {
     for(j=0;j<lambdavec[i].size();j++)
        genoptions_mult[0]->out(ST::doubletostring(lambdavec[i][j]) + "\n");
     }
  */

  unsigned k;
  unsigned zaehler = 0;
  bool fehler_randomslope = false;                        // überprüfen, dass Randomslopes
  for(i=0;i<fullcondp.size();i++)
     {
     if(fullcondp[i]->get_fctype() == MCMC::randomslopes)
       {
       for(j=1;j<names_fixed.size();j++)
          {
          if(names_fixed[j]==names_nonparametric[zaehler][0])
             {
             k = j;
             fehler_randomslope = true;
             }
          }
       }
     if(fullcondp[i]->get_fctype() != MCMC::variance &&
                    fullcondp[i]->get_fctype() != MCMC::fixed)
        zaehler = zaehler + 1;
     }
  if(fehler_randomslope==true)                 // nicht auch als fixe Effekte angegeben werden!
    {
    genoptions_mult[0]->outerror("\n\n ERROR: You must not put fixed effect " +
       names_fixed[k] + " in the model! \n");
    return true;
    }

  fullcond_alle = fullcondp;

  vector<vector<unsigned> > startindex;
  vector<vector<double> > startfix;
  startwerte(startmodel,lambdavec,names_fixed,startindex,startfix);

  if(likep_mult[0]->get_family() != "Gaussian")
    {
    vector<ST::string> title;
    title.push_back("");
    vector<double> modell_init;
    /*
    fixed_entfernen(fullcondp[0],modell_init,names_fixed,names_nonp,D,modelv);
    fullcond_entfernen(fullcondp,modell_init,anfang,ende,names_fixed,
          names_nonp,D,modelv);
    */
    for(i=1;i<names_fixed.size();i++)       // nur Kontrolle!
       modell_init.push_back(-1);
    for(i=0;i<lambdavec.size();i++)
      modell_init.push_back(lambdavec[i][floor(3*lambdavec[i].size()/4)]);

   for(i=0;i<modell_init.size();i++)
       genoptions_mult[0]->out(ST::doubletostring(modell_init[i]) + "\n");
    lambdas_update(modell_init);
    posteriormode(title,true);
    }

  options_text(number,lambdavec,names_fixed,names_nonparametric,startfix,startindex,name);

  ST::string path_tex = path + "_model_summary.tex";
  outtex.open(path_tex.strtochar());
  make_graphics(name,lambdavec,startindex);

  bool first = true;
  bool abbruch = false;

  ST::string pathmodels = path + "_models.raw";
  ST::string pathcrit = path + "_criterium.raw";
  outcriterium.open(pathcrit.strtochar());
  outcriterium << "step   " << criterion << endl;
  outmodels.open(pathmodels.strtochar());
  outmodels << "step   " << criterion << "   model" << endl;

  for(i=0;i<startindex.size();i++)
     {
     abbruch = single_stepwise(kriterium_alt,modell_alt,text_alt,
         anfang,ende,lambdavec,startindex[i],startfix[i],fullcond_alle,names_fixed,
         names_nonparametric,D,modelv,true);

     if(abbruch==true)
        return true;

     if(first)
       {
       first = false;
       kriterium_final = kriterium_alt;
       modell_final = modell_alt;
       text_final = text_alt;
       }

     else
       {
       if(kriterium_final>kriterium_alt)
         {
         kriterium_final = kriterium_alt;
         modell_final = modell_alt;
         text_final = text_alt;
         }
       }
     }

  vector<ST::string> title;
  title.push_back("");
  ST::string header = "  Final Model:";
  fullcondp = fullcond_alle;

  fixed_entfernen(fullcondp[0],modell_final,names_fixed,names_nonparametric,
                  D,modelv);
  fullcond_entfernen(fullcondp,modell_final,anfang,ende,names_fixed,
          names_nonparametric,D,modelv);
  lambdas_update(modell_final);
  tr_akt = "trace_on";
  maketext(header,modell_final,kriterium_final,text_final,false,tr_akt,false);
  genoptions_mult[0]->out("\n\n");
  kriterium_tex = kriterium_final;

  if(fine_tuning == true)
     abbruch = finetuning(modell_final,lambdavec,names_fixed,anfang,ende,D,modelv);

  if(abbruch==true)
    return true;

  fullcond_z = fullcondp;
  for(i=0;i<fullcond_z.size();i++)
     {
     if(fullcond_z[i]->get_fctype() != MCMC::variance)
        fullcond_z[i]->set_fcnumber(i);
     }

  posteriormode(title,false);  // Problem: linearer Prädiktor bei "true" standardisiert! Hier wird zurückgerechnet!
                               // danach nicht mehr compute_criterion() aufrufen!!!
  make_tex_end(path,modell_final,names_fixed);

  //            Files müssen wieder geschlossen werden!!!
  outtex.close();
  outcriterium.close();
  outmodels.close();

  // gibt Lambdas aus, damit man die richtig bestimmten Variablen zählen kann!
  /*
  ST::string zaehlername = "c:\\cprog\\stepwise\\simulation\\df_equidist\\zaehler_";
  zaehlername = zaehlername + likep_mult[0]->get_responsename() + ".ascii";
  //zaehlername = zaehlername + "_" + ST::inttostring(increment) + ".ascii";
  ofstream out(zaehlername.strtochar());
  ST::string beschriftung = "  krit   ";
  ST::string eintrag = "  " + ST::doubletostring(kriterium_final) + "   ";
  for(i=1;i<names_fixed.size();i++)
     beschriftung = beschriftung + names_fixed[i] + "   ";
  for(i=0;i<names_nonparametric.size();i++)
     beschriftung = beschriftung + names_nonparametric[i][0] + "   ";
  for(i=0;i<modell_final.size();i++)
     eintrag = eintrag + ST::doubletostring(modell_final[i]) + "   ";
  out << beschriftung << endl;
  out << eintrag << endl;
  */

  return false;
  }


bool STEPWISErun::single_stepwise(double & kriterium_alt,
         vector<double> & modell_alt, ST::string & text_alt,
         const vector<unsigned> & anfang, const vector<unsigned> & ende,
         const vector<vector<double> > & lambdavec,
         const vector<unsigned> & start, const vector<double> & startfix,
         const vector<FULLCOND*> & fullcond_fest,
         const vector<ST::string> & names_fixed,
         const vector<vector<ST::string> > & names_nonp, const datamatrix & D,
         const vector<ST::string> & modelv, const bool & tex)
  {

  vector<vector<vector<double> > > modellematrix;
  vector<double> modell_neu;
  vector<vector<double> > kriteriummatrix;
  double kriterium_neu;
  int steps_aktuell = 0;
  ST::string tr_akt = "trace_on";
  vector<vector<double> > startiteration;
  vector<double> kriteriumiteration;
//  vector<double> freiheitsgrade;                        //**

  unsigned i;
  fullcondp = fullcond_fest;
  unsigned z = 0;
  for(i=0;i<fullcondp.size();i++)
     {
     if(fullcondp[i]->get_fctype() == MCMC::fixed)
       {
//       freiheitsgrade.push_back(0);                     //**
       unsigned j;
       for(j=0;j<names_fixed.size()-1;j++)
          modell_neu.push_back(startfix[j]);
       }
     else if(fullcondp[i]->get_fctype() != MCMC::fixed
                      && fullcondp[i]->get_fctype() != MCMC::variance)
       {
//       freiheitsgrade.push_back(0);                     //**
       double lambda = lambdavec[z][start[z]];
       modell_neu.push_back(lambda);
       z = z + 1;
       }
     }
  modell_alt = modell_neu;
  startiteration.push_back(modell_alt);
  modellematrix.push_back(startiteration);
  vector<ST::string> title;
  title.push_back("");
  fixed_entfernen(fullcondp[0],modell_alt,names_fixed,names_nonp,D,modelv);
  fullcond_entfernen(fullcondp,modell_alt,anfang,ende,names_fixed,
          names_nonp,D,modelv);
  lambdas_update(modell_alt);
  posteriormode(title,true);

//  fgvector_besetzen(freiheitsgrade,modell_alt,names_fixed);       //**

  kriterium_alt = compute_criterion();

  if(tex==true)
     {
     kriterium_tex = kriterium_alt;
     make_predictor();
     }

  kriteriumiteration.push_back(kriterium_alt);
  kriteriummatrix.push_back(kriteriumiteration);
  // ENDE: startmodel

  kriterium_neu = kriterium_alt;
  outcriterium << steps_aktuell << "   " << kriterium_neu << endl;
  outmodels << steps_aktuell << "   " << kriterium_neu << "   ";
  ST::string header;
  bool fertig = false;    // überprüft, ob es noch nicht gerechnete Modelle gibt
  ST::string text_neu;
  bool eins = true;

  while(kriterium_neu<=kriterium_alt && fertig==false && steps_aktuell<steps)
       {
       //outcriterium << steps_aktuell << "   " << kriterium_neu << endl;

       steps_aktuell = steps_aktuell + 1;
       vector<ST::string> textiteration;
       fertig = true;
       kriterium_alt = kriterium_neu;
       modell_alt = modell_neu;
       header = "  Startmodel:";
       bool neutext;
       if(eins == true)
         {
         eins = false;
         neutext = true;  // Modell rausschreiben !!
         }
       else
         {
         if(trace == "trace_off")
           tr_akt = "trace_off";
         neutext = false;
         }
       maketext(header, modell_alt,kriterium_alt,text_neu,neutext,          //criterion,
                tr_akt,neutext); //,outmodels);
       text_alt = text_neu;

       vector<vector<double> > modeliteration;
       double kriterium;
       vector<double> kriteriumiteration2;
       unsigned zaehlen = names_fixed.size()-1;
       for(i=0;i<fullcond_fest.size();i++)
          {
          if(fullcond_fest[i]->get_fctype() == MCMC::fixed)
            {
            unsigned f;
            for(f=1;f<names_fixed.size();f++)
               {
               modell_neu = modell_alt;
               if(modell_alt[f-1]==-1)
                 modell_neu[f-1]= 0;
               else if(modell_alt[f-1]==0)
                 modell_neu[f-1] = -1;
                 newmodel_rechnen(fertig,modell_neu,kriteriumiteration2,
                        modeliteration,modellematrix,fullcond_fest,anfang,ende,
                        textiteration,        //criterion,trace,
                        names_fixed,names_nonp,D,modelv); //,outmodels);
               }
            }
          else if(fullcond_fest[i]->get_fctype() != MCMC::fixed
                           && fullcond_fest[i]->get_fctype() != MCMC::variance)
            {
            unsigned sch;
            for(sch=1;sch<=increment;sch++)
               {
               modell_neu = modell_alt;
               fullcondp = fullcond_fest;
               fixed_entfernen(fullcondp[0],modell_neu,names_fixed,names_nonp,
                            D,modelv);
               fullcond_entfernen(fullcondp,modell_neu,anfang,ende,names_fixed,
                               names_nonp,D,modelv);
               lambdas_update(modell_neu);
               bool lambda_exist;    // zum Überprüfen, ob neues Lambda im Vektor enthalten ist
               unsigned index = search_lambdaindex(modell_alt[zaehlen],
                                lambdavec[zaehlen-names_fixed.size()+1],lambda_exist);
               lambda_exist = false;
               if(index < lambdavec[zaehlen-names_fixed.size()+1].size()-sch)
                 lambda_exist = true;
               if(lambda_exist==true)
                 {
                 modell_neu[zaehlen] = lambdavec[zaehlen-names_fixed.size()+1][index+sch];
                 newmodel_rechnen(fertig,modell_neu,kriteriumiteration2,
                    modeliteration,modellematrix,fullcond_fest,anfang,ende,
                    textiteration,//criterion,trace,
                    names_fixed,names_nonp,D,modelv);  //,outmodels);
                 }

               lambda_exist = false;
               modell_neu = modell_alt;
               fullcondp = fullcond_fest;
               fixed_entfernen(fullcondp[0],modell_neu,names_fixed,names_nonp,
                            D,modelv);
               fullcond_entfernen(fullcondp,modell_neu,anfang,ende,names_fixed,
                               names_nonp,D,modelv);
               lambdas_update(modell_neu);
               if(index >= sch)
                 lambda_exist = true;
               if(lambda_exist==true)
                 {
                 modell_neu[zaehlen] = lambdavec[zaehlen-names_fixed.size()+1][index-sch];
                 newmodel_rechnen(fertig,modell_neu,kriteriumiteration2,
                    modeliteration,modellematrix,fullcond_fest,anfang,ende,
                    textiteration,//criterion,trace,
                    names_fixed,names_nonp,D,modelv);  //,outmodels);
                 }
               }
            zaehlen = zaehlen + 1;
            }
          }
       if(fertig==false)
         {
         kriterium = kriteriumiteration2[0];
         unsigned j;
         for(j=0;j<kriteriumiteration2.size();j++)  //berechnet den besten Wert
            {
            if(kriteriumiteration2[j]<=kriterium)
              {
              kriterium = kriteriumiteration2[j];
              modell_neu = modeliteration[j];
              text_neu = textiteration[j];
              }
            }
         kriterium_neu = kriterium;
         outcriterium << steps_aktuell << "   " << kriterium_neu << endl;
         outmodels << steps_aktuell << "   " << kriterium_neu << "   ";
         modellematrix.push_back(modeliteration);
         kriteriummatrix.push_back(kriteriumiteration2);

         header = "\n\nBest Model of this iteration:";
         fullcondp = fullcond_fest;
         fixed_entfernen(fullcondp[0],modell_neu,names_fixed,names_nonp,
                         D,modelv);
         fullcond_entfernen(fullcondp,modell_neu,anfang,ende,names_fixed,
                            names_nonp,D,modelv);
         lambdas_update(modell_neu);
         neutext = false;
         maketext(header,modell_neu,kriterium_neu,text_neu,neutext,         //criterion,
                  trace,true); //outmodels);                   // Modell rausschreiben!!

         if(steps_aktuell==steps)
           {
           if(kriterium_alt>kriterium_neu)
             {
             kriterium_alt = kriterium_neu;
             modell_alt = modell_neu;
             text_alt = text_neu;
             }
           }
         }
       else
         {
         if(trace == "trace_on")
           {
           genoptions_mult[0]->out("\n\n");
           genoptions_mult[0]->out("\n\n");
           genoptions_mult[0]->out("  There are no new models for this iteration! \n");
           }
         }
       if(trace == "trace_on")
         {
         genoptions_mult[0]->out("\n\n");
         genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
         genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
         }

       #if defined(BORLAND_OUTPUT_WINDOW)
       Application->ProcessMessages();

       if (Frame->stop)
         {
         //break;
         genoptions_mult[0]->out("\n");
         genoptions_mult[0]->out(
         "STEPWISE PROCEDURE TERMINATED BY USER BREAK\n");
         genoptions_mult[0]->out("\n");
         genoptions_mult[0]->out("Estimation results: none\n");
         genoptions_mult[0]->out("\n");
         return true;
         }

       if (Frame->pause)
         {
         genoptions_mult[0]->out("\n");
         genoptions_mult[0]->out("STEPWISE PROCEDURE PAUSED\n");
         genoptions_mult[0]->out("Click CONTINUE to proceed\n");
         genoptions_mult[0]->out("\n");

         while (Frame->pause)
           {
           Application->ProcessMessages();
           }

         genoptions_mult[0]->out("STEPWISE PROCEDURE CONTINUED\n");
         genoptions_mult[0]->out("\n");
         }
         
       #elif defined(JAVA_OUTPUT_WINDOW)
       bool stop = genoptions_mult[0]->adminb_p->breakcommand();
       if(stop)
         {
         genoptions_mult[0]->out("\n");
         genoptions_mult[0]->out(
         "STEPWISE PROCEDURE TERMINATED BY USER BREAK\n");
         genoptions_mult[0]->out("\n");
         genoptions_mult[0]->out("Estimation results: none\n");
         genoptions_mult[0]->out("\n");
         //break;
         return true;
         }
       #endif

       }

  header = "  Final Model:";
  //bool neutext = false;
  tr_akt = "trace_on";
  maketext(header,modell_alt,kriterium_alt,text_alt,false,tr_akt,false); //,outmodels);
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("  Used number of iterations: " + ST::inttostring(steps_aktuell));
  genoptions_mult[0]->out("\n\n");
  genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
  genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");

  return false;
  }


void STEPWISErun::modelcomparison(bool & s, const vector<double> & m,
                   const vector<vector<vector<double> > > & mmatrix)
  {
  int x = mmatrix.size()-1;
  while(x>=0 && s==false)
       {
       int y = mmatrix[x].size()-1;
       while(y>=0 && s==false)
            {
            if(mmatrix[x][y]==m)
              s = true;
            y = y - 1;
            }
       x = x - 1;
       }
  }

void STEPWISErun::maketext(const ST::string & h, const vector<double> & m,
                          const double & a, ST::string & text,
                          const bool & neutext, const ST::string & tr,
                          const bool & datei) //, ofstream & outmodels)
  {
  if(tr == "trace_on")
    {
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out(h);
    }
  ST::string modeltext;
  if(neutext==true)
    {
    modeltext = "  " + likep_mult[0]->get_responsename() + " = ";
    unsigned i;
    for(i=0;i<fullcondp.size();i++)
       {
       if(fullcondp[i]->get_fctype() != MCMC::fixed
                   && fullcondp[i]->get_fctype() != MCMC::variance)
         modeltext = modeltext + " + " + fullcondp[i]->get_effect();
       else if(fullcondp[i]->get_fctype() == MCMC::fixed)
          {
          //unsigned z;
          //for(z=1;z<fullcondp[i]->get_datanames().size();z++)
           //   modeltext = modeltext + " + " + fullcondp[i]->get_datanames()[z];
          modeltext = modeltext + fullcondp[i]->get_effect();
          }
       }
    text = modeltext;
    }
  else
    {
    modeltext = text;
    }
  if(tr == "trace_on")
    {
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out(modeltext);
    genoptions_mult[0]->out("\n " + criterion + " = " + ST::doubletostring(a,6));
    }
  if(datei==true)
    outmodels << modeltext << endl;
//  ST::string lambdatext = "const";
//  for(unsigned z=0;z<m.size();z++)
//     {
//     lambdatext = lambdatext + " + " + ST::doubletostring(m[z]);
//     }
//  genoptions_mult[0]->out(lambdatext);
  }

void STEPWISErun::lambdas_update(const vector<double> & m)
  {
  vector<double> modell_kurz;
  unsigned k;
  for(k=0;k<m.size();k++)
     {
     if(m[k]>0)
       modell_kurz.push_back(m[k]);
     }
  unsigned j = 0;
  unsigned i;
  for(i=0;i<fullcondp.size();i++)
     {
     if(fullcondp[i]->get_fctype()!=MCMC::fixed
                 && fullcondp[i]->get_fctype()!=MCMC::variance)
       {
       fullcondp[i]->update_stepwise(modell_kurz[j]);
       j = j + 1;
       }
     }
  }

void STEPWISErun::newmodel_rechnen(bool & f, const vector<double> & m,
    vector<double> & a, vector<vector<double> > & mi,
    const vector<vector<vector<double> > > & mmatrix,
    const vector<FULLCOND*> & alle, const vector<unsigned> & beg,
    const vector<unsigned> & ende, vector<ST::string> & textit,
    //const ST::string & criterion, const ST::string & tr,
    const vector<ST::string> & names_fixed,
    const vector<vector<ST::string> > & names_nonp, const datamatrix & D,
    const vector<ST::string> & modelv) //, ofstream & outmodels)
  {
  bool schon_gerechnet = false;  // überprüft, ob Modell schon mal gerechnet wurde

  modelcomparison(schon_gerechnet, m, mmatrix);

  if(schon_gerechnet==false)
    {
    f = false;
    mi.push_back(m);
    fullcondp = alle;
    fixed_entfernen(fullcondp[0], m, names_fixed, names_nonp, D, modelv);
    fullcond_entfernen(fullcondp, m, beg, ende, names_fixed,
          names_nonp, D, modelv);
    lambdas_update(m);
    vector<ST::string> title;
    title.push_back("");
    posteriormode(title,true);      //keine Ausgabe
    double kriterium = compute_criterion();
    ST::string header = "  Trial: ";
    ST::string text;
    //bool neutext = true;
    maketext(header,m,kriterium,text,true,trace,false); //,outmodels);
    textit.push_back(text);
    a.push_back(kriterium);
    }
  }

unsigned STEPWISErun::search_lambdaindex(const double & m,
                       const vector<double> lam, bool & b) const
  {
  unsigned index = 0;
  double lambda = m;
  unsigned j=0;
  b = false;
  while(j<lam.size() && b==false)
       {
       if(lambda==lam[j])
         {
         index = j;
         b = true;
         }
       j++;
       }
  return index;
  }

unsigned STEPWISErun::search_lambdastartindex(const double & start,
                       const vector<double> & lambdas) const
  {
  bool gefunden = false;
  unsigned index = search_lambdaindex(start, lambdas, gefunden);

  if(gefunden==false)
    {
    vector<double> diff;
    unsigned i;
    for(i=0;i<lambdas.size();i++)
       diff.push_back(fabs(lambdas[i]-start));
    double irgendwas = diff[0];
    for(i=1;i<diff.size();i++)
       {
       if(diff[i]<=irgendwas)
         {
         irgendwas = diff[i];
         index = i;
         }
       }
    }

  return index;
  }

void STEPWISErun::startwerte(const ST::string & startmodel,
          const vector<vector<double> > & lambdavec,
          const vector<ST::string> & names_fixed,
          vector<vector<unsigned> > & startindex,
          vector<vector<double> > & startfix)
  {
  unsigned i;
  unsigned j;
  if(startmodel == "empty" || startmodel == "both")
    {
    vector<unsigned> indexhelp;
    vector<double> fixhelp;
    unsigned z = 0;
    for(i=0;i<fullcondp.size();i++)
       {
       if(fullcondp[i]->get_fctype() != MCMC::fixed
                   && fullcondp[i]->get_fctype() != MCMC::variance)
         {
         indexhelp.push_back(lambdavec[z].size()-1);
         z = z + 1;
         }
       else if(fullcondp[i]->get_fctype() == MCMC::fixed)
         {
         for(j=1;j<names_fixed.size();j++)
            fixhelp.push_back(0);
         }
       }
    startindex.push_back(indexhelp);
    startfix.push_back(fixhelp);
    }
  if(startmodel == "both" || startmodel == "full")
    {
    vector<unsigned> indexhelp;
    vector<double> fixhelp;
    for(i=0;i<fullcondp.size();i++)
       {
       if(fullcondp[i]->get_fctype() != MCMC::fixed
                   && fullcondp[i]->get_fctype() != MCMC::variance)
         indexhelp.push_back(0);
       else if(fullcondp[i]->get_fctype() == MCMC::fixed)
         {
         for(j=1;j<names_fixed.size();j++)
            fixhelp.push_back(-1);
         }
       }
    startindex.push_back(indexhelp);
    startfix.push_back(fixhelp);
    }
  else if(startmodel == "userdefined")
    {
    vector<unsigned> indexhelp;
    vector<double> fixhelp;
    unsigned z = 0;
    for(i=0;i<fullcondp.size();i++)
       {
       if(fullcondp[i]->get_fctype() != MCMC::fixed
                   && fullcondp[i]->get_fctype() != MCMC::variance)
         {
         double start = fullcondp[i]->get_lambdastart();
         unsigned index = search_lambdastartindex(start, lambdavec[z]);
         indexhelp.push_back(index);
         z = z + 1;
         }
       else if(fullcondp[i]->get_fctype() == MCMC::fixed)
         {
         for(j=1;j<names_fixed.size();j++)
            fixhelp.push_back(-1);                       //fehlt: vom Benutzer angeben lassen!!!
         }
       }
    startindex.push_back(indexhelp);
    startfix.push_back(fixhelp);
    }
  }

void STEPWISErun::fullcond_entfernen(vector<FULLCOND*> & fullc,
         const vector<double> & m, const vector<unsigned> & beg,
         const vector<unsigned> & ende, const vector<ST::string> & names_fixed,
         const vector<vector<ST::string> > & names_nonp, const datamatrix & D,
         const vector<ST::string> & modelv)
  {
  vector<FULLCOND*> fullcond_neu;
  unsigned j = names_fixed.size()-1;
  unsigned i;
  for(i=0;i<fullc.size();i++)
     {
     if(fullc[i]->get_fctype()!=MCMC::fixed
             && fullc[i]->get_fctype()!=MCMC::variance)
       {
       if(m[j]>0)
         {
         fullcond_neu.push_back(fullc[beg[j-names_fixed.size()+1]]);
         fullcond_neu.push_back(fullc[ende[j-names_fixed.size()+1]]);
         }
       else if(m[j]==0)
       
         {
         fullc[i]->reset_effect(0);
         }
       else if(m[j] == -1)
         {
         fullc[i]->reset_effect(0);
         fullc[0]->include_effect(names_nonp[j-names_fixed.size()+1],
                                  fullc[i]->get_data_forfixedeffects());
         }
       j = j + 1;
       }
     else if(fullc[i]->get_fctype() == MCMC::fixed)
       fullcond_neu.push_back(fullc[i]);
     }
  fullc = fullcond_neu;
  end[0] = fullc.size()-1;
  }

void STEPWISErun::fixed_entfernen(FULLCOND* & fullc,
         const vector<double> &  modell, const vector<ST::string> & names_fixed,
         const vector<vector<ST::string> > & names_nonp, const datamatrix & D,
         const vector<ST::string> & modelv)
  {
  unsigned z;
  for(z=0;z<names_fixed.size()-1;z++)
     {
     if(modell[z]==0)
       {
       bool raus = false;
       unsigned j = 1;
       while(j<fullc->get_datanames().size() && raus==false)
            {
            if(fullc->get_datanames()[j]==names_fixed[z+1])
              {
              raus = true;

              fullc->reset_effect(j);
              }
            j = j + 1;
            }
       }
     else if(modell[z]==-1)
       {
       unsigned i = 1;
       bool rein = false;
       while(i<fullc->get_datanames().size() && rein==false)
            {
            if(fullc->get_datanames()[i]==names_fixed[z+1])
              rein = true;
            i = i + 1;
            }
       if(rein==false)
         {
         bool gefunden = false;
         i = 0;
         while(i<modelv.size() && gefunden==false)
              {
              if(names_fixed[z+1]==modelv[i])
                gefunden = true;
              i = i + 1;
              }

         vector<ST::string> help_name;
         help_name.push_back(names_fixed[z+1]);
         fullc->include_effect(help_name,datamatrix(D.getCol(i-1)));
         }
       }
     }

  for(z=names_fixed.size()-1;z<modell.size();z++)
     {
     bool gefunden = false;
     unsigned i = 1;
     while(i<fullc->get_datanames().size() && gefunden==false)
        {
        if(fullc->get_datanames()[i]==names_nonp[z-names_fixed.size()+1][0])
           {
           gefunden = true;
           fullc->reset_effect(i);
           }
           i = i + 1;
        }
     if(gefunden==true && names_nonp[z-names_fixed.size()+1].size()>1)
        {
        unsigned j;
        for(j=1;j<names_nonp[z-names_fixed.size()+1].size();j++)
           {
           bool gefunden2 = false;
           i = 1;
           while(i<fullc->get_datanames().size() && gefunden2==false)
              {
              if(fullc->get_datanames()[i]==names_nonp[z-names_fixed.size()+1][j])
                {
                gefunden2 = true;
                fullc->reset_effect(i);
                }
              i = i + 1;
              }
           }
        }
     }
  }

void STEPWISErun::options_text(//const ST::string & criterion, const int & steps,
         const int & number, //const ST::string & startmodel,
         const vector<vector<double> > & lambdavec,
         const vector<ST::string> & names_fixed,
         const vector<vector<ST::string> > & names_nonp,
         const vector<vector<double> > & startfix,
         const vector<vector<unsigned> > & startindex, const ST::string & name)
  {
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("STEPWISE OBJECT " + name + ": stepwise procedure \n");
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("GENERAL OPTIONS: \n");
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("  Performance criterion: " + criterion + " \n");
  genoptions_mult[0]->out("  Maximal number of iterations: " + ST::inttostring(steps) + "\n");
  //genoptions_mult[0]->out("  Number of different smoothing parameters: "
  //   + ST::inttostring(number) + "\n");
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("  RESPONSE DISTRIBUTION: \n");
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("  Family: " + likep_mult[0]->get_family() + "\n");
  genoptions_mult[0]->out("  Number of observations: "
     + ST::inttostring(likep_mult[0]->get_nrobs()) + "\n");
  //genoptions_mult[0]->out("  Number of observations with positive weights: " + );
  //genoptions_mult[0]->out("  Response function: " + );
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("OPTIONS FOR STEPWISE PROCEDURE: \n");
  unsigned i;
  for(i=1;i<names_fixed.size();i++)
     {
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  OPTIONS FOR FIXED EFFECTS TERM: "
        + names_fixed[i] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Prior: diffuse prior \n");
     unsigned j;
     for(j=0;j<startfix.size();j++)
        if(startfix[j][i-1] == 0)
           genoptions_mult[0]->out("  Startvalue of the "
              + ST::doubletostring(j+1) + ". startmodel is \"effect excluded\" \n");
        else
           genoptions_mult[0]->out("  Startvalue of the "
              + ST::doubletostring(j+1) + ". startmodel is the fixed effect \n");
     }
  unsigned z = 0;
  for(i=1;i<fullcondp.size();i++)
     {
     if(fullcondp[i]->get_fctype() != MCMC::fixed
                 && fullcondp[i]->get_fctype() != MCMC::variance)
       {
       genoptions_mult[0]->out("\n");
       genoptions_mult[0]->out("  OPTIONS FOR NONPARAMETRIC TERM: "
          + names_nonp[z][0] + "\n");
       genoptions_mult[0]->out("\n");
       //genoptions_mult[0]->out("Prior: " + fullcondp[i]->get_priorassumptions());
       if(fullcondp[i]->get_lambdamin()!=0 && fullcondp[i]->get_lambdamin()!=-1)
          {
          genoptions_mult[0]->out("  Minimal value for the smoothing parameter: "
             + ST::doubletostring(fullcondp[i]->get_lambdamin()) + "\n");
          fullcondp[i]->update_stepwise(fullcondp[i]->get_lambdamin());
          if(fullcondp[i]->get_df_equidist()==false && fullcondp[i]->get_lambdamin_opt()==false)
            genoptions_mult[0]->out("  This is equivalent to degrees of freedom: "
               + ST::doubletostring(fullcondp[i]->compute_df(),6) + "\n");
          else
            genoptions_mult[0]->out("  This is equivalent to degrees of freedom: approximately "
               + ST::doubletostring(fullcondp[i]->get_df_lambdamin()) + ", exact "
               + ST::doubletostring(fullcondp[i]->compute_df(),6) + "\n");
          }
       if(fullcondp[i]->get_lambdamax()!=0 && fullcondp[i]->get_lambdamax()!=-1)
          {
          genoptions_mult[0]->out("  Maximal value for the smoothing parameter: "
             + ST::doubletostring(fullcondp[i]->get_lambdamax(),6) + "\n");
          fullcondp[i]->update_stepwise(fullcondp[i]->get_lambdamax());
          if(fullcondp[i]->get_df_equidist()==false && fullcondp[i]->get_lambdamax_opt()==false)
            genoptions_mult[0]->out("  This is equivalent to degrees of freedom: "
               + ST::doubletostring(fullcondp[i]->compute_df(),6) + "\n");
          else
            genoptions_mult[0]->out("  This is equivalent to degrees of freedom: approximately "
               + ST::doubletostring(fullcondp[i]->get_df_lambdamax()) + ", exact "
               + ST::doubletostring(fullcondp[i]->compute_df(),6) + "\n");
          }
       if(fullcondp[i]->get_df_equidist()==true)
         genoptions_mult[0]->out("  Number of different smoothing parameters with equidistant degrees of freedom: "
            + ST::doubletostring(fullcondp[i]->get_number()) + "\n");
       else
         genoptions_mult[0]->out("  Number of different smoothing parameters on a logarithmic scale: "
            + ST::doubletostring(fullcondp[i]->get_number()) + "\n");
       unsigned j;
       for(j=0;j<startindex.size();j++)
          {
          if(lambdavec[z][startindex[j][z]]==0)
             genoptions_mult[0]->out("  Startvalue of the "
                + ST::doubletostring(j+1) + ". startmodel is \"effect excluded\" \n");
          else if(lambdavec[z][startindex[j][z]]==-1)
             genoptions_mult[0]->out("  Startvalue of the "
                + ST::doubletostring(j+1) + ". startmodel is the fixed effect \n");
          else
             {
             genoptions_mult[0]->out("  Startvalue of the smoothing parameter for the "
                + ST::doubletostring(j+1) + ". startmodel: "
                + ST::doubletostring(lambdavec[z][startindex[j][z]],6) + "\n");
             fullcondp[i]->update_stepwise(lambdavec[z][startindex[j][z]]);
             genoptions_mult[0]->out("  This is equivalent to degrees of freedom: "
                + ST::doubletostring(fullcondp[i]->compute_df(),6) + "\n");
             }
          }
       z = z + 1;
       }
     }
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("STEPWISE PROCEDURE STARTED \n");
  genoptions_mult[0]->out("\n");
  }


bool STEPWISErun::finetuning(vector<double> & modell, vector<vector<double> > & lambdavec,
                 vector<ST::string> & names_fixed, const vector<unsigned> & beg,
                 const vector<unsigned> & end, const datamatrix & D,
                 const vector<ST::string> & modelv)
  {
  genoptions_mult[0]->out("  BEGIN FINE-TUNING:");
  genoptions_mult[0]->out("\n\n");
  vector<FULLCOND*> fullcond_neu;
  double number = 11;
  unsigned i = 0;
  unsigned z = names_fixed.size()-1;
  unsigned anzahl = 0;
  vector<double> lambdamaxe;
  for(i=0;i<fullcond_alle.size();i++)
     {
     if(fullcond_alle[i]->get_fctype() != MCMC::fixed
                   && fullcond_alle[i]->get_fctype() != MCMC::variance)
        {
        if(modell[z]!=0)
           {
           if(modell[z]!=-1)
              anzahl = anzahl + 1;
           fullcond_neu.push_back(fullcond_alle[beg[z-names_fixed.size()+1]]);
           if(i<fullcond_alle.size()-1 && fullcond_alle[i+1]->get_fctype() == MCMC::variance)
              fullcond_neu.push_back(fullcond_alle[end[z-names_fixed.size()+1]]);

           bool lambda_exist;
           unsigned index = search_lambdaindex(modell[z],lambdavec[z-names_fixed.size()+1],lambda_exist);
           double lambdastart = lambdavec[z-names_fixed.size()+1][index];
           double lambdamin;
           double lambdamax;
           if(lambdavec[z-names_fixed.size()+1].size()>2)
              {
              lambda_exist = false;
              int j = 2;
              while(lambda_exist==false && j>0)
                 {
		 if(index>=j)
                    {
                    if(lambdavec[z-names_fixed.size()+1][index-j]!=0 &&
                                   lambdavec[z-names_fixed.size()+1][index-j]!=-1)
                       lambda_exist = true;
                    }
                 j = j - 1;
		 }
              if(lambda_exist==true)
                lambdamin = lambdavec[z-names_fixed.size()+1][index-j-1];
              else
                lambdamin = fullcond_alle[i]->get_lambdamin()/10;

              lambda_exist = false;
              j = 2;
              while(lambda_exist==false && j>0)
                 {
		 if(index<lambdavec[z-names_fixed.size()+1].size()-j)
		    {
		    if(lambdavec[z-names_fixed.size()+1][index+j]!=0 &&
                                 lambdavec[z-names_fixed.size()+1][index+j]!=-1)
                       lambda_exist = true;
		    }
		 j = j - 1;
 		 }
	      if(lambda_exist==true)
                lambdamax = lambdavec[z-names_fixed.size()+1][index+j+1];     //**!!
              else if(lambdavec[z-names_fixed.size()+1][index]==-1)
                lambdamax = 1.5*lambdavec[z-names_fixed.size()+1][index-1] -
                              0.5*lambdavec[z-names_fixed.size()+1][index-2];
              else
                lambdamax = 1.5*lambdavec[z-names_fixed.size()+1][index] -
                              0.5*lambdavec[z-names_fixed.size()+1][index-1];

              lambdamaxe.push_back(lambdamax);
              if(lambdastart!=-1)
                fullcond_alle[i]->set_stepwise_options(lambdastart,lambdastart,lambdamin,
                                fullcond_alle[i]->get_forced(),0,0,false,false,number,false);
              else
                fullcond_alle[i]->set_stepwise_options(lambdastart,lambdamax,lambdamin,
                                fullcond_alle[i]->get_forced(),0,0,false,false,21,false);
              }
           }
	z = z + 1;
        }
     else if(fullcond_alle[i]->get_fctype() == MCMC::fixed)
        fullcond_neu.push_back(fullcond_alle[i]);
     }
  //fullcond_alle = fullcond_neu;
  fullcondp = fullcond_neu;
  fullcond_alle = fullcondp; // hier anders, wenn es für als fix ausgew. Var. wieder mehr Mögl. geben soll

  vector<unsigned> anfang;
  vector<unsigned> ende;
  vector<vector<ST::string> > names_nonp;
  vector<vector<double> > lambdavec_ursp;
  vector<vector<double> > lambdavec_fine;
  vector<ST::string> namesfixed_neu;
  initialise_lambdas(names_nonp,namesfixed_neu,lambdavec_ursp,anfang,ende,number);

  if(anzahl!=names_nonp.size())
     {
     vector<ST::string> namesfixed_fine;
     unsigned j;
     namesfixed_fine.push_back(namesfixed_neu[0]);
     for(i=1;i<namesfixed_neu.size();i++)
        {
        bool gefunden = false;
        j = 0;
        while(j<names_nonp.size() && gefunden==false)
           {
           if(names_nonp[j][0]==namesfixed_neu[i])
              gefunden = true;
           else if(names_nonp[j].size()>1)
              gefunden = true;
           j = j + 1;
           }
        if(gefunden==false)
           namesfixed_fine.push_back(namesfixed_neu[i]);
        }
     names_fixed = namesfixed_fine;
     }
  else
     names_fixed = namesfixed_neu;

  z = 0;
  unsigned w = 0;
  for(i=0;i<fullcondp.size();i++)
     {
     if(fullcondp[i]->get_fctype() != MCMC::fixed
                   && fullcondp[i]->get_fctype() != MCMC::variance)
        {
        if(i<fullcondp.size()-1 && fullcondp[i+1]->get_fctype() == MCMC::variance)
          {
          if(fullcondp[i]->get_lambdastart()!=-1)
            {
            fullcondp[i]->set_stepwise_options(fullcondp[i]->get_lambdastart(),lambdamaxe[w],
                         fullcondp[i]->get_lambdastart(),fullcondp[i]->get_forced(),0,0,false,false,number,false);
            vector<double> untervector;
            unsigned j = 0;
            for(j=0;j<lambdavec_ursp[z].size()-3;j++)
              untervector.push_back(lambdavec_ursp[z][j]);
            if(fullcondp[i]->get_forced()==true)
              untervector.push_back(lambdavec_ursp[z][lambdavec_ursp[z].size()-3]);
            fullcondp[i]->compute_lambdavec(untervector,number);
            lambdavec_fine.push_back(untervector);
            }
          else
            lambdavec_fine.push_back(lambdavec_ursp[z]);
          w = w + 1;
          }
        else
          lambdavec_fine.push_back(lambdavec_ursp[z]);
        z = z + 1;
        }
     }

  vector<vector<unsigned> > startindex;
  vector<vector<double> > startfix;
  startwerte("userdefined",lambdavec_fine,names_fixed,startindex,startfix);

  ST::string text;
  double kriterium;

  bool abbruch = false;

  abbruch = single_stepwise(kriterium,modell,text,anfang,ende,lambdavec_fine,
              startindex[0],startfix[0],fullcond_neu,names_fixed,names_nonp,D,
              modelv,false);

  if(abbruch==true)
    return true;

  ST::string header = "  Final Model after Fine-Tuning:";
  fullcondp = fullcond_neu;
  fixed_entfernen(fullcondp[0],modell,names_fixed,names_nonp,
                  D,modelv);
  fullcond_entfernen(fullcondp,modell,anfang,ende,names_fixed,
                     names_nonp,D,modelv);
  lambdas_update(modell);
  ST::string tr_akt = "trace_on";
  maketext(header,modell,kriterium,text,false,tr_akt,false);
  kriterium_tex = kriterium;
  genoptions_mult[0]->out("\n\n");

  return false;
  }


void STEPWISErun::initialise_lambdas(vector<vector<ST::string> > & names_nonp,
                 vector<ST::string> & names_fixed, vector<vector<double> > & lambdavec,
                 vector<unsigned> & anfang, vector<unsigned> & ende, const int & number)
  {
  unsigned i;
  for(i=0;i<fullcondp.size();i++)
     {
     if(fullcondp[i]->get_fctype() != MCMC::fixed
                   && fullcondp[i]->get_fctype() != MCMC::variance)
       {
       anfang.push_back(i);
       ende.push_back(i+1);
       int nummer = number;
       if(fullcondp[i]->get_data_forfixedeffects().cols() > 1)  //bei Faktor-Variablen
           names_nonp.push_back(fullcondp[i]->get_datanames());
       else
          {
          vector<ST::string> names_help;
          names_help.push_back(fullcondp[i]->get_datanames()[fullcondp[i]->get_datanames().size()-1]);
          names_nonp.push_back(names_help);
          nummer = fullcondp[i]->get_number();
          if(nummer==0)
             nummer = number;
          fullcondp[i]->set_stepwise_options(fullcondp[i]->get_lambdastart(),fullcondp[i]->get_lambdamax(),
                fullcondp[i]->get_lambdamin(),fullcondp[i]->get_forced(),fullcondp[i]->get_df_lambdamax(),
                fullcondp[i]->get_df_lambdamin(),fullcondp[i]->get_lambdamax_opt(),fullcondp[i]->get_lambdamin_opt(),
                nummer,fullcondp[i]->get_df_equidist());
          }
       vector<double> untervector;
       fullcondp[i]->compute_lambdavec(untervector,nummer);
       lambdavec.push_back(untervector);
       }
     else if(fullcondp[i]->get_fctype() == MCMC::fixed)
       names_fixed = fullcondp[i]->get_datanames();
     }
  }

//------------------------------------------------------------------------------
//-------TEX-File---------------------------------------------------------------
//------------------------------------------------------------------------------

void STEPWISErun::make_graphics(const ST::string & name,
       vector<vector<double> > & lambdavec, vector<vector<unsigned> > & startindex)
  {
  ST::string title = "STEPWISEREG OBJECT " + name + ": stepwise procedure";
    //erzeugt den Kopf des Tex-Files
  outtex << "\\documentclass[a4paper, 12pt]{article}" << endl
         << "\n" << "\\usepackage{graphicx}" << endl
         << "\\parindent0em" << endl
         << "\\textheight22cm \\textwidth15cm \\oddsidemargin0.5cm" << endl
         << "\n\\begin{document}" << endl
         << "\\begin{center}" << endl
         << "\\LARGE{\\bf " << title << "}"
         << endl << "\\end{center} \n\\vspace{1cm}" << endl;

  make_model();

  make_options();

  make_prior(lambdavec,startindex);

  outtex << "\n\\noindent {\\bf \\large Start Predictor";
  if(startindex.size()>1)
     outtex << "s";
  outtex << ":}\\\\" << endl;

  }

void STEPWISErun::make_tex_end(ST::string & path, const vector<double> & modell,
                               const vector<ST::string> & names_f)
  {
  ST::string path_batch = path + "_graphics.prg";
  ST::string path_splus = path +  "_splus.txt";
  //ST::string path_stata = path +  "_stata.do";
  ST::string path_tex = path + "_model_summary.tex";

  outtex << "\n\\noindent {\\bf \\large Final Predictor:}\\\\" << endl;
  make_predictor();

  unsigned j;
  unsigned z = names_f.size()-1;
  outtex << "\n\\noindent {\\bf \\large Final Properties:}\\\\ \n\\\\" << endl;
  for(j=1;j<fullcond_alle.size();j++)
     {
     if(fullcond_alle[j]->get_fctype() != MCMC::fixed
                 && fullcond_alle[j]->get_fctype() != MCMC::variance)
       {
       if(modell[z]!=0 && modell[z]!=-1)
         {
         vector<ST::string> prior = fullcond_alle[j]->get_priorassumptions();
         outtex << prior[0] << "\\\\" << endl
                << "smoothing parameter: $\\lambda = "
                << ST::doubletostring(modell[z],6) << " \\,\\, \\hat{=} \\,\\, \\mbox{df} = "
                << ST::doubletostring(fullcond_alle[j]->compute_df(),6)
                << "$ \\\\ \n\\\\" << endl;
         }
       z = z + 1;
       }
     }

  vector<ST::string> distr_results;
  distr_results = likep_mult[0]->get_results_latex();
  unsigned i;
  for (i=0;i<distr_results.size();i++)
     {
     outtex << distr_results[i] << endl;
     }

  make_fixed_table();

    // Pfade der Files
    //werden im BayesX-Output angegeben
  genoptions_mult[0]->out("  Files of model summary: \n" , true);
  genoptions_mult[0]->out("\n");

  make_plots(path_batch,path_splus);

  genoptions_mult[0]->out("  --------------------------------------------------------------------------- \n");
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("  Latex file of model summary is stored in file \n");
  genoptions_mult[0]->out("  " + path_tex + "\n");
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("  --------------------------------------------------------------------------- \n");
  genoptions_mult[0]->out("\n");

  outtex << "\\end{document}" << endl;
  }

void STEPWISErun::make_options(void)
  {
  //double l1 = genoptions_mult[0]->get_level1();
  //double l2 = genoptions_mult[0]->get_level2();

  //schreibt Stepwise options ins Tex-File
  outtex << "\n\\noindent {\\bf \\large Stepwise Options:}" << endl
         << "\\begin{tabbing}" << endl
         //<< "Levels for credible intervals: \\= \\\\" << endl
         //<< "Level 1: \\> " << l1 << "\\\\" << endl
         //<< "Level 2: \\> " << l2 << "\\\\" << endl
         << "Maximum number of Iterations: \\= " << steps << " \\\\" << endl
         << "Performance criterion: \\> " << criterion << " \\\\" << endl
         << "Startmodel: \\> " << startmodel << " \\\\" << endl
         //<< "Used number of Iterations: \\> blabla \\\\" << endl
         << "Increment: \\> " << increment << " \\\\" << endl;
         if(fine_tuning==true)
           outtex << "With fine\\_tuning \\\\ " << endl;
  outtex << "\\end{tabbing}\n"  << "\\vspace{0.5cm}" <<  endl;
  }

void STEPWISErun::make_predictor(void)
  {
  unsigned i;
  ST::string term = "$\\eta$ & $=$ & $\\gamma_0";
  for(i=0;i<fullcondp.size();i++)
     {
     if(fullcondp[i]->get_fctype() != MCMC::fixed
                   && fullcondp[i]->get_fctype() != MCMC::variance)
        {
        term = term + " + " + fullcondp[i]->get_term_symbolic();
        }
     else if(fullcondp[i]->get_fctype() == MCMC::fixed)
        {
        char hcharu = '_';
        ST::string hstringu = "\\_";
        unsigned z;
        for(z=1;z<fullcondp[i]->get_datanames().size();z++)
              term = term + " + " + fullcondp[i]->get_datanames()[z].insert_string_char(hcharu,hstringu);
        }
     }
  outtex << endl << "\n\\begin{tabular}{ccp{12cm}}\n" << term
         << "$\n\\end{tabular}\n\\\\ \n\\\\" << endl;
  outtex << criterion << " = " << ST::doubletostring(kriterium_tex,6) << " \\\\ \n\\\\" << endl;
  }

void STEPWISErun::make_model(void)
  {
  //Vert-Fam wird übergeben
  ST::string fam = likep_mult[0]->get_family();
  fam = fam.replaceallsigns('_', ' ');

  //Anz. Beob. wird übergeben
  unsigned obs = likep_mult[0]->get_nrobs();

  //Name der Resp.-Var. übergeben
  ST::string resp = likep_mult[0]->get_responsename();

  //schreibt das Modell und die Prior-annahmen ins Tex-File
  outtex << "\n\\noindent {\\bf \\large Response:}" << endl
         << "\\begin{tabbing}\n"
         << "Number of observations: \\= " << obs << "\\\\" << endl
         << "Response Variable: \\> " << resp << "\\\\" << endl
         << "Family: \\> " << fam << "\\\\" << endl
         << "\\end{tabbing} \n" << endl;
  }

void STEPWISErun::make_prior(vector<vector<double> > & lambdavec,
                 vector<vector<unsigned> > & startindex)
  {
  vector<ST::string> names_fixed = fullcond_alle[0]->get_datanames();
  outtex << "\n\\noindent {\\bf \\large Priors:}\\\\" << endl << "\\\\" << endl;
  unsigned i;
  if(names_fixed.size()>1)
     {
     outtex << "Fixed Effects:\\\\";
     ST::string term = "";
     for(i=1;i<names_fixed.size()-1;i++)
        term = term + "$" + names_fixed[i] + "$, ";
     term = term +  "$" + names_fixed[names_fixed.size()-1];
     outtex << endl << "\\begin{tabular}{ccp{12cm}}\n" << term
            << "$\n\\end{tabular}\n" << endl;
     if(startmodel == "full" || startmodel == "userdefined")
        outtex << "Startvalue is the fixed effect \\\\ \n\\\\" << endl;
     else if(startmodel == "empty")
        outtex << "Startvalue is \\glqq effect excluded\\grqq \\\\ \n\\\\" << endl;
     else
        outtex << "1. Startvalue is \\glqq effect excluded\\grqq \\\\" << endl
               << "2. Startvalue is the fixed effect \\\\ \n\\\\" << endl;
    }
  unsigned j;
  unsigned z = 0;
  for(i=1;i<fullcond_alle.size();i++)
     {
     if(fullcond_alle[i]->get_fctype() != MCMC::fixed
                 && fullcond_alle[i]->get_fctype() != MCMC::variance)
       {
       j = 0;
       vector<ST::string> prior = fullcond_alle[i]->get_priorassumptions();
       for(j=0;j<prior.size()-1;j++)
          outtex << prior[j] << "\\\\" << endl;
       if(fullcond_alle[i]->get_lambdamin()!=0 && fullcond_alle[i]->get_lambdamin()!=-1)
          {
          outtex << "Minimal value for the smoothing parameter: $\\lambda = "
                 << ST::doubletostring(fullcond_alle[i]->get_lambdamin())
                 << " \\,\\, \\hat{=} \\,\\, \\mbox{df} = ";
          fullcond_alle[i]->update_stepwise(fullcond_alle[i]->get_lambdamin());
          if(fullcond_alle[i]->get_df_equidist()==false && fullcond_alle[i]->get_lambdamin_opt()==false)
            outtex << ST::doubletostring(fullcond_alle[i]->compute_df(),6);
          else
            outtex << ST::doubletostring(fullcond_alle[i]->compute_df(),6)
                   << " \\approx " << ST::doubletostring(fullcond_alle[i]->get_df_lambdamin());
          outtex << "$ \\\\ \n";

          outtex << "Maximal value for the smoothing parameter: $\\lambda = "
                 << ST::doubletostring(fullcondp[i]->get_lambdamax(),6)
                 << " \\,\\, \\hat{=} \\,\\, \\mbox{df} = ";
          fullcond_alle[i]->update_stepwise(fullcond_alle[i]->get_lambdamax());
          if(fullcond_alle[i]->get_df_equidist()==false && fullcond_alle[i]->get_lambdamax_opt()==false)
            outtex << ST::doubletostring(fullcond_alle[i]->compute_df(),6);
          else
            outtex << ST::doubletostring(fullcond_alle[i]->compute_df(),6)
                   << " \\approx " << ST::doubletostring(fullcond_alle[i]->get_df_lambdamax());
          outtex << "$ \\\\ \n";

          if(fullcond_alle[i]->get_df_equidist()==true)
            outtex << "Number of different smoothing parameters with equidistant degrees of freedom: "
                   << ST::doubletostring(fullcond_alle[i]->get_number()) << " \\\\ \n";
          else
            outtex << "Number of different smoothing parameters on a logarithmic scale: "
                   << ST::doubletostring(fullcondp[i]->get_number()) << " \\\\ \n";
          }

       if(fullcond_alle[i]->get_forced()==true)
         outtex << "Without the excluded effect" << " \\\\ \n";

       j = 0;
       for(j=0;j<startindex.size();j++)
          {
          if(lambdavec[z][startindex[j][z]]==0)
             {
             if(startindex.size()>1)
                outtex << ST::doubletostring(j+1) << ". ";
             outtex << "Startvalue is \\glqq effect excluded\\grqq \\\\ \n";
             }
          else if(lambdavec[z][startindex[j][z]]==-1)
             {
             if(startindex.size()>1)
                outtex << ST::doubletostring(j+1) << ". ";
             outtex << "Startvalue is the fixed effect \\\\ \n";
             }
          else
             {
             if(startindex.size()>1)
                outtex << ST::doubletostring(j+1) << ". ";
             outtex << "Startvalue of the smoothing parameter: $\\lambda = "
                    << ST::doubletostring(lambdavec[z][startindex[j][z]],6)
                    << " \\,\\, \\hat{=} \\,\\, \\mbox{df} = ";
             fullcond_alle[i]->update_stepwise(lambdavec[z][startindex[j][z]]);
             outtex << ST::doubletostring(fullcond_alle[i]->compute_df(),6) << "$ \\\\ \n";
             }
          }
       outtex << "\\\\" << endl;
       z = z + 1;
       }
     }
  }

void STEPWISErun::make_fixed_table(void)
  {
  /*
  // falls andere Quantile gewünscht werden
  double u = fullcondp[begin[0]]->get_level1();
  double o = fullcondp[begin[0]]->get_level2();
  double u1 = fullcondp[begin[0]]->get_lower1();
  double u2 = fullcondp[begin[0]]->get_upper2();
  double o1 = fullcondp[begin[0]]->get_lower2();
  double o2 = fullcondp[begin[0]]->get_upper1();
  ST::string u_str = ST::doubletostring(u,0);
  ST::string o_str = ST::doubletostring(o,0);
  ST::string u1_str = ST::doubletostring(u1,5);
  ST::string u2_str = ST::doubletostring(u2,5);
  ST::string o1_str = ST::doubletostring(o1,5);
  ST::string o2_str = ST::doubletostring(o2,5);
  */
  vector<ST::string> h;
  unsigned j;
  unsigned r = 2;
  outtex << "\n\\newpage \n" << endl << "\n\\noindent {\\bf \\large Fixed Effects:}\\\\"
           << endl << "\\\\" << endl;
  outtex << "\\begin{tabular}{|r|r|}" << endl << "\\hline" << endl
         << "Variable & Mean \\\\" << endl << "\\hline" << endl;
  h = fullcondp[0]->get_results_latex();
  for(j=0;j<h.size();j++)
     {
     r++;
     if (r < 39)
        outtex << h[j].substr(0,h[j].length()-17) << "\\\\" << endl;
     else
        {
        r=1;
        outtex << "\\hline \n\\end{tabular}" << endl;
        outtex << "\n\\newpage \n" << endl
               << "\n\\noindent {\\bf \\large Fixed Effects (continued):}\\\\"
               << endl << "\\\\" << endl;
        outtex << "\\begin{tabular}{|r|r|}" << endl << "\\hline" << endl
               << "Variable & Mean\\\\" << endl << "\\hline" << endl;
        outtex << h[j] << endl;
        }
     }
  outtex << "\\hline \n\\end{tabular}" << endl;
  }

void STEPWISErun::make_plots(ST::string & path_batch,
                             ST::string & path_splus)      //,ST::string & path_stata)
  {

  char hcharu = '_';
  ST::string hstringu = "\\_";

  unsigned j;

  ST::string pathresult;

  bool stil = false;

  // Schleife überprüft, ob es ein fullcond-Object
  // gibt, bei dem Effekt gezeichnet werden kann
  MCMC::plotstyles plst;
  for(j=0;j<fullcondp.size();j++)
    {
    plst = fullcondp[j]->get_plotstyle();
    if(plst != MCMC::noplot)
      stil = true;
    }

  if(stil == true)
    {
    //erzeugt File, das Plot-Befehle für Java-Version enthält
    ofstream outbatch(path_batch.strtochar());

    //erzeugt File, das SPlus-Befehle zum Plotten enthält
    ofstream outsplus(path_splus.strtochar());

    outtex << "\n\\newpage" << "\n\\noindent {\\bf \\large Plots:}" << endl;
    outsplus << "# NOTE: 'directory' must be substituted by the directory"
             << " where the sfunctions are stored \n"
             << endl
    // einlesen der Source-Files für S-Plus
             << "source(\"'directory'\\\\sfunctions\\\\plotsample.s\")" << endl
             << "source(\"'directory'\\\\sfunctions\\\\plotnonp.s\")" << endl
             << "source(\"'directory'\\\\sfunctions\\\\plotsurf.s\")" << endl
             << "source(\"'directory'\\\\sfunctions\\\\drawmap.s\")" << endl
             << "source(\"'directory'\\\\sfunctions\\\\readbndfile.s\")\n" << endl;

    genoptions_mult[0]->out("  --------------------------------------------------------------------------- \n");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  Batch file for visualizing effects of nonlinear functions is stored in file \n");
    genoptions_mult[0]->out("  " + path_batch + "\n");
    genoptions_mult[0]->out("\n");

    bool stil2 = true;
    for(j=begin[0];j<=end[0];j++)  //Schleife überprüft, ob es map-Objekt gibt
      {
      plst = fullcondp[j]->get_plotstyle();
      if(plst == MCMC::drawmap)
        stil2 = false;
      }

    if(stil2 == true)
      {
      genoptions_mult[0]->out("  --------------------------------------------------------------------------- \n");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("  Batch file for visualizing effects of nonlinear functions ");
      genoptions_mult[0]->out("  in S-Plus is stored in file \n");
      genoptions_mult[0]->out("  " + path_splus + "\n");
      genoptions_mult[0]->out("\n");
      }

    if(stil2 == false)
      {
      genoptions_mult[0]->out("  NOTE: 'input filename' must be substituted by the filename of the boundary-file \n");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("  --------------------------------------------------------------------------- \n");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("  Batch file for visualizing effects of nonlinear functions ");
      genoptions_mult[0]->out("  in S-Plus is stored in file \n");
      genoptions_mult[0]->out("  " + path_splus + "\n");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("  NOTE: 'input filename' must be substituted by the filename of the boundary-file \n");
      genoptions_mult[0]->out("\n");
      }

    outbatch << "% usefile " << path_batch << endl;

    // falls andere Quantile gewünscht werden
    //double u = fullcondp[begin[0]]->get_level1();
    //double o = fullcondp[begin[0]]->get_level2();
    double u1 = fullcondp[begin[0]]->get_lower1();
    //double u2 = fullcondp[begin[0]]->get_upper2();
    //double o1 = fullcondp[begin[0]]->get_lower2();
    //double o2 = fullcondp[begin[0]]->get_upper1();
    //ST::string u_str = ST::doubletostring(u,0);
    //ST::string o_str = ST::doubletostring(o,0);
    ST::string u1_str = ST::doubletostring(u1,5);
    //ST::string u2_str = ST::doubletostring(u2,5);
    //ST::string o1_str = ST::doubletostring(o1,5);
    //ST::string o2_str = ST::doubletostring(o2,5);

    // durchlaufen der Fullconditionals
    for(j=0;j<fullcondp.size();j++)
      {
      // Pfad der Regr.-Ergebnisse
      pathresult = fullcondp[j]->get_pathresult();

      // Plotstyle: noplot, plotnonp, drawmap
      plst = fullcondp[j]->get_plotstyle();

      if (plst != MCMC::noplot)
        {
        // Pfade für ps-, tex-, SPlus-files
        ST::string pathps = pathresult.substr(0, pathresult.length()-4);
        ST::string pathgr = pathps.replaceallsigns('\\', '/');

        char hchar = '\\';
        ST::string hstring = "\\\\";

        ST::string pathps_spl = pathps.insert_string_char(hchar,hstring);
        ST::string pathres_spl = pathresult.insert_string_char(hchar,hstring);
        if(plst == MCMC::plotnonp)
           {
           outbatch << "\n";                // Befehle f. d. batch-file
           outbatch << "dataset _dat" << endl;
           outbatch << "_dat.infile using " << pathresult << endl;
           outbatch << "graph _g" << endl;
           vector<ST::string> varnames = fullcondp[j]->get_datanames();
           ST::string xvar = varnames[0];
           outbatch << "_g.plot " << xvar
                    << " pmean pqu" << u1_str.replaceallsigns('.','p')
                    // << " pqu"
                    // << o1_str.replaceallsigns('.','p') << " pqu"
                    // << o2_str.replaceallsigns('.','p') << " pqu"
                    // << u2_str.replaceallsigns('.','p')
                    << ", " << "title = \"Effect of " << xvar << "\" xlab = " << xvar
                    << " ylab = \" \" " << "outfile = " << pathps
                    << ".ps replace using _dat" << endl;
           outbatch << "drop _dat" << endl;
           outbatch << "drop _g" << endl;
            // Plot-Befehle f. d. SPlus-file
           outsplus << "plotnonp(\"" << pathres_spl << "\", psname = \""
                    << pathps_spl << ".ps\")" << endl;
            // Plot-Befehle f. d. tex-file
           outtex << "\n\\begin{figure}[h!]" << endl
                  << "\\centering" << endl
                  << "\\includegraphics[scale=0.6]{" << pathgr << ".ps}" << endl
                  << "\\caption{Non--linear Effect of '"
                  << xvar.insert_string_char(hcharu,hstringu) << "'";
           outtex << "." << endl << "Shown are the posterior means.}" << endl
                  << "\\end{figure}" << endl;
           }
          // für map-Funktionen
        else if (plst == MCMC::drawmap)
           {
           outbatch << "\n";                 // Befehle f. d. batch-file
           outbatch << "dataset _dat" << endl;
           outbatch << "_dat.infile using " << pathresult << endl;
           outbatch << "map _map" << endl;
           outbatch << "_map.infile using input_filename" << endl;
           outbatch << "graph _g" << endl;
           vector<ST::string> varnames = fullcondp[j]->get_datanames();
           ST::string regionvar = varnames[0];
           outbatch << "_g.drawmap " << "pmean" << " " << regionvar
                    << ", map = _map color outfile = " << pathps
                    << "_pmean.ps replace using _dat" << endl;
           /*
           outbatch << "_g.drawmap " << "pcat" << u_str << " " << regionvar
                      << ", map = _map nolegend pcat outfile = " << pathps
                      << "_pcat" << u_str << ".ps replace using _dat" << endl;
           outbatch << "_g.drawmap " << "pcat" << o_str << " " << regionvar
                      << ", map = _map nolegend pcat outfile = " << pathps
                      << "_pcat" << o_str << ".ps replace using _dat" << endl;
           */
           outbatch << "drop _dat" << endl;
           outbatch << "drop _g" << endl;
           outbatch << "drop _map" << endl;
            // Plot-Befehle f. d. SPlus-file
           outsplus << "# NOTE: 'input_filename' must be substituted by the "
                    << "filename of the boundary-file \n"
                    << "# NOTE: choose a 'name' for the map \n" << endl
                    << "readbndfile(\"'input_filename'\", \"'name'\")" << endl
                    << "drawmap(map = 'name', outfile = \"" << pathps_spl
                    <<"_pmean.ps\", dfile = \"" << pathres_spl
                    << "\" ,plotvar = \"pmean\", regionvar = \""
                    << regionvar << "\", color = T)" << endl;
           /*
           outsplus << "drawmap(map = 'name', outfile = \"" << pathps_spl
                    <<"_pcat" << u_str << ".ps\", dfile = \"" << pathres_spl
                    << "\" ,plotvar = \"pcat" << u_str << "\", regionvar = \""
                    << regionvar << "\", legend = F, pcat = T)" << endl;
           outsplus << "drawmap(map = 'name', outfile = \"" << pathps_spl
                    <<"_pcat" << o_str << ".ps\", dfile = \"" << pathres_spl
                    << "\",plotvar = \"pcat" << o_str << "\", regionvar = \""
                    << regionvar << "\", legend = F, pcat = T)" << endl;
           */
            // Plot-Befehle f. d. tex-file
           outtex << "\n\\begin{figure}[h!]" << endl
                  << "\\centering" << endl
                  << "\\includegraphics[scale=0.6]{" << pathgr << "_pmean.ps}"
                  << endl
                  << "\\caption{Non--linear Effect of '" <<
                  regionvar.insert_string_char(hcharu,hstringu) << "'";
           outtex << ". Shown are the posterior means.}" << endl
                  << "\\end{figure}" << endl;
           /*
           outtex << "\n\\begin{figure}[htb]" << endl
                  << "\\centering" << endl
                  << "\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                  << u_str << ".ps}" << endl
                  << "\\caption{Non--linear Effect of '" << regionvar << "'";
           outtex << ". Posterior probabilities for a nominal level of "
                  << u_str << "\\%." << endl
                  << "Black denotes regions with strictly negative credible intervals,"
                  << endl
                  << "white denotes regions with strictly positive credible intervals.}"
                  << endl << "\\end{figure}" << endl;
           outtex << "\n\\begin{figure}[htb]" << endl
                  << "\\centering" << endl
                  << "\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                  << o_str << ".ps}" << endl
                  << "\\caption{Non--linear Effect of '" << regionvar << "'";
           outtex << ". Posterior probabilities for a nominal level of "
                  << o_str << "\\%." << endl
                  << "Black denotes regions with strictly negative credible intervals,"
                  << endl
                  << "white denotes regions with strictly positive credible intervals.}"
                  << endl << "\\end{figure}" << endl;
           */

           } // end: else if (plst == MCMC::drawmap && fullcondp[j]->get_col()==i)
        } // end: if (plst != MCMC::noplot)
      } // end: for(j=begin[0];j<=end[nr];j++)
    }
  }


}


