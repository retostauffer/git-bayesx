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


// -----------------------------------------------------------------------------
// ------- die grundlegenden Funktionen ----------------------------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::stepwise(const ST::string & procedure, const ST::string & crit,
        const int & stp, const ST::string & trac, const int & number,
        const ST::string & stam, const int & inc, const bool & finet,
        const datamatrix & Da, const vector<ST::string> & modelvar,
        const ST::string & name, vector<FULLCOND*> & fullcond_z, ST::string & path)
  {

  D = Da;
  modelv = modelvar;
  algorithm = procedure;
  criterion = crit;
  increment = inc;
  steps = stp;
  startmodel = stam;
  fine_tuning = finet;
  trace = trac;

  //vector<double> modell_alt;
  //double kriterium_alt;
  //ST::string text_alt;
  modell_alt.erase(modell_alt.begin(),modell_alt.end());
  vector<double> modell_final;
  double kriterium_final;
  ST::string text_final;

  ST::string tr_akt = "trace_on";

  lambdavec.erase(lambdavec.begin(),lambdavec.end());
  names_fixed.erase(names_fixed.begin(),names_fixed.end());
  names_nonp.erase(names_nonp.begin(),names_nonp.end());
  //vector<vector<double> > lambdavec;
  //vector<ST::string> names_fixed;
  //vector<vector<ST::string> > names_nonparametric;

  fullcond_alle = fullcondp;

  bool gewichte = false;
  if(likep_mult[0]->get_family() != "Gaussian")
    gewichte = true;
  initialise_lambdas(names_nonp,names_fixed,lambdavec,number,gewichte);

  unsigned i;
  
  // Fehlerabfragen
  unsigned j;
  for(j=0;j<fullcond_alle.size();j++)
     {
     if(fullcond_alle[j]->geterrors().size()>0)
      {
      for(i=0;i<fullcond_alle[j]->geterrors().size();i++)
         genoptions_mult[0]->out(fullcond_alle[j]->geterrors()[i]);
      return true;
      }
    }
  // überprüfen, dass Randomslopes nicht auch als fixe Effekte angegeben werden!
  if( vcm_doppelt(names_fixed,names_nonp) == true)
     return true;

  vector<vector<unsigned> > startindex;
  vector<vector<double> > startfix;
  startwerte(startmodel,lambdavec,names_fixed,startindex,startfix);

  options_text(number,lambdavec,names_fixed,names_nonp,startfix,startindex,name);

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
         lambdavec,startindex[i],startfix[i],names_fixed,
         names_nonp,D,modelv,true);

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
  fix_komplett(modell_final,names_fixed,names_nonp,D,modelv);
  fullcond_komplett(modell_final,names_fixed.size(),names_nonp,D,modelv);
  tr_akt = "trace_on";
  maketext(header,modell_final,kriterium_final,text_final,false,tr_akt,false);
  genoptions_mult[0]->out("\n\n");
  kriterium_tex = kriterium_final;

  if(fine_tuning == true)
     abbruch = finetuning(modell_final,lambdavec,names_fixed,D,modelv);

  if(abbruch==true)
    return true;

  fullcond_z = fullcondp;
  for(i=0;i<fullcond_z.size();i++)
     fullcond_z[i]->set_fcnumber(i);

  //likep_mult[0]->set_linpred_null();

  posteriormode(title,false);  // Problem: linearer Prädiktor bei "true" standardisiert! Hier wird zurückgerechnet!
                               // danach nicht mehr compute_criterion() aufrufen!!!
  make_tex_end(path,modell_final,names_fixed);

  // Files müssen wieder geschlossen werden!!!
  outtex.close();
  outcriterium.close();
  outmodels.close();


  // gibt Lambdas aus, damit man die richtig bestimmten Variablen zählen kann!
  ST::string zaehlername = path + "_lambdas_" + likep_mult[0]->get_responsename() + ".ascii";
  //zaehlername = zaehlername + "_" + ST::inttostring(increment) + ".ascii";
  ofstream out(zaehlername.strtochar());
  ST::string beschriftung = "  krit   ";
  ST::string eintrag = "  " + ST::doubletostring(kriterium_final) + "   ";
  for(i=1;i<names_fixed.size();i++)
     beschriftung = beschriftung + names_fixed[i] + "   ";
  for(i=0;i<names_nonp.size();i++)
     beschriftung = beschriftung + names_nonp[i][0] + "   ";
  for(i=0;i<modell_final.size();i++)
     eintrag = eintrag + ST::doubletostring(modell_final[i]) + "   ";
  out << beschriftung << endl;
  out << eintrag << endl;
  
  
  return false;
  }


bool STEPWISErun::single_stepwise(double & kriterium_alt,
         vector<double> & modell_alt, ST::string & text_alt,
         const vector<vector<double> > & lambdavec,
         const vector<unsigned> & start, const vector<double> & startfix,
         const vector<ST::string> & names_fixed,
         const vector<vector<ST::string> > & names_nonp, const datamatrix & D,
         const vector<ST::string> & modelv, const bool & tex)
  {
  modell_neu.erase(modell_neu.begin(),modell_neu.end());
  modellematrix.erase(modellematrix.begin(),modellematrix.end());
  //vector<vector<vector<double> > > modellematrix;
  //vector<double> modell_neu;
  //vector<vector<double> > kriteriummatrix;       UNNÖTIG
  //double kriterium_neu;
  int steps_aktuell = 0;
  ST::string tr_akt = "trace_on";
  vector<vector<double> > startiteration;
  //vector<double> kriteriumiteration;            UNNÖTIG

  unsigned i;
  for(i=0;i<names_fixed.size()-1;i++)
     modell_neu.push_back(startfix[i]);
  for(i=1;i<fullcond_alle.size();i++)
     {
     double lambda = lambdavec[i-1][start[i-1]];
     modell_neu.push_back(lambda);
     }
  modell_alt = modell_neu;
  startiteration.push_back(modell_alt);
  modellematrix.push_back(startiteration);
  vector<ST::string> title;
  title.push_back("");
  fix_komplett(modell_alt,names_fixed,names_nonp,D,modelv);          
  fullcond_komplett(modell_alt,names_fixed.size(),names_nonp,D,modelv);
  posteriormode(title,true);

  kriterium_alt = compute_criterion();

  /*
  // Versuch!!

  for(i=0;i<lambdavec[0].size()-2;i++)
    {
    if(lambdavec[0][i]!=modell_alt[names_fixed.size()-2])
      {
      fullcond_alle[1]->update_stepwise(lambdavec[0][i]);
      fullcond_alle[1]->posteriormode();
      double kriterium_versuch = compute_criterion();
      genoptions_mult[0]->out("AIC(" + ST::doubletostring(fullcond_alle[1]->compute_df(),6) + ") = " + ST::doubletostring(kriterium_versuch) + "\n");
      }
    }
  for(i=lambdavec[0].size()-2;i<lambdavec[0].size();i++)
    {
    if(lambdavec[0][i]!=modell_alt[names_fixed.size()-2])
      {
      fullcond_alle[1]->reset_effect(0);
      double kriterium_versuch;
      if(i==lambdavec[0].size()-2)
        {
        fullcondp[0]->include_effect(names_nonp[0],
                                  fullcond_alle[1]->get_data_forfixedeffects());
        fullcond_alle[0]->posteriormode();
        kriterium_versuch = compute_criterion()-2*(fullcond_alle[1]->compute_df()) + 2;
        reset_fix(names_nonp[0][0]);
        }
      else
        {
        fullcond_alle[0]->posteriormode();
        kriterium_versuch = compute_criterion()-2*(fullcond_alle[1]->compute_df());
        }
      genoptions_mult[0]->out("AIC(" + ST::doubletostring(fullcond_alle[1]->compute_df(),6) + ") = " + ST::doubletostring(kriterium_versuch) + "\n");
      }
    }

  // Versuch ENDE
  */
  
    if(tex==true)
     {
     kriterium_tex = kriterium_alt;
     make_predictor();
     }

  //kriteriumiteration.push_back(kriterium_alt);
  //kriteriummatrix.push_back(kriteriumiteration);
  // ENDE: startmodel

  kriterium_neu = kriterium_alt;
  outcriterium << steps_aktuell << "   " << kriterium_neu << endl;
  outmodels << steps_aktuell << "   " << kriterium_neu << "   ";
  ST::string header;
  fertig = false;    // überprüft, ob es noch nicht gerechnete Modelle gibt
  ST::string text_neu;
  bool eins = true;

       // Schleife für Minimierung
  while(kriterium_neu<=kriterium_alt && fertig==false && steps_aktuell<steps)
       {
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
       maketext(header, modell_alt,kriterium_alt,text_neu,neutext,
                tr_akt,neutext);
       text_alt = text_neu;

       vector<vector<double> > modeliteration;
       double kriterium;
       vector<double> kriteriumiteration2;

       for(i=1;i<names_fixed.size();i++)
          {
          modell_neu = modell_alt;
          if(modell_alt[i-1]==-1)
            modell_neu[i-1]= 0;
          else if(modell_alt[i-1]==0)
            modell_neu[i-1] = -1;
          if(modelcomparison(modell_neu,modellematrix)==false)
            newmodel_fix(fertig,modell_neu[i-1],kriteriumiteration2,
                 modeliteration,modell_neu,textiteration,names_fixed[i],
                 D,modelv);
          }

       unsigned z = 1;
       while(z<fullcond_alle.size() && fullcond_alle[z]->get_fctype()==factor)
          {
          modell_neu = modell_alt;
          if(modell_alt[z+names_fixed.size()-2]==-1)
            modell_neu[z+names_fixed.size()-2]= 0;
          else if(modell_alt[z+names_fixed.size()-2]==0)
            modell_neu[z+names_fixed.size()-2] = -1;
          if(modelcomparison(modell_neu,modellematrix)==false)
            newmodel_factor(fertig,modell_neu[z+names_fixed.size()-2],z,kriteriumiteration2,
                 modeliteration,modell_neu,textiteration,names_nonp[z-1],
                 D,modelv);
          z = z + 1;
          }

       if(algorithm == "stepwise")
         step_nonpfkt(kriteriumiteration2,modeliteration,textiteration,z);
       else // if(algorithm == "stepmin")
         stepmin(kriteriumiteration2,modeliteration,textiteration,z,kriterium_alt);

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
         //kriteriummatrix.push_back(kriteriumiteration2);

         header = "\n\nBest Model of this iteration:";
         fix_komplett(modell_neu,names_fixed,names_nonp,D,modelv);
         fullcond_komplett(modell_neu,names_fixed.size(),names_nonp,D,modelv);
         neutext = false;
         maketext(header,modell_neu,kriterium_neu,text_neu,neutext,
                  trace,true);                    // Modell rausschreiben!!

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

       if(make_pause() == true)
          return true;
       }

  header = "  Final Model:";
  tr_akt = "trace_on";
  maketext(header,modell_alt,kriterium_alt,text_alt,false,tr_akt,false);
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("  Used number of iterations: " + ST::inttostring(steps_aktuell));
  genoptions_mult[0]->out("\n\n");
  genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
  genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");

  return false;
  }

  
// -----------------------------------------------------------------------------
// ------- Funktionen, die sich bei Stepwise / Stepmin unterscheiden -----------
// -----------------------------------------------------------------------------

void STEPWISErun::step_nonpfkt(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)         //!!!
  {
  unsigned i;
  for(i=z;i<fullcond_alle.size();i++)
    {
    unsigned sch;
    for(sch=1;sch<=increment;sch++)
       {
       modell_neu = modell_alt;
       bool lambda_exist;    // zum Überprüfen, ob neues Lambda im Vektor enthalten ist
       unsigned index = search_lambdaindex(modell_alt[names_fixed.size()-2+i],
                                lambdavec[i-1],lambda_exist);
       lambda_exist = false;
       if(index < lambdavec[i-1].size()-sch)
          lambda_exist = true;
       if(lambda_exist==true)
          {
          modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][index+sch];
          if(modelcomparison(modell_neu,modellematrix)==false)
          newmodel_nonp(fertig,i,modell_neu,modell_alt,kriteriumiteration2,
                       modeliteration,textiteration,names_fixed,names_nonp,D,modelv);
          }

       lambda_exist = false;
       modell_neu = modell_alt;
       if(index >= sch)
          lambda_exist = true;
       if(lambda_exist==true)
          {
          modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][index-sch];
          if(modelcomparison(modell_neu,modellematrix)==false)
          newmodel_nonp(fertig,i,modell_neu,modell_alt,kriteriumiteration2,
                       modeliteration,textiteration,names_fixed,names_nonp,D,modelv);
          }
       }
    }
  }

// Fehler:
// 1) fixe Effekte -> Absturz
// 2) Es kann sein, daß das Startmodell nicht wieder korrekt hergestellt wird!  

void STEPWISErun::stepmin(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z, double & kriterium)        
  {
  // Stellt linearen Prädiktor auf das Startmodell ein. Besser wäre, den lin. Prädiktor zu speichern!!!
  fullcondp[0]->set_effect_zero();
  vector<ST::string> title;
  title.push_back("");
  posteriormode(title,true);

  unsigned i;
  for(i=z;i<fullcond_alle.size();i++)
    {
    modell_neu = modell_alt;

//unsigned y;

    vector<double> krit_fkt;
    if(modell_alt[i+names_fixed.size()-2]==0)
      stepmin_leer(i,krit_fkt,kriterium);
    else if(modell_alt[i+names_fixed.size()-2]==-1)
      stepmin_fix(i,krit_fkt,kriterium);
    else
      stepmin_nonp(i,krit_fkt,kriterium);


//for(y=0;y<krit_fkt.size();y++)
//genoptions_mult[0]->out("AIC(" + ST::doubletostring(krit_fkt.size()-1-y) + ") = " + ST::doubletostring(krit_fkt[y]) + "\n");

    double kriterium_min = krit_fkt[0];
    unsigned j;
    unsigned lambda_ind;
    for(j=1;j<krit_fkt.size();j++)  //berechnet den besten Wert
       {
       if(krit_fkt[j]<=kriterium_min)
         {
         kriterium_min = krit_fkt[j];
         lambda_ind = j;
         }
       }

    modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][lambda_ind];
    if(modell_neu[names_fixed.size()-2+i] != modell_alt[names_fixed.size()-2+i])
      {
      if(modelcomparison(modell_neu,modellematrix)==false)
        {
        newmodel_nonp(fertig,i,modell_neu,modell_alt,kriteriumiteration2,
               modeliteration,textiteration,names_fixed,names_nonp,D,modelv);

        // Stellt linearen Prädiktor wieder her. Besser wäre, den lin. Prädiktor zu speichern!!!
        fullcondp[0]->set_effect_zero();
        vector<ST::string> title;
        title.push_back("");
        posteriormode(title,true);
        }
      }
    }
 }

  
// -----------------------------------------------------------------------------
// ------------------ Funktionen für Stepmin -----------------------------------
// -----------------------------------------------------------------------------

void STEPWISErun::stepmin_nonp(unsigned & z, vector<double> & krit_fkt,
                                double & kriterium)
  {
  unsigned i;

  double df = 0;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();
  df = df - fullcond_alle[z]->compute_df();

  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=modell_alt[z+names_fixed.size()-2])
      {
      double kriterium_versuch;
      if(lambdavec[z-1][i]!=-1 && lambdavec[z-1][i]!=0)
        {
        fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
        fullcond_alle[z]->posteriormode();
        kriterium_versuch = criterion_min(df + fullcond_alle[z]->compute_df());
        }
      else if(lambdavec[z-1][i]==-1)
        {
        fullcond_alle[z]->reset_effect(0);
        fullcondp[0]->include_effect(names_nonp[z-1],
                                  fullcond_alle[z]->get_data_forfixedeffects());
        fullcond_alle[0]->posteriormode();
        kriterium_versuch = criterion_min(df + 1);
        reset_fix(names_nonp[z-1][0]);
        fullcond_alle[0]->posteriormode();
        }
      else
        {
        fullcond_alle[z]->reset_effect(0);
        kriterium_versuch = criterion_min(df);
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }
  double help = modell_alt[z+names_fixed.size()-2];
  fullcond_alle[z]->update_stepwise(modell_alt[z+names_fixed.size()-2]);
  fullcond_alle[z]->posteriormode();
  }


void STEPWISErun::stepmin_fix(unsigned & z, vector<double> & krit_fkt,
                                double & kriterium)
  {
  unsigned i;
  double df = 0;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();
  df = df - 1;

  reset_fix(names_nonp[z-1][0]);
  fullcond_alle[0]->posteriormode();
  fullcondp.push_back(fullcond_alle[z]);

  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=-1)
      {
      double kriterium_versuch;
      if(lambdavec[z-1][i]!=0)
        {
        fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
        fullcond_alle[z]->posteriormode();
        kriterium_versuch = criterion_min(df + fullcond_alle[z]->compute_df());
        }
      else
        {
        fullcond_alle[z]->reset_effect(0);
        kriterium_versuch = criterion_min(df);
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }
  fullcond_alle[z]->reset_effect(0);
  fullcondp[0]->include_effect(names_nonp[z-1],
                                  fullcond_alle[z]->get_data_forfixedeffects());
  fullcondp.erase(fullcondp.end()-1,fullcondp.end());
  fullcond_alle[0]->posteriormode();
  }


void STEPWISErun::stepmin_leer(unsigned & z, vector<double> & krit_fkt,
                                double & kriterium)
  {
  unsigned i;

  double df = 0;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();

  fullcondp.push_back(fullcond_alle[z]);

  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=0)
      {
      double kriterium_versuch;
      if(lambdavec[z-1][i]!=-1)
        {
        fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
        fullcond_alle[z]->posteriormode();
        kriterium_versuch = criterion_min(df + fullcond_alle[z]->compute_df());
        }
      else
        {
        fullcond_alle[z]->reset_effect(0);
        fullcondp[0]->include_effect(names_nonp[z-1],
                                  fullcond_alle[z]->get_data_forfixedeffects());
        fullcond_alle[0]->posteriormode();
        kriterium_versuch = criterion_min(df + 1);
        reset_fix(names_nonp[z-1][0]);
        fullcond_alle[0]->posteriormode();
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }
  fullcond_alle[z]->reset_effect(0);
  fullcondp.erase(fullcondp.end()-1,fullcondp.end());
  }


double STEPWISErun::criterion_min(double & df)
  {
  double kriterium;

  if(criterion=="GCV")
    kriterium = likep_mult[0]->compute_gcv(df);
  else if(criterion=="AIC")
    kriterium = likep_mult[0]->compute_aic(df);
  else if(criterion=="BIC")
    kriterium = likep_mult[0]->compute_bic(df);
  else //if(criterion=="AIC_imp")
    kriterium = likep_mult[0]->compute_improvedaic(df);
  //else if(criterion=="MSEP")
  //  kriterium = likep_mult[0]->compute_msep();
  //else //if(criterion=="AUC")
  //  kriterium = -1 * likep_mult[0]->compute_auc();

  return kriterium;
  }


// -----------------------------------------------------------------------------
// --------------------------- Fine-Tuning -------------------------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::finetuning(vector<double> & modell, vector<vector<double> > & lambdavec,
                 vector<ST::string> & names_fixed, const datamatrix & D,
                 const vector<ST::string> & modelv)
  {
  genoptions_mult[0]->out("  BEGIN FINE-TUNING:");
  genoptions_mult[0]->out("\n\n");
  vector<FULLCOND*> fullcond_neu;
  double number = 11;
  unsigned i = 0;
  unsigned anzahl = 0;
  vector<double> lambdamaxe;

  fullcond_neu.push_back(fullcondp[0]);
  for(i=1;i<fullcond_alle.size();i++)
     {
     if(modell[names_fixed.size()-2+i]!=0)
           {
           if(modell[names_fixed.size()-2+i]!=-1)
              anzahl = anzahl + 1;
           fullcond_neu.push_back(fullcond_alle[i]);

           bool lambda_exist;
           unsigned index = search_lambdaindex(modell[names_fixed.size()-2+i],lambdavec[i-1],lambda_exist);
           double lambdastart = lambdavec[i-1][index];
           double lambdamin;
           double lambdamax;
           if(fullcond_alle[i]->get_fctype() != MCMC::factor)
              {
              lambda_exist = false;
              int j = 2;
              while(lambda_exist==false && j>0)
                 {
                 if(index>=j)
                    {
                    if(lambdavec[i-1][index-j]!=0 && lambdavec[i-1][index-j]!=-1)
                       lambda_exist = true;
                    }
                 j = j - 1;
		         }
              if(lambda_exist==true)
                lambdamin = lambdavec[i-1][index-j-1];
              else
                lambdamin = fullcond_alle[i]->get_lambdamin()/10;

              lambda_exist = false;
              j = 2;
              while(lambda_exist==false && j>0)
                 {
		         if(index<lambdavec[i-1].size()-j)
		            {
		            if(lambdavec[i-1][index+j]!=0 && lambdavec[i-1][index+j]!=-1)
                    lambda_exist = true;
		            }
		        j = j - 1;
 		        }
              if(lambda_exist==true)
                lambdamax = lambdavec[i-1][index+j+1];
              else if(lambdavec[i-1][index]==-1)
                lambdamax = 1.5*lambdavec[i-1][index-1] - 0.5*lambdavec[i-1][index-2];
              else
                lambdamax = 1.5*lambdavec[i-1][index] - 0.5*lambdavec[i-1][index-1];

              lambdamaxe.push_back(lambdamax);
              if(lambdastart!=-1)
                fullcond_alle[i]->set_stepwise_options(lambdastart,lambdastart,lambdamin,
                                fullcond_alle[i]->get_forced(),0,0,false,false,number,false);
              else
                fullcond_alle[i]->set_stepwise_options(lambdastart,lambdamax,lambdamin,
                                fullcond_alle[i]->get_forced(),0,0,false,false,21,false);
              }
           }
     }

  fullcondp = fullcond_neu;
  fullcond_alle = fullcondp; // hier anders, wenn es für als fix ausgew. Var. wieder mehr Mögl. geben soll

  vector<vector<ST::string> > names_nonp;
  vector<vector<double> > lambdavec_ursp;
  vector<vector<double> > lambdavec_fine;
  vector<ST::string> namesfixed_neu;
  initialise_lambdas(names_nonp,namesfixed_neu,lambdavec_ursp,number,false);

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

  unsigned w = 0;
  for(i=1;i<fullcond_alle.size();i++)      // hier wird die zweite Hälfte (Wert > Startwert) der Lambdas
     {                                     // bestimmt (damit gleich viele Lambdas links und rechts vom Startwert)
     if(fullcond_alle[i]->get_fctype() != MCMC::factor)
          {
          if(fullcond_alle[i]->get_lambdastart()!=-1)
            {
            fullcond_alle[i]->set_stepwise_options(fullcond_alle[i]->get_lambdastart(),
                  lambdamaxe[w],fullcond_alle[i]->get_lambdastart(),
                  fullcond_alle[i]->get_forced(),0,0,false,false,number,false);
            vector<double> untervector;
            unsigned j = 0;
            for(j=0;j<lambdavec_ursp[i-1].size()-3;j++)
              untervector.push_back(lambdavec_ursp[i-1][j]);
            if(fullcond_alle[i]->get_forced()==true)
              untervector.push_back(lambdavec_ursp[i-1][lambdavec_ursp[i-1].size()-3]);
            fullcond_alle[i]->compute_lambdavec(untervector,number);
            lambdavec_fine.push_back(untervector);
            }
          else
            lambdavec_fine.push_back(lambdavec_ursp[i-1]);
          w = w + 1;
          }
        else
          lambdavec_fine.push_back(lambdavec_ursp[i-1]);
     }

  vector<vector<unsigned> > startindex;
  vector<vector<double> > startfix;
  startwerte("userdefined",lambdavec_fine,names_fixed,startindex,startfix);

  ST::string text;
  double kriterium;

  bool abbruch = false;

  abbruch = single_stepwise(kriterium,modell,text,lambdavec_fine,
              startindex[0],startfix[0],names_fixed,names_nonp,D,
              modelv,false);

  if(abbruch==true)
    return true;

  ST::string header = "  Final Model after Fine-Tuning:";
  fix_komplett(modell,names_fixed,names_nonp,D,modelv);
  fullcond_komplett(modell,names_fixed.size(),names_nonp,D,modelv);
  ST::string tr_akt = "trace_on";
  maketext(header,modell,kriterium,text,false,tr_akt,false);
  kriterium_tex = kriterium;
  genoptions_mult[0]->out("\n\n");

  return false;
  }


// -----------------------------------------------------------------------------
// ------- Funktionen für die Erstellung des Startmodels -----------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::vcm_doppelt(const vector<ST::string> & names_fixed,
      const vector<vector<ST::string> > & names_nonp)
  {
  unsigned i;
  unsigned j;
  unsigned k;
  unsigned zaehler = 0;
  bool fehler_randomslope = false;
  for(i=1;i<fullcond_alle.size();i++)
     {
     if(fullcond_alle[i]->get_fctype() == MCMC::randomslopes)
       {
       for(j=1;j<names_fixed.size();j++)
          {
          if(names_fixed[j]==names_nonp[zaehler][0])
             {
             k = j;
             fehler_randomslope = true;
             }
          }
       }
     zaehler = zaehler + 1;
     }
  if(fehler_randomslope==true)                
    {
    genoptions_mult[0]->outerror("\n\n ERROR: You must not put fixed effect " +
       names_fixed[k] + " in the model! \n");
    return true;
    }
  else
    return false;
  }


void STEPWISErun::initialise_lambdas(vector<vector<ST::string> > & names_nonp,
                 vector<ST::string> & names_fixed, vector<vector<double> > & lambdavec,
                 const int & number, const bool & gewichte)
  {
  names_fixed = fullcondp[0]->get_datanames();
  unsigned i;

       // berechnet ein Modell, um Gewichte für Augabe usw. zu erhalten!!!
  if(gewichte == true)
    {
    vector<ST::string> title;
    title.push_back("");
    vector<double> modell_init;
    for(i=1;i<names_fixed.size();i++)
       modell_init.push_back(-1);
    for(i=1;i<fullcond_alle.size();i++)
      {
      if(fullcond_alle[i]->get_fctype() != MCMC::factor)
         fullcond_alle[i]->update_stepwise(100);
      else if(fullcond_alle[i]->get_fctype() == MCMC::factor)
         {
         fullcondp.erase(fullcondp.begin() + 1);   //löscht das Element an Pos.1 aus fullcondp-Vektor
         fullcondp[0]->include_effect(fullcond_alle[i]->get_datanames(),
                                   fullcond_alle[i]->get_data_forfixedeffects());
         }
      }
    end[0] = fullcondp.size()-1;

    posteriormode(title,true);
    }

  for(i=1;i<fullcond_alle.size();i++)
     {
     int nummer = number;
     if(fullcond_alle[i]->get_data_forfixedeffects().cols() > 1)  //bei Faktor-Variablen
        names_nonp.push_back(fullcond_alle[i]->get_datanames());
     else
        {
        vector<ST::string> names_help;
        names_help.push_back(fullcond_alle[i]->get_datanames()[fullcond_alle[i]->get_datanames().size()-1]);
        names_nonp.push_back(names_help);
        nummer = fullcond_alle[i]->get_number();
        if(nummer==0)
           nummer = number;
        fullcond_alle[i]->set_stepwise_options(fullcond_alle[i]->get_lambdastart(),
                       fullcond_alle[i]->get_lambdamax(),fullcond_alle[i]->get_lambdamin(),
                       fullcond_alle[i]->get_forced(),fullcond_alle[i]->get_df_lambdamax(),
                       fullcond_alle[i]->get_df_lambdamin(),fullcond_alle[i]->get_lambdamax_opt(),
                       fullcond_alle[i]->get_lambdamin_opt(),nummer,
                       fullcond_alle[i]->get_df_equidist());
        }
     vector<double> untervector;
     fullcond_alle[i]->compute_lambdavec(untervector,nummer);
     lambdavec.push_back(untervector);
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
  if(startmodel == "empty" || startmodel == "both")
    {
    vector<unsigned> indexhelp;
    vector<double> fixhelp;

    for(i=1;i<names_fixed.size();i++)
       fixhelp.push_back(0);

    for(i=1;i<fullcond_alle.size();i++)
       indexhelp.push_back(lambdavec[i-1].size()-1);

    startindex.push_back(indexhelp);
    startfix.push_back(fixhelp);
    }
  if(startmodel == "both" || startmodel == "full")
    {
    vector<unsigned> indexhelp;
    vector<double> fixhelp;

    for(i=1;i<names_fixed.size();i++)
       fixhelp.push_back(-1);
    for(i=1;i<fullcond_alle.size();i++)
       indexhelp.push_back(0);

    startindex.push_back(indexhelp);
    startfix.push_back(fixhelp);
    }
  else if(startmodel == "userdefined")
    {
    vector<unsigned> indexhelp;
    vector<double> fixhelp;

    for(i=1;i<names_fixed.size();i++)
       fixhelp.push_back(-1);                       //fehlt: vom Benutzer angeben lassen!!!
    for(i=1;i<fullcond_alle.size();i++)
       {
       double start = fullcondp[i]->get_lambdastart();
       unsigned index = search_lambdastartindex(start,lambdavec[i-1]);
       indexhelp.push_back(index);
       }

    startindex.push_back(indexhelp);
    startfix.push_back(fixhelp);
    }
  }


// -----------------------------------------------------------------------------
// ------- Funktionen für die Berechnung neuer Modelle -------------------------
// -----------------------------------------------------------------------------

double STEPWISErun::compute_criterion(void)
  {
  double df = 0;
  if(criterion != "MSEP" && criterion != "AUC")
    {
    unsigned i;
    for(i=0;i<fullcondp.size();i++)
      df = df + fullcondp[i]->compute_df();
    }
  double kriterium;

  if(criterion=="GCV")
    kriterium = likep_mult[0]->compute_gcv(df);
  else if(criterion=="AIC")
    kriterium = likep_mult[0]->compute_aic(df);
  else if(criterion=="BIC")
    kriterium = likep_mult[0]->compute_bic(df);
  else if(criterion=="AIC_imp")
    kriterium = likep_mult[0]->compute_improvedaic(df);
  else if(criterion=="MSEP")
    kriterium = likep_mult[0]->compute_msep();
  else //if(criterion=="AUC")
    kriterium = -1 * likep_mult[0]->compute_auc();

  return kriterium;
  }


void STEPWISErun::newmodel(bool & fertig, const vector<double> & modell,
    vector<double> & krit, vector<vector<double> > & mi,
    vector<ST::string> & textit)
  {
  fertig = false;
  mi.push_back(modell);
  vector<ST::string> title;
  title.push_back("");
  posteriormode(title,true);      //keine Ausgabe
  double kriterium = compute_criterion();
  ST::string header = "  Trial: ";
  ST::string text;
  maketext(header,modell,kriterium,text,true,trace,false);
  textit.push_back(text);
  krit.push_back(kriterium);
  }


void STEPWISErun::newmodel_fix(bool & fertig, const double & mo,
    vector<double> & krit, vector<vector<double> > & mi,
    const vector<double> & modell, vector<ST::string> & textit,
    const ST::string & name, const datamatrix & D, const vector<ST::string> & modelv)
  {
  if(mo==0)
    reset_fix(name);
  else
    include_fix(name,D,modelv);
  fullcondp[0]->set_effect_zero();
  newmodel(fertig,modell,krit,mi,textit);
  if(mo==0)
    include_fix(name,D,modelv);
  else
    reset_fix(name);
  //fullcondp[0]->set_effect_zero();
  }


void STEPWISErun::newmodel_factor(bool & fertig, const double & mo, const unsigned & index,
    vector<double> & krit, vector<vector<double> > & mi,
    const vector<double> & modell, vector<ST::string> & textit,
    const vector<ST::string> & name, const datamatrix & D, const vector<ST::string> & modelv)
  {
  unsigned i;
  if(mo==0)
    {
    for(i=0;i<name.size();i++)
      reset_fix(name[i]);
    }
  else
    fullcondp[0]->include_effect(name,fullcond_alle[index]->get_data_forfixedeffects());
  fullcondp[0]->set_effect_zero();
  newmodel(fertig,modell,krit,mi,textit);
  if(mo==0)
    fullcondp[0]->include_effect(name,fullcond_alle[index]->get_data_forfixedeffects());
  else
    {
    for(i=0;i<name.size();i++)
      reset_fix(name[i]);
    }
  //fullcondp[0]->set_effect_zero();
  }


void STEPWISErun::newmodel_nonp(bool & fertig, const unsigned & index,
    const vector<double> & modell, const vector<double> & modell_alt,
    vector<double> & krit, vector<vector<double> > & mi,
    vector<ST::string> & textit, const vector<ST::string> & names_fixed,
    const vector<vector<ST::string> > & names_nonp, const datamatrix & D,
    const vector<ST::string> & modelv)
  {
  fullcond_einzeln(modell,modell_alt,index,names_fixed.size(),names_nonp,D,modelv);
  fullcondp[0]->set_effect_zero();
  newmodel(fertig,modell,krit,mi,textit);
  fullcond_einzeln(modell_alt,modell,index,names_fixed.size(),names_nonp,D,modelv);
  //fullcondp[0]->set_effect_zero();
  }


bool STEPWISErun::modelcomparison(const vector<double> & m,
                   const vector<vector<vector<double> > > & mmatrix)
  {
  bool s = false;
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
  return s;
  }


// -----------------------------------------------------------------------------
// ------- Funktionen für die Erstellung des fullcondp-Vektors -----------------
// -----------------------------------------------------------------------------

void STEPWISErun::fullcond_einzeln(const vector<double> & modell_neu,
         const vector<double> & modell_alt, const unsigned & index,
         const unsigned & nf_size, const vector<vector<ST::string> > & names_nonp,
         const datamatrix & D, const vector<ST::string> & modelv)
  {
  vector<FULLCOND*> fullcond_neu;
  unsigned i;
  fullcond_neu.push_back(fullcondp[0]);

  unsigned z = 1;
  while(z<fullcond_alle.size() && fullcond_alle[z]->get_fctype()==factor)
    {
    z = z + 1;
    }

  for(i=z;i<fullcond_alle.size();i++)
     {
     if(modell_alt[nf_size-2+i]==-1 && index==i)
        reset_fix(names_nonp[i-1][0]);
     if(modell_neu[nf_size-2+i]>0)
        {
        fullcond_neu.push_back(fullcond_alle[i]);
        if(i == index)
          fullcond_alle[i]->update_stepwise(modell_neu[nf_size-2+i]);
        }
     else if(modell_neu[nf_size-2+i]==0)    // && i == index)
        fullcond_alle[i]->reset_effect(0);
     else if(modell_neu[nf_size-2+i] == -1) // && i == index)
        {
        fullcond_alle[i]->reset_effect(0);
        if(i == index)
           fullcond_neu[0]->include_effect(names_nonp[i-1],
                                  fullcond_alle[i]->get_data_forfixedeffects());
        }
     }

  fullcondp = fullcond_neu;
  end[0] = fullcondp.size()-1;
  }


void STEPWISErun::fullcond_komplett(const vector<double> & m,
         const unsigned & nf_size,
         const vector<vector<ST::string> > & names_nonp, const datamatrix & D,
         const vector<ST::string> & modelv)
  {
  vector<FULLCOND*> fullcond_neu;
  unsigned i;

  fullcond_neu.push_back(fullcondp[0]);
  for(i=1;i<fullcond_alle.size();i++)
     {
     if(m[nf_size-2+i]>0)
        {
        fullcond_alle[i]->update_stepwise(m[nf_size-2+i]);
        fullcond_neu.push_back(fullcond_alle[i]);
        }
     else if(m[nf_size-2+i]==0)
        fullcond_alle[i]->reset_effect(0);
     else if(m[nf_size-2+i] == -1)
        {
        fullcond_alle[i]->reset_effect(0);
        fullcond_neu[0]->include_effect(names_nonp[i-1],
                                  fullcond_alle[i]->get_data_forfixedeffects());
        }
     }

  fullcondp = fullcond_neu;
  end[0] = fullcondp.size()-1;
  fullcondp[0]->set_effect_zero();
  }
  

void STEPWISErun::fix_komplett(const vector<double> &  modell,
         const vector<ST::string> & names_fixed,
         const vector<vector<ST::string> > & names_nonp, const datamatrix & D,
         const vector<ST::string> & modelv)
  {
  unsigned z;
  for(z=0;z<names_fixed.size()-1;z++)
     {
     if(modell[z]==0)
        reset_fix(names_fixed[z+1]);
     else if(modell[z]==-1)
       {
       unsigned i = 1;
       bool rein = false;
       while(i<fullcondp[0]->get_datanames().size() && rein==false)
            {
            if(fullcondp[0]->get_datanames()[i]==names_fixed[z+1])
              rein = true;
            i = i + 1;
            }
       if(rein==false)
          include_fix(names_fixed[z+1],D,modelv);
       }
     }

  for(z=names_fixed.size()-1;z<modell.size();z++)
     {
     bool gefunden = false;
     unsigned i = 1;
     while(i<fullcondp[0]->get_datanames().size() && gefunden==false)
        {
        if(fullcondp[0]->get_datanames()[i]==names_nonp[z-names_fixed.size()+1][0])    // ersetzen durch reset_fix?
           {
           gefunden = true;
           fullcondp[0]->reset_effect(i);
           }
           i = i + 1;
        }
     if(gefunden==true && names_nonp[z-names_fixed.size()+1].size()>1)
        {
        unsigned j;
        for(j=1;j<names_nonp[z-names_fixed.size()+1].size();j++)
           reset_fix(names_nonp[z-names_fixed.size()+1][j]);
        }
     }
  }


void STEPWISErun::reset_fix(const ST::string & name)
  {
  bool raus = false;
  unsigned j = 1;
  while(j<fullcondp[0]->get_datanames().size() && raus==false)
     {
     if(fullcondp[0]->get_datanames()[j]==name)
        {
        raus = true;
        fullcondp[0]->reset_effect(j);
        }
     j = j + 1;
     }
  }


void STEPWISErun::include_fix(const ST::string & name, const datamatrix & D,
         const vector<ST::string> & modelv)
  {
  bool gefunden = false;
  unsigned i = 0;
  while(i<modelv.size() && gefunden==false)
     {
     if(name==modelv[i])
        gefunden = true;
     i = i + 1;
     }
  vector<ST::string> help_name;
  help_name.push_back(name);
  fullcondp[0]->include_effect(help_name,datamatrix(D.getCol(i-1)));
  }


// -----------------------------------------------------------------------------
// ------- Funktionen für die Ausgabe im Output-Fenster ------------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::make_pause(void)
  { 
  #if defined(BORLAND_OUTPUT_WINDOW)
  Application->ProcessMessages();

  if (Frame->stop)
    {
    //break;
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("STEPWISE PROCEDURE TERMINATED BY USER BREAK\n");
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
    genoptions_mult[0]->out("STEPWISE PROCEDURE TERMINATED BY USER BREAK\n");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("Estimation results: none\n");
    genoptions_mult[0]->out("\n");
    //break;
    return true;
    }
  #endif

  return false; 
  }


void STEPWISErun::maketext(const ST::string & h, const vector<double> & m,
                          const double & a, ST::string & text,
                          const bool & neutext, const ST::string & tr,
                          const bool & datei)
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
    modeltext = modeltext + fullcondp[0]->get_effect();
    for(i=1;i<fullcondp.size();i++)
       modeltext = modeltext + " + " + fullcondp[i]->get_effect();
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
  }


void STEPWISErun::options_text(const int & number,
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
     
  for(i=1;i<fullcondp.size();i++)
     {
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  OPTIONS FOR NONPARAMETRIC TERM: "
          + names_nonp[i-1][0] + "\n");
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
     else if(fullcondp[i]->get_fctype()!=MCMC::factor)
         genoptions_mult[0]->out("  Number of different smoothing parameters on a logarithmic scale: "
            + ST::doubletostring(fullcondp[i]->get_number()) + "\n");
     unsigned j;
     for(j=0;j<startindex.size();j++)
          {
          if(lambdavec[i-1][startindex[j][i-1]]==0)
             genoptions_mult[0]->out("  Startvalue of the "
                + ST::doubletostring(j+1) + ". startmodel is \"effect excluded\" \n");
          else if(lambdavec[i-1][startindex[j][i-1]]==-1)
             genoptions_mult[0]->out("  Startvalue of the "
                + ST::doubletostring(j+1) + ". startmodel is the fixed effect \n");
          else
             {
             genoptions_mult[0]->out("  Startvalue of the smoothing parameter for the "
                + ST::doubletostring(j+1) + ". startmodel: "
                + ST::doubletostring(lambdavec[i-1][startindex[j][i-1]],6) + "\n");
             fullcondp[i]->update_stepwise(lambdavec[i-1][startindex[j][i-1]]);
             genoptions_mult[0]->out("  This is equivalent to degrees of freedom: "
                + ST::doubletostring(fullcondp[i]->compute_df(),6) + "\n");
             }
          }
     }
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("STEPWISE PROCEDURE STARTED \n");
  genoptions_mult[0]->out("\n");
  }


// -----------------------------------------------------------------------------
// ------- Funktionen für die Ausgabe im Tex-File ------------------------------
// -----------------------------------------------------------------------------

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
  outtex << "\n\\noindent {\\bf \\large Final Properties:}\\\\ \n\\\\" << endl;
  for(j=1;j<fullcond_alle.size();j++)
     {
     if(modell[names_f.size()-2+j]!=0 && modell[names_f.size()-2+j]!=-1)
         {
         vector<ST::string> prior = fullcond_alle[j]->get_priorassumptions();
         outtex << prior[0] << "\\\\" << endl
                << "smoothing parameter: $\\lambda = "
                << ST::doubletostring(modell[names_f.size()-2+j],6)
                << " \\,\\, \\hat{=} \\,\\, \\mbox{df} = "
                << ST::doubletostring(fullcond_alle[j]->compute_df(),6)
                << "$ \\\\ \n\\\\" << endl;
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

  char hcharu = '_';
  ST::string hstringu = "\\_";
  for(i=1;i<fullcondp[0]->get_datanames().size();i++)
     term = term + " + " + fullcondp[0]->get_datanames()[i].insert_string_char(hcharu,hstringu);
  for(i=1;i<fullcondp.size();i++)
     term = term + " + " + fullcondp[i]->get_term_symbolic();

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
  for(i=1;i<fullcond_alle.size();i++)
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
        if(lambdavec[i-1][startindex[j][i-1]]==0)
           {
           if(startindex.size()>1)
              outtex << ST::doubletostring(j+1) << ". ";
           outtex << "Startvalue is \\glqq effect excluded\\grqq \\\\ \n";
           }
        else if(lambdavec[i-1][startindex[j][i-1]]==-1)
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
                  << ST::doubletostring(lambdavec[i-1][startindex[j][i-1]],6)
                  << " \\,\\, \\hat{=} \\,\\, \\mbox{df} = ";
           fullcond_alle[i]->update_stepwise(lambdavec[i-1][startindex[j][i-1]]);
           outtex << ST::doubletostring(fullcond_alle[i]->compute_df(),6) << "$ \\\\ \n";
           }
        }
     outtex << "\\\\" << endl;
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


