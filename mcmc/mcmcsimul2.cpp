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

bool STEPWISErun::stepwise(const ST::string & procedure, const ST::string & minimum,
        const ST::string & crit, const int & stp, const ST::string & trac,
        const int & number, const ST::string & stam, const int & inc, const bool & finet,
        const bool & maveraging, int & fenster,
        const datamatrix & Da, const vector<ST::string> & modelvar,
        const ST::string & name, vector<FULLCOND*> & fullcond_z, ST::string & path)
  {

  D = Da;
  modelv = modelvar;
  algorithm = procedure;
  minim = minimum;
  criterion = crit;
  increment = inc;
  steps = stp;
  startmodel = stam;
  fine_tuning = finet;
  trace = trac;
  bool modelaveraging = maveraging;
  window = fenster;
  smoothing = "global";

  modell_alt.erase(modell_alt.begin(),modell_alt.end());
  vector<double> modell_final;
  double kriterium_final;
  ST::string text_final;

  ST::string tr_akt = "trace_on";
  posttitle.push_back("");

  lambdavec.erase(lambdavec.begin(),lambdavec.end());
  names_fixed.erase(names_fixed.begin(),names_fixed.end());
  names_nonp.erase(names_nonp.begin(),names_nonp.end());

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
  if( vcm_doppelt() == true)
     return true;

  vector<vector<unsigned> > startindex;
  vector<vector<double> > startfix;
  startwerte(startmodel,startindex,startfix);

  options_text(number,startfix,startindex,name);

  ST::string path_tex = path + "_model_summary.tex";
  outtex.open(path_tex.strtochar());
  make_graphics(name,startindex);

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
     abbruch = single_stepwise(startindex[i],startfix[i],true);

     if(abbruch==true)
        return true;

     if(minim=="adaptiv" || minim=="adaptiv_golden")
       {
       posteriormode(posttitle,true);
       kriterium_alt = compute_criterion();
       outcriterium << "B   " << ST::doubletostring(kriterium_alt,8) << endl;
       outmodels << endl << "B   " << ST::doubletostring(kriterium_alt,8) << endl;
       }

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

  ST::string header = "  Final Model:";
  fix_komplett(modell_final);
  fullcond_komplett(modell_final);
  tr_akt = "trace_on";
  maketext(header,modell_final,kriterium_final,text_final,false,tr_akt,false);
  genoptions_mult[0]->out("\n\n");
  kriterium_tex = kriterium_final;

  if(fine_tuning == true)
     abbruch = finetuning(modell_final);
     //abbruch = fine_local(modell_final);

  if(abbruch==true)
    return true;

  //likep_mult[0]->set_linpred_null();
  if(modelaveraging == true)
    {
    fullcond_z = fullcond_alle;
    modell_alt = modell_final;
    compute_average();
    }
  else
    {
    fullcond_z = fullcondp;
    for(i=0;i<fullcond_z.size();i++)
      fullcond_z[i]->set_fcnumber(i);
//posteriormode(posttitle,true);
//double mse = likep_mult[0]->compute_msep();
//genoptions_mult[0]->out("MSE(absatz) = " + ST::doubletostring(mse) + "\n");
//ST::string hilf = path + "_mse.res";
//ofstream out(hilf.strtochar());
//out << "MSE(absatz) = " << ST::doubletostring(mse) << endl;
    posteriormode(posttitle,false); // Problem: linearer Prädiktor bei "true" standardisiert! Hier wird zurückgerechnet!
    }
                                    // danach nicht mehr compute_criterion() aufrufen!!!
  make_tex_end(path,modell_final);

  // Files müssen wieder geschlossen werden!!!
  outtex.close();
  outcriterium.close();
  outmodels.close();  

  /*
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
  */
  
  return false;
  }


bool STEPWISErun::single_stepwise(const vector<unsigned> & start,
                            const vector<double> & startfix, const bool & tex)
  {
  modell_neu.erase(modell_neu.begin(),modell_neu.end());
  modellematrix.erase(modellematrix.begin(),modellematrix.end());
  steps_aktuell = 0;
  ST::string tr_akt = "trace_on";
  vector<vector<double> > startiteration;

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
  fix_komplett(modell_alt);
  fullcond_komplett(modell_alt);
  posteriormode(posttitle,true);

  kriterium_alt = compute_criterion();

  if(tex==true)
     {
     kriterium_tex = kriterium_alt;
     make_predictor();
     }

  kriterium_neu = kriterium_alt;
  outcriterium << steps_aktuell << "   " << ST::doubletostring(kriterium_neu,8) << endl;
  outmodels << steps_aktuell << "   " << ST::doubletostring(kriterium_neu,8) << "   ";
  ST::string header;
  fertig = false;    // überprüft, ob es noch nicht gerechnete Modelle gibt
  ST::string text_neu;

  bool abbruch = false;
  if(algorithm != "coorddescent")
    abbruch = stepfunctions();
  else
    abbruch = koordabstieg();

  if(abbruch==true)
    return true;

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
// ------------------- Funktionen, für Stepwise / Stepmin ----------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::stepfunctions(void)
  {
  ST::string tr_akt = "trace_on";
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
       ST::string header = "  Startmodel:";
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
       maketext(header,modell_alt,kriterium_alt,text_neu,neutext,
                tr_akt,neutext);
       text_alt = text_neu;

       vector<vector<double> > modeliteration;
       double kriterium;
       vector<double> kriteriumiteration2;
       unsigned j;

       if(algorithm == "stepwise")
         {
         unsigned z = stepwise_fixfactor(kriteriumiteration2,modeliteration,textiteration);
         stepwise_nonp(kriteriumiteration2,modeliteration,textiteration,z);
         }
       else // if(algorithm == "stepmin")
         {
         if(minim == "exact" || minim == "exact_golden")
           {
           unsigned z = stepwise_fixfactor(kriteriumiteration2,modeliteration,textiteration);
           minexact_nonp(kriteriumiteration2,modeliteration,textiteration,z);
           }
         else
           {
           fullcondp[0]->set_effect_zero();
           posteriormode(posttitle,true);
           step_minfix(kriteriumiteration2,modeliteration,textiteration);
           unsigned z = step_minfactor(kriteriumiteration2,modeliteration,textiteration);
           stepmin_nonp(kriteriumiteration2,modeliteration,textiteration,z);
           }
         }

       if(fertig==false)
         {
         kriterium = kriteriumiteration2[0];
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
         outcriterium << steps_aktuell << "   " << ST::doubletostring(kriterium_neu,8) << endl;
         outmodels << steps_aktuell << "   " << ST::doubletostring(kriterium_neu,8) << "   ";
         modellematrix.push_back(modeliteration);

         header = "\n\nBest Model of this iteration:";
         fix_komplett(modell_neu);
         fullcond_komplett(modell_neu);
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
         if(trace == "trace_on" || trace == "trace_minim")
           {
           genoptions_mult[0]->out("\n\n");
           genoptions_mult[0]->out("\n\n");
           genoptions_mult[0]->out("  There are no new models for this iteration! \n");
           }
         outcriterium << ST::inttostring(steps_aktuell) << "   " << ST::doubletostring(kriterium_neu,8) << endl;
         outmodels << ST::inttostring(steps_aktuell) << "   " << ST::doubletostring(kriterium_neu,8) << "   ";
         }
       if(trace == "trace_on" || trace == "trace_minim")
         {
         genoptions_mult[0]->out("\n\n");
         genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
         genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
         }

       if(make_pause() == true)
          return true;
       }

  if(minim == "apprexact")
    {
    minim = "exact";
    if(fertig==false)
      {
      kriterium_neu = kriterium_alt;
      modell_neu = modell_alt;
      fix_komplett(modell_alt);
      fullcond_komplett(modell_alt);
      }
    //if(trace == "trace_on" || trace == "trace_minim")
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out("  Beginning of the exact minimization! \n");
    fertig = false;
    stepfunctions();
    }

  return false;
  }


unsigned STEPWISErun::stepwise_fixfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration)
  {
  unsigned i;
  for(i=1;i<names_fixed.size();i++)
    {
    modell_neu = modell_alt;
    if(modell_alt[i-1]==-1)
      modell_neu[i-1]= 0;
    else if(modell_alt[i-1]==0)
      modell_neu[i-1] = -1;
    if(modelcomparison(modell_neu,modellematrix)==false)
      newmodel_fix(modell_neu[i-1],kriteriumiteration2,modeliteration,
                            textiteration,names_fixed[i]);
    }

  unsigned z = 1;
  while(z<fullcond_alle.size() && fullcond_alle[z]->get_fctype()==factor)
     {
     modell_neu = modell_alt;
     if(modell_alt[z+names_fixed.size()-2]==-1 && fullcond_alle[z]->get_forced()==false)
       modell_neu[z+names_fixed.size()-2]= 0;
     else if(modell_alt[z+names_fixed.size()-2]==0)
       modell_neu[z+names_fixed.size()-2] = -1;
     if(modelcomparison(modell_neu,modellematrix)==false)
       newmodel_factor(modell_neu[z+names_fixed.size()-2],z,kriteriumiteration2,
                 modeliteration,textiteration,names_nonp[z-1]);
     z = z + 1;
     }
  return z;
  }

void STEPWISErun::stepwise_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  unsigned i;
  for(i=z;i<fullcond_alle.size();i++)
    {
    unsigned sch;
    for(sch=1;sch<=unsigned(increment);sch++)
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
            newmodel_nonp(i,kriteriumiteration2,modeliteration,textiteration);
          }

       lambda_exist = false;
       modell_neu = modell_alt;
       if(index >= sch)
          lambda_exist = true;
       if(lambda_exist==true)
          {
          modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][index-sch];
          if(modelcomparison(modell_neu,modellematrix)==false)
            newmodel_nonp(i,kriteriumiteration2,modeliteration,textiteration);
          }
       }
    }
  }
  

void STEPWISErun::stepmin_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  unsigned i;
  for(i=z;i<fullcond_alle.size();i++)
    {
    modell_neu = modell_alt;
    unsigned lambda_ind;
    if((minim != "approx_golden" && minim != "adaptiv_golden") || lambdavec[z-1].size()<7)
      {
      vector<double> krit_fkt;
      if(modell_alt[i+names_fixed.size()-2]==0)
        stepmin_nonp_leer(i,krit_fkt,kriterium_alt);
      else if(modell_alt[i+names_fixed.size()-2]==-1)
        stepmin_nonp_fix(i,krit_fkt,kriterium_alt);
      else
        stepmin_nonp_nonp(i,krit_fkt,kriterium_alt);

      double kriterium_min = krit_fkt[0];
      unsigned j;
      lambda_ind = 0;
      for(j=1;j<krit_fkt.size();j++)  //berechnet den besten Wert
         {
         if(krit_fkt[j]<=kriterium_min)
           {
           kriterium_min = krit_fkt[j];
           lambda_ind = j;
           }
         }
      }
    else
      lambda_ind = golden_section(i,kriterium_alt);

    modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][lambda_ind];
    if(modell_neu[names_fixed.size()-2+i] != modell_alt[names_fixed.size()-2+i])
      {
      if(modelcomparison(modell_neu,modellematrix)==false)
        {
        newmodel_nonp(i,kriteriumiteration2,modeliteration,textiteration);
        // Stellt linearen Prädiktor wieder her. Besser wäre, den lin. Prädiktor zu speichern!!!
        fullcondp[0]->set_effect_zero();
        posteriormode(posttitle,true);
        }
      }
    }
  }

void STEPWISErun::minexact_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  unsigned i;
  unsigned j;
  for(i=z;i<fullcond_alle.size();i++)
    {
    modell_neu = modell_alt;
    unsigned lambda_ind;
    if(minim != "exact_golden" || lambdavec[z-1].size()<7)
      {
      vector<double> krit_fkt;
      double df = 0;
      if(modell_alt[i+names_fixed.size()-2]==0)
        {
        for(j=0;j<fullcondp.size();j++)
          df = df + fullcondp[j]->compute_df();
        fullcondp.push_back(fullcond_alle[i]);
        minexact_nonp_leer(i,krit_fkt,kriterium_alt,df);
        }
      else if(modell_alt[i+names_fixed.size()-2]==-1)
        {
        for(j=0;j<fullcondp.size();j++)
          df = df + fullcondp[j]->compute_df();
        df = df - 1;
        reset_fix(names_nonp[i-1][0]);
        fullcondp.push_back(fullcond_alle[i]);
        minexact_nonp_fix(i,krit_fkt,kriterium_alt,df);
        }
      else
        {
        for(j=0;j<fullcondp.size();j++)
          df = df + fullcondp[j]->compute_df();
        df = df - fullcond_alle[i]->compute_df();
        minexact_nonp_nonp(i,krit_fkt,kriterium_alt,df);
        }

      double kriterium_min = krit_fkt[0];
      lambda_ind = 0;
      for(j=1;j<krit_fkt.size();j++)  //berechnet den besten Wert
         {
         if(krit_fkt[j]<=kriterium_min)
           {
           kriterium_min = krit_fkt[j];
           lambda_ind = j;
           }
         }
      }
    else
      lambda_ind = golden_section(i,kriterium_alt);

    modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][lambda_ind];
    if(modell_neu[names_fixed.size()-2+i] != modell_alt[names_fixed.size()-2+i])
      {
      if(modelcomparison(modell_neu,modellematrix)==false)
        {
        newmodel_nonp(i,kriteriumiteration2,modeliteration,textiteration);
        // Stellt linearen Prädiktor wieder her.
        fullcondp[0]->set_effect_zero();
        posteriormode(posttitle,true);
        }
      }
    }
  }


// -----------------------------------------------------------------------------
// ------------------ Funktionen für Stepmin -----------------------------------
// -----------------------------------------------------------------------------

void STEPWISErun::step_minfix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration)
  {
  unsigned i;
  for(i=1;i<names_fixed.size();i++)
    {
    if(modell_alt[i-1]==-1)
      stepmin_fix_leer(kriteriumiteration2,modeliteration,textiteration,i);
    else if(modell_alt[i-1]==0)
      stepmin_leer_fix(kriteriumiteration2,modeliteration,textiteration,i);
    }
  }

void STEPWISErun::stepmin_fix_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration, unsigned & i)
  {
  double df = 0;
  unsigned y;
  for(y=0;y<fullcondp.size();y++)
    df = df + fullcondp[y]->compute_df();
  df = df - 1;

  fullcond_alle[0]->safe_const();
  reset_fix(names_fixed[i]);
  fullcond_alle[0]->posteriormode_const();
  kriterium_neu = criterion_min(df);
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    posteriormode(posttitle,true);
    double kriterium_test = criterion_min(df);
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_fixed[i] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    include_fix(names_fixed[i]);
    fullcondp[0]->set_effect_zero();
    posteriormode(posttitle,true);
    reset_fix(names_fixed[i]);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  " + names_fixed[i] + " Testvalue: approx \n");
     genoptions_mult[0]->out("        " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(kriterium_neu < kriterium_alt)
    {
    modell_neu[i-1] = 0;
    if(modelcomparison(modell_neu,modellematrix)==false)
      {
      newmodel(kriteriumiteration2,modeliteration,textiteration);
      include_fix(names_fixed[i]);
      fullcondp[0]->set_effect_zero();
      posteriormode(posttitle,true);
      }
    else
      {
      int c = column_for_fix(names_fixed[i]);
      vector<ST::string> name_help;
      name_help.push_back(names_fixed[i]);
      fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)));
      }
    modell_neu[i-1] = -1;
    }
  else
    {
    int c = column_for_fix(names_fixed[i]);
    vector<ST::string> name_help;
    name_help.push_back(names_fixed[i]);
    fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)));
    }
  }

void STEPWISErun::stepmin_leer_fix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration, unsigned & i)
  {
  double df = 0;
  unsigned y;
  for(y=0;y<fullcondp.size();y++)
    df = df + fullcondp[y]->compute_df();
  df = df + 1;

  fullcond_alle[0]->safe_const();
  int c = column_for_fix(names_fixed[i]);
  vector<ST::string> name_help;
  name_help.push_back(names_fixed[i]);
  fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)));
  kriterium_neu = criterion_min(df);
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    posteriormode(posttitle,true);
    double kriterium_test = criterion_min(df);
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_fixed[i] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    reset_fix(names_fixed[i]);
    fullcondp[0]->set_effect_zero();
    posteriormode(posttitle,true);
    fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)));    
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  " + names_fixed[i] + " Testvalue: approx \n");
     genoptions_mult[0]->out("        " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(kriterium_neu < kriterium_alt)
     {
     modell_neu[i-1] = -1;
     if(modelcomparison(modell_neu,modellematrix)==false)
       {
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       reset_fix(names_fixed[i]);
       fullcondp[0]->set_effect_zero();
       posteriormode(posttitle,true);
       }
     else
       reset_fix(names_fixed[i]);
     modell_neu[i-1] = 0;
     }
  else
     reset_fix(names_fixed[i]);
  }

unsigned STEPWISErun::step_minfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration)
  {
  unsigned z = 1;
  while(z<fullcond_alle.size() && fullcond_alle[z]->get_fctype()==factor)
     {
     if(modell_alt[z+names_fixed.size()-2]==-1 && fullcond_alle[z]->get_forced()==false)
       stepmin_factor_leer(kriteriumiteration2,modeliteration,textiteration,z);
     else if(modell_alt[z+names_fixed.size()-2]==0)
       stepmin_leer_factor(kriteriumiteration2,modeliteration,textiteration,z);
     z = z + 1;
     }
  return z;
  }

void STEPWISErun::stepmin_factor_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  double df = 0;
  unsigned i;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();
  df = df - fullcond_alle[z]->get_data_forfixedeffects().cols();

  fullcond_alle[0]->safe_const();
  for(i=0;i<names_nonp[z-1].size();i++)
    reset_fix(names_nonp[z-1][i]);
  fullcond_alle[0]->posteriormode_const();
  kriterium_neu = criterion_min(df);
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    posteriormode(posttitle,true);
    double kriterium_test = criterion_min(df);
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    fullcondp[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
    fullcondp[0]->set_effect_zero();
    posteriormode(posttitle,true);
    for(i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalue: approx \n");
     genoptions_mult[0]->out("        " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(kriterium_neu < kriterium_alt)
     {
     modell_neu[z+names_fixed.size()-2] = 0;
     if(modelcomparison(modell_neu,modellematrix)==false)
       {
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       fullcondp[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
       fullcondp[0]->set_effect_zero();
       posteriormode(posttitle,true);
       }
     else
       fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                  fullcond_alle[z]->get_data_forfixedeffects());
     modell_neu[z+names_fixed.size()-2] = -1;
     }
  else
     fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                 fullcond_alle[z]->get_data_forfixedeffects());
  }

void STEPWISErun::stepmin_leer_factor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  double df = 0;
  unsigned i;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();
  df = df + fullcond_alle[z]->get_data_forfixedeffects().cols();

  fullcond_alle[0]->safe_const();
  fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                 fullcond_alle[z]->get_data_forfixedeffects());
  kriterium_neu = criterion_min(df);
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    posteriormode(posttitle,true);
    double kriterium_test = criterion_min(df);
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    for(i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    fullcondp[0]->set_effect_zero();
    posteriormode(posttitle,true);
    fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                 fullcond_alle[z]->get_data_forfixedeffects());
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalue: approx \n");
     genoptions_mult[0]->out("        " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(kriterium_neu < kriterium_alt)
     {
     modell_neu[z+names_fixed.size()-2] = -1;
     if(modelcomparison(modell_neu,modellematrix)==false)
       {
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       for(i=0;i<names_nonp[z-1].size();i++)
         reset_fix(names_nonp[z-1][i]);
       fullcondp[0]->set_effect_zero();
       posteriormode(posttitle,true);
       }
     else
       {
       for(i=0;i<names_nonp[z-1].size();i++)
         reset_fix(names_nonp[z-1][i]);
       }
     modell_neu[z+names_fixed.size()-2] = 0;
     }
  else
     {
     for(i=0;i<names_nonp[z-1].size();i++)
       reset_fix(names_nonp[z-1][i]);
     }
  }

void STEPWISErun::stepmin_nonp_nonp(unsigned & z, vector<double> & krit_fkt,double & kriterium)
  {
  unsigned i;
  double df = 0;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();

  if(minim == "adaptiv" || minim == "adap_exact")
    {
    fullcond_alle[z]->update_stepwise(modell_alt[z+names_fixed.size()-2]);
    fullcond_alle[z]->posteriormode();
    kriterium = criterion_min(df);
    }

  df = df - fullcond_alle[z]->compute_df();
  fullcondp[0]->safe_const();
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
        fullcond_alle[0]->set_const_old();
        }
      else if(lambdavec[z-1][i]==-1)
        {
        fullcond_alle[z]->set_inthemodel(0);
        fullcond_alle[z]->reset_effect(0);
        fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                  fullcond_alle[z]->get_data_forfixedeffects());
        kriterium_versuch = criterion_min(df + 1);
        fullcond_alle[0]->set_const_old();
        reset_fix(names_nonp[z-1][0]);
        }
      else
        {
        fullcond_alle[z]->set_inthemodel(0);
        fullcond_alle[z]->reset_effect(0);
        fullcond_alle[0]->posteriormode_const();
        kriterium_versuch = criterion_min(df);
        fullcond_alle[0]->set_const_old();
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalues: approx \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(minim == "approx_control")
    {
    vector<double> kriterium_control;
    genoptions_mult[0]->out("\n\n");
    fullcond_alle[z]->set_inthemodel(1);
    minexact_nonp_nonp(z,kriterium_control,kriterium,df);
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalues: approx    exact \n");
    for(i=0;i<kriterium_control.size();i++)
      genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],6) + "   "
             + ST::doubletostring(kriterium_control[i],6) + "\n");
    }
  else
    {
    fullcond_alle[z]->set_inthemodel(1);
    fullcond_alle[z]->update_stepwise(modell_alt[z+names_fixed.size()-2]);
    fullcond_alle[z]->posteriormode();
    }
  }


void STEPWISErun::stepmin_nonp_fix(unsigned & z, vector<double> & krit_fkt, double & kriterium)
  {
  unsigned i;
  double df = 0;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();

  if(minim == "adaptiv" || minim == "adap_exact")
    {
    reset_fix(names_nonp[z-1][0]);
    fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                  fullcond_alle[z]->get_data_forfixedeffects());
    kriterium = criterion_min(df);
    }

  df = df - 1;
  fullcond_alle[0]->safe_const();
  reset_fix(names_nonp[z-1][0]);
  fullcondp.push_back(fullcond_alle[z]);
  fullcond_alle[z]->set_inthemodel(1);

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
        fullcond_alle[0]->set_const_old();
        }
      else
        {
        fullcond_alle[z]->set_inthemodel(0);
        fullcond_alle[z]->reset_effect(0);
        fullcond_alle[0]->posteriormode_const();
        kriterium_versuch = criterion_min(df);
        fullcond_alle[0]->set_const_old();
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalues: approx \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(minim == "approx_control")
    {
    vector<double> kriterium_control;
    fullcond_alle[z]->set_inthemodel(0);
    minexact_nonp_fix(z,kriterium_control,kriterium,df);
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalues: approx    exact \n");
    for(i=0;i<kriterium_control.size();i++)
      genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],6) + "   "
             + ST::doubletostring(kriterium_control[i],6) + "\n");
    }
  else
    {
    fullcond_alle[z]->set_inthemodel(0);
    fullcond_alle[z]->reset_effect(0);
    fullcondp.erase(fullcondp.end()-1,fullcondp.end());
    fullcondp[0]->posteriormode_single(names_nonp[z-1],
                                fullcond_alle[z]->get_data_forfixedeffects());
    }
  }


void STEPWISErun::stepmin_nonp_leer(unsigned & z, vector<double> & krit_fkt ,double & kriterium)
  {
  unsigned i;

  double df = 0;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();
  // if(minim == "adaptiv" || minim == "adap_exact")  hier nicht nötig (siehe koordmin_leer_fix) 

  fullcond_alle[0]->safe_const();
  fullcondp.push_back(fullcond_alle[z]);
  fullcond_alle[z]->set_inthemodel(1);

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
        fullcond_alle[0]->set_const_old();
        }
      else
        {
        fullcond_alle[z]->set_inthemodel(0);
        fullcond_alle[z]->reset_effect(0);
        fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                   fullcond_alle[z]->get_data_forfixedeffects());
        kriterium_versuch = criterion_min(df + 1);
        reset_fix(names_nonp[z-1][0]);
        fullcond_alle[0]->set_const_old();          
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalues: approx \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(minim == "approx_control")
    {
    vector<double> kriterium_control;
    fullcond_alle[z]->set_inthemodel(0);
    minexact_nonp_leer(z,kriterium_control,kriterium,df);
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalues: approx    exact \n");
    for(i=0;i<kriterium_control.size();i++)
      genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],6) + "   "
             + ST::doubletostring(kriterium_control[i],6) + "\n");
    }
  else
    {
    fullcond_alle[z]->set_inthemodel(0);
    fullcond_alle[z]->reset_effect(0);
    fullcondp.erase(fullcondp.end()-1,fullcondp.end());
    }
  }

//------------------------------------------------------------------------------------

void STEPWISErun::minexact_nonp_nonp(unsigned & z, vector<double> & krit_fkt,
                  double & kriterium, double & df)
  {
  unsigned i;
  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=modell_alt[z+names_fixed.size()-2])
      {
      double kriterium_versuch;
      if(lambdavec[z-1][i]!=-1 && lambdavec[z-1][i]!=0)
        {
        fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
        fullcondp[0]->set_effect_zero();
        posteriormode(posttitle,true);
        kriterium_versuch = criterion_min(df + fullcond_alle[z]->compute_df());
        }
      else if(lambdavec[z-1][i]==-1)
        {
        fullcond_alle[z]->set_inthemodel(0);
        vector<FULLCOND*> fullcond_start = fullcondp;
        vector<double> modell1 = modell_alt;
        modell1[z+names_fixed.size()-2] = -1;
        fullcond_einzeln(modell1,modell_alt,z);  // hier muß der Fullcond-Vekor angepaßt werden!!!
        fullcondp[0]->set_effect_zero();
        posteriormode(posttitle,true);
        kriterium_versuch = criterion_min(df + 1);
        fullcondp = fullcond_start;
        end[0] = fullcondp.size()-1;
        reset_fix(names_nonp[z-1][0]);
        }
      else
        {
        fullcond_alle[z]->set_inthemodel(0);
        vector<FULLCOND*> fullcond_start = fullcondp;
        vector<double> modell1 = modell_alt;
        modell1[z+names_fixed.size()-2] = 0;
        fullcond_einzeln(modell1,modell_alt,z);  // hier muß der Fullcond-Vekor angepaßt werden!!!
        fullcondp[0]->set_effect_zero();
        posteriormode(posttitle,true);
        kriterium_versuch = criterion_min(df);
        fullcondp = fullcond_start;
        end[0] = fullcondp.size()-1;
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }
  fullcond_alle[z]->set_inthemodel(1);
  fullcond_alle[z]->update_stepwise(modell_alt[z+names_fixed.size()-2]);
  fullcondp[0]->set_effect_zero();
  posteriormode(posttitle,true);

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalues: exact \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }
  }


void STEPWISErun::minexact_nonp_fix(unsigned & z, vector<double> & krit_fkt,
          double & kriterium, double & df)
  {
  unsigned i;
  end[0] = fullcondp.size()-1;
  fullcond_alle[z]->set_inthemodel(1);
  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=-1)
      {
      double kriterium_versuch;
      if(lambdavec[z-1][i]!=0)
        {
        fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
        fullcondp[0]->set_effect_zero();
        posteriormode(posttitle,true);
        kriterium_versuch = criterion_min(df + fullcond_alle[z]->compute_df());
        }
      else
        {
        fullcond_alle[z]->set_inthemodel(0);
        vector<FULLCOND*> fullcond_start = fullcondp;
        fullcondp.erase(fullcondp.end()-1,fullcondp.end());
        end[0] = fullcondp.size()-1;
        fullcond_alle[z]->reset_effect(0);
        fullcondp[0]->set_effect_zero();
        posteriormode(posttitle,true);
        kriterium_versuch = criterion_min(df);
        fullcondp = fullcond_start;
        end[0] = fullcondp.size()-1;
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }
  fullcond_alle[z]->set_inthemodel(0);
  fullcond_alle[z]->reset_effect(0);
  fullcondp[0]->include_effect(names_nonp[z-1],
                                fullcond_alle[z]->get_data_forfixedeffects());
  fullcondp.erase(fullcondp.end()-1,fullcondp.end());
  end[0] = fullcondp.size()-1;
  fullcondp[0]->set_effect_zero();
  posteriormode(posttitle,true);

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalues: exact \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }
  }


void STEPWISErun::minexact_nonp_leer(unsigned & z, vector<double> & krit_fkt,
                  double & kriterium, double & df)
  {
  unsigned i;
  end[0] = fullcondp.size()-1;
  fullcond_alle[z]->set_inthemodel(1);
  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=0)
      {
      double kriterium_versuch;
      if(lambdavec[z-1][i]!=-1)
        {
        fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
        fullcondp[0]->set_effect_zero();
        posteriormode(posttitle,true);
        kriterium_versuch = criterion_min(df + fullcond_alle[z]->compute_df());
        }
      else
        {
        fullcond_alle[z]->set_inthemodel(0);
        vector<FULLCOND*> fullcond_start = fullcondp;
        fullcondp.erase(fullcondp.end()-1,fullcondp.end());
        end[0] = fullcondp.size()-1;
        fullcond_alle[z]->reset_effect(0);
        fullcondp[0]->include_effect(names_nonp[z-1],
                                fullcond_alle[z]->get_data_forfixedeffects());
        fullcondp[0]->set_effect_zero();
        posteriormode(posttitle,true);
        kriterium_versuch = criterion_min(df + 1);
        fullcondp = fullcond_start;
        end[0] = fullcondp.size()-1;
        reset_fix(names_nonp[z-1][0]);
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }
  fullcond_alle[z]->set_inthemodel(0);
  fullcond_alle[z]->reset_effect(0);
  fullcondp.erase(fullcondp.end()-1,fullcondp.end());
  fullcondp[0]->set_effect_zero();
  end[0] = fullcondp.size()-1;
  posteriormode(posttitle,true);

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalues: exact \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }
  }

//-----------------------------------------------------------------------------------

double STEPWISErun::criterion_min(const double & df)
  {
  double kriterium;

  /*double df2 = df;
  if(df2>=499)      // Versuch
    df2 = 498;*/

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

//     Versuch, ob penalisierte statt normaler Log-Likelihood verwendet werden sollte
/*
unsigned j;
double sum = 0;
for(j=1;j<fullcondp.size();j++)
  {
  sum += fullcondp[j]->compute_penalterm();
  }
sum /= likep_mult[0]->compute_rss();
sum *= likep_mult[0]->get_nrobs_wpw();
kriterium += sum;
*/

  return kriterium;
  }


// -----------------------------------------------------------------------------
// ------------------ Funktionen für Koordinatenmethode ------------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::koordabstieg(void)
  {
      // Schleife für Minimierung
  ST::string tr_akt = "trace_on";
  ST::string text_neu;
  bool eins = true;
  double kriterium_aktuell;
//  modell_uralt = modell_neu;
  while(kriterium_neu <= kriterium_alt && fertig==false && steps_aktuell<steps)
       {
       steps_aktuell = steps_aktuell + 1;
       vector<ST::string> textiteration;
       fertig = true;
       kriterium_aktuell = kriterium_neu;
       kriterium_alt = kriterium_neu;
       modell_alt = modell_neu;

       ST::string header = "  Startmodel:";
       maketext(header,modell_alt,kriterium_alt,text_neu,true,tr_akt,true);
       text_alt = text_neu;

       if(eins == true)
         {
         eins = false;
         if(trace == "trace_off")
           tr_akt = "trace_off";
         }

       vector<vector<double> > modeliteration;
       vector<double> kriteriumiteration2;

       if(minim == "exact" || minim == "exact_golden")
         {
         unsigned z = koordexact_fixfactor(kriteriumiteration2,modeliteration,
                            textiteration,kriterium_aktuell);
         koordexact_nonp(kriteriumiteration2,modeliteration,textiteration,z,kriterium_aktuell);
         }

       else
         {
         koord_minfix(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell);
         unsigned z = koord_minfactor(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell);
         koord_minnonp(kriteriumiteration2,modeliteration,textiteration,z,kriterium_aktuell);

         if((minim == "adaptiv" || minim == "adap_exact" || minim == "adaptiv_golden")
                                      && likep_mult[0]->get_family() != "Gaussian")
           likep_mult[0]->compute_iwls();
           
         if((minim == "adaptiv" || minim == "adap_exact" || minim == "adaptiv_golden")
                                     && steps_aktuell>=3)
           {
           unsigned hilfe = modellematrix.size()-3;
           if(modell_alt == modellematrix[hilfe][modellematrix[hilfe].size()-1])
             fertig = true;
           } 
         }

       kriterium_neu = kriterium_aktuell;

       if(fertig==false)
         {
         outcriterium << steps_aktuell << "   " << ST::doubletostring(kriterium_neu,8) << endl;
         outmodels << steps_aktuell << "   " << ST::doubletostring(kriterium_neu,8) << "   ";
         modellematrix.push_back(modeliteration);

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
         if(trace == "trace_on" || trace == "trace_minim")
           {
           genoptions_mult[0]->out("\n\n");
           genoptions_mult[0]->out("\n\n");
           genoptions_mult[0]->out("  There are no new models for this iteration! \n");
           }
         outcriterium << ST::inttostring(steps_aktuell) << "   " << ST::doubletostring(kriterium_neu,8) << endl;
         outmodels << ST::inttostring(steps_aktuell) << "   " << ST::doubletostring(kriterium_neu,8) << "   ";
         }
       if(trace == "trace_on" || trace == "trace_minim")
         {
         genoptions_mult[0]->out("\n\n");
         genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
         genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
         }

       if(make_pause() == true)
          return true;
       }

  if(minim == "apprexact" || minim == "adap_exact")
    {
    if(fertig==false)
      {
      modell_neu = modell_alt;
      fix_komplett(modell_alt);
      fullcond_komplett(modell_alt);
      if(minim == "adap_exact")
        {
        posteriormode(posttitle,true);
        kriterium_alt = compute_criterion();
        }
      kriterium_neu = kriterium_alt;
      }
    else if(fertig == true && minim == "adap_exact")
      {
      posteriormode(posttitle,true);
      kriterium_alt = compute_criterion();
      kriterium_neu = kriterium_alt;
      }
    minim = "exact";
    //if(trace == "trace_on" || trace == "trace_minim")
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out("  Beginning of the exact minimization! \n");
    fertig = false;
    koordabstieg();
    }

  return false;
  }


void STEPWISErun::koord_minfix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell)
  {
  unsigned i;
  for(i=1;i<names_fixed.size();i++)
    {
    if(modell_alt[i-1]==-1)
      koord_fix_leer(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell,i);
    else if(modell_alt[i-1]==0)
      koord_leer_fix(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell,i);
    modell_alt = modell_neu;
    }
  }

void STEPWISErun::koord_fix_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & i)
  {
  double df = 0;
  unsigned y;
  for(y=0;y<fullcondp.size();y++)
    df = df + fullcondp[y]->compute_df();

  double kriterium_adaptiv = kriterium_aktuell;
  if(minim == "adaptiv" || minim == "adap_exact")
    {
    reset_fix(names_fixed[i]);          // "posteriormode_single" ändern: ohne include_effect! (dann ist raus- + reinnehmen überflüssig!)
    int c = column_for_fix(names_fixed[i]);
    vector<ST::string> name_help;
    name_help.push_back(names_fixed[i]);
    fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)));
    kriterium_aktuell = criterion_min(df);
    }

  df = df - 1;
  modell_neu[i-1] = 0;
  fullcond_alle[0]->safe_const();
  reset_fix(names_fixed[i]);
  fullcond_alle[0]->posteriormode_const();
  kriterium_neu = criterion_min(df);
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    posteriormode(posttitle,true);
    double kriterium_test = criterion_min(df);
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_fixed[i] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    include_fix(names_fixed[i]);
    posteriormode(posttitle,true);
    reset_fix(names_fixed[i]);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  " + names_fixed[i] + " Testvalue: approx \n");
     if(minim=="adaptiv" || minim == "adap_exact")
       genoptions_mult[0]->out("        " + ST::doubletostring(kriterium_aktuell,6) + "\n");
     genoptions_mult[0]->out("        " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

 if(kriterium_neu < kriterium_aktuell && minim != "adaptiv" && minim != "adap_exact")
    {
    bool neu = modelcomparison(modell_neu,modellematrix);
    if(neu==false)
      {
      newmodel(kriteriumiteration2,modeliteration,textiteration);
      kriterium_neu = kriteriumiteration2[kriteriumiteration2.size()-1];
      }
    if(neu==true || kriterium_aktuell < kriterium_neu)
      {
      int c = column_for_fix(names_fixed[i]);
      vector<ST::string> name_help;
      name_help.push_back(names_fixed[i]);
      fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)));
      modell_neu[i-1] = -1;
      if(kriterium_aktuell < kriterium_neu) // verhindert, daß "approx" schlechter wird!
        posteriormode(posttitle,true);
      }
    }

 if(kriterium_aktuell >= kriterium_neu)
   kriterium_aktuell = kriterium_neu;
 else //if(kriterium_neu > kriterium_aktuell)
    {
    int c = column_for_fix(names_fixed[i]);
    vector<ST::string> name_help;
    name_help.push_back(names_fixed[i]);
    fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)));   
    modell_neu[i-1] = -1;
    }

  if(minim == "adaptiv" || minim == "adap_exact")
    {
    if(fabs((kriterium_adaptiv - kriterium_aktuell)/kriterium_adaptiv) >= pow10(-6))
      fertig = false;
    if(modell_alt[i-1] != modell_neu[i-1] && (trace == "trace_on" || trace == "trace_minim"))
      {
      ST::string text;
      maketext("  Trial:",modell_neu,kriterium_aktuell,text,true,trace,false);
      }
    kriterium_alt = kriterium_aktuell;
    modell_alt[i-1] = modell_neu[i-1];
    modeliteration.push_back(modell_alt);
    }
  }

void STEPWISErun::koord_leer_fix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & i)
  {
  double df = 0;
  unsigned y;
  for(y=0;y<fullcondp.size();y++)
    df = df + fullcondp[y]->compute_df();
  double kriterium_adaptiv = kriterium_aktuell;
  // if(minim == "adaptiv" || minim == "adap_exact")   hier nicht nötig, weil Intercept immer erneuert wird!

  df = df + 1;
  modell_neu[i-1] = -1;
  fullcond_alle[0]->safe_const();
  int c = column_for_fix(names_fixed[i]);
  vector<ST::string> name_help;
  name_help.push_back(names_fixed[i]);
  fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)));
  kriterium_neu = criterion_min(df);
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    posteriormode(posttitle,true);
    double kriterium_test = criterion_min(df);
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_fixed[i] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    reset_fix(names_fixed[i]);
    fullcondp[0]->set_effect_zero();
    posteriormode(posttitle,true);
    fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)));
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  " + names_fixed[i] + " Testvalue: approx \n");
     genoptions_mult[0]->out("        " + ST::doubletostring(kriterium_neu,6) + "\n");
     if(minim=="adaptiv" || minim == "adap_exact")
       genoptions_mult[0]->out("        " + ST::doubletostring(kriterium_aktuell,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(kriterium_neu < kriterium_aktuell && minim != "adaptiv" && minim != "adap_exact")
     {
     bool neu = modelcomparison(modell_neu,modellematrix);
     if(neu==false)
       {
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       kriterium_neu = kriteriumiteration2[kriteriumiteration2.size()-1];
       }
     if(neu==true || kriterium_aktuell < kriterium_neu)
       {
       reset_fix(names_fixed[i]);
       modell_neu[i-1] = 0;
       if(kriterium_aktuell < kriterium_neu)
         posteriormode(posttitle,true);
       }
     }

  if(kriterium_aktuell >= kriterium_neu)
    kriterium_aktuell = kriterium_neu;
  else //if(kriterium_neu >= kriterium_aktuell)
    {
    reset_fix(names_fixed[i]);
    modell_neu[i-1] = 0;
    }

  if(minim == "adaptiv" || minim == "adap_exact")
    {
    if(fabs((kriterium_adaptiv - kriterium_aktuell)/kriterium_adaptiv) >= pow10(-6))
      fertig = false;
    if(modell_alt[i-1] != modell_neu[i-1] && (trace == "trace_on" || trace == "trace_minim"))
      {
      ST::string text;
      maketext("  Trial:",modell_neu,kriterium_aktuell,text,true,trace,false);
      }
    kriterium_alt = kriterium_aktuell;
    modell_alt[i-1] = modell_neu[i-1];
    modeliteration.push_back(modell_alt);
    }
  }

unsigned STEPWISErun::koord_minfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell)
  {
  unsigned z = 1;
  while(z<fullcond_alle.size() && fullcond_alle[z]->get_fctype()==factor)
     {
     if(modell_alt[z+names_fixed.size()-2]==-1 && fullcond_alle[z]->get_forced()==false)
       koord_factor_leer(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell,z);
     else if(modell_alt[z+names_fixed.size()-2]==0)
       koord_leer_factor(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell,z);
     modell_alt = modell_neu;
     z = z + 1;
     }
  return z;
  }

void STEPWISErun::koord_factor_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & z)
  {
  double df = 0;
  unsigned i;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();

  double kriterium_adaptiv = kriterium_aktuell;
  if(minim == "adaptiv" || minim == "adap_exact")
    {
    for(i=0;i<names_nonp[z-1].size();i++)    // "posteriormode_single" ändern: ohne include_effect! (dann ist raus- + reinnehmen überflüssig!)
      reset_fix(names_nonp[z-1][i]);
    fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                        fullcond_alle[z]->get_data_forfixedeffects());
    kriterium_aktuell = criterion_min(df);
    }

  df = df - fullcond_alle[z]->get_data_forfixedeffects().cols();
  modell_neu[z+names_fixed.size()-2] = 0;
  fullcond_alle[0]->safe_const();
  for(i=0;i<names_nonp[z-1].size();i++)
    reset_fix(names_nonp[z-1][i]);
  fullcond_alle[0]->posteriormode_const();
  kriterium_neu = criterion_min(df);
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    posteriormode(posttitle,true);
    double kriterium_test = criterion_min(df);
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    fullcondp[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
    fullcondp[0]->set_effect_zero();
    posteriormode(posttitle,true);
    for(i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalue: approx \n");
     if(minim=="adaptiv" || minim == "adap_exact")
       genoptions_mult[0]->out("        " + ST::doubletostring(kriterium_aktuell,6) + "\n");
     genoptions_mult[0]->out("        " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(kriterium_neu < kriterium_aktuell && minim != "adaptiv" && minim != "adap_exact")
     {
     bool neu = modelcomparison(modell_neu,modellematrix);
     if(neu==false)
       {
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       kriterium_neu = kriteriumiteration2[kriteriumiteration2.size()-1];
       }
     if(neu==true || kriterium_aktuell < kriterium_neu)
       {
       fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                        fullcond_alle[z]->get_data_forfixedeffects());
       modell_neu[z+names_fixed.size()-2] = -1;
       if(kriterium_aktuell < kriterium_neu)
         posteriormode(posttitle,true);
       }
     }

  if(kriterium_aktuell >= kriterium_neu)
     kriterium_aktuell = kriterium_neu;
  else //if(kriterium_neu >= kriterium_aktuell)
    {
    fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                     fullcond_alle[z]->get_data_forfixedeffects());
    modell_neu[z+names_fixed.size()-2] = -1;
    }

  if(minim == "adaptiv" || minim == "adap_exact")
    {
    if(fabs((kriterium_adaptiv - kriterium_aktuell)/kriterium_adaptiv) >= pow10(-6))
      fertig = false;
    if(modell_alt[z+names_fixed.size()-2] != modell_neu[z+names_fixed.size()-2]
                    && (trace == "trace_on" || trace == "trace_minim"))
      {
      ST::string text;
      maketext("  Trial:",modell_neu,kriterium_aktuell,text,true,trace,false);
      }
    kriterium_alt = kriterium_aktuell;
    modell_alt[z+names_fixed.size()-2] = modell_neu[z+names_fixed.size()-2];
    modeliteration.push_back(modell_alt);
    }
  }

void STEPWISErun::koord_leer_factor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & z)
  {
  double df = 0;
  unsigned i;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();
  double kriterium_adaptiv = kriterium_aktuell;
  // if(minim == "adaptiv" || minim == "adap_exact")   hier überflüssig (siehe oben)!

  df = df + fullcond_alle[z]->get_data_forfixedeffects().cols();
  modell_neu[z+names_fixed.size()-2] = -1;
  fullcond_alle[0]->safe_const();
  fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                 fullcond_alle[z]->get_data_forfixedeffects());
  kriterium_neu = criterion_min(df);
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    posteriormode(posttitle,true);
    double kriterium_test = criterion_min(df);
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    for(i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    fullcondp[0]->set_effect_zero();
    posteriormode(posttitle,true);
    fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                  fullcond_alle[z]->get_data_forfixedeffects());
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalue: approx \n");
     genoptions_mult[0]->out("        " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out("\n");
     if(minim=="adaptiv" || minim == "adap_exact")
       genoptions_mult[0]->out("        " + ST::doubletostring(kriterium_aktuell,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(kriterium_neu < kriterium_aktuell && minim != "adaptiv" && minim != "adap_exact")
     {
     bool neu = modelcomparison(modell_neu,modellematrix);
     if(neu==false)
       {
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       kriterium_neu = kriteriumiteration2[kriteriumiteration2.size()-1];
       }
     if(neu==true || kriterium_aktuell < kriterium_neu)
       {
       for(i=0;i<names_nonp[z-1].size();i++)
         reset_fix(names_nonp[z-1][i]);
       modell_neu[z+names_fixed.size()-2] = 0;
       if(kriterium_aktuell < kriterium_neu)
         posteriormode(posttitle,true);
       }
     }

  if(kriterium_aktuell >= kriterium_neu)
    kriterium_aktuell = kriterium_neu;
  else if(kriterium_neu >= kriterium_aktuell)
    {
    for(i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    modell_neu[z+names_fixed.size()-2] = 0;
    }

  if(minim == "adaptiv" || minim == "adap_exact")
    {
    if(fabs((kriterium_adaptiv - kriterium_aktuell)/kriterium_adaptiv) >= pow10(-6))
      fertig = false;
    if(modell_alt[z+names_fixed.size()-2] != modell_neu[z+names_fixed.size()-2]
                    && (trace == "trace_on" || trace == "trace_minim"))
      {
      ST::string text;
      maketext("  Trial:",modell_neu,kriterium_aktuell,text,true,trace,false);
      }
    kriterium_alt = kriterium_aktuell;
    modell_alt[z+names_fixed.size()-2] = modell_neu[z+names_fixed.size()-2];
    modeliteration.push_back(modell_alt);
    }
  }

void STEPWISErun::koord_minnonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z, double & kriterium_aktuell)
  {
  unsigned i;
  for(i=z;i<fullcond_alle.size();i++)
    {
    unsigned lambda_ind;
    double kriterium_min;
    double kriterium_test = kriterium_aktuell;
    if((minim != "approx_golden" && minim != "adaptiv_golden") || lambdavec[z-1].size()<7)
      {
      vector<double> krit_fkt;
      if(modell_alt[i+names_fixed.size()-2]==0)
        stepmin_nonp_leer(i,krit_fkt,kriterium_aktuell);
      else if(modell_alt[i+names_fixed.size()-2]==-1)
        stepmin_nonp_fix(i,krit_fkt,kriterium_aktuell);
      else
        stepmin_nonp_nonp(i,krit_fkt,kriterium_aktuell);

      kriterium_min = krit_fkt[0];
      unsigned j;
      lambda_ind = 0;
      for(j=1;j<krit_fkt.size();j++)  //berechnet den besten Wert
         {
         if(krit_fkt[j]<=kriterium_min)
           {
           kriterium_min = krit_fkt[j];
           lambda_ind = j;
           }
         }
      }
    else
      {
      kriterium_min = kriterium_aktuell;
      lambda_ind = golden_section(i,kriterium_min);
      }

    modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][lambda_ind];
    if(minim != "adaptiv" && minim != "adap_exact" && minim != "adaptiv_golden")
      {
      if(modell_neu[names_fixed.size()-2+i] != modell_alt[names_fixed.size()-2+i])
        {
        bool neu = modelcomparison(modell_neu,modellematrix);
        fullcond_einzeln(modell_neu,modell_alt,i);
        if(neu==false)
          {
          fullcondp[0]->set_effect_zero();
          newmodel(kriteriumiteration2,modeliteration,textiteration);
          kriterium_test = kriteriumiteration2[kriteriumiteration2.size()-1];
          }
        if(kriterium_aktuell > kriterium_test)
          {
          modell_alt = modell_neu;
          kriterium_aktuell = kriterium_test;
          }
        else //if(kriterium_aktuell <= kriterium_test)
          {
          if( (trace == "trace_minim" || trace == "trace_on") && neu == false)
            genoptions_mult[0]->out("\n\n  Trial won't become the new model! \n");
          fullcond_einzeln(modell_alt,modell_neu,i);
          modell_neu = modell_alt;
          posteriormode(posttitle,true);
          }
        }
      }
    else //if(minim == "adaptiv" || minim == "adap_exact")
      {
      if(modell_alt[names_fixed.size()-2+i] != modell_neu[names_fixed.size()-2+i])
        {
        fullcond_einzeln(modell_neu,modell_alt,i);
        if(modell_neu[names_fixed.size()-2+i] == 0)
          fullcond_alle[0]->posteriormode_const();
        else if(modell_neu[names_fixed.size()-2+i] == -1)
          {
          fullcond_alle[i]->reset_effect(0);
          reset_fix(names_nonp[i-1][0]);
          fullcond_alle[0]->posteriormode_single(names_nonp[i-1],
                                 fullcond_alle[i]->get_data_forfixedeffects());
          }
        else
          {
          fullcond_alle[i]->update_stepwise(modell_neu[names_fixed.size()-2+i]);
          fullcond_alle[i]->posteriormode();
          }
        if(trace == "trace_on" || trace == "trace_minim")
          {
          ST::string text;
          maketext("  Trial:",modell_neu,kriterium_min,text,true,trace,false);
          }
        }
      //else
      //  maketext(header,modell_neu,kriterium_alt,text_neu,true,tr_akt,true);
      modell_alt = modell_neu;
      kriterium_aktuell = kriterium_min;
      kriterium_alt = kriterium_aktuell;
      if(fabs((kriterium_test - kriterium_aktuell)/kriterium_test) >= pow10(-6))
        fertig = false;
      modeliteration.push_back(modell_alt);
      }
    }
  }


//------------------------------------------------------------------------------

unsigned STEPWISErun::koordexact_fixfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell)
  {
  bool help_fertig = true;
  unsigned i;
  for(i=1;i<names_fixed.size();i++)
     {
     if(modell_alt[i-1]==-1)
       modell_neu[i-1]= 0;
     else if(modell_alt[i-1]==0)
       modell_neu[i-1] = -1;
     if(modelcomparison(modell_neu,modellematrix)==false)
       {
       if(modell_neu[i-1]==0)
         reset_fix(names_fixed[i]);
       else
         include_fix(names_fixed[i]);
       fullcondp[0]->set_effect_zero();
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       if(kriteriumiteration2[kriteriumiteration2.size()-1] > kriterium_aktuell)
          {
          if(modell_neu[i-1]==0)
            include_fix(names_fixed[i]);
          else
            reset_fix(names_fixed[i]);
          modell_neu = modell_alt;
          }
       else
          {
          help_fertig = false;
          modell_alt = modell_neu;
          kriterium_aktuell = kriteriumiteration2[kriteriumiteration2.size()-1];
          }
       }
     else
       modell_neu = modell_alt;
     }

  unsigned z = 1;
  while(z<fullcond_alle.size() && fullcond_alle[z]->get_fctype()==factor)
     {
     if(modell_alt[z+names_fixed.size()-2]==-1 && fullcond_alle[z]->get_forced()==false)
       modell_neu[z+names_fixed.size()-2]= 0;
     else if(modell_alt[z+names_fixed.size()-2]==0)
       modell_neu[z+names_fixed.size()-2] = -1;
     if(modelcomparison(modell_neu,modellematrix)==false)
        {
        if(modell_neu[z+names_fixed.size()-2]==0)
          {
          for(i=0;i<names_nonp[z-1].size();i++)
            reset_fix(names_nonp[z-1][i]);
          }
        else
          fullcondp[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
        fullcondp[0]->set_effect_zero();
        newmodel(kriteriumiteration2,modeliteration,textiteration);
        if(kriteriumiteration2[kriteriumiteration2.size()-1] > kriterium_aktuell)
           {
           if(modell_neu[z+names_fixed.size()-2]==0)
             fullcondp[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
           else
             {
             for(i=0;i<names_nonp[z-1].size();i++)
               reset_fix(names_nonp[z-1][i]);
             }
           modell_neu = modell_alt;
           }
        else
           {
           help_fertig = false;
           modell_alt = modell_neu;
           kriterium_aktuell = kriteriumiteration2[kriteriumiteration2.size()-1];
           }
        }
     else
        modell_neu = modell_alt;
     z = z + 1;
     }
  fertig = help_fertig;
  return z;
  }

void STEPWISErun::koordexact_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z, double & kriterium_aktuell)
  {
  unsigned i;
  unsigned j;
  for(i=z;i<fullcond_alle.size();i++)
    {
    modell_neu = modell_alt;
    unsigned lambda_ind;
    if(minim != "exact_golden" || lambdavec[z-1].size()<7)
      {
      vector<double> krit_fkt;
      double df = 0;
      if(modell_alt[i+names_fixed.size()-2]==0)
        {
        for(j=0;j<fullcondp.size();j++)
          df = df + fullcondp[j]->compute_df();
        fullcondp.push_back(fullcond_alle[i]);
        minexact_nonp_leer(i,krit_fkt,kriterium_aktuell,df);
        }
      else if(modell_alt[i+names_fixed.size()-2]==-1)
        {
        for(j=0;j<fullcondp.size();j++)
          df = df + fullcondp[j]->compute_df();
        df = df - 1;
        reset_fix(names_nonp[i-1][0]);
        fullcondp.push_back(fullcond_alle[i]);
        minexact_nonp_fix(i,krit_fkt,kriterium_aktuell,df);
        }
      else
        {
        for(j=0;j<fullcondp.size();j++)
          df = df + fullcondp[j]->compute_df();
        df = df - fullcond_alle[i]->compute_df();
        minexact_nonp_nonp(i,krit_fkt,kriterium_aktuell,df);
        }

      double kriterium_min = krit_fkt[0];
      lambda_ind = 0;
      for(j=1;j<krit_fkt.size();j++)  //berechnet den besten Wert
         {
         if(krit_fkt[j]<=kriterium_min)
           {
           kriterium_min = krit_fkt[j];
           lambda_ind = j;
           }
         }
      }
    else
      lambda_ind = golden_section(i,kriterium_aktuell);

    modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][lambda_ind];
    if(modell_neu[names_fixed.size()-2+i] != modell_alt[names_fixed.size()-2+i])
      {
      if(modelcomparison(modell_neu,modellematrix)==false)
        {
        fullcond_einzeln(modell_neu,modell_alt,i);
        fullcondp[0]->set_effect_zero();
        newmodel(kriteriumiteration2,modeliteration,textiteration);
        kriterium_aktuell = kriteriumiteration2[kriteriumiteration2.size()-1];
        }
      }
    modell_alt = modell_neu;
    }
  }


// -----------------------------------------------------------------------------
// --------------------------- Fine-Tuning -------------------------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::finetuning(vector<double> & modell)
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
                 if(int(index)>=j)
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

  vector<vector<ST::string> > names_nonp_fine;
  vector<vector<double> > lambdavec_ursp;
  vector<vector<double> > lambdavec_fine;
  vector<ST::string> namesfixed_neu;
  initialise_lambdas(names_nonp_fine,namesfixed_neu,lambdavec_ursp,number,false);
  names_fixed.erase(names_fixed.begin(),names_fixed.end());

  if(anzahl!=names_nonp_fine.size())
     {
     vector<ST::string> namesfixed_fine;
     unsigned j;
     namesfixed_fine.push_back(namesfixed_neu[0]);
     for(i=1;i<namesfixed_neu.size();i++)
        {
        bool gefunden = false;
        j = 0;
        while(j<names_nonp_fine.size() && gefunden==false)
           {
           if(names_nonp_fine[j][0]==namesfixed_neu[i])
              gefunden = true;
           else if(names_nonp_fine[j].size()>1)
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

  names_nonp.erase(names_nonp.begin(),names_nonp.end());
  names_nonp = names_nonp_fine;
  lambdavec.erase(lambdavec.begin(),lambdavec.end());
  lambdavec = lambdavec_fine;

  vector<vector<unsigned> > startindex;
  vector<vector<double> > startfix;
  startwerte("userdefined",startindex,startfix);

  ST::string text;
  bool abbruch = false;

  abbruch = single_stepwise(startindex[0],startfix[0],false);

  if(abbruch==true)
    return true;

  ST::string header = "  Final Model after Fine-Tuning:";
  fix_komplett(modell_alt);
  fullcond_komplett(modell_alt);
  ST::string tr_akt = "trace_on";
  maketext(header,modell_alt,kriterium_alt,text_alt,false,tr_akt,false);
  kriterium_tex = kriterium_alt;
  genoptions_mult[0]->out("\n\n");
  modell = modell_alt;

  return false;
  }


bool STEPWISErun::fine_local(vector<double> & modell)
  {
  smoothing = "local";
  genoptions_mult[0]->out("  BEGIN FINE-TUNING with local smoothing parameters:");
  genoptions_mult[0]->out("\n\n");
  vector<FULLCOND*> fullcond_local;
  vector<vector<double> > lambdavec_local;
  fullcond_local.push_back(fullcondp[0]);
  vector<double> modell_local;
  vector<vector<ST::string> > names_nonp_local;

  double nummer = 5;

  unsigned i,j;
  unsigned anzahl;
  for(i=1;i<fullcond_alle.size();i++)
     {
     if(modell[names_fixed.size()-2+i]!=0 && modell[names_fixed.size()-2+i]!=-1)
       {
       fullcond_alle[i]->set_smoothing("local");

       vector<double> lambdas;
       bool lambda_exist;
       unsigned index = search_lambdaindex(modell[names_fixed.size()-2+i],lambdavec[i-1],lambda_exist);
       double lambdastart = lambdavec[i-1][index];
       double lambdamin;
       double lambdamax; 

       lambda_exist = false;                       // für die Lambdas < Startwert
       int j = 2;
       while(lambda_exist==false && j>0)
          {
          if(int(index)>=j)
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
       double l = log10(lambdamin);
       double u = log10(lambdastart);
       for(j=0;j<(floor(nummer/2));j++)
         lambdas.push_back(pow(10,l+double(j)*((u-l)/(double(floor(nummer/2))-1))));

       lambda_exist = false;                        // für die Lambdas >= Startwert
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
       else
         lambdamax = 1.5*lambdavec[i-1][index] - 0.5*lambdavec[i-1][index-1];
       l = log10(lambdastart);
       u = log10(lambdamax);
       for(j=0;j<(nummer-floor(nummer/2));j++)
         lambdas.push_back(pow(10,l+double(j)*((u-l)/(double(nummer-floor(nummer/2))-1))));

       fullcond_alle[i]->set_stepwise_options(lambdastart,lambdamax,lambdamin,true,
               fullcond_alle[i]->get_df_lambdamax(),fullcond_alle[i]->get_df_lambdamin(),false,false,nummer,false);

       //anzahl = fullcond_alle[i]->get_rankK();
       anzahl = 10;  // falsch!!!, nur zum Kompilieren!
       for(j=1;j<=anzahl;j++)
          {
          lambdavec_local.push_back(lambdas);
          modell_local.push_back(lambdastart); 
          fullcond_local.push_back(fullcond_alle[i]);
          }

       vector<ST::string> names_help;
       names_help.push_back(fullcond_alle[i]->get_datanames()[fullcond_alle[i]->get_datanames().size()-1]);
       names_nonp_local.push_back(names_help);
       }
     }

  names_fixed.erase(names_fixed.begin(),names_fixed.end());
  names_fixed.push_back("const");
  names_nonp = names_nonp_local;

  fullcond_alle = fullcond_local;
  lambdavec.erase(lambdavec.begin(),lambdavec.end());
  lambdavec = lambdavec_local;

  vector<vector<unsigned> > startindex;
  vector<vector<double> > startfix;
  startwerte("userdefined",startindex,startfix);

  ST::string text;
  bool abbruch = false;

  abbruch = single_stepwise(startindex[0],startfix[0],false);

  if(abbruch==true)
    return true;

  ST::string header = "  Final Model after Fine-Tuning with local smoothing parameters:";
  fix_komplett(modell_alt);
  fullcond_komplett(modell_alt);
  ST::string tr_akt = "trace_on";
  maketext(header,modell_alt,kriterium_alt,text_alt,false,tr_akt,false);
  kriterium_tex = kriterium_alt;
  genoptions_mult[0]->out("\n\n");
  modell = modell_alt;

  return false;
  }


// -----------------------------------------------------------------------------
// ------- Funktionen für die Erstellung des Startmodels -----------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::vcm_doppelt(void)
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


void STEPWISErun::initialise_lambdas(vector<vector<ST::string> > & namen_nonp,
       vector<ST::string> & namen_fix, vector<vector<double> > & lambdavector,
       const int & number, const bool & gewichte)
  {
  namen_fix = fullcondp[0]->get_datanames();
  unsigned i;

       // berechnet ein Modell, um Gewichte für Augabe usw. zu erhalten!!!
  if(gewichte == true)
    {
    vector<double> modell_init;
    for(i=1;i<namen_fix.size();i++)
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

    posteriormode(posttitle,true);
    }

  for(i=1;i<fullcond_alle.size();i++)
     {
     int nummer = number;
     if(fullcond_alle[i]->get_data_forfixedeffects().cols() > 1)  //bei Faktor-Variablen
        namen_nonp.push_back(fullcond_alle[i]->get_datanames());
     else
        {
        vector<ST::string> names_help;
        names_help.push_back(fullcond_alle[i]->get_datanames()[fullcond_alle[i]->get_datanames().size()-1]);
        namen_nonp.push_back(names_help);
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
     lambdavector.push_back(untervector);
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
  unsigned index = search_lambdaindex(start,lambdas,gefunden);

  if(gefunden==false)
    {
    vector<double> diff;
    unsigned i;
    for(i=0;i<lambdas.size();i++)
       {
       if(lambdas[i]!=0 && lambdas[i]!=-1)
         diff.push_back(fabs(lambdas[i]-start));
       }
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
          vector<vector<unsigned> > & startindex,
          vector<vector<double> > & startfix)
  {
  unsigned i;

  if(startmodel == "empty" || startmodel == "both" || startmodel == "emplin")
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

  if(startmodel == "emplin")
    {
    vector<unsigned> indexhelp;
    vector<double> fixhelp;

    for(i=1;i<names_fixed.size();i++)
       fixhelp.push_back(-1);                      
    for(i=1;i<fullcond_alle.size();i++)
       {
       unsigned index = search_lambdastartindex(-1,lambdavec[i-1]);
       indexhelp.push_back(index);
       }

    startindex.push_back(indexhelp);
    startfix.push_back(fixhelp);
    }

  if(startmodel == "userdefined")
    {
    vector<unsigned> indexhelp;
    vector<double> fixhelp;

    for(i=1;i<names_fixed.size();i++)
       fixhelp.push_back(-1);         //fehlt: vom Benutzer angeben lassen!!!
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

  /*if(df>=499)      // Versuch
    df = 498;*/

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

// Versuch, ob penalisierte statt normaler Log-Likelihood verwendet werden sollte
/*
unsigned j;
double sum = 0;
for(j=1;j<fullcondp.size();j++)
  {
  sum += fullcondp[j]->compute_penalterm();
  }
sum /= likep_mult[0]->compute_rss();
sum *= likep_mult[0]->get_nrobs_wpw();
kriterium += sum;
*/

  return kriterium;
  }


void STEPWISErun::newmodel(vector<double> & krit,
  vector<vector<double> > & mi, vector<ST::string> & textit)
  {
  fertig = false;
  mi.push_back(modell_neu);
  posteriormode(posttitle,true);      //keine Ausgabe
  double kriterium = compute_criterion();
  ST::string header = "  Trial: ";
  ST::string text;
  maketext(header,modell_neu,kriterium,text,true,trace,false);
  textit.push_back(text);
  krit.push_back(kriterium);
  }


void STEPWISErun::newmodel_fix(const double & mo, vector<double> & krit,
   vector<vector<double> > & mi, vector<ST::string> & textit,
   const ST::string & name)
  {
  if(mo==0)
    reset_fix(name);
  else
    include_fix(name);
  fullcondp[0]->set_effect_zero();
  newmodel(krit,mi,textit);
  if(mo==0)
    include_fix(name);
  else
    reset_fix(name);
  //fullcondp[0]->set_effect_zero();
  }


void STEPWISErun::newmodel_factor(const double & mo, const unsigned & index,
    vector<double> & krit, vector<vector<double> > & mi,
    vector<ST::string> & textit, const vector<ST::string> & name)
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
  newmodel(krit,mi,textit);
  if(mo==0)
    fullcondp[0]->include_effect(name,fullcond_alle[index]->get_data_forfixedeffects());
  else
    {
    for(i=0;i<name.size();i++)
      reset_fix(name[i]);
    }
  //fullcondp[0]->set_effect_zero();
  }


void STEPWISErun::newmodel_nonp(const unsigned & index,
    vector<double> & krit, vector<vector<double> > & mi,
    vector<ST::string> & textit)
  {
  fullcond_einzeln(modell_neu,modell_alt,index);
  fullcondp[0]->set_effect_zero();
  newmodel(krit,mi,textit);
  fullcond_einzeln(modell_alt,modell_neu,index);
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

void STEPWISErun::fullcond_einzeln(const vector<double> & modell1,
         const vector<double> & modell2, const unsigned & index)
  {

if(smoothing == "global")    // Versuch!!!
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
     fullcond_alle[i]->set_inthemodel(modell1[names_fixed.size()-2+i]);
     if(modell2[names_fixed.size()-2+i]==-1 && index==i)
        reset_fix(names_nonp[i-1][0]);
     if(modell1[names_fixed.size()-2+i]>0)
        {
        fullcond_neu.push_back(fullcond_alle[i]);
        if(i == index)
          fullcond_alle[i]->update_stepwise(modell1[names_fixed.size()-2+i]);
        }
     else if(modell1[names_fixed.size()-2+i]==0)    // && i == index)
        fullcond_alle[i]->reset_effect(0);
     else if(modell1[names_fixed.size()-2+i] == -1) // && i == index)
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
else //if(smoothing == "local")
  fullcond_alle[index]->update_stepwise(modell1[names_fixed.size()-2+index]);

  }


void STEPWISErun::fullcond_komplett(const vector<double> & m)
  {

if(smoothing == "global")
  {
  vector<FULLCOND*> fullcond_neu;
  unsigned i;

  fullcond_neu.push_back(fullcondp[0]);
  for(i=1;i<fullcond_alle.size();i++)
     {
     fullcond_alle[i]->set_inthemodel(m[names_fixed.size()-2+i]);
     if(m[names_fixed.size()-2+i]>0)
        {
        fullcond_alle[i]->update_stepwise(m[names_fixed.size()-2+i]);
        fullcond_neu.push_back(fullcond_alle[i]);
        }
     else if(m[names_fixed.size()-2+i]==0)
        fullcond_alle[i]->reset_effect(0);
     else if(m[names_fixed.size()-2+i] == -1)
        {
        fullcond_alle[i]->reset_effect(0);
        fullcond_neu[0]->include_effect(names_nonp[i-1],
                                  fullcond_alle[i]->get_data_forfixedeffects());
        }
     }

  fullcondp = fullcond_neu;
  end[0] = fullcondp.size()-1;
  if(minim!="adaptiv")   // stürzt sonst bei "binomial" beim letzten "posteriormode" ab, warum???
    fullcondp[0]->set_effect_zero();
  }
else //if(smoothing == "local")
  {
  unsigned i;
  for(i=1;i<fullcond_alle.size();i++)
     fullcond_alle[i]->update_stepwise(m[names_fixed.size()-2+i]);
  }

  }
  

void STEPWISErun::fix_komplett(const vector<double> &  modell)
  {

if(smoothing == "global")
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
          include_fix(names_fixed[z+1]);
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
  }  // END: if(smoothing == "global")
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


void STEPWISErun::include_fix(const ST::string & name)
  {
  int i = column_for_fix(name);
  vector<ST::string> help_name;
  help_name.push_back(name);
  fullcondp[0]->include_effect(help_name,datamatrix(D.getCol(i)));
  }

int STEPWISErun::column_for_fix(const ST::string & name)
  {
  bool gefunden = false;
  unsigned i = 0;
  while(i<modelv.size() && gefunden==false)
     {
     if(name==modelv[i])
        gefunden = true;
     i = i + 1;
     }
  return i-1;
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
  if(tr == "trace_on" || trace == "trace_minim")
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
  if(tr == "trace_on" || trace == "trace_minim")
    {
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out(modeltext);
    genoptions_mult[0]->out("\n " + criterion + " = " + ST::doubletostring(a,8));
    }
  if(datei==true)
    outmodels << modeltext << endl;
  }


void STEPWISErun::options_text(const int & number,
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
                  vector<vector<unsigned> > & startindex)
  {
  ST::string title = "STEPWISEREG OBJECT " + name + ": " + algorithm + " procedure";
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

  make_prior(startindex);

  outtex << "\n\\noindent {\\bf \\large Start Predictor";
  if(startindex.size()>1)
     outtex << "s";
  outtex << ":}\\\\" << endl;

  }


void STEPWISErun::make_tex_end(ST::string & path, const vector<double> & modell)
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
     if(modell[names_fixed.size()-2+j]!=0 && modell[names_fixed.size()-2+j]!=-1)
         {
         vector<ST::string> prior = fullcond_alle[j]->get_priorassumptions();
         outtex << prior[0] << "\\\\" << endl
                << "smoothing parameter: $\\lambda = "
                << ST::doubletostring(modell[names_fixed.size()-2+j],6)
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
  char hcharu = '_';
  ST::string hstringu = "\\_";

  //schreibt Stepwise options ins Tex-File
  outtex << "\n\\noindent {\\bf \\large Stepwise Options:}" << endl
         << "\\begin{tabbing}" << endl
         //<< "Levels for credible intervals: \\= \\\\" << endl
         //<< "Level 1: \\> " << l1 << "\\\\" << endl
         //<< "Level 2: \\> " << l2 << "\\\\" << endl
         << "Maximum number of Iterations: \\= " << steps << " \\\\" << endl
         << "Performance criterion: \\> " << criterion.insert_string_char(hcharu,hstringu) << " \\\\" << endl
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


void STEPWISErun::make_prior(vector<vector<unsigned> > & startindex)
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
     outtex << endl << "\\begin{tabular}{p{12cm}}\n" << term
            << "$\n\\end{tabular}\n" << endl;
     if(startmodel == "full" || startmodel == "userdefined")
        outtex << "Startvalue is the fixed effect \\\\ \n\\\\" << endl;
     else if(startmodel == "empty")
        outtex << "Startvalue is 'effect excluded' \\\\ \n\\\\" << endl;
     else
        outtex << "1. Startvalue is 'effect excluded' \\\\" << endl
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
           outtex << "Startvalue is 'effect excluded' \\\\ \n";
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

// -----------------------------------------------------------------------------
// ------------- Golden Section Search -----------------------------------------
// -----------------------------------------------------------------------------

unsigned STEPWISErun::golden_section(unsigned & z, double & kriterium)
  {
  // Algorithmus macht nur Sinn, wenn genug verschiedene Lambdas da sind! Wieviele? >=7

  vector<unsigned> index_vec;
  vector<double> krit_vec;
  unsigned x0, b, x3, x1, x2;
  double f0, fb, f3, f1, f2;
  double fraus = 100000;
  double xraus = -1;             
  double gold = (3-sqrt(5))/2;
  unsigned i;

  double df = startbedingungen(z,kriterium);

  bool falsch = false;
  unsigned index = search_lambdaindex(modell_alt[z+names_fixed.size()-2],
                                lambdavec[z-1],falsch);

  // für FG, die nicht durch ein lambda gebildet werden können:
  int grenzfall = fullcond_alle[z]->get_grenzfall();
  if(grenzfall == 0)
    index_vec.push_back(lambdavec[z-1].size()-1);
  else if(grenzfall == 1 && fullcond_alle[z]->get_forced()==false)
    {
    index_vec.push_back(lambdavec[z-1].size()-2);
    if(index != lambdavec[z-1].size()-1)
      fraus = wert_einzeln(z,lambdavec[z-1].size()-1,df);
    else
      fraus = kriterium;
    xraus = lambdavec[z-1].size()-1;
    }
  else if(grenzfall == 1 && fullcond_alle[z]->get_forced()==true)
    index_vec.push_back(lambdavec[z-1].size()-1);
  else if(grenzfall >= 2 && fullcond_alle[z]->get_forced()==false)
    {
    if(lambdavec[z-1].size()-2 == -1)
      {
      index_vec.push_back(lambdavec[z-1].size()-3);
      if(index != lambdavec[z-1].size()-1)
        fraus = wert_einzeln(z,lambdavec[z-1].size()-1,df);
      else
        fraus = kriterium;
      double fraus2;
      if(index != lambdavec[z-1].size()-2)
        fraus2 = wert_einzeln(z,lambdavec[z-1].size()-2,df);
      else
        fraus2 = kriterium;
      if(fraus2 < fraus)
        xraus = lambdavec[z-1].size()-2;
      else
        xraus = lambdavec[z-1].size()-1;
      }
    else      // trifft bei "season" zu!
      {
      index_vec.push_back(lambdavec[z-1].size()-2);
      if(index != lambdavec[z-1].size()-1)
        fraus = wert_einzeln(z,lambdavec[z-1].size()-1,df);
      else
        fraus = kriterium;
      xraus = lambdavec[z-1].size()-1;
      }
    }
  else if(grenzfall >= 2 && fullcond_alle[z]->get_forced()==true)
    {
    if(lambdavec[z-1].size()-1 == -1)
      {
      index_vec.push_back(lambdavec[z-1].size()-2);
      if(index != lambdavec[z-1].size()-1)
        fraus = wert_einzeln(z,lambdavec[z-1].size()-1,df);
      else
        fraus = kriterium;
      xraus = lambdavec[z-1].size()-1;
      }
    else   // kann bei "season" zutreffen!
      index_vec.push_back(lambdavec[z-1].size()-1);
    }

  /*if(fullcond_alle[z]->get_forced()==false)
    {
    index_vec.push_back(lambdavec[z-1].size()-2);
    if(index != lambdavec[z-1].size()-1)
      fraus = wert_einzeln(z,lambdavec[z-1].size()-1,df);
    else
      fraus = kriterium;
    }
  else
    index_vec.push_back(lambdavec[z-1].size()-1);   */

  index_vec.push_back(0);
  x0 = index_vec[0];
  x3 = index_vec[1];

  if(index!=x0)
    f0 = wert_einzeln(z,x0,df);
  else
    f0 = kriterium;
  krit_vec.push_back(f0);
  if(index!=x3)
    f3 = wert_einzeln(z,x3,df);
  else
    f3 = kriterium;
  krit_vec.push_back(f3);
  index_vec.push_back(index);
  krit_vec.push_back(kriterium);

  if(index!=x0 && index!=x3 && index!=lambdavec[z-1].size()-1)
    b = index;
  else
    b = x0 - floor((x0-x3)/2*gold + 0.5);
  int indexb = index_suchen(b,index_vec);
  if(indexb < 0)
    fb = wert_einzeln(z,b,df);
  else
    fb = krit_vec[indexb];
  index_vec.push_back(b);
  krit_vec.push_back(fb);

  b = start_b_suchen(index_vec,krit_vec,df,b,z);     // sucht passenden Startwert für b!
  if(b >= lambdavec[z-1].size())
     {
     unsigned xmin;
     if(f3 < f0 && f3 < fraus && f3 < kriterium)
       xmin = x3;
     else if(f0 <= f3 && f0 < fraus && f0 < kriterium)
       xmin = x0;
     else if(fraus < kriterium && fraus <= f3 && fraus <= f0)
       xmin = lambdavec[z-1].size()-1;
     else
       xmin = index;

     if(trace == "trace_minim")
       {
       genoptions_mult[0]->out("\n\n");
       genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "  Testlambdas:  " + " Testvalues: \n");
       for(i=0;i<krit_vec.size();i++)
          genoptions_mult[0]->out("      " + ST::doubletostring(lambdavec[z-1][index_vec[i]],6)
                               + "    " + ST::doubletostring(krit_vec[i],6) + "\n");
       genoptions_mult[0]->out("\n");
       }

     if(minim != "exact_golden")
       approx_zurueck(z);
     else
       exact_zurueck(z);
     return xmin;
     }
  else
     {
     if(abs(x3-b) > abs(b-x0))
        {
        x1 = b;
        f1 = fb;
        x2 = floor(b - gold*(b-x3) + 0.5);
        f2 = compute_findex(index_vec,krit_vec,x2,z,df);
        }
     else // if(abs(x3-b) <= abs(b-x0))
        {
        x2 = b;
        f2 = fb;
        x1 = floor(b + gold*(x0-b) + 0.5);
        f1 = compute_findex(index_vec,krit_vec,x1,z,df);
        }

     while(abs(int(x3-x0)) > 1)
        {
        if(f2 < f1)  // && x1 != x2)
          {
          x0 = x1;
          f0 = f1;
          x1 = x2;
          f1 = f2;
          x2 = floor((1-gold)*x1 + gold*x3 + 0.5);
          if(x1 != x2)
            f2 = compute_findex(index_vec,krit_vec,x2,z,df);
          else
            {
            if(x3<x2-1)
              {
              x2 = x2 - 1;
              f2 = compute_findex(index_vec,krit_vec,x2,z,df);
              }
            else
              f2 = f3;
            }
          }
        else // if (f2 >= f1 && x1 != x2)
          {
          x3 = x2;
          f3 = f2;
          x2 = x1;
          f2 = f1;
          x1 = floor((1-gold)*x2 + gold*x0 + 0.5);
          if(x1 != x2)
            f1 = compute_findex(index_vec,krit_vec,x1,z,df);
          else
            {
            if(x0>x1+1)
              {
              x1 = x1 + 1;
              f1 = compute_findex(index_vec,krit_vec,x1,z,df);
              }
            else
              f1 = f0;
            }
          }
        /*
        else if (f2 < f1 && x1 == x2)
          {
          x0 = x1;
          x1 = x2;
          x2 = floor((1-gold)*x1 + gold*x3 + 0.5);
          f1 = f2;
          f2 = compute_findex(index_vec,krit_vec,x2,z,df);
          }
        else // if (f2 >= f1 && x1 == x2)  //Problem: Wenn x2 = x1, dann
          {
          x3 = x2;
          x2 = x1;
          x1 = floor((1-gold)*x2 + gold*x0 + 0.5);
          f2 = f1;
          f1 = compute_findex(index_vec,krit_vec,x1,z,df);
          }
        */
        }

     unsigned xmin;
     if(f1 <= f2 && f1 < fraus && f1 < kriterium)
       {
       xmin = x1;
       kriterium = f1;
       }
     else if(f2 < f1 && f2 < fraus && f2 < kriterium)
       {
       xmin = x2;
       kriterium = f2;
       }
     else if(fraus < kriterium && fraus <= f1 && fraus <= f2)
       {
       xmin = xraus;
       kriterium = fraus;
       }
     else
       xmin = index;

     if(trace == "trace_minim")
       {
       genoptions_mult[0]->out("\n\n");
       genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "  Testlambdas:  " + " Testvalues: \n");
       for(i=0;i<krit_vec.size();i++)
          genoptions_mult[0]->out("      " + ST::doubletostring(lambdavec[z-1][index_vec[i]],6)
                               + "    " + ST::doubletostring(krit_vec[i],6) + "\n");
       genoptions_mult[0]->out("\n");
       }
     
     if(minim != "exact_golden")
       approx_zurueck(z);
     else
       exact_zurueck(z);
     return xmin;
     }
  }

double STEPWISErun::startbedingungen(unsigned & z, double & kriterium)
  {
  unsigned i;
  double df = 0;
  for(i=0;i<fullcondp.size();i++)               // Startbedingungen herstellen: FGe ohne Variable, Fullcondp-Vektor anpassen, ...
    df = df + fullcondp[i]->compute_df();
  if(minim != "exact_golden")
    fullcond_alle[0]->safe_const();
  if(modell_alt[z+names_fixed.size()-2]!=-1 && modell_alt[z+names_fixed.size()-2]!=0)
    {
    if(minim == "adaptiv_golden")
      {
      fullcond_alle[z]->update_stepwise(modell_alt[z+names_fixed.size()-2]);
      fullcond_alle[z]->posteriormode();
      kriterium = criterion_min(df);
      }
    df = df - fullcond_alle[z]->compute_df();
    fullcond_alle[z]->reset_effect(0);
    }
  else if(modell_alt[z+names_fixed.size()-2] == -1)
    {
    if(minim == "adaptiv_golden")
      {
      reset_fix(names_nonp[z-1][0]);
      fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                  fullcond_alle[z]->get_data_forfixedeffects());
      kriterium = criterion_min(df);
      }
    df = df - 1;
    fullcondp.push_back(fullcond_alle[z]);
    reset_fix(names_nonp[z-1][0]);
    }
  else
    fullcondp.push_back(fullcond_alle[z]);

  end[0] = fullcondp.size()-1;
  return df;
  }


void STEPWISErun::approx_zurueck(unsigned & z)
  {
  if(modell_alt[z+names_fixed.size()-2]!=-1 && modell_alt[z+names_fixed.size()-2]!=0)
    {
    fullcond_alle[z]->set_inthemodel(1);
    fullcond_alle[z]->update_stepwise(modell_alt[z+names_fixed.size()-2]);
    fullcond_alle[z]->posteriormode();
    }
  else if(modell_alt[z+names_fixed.size()-2] == -1)
    {
    fullcond_alle[z]->set_inthemodel(0);
    fullcond_alle[z]->reset_effect(0);
    fullcondp.erase(fullcondp.end()-1,fullcondp.end());
    fullcondp[0]->posteriormode_single(names_nonp[z-1],
                                fullcond_alle[z]->get_data_forfixedeffects());
    }
  else
    {
    fullcond_alle[z]->set_inthemodel(0);
    fullcond_alle[z]->reset_effect(0);
    fullcondp.erase(fullcondp.end()-1,fullcondp.end());
    }
  end[0] = fullcondp.size()-1;
  }


void STEPWISErun::exact_zurueck(unsigned & z)
  {
  if(modell_alt[z+names_fixed.size()-2]!=-1 && modell_alt[z+names_fixed.size()-2]!=0)
    {
    fullcond_alle[z]->update_stepwise(modell_alt[z+names_fixed.size()-2]);
    fullcond_alle[z]->set_inthemodel(1);
    }
  else if(modell_alt[z+names_fixed.size()-2] == -1)
    {
    fullcond_alle[z]->set_inthemodel(0);
    fullcond_alle[z]->reset_effect(0);
    fullcondp[0]->include_effect(names_nonp[z-1],
                                fullcond_alle[z]->get_data_forfixedeffects());
    fullcondp.erase(fullcondp.end()-1,fullcondp.end());
    }
  else
    {
    fullcond_alle[z]->set_inthemodel(0);
    fullcond_alle[z]->reset_effect(0);
    fullcondp.erase(fullcondp.end()-1,fullcondp.end());
    }
  end[0] = fullcondp.size()-1;
  posteriormode(posttitle,true);
  }


double STEPWISErun::wert_einzeln(unsigned & z, unsigned i, double & df)
  {
  double kriterium;
  if(minim == "exact_golden")
    kriterium = exact_einzeln(z,i,df);
  else    //if(minim == "approx_golden" || minim == "adaptiv_golden)
    kriterium = approx_einzeln(z,i,df);
  return kriterium;
  }

double STEPWISErun::approx_einzeln(unsigned & z, unsigned & i, double & df)
  {
  double kriterium_versuch;
  if(lambdavec[z-1][i]!=-1 && lambdavec[z-1][i]!=0)
    {
    fullcond_alle[z]->set_inthemodel(1);
    fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
    fullcond_alle[z]->posteriormode();
    kriterium_versuch = criterion_min(df + fullcond_alle[z]->compute_df());
    fullcond_alle[0]->set_const_old();
    fullcond_alle[z]->reset_effect(0);
    }
  else if(lambdavec[z-1][i]==-1)
    {
    fullcond_alle[z]->set_inthemodel(0);
    fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                  fullcond_alle[z]->get_data_forfixedeffects());
    kriterium_versuch = criterion_min(df + 1);
    fullcond_alle[0]->set_const_old();
    reset_fix(names_nonp[z-1][0]);
    }
  else
    {
    fullcond_alle[z]->set_inthemodel(0);
    fullcond_alle[0]->posteriormode_const();
    kriterium_versuch = criterion_min(df);
    fullcond_alle[0]->set_const_old();
    }
  return kriterium_versuch;
  }


double STEPWISErun::exact_einzeln(unsigned & z, unsigned & i, double & df)
  {
  double kriterium_versuch;
  if(lambdavec[z-1][i]!=-1 && lambdavec[z-1][i]!=0)
    {
    fullcond_alle[z]->set_inthemodel(1);
    fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
    posteriormode(posttitle,true);
    kriterium_versuch = criterion_min(df + fullcond_alle[z]->compute_df());
    }
  else if(lambdavec[z-1][i]==-1)
    {
    fullcond_alle[z]->set_inthemodel(0);
    vector<FULLCOND*> fullcond_start = fullcondp;
    vector<double> modell1 = modell_alt;
    modell1[z+names_fixed.size()-2] = -1;
    fullcond_einzeln(modell1,modell_alt,z);  // hier muß der Fullcond-Vekor angepaßt werden!!!
    posteriormode(posttitle,true);
    kriterium_versuch = criterion_min(df + 1);
    fullcondp = fullcond_start;
    end[0] = fullcondp.size()-1;
    reset_fix(names_nonp[z-1][0]);
    }
  else
    {
    fullcond_alle[z]->set_inthemodel(0);
    vector<FULLCOND*> fullcond_start = fullcondp;
    vector<double> modell1 = modell_alt;
    modell1[z+names_fixed.size()-2] = 0;
    fullcond_einzeln(modell1,modell_alt,z);  // hier muß der Fullcond-Vekor angepaßt werden!!!
    posteriormode(posttitle,true);
    kriterium_versuch = criterion_min(df);
    fullcondp = fullcond_start;
    end[0] = fullcondp.size()-1;
    }
  return kriterium_versuch;
  }


int STEPWISErun::index_suchen(const unsigned & index, const vector<unsigned> & index_vec)
  {
  bool s = false;
  int x = index_vec.size()-1;
  while(x>=0 && s==false)
       {
       if(index_vec[x]==index)
          s = true;
       x = x - 1;
       }
  if(s==true)
    return x+1;
  else
    return x;
  }


int STEPWISErun::start_b_suchen(vector<unsigned> & index_vec, vector<double> & krit_vec,
                  double & df, unsigned & b, unsigned & z)
  {
  double fb = krit_vec[krit_vec.size()-1];
  // unsigned b_alt = b;
  int b_neu = index_vec[0]-1;
  bool sofort = true;
  bool found = true;
  while((krit_vec[0] < fb || krit_vec[1] < fb) && found == true)
     {
     sofort = false;
     found = false;
     while(b_neu > 0)
       {
       found = true;
       if(b_neu != int(b))
         {
         if(b_neu != int(index_vec[2]))
           {
           fb = wert_einzeln(z,unsigned(b_neu),df);
           index_vec.push_back(b_neu);
           krit_vec.push_back(fb);
           }
         else
           fb = krit_vec[2];
         }
       b_neu = b_neu - 2;
       }
     }

  if(sofort == false)
    b = b_neu + 2;

  if(found == true)
    {
    return b;
    }
  else
    return lambdavec[z-1].size()+1;
  }

double STEPWISErun::compute_findex(vector<unsigned> & index_vec, vector<double> & krit_vec,
                  unsigned & index, unsigned & z, double & df)
  {
  int indexb = index_suchen(index,index_vec);
  double f;
  if(indexb < 0)
    {
    f = wert_einzeln(z,index,df);
    index_vec.push_back(index);
    krit_vec.push_back(f);
    }
  else
    f = krit_vec[indexb];
  return f;
  }


// -----------------------------------------------------------------------------
// ------------- Model Averaging -----------------------------------------------
// -----------------------------------------------------------------------------

void STEPWISErun::compute_average(void)
  {
  vector<double> kriterien_alle;
  vector<vector<double> > modelle;
  vector<ST::string> ausgabe;

  unsigned i;
  unsigned j;
  for(i=0;i<modell_neu.size();i++)    // fixe Effekte im Endmodell u.U. in falscher Reihenfolge!
    modell_neu[i] = 0;                // deswegen: zuerts Modell leer machen,
  fix_komplett(modell_neu);
  fullcond_komplett(modell_neu);

  fix_komplett(modell_alt);          // dann Endmodell mit richtiger Reihenfolge erzeugen!
  fullcond_komplett(modell_alt);
  posteriormode(posttitle,true);
  double kriterium = compute_criterion();
  genoptions_mult[0]->out("\n\n\n");  
  genoptions_mult[0]->out("  Beginning of the Model_Averaging! \n");
  ST::string header = "  Final Model: ";
  ST::string text;
  maketext(header,modell_alt,kriterium,text,true,trace,false);

  kriterien_alle.push_back(kriterium);
  save_alle_betas(modell_alt);

  for(i=1;i<names_fixed.size();i++)
    {
    modell_neu = modell_alt;
    if(modell_alt[i-1]==-1)
      {
      modell_neu[i-1]= 0;
      fix_komplett(modell_neu);
      fullcond_komplett(modell_neu);
      }
    else if(modell_alt[i-1]==0)
      {
      modell_neu[i-1] = -1;
      fix_komplett(modell_neu);
      fullcond_komplett(modell_neu);
      }
    fullcondp[0]->set_effect_zero();
    newmodel(kriterien_alle,modelle,ausgabe);

    save_alle_betas(modell_neu);
    }

  unsigned z = 1;
  while(z<fullcond_alle.size() && fullcond_alle[z]->get_fctype()==factor)
     {
     modell_neu = modell_alt;
     if(modell_alt[z+names_fixed.size()-2]==-1 && fullcond_alle[z]->get_forced()==false)
       {
       modell_neu[z+names_fixed.size()-2]= 0;
       fix_komplett(modell_neu);
       fullcond_komplett(modell_neu);
       }
     else if(modell_alt[z+names_fixed.size()-2]==0)
       {
       modell_neu[z+names_fixed.size()-2] = -1;
       fix_komplett(modell_neu);
       fullcond_komplett(modell_neu);
       }
     fullcondp[0]->set_effect_zero();
     newmodel(kriterien_alle,modelle,ausgabe);

     save_alle_betas(modell_neu);
     z = z + 1;
     }

  for(i=z;i<fullcond_alle.size();i++)
    {
    int sch;
    for(sch=1;sch<=window;sch++)
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
          fix_komplett(modell_neu);
          fullcond_komplett(modell_neu);
          fullcondp[0]->set_effect_zero();
          newmodel(kriterien_alle,modelle,ausgabe);

          save_alle_betas(modell_neu);
          }

       lambda_exist = false;
       modell_neu = modell_alt;
       if(int(index) >= sch)
          lambda_exist = true;
       if(lambda_exist==true)
          {
          modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][index-sch];
          fix_komplett(modell_neu);
          fullcond_komplett(modell_neu);
          fullcondp[0]->set_effect_zero();
          newmodel(kriterien_alle,modelle,ausgabe);

          save_alle_betas(modell_neu);
          }
       }
    }

  int l;         // Gewichte berechnen!
  double kriterium_summe = 0;
  for(l=kriterien_alle.size()-1;l>=0;l--)
     {
     kriterien_alle[l] = exp(-0.5*(kriterien_alle[l]-kriterien_alle[0]));   // es ist egal, ob kriterien_alle[0] wirklich das Minimum ist!!!
     kriterium_summe += kriterien_alle[l];
     }
  for(i=0;i<kriterien_alle.size();i++)
     kriterien_alle[i] /= kriterium_summe;
/*
vector<double> versuch;
versuch.push_back(0.15);
versuch.push_back(0.15);
versuch.push_back(0.15);
versuch.push_back(0.15);
versuch.push_back(0.1);
versuch.push_back(0.3);
*/
  for(i=0;i<names_fixed.size()-1;i++)    // sorgt dafür, dass alle fixen Effekte und fullcond-Vektor stimmen,
    modell_neu[i] = -1;
  for(j=i;j<modell_neu.size();j++)       // d.h. alle fullcondp = fullcond_alle und bei den fixen Effekten nur noch aus fullcondp[0]!
    {
    if(fullcond_alle[j-names_fixed.size()+2]->get_fctype()==factor)
      modell_neu[j] = -1;
    else
      modell_neu[j] = 1;
    }
  fix_komplett(modell_neu);
  fullcond_komplett(modell_neu);
  for(i=1;i<fullcondp.size();i++)
    fullcondp[i]->average_posteriormode(kriterien_alle);
    //fullcond_alle[i]->average_posteriormode(versuch);
  fullcond_alle[0]->average_posteriormode(kriterien_alle);
  //fullcond_alle[0]->average_posteriormode(versuch);

  likep_mult[0]->posteriormode();       // berechnet "sigma2"
  genoptions_mult[0]->out("\n\n");
  likep_mult[0]->outresults();          // Ausgabe von "predictmean.raw", mit nicht standardidiertem Prädiktor
  for(i=0;i<fullcondp.size();i++)
    fullcondp[i]->outresults();     // Ausgabe der .res-Dateien
  }


void STEPWISErun::save_alle_betas(vector<double> & modell)
  {
  unsigned anzahl = names_fixed.size();
  fullcond_alle[0]->save_betas(modell,anzahl);
  anzahl = 1;
  unsigned j;
  for(j=0;j<names_fixed.size()-1;j++)
    {
    if(modell[j] == -1)
      anzahl++;
    }
  for(j=1;j<fullcond_alle.size();j++)
    {
    if(modell[j+names_fixed.size()-2] == -1)
      {
      fullcond_alle[j]->save_betas(modell,anzahl);
      anzahl += fullcond_alle[j]->get_data_forfixedeffects().cols();
      }
    else if(modell[j+names_fixed.size()-2] == 0)
      fullcond_alle[j]->save_betas(modell,0);
    else
      fullcond_alle[j]->save_betas(modell,-1);   // ACHTUNG: hier steht -1 für nichtlinear!!!
    }
  fullcond_alle[0]->save_betas2();
  }



}



