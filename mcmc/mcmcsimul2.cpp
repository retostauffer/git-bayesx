
#include "first.h"

#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include<StatwinFrame.h>
#endif

#include"mcmcsimul2.h"
#include<time.h>
#include"clstring.h"
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
        const bool & fineloc, const bool & maveraging, int & fenster,
        const datamatrix & Da, const vector<ST::string> & modelvar,
        const ST::string & name, vector<FULLCOND*> & fullcond_z, ST::string & path,
        const bool & CI, bool & hier, bool & gm, const double & prop, const bool & minib)
  {

  D = Da;
  modelv = modelvar;
  algorithm = procedure;
  minim = minimum;
  miniback_off = minib;
  criterion = crit;
  increment = inc;
  steps = stp;
  startmodel = stam;
  fine_tuning = finet;
  finelocal = fineloc;
  trace = trac;
  bool modelaveraging = maveraging;
  window = fenster;
  smoothing = "global";
  hierarchical = hier;
  ganze_matrix = gm;
  df_exact = false;

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

  set_center(likep_mult[0],fullcond_alle,begin[0],end[0]);  // sorgt daf�r, da� Funktionen zentriert werden!

  bool gewichte = false;
  if(likep_mult[0]->get_family() != "Gaussian")
    gewichte = true;
  initialise_lambdas(names_nonp,names_fixed,lambdavec,number,gewichte);

  if(criterion == "MSEP" || criterion == "AUC" ||criterion == "CV5" || criterion == "CV10")
    initialise_weights(prop);

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
  // �berpr�fen, dass Randomslopes nicht auch als fixe Effekte angegeben werden!
  if( vcm_doppelt() == true)      // F�r VCM-Modelle!!! --> mu� die Funktion raus?
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
  outmodels << "step   " << criterion << "   model" << endl << endl;

  for(i=0;i<startindex.size();i++)
     {
     abbruch = single_stepwise(startindex[i],startfix[i],true);

     if(abbruch==true)
        return true;

     if(minim=="adaptiv" || minim=="adaptiv_golden")
       {
       schaetzen(0,kriterium_alt,true,"backfitting");
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
  if(finelocal == true)
     abbruch = fine_local(modell_final);

  if(abbruch==true)
    return true;

  if(modelaveraging == true)
    {
    ST::string hilf = path + "_weights.raw";
    ofstream outweight(hilf.strtochar());
    fullcond_z = fullcond_alle;
    modell_alt = modell_final;
    compute_average(outweight);
    fullcond_z = fullcondp;
    for(i=0;i<fullcond_z.size();i++)
      fullcond_z[i]->set_fcnumber(i);
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

    if(CI == false || finelocal == true)  // bei lokalem Gl�ttungsparameter erst mal keine Konfidenzintervalle!!!
      {
      if(criterion == "MSEP" || criterion == "AUC")
        likep_mult[0]->weight_for_all();
      posteriormode(posttitle,false); // Problem: linearer Pr�diktor bei "true" standardisiert! Hier wird zur�ckgerechnet!
                                      // danach nicht mehr compute_criterion() aufrufen!!!

/*double df = df_ganzehatmatrix();
df = 0;
for(unsigned k=0;k<fullcondp.size();k++)
  df += fullcondp[k]->compute_df();
df = df; */
      }
    else
      {
  #if defined(BORLAND_OUTPUT_WINDOW)
      genoptions_mult[0]->out("  Calculating of confidence intervals by MCMC-techniques is not possible with this version!");
      posteriormode(posttitle,false);
  #elif defined(JAVA_OUTPUT_WINDOW)
      ST::string befehl = name + ".regress " + likep_mult[0]->get_responsename() + " = ";
      ST::string textfix = fullcond_alle[0]->get_befehl();
      if(textfix != "" && fullcondp.size()>1)
        befehl = befehl + textfix + " + ";
      for(i=1;i<fullcondp.size()-1;i++)
        befehl = befehl + fullcondp[i]->get_befehl() + " + ";
      befehl = befehl + fullcondp[fullcondp.size()-1]->get_befehl();
      befehl = befehl + ", constlambda ";
      ST::string family;
      if(likep_mult[0]->get_family() == "Gaussian")
        {
        double scale = likep_mult[0]->get_scale();
        double transform = likep_mult[0]->get_trmult(0);
        scale *= transform*transform;
        befehl = befehl + "constscale scale=" + ST::doubletostring(scale,6);
        family = "gaussian";
        }
      else if(likep_mult[0]->get_family() == "Binomial (logit link)")
        family = "binomial";
      else if(likep_mult[0]->get_family() == "Binomial (probit link)")
        family = "binomialprobit";
      else if(likep_mult[0]->get_family() == "Poisson")
        family = "poisson";
      befehl = befehl + " burnin=2000 step=10 iterations=12000 predict family=" + family;
      path = befehl;
  #endif
      }
    }

  make_tex_end(path,modell_final);

  // Files m�ssen wieder geschlossen werden!!!
  outtex.close();
  outcriterium.close();
  outmodels.close();

  // gibt Lambdas aus, damit man die richtig bestimmten Variablen z�hlen kann!
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

  if(hierarchical == true)
    {
    for(i=fullcond_alle.size()-1;i>=1;i--)   // Abfrage, ob Startmodell hierarchisch ist!
       {
       ST::string possible = "alles";
       fullcond_alle[i]->hierarchical(possible);
       bool falsch = true;

       if((modell_alt[names_fixed.size()-2+i] == -1 || modell_alt[names_fixed.size()-2+i] > 0) && possible == "raus")
         falsch = false;
       if(modell_alt[names_fixed.size()-2+i] > 0 && possible == "rfix")
         falsch = false;
       if((modell_alt[names_fixed.size()-2+i] == -1 || modell_alt[names_fixed.size()-2+i] == 0) && possible == "spline")
         falsch = false;
       if(modell_alt[names_fixed.size()-2+i] == 0 && (possible == "spfix" || possible == "vfix"))
         falsch = false;

       if(falsch == false)
         {
         genoptions_mult[0]->out("  NOTE: The startmodel is no hierarchical model! Choose another one.");
         return true;
         }
       }
    }

  if(likep_mult[0]->get_family() == "Gamma")
    {
    //likep_mult[0]->set_scale(0.1);
    fullcondp[0]->set_effect_zero();
    }
  schaetzen(0,kriterium_alt,true,"backfitting");

  if(tex==true)
     {
     kriterium_tex = kriterium_alt;
     make_predictor();
     }

  kriterium_neu = kriterium_alt;
  outcriterium << steps_aktuell << "   " << ST::doubletostring(kriterium_neu,8) << endl;
  outmodels << steps_aktuell << "   " << ST::doubletostring(kriterium_neu,8) << "   ";
  ST::string header;
  fertig = false;    // �berpr�ft, ob es noch nicht gerechnete Modelle gibt
  ST::string text_neu;

  bool abbruch = false;
  if(algorithm != "coorddescent")
    abbruch = stepfunctions();
  else
    {
    abbruch = koordabstieg();

    /*if(likep_mult[0]->get_family() != "Gaussian" && minim=="adaptiv")
      {
      while(modell_neu != modellematrix[0][0])
        {
        fertig = false;
        modell_neu = modell_alt;
        startiteration.erase(startiteration.begin(),startiteration.end());
        startiteration.push_back(modell_alt);
        modellematrix.erase(modellematrix.begin(),modellematrix.end());
        modellematrix.push_back(startiteration);
//        posteriormode(posttitle,true);
//        kriterium_alt = compute_criterion(); 
schaetzen(0,kriterium_alt,0,true,"backfitting");
        kriterium_neu = kriterium_alt;
        abbruch = koordabstieg();
        if(abbruch==true)
         return true;
        }
      }*/
    }

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


void STEPWISErun::schaetzen(int z, double & kriterium, bool neu, ST::string variante)
  {
  if(variante == "backfitting")
    {
    if(criterion != "CV5" && criterion != "CV10")
      {
      posteriormode(posttitle,true);
      kriterium = compute_criterion();
      }
    else
      {
      likep_mult[0]->save_weightiwls();
      kriterium = 0;
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        posteriormode(posttitle,true);
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        }
      }
    }

  if(variante == "minibackfitting")
    {
    if(criterion != "CV5" && criterion != "CV10")
      {
      minibackfitting(fullcondp);
      kriterium = compute_criterion();
      }
    else
      {
      likep_mult[0]->save_weightiwls();
      kriterium = 0;
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        minibackfitting(fullcondp);
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        }
      }
    }

  else if(variante == "nonp")
    {
    if(criterion != "CV5" && criterion != "CV10")
      {
      fullcond_alle[z]->posteriormode();
      kriterium = compute_criterion();
      }
    else
      {
      kriterium = 0;
      likep_mult[0]->save_weightiwls();
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        fullcond_alle[z]->posteriormode();
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        }
      }
    }

  else if(variante == "nonpnonp")
    {
    if(criterion != "CV5" && criterion != "CV10")
      {
      if(neu)
        fullcond_alle[z]->const_varcoeff();
      fullcond_alle[z]->posteriormode();
      fullcond_alle[z]->wiederholen(fullcond_alle[z],false);
      kriterium = compute_criterion();
      }
    else
      {
      kriterium = 0;
      likep_mult[0]->save_weightiwls();
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        if(neu)
          fullcond_alle[z]->const_varcoeff();
        fullcond_alle[z]->posteriormode();
        fullcond_alle[z]->wiederholen(fullcond_alle[z],false);
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        }
      }
    }

  else if(variante == "fixnonp")
    {
    if(criterion != "CV5" && criterion != "CV10")
      {
      fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],-1,false);
      fullcond_alle[z]->posteriormode();
      fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],1,false);
      kriterium = compute_criterion();
      }
    else
      {
      kriterium = 0;
      likep_mult[0]->save_weightiwls();
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],-1,false);
        fullcond_alle[z]->posteriormode();
        fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],1,false);
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        }
      }
    }

  else if(variante == "factor")
    {
    if(criterion != "CV5" && criterion != "CV10")
      {
      if(!neu) 
        fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                               fullcond_alle[z]->get_data_forfixedeffects(),true);
      else
        fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                               fullcond_alle[z]->get_data_forfixedeffects(),false);
      kriterium = compute_criterion();
      }
    else
      {
      kriterium = 0;
      likep_mult[0]->save_weightiwls();
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      if(!neu)
        fullcond_alle[0]->include_effect(names_nonp[z-1],
                               fullcond_alle[z]->get_data_forfixedeffects());
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                               fullcond_alle[z]->get_data_forfixedeffects(),false);
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        }
      }
    }

  else if(variante == "fix")
    {
    int col = column_for_fix(names_fixed[z]);
    vector<ST::string> name_help;
    name_help.push_back(names_fixed[z]);

    if(criterion != "CV5" && criterion != "CV10")
      {
      if(neu)
        fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(col)),false);
      else
        fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(col)),true);
      kriterium = compute_criterion();
      }
    else
      {
      kriterium = 0;
      likep_mult[0]->save_weightiwls();
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      if(!neu)
        fullcond_alle[0]->include_effect(name_help,datamatrix(D.getCol(col)));
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(col)),false);
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        }
      }
    }

  else if(variante == "nonpfix")
    {
    if(criterion != "CV5" && criterion != "CV10")
      {
      fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                               fullcond_alle[z]->get_data_forfixedeffects(),true);
      fullcond_alle[z]->wiederholen(fullcond_alle[z],false);
      kriterium = compute_criterion();
      }
    else
      {
      kriterium = 0;
      likep_mult[0]->save_weightiwls();
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      fullcond_alle[0]->include_effect(names_nonp[z-1],
                               fullcond_alle[z]->get_data_forfixedeffects());
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                               fullcond_alle[z]->get_data_forfixedeffects(),false);
        fullcond_alle[z]->wiederholen(fullcond_alle[z],false);
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        }
      }
    }

  else if(variante == "fixfix")
    {
    if(criterion != "CV5" && criterion != "CV10")
      {
      fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                               fullcond_alle[z]->get_data_forfixedeffects(),false);
      fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],1,false);
      kriterium = compute_criterion();
      }
    else
      {
      kriterium = 0;
      likep_mult[0]->save_weightiwls();
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                               fullcond_alle[z]->get_data_forfixedeffects(),false);
        fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],1,false);
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        }
      }
    }

  else if(variante == "leer")
    {
    if(criterion != "CV5" && criterion != "CV10")
      {
      fullcond_alle[0]->posteriormode_const();
      kriterium = compute_criterion();
      }
    else
      {
      kriterium = 0;
      likep_mult[0]->save_weightiwls();
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        fullcond_alle[0]->posteriormode_const();
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        }
      }
    }

  else if(variante == "nonpleer")
    {
    ST::string possible = "alles";
    if(hierarchical == true)
      fullcond_alle[z]->hierarchical(possible);
  
    if(criterion != "CV5" && criterion != "CV10")
      {
      if(possible != "valles")
        fullcond_alle[0]->posteriormode_const();
      else if(possible == "valles")  // bei Rauslassen von VC mu� zugeh�riger fixer Effekt upgedatet werden!
        fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                               fullcond_alle[z]->get_data_forfixedeffects(),false);
      fullcond_alle[z]->wiederholen(fullcond_alle[z],false);
      kriterium = compute_criterion();
      }
    else
      {
      kriterium = 0;
      likep_mult[0]->save_weightiwls();
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        if(possible != "valles")
          fullcond_alle[0]->posteriormode_const();
        else if(possible == "valles")  // bei Rauslassen von VC mu� zugeh�riger fixer Effekt upgedatet werden!
          fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                fullcond_alle[z]->get_data_forfixedeffects(),false);
        fullcond_alle[z]->wiederholen(fullcond_alle[z],false);
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        }
      }
    }
 
  } 


// -----------------------------------------------------------------------------
// ------------------- Funktionen, f�r Stepwise / Stepmin ----------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::stepfunctions(void)
  {
  ST::string tr_akt = "trace_on";
  ST::string text_neu;
  bool eins = true;
       // Schleife f�r Minimierung
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
           korrektur();  //fullcondp[0]->posteriormode_const();
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
     ST::string possible = "alles";
     // VCM
     if(hierarchical)
       fullcond_alle[z]->hierarchical(possible);

     if(modell_alt[z+names_fixed.size()-2]==-1 && fullcond_alle[z]->get_forced()==false
                                               && possible == "alles")
       modell_neu[z+names_fixed.size()-2]= 0;
     else if(modell_alt[z+names_fixed.size()-2]==0)
       {
       modell_neu[z+names_fixed.size()-2] = -1;
       // VCM
       if(possible == "vfix")
         {
         for(i=0;i<names_nonp[z-1].size();i++)
           reset_fix(names_nonp[z-1][i]);
         kriterium_alt = MAXDOUBLE;
         }
       }
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

    ST::string possible = "alles";
    if(hierarchical == true)
      fullcond_alle[i]->hierarchical(possible);

    unsigned sch;
    for(sch=1;sch<=unsigned(increment);sch++)
       {
       modell_neu = modell_alt;
       bool lambda_exist;    // zum �berpr�fen, ob neues Lambda im Vektor enthalten ist
       unsigned index = search_lambdaindex(modell_alt[names_fixed.size()-2+i],
                                lambdavec[i-1],lambda_exist);
       lambda_exist = false;
       if(index < lambdavec[i-1].size()-sch)
          lambda_exist = true;
       if(lambda_exist==true && hierarchical == true)
         {
         // f�r 2d-Spline
         if(lambdavec[i-1][index+sch] == 0 && (possible == "spline" || possible == "spfix"))
           lambda_exist = false;
         if(lambdavec[i-1][index+sch] == -1 && (possible == "spline" || possible == "raus"))
           lambda_exist = false;
         if(lambdavec[i-1][index+sch] > 0 && (possible == "rfix" || possible == "raus"))
           lambda_exist = false;

         // f�r VCM
         if(lambdavec[i-1][index+sch] == 0 && possible == "vfix")
           lambda_exist = false;
         if(lambdavec[i-1][index+sch] == -1 && possible == "vfix")
           {   // bedeutet, da� fixer Effekt vorher nicht im Modell war, aber durch VC aufgenommen wurde.
           for(i=0;i<names_nonp[z-1].size();i++)
            reset_fix(names_nonp[z-1][i]);
           }
         }
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
       if(lambda_exist==true && hierarchical == true)
         {
         // f�r 2d-Spline
         if(lambdavec[i-1][index-sch] == 0 && (possible == "spline" || possible == "spfix"))
           lambda_exist = false;
         if(lambdavec[i-1][index-sch] == -1 && (possible == "spline" || possible == "raus"))
           lambda_exist = false;
         if(lambdavec[i-1][index-sch] > 0 && (possible == "rfix" || possible == "raus"))
           lambda_exist = false;

         // f�r VCM
         if(lambdavec[i-1][index-sch] == 0 && possible == "vfix")
           lambda_exist = false;
         if(lambdavec[i-1][index-sch] == -1 && possible == "vfix")
           {   // bedeutet, da� fixer Effekt vorher nicht im Modell war, aber durch VC aufgenommen wurde.
           for(i=0;i<names_nonp[z-1].size();i++)
             reset_fix(names_nonp[z-1][i]);
           }
         }
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

    unsigned y;
    for(y=1;y<fullcondp.size();y++)
      {
      if(fullcondp[y] == fullcond_alle[i])
        fullcondp[y]->remove_centering();
      }
    for(y=1;y<fullcond_alle.size();y++)
      fullcond_alle[y]->set_center(false);    // sorgt daf�r, da� Funktionen nicht zentriert werden!

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

    for(y=1;y<fullcond_alle.size();y++)
     {
     if(fullcond_alle[y]->is_identifiable() == false)
       fullcond_alle[y]->set_center(true);    // sorgt daf�r, da� Funktionen zentriert werden!
     }

    modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][lambda_ind];
    if(modell_neu[names_fixed.size()-2+i] != modell_alt[names_fixed.size()-2+i])
      {
      if(modelcomparison(modell_neu,modellematrix)==false)
        {
        newmodel_nonp(i,kriteriumiteration2,modeliteration,textiteration);
        // Stellt linearen Pr�diktor wieder her. Besser w�re, den lin. Pr�diktor zu speichern!!!
        korrektur(); //fullcondp[0]->posteriormode_const();
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
      if(modell_alt[i+names_fixed.size()-2]==0)
        minexact_nonp_leer(i,krit_fkt,kriterium_alt);
      else if(modell_alt[i+names_fixed.size()-2]==-1)
        {
        reset_fix(names_nonp[i-1][0]);
        minexact_nonp_fix(i,krit_fkt,kriterium_alt);
        }
      else
        minexact_nonp_nonp(i,krit_fkt,kriterium_alt);

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
        // Stellt linearen Pr�diktor wieder her.
        korrektur();  // fullcondp[0]->posteriormode_const();
        posteriormode(posttitle,true);
        }
      }
    }
  }


// -----------------------------------------------------------------------------
// ------------------ Funktionen f�r Stepmin -----------------------------------
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
  fullcond_alle[0]->safe_const();
  reset_fix(names_fixed[i]);
  schaetzen(0,kriterium_neu,true,"leer");

  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    double kriterium_test;
    schaetzen(-1,kriterium_test,false,"backfitting");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_fixed[i] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    include_fix(names_fixed[i]);
    korrektur();  // fullcondp[0]->posteriormode_const();
    posteriormode(posttitle,true);
    reset_fix(names_fixed[i]);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_fixed[i] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_alt,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(kriterium_neu < kriterium_alt)
    {
    modell_neu[i-1] = 0;
    if(modelcomparison(modell_neu,modellematrix)==false)
      {
      newmodel(kriteriumiteration2,modeliteration,textiteration);
      include_fix(names_fixed[i]);
      korrektur();  // fullcondp[0]->posteriormode_const();
      posteriormode(posttitle,true);
      }
    else
      {
      int c = column_for_fix(names_fixed[i]);
      vector<ST::string> name_help;
      name_help.push_back(names_fixed[i]);
      fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
      }
    modell_neu[i-1] = -1;
    }
  else
    {
    int c = column_for_fix(names_fixed[i]);
    vector<ST::string> name_help;
    name_help.push_back(names_fixed[i]);
    fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
    }
  }

void STEPWISErun::stepmin_leer_fix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration, unsigned & i)
  {
  fullcond_alle[0]->safe_const();
  schaetzen(i,kriterium_neu,false,"fix");
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    double kriterium_test;
    schaetzen(-1,kriterium_test,false,"backfitting");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_fixed[i] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    reset_fix(names_fixed[i]);
    int c = column_for_fix(names_fixed[i]);
    vector<ST::string> name_help;
    name_help.push_back(names_fixed[i]);
    korrektur();  // fullcondp[0]->posteriormode_const();
    posteriormode(posttitle,true);
    fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_fixed[i] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_alt,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(kriterium_neu < kriterium_alt)
     {
     modell_neu[i-1] = -1;
     if(modelcomparison(modell_neu,modellematrix)==false)
       {
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       reset_fix(names_fixed[i]);
       korrektur();  // fullcondp[0]->posteriormode_const();
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
     ST::string possible = "alles";
     // VCM
     if(hierarchical)
       fullcond_alle[z]->hierarchical(possible);

     if(modell_alt[z+names_fixed.size()-2]==-1 && fullcond_alle[z]->get_forced()==false
                                               && possible == "alles")
       stepmin_factor_leer(kriteriumiteration2,modeliteration,textiteration,z);
     else if(modell_alt[z+names_fixed.size()-2]==0)
       {
       // VCM
       if(possible == "vfix")
         {
         for(unsigned i=0;i<names_nonp[z-1].size();i++)
           reset_fix(names_nonp[z-1][i]);
         kriterium_alt = MAXDOUBLE;
         }
       stepmin_leer_factor(kriteriumiteration2,modeliteration,textiteration,z);
       }
     z = z + 1;
     }
  return z;
  }

void STEPWISErun::stepmin_factor_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  unsigned i;
  fullcond_alle[0]->safe_const();
  for(i=0;i<names_nonp[z-1].size();i++)
    reset_fix(names_nonp[z-1][i]);
  schaetzen(0,kriterium_neu,true,"leer");
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    double kriterium_test;
    schaetzen(-1,kriterium_test,false,"backfitting");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    fullcondp[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
    korrektur();  // fullcondp[0]->posteriormode_const();
    posteriormode(posttitle,true);
    for(i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_alt,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(kriterium_neu < kriterium_alt)
     {
     modell_neu[z+names_fixed.size()-2] = 0;
     if(modelcomparison(modell_neu,modellematrix)==false)
       {
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       fullcondp[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
       korrektur();  // fullcondp[0]->posteriormode_const();
       posteriormode(posttitle,true);
       }
     else
       fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                  fullcond_alle[z]->get_data_forfixedeffects(),true);
     modell_neu[z+names_fixed.size()-2] = -1;
     }
  else
     fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                 fullcond_alle[z]->get_data_forfixedeffects(),true);
  }

void STEPWISErun::stepmin_leer_factor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  unsigned i;
  fullcond_alle[0]->safe_const();
  schaetzen(z,kriterium_neu,false,"factor");
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    double kriterium_test;
    schaetzen(-1,kriterium_test,false,"backfitting");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    for(i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    korrektur();  // fullcondp[0]->posteriormode_const();
    posteriormode(posttitle,true);
    fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                 fullcond_alle[z]->get_data_forfixedeffects(),true);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_alt,6) + "\n");
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
       korrektur();  // fullcondp[0]->posteriormode_const();
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
  if(smoothing == "local")
    fullcond_alle[z]->set_lambda_nr();

  unsigned i;
/*double versuch = criterion_min(df);
versuch = 3;                                       // bis hier
for(i=0;i<fullcondp.size();i++)
  fullcondp[i]->hilfeee(); */

  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);

  vector<FULLCOND*> fullcond_ori = fullcondp;
  unsigned pos;
  bool miniback = false;
  if(!miniback_off)
    miniback = blockbilden(fullcondp,z,pos);
  if(miniback==false)
    fullcondp = fullcond_ori;

  if(minim == "adaptiv" || minim == "adap_exact" || criterion == "CV5" || criterion == "CV10")
    {
    fullcond_alle[z]->update_stepwise(modell_alt[z+names_fixed.size()-2]);
    if(miniback == false)
      schaetzen(z,kriterium,true,"nonpnonp");
    else
      schaetzen(z,kriterium,true,"minibackfitting");
    }

  fullcondp[0]->safe_const();
  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=modell_alt[z+names_fixed.size()-2])
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[z-1][i]!=-1 && lambdavec[z-1][i]!=0)
        {
        if(possible == "alles" || possible == "valles" || possible == "spline" || possible == "spfix")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
          if(miniback == false)
            schaetzen(z,kriterium_versuch,false,"nonpnonp");
          else
            schaetzen(z,kriterium_versuch,false,"minibackfitting");
          fullcond_alle[0]->set_const_old();
          }
        }
      else if(lambdavec[z-1][i]==-1)
        {
        if(possible == "alles" || possible == "spfix")
          {
          fullcond_alle[z]->set_inthemodel(-1);
          fullcond_alle[z]->reset_effect(0);

          if(df_exact == false)
            {
            if(miniback == false)
              schaetzen(z,kriterium_versuch,false,"nonpfix");
            else
              {
              minifullcond_aendern(fullcond_alle[z],fullcondp,pos);
              fullcond_alle[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
              schaetzen(z,kriterium_versuch,false,"minibackfitting");
              }
            }
          else // if(df_exact == true)
            {
            vector<FULLCOND*> fullcond_start = fullcondp;
            vector<double> modell1 = modell_alt;
            modell1[z+names_fixed.size()-2] = -1;
            fullcond_einzeln(modell1,modell_alt,z);
            fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                               fullcond_alle[z]->get_data_forfixedeffects(),false);
            fullcond_alle[z]->wiederholen(fullcond_alle[z],false);
            kriterium_versuch = compute_criterion();
            fullcondp = fullcond_start;
            end[0] = fullcondp.size()-1;
            }

          fullcond_alle[0]->set_const_old();
          reset_fix(names_nonp[z-1][0]);
          }
        }
      else
        {
        if(possible == "alles" || possible == "valles")
          {
          fullcond_alle[z]->set_inthemodel(0);
          fullcond_alle[z]->reset_effect(0);

          vector<FULLCOND*> fullcond_start;
          if(df_exact == true)
            {
            fullcond_start = fullcondp;
            vector<double> modell1 = modell_alt;
            modell1[z+names_fixed.size()-2] = 0;
            fullcond_einzeln(modell1,modell_alt,z);
            }

          if(miniback == false)
            schaetzen(z,kriterium_versuch,false,"nonpleer");
          else
            {
            minifullcond_aendern(fullcond_alle[z],fullcondp,pos);
            schaetzen(z,kriterium_versuch,false,"minibackfitting");
            }

          if(df_exact == true)
            {
            fullcondp = fullcond_start;
            end[0] = fullcondp.size()-1;
            }

          fullcond_alle[0]->set_const_old();
          }
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }

  fullcond_alle[z]->set_inthemodel(1);
  fullcond_alle[z]->update_stepwise(modell_alt[z+names_fixed.size()-2]);
  if(miniback == true)
    {
    minifullcond_aendern(fullcond_alle[z],fullcondp,pos);
    schaetzen(z,kriterium,false,"minibackfitting");
    fullcondp = fullcond_ori;
    fullcond_alle[0]->merge_data();
    }
  else
    {
    if(fullcond_alle[z]->is_identifiable() == false)
      fullcond_alle[z]->set_center(true);
    fullcond_alle[z]->posteriormode();
    fullcond_alle[z]->wiederholen(fullcond_alle[z],true);          // f�r Haupteffekt
    fullcond_alle[z]->remove_centering_fix();                 // f�r Interaktion mit fixem Haupteffekt
    fullcond_alle[0]->update_linold();
    }

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[z-1][i],6).helpfill(8)
                               + "   " + ST::doubletostring(krit_fkt[i],12) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(minim == "approx_control")
    {
    vector<double> kriterium_control;
    genoptions_mult[0]->out("\n\n");
    unsigned y;
    for(y=1;y<fullcond_alle.size();y++)
      {
      if(fullcond_alle[y]->is_identifiable() == false)
        fullcond_alle[y]->set_center(true);
      }
    minexact_nonp_nonp(z,kriterium_control,kriterium);
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalues: approx    exact \n");
    for(i=0;i<kriterium_control.size();i++)
      genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],12) + "   "
                        + ST::doubletostring(kriterium_control[i],6) + "\n");
    }
    //  else
    //    {

/*versuch = criterion_min(df + fullcond_alle[z]->compute_df());       // nur Kontrolle
versuch = 3;
for(i=0;i<fullcondp.size();i++)
  fullcondp[i]->hilfeee(); */
    //    }
  }


void STEPWISErun::stepmin_nonp_fix(unsigned & z, vector<double> & krit_fkt, double & kriterium)
  {
  unsigned i;
/*double versuch = criterion_min(df);
versuch = 3;                                       // bis hier
for(i=0;i<fullcondp.size();i++)
  fullcondp[i]->hilfeee(); */

  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);

  vector<FULLCOND*> fullcond_ori = fullcondp;
  unsigned pos;
  bool miniback = false;
  if(!miniback_off)
    miniback = blockbilden(fullcondp,z,pos);
  if(miniback==false)
    fullcondp = fullcond_ori;

  if(miniback==false)
    fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],-1,true);
  if(minim == "adaptiv" || minim == "adap_exact" || criterion == "CV5" || criterion == "CV10")
    {
    if(miniback == false)
      schaetzen(z,kriterium,true,"fixfix");
    else
      schaetzen(z,kriterium,true,"minibackfitting");
    }
  else if(miniback==false)
    fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],1,false);

  fullcond_alle[0]->safe_const();
  reset_fix(names_nonp[z-1][0]);
  if(miniback==false)
    fullcondp.push_back(fullcond_alle[z]);
  fullcond_alle[z]->set_inthemodel(1);

  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=-1)
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[z-1][i]!=0)
        {
        if(possible == "alles" || possible == "spfix" || possible == "vfix")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
          if(miniback == false)
            schaetzen(z,kriterium_versuch,false,"fixnonp");
          else
            {
            minifullcond_aendern(fullcond_alle[z],fullcondp,pos);
            schaetzen(z,kriterium_versuch,false,"minibackfitting");
            }
          fullcond_alle[0]->set_const_old();
          }
        }
      else
        {
        if(possible == "alles" || possible == "rfix")
          {
          fullcond_alle[z]->set_inthemodel(0);
          fullcond_alle[z]->reset_effect(0);

          if(df_exact)
            fullcondp.erase(fullcondp.end()-1,fullcondp.end());    // wegen df_exact

          if(miniback == false)
            schaetzen(0,kriterium_versuch,true,"leer");
          else
            {
            minifullcond_aendern(fullcond_alle[z],fullcondp,pos);
            schaetzen(z,kriterium_versuch,false,"minibackfitting");
            }

          if(df_exact)
            fullcondp.push_back(fullcond_alle[z]);                 // wegen df_exact
            
          fullcond_alle[0]->set_const_old();
          }
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }

  fullcond_alle[z]->set_inthemodel(-1);
  fullcond_alle[z]->reset_effect(0);
  if(miniback == true)
    {
    minifullcond_aendern(fullcond_alle[z],fullcondp,pos);
    fullcond_alle[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
    schaetzen(z,kriterium,false,"minibackfitting");
    fullcondp = fullcond_ori;
    fullcond_alle[0]->merge_data();
    }
  else
    {
    fullcondp.erase(fullcondp.end()-1,fullcondp.end());
    fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],-1,false);
    fullcondp[0]->posteriormode_single(names_nonp[z-1],
                                fullcond_alle[z]->get_data_forfixedeffects(),true);
    fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],1,true);
    fullcond_alle[0]->update_linold();
    }

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[z-1][i],6).helpfill(8)
                              + "   " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(minim == "approx_control")
    {
    reset_fix(names_nonp[z-1][0]);
    vector<double> kriterium_control;
    //if(df_exact)
    //  fullcondp.erase(fullcondp.end()-1,fullcondp.end());

    unsigned y;
    for(y=1;y<fullcond_alle.size();y++)
      {
      if(fullcond_alle[y]->is_identifiable() == false)
        fullcond_alle[y]->set_center(true);
      }
    minexact_nonp_fix(z,kriterium_control,kriterium);
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalues: approx    exact \n");
    for(i=0;i<kriterium_control.size();i++)
      genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],6) + "   "
             + ST::doubletostring(kriterium_control[i],6) + "\n");
    }
//  else
//    {
    //if(df_exact)
    //  fullcondp.erase(fullcondp.end()-1,fullcondp.end());

/*for(i=0;i<fullcondp.size();i++)
  fullcondp[i]->hilfeee();
double versuch = criterion_min(df + 1);       // nur Kontrolle
versuch = 3;*/
//    }
  }


void STEPWISErun::stepmin_nonp_leer(unsigned & z, vector<double> & krit_fkt, double & kriterium)
  {
  unsigned i;
      // nur Kontrolle
/*fullcond_alle[0]->posteriormode_const();
//fullcond_alle[z]->wiederholen(fullcond_alle[z]);
//fullcond_alle[0]->posteriormode_const();
double versuch = criterion_min(df);
versuch = 3;*/                                       // bis hier

  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);
  if(possible == "valles")
    possible = "alles";
    
  // VCM
  if(possible == "vfix")
    {
    for(i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    kriterium = MAXDOUBLE;
    fullcond_alle[z]->set_inthemodel(0);
    modell_alt[names_fixed.size()-2+z] = -1;
    }

   vector<FULLCOND*> fullcond_ori = fullcondp;
   unsigned pos;
   bool miniback = false;
   if(!miniback_off)
     miniback = blockbilden(fullcondp,z,pos);
   if(miniback==false)
     fullcondp = fullcond_ori;

  // if(minim == "adaptiv" || minim == "adap_exact")
           // ---> hier nicht n�tig (siehe koordmin_leer_fix)
  if((criterion == "CV5" || criterion == "CV10") && possible != "vfix")
    {
    if(miniback == false)
      schaetzen(i,kriterium,true,"leer");
    else
      schaetzen(z,kriterium,true,"minibackfitting");
    }

  if(miniback==false)
    {
    fullcondp = fullcond_ori;
    fullcond_alle[z]->const_varcoeff();   // Konstante anpassen f�r varrierenden Koeffizienten
    }

/*fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],-1,true);
fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],1,false);*/

  fullcond_alle[0]->safe_const();
  if(miniback == false)
    fullcondp.push_back(fullcond_alle[z]);
  fullcond_alle[z]->set_inthemodel(1);

  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=0)
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[z-1][i]!=-1)
        {
        if(possible == "alles")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
          if(miniback == false)
            schaetzen(z,kriterium_versuch,false,"nonp");
          else
            {
            minifullcond_aendern(fullcond_alle[z],fullcondp,pos);
            schaetzen(z,kriterium_versuch,false,"minibackfitting");
            }
          fullcond_alle[0]->set_const_old();
          }
        }
      else
        {
        if(possible == "rfix" || possible == "alles" || possible == "vfix")
          {
          fullcond_alle[z]->set_inthemodel(-1);
          fullcond_alle[z]->reset_effect(0);

          if(df_exact)
            fullcondp.erase(fullcondp.end()-1,fullcondp.end());                               // wegen df_exact

          if(miniback == false)
            schaetzen(z,kriterium_versuch,false,"factor");
          else
            {
            minifullcond_aendern(fullcond_alle[z],fullcondp,pos);
            fullcond_alle[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
            schaetzen(z,kriterium_versuch,false,"minibackfitting");
            }

          if(df_exact)
            fullcondp.push_back(fullcond_alle[z]);                                             // wegen df_exact

          reset_fix(names_nonp[z-1][0]);
          fullcond_alle[0]->set_const_old();
          }
        }
      krit_fkt.push_back(kriterium_versuch); // L�nge des Vektors mu� zu Anzahl M�glichkeiten passen!
      }
    else
      krit_fkt.push_back(kriterium);
    }

  fullcond_alle[z]->set_inthemodel(0);
  fullcond_alle[z]->reset_effect(0);
  if(miniback == true)
    {
    minifullcond_aendern(fullcond_alle[z],fullcondp,pos);
    schaetzen(z,kriterium,false,"minibackfitting");
    fullcondp = fullcond_ori;
    fullcond_alle[0]->merge_data();
    }
  else
    {
    fullcondp.erase(fullcondp.end()-1,fullcondp.end());
    fullcond_alle[0]->posteriormode_const();
    }

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[z-1][i],6).helpfill(8)
                              + "   " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(minim == "approx_control")
    {
/*fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],-1,false);
fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],1,true);
fullcond_alle[0]->update_linold();*/

    vector<double> kriterium_control;
    unsigned y;
    for(y=1;y<fullcond_alle.size();y++)
      {
      if(fullcond_alle[y]->is_identifiable() == false)
        fullcond_alle[y]->set_center(true);
      }
    minexact_nonp_leer(z,kriterium_control,kriterium);
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalues: approx    exact \n");
    for(i=0;i<kriterium_control.size();i++)
      genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],6) + "   "
             + ST::doubletostring(kriterium_control[i],6) + "\n");
    }
//  else
//    {

/*unsigned y;
for(y=0;y<fullcondp.size();y++)
  fullcondp[y]->hilfeee();
double versuch = criterion_min(df);       // nur Kontrolle
versuch = 3;*/
//    }
  }

//------------------------------------------------------------------------------------

void STEPWISErun::minexact_nonp_nonp(unsigned & z, vector<double> & krit_fkt,
                  double & kriterium)
  {
  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);
  if(possible == "valles")
    possible = "alles";

  unsigned i;
  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=modell_alt[z+names_fixed.size()-2])
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[z-1][i]!=-1 && lambdavec[z-1][i]!=0)
        {
        if(possible == "alles" || possible == "spline" || possible == "spfix")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
          korrektur();  // fullcondp[0]->posteriormode_const();
          schaetzen(z,kriterium_versuch,false,"backfitting");
          }
        }
      else if(lambdavec[z-1][i]==-1)
        {
        if(possible == "alles" || possible == "spfix")
          {
          fullcond_alle[z]->set_inthemodel(-1);
          vector<FULLCOND*> fullcond_start = fullcondp;
          vector<double> modell1 = modell_alt;
          modell1[z+names_fixed.size()-2] = -1;
          fullcond_einzeln(modell1,modell_alt,z);  // hier mu� der Fullcond-Vekor angepa�t werden!!!
          korrektur();  // fullcondp[0]->posteriormode_const();
          schaetzen(z,kriterium_versuch,false,"backfitting");
          fullcondp = fullcond_start;
          end[0] = fullcondp.size()-1;
          reset_fix(names_nonp[z-1][0]);
          }
        }
      else
        {
        if(possible == "alles")
          {
          fullcond_alle[z]->set_inthemodel(0);
          vector<FULLCOND*> fullcond_start = fullcondp;
          vector<double> modell1 = modell_alt;
          modell1[z+names_fixed.size()-2] = 0;
          fullcond_einzeln(modell1,modell_alt,z);  // hier mu� der Fullcond-Vekor angepa�t werden!!!
          korrektur();  // fullcondp[0]->posteriormode_const();
          schaetzen(z,kriterium_versuch,false,"backfitting");
          fullcondp = fullcond_start;
          end[0] = fullcondp.size()-1;
          }
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }
  fullcond_alle[z]->set_inthemodel(1);
  fullcond_alle[z]->update_stepwise(modell_alt[z+names_fixed.size()-2]);
  korrektur();  // fullcondp[0]->posteriormode_const();
  posteriormode(posttitle,true);

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (exact): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[z-1][i],6).helpfill(8) + "   " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }
  }


void STEPWISErun::minexact_nonp_fix(unsigned & z, vector<double> & krit_fkt,
          double & kriterium)
  {
  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);

  unsigned i;
  vector<FULLCOND*> fullcond_begin = fullcondp;
  vector<double> modell1 = modell_alt;
  modell1[z+names_fixed.size()-2] = 1;
  fullcond_einzeln(modell1,modell_alt,z);  // hier mu� der Fullcond-Vekor angepa�t werden!!!
  fullcond_alle[z]->set_inthemodel(1);
  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=-1)
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[z-1][i]!=0)
        {
        if(possible == "alles" || possible == "spfix" || possible == "vfix")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
          korrektur();  // fullcondp[0]->posteriormode_const();
          schaetzen(z,kriterium_versuch,false,"backfitting");
          }
        }
      else
        {
        if(possible == "alles" || possible == "rfix")
          {
          fullcond_alle[z]->set_inthemodel(0);
          fullcondp = fullcond_begin;
          end[0] = fullcondp.size()-1;
          fullcond_alle[z]->reset_effect(0);
          korrektur();  // fullcondp[0]->posteriormode_const();
          schaetzen(z,kriterium_versuch,false,"backfitting");
          }
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }
  fullcond_alle[z]->set_inthemodel(-1);
  fullcond_alle[z]->reset_effect(0);
  fullcond_alle[0]->include_effect(names_nonp[z-1],
                                fullcond_alle[z]->get_data_forfixedeffects());
  fullcondp = fullcond_begin;
  end[0] = fullcondp.size()-1;
  korrektur();  // fullcondp[0]->posteriormode_const();
  posteriormode(posttitle,true);

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (exact): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[z-1][i],6).helpfill(8) + "   " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }
  }


void STEPWISErun::minexact_nonp_leer(unsigned & z, vector<double> & krit_fkt,
                  double & kriterium)
  {
  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);
  if(possible == "valles")
    possible = "alles";

  // VCM
  if(possible == "vfix")
    {
    for(unsigned i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    kriterium = MAXDOUBLE;
    }

  unsigned i;
  vector<FULLCOND*> fullcond_begin = fullcondp;
  vector<double> modell1 = modell_alt;
  modell1[z+names_fixed.size()-2] = 1;
  fullcond_einzeln(modell1,modell_alt,z);  // hier mu� der Fullcond-Vekor angepa�t werden!!!
  fullcond_alle[z]->set_inthemodel(1);
  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=0)
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[z-1][i]!=-1)
        {
        if(possible == "alles")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
          korrektur();  // fullcondp[0]->posteriormode_const();
          schaetzen(z,kriterium_versuch,false,"backfitting");
          }
        }
      else
        {
        if(possible == "rfix" || possible == "alles" || possible == "vfix")
          {
          fullcond_alle[z]->set_inthemodel(-1);
          fullcondp = fullcond_begin;
          end[0] = fullcondp.size()-1;
          fullcond_alle[z]->reset_effect(0);
          fullcondp[0]->include_effect(names_nonp[z-1],
                                fullcond_alle[z]->get_data_forfixedeffects());
          korrektur();  // fullcondp[0]->posteriormode_const();
          schaetzen(z,kriterium_versuch,false,"backfitting");
          reset_fix(names_nonp[z-1][0]);
          }
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }
  fullcond_alle[z]->set_inthemodel(0);
  fullcond_alle[z]->reset_effect(0);
  fullcondp = fullcond_begin;
  end[0] = fullcondp.size()-1;
  korrektur();  // fullcondp[0]->posteriormode_const();
  posteriormode(posttitle,true);

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (exact): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[z-1][i],6).helpfill(8) + "   " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }
  }

//-----------------------------------------------------------------------------------

double STEPWISErun::criterion_min(const double & df)
  {
  double df1;
  if(df_exact == true)
    df1 = df_ganzehatmatrix();
  else df1 = df;

  double kriterium;
  if(criterion=="GCV")
    kriterium = likep_mult[0]->compute_gcv(df1);
  else if(criterion=="GCV2")
    kriterium = likep_mult[0]->compute_gcv2(df1);
  else if(criterion=="AIC")
    kriterium = likep_mult[0]->compute_aic(df1);
  else if(criterion=="BIC")
    kriterium = likep_mult[0]->compute_bic(df1);
  else if(criterion=="AIC_imp")
    kriterium = likep_mult[0]->compute_improvedaic(df1);
  else if(criterion=="MSEP" || criterion=="CV5" || criterion=="CV10")
    kriterium = likep_mult[0]->compute_msep();
  else //if(criterion=="AUC")
    kriterium = -1 * likep_mult[0]->compute_auc();

if(criterion=="CV5" || criterion=="CV10")
  kriterium = kriterium / likep_mult[0]->get_nrobs();

  return kriterium;
  }

double STEPWISErun::criterion_min(const double & df, const ST::string & auswahl)
  {
  double df1;
  if(df_exact == true)
    df1 = df_ganzehatmatrix();
  else
    df1 = df;

  double kriterium;
  if(auswahl=="GCV")
    kriterium = likep_mult[0]->compute_gcv(df1);
  else if(auswahl=="GCV2")
    kriterium = likep_mult[0]->compute_gcv2(df1);
  else if(auswahl=="AIC")
    kriterium = likep_mult[0]->compute_aic(df1);
  else if(auswahl=="BIC")
    kriterium = likep_mult[0]->compute_bic(df1);
  else if(auswahl=="AIC_imp")
    kriterium = likep_mult[0]->compute_improvedaic(df1);
  else if(auswahl=="MSEP")
    kriterium = likep_mult[0]->compute_msep();
  //else //if(auswahl=="AUC")
  //  kriterium = -1 * likep_mult[0]->compute_auc();

  return kriterium;
  }

// -----------------------------------------------------------------------------
// ------------------ Funktionen f�r Koordinatenmethode ------------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::koordabstieg(void)
  {
      // Schleife f�r Minimierung
  ST::string tr_akt = "trace_on";
  ST::string text_neu;
  bool eins = true;
  double kriterium_aktuell;
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
                                     && modellematrix.size()>=3 && ganze_matrix == false)
           {
           unsigned hilfe = modellematrix.size()-3;
           if(modell_alt == modellematrix[hilfe][modellematrix[hilfe].size()-1])
              fertig = true;
           }
         else if(minim == "adaptiv" && ganze_matrix == true && modellematrix.size()>=3)
           {
           unsigned hilfe = modellematrix.size()-3;
           if(modell_alt == modellematrix[hilfe][modellematrix[hilfe].size()-1] && df_exact == false)
               {
               df_exact = true;
               genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
               genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
               genoptions_mult[0]->out("\n\n");
               genoptions_mult[0]->out("\n   From now on, degrees of freedom are based on the overall hatmatrix! \n");
               }
           else if(modell_alt == modellematrix[hilfe][modellematrix[hilfe].size()-1] && df_exact == true)
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
        schaetzen(0,kriterium_alt,true,"backfitting");
        }
      kriterium_neu = kriterium_alt;
      }
    else if(fertig == true && minim == "adap_exact")
      {
      schaetzen(0,kriterium_alt,true,"backfitting");
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
  double kriterium_adaptiv = kriterium_aktuell;
  if(minim == "adaptiv" || minim == "adap_exact" || minim == "adaptiv_golden"
                        || criterion == "CV5" || criterion == "CV10")
    {
    schaetzen(i,kriterium_aktuell,true,"fix");
    }

  modell_neu[i-1] = 0;
  fullcond_alle[0]->safe_const();
  reset_fix(names_fixed[i]);
  schaetzen(0,kriterium_neu,false,"leer");
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    double kriterium_test;
    schaetzen(-1,kriterium_test,false,"backfitting");
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
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_fixed[i] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_aktuell,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

 if(minim != "adaptiv" && minim != "adap_exact" && minim!= "adaptiv_golden")
   {
   if(kriterium_neu < kriterium_aktuell)
     {
     kriterium_aktuell = kriterium_adaptiv;
     bool neu = modelcomparison(modell_neu,modellematrix);
     if(neu==false)
       {
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       kriterium_neu = kriteriumiteration2[kriteriumiteration2.size()-1];
       }
     else
       kriterium_neu = kriterium_aktuell;

     if(neu==true || kriterium_aktuell < kriterium_neu)
       {
       int c = column_for_fix(names_fixed[i]);
       vector<ST::string> name_help;
       name_help.push_back(names_fixed[i]);
       fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
       modell_neu[i-1] = -1;
       if(kriterium_aktuell < kriterium_neu) // verhindert, da� "approx" schlechter wird!
         {
         posteriormode(posttitle,true);
         if( (trace == "trace_minim" || trace == "trace_on") && neu == false)
             genoptions_mult[0]->out("\n\n  Trial won't become the new model! \n");
         }
       }
     else
       kriterium_aktuell = kriterium_neu;
     }
   else //if(kriterium_neu > kriterium_aktuell)
     {
     kriterium_neu = kriterium_adaptiv;
     kriterium_aktuell = kriterium_adaptiv;
     int c = column_for_fix(names_fixed[i]);
     vector<ST::string> name_help;
     name_help.push_back(names_fixed[i]);
     fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
     modell_neu[i-1] = -1;
     }
   }

  if(minim == "adaptiv" || minim == "adap_exact" || minim == "adaptiv_golden")
    {
    if(kriterium_aktuell >= kriterium_neu)
      kriterium_aktuell = kriterium_neu;
   else //if(kriterium_neu > kriterium_aktuell)
     {
     kriterium_neu = kriterium_adaptiv;
     kriterium_aktuell = kriterium_adaptiv;
     int c = column_for_fix(names_fixed[i]);
     vector<ST::string> name_help;
     name_help.push_back(names_fixed[i]);
     fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
     modell_neu[i-1] = -1;
     }
      
    if(fabs((kriterium_adaptiv - kriterium_aktuell)/kriterium_adaptiv) >= std::pow(10,-6.0))
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
  double kriterium_adaptiv = kriterium_aktuell;
  // if(minim == "adaptiv" || minim == "adap_exact")
     // ---> hier nicht n�tig, weil Intercept immer erneuert wird!
  if(criterion == "CV5" || criterion == "CV10")
    {
    schaetzen(i,kriterium_aktuell,true,"leer");
    }

  modell_neu[i-1] = -1;
  fullcond_alle[0]->safe_const();
  schaetzen(i,kriterium_neu,false,"fix");
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    double kriterium_test;
    schaetzen(-1,kriterium_test,false,"backfitting");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_fixed[i] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    reset_fix(names_fixed[i]);
    int c = column_for_fix(names_fixed[i]);
    vector<ST::string> name_help;
    name_help.push_back(names_fixed[i]);
    korrektur();  // fullcondp[0]->posteriormode_const();
    posteriormode(posttitle,true);
    fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_fixed[i] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_aktuell,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(minim != "adaptiv" && minim != "adap_exact" && minim != "adaptiv_golden")
    {
    if(kriterium_neu < kriterium_aktuell)
      {
      kriterium_aktuell = kriterium_adaptiv;
      bool neu = modelcomparison(modell_neu,modellematrix);
      if(neu==false)
        {
        newmodel(kriteriumiteration2,modeliteration,textiteration);
        kriterium_neu = kriteriumiteration2[kriteriumiteration2.size()-1];
        }
     else
       kriterium_neu = kriterium_aktuell;

      if(neu==true || kriterium_aktuell < kriterium_neu)
        {
        reset_fix(names_fixed[i]);
        modell_neu[i-1] = 0;
        if(kriterium_aktuell < kriterium_neu)
          {
          posteriormode(posttitle,true);
          if( (trace == "trace_minim" || trace == "trace_on") && neu == false)
             genoptions_mult[0]->out("\n\n  Trial won't become the new model! \n");
          }
        }
     else
       kriterium_aktuell = kriterium_neu;
      }
    else //if(kriterium_neu >= kriterium_aktuell)
      {
      reset_fix(names_fixed[i]);
      modell_neu[i-1] = 0;
      }
    }

  if(minim == "adaptiv" || minim == "adap_exact" || minim == "adaptiv_golden")
    {
    if(kriterium_aktuell >= kriterium_neu)
      kriterium_aktuell = kriterium_neu;
    else  // if(kriterium_neu >= kriterium_aktuell)
      {
      reset_fix(names_fixed[i]);
      modell_neu[i-1] = 0;
      }

    if(fabs((kriterium_adaptiv - kriterium_aktuell)/kriterium_adaptiv) >= std::pow(10,-6.0))
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
     ST::string possible = "alles";
     // VCM
     if(hierarchical)
       fullcond_alle[z]->hierarchical(possible);

     if(modell_alt[z+names_fixed.size()-2]==-1 && fullcond_alle[z]->get_forced()==false
                                               && possible == "alles")
       koord_factor_leer(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell,z);
     else if(modell_alt[z+names_fixed.size()-2]==-1 && (fullcond_alle[z]->get_forced()==true
                                               || possible == "vfix"))
       {
       if(minim == "adaptiv" || minim == "adap_exact" || minim == "adaptiv_golden")
                 // || criterion == "CV5" || criterion == "CV10")  --> hier nicht n�tig!
         {
         kriterium_aktuell = MAXDOUBLE;
         koord_factor_leer(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell,z);
         }
       }
     else if(modell_alt[z+names_fixed.size()-2]==0)
       {
       // VCM
       if(possible == "vfix")
         {
         for(unsigned i=0;i<names_nonp[z-1].size();i++)
           reset_fix(names_nonp[z-1][i]);
         kriterium_aktuell = MAXDOUBLE;
         fullcond_alle[z]->set_inthemodel(0);
         }
       koord_leer_factor(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell,z);
       }
     modell_alt = modell_neu;
     z = z + 1;
     }
  return z;
  }

void STEPWISErun::koord_factor_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & z)
  {
  unsigned i;
  vector<FULLCOND*> fullcond_ori = fullcondp;
  unsigned pos;
  bool miniback = false;
  if(!miniback_off)
    miniback = blockbilden(fullcondp,z,pos);
  if(miniback==false)
    fullcondp = fullcond_ori;

  double kriterium_adaptiv = kriterium_aktuell;
  if(minim == "adaptiv" || minim == "adap_exact" || minim == "adaptiv_golden"
                        || criterion == "CV5" || criterion == "CV10")
    {
    if(miniback==false)
      schaetzen(z,kriterium_aktuell,true,"factor");
    else
      {
      //fullcond_alle[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
      schaetzen(z,kriterium_aktuell,true,"minibackfitting");
      }
    }

  if(kriterium_adaptiv < MAXDOUBLE)
    {
    modell_neu[z+names_fixed.size()-2] = 0;
    fullcond_alle[z]->set_inthemodel(0);
    fullcond_alle[0]->safe_const();
    for(i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    if(miniback==false)
      {
      schaetzen(0,kriterium_neu,true,"leer");
      }
    else
      {
      schaetzen(0,kriterium_neu,false,"minibackfitting");
      fullcondp = fullcond_ori;
      fullcond_alle[0]->merge_data();
      }
    fullcond_alle[0]->set_const_old();

    if(minim == "approx_control")
      {
      double kriterium_test;
      schaetzen(-1,kriterium_test,false,"backfitting");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalue: approx = "
            + ST::doubletostring(kriterium_neu,6) + " exact = "
            + ST::doubletostring(kriterium_test,6) + "\n");
      fullcondp[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
      korrektur();  // fullcondp[0]->posteriormode_const();
      posteriormode(posttitle,true);
      for(i=0;i<names_nonp[z-1].size();i++)
        reset_fix(names_nonp[z-1][i]);
      }
    if(trace == "trace_minim" && minim != "approx_control")
      {
      genoptions_mult[0]->out("\n\n");
      genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
      genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_aktuell,6) + "\n");
      genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
      genoptions_mult[0]->out("\n");
      }

    if(minim != "adaptiv" && minim != "adap_exact" && minim != "adaptiv_golden")
      {
      if(kriterium_neu < kriterium_aktuell)
        {
        kriterium_aktuell = kriterium_adaptiv;
        bool neu = modelcomparison(modell_neu,modellematrix);
        if(neu==false)
          {
          newmodel(kriteriumiteration2,modeliteration,textiteration);
          kriterium_neu = kriteriumiteration2[kriteriumiteration2.size()-1];
          }
        else
          kriterium_neu = kriterium_aktuell;

        if(neu==true || kriterium_aktuell < kriterium_neu)
          {
          fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                         fullcond_alle[z]->get_data_forfixedeffects(),true);
          modell_neu[z+names_fixed.size()-2] = -1;
          fullcond_alle[z]->set_inthemodel(-1);
          if(kriterium_aktuell < kriterium_neu)
            {
            posteriormode(posttitle,true);
            if( (trace == "trace_minim" || trace == "trace_on") && neu == false)
              genoptions_mult[0]->out("\n\n  Trial won't become the new model! \n");
            }
          }
        else
          kriterium_aktuell = kriterium_neu;
        }
      else //if(kriterium_neu >= kriterium_aktuell)
        {
        kriterium_aktuell = kriterium_adaptiv;
        fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                     fullcond_alle[z]->get_data_forfixedeffects(),true);
        modell_neu[z+names_fixed.size()-2] = -1;
        fullcond_alle[z]->set_inthemodel(-1);
        }
      }

    if(minim == "adaptiv" || minim == "adap_exact" || minim == "adaptiv_golden")
      {
      if(kriterium_aktuell >= kriterium_neu)
        kriterium_aktuell = kriterium_neu;
      else //if(kriterium_neu >= kriterium_aktuell)
        {
        kriterium_aktuell = kriterium_adaptiv;
        fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                     fullcond_alle[z]->get_data_forfixedeffects(),true);
        modell_neu[z+names_fixed.size()-2] = -1;
        fullcond_alle[z]->set_inthemodel(-1);        
        }

      if(fabs((kriterium_adaptiv - kriterium_aktuell)/kriterium_adaptiv) >= std::pow(10,-6.0))
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
  else
    {
    if(miniback==true)
      {
      fullcondp = fullcond_ori;
      fullcond_alle[0]->merge_data();
      }
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
    genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_aktuell,6) + "\n");
    genoptions_mult[0]->out("\n");
    }
  }

void STEPWISErun::koord_leer_factor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & z)
  {
  unsigned i;
  vector<FULLCOND*> fullcond_ori = fullcondp;
  unsigned pos;
  bool miniback = false;
  if(!miniback_off)
    miniback = blockbilden(fullcondp,z,pos);
  if(miniback==false)
    fullcondp = fullcond_ori;

  double kriterium_adaptiv = kriterium_aktuell;
  // if(minim == "adaptiv" || minim == "adap_exact")
      // ---> hier �berfl�ssig (siehe oben)!
  if(criterion == "CV5" || criterion == "CV10")
    {
    if(miniback==false)
      schaetzen(z,kriterium_aktuell,true,"leer");
    else
      {
      schaetzen(z,kriterium_aktuell,true,"minibackfitting");
      }
    }

  modell_neu[z+names_fixed.size()-2] = -1;
  fullcond_alle[z]->set_inthemodel(-1);
  fullcond_alle[0]->safe_const();
  if(miniback==false)
    schaetzen(z,kriterium_neu,false,"factor");
  else
    {
    fullcond_alle[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
    schaetzen(z,kriterium_neu,false,"minibackfitting");
    fullcondp = fullcond_ori;
    fullcond_alle[0]->merge_data();
    }
  fullcond_alle[0]->set_const_old();

  if(minim == "approx_control" && kriterium_aktuell < MAXDOUBLE)
    {
    double kriterium_test;
    schaetzen(-1,kriterium_test,false,"backfitting");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    for(i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    korrektur();  // fullcondp[0]->posteriormode_const();
    posteriormode(posttitle,true);
    fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                  fullcond_alle[z]->get_data_forfixedeffects(),true);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_aktuell,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(minim != "adaptiv" && minim != "adap_exact" && minim != "adaptiv_golden")
    {
    if(kriterium_neu < kriterium_aktuell)
      {
      kriterium_aktuell = kriterium_adaptiv;
      bool neu = modelcomparison(modell_neu,modellematrix);
      if(neu==false)
        {
        newmodel(kriteriumiteration2,modeliteration,textiteration);
        kriterium_neu = kriteriumiteration2[kriteriumiteration2.size()-1];
        }
      else
        kriterium_neu = kriterium_aktuell;

      if(neu==true || kriterium_aktuell < kriterium_neu)
        {
        for(i=0;i<names_nonp[z-1].size();i++)
          reset_fix(names_nonp[z-1][i]);
        modell_neu[z+names_fixed.size()-2] = 0;
        if(kriterium_aktuell < kriterium_neu)
          {
          posteriormode(posttitle,true);
          if( (trace == "trace_minim" || trace == "trace_on") && neu == false)
             genoptions_mult[0]->out("\n\n  Trial won't become the new model! \n");
          }
        }
      else
        kriterium_aktuell = kriterium_neu;
      }
    else // if(kriterium_neu >= kriterium_aktuell)
      {
      for(i=0;i<names_nonp[z-1].size();i++)
        reset_fix(names_nonp[z-1][i]);
      modell_neu[z+names_fixed.size()-2] = 0;
      fullcond_alle[z]->set_inthemodel(0);      
      }
    }

  if(minim == "adaptiv" || minim == "adap_exact" || minim == "adaptiv_golden")
    {
    if(kriterium_aktuell >= kriterium_neu)
      kriterium_aktuell = kriterium_neu;
    else  // if(kriterium_neu >= kriterium_aktuell)
      {
      for(i=0;i<names_nonp[z-1].size();i++)
        reset_fix(names_nonp[z-1][i]);
      modell_neu[z+names_fixed.size()-2] = 0;
      fullcond_alle[z]->set_inthemodel(0);      
      }

    if(fabs((kriterium_adaptiv - kriterium_aktuell)/kriterium_adaptiv) >= std::pow(10,-6.0))
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

    unsigned y;
    for(y=1;y<fullcondp.size();y++)
      {
      if(fullcondp[y] == fullcond_alle[i])
        fullcondp[y]->remove_centering();
      }
    for(y=1;y<fullcond_alle.size();y++)
      fullcond_alle[y]->set_center(false);    // sorgt daf�r, da� Funktionen nicht zentriert werden!

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

    for(y=1;y<fullcond_alle.size();y++)
     {
     if(fullcond_alle[y]->is_identifiable() == false)
       fullcond_alle[y]->set_center(true);    // sorgt daf�r, da� Funktionen zentriert werden!
     }

    modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][lambda_ind];
    if(minim != "adaptiv" && minim != "adap_exact" && minim != "adaptiv_golden")
      {
      kriterium_aktuell = kriterium_test;
      if(modell_neu[names_fixed.size()-2+i] != modell_alt[names_fixed.size()-2+i])
        {
        bool neu = modelcomparison(modell_neu,modellematrix);
        fullcond_einzeln(modell_neu,modell_alt,i);
        if(neu==false)
          {
          korrektur();  // fullcondp[0]->posteriormode_const();
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
      /*else
        {
        posteriormode(posttitle,true);  // nur Versuch!!!
        }*/
      }
    else //if(minim == "adaptiv" || minim == "adap_exact")
      {
      if(modell_alt[names_fixed.size()-2+i] != modell_neu[names_fixed.size()-2+i])
        {
        fullcond_einzeln(modell_neu,modell_alt,i);

        vector<FULLCOND*> fullcond_ori = fullcondp;
        unsigned pos;
        bool miniback = false;
        if(!miniback_off)
          miniback = blockbilden(fullcondp,z,pos);
        if(miniback==false)
          {
          fullcondp = fullcond_ori;

          if(modell_neu[names_fixed.size()-2+i] == 0)      // noch mal �berpr�fen!!!
            {
            //fullcond_alle[i]->reset_effect(0);   // Nicht n�tig, wegen "fullcond_einzeln"-> fullcond_alle[i] nicht in fullcondp!!!
            if(modell_alt[names_fixed.size()-2+i] > 0)
              {
              fullcond_alle[i]->wiederholen_fix(fullcond_alle[i],1,true);
              if(hierarchical)                        // neu f�r VCM!  Versuch!!!
                {
                ST::string possible = "alles";
                fullcond_alle[i]->hierarchical(possible);
                if(possible == "valles")
                  fullcond_alle[i]->posteriormode_single(names_nonp[i-1],
                                   fullcond_alle[i]->get_data_forfixedeffects(),false);
                }
              }
            fullcond_alle[0]->posteriormode_const();
            }
          else if(modell_neu[names_fixed.size()-2+i] == -1)
            {
            if(modell_alt[names_fixed.size()-2+i] == 0)
              fullcond_alle[i]->wiederholen_fix(fullcond_alle[i],-1,true);

            fullcond_alle[i]->reset_effect(0);
            fullcond_alle[0]->posteriormode_single(names_nonp[i-1],
                                 fullcond_alle[i]->get_data_forfixedeffects(),false);

            //if(modell_alt[names_fixed.size()-2+i] > 0)     // Warum???
            fullcond_alle[i]->wiederholen_fix(fullcond_alle[i],1,true);
            }
          else
            {
            if(modell_alt[names_fixed.size()-2+i] > 0)
              fullcond_alle[i]->remove_centering();
            else
              fullcond_alle[i]->wiederholen_fix(fullcond_alle[i],-1,true);

            if(modell_alt[names_fixed.size()-2+i] == 0)
              fullcond_alle[i]->const_varcoeff();

            fullcond_alle[i]->update_stepwise(modell_neu[names_fixed.size()-2+i]);
            fullcond_alle[i]->posteriormode();
            fullcond_alle[0]->update_linold();
            fullcond_alle[i]->wiederholen(fullcond_alle[i],true);
            fullcond_alle[i]->remove_centering_fix();
            }
          }
        else
          {
          minibackfitting(fullcondp);
          fullcond_alle[0]->merge_data();
          fullcondp = fullcond_ori;
          }

        if(trace == "trace_on" || trace == "trace_minim")
          {
          ST::string text;
          maketext("  Trial:",modell_neu,kriterium_min,text,true,trace,false);
          }
/*for(y=0;y<fullcondp.size();y++)    // nur Kontrolle!
  fullcondp[y]->hilfeee();
double test = compute_criterion();
test = test;*/
        }
      //else
      //  maketext(header,modell_neu,kriterium_alt,text_neu,true,tr_akt,true);
      modell_alt = modell_neu;
      kriterium_aktuell = kriterium_min;
      kriterium_alt = kriterium_aktuell;
      if(fabs((kriterium_test - kriterium_aktuell)/kriterium_test) >= std::pow(10,-6.0))
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
       korrektur();  // fullcondp[0]->posteriormode_const();
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
     ST::string possible = "alles";
     // VCM
     if(hierarchical)
       fullcond_alle[z]->hierarchical(possible);

     if(modell_alt[z+names_fixed.size()-2]==-1 && fullcond_alle[z]->get_forced()==false
                                               && possible == "alles")
       modell_neu[z+names_fixed.size()-2]= 0;
     else if(modell_alt[z+names_fixed.size()-2]==0)
       {
       if(possible == "vfix")         // VCM
         {
         for(i=0;i<names_nonp[z-1].size();i++)
           reset_fix(names_nonp[z-1][i]);
         kriterium_aktuell = MAXDOUBLE;
         }
       modell_neu[z+names_fixed.size()-2] = -1;
       }
     if(modelcomparison(modell_neu,modellematrix)==false)
        {
        if(modell_neu[z+names_fixed.size()-2]==0)
          {
          for(i=0;i<names_nonp[z-1].size();i++)
            reset_fix(names_nonp[z-1][i]);
          }
        else
          fullcondp[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
        korrektur();  // fullcondp[0]->posteriormode_const();
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
      if(modell_alt[i+names_fixed.size()-2]==0)
        minexact_nonp_leer(i,krit_fkt,kriterium_aktuell);
      else if(modell_alt[i+names_fixed.size()-2]==-1)
        {
        reset_fix(names_nonp[i-1][0]);
        minexact_nonp_fix(i,krit_fkt,kriterium_aktuell);
        }
      else
        minexact_nonp_nonp(i,krit_fkt,kriterium_aktuell);

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
        korrektur();  // fullcondp[0]->posteriormode_const();
        newmodel(kriteriumiteration2,modeliteration,textiteration);
        kriterium_aktuell = kriteriumiteration2[kriteriumiteration2.size()-1];
        }
      }
    modell_alt = modell_neu;
    }
  }


// -----------------------------------------------------------------------------
// --------------------------- Mini-Backfitting --------------------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::blockbilden(vector<FULLCOND*> & fullcond_block, unsigned & z, unsigned & pos)
  {
  vector<FULLCOND*> interactionspointer;
  fullcond_alle[z]->get_interactionspointer(interactionspointer);
  if(interactionspointer.size()==0)
    return false;
  else
    {
    fullcond_block.erase(fullcond_block.begin(),fullcond_block.end());
    fullcond_block.push_back(fullcond_alle[0]);

    vector<FULLCOND*> fullcond_help;
    fullcond_help.push_back(fullcond_alle[z]);

    unsigned i,j;
    bool drin, fix;
    bool gibts = false;

    i = 0;
    while(i<fullcond_help.size())
      {
      fullcond_help[i]->get_interactionspointer(interactionspointer);
      for(j=0;j<interactionspointer.size();j++)
        {
        interactionspointer[j]->get_inthemodel(drin,fix);
        if(drin == true || fix == true)
          {
          unsigned z = 0;
          bool gefunden = false;
          while(z<fullcond_help.size() && gefunden == false)
            {
            if(fullcond_help[z] == interactionspointer[j])
              gefunden = true;
            z++;
            }
          if(gefunden == false)
            {
            fullcond_help.push_back(interactionspointer[j]);
            gibts = true;
            }
          }
        }
      i++;
      }

//for(i=0;i<fullcond_help.size();i++)
//  genoptions_mult[0]->out(fullcond_help[i]->get_datanames()[0] + "  \n");

    if(gibts == true)
      {
// in richtige Reihenfolge bringen (wie bei "fullcond_alle")
      vector<ST::string> fixenames;
      unsigned p = 0;
      for(i=1;i<fullcond_alle.size();i++)
        {
        bool gefunden = false;
        j = 0;
        while(gefunden == false && j < fullcond_help.size())
          {
          if(fullcond_alle[i] == fullcond_help[j])
            {
            gefunden = true;
            fullcond_help[j]->get_inthemodel(drin,fix);   // entscheiden, ob fix oder nonp
            if(fix==true)
              {
              fixenames.push_back(names_nonp[i-1][0]);
              if(i == z)
                pos = p+1;
              }
            else if(drin==true)
              {
              p++;
              fullcond_block.push_back(fullcond_alle[i]);
              if(i == z)
                pos = p;
              }
            else
              {
              if(i == z)
                pos = p+1;
              }
            }
          j++;
          }
        }
      fullcond_alle[0]->split_data(fixenames);
      for(i=1;i<fullcond_block.size();i++)
        {
        if(fullcond_block[i]->is_identifiable() == false)
          fullcond_block[i]->set_center(true);
        }
      if(fullcond_alle[z]->is_identifiable() == false)
        fullcond_alle[z]->set_center(true);
      }

    return gibts;
    }
  }


void STEPWISErun::minifullcond_aendern(FULLCOND* & fullcondz, vector<FULLCOND*> & fullcond, unsigned & pos)
  {
  bool drin, fix;
  fullcondz->get_inthemodel(drin,fix);
  if(drin==false)
    {
    if(pos < fullcond.size() && fullcond[pos] == fullcondz)
      fullcond.erase(fullcond.begin()+pos,fullcond.begin()+pos+1);
    }
  else if(drin==true)
    {
    if(pos >= fullcond.size())
      fullcond.push_back(fullcondz);
    else if(fullcond[pos] != fullcondz)
      //vector<FULLCOND*>::iterator fit = fullcond.begin() + pos;
      fullcond.insert(fullcond.begin()+pos,fullcondz);
    }
  }


void STEPWISErun::minibackfitting(vector<FULLCOND*> & fullcond)
  {
  unsigned it=1;
  bool converged = false;
  unsigned j;

  while ((!converged) && (it <= 100))
    {
    converged = true;

    if (likepexisting)
      if (likep_mult[0]->posteriormode() == false)
         converged = false;

    for(j=0;j<fullcond.size();j++)
      {
      if(fullcond[j]->posteriormode() == false)
        converged = false;
      } // end: for(j=0;j<fullcond.size();j++)

    #if defined(BORLAND_OUTPUT_WINDOW)
    Application->ProcessMessages();
    if(Frame->stop)
      {
      break;
      }
    if(Frame->pause)
      {
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("SIMULATION PAUSED\n");
      genoptions_mult[0]->out("Click CONTINUE to proceed\n");
      genoptions_mult[0]->out("\n");
      while (Frame->pause)
        {
        Application->ProcessMessages();
        }
      genoptions_mult[0]->out("SIMULATION CONTINUED\n");
      genoptions_mult[0]->out("\n");
      }
    #elif defined(JAVA_OUTPUT_WINDOW)
    bool stop = genoptions_mult[0]->adminb_p->breakcommand();
    if(stop)
      break;
    #endif

    it++;

    } // end:   while ((!converged) && (it <= 100))

  for(j=0;j<fullcond.size();j++)    
    fullcondp[j]->remove_centering_fix();
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
  fullcond_alle = fullcondp; // hier anders, wenn es f�r als fix ausgew. Var. wieder mehr M�gl. geben soll

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
  for(i=1;i<fullcond_alle.size();i++)      // hier wird die zweite H�lfte (Wert > Startwert) der Lambdas
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
            int intnumber = (int)number;
            fullcond_alle[i]->compute_lambdavec(untervector,intnumber);
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
  algorithm = "coorddescent";
  minim = "approx";          // bei adaptiv springen die Werte zu stark
                             // bei stepwise mu� erst noch �berpr�ft werden, ob "lambda_nr" immer richtig gesetzt wird!
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
// Vorschlag:
//       fullcond_alle[i]->set_smoothing("local");
       ST::string helpstring = "local";
       fullcond_alle[i]->set_smoothing(helpstring);

       vector<double> lambdas;
       double lambdastart = modell[names_fixed.size()-2+i];

       fullcond_alle[i]->set_lambdas_vector(lambdastart);

       double lambdamin,lambdamax;
       lambdas = lambdavec[i-1];
       if(lambdas[lambdas.size()-1] == 0)
       lambdas.erase(lambdas.end()-1,lambdas.end());
       if(lambdas[lambdas.size()-1] == -1)
         lambdas.erase(lambdas.end()-1,lambdas.end());
       lambdamin = lambdas[0];
       lambdamax = lambdas[lambdas.size()-1];

       fullcond_alle[i]->set_stepwise_options(lambdastart,lambdamax,lambdamin,true,
               fullcond_alle[i]->get_df_lambdamax(),fullcond_alle[i]->get_df_lambdamin(),false,false,nummer,false);

       vector<ST::string> names_help;
       names_help.push_back(fullcond_alle[i]->get_datanames()[fullcond_alle[i]->get_datanames().size()-1]);

       anzahl = fullcond_alle[i]->get_rankK2();
       for(j=1;j<=anzahl;j++)       
          {
          lambdavec_local.push_back(lambdas);
          modell_local.push_back(lambdastart);
          fullcond_local.push_back(fullcond_alle[i]);
          names_nonp_local.push_back(names_help);
          }
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

  outmodels << endl;
  ST::string text;
  bool abbruch = false;
  abbruch = single_stepwise(startindex[0],startfix[0],false);

  if(abbruch==true)
    return true;

  ST::string header = "  Final Model after Fine-Tuning with local smoothing parameters:";
  fix_komplett(modell_alt);
  fullcond_komplett(modell_alt);

  smoothing = "global";    // Versuch!!!
  kriterium_alt = compute_criterion();

  ST::string tr_akt = "trace_on";
  maketext(header,modell_alt,kriterium_alt,text_alt,false,tr_akt,false);
  kriterium_tex = kriterium_alt;
  genoptions_mult[0]->out("\n\n");
  modell = modell_alt;

  return false;
  }


// -----------------------------------------------------------------------------
// ------- Funktionen f�r die Erstellung des Startmodels -----------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::vcm_doppelt(void)
  {
  unsigned i = 1;
  unsigned j;
  unsigned k;
  bool fehler_randomslope = false;
  while(i<fullcond_alle.size() && fehler_randomslope == false)
     {
     j = 1;
     while(j<names_fixed.size() && fehler_randomslope == false)
       {
       if(names_fixed[j]==names_nonp[i-1][0])
         {
         k = j;
         fehler_randomslope = true;
         }
       j++;
       }
     i++;
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

       // berechnet ein Modell, um Gewichte f�r Augabe usw. zu erhalten!!!
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
         fullcondp.erase(fullcondp.begin() + 1);   //l�scht das Element an Pos.1 aus fullcondp-Vektor
         fullcondp[0]->include_effect(fullcond_alle[i]->get_datanames(),
                                   fullcond_alle[i]->get_data_forfixedeffects());
         }
      }
    end[0] = fullcondp.size()-1;

//likep_mult[0]->compute_iwls();
    posteriormode(posttitle,true);
    }
  else
   likep_mult[0]->compute_iwls();

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
     fullcond_alle[i]->set_inthemodel(1);
     fullcond_alle[i]->compute_lambdavec(untervector,nummer);
     fullcond_alle[i]->set_inthemodel(0);
     lambdavector.push_back(untervector);
     }
  }


void STEPWISErun::initialise_weights(double prop)
  {
  if(criterion == "MSEP" || criterion == "AUC")
    {
    datamatrix weight = likep_mult[0]->get_weight();
    bool gewicht = false;
    if(weight.min(0) > 0)
      {
      weight = datamatrix(weight.rows(),1,0);
      unsigned i;
      for(i=1;i<fullcond_alle.size();i++)
        fullcond_alle[i]->create_weight(weight);
      gewicht = true;
      }
    likep_mult[0]->create_weight(weight,prop,gewicht,false);
    }
  else
    {
    datamatrix w = datamatrix(1,1,0);
    bool gewicht = false;
    if(criterion== "CV5")
      likep_mult[0]->create_weight(w,5,gewicht,true);
    else
      likep_mult[0]->create_weight(w,10,gewicht,true);
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
       else
         diff.push_back(MAXDOUBLE);
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
       double start = fullcond_alle[i]->get_lambdastart();
       unsigned index = search_lambdastartindex(start,lambdavec[i-1]);
       indexhelp.push_back(index);
       }

    startindex.push_back(indexhelp);
    startfix.push_back(fixhelp);
    }
  }


// -----------------------------------------------------------------------------
// ------- Funktionen f�r die Berechnung neuer Modelle -------------------------
// -----------------------------------------------------------------------------

double STEPWISErun::compute_criterion(void)
  {
  double df = 0;
  if(criterion != "MSEP" && criterion != "AUC" && criterion != "CV5" && criterion != "CV10")
    {
    if(df_exact == false)
      {
      unsigned i;
      for(i=0;i<fullcond_alle.size();i++)
        df = df + fullcond_alle[i]->compute_df();
      likep_mult[0]->set_iwlsweights_notchanged(true);
      }
    else
      {
      df = df_ganzehatmatrix();
      }
    }
  double kriterium;

  if(criterion=="GCV")
    kriterium = likep_mult[0]->compute_gcv(df);
  else if(criterion=="GCV2")
    kriterium = likep_mult[0]->compute_gcv2(df);
  else if(criterion=="AIC")
    kriterium = likep_mult[0]->compute_aic(df);
  else if(criterion=="BIC")
    kriterium = likep_mult[0]->compute_bic(df);
  else if(criterion=="AIC_imp")
    kriterium = likep_mult[0]->compute_improvedaic(df);
  else if(criterion=="MSEP" || criterion=="CV5" || criterion=="CV10")
    kriterium = likep_mult[0]->compute_msep();
  else //if(criterion=="AUC")
    kriterium = -1 * likep_mult[0]->compute_auc();

  if(criterion=="CV5" || criterion=="CV10")
    kriterium = kriterium / likep_mult[0]->get_nrobs();
  else if(criterion=="MSEP")
    kriterium = kriterium /(likep_mult[0]->get_nrobs() - likep_mult[0]->get_nrobs_wpw());

  return kriterium;
  }


void STEPWISErun::newmodel(vector<double> & krit,
  vector<vector<double> > & mi, vector<ST::string> & textit)
  {
  fertig = false;
  mi.push_back(modell_neu);

/*BIC_min += 1;
if(BIC_min == 1500) //1698)
  {
  //ofstream out("c:\\bayesx\\output\\nr.txt");
  //out << modell_neu[0] << "  " << modell_neu[1] << "  " << modell_neu[2] << "  " << modell_neu[3] << "  " << modell_neu[4] << endl;
  ST::string scheisse = "Scheisse";
  } */

  double kriterium;
  schaetzen(0,kriterium,true,"backfitting");
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
  korrektur();  // fullcondp[0]->posteriormode_const();
  newmodel(krit,mi,textit);
  if(mo==0)
    include_fix(name);
  else
    reset_fix(name);
  //fullcondp[0]->posteriormode_const();
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
  korrektur();  // fullcondp[0]->posteriormode_const();
  newmodel(krit,mi,textit);
  if(mo==0)
    fullcondp[0]->include_effect(name,fullcond_alle[index]->get_data_forfixedeffects());
  else
    {
    for(i=0;i<name.size();i++)
      reset_fix(name[i]);
    }
  //fullcondp[0]->posteriormode_const();
  }


void STEPWISErun::newmodel_nonp(const unsigned & index,
    vector<double> & krit, vector<vector<double> > & mi,
    vector<ST::string> & textit)
  {
  fullcond_einzeln(modell_neu,modell_alt,index);
  korrektur();  // fullcondp[0]->posteriormode_const();
  newmodel(krit,mi,textit);
  fullcond_einzeln(modell_alt,modell_neu,index);
  //fullcondp[0]->posteriormode_const();
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


double STEPWISErun::df_ganzehatmatrix(void)
  {
  double df;
  if(fullcondp.size()>1)
    {
    unsigned i;
    if(fullcondp.size() != zcut.size() || fullcondp[0]->get_datanames().size() != xcut[1])
      {
      xcut.erase(xcut.begin(),xcut.end());
      zcut.erase(zcut.begin(),zcut.end());
      xcut.push_back(0);                          // enth�lt Struktur von Matrix X (unbestr. Anteil)
      zcut.push_back(0);                          // enth�lt Struktur von Matrix Z (bestr. Anteil)
      xcut.push_back(fullcondp[0]->get_datanames().size());
      for(i=1;i<fullcondp.size();i++)
        {
        xcut.push_back(xcut[i] + fullcondp[i]->get_dimX());  // gibt die h�chste Spalte von X an, die Eintr�ge von fullcond[i] enth�lt
        zcut.push_back(zcut[i-1] + fullcondp[i]->get_dimZ());
        }                                                 // Bem.: Dimension von "zcut" hier immer eins kleiner als bei "xcut"

      X = datamatrix(likep_mult[0]->get_nrobs(),xcut[xcut.size()-1],0);    // Initialisierung von X
      Z = datamatrix(likep_mult[0]->get_nrobs(),zcut[zcut.size()-1],0);    // Initialisierung von Z

      fullcond_alle[0]->createreml(X,Z,xcut[0],0);             // tr�gt fixe Effekte in X ein
      for(i=1;i<fullcondp.size();i++)
        {
        fullcondp[i]->createreml(X,Z,xcut[i],zcut[i-1]);   // schreibt Matrizen X und Z voll
        }
      }

    statmatrix<double>theta(zcut.size()-1,1,0);        // speichert Gl�ttungsparameter als Spalte ab, jeden Wert 1x
    for(i=1; i<fullcondp.size(); i++)
      {
      theta(i-1,0) = fullcondp[i]->get_lambda();
      }
    unsigned l,k;
    statmatrix<double>Qinv(Z.cols(),1,0);
    for(i=0, l=0; i<theta.rows(); i++)                  // Strafmatrix ist hier Einheitsmatrix
      {
      for(k=zcut[i]; k<zcut[i+1]; k++, l++)
        {
        Qinv(l,0) = theta(i,0);                         // Diagonale enth�lt dann dimZ * das jeweilige Lambda
        }
      }

    statmatrix<double>H(X.cols()+Z.cols(),X.cols()+Z.cols(),0);
    statmatrix<double>Hinv(H.rows(),H.rows(),0);
    H.weightedsscp2(X,Z,likep_mult[0]->get_weightiwls());       // berechnet (X,Z)'W(X,Z) und speichert es in Matrix H
    H.addtodiag(Qinv,X.cols(),X.cols()+Z.cols());             // berechnet H + Qinv und speichert es wieder in H
    Hinv = H.inverse();
    H.addtodiag(-Qinv,X.cols(),X.cols()+Z.cols());
    df = (H*Hinv).trace();
    }
  else
    df = fullcondp[0]->get_datanames().size();

  return df;
  }

// -----------------------------------------------------------------------------
// ------- Funktionen f�r die Erstellung des fullcondp-Vektors -----------------
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
     //fullcond_alle[i]->set_inthemodel(modell1[names_fixed.size()-2+i]);
     if(modell2[names_fixed.size()-2+i]==-1 && index==i)
        reset_fix(names_nonp[i-1][0]);
     if(modell1[names_fixed.size()-2+i] != -1 && modell1[names_fixed.size()-2+i] != 0)
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
fullcond_alle[index]->set_inthemodel(modell1[names_fixed.size()-2+index]);     

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
     if(m[names_fixed.size()-2+i] != 0 && m[names_fixed.size()-2+i] != -1)
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
  korrektur();  // fullcondp[0]->posteriormode_const();
  }
else //if(smoothing == "local")
  {
  unsigned i;
  for(i=1;i<fullcond_alle.size();i++)
    {
    fullcond_alle[i]->set_lambda_nr();
    fullcond_alle[i]->update_stepwise(m[names_fixed.size()-2+i]);
    }
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


void STEPWISErun::fix_ganz_komplett(const vector<double> &  modell)
  {
  if(smoothing == "global")
    {
    unsigned z;
    for(z=0;z<names_fixed.size()-1;z++)
       reset_fix(names_fixed[z+1]);

    for(z=0;z<names_fixed.size()-1;z++)
       {
       //if(modell[z]==0)
       //   reset_fix(names_fixed[z+1]);
       if(modell[z]==-1)
         {
         /*unsigned i = 1;
         bool rein = false;
         while(i<fullcondp[0]->get_datanames().size() && rein==false)
              {
              if(fullcondp[0]->get_datanames()[i]==names_fixed[z+1])
                rein = true;
              i = i + 1;
              }
         if(rein==false)*/
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
              // || fullcondp[0]->get_datanames()[j]== name+"_1")    // neu: Zusatz f�r VCM???
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


void STEPWISErun::korrektur(void)
  {
  //if(likep_mult[0]->get_family() == "Gamma")
  //  fullcond_alle[0]->set_effect_zero();
  //else
  //if(likep_mult[0]->get_family() != "Gaussian")
    fullcond_alle[0]->posteriormode_const();

/*unsigned i;
for(i=0;i<fullcondp.size();i++)
  fullcondp[i]->setbeta(0);*/
  }

// -----------------------------------------------------------------------------
// ------- Funktionen f�r die Ausgabe im Output-Fenster ------------------------
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
    outmodels << modeltext << endl << endl;
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
     fullcondp[i]->set_inthemodel(1);
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  OPTIONS FOR NONPARAMETRIC TERM: "
          + names_nonp[i-1][0] + "\n");
     genoptions_mult[0]->out("\n");
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
     fullcondp[i]->set_inthemodel(0);          
     }
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("STEPWISE PROCEDURE STARTED \n");
  genoptions_mult[0]->out("\n");
  }


// -----------------------------------------------------------------------------
// ------- Funktionen f�r die Ausgabe im Tex-File ------------------------------
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
  outtex << criterion.insert_string_char(hcharu,hstringu) << " = " << ST::doubletostring(kriterium_tex,6) << " \\\\ \n\\\\" << endl;
  }


void STEPWISErun::make_model(void)
  {
  //Vert-Fam wird �bergeben
  ST::string fam = likep_mult[0]->get_family();
  fam = fam.replaceallsigns('_', ' ');

  //Anz. Beob. wird �bergeben
  unsigned obs = likep_mult[0]->get_nrobs();

  //Name der Resp.-Var. �bergeben
  ST::string resp = likep_mult[0]->get_responsename();
  char hcharu = '_';
  ST::string hstringu = "\\_";
  resp = resp.insert_string_char(hcharu,hstringu);

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
  // falls andere Quantile gew�nscht werden
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

  // Schleife �berpr�ft, ob es ein fullcond-Object
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
    //erzeugt File, das Plot-Befehle f�r Java-Version enth�lt
    ofstream outbatch(path_batch.strtochar());

    //erzeugt File, das SPlus-Befehle zum Plotten enth�lt
    ofstream outsplus(path_splus.strtochar());

    outtex << "\n\\newpage" << "\n\\noindent {\\bf \\large Plots:}" << endl;
    outsplus << "# NOTE: 'directory' must be substituted by the directory"
             << " where the sfunctions are stored \n"
             << endl
    // einlesen der Source-Files f�r S-Plus
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
    for(j=begin[0];j<=end[0];j++)  //Schleife �berpr�ft, ob es map-Objekt gibt
      {
      plst = fullcondp[j]->get_plotstyle();
      if(plst == MCMC::drawmap || plst == MCMC::drawmapgraph)
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

    // falls andere Quantile gew�nscht werden
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

      // Plotstyle: noplot, plotnonp, drawmap, drawmapgraph
      plst = fullcondp[j]->get_plotstyle();

      if (plst != MCMC::noplot)
        {
        // Pfade f�r ps-, tex-, SPlus-files
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
          // f�r map-Funktionen
        else if (plst == MCMC::drawmap || plst == MCMC::drawmapgraph)
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
           if(plst == MCMC::drawmap)
             {
             outtex << "\n\\begin{figure}[h!]" << endl
                    << "\\centering" << endl
                    << "\\includegraphics[scale=0.6]{" << pathgr << "_pmean.ps}"
                    << endl
                    << "\\caption{Non--linear Effect of '" <<
                    regionvar.insert_string_char(hcharu,hstringu) << "'";
             outtex << ". Shown are the posterior means.}" << endl
                    << "\\end{figure}" << endl;
             }
           else if(plst == MCMC::drawmapgraph)
             {
             outtex << "\n%\\begin{figure}[h!]" << endl
                    << "%\\centering" << endl
                    << "%\\includegraphics[scale=0.6]{" << pathgr << "_pmean.ps}"
                    << endl
                    << "%\\caption{Non--linear Effect of '" <<
                    regionvar.insert_string_char(hcharu,hstringu) << "'";
             outtex << ". Shown are the posterior means.}" << endl
                    << "%\\end{figure}" << endl;
             }
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

  // f�r FG, die nicht durch ein lambda gebildet werden k�nnen:
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

  b = start_b_suchen(index_vec,krit_vec,df,b,z);     // sucht passenden Startwert f�r b!
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
      fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                  fullcond_alle[z]->get_data_forfixedeffects(),false);
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
    fullcond_alle[z]->set_inthemodel(-1);
    fullcond_alle[z]->reset_effect(0);
    fullcondp.erase(fullcondp.end()-1,fullcondp.end());
    fullcondp[0]->posteriormode_single(names_nonp[z-1],
                                fullcond_alle[z]->get_data_forfixedeffects(),true);
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
    fullcond_alle[z]->set_inthemodel(-1);
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
    fullcond_alle[z]->set_inthemodel(-1);
    fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                  fullcond_alle[z]->get_data_forfixedeffects(),true);
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
    fullcond_alle[z]->set_inthemodel(-1);
    vector<FULLCOND*> fullcond_start = fullcondp;
    vector<double> modell1 = modell_alt;
    modell1[z+names_fixed.size()-2] = -1;
    fullcond_einzeln(modell1,modell_alt,z);  // hier mu� der Fullcond-Vekor angepa�t werden!!!
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
    fullcond_einzeln(modell1,modell_alt,z);  // hier mu� der Fullcond-Vekor angepa�t werden!!!
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

void STEPWISErun::compute_average(ofstream & outweight)
  {
  vector<double> kriterien_alle;
  vector<double> priori;         // Versuch
  double prior = 0;
  vector<vector<double> > modelle;
  vector<ST::string> ausgabe;
  vector<int> vorgekommen;

  unsigned i;
  unsigned j;
  for(i=0;i<modell_neu.size();i++)    // fixe Effekte im Endmodell u.U. in falscher Reihenfolge!
    modell_neu[i] = 0;                // deswegen: zuerts Modell leer machen,
  fix_ganz_komplett(modell_neu);
  fullcond_komplett(modell_neu);

  fix_ganz_komplett(modell_alt);          // dann Endmodell mit richtiger Reihenfolge erzeugen!
  fullcond_komplett(modell_alt);
  posteriormode(posttitle,true);
  kriterium_alt = compute_criterion();
  genoptions_mult[0]->out("\n\n\n");
  genoptions_mult[0]->out("  Beginning of the Model_Averaging! \n");
  ST::string header = "  Final Model: ";
  ST::string text;
  maketext(header,modell_alt,kriterium_alt,text,true,trace,false);
  kriterien_alle.push_back(kriterium_alt);
  ausgabe.push_back(text);

BIC_min=1;
  //prior = pow(fullcond_alle[0]->compute_df(),0.25);
  //for(j=1;j<fullcondp.size();j++)
  //  prior *= 1/((fullcondp[j]->compute_df()+1)*(fullcondp[j]->compute_df()+1));
  prior=1;
  priori.push_back(prior);
  save_alle_betas(modell_alt);

  for(i=0;i<modell_alt.size();i++)
    {
    if(modell_alt[i] != 0)
      vorgekommen.push_back(1);
    else
      vorgekommen.push_back(0);
    }

  bool mcmcmc = false;      // sp�ter als Option!!!
  if(mcmcmc == true)
    mc3(kriterien_alle,priori,modelle,ausgabe,outweight,vorgekommen);
  else
    occam(kriterien_alle,priori,modelle,ausgabe,outweight,vorgekommen);

  for(i=0;i<names_fixed.size()-1;i++)    // sorgt daf�r, dass alle fixen Effekte und fullcond-Vektor stimmen,
    {
    if(vorgekommen[i] == 1)
      modell_neu[i] = -1;
    else
      modell_neu[i] = 0;
    }
  for(j=i;j<modell_neu.size();j++)       // d.h. fullcondp mu� alle Variablen enthalten, die in einem der Modelle vorgekommen sind!
    {
    if(vorgekommen[j] == 0)
      modell_neu[j] = 0;
    else
      {
      if(fullcond_alle[j-names_fixed.size()+2]->get_fctype()==factor)
        modell_neu[j] = -1;
      else
        modell_neu[j] = 1;
      }
    }
  fix_ganz_komplett(modell_neu);
  fullcond_komplett(modell_neu);
  for(i=1;i<fullcondp.size();i++)
    fullcondp[i]->average_posteriormode(kriterien_alle);
  fullcond_alle[0]->average_posteriormode(kriterien_alle);

  likep_mult[0]->posteriormode();       // berechnet "sigma2"
  genoptions_mult[0]->out("\n\n");
  likep_mult[0]->outresults();          // Ausgabe von "predictmean.raw", mit nicht standardidiertem Pr�diktor
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
    int helpint;
    if(modell[j+names_fixed.size()-2] == -1)
      {
      helpint = int(anzahl);
      fullcond_alle[j]->save_betas(modell,helpint);
      anzahl += fullcond_alle[j]->get_data_forfixedeffects().cols();
      }
    else if(modell[j+names_fixed.size()-2] == 0)
      {
      helpint = 0;
      fullcond_alle[j]->save_betas(modell,helpint);
      }
    else
      {
      helpint = -1;
      fullcond_alle[j]->save_betas(modell,helpint);   // ACHTUNG: hier steht -1 f�r nichtlinear!!!
      }
    }
  fullcond_alle[0]->save_betas2();
  }


void STEPWISErun::alle_modelle_berechnen(double z, vector<double> & hilf,
                      const vector<double> & best,
                      const vector<double> & unten, const vector<double> & oben,
                      vector<double> & kriterien_alle, vector<double> & priori,
                      vector<vector<double> > & modelle, vector<ST::string> & ausgabe)

  {
if(BIC_min == 30000)// 9227)
  fertig = true;

  #if defined(BORLAND_OUTPUT_WINDOW)
  Application->ProcessMessages();
  if (Frame->stop)
    {
    fertig = true;
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("MODEL AVERAGING TERMINATED BY USER BREAK\n");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("The remaining models won't be estimated!");
    genoptions_mult[0]->out("\n");

genoptions_mult[0]->out(" Anzahl gesch�tzter Modelle: " + ST::inttostring(BIC_min));
genoptions_mult[0]->out("\n");

    }
  #elif defined(JAVA_OUTPUT_WINDOW)
  bool stop = genoptions_mult[0]->adminb_p->breakcommand();
  if(stop)
    {
    fertig = true;
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("MODEL AVERAGING TERMINATED BY USER BREAK\n");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("The remaining models won't be estimated!");
    genoptions_mult[0]->out("\n");
    }
  #endif

  bool weitermachen = true;  //***

  while(hilf.size()>1 && fertig == false)
     {
     modell_neu[z+names_fixed.size()+hilf.size()-3] = lambdavec[hilf.size()-1][hilf[hilf.size()-1]];
     hilf.erase(hilf.end()-1,hilf.end());
     }

  if(hilf[0]<=oben[0] && fertig == false)
    {
    bool gefunden;
    double quotient, prior;
    unsigned index = best[0];                   
    gefunden = false;
    hilf[0] = index+1;
    while(hilf[0]>unten[0] && gefunden == false)
      {
      hilf[0] = hilf[0] - 1;
      modell_neu[z+names_fixed.size()-2] = lambdavec[0][hilf[0]];
      int ok = 0;
      unsigned j;
      for(j=0;j<modell_neu.size();j++)
        {
        if(modell_neu[j] != modell_alt[j])
           ok += 1;
        }
      if(ok>=2)     // bei ok==0 -> Startmodell, bei ok==1 -> Modelle entlang der Variablen!
        {
        fix_komplett(modell_neu);    // _ganz ?
        fullcond_komplett(modell_neu);
        fullcondp[0]->posteriormode_const();
        newmodel(kriterien_alle,modelle,ausgabe);

BIC_min+=1;        

        //prior = pow(fullcond_alle[0]->compute_df(),0.25);
        //for(j=1;j<fullcondp.size();j++)
        //  prior *= 1/((fullcondp[j]->compute_df()+1)*(fullcondp[j]->compute_df()+1));
        prior=1;

        quotient = priori[0]/(exp(-0.5*(kriterien_alle[kriterien_alle.size()-1]-kriterien_alle[0]))*prior);
        if(quotient > window)
          {
          gefunden = true;
          kriterien_alle.erase(kriterien_alle.end()-1,kriterien_alle.end());
          modelle.erase(modelle.end()-1,modelle.end());
          ausgabe.erase(ausgabe.end()-1,ausgabe.end());
          }
        else
          {
          priori.push_back(prior);
          save_alle_betas(modell_neu);
          }
        }
      }

    if(gefunden == true && hilf[0] == index)
      {
      hilf[0] = oben[0];
      weitermachen = false;
      }
    else if((gefunden == false || hilf[0] != index) && fertig == false)
      {
      hilf[0] = index;
      gefunden = false;
      }
    while(hilf[0]<oben[0] && gefunden == false && fertig == false)
      {
      hilf[0] = hilf[0] + 1;
      modell_neu[z+names_fixed.size()-2] = lambdavec[0][hilf[0]];
      int ok = 0;
      unsigned j;
      for(j=0;j<modell_neu.size();j++)
        {
        if(modell_neu[j] != modell_alt[j])
           ok += 1;
        }
      if(ok>=2)
        {
        fix_komplett(modell_neu);     // _ganz ?
        fullcond_komplett(modell_neu);
        fullcondp[0]->posteriormode_const();
        newmodel(kriterien_alle,modelle,ausgabe);

BIC_min+=1;        

        //prior = pow(fullcond_alle[0]->compute_df(),0.25);
        //for(j=1;j<fullcondp.size();j++)
        //  prior *= 1/((fullcondp[j]->compute_df()+1)*(fullcondp[j]->compute_df()+1));
        prior=1;

        quotient = priori[0]/(exp(-0.5*(kriterien_alle[kriterien_alle.size()-1]-kriterien_alle[0]))*prior);
        if(quotient > window)
          {
          gefunden = true;
          kriterien_alle.erase(kriterien_alle.end()-1,kriterien_alle.end());
          modelle.erase(modelle.end()-1,modelle.end());
          ausgabe.erase(ausgabe.end()-1,ausgabe.end());
          }
        else
          {
          priori.push_back(prior);
          save_alle_betas(modell_neu);
          }
        }
      }
    hilf[0] = oben[0];
    }

  if(hilf[0]>=oben[0] && fertig == false)
    {
    bool geaendert = false;
    if(hilf.size()==oben.size())
      {
      fertig = true;
      geaendert = true;
      }
    while(hilf.size()<oben.size() && geaendert==false && fertig == false)
      {
      bool lambda_exist;
      unsigned index = search_lambdaindex(modell_neu[z+names_fixed.size()+hilf.size()-2],
                                          lambdavec[hilf.size()],lambda_exist);

      hilf.push_back(index);                                          
      if(hilf.size()==oben.size())
        {
        if(hilf[hilf.size()-1] == oben[hilf.size()-1] && best[hilf.size()-1] != oben[hilf.size()-1])
          fertig = true;
        else if(hilf[hilf.size()-1] == unten[hilf.size()-1] && best[hilf.size()-1] == oben[hilf.size()-1])
          fertig = true;
        }
      if(hilf[hilf.size()-1]<oben[hilf.size()-1] ||
            (hilf[hilf.size()-1]==oben[hilf.size()-1] &&
            oben[hilf.size()-1] != unten[hilf.size()-1] &&
            oben[hilf.size()-1] == best[hilf.size()-1]))
        {
        geaendert = true;

        if(index == unten[hilf.size()-1])
          {
          if(best[hilf.size()-1] < oben[hilf.size()-1])
            {
            index = best[hilf.size()-1] + 1;
            modell_neu[z+names_fixed.size()+hilf.size()-3] = lambdavec[hilf.size()-1][index];
            }
          else
            {
            geaendert = false;
            hilf[hilf.size()-1] = best[hilf.size()-1];
            }
          }
        else
          {
          if(index > best[hilf.size()-1])
            {
            if(weitermachen == true)
              modell_neu[z+names_fixed.size()+hilf.size()-3] = lambdavec[hilf.size()-1][index+1];
            else
              {
              hilf[hilf.size()-1] = oben[hilf.size()-1];
              geaendert = false;
              weitermachen = true;
              }
            }
          else if(index <= best[hilf.size()-1] && best[hilf.size()-1] > unten[hilf.size()-1])
            {
            if(weitermachen == true)
              modell_neu[z+names_fixed.size()+hilf.size()-3] = lambdavec[hilf.size()-1][index-1];
            else
              {
              if(index != best[hilf.size()-1])
                {
                hilf[hilf.size()-1] = best[hilf.size()-1] + 1;
                modell_neu[z+names_fixed.size()+hilf.size()-3] = lambdavec[hilf.size()-1][best[hilf.size()-1] + 1];
                }
              else //if(hilf[hilf.size()-1] == best[hilf.size()-1])
                {
                hilf[hilf.size()-1] = oben[hilf.size()-1];
                geaendert = false;
                weitermachen = true;
                }
              }
            }
          else
            geaendert = false;
          }
        if(geaendert == true)
          {
          hilf.erase(hilf.end()-1,hilf.end());
          int j;
          for(j=hilf.size()-1;j>=1;j--)
            hilf[j] = best[j];
          hilf[0] = unten[0]-1;
          }
        }
      }
    }

  }


void STEPWISErun::occam(vector<double> & kriterien_alle, vector<double> & priori,
                      vector<vector<double> > & modelle, vector<ST::string> & ausgabe,
                      ofstream & outweight, vector<int> & vorgekommen)
  {
  unsigned i;
  double prior;
  for(i=1;i<names_fixed.size();i++)
    {
    modell_neu = modell_alt;
    if(modell_alt[i-1]==-1)
       modell_neu[i-1]= 0;
    else if(modell_alt[i-1]==0)
       {
       modell_neu[i-1] = -1;
       vorgekommen[i-1] = 1;
       }
    fix_ganz_komplett(modell_neu);
    fullcond_komplett(modell_neu);
    fullcondp[0]->posteriormode_const();
    newmodel(kriterien_alle,modelle,ausgabe);

    //prior = pow(fullcond_alle[0]->compute_df(),0.25);
    //for(j=1;j<fullcondp.size();j++)
    //  prior *= 1/((fullcondp[j]->compute_df()+1)*(fullcondp[j]->compute_df()+1));
    prior=1;
    priori.push_back(prior);
    save_alle_betas(modell_neu);
    }

  fix_ganz_komplett(modell_alt);      // macht das '_ganz' unten �berfl�ssig? 
  fullcond_komplett(modell_alt);

  unsigned z = 1;
  while(z<fullcond_alle.size() && fullcond_alle[z]->get_fctype()==factor)
     {
     modell_neu = modell_alt;
     if(modell_alt[z+names_fixed.size()-2]==-1 && fullcond_alle[z]->get_forced()==false)
        modell_neu[z+names_fixed.size()-2]= 0;
     else if(modell_alt[z+names_fixed.size()-2]==0)
        {
        modell_neu[z+names_fixed.size()-2] = -1;
        vorgekommen[z+names_fixed.size()-2] = 1;
        }
     if(modell_alt[z+names_fixed.size()-2]!=-1 || modell_neu[z+names_fixed.size()-2]!=-1)
        {
        fix_komplett(modell_neu);   // _ganz ?
        fullcond_komplett(modell_neu);
        fullcondp[0]->posteriormode_const();
        newmodel(kriterien_alle,modelle,ausgabe);

        //prior = pow(fullcond_alle[0]->compute_df(),0.25);
        //for(j=1;j<fullcondp.size();j++)
        //  prior *= 1/((fullcondp[j]->compute_df()+1)*(fullcondp[j]->compute_df()+1));
        prior=1;
        priori.push_back(prior);
        save_alle_betas(modell_neu);
        }
     z = z + 1;
     }

  vector<double> unten;
  vector<double> oben;
  vector<double> best;
  for(i=z;i<fullcond_alle.size();i++)
    {
    bool gefunden;
    double quotient;
    unsigned index = search_lambdaindex(modell_alt[names_fixed.size()-2+i],
                                        lambdavec[i-1],gefunden);
    best.push_back(index);
    gefunden = false;    // zum �berpr�fen, ob n�chstes Modell gut genug ist

    modell_neu = modell_alt;
    unsigned sch = 1;
    if(index == lambdavec[i-1].size()-1)
       oben.push_back(index);
    while(gefunden == false && index+sch < lambdavec[i-1].size())
      {
      modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][index+sch];
      fix_komplett(modell_neu);
      fullcond_komplett(modell_neu);
      fullcondp[0]->posteriormode_const();
      newmodel(kriterien_alle,modelle,ausgabe);

BIC_min+=1;

      //prior = pow(fullcond_alle[0]->compute_df(),0.25);
      //for(j=1;j<fullcondp.size();j++)
      //  prior *= 1/((fullcondp[j]->compute_df()+1)*(fullcondp[j]->compute_df()+1));
      prior=1;

      quotient = priori[0]/(exp(-0.5*(kriterien_alle[kriterien_alle.size()-1]-kriterien_alle[0]))*prior);
      if(quotient > window)    
        {
        gefunden = true;
        oben.push_back(index+sch-1);
        kriterien_alle.erase(kriterien_alle.end()-1,kriterien_alle.end());
        modelle.erase(modelle.end()-1,modelle.end());
        ausgabe.erase(ausgabe.end()-1,ausgabe.end());

BIC_min-=1;
        }
      else
        {
        priori.push_back(prior);
        save_alle_betas(modell_neu);
        if(index+sch == lambdavec[i-1].size()-1)
          {
          gefunden = true;
          oben.push_back(index+sch);
          }
        }
      sch++;
      }

    modell_neu = modell_alt;
    sch = 1;
    gefunden = false;
    if(index == 0)
       unten.push_back(index);
    while(gefunden== false && int(index-sch) >= 0)
      {
      modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][index-sch];
      fix_komplett(modell_neu);
      fullcond_komplett(modell_neu);
      fullcondp[0]->posteriormode_const();
      newmodel(kriterien_alle,modelle,ausgabe);

BIC_min+=1;

      //prior = pow(fullcond_alle[0]->compute_df(),0.25);
      //for(j=1;j<fullcondp.size();j++)
      //  prior *= 1/((fullcondp[j]->compute_df()+1)*(fullcondp[j]->compute_df()+1));
      prior=1;

      quotient = priori[0]/(exp(-0.5*(kriterien_alle[kriterien_alle.size()-1]-kriterien_alle[0]))*prior);
      if(quotient > window)    
        {
        gefunden = true;
        unten.push_back(index-sch+1);
        kriterien_alle.erase(kriterien_alle.end()-1,kriterien_alle.end());
        modelle.erase(modelle.end()-1,modelle.end());
        ausgabe.erase(ausgabe.end()-1,ausgabe.end());

BIC_min-=1;
        }
      else
        {
        vorgekommen[names_fixed.size()-2+i] = 1;  // df = 0 ist der h�chste Index. Deshalb kann �nderung "raus" -> "rein" nur hier vorkommen!
        priori.push_back(prior);
        save_alle_betas(modell_neu);
        if(index-sch == 0)
          {
          gefunden = true;
          unten.push_back(index-sch);
          }
        }
      sch++;
      }
    }

  unsigned anzahl_modelle = 1;
  fertig = false;
  unsigned j;
  for(j=0;j<unten.size();j++)
    {
genoptions_mult[0]->out(ST::inttostring(oben[j]) + "   " + ST::inttostring(unten[j]));
    anzahl_modelle *= (oben[j]-unten[j]+1);
    }
  if(anzahl_modelle == 1)
    {
    anzahl_modelle = 0;
    fertig = true;
    }
  genoptions_mult[0]->out("\n\n");
  genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
  genoptions_mult[0]->out("  Number of models for model_averaging: " + ST::inttostring(anzahl_modelle) + "\n");
  genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
  genoptions_mult[0]->out("\n\n");

  modell_neu = modell_alt;

  //vector<double> hilf = unten;  
  vector<double> hilf = best;
  // bool nachoben = true;     // f�r "zweitesmal = true"
  // bool nachunten = true;    // f�r "zweitesmal = false"
  while(fertig == false)
    alle_modelle_berechnen(z,hilf,best,unten,oben,kriterien_alle,priori,modelle,ausgabe);

  int l;         // Gewichte berechnen!
  double kriterium_summe = 0;
  vector<double> weights_alle = kriterien_alle;

  outweight << "no  " << criterion << "  weight  model" << endl << endl;
  for(l=kriterien_alle.size()-1;l>=0;l--)
     {
     kriterien_alle[l] = exp(-0.5*(kriterien_alle[l]-kriterien_alle[0])) * priori[l];   // es ist egal, ob kriterien_alle[0] wirklich das Minimum ist!!!
     kriterium_summe += kriterien_alle[l];
     }
  for(i=0;i<kriterien_alle.size();i++)
     {
     kriterien_alle[i] /= kriterium_summe;
     outweight << ST::inttostring(i) << "  " << ST::doubletostring(weights_alle[i],6)
               << "  " << ST::doubletostring(kriterien_alle[i],6)
               << "  " << ausgabe[i] << endl << endl;

     }

  }


void STEPWISErun::mc3(vector<double> & kriterien_alle, vector<double> & priori,
                      vector<vector<double> > & modelle, vector<ST::string> & ausgabe,
                      ofstream & outweight, vector<int> & vorgekommen)
  {

  // kriterium_alt, modell_alt enthalten immer die Sachen vom alten Modell,
  // kriterium_neu, modell_neu immer die Sachen vom neuen Modell
  modell_neu = modell_alt;
  kriterium_neu = kriterium_alt;

  unsigned i;
  unsigned iterations = 10000;  // sp�ter als Option!
  double var, u, alpha;
  double prior, quotient;
  for(i=0;i<iterations;i++)
    {
    var = floor(uniform()*modell_alt.size());
    if(var == modell_alt.size())  // hat Wkt = 0!
      var = var - 1;

    // Abfrage, ob fixer Effekt, Faktorvariable oder nichtparametrisch
    if(var < names_fixed.size()-1)          // Fixer Effekt
      {
      u = uniform();            // u <= 0.5 hei�t: Modell bleibt gleich
                                // u > 0.5 hei�t: Modell wird ver�ndert
      if(u > 0.5)
        {
        if(modell_alt[var]==-1)
          modell_neu[var]= 0;
        else if(modell_alt[var]==0)
          modell_neu[var] = -1;
        }
      else //if(u <= 0.5)
        u = 0;
      }
    else if(fullcond_alle[var-names_fixed.size()+2]->get_fctype()==factor)  // Faktor
      {
      if(fullcond_alle[var-names_fixed.size()+2]->get_forced()==false)
        u = uniform();            // u <= 0.5 hei�t: Modell bleibt gleich
      else                        // u > 0.5 hei�t: Modell wird ver�ndert
        u = 0;

      if(u > 0.5)
        {
        if(modell_alt[var]==-1)
          modell_neu[var]= 0;
        else if(modell_alt[var]==0)
          modell_neu[var] = -1;
        }
      else //if(u <= 0.5)
        u = 0;
      }
    else  // nichtparametrisch
      {
      bool gefunden;
      unsigned index = search_lambdaindex(modell_alt[var],
                                        lambdavec[var-names_fixed.size()+1],gefunden);
      if(index+1 < lambdavec[var-names_fixed.size()+1].size() && int(index)-1 >= 0)
        u = floor(uniform()*3) - 1;
      else if(index+1 >= lambdavec[var-names_fixed.size()+1].size() && int(index)-1 >= 0)
        u = floor(uniform()*2) - 1;
      else if(index+1 < lambdavec[var-names_fixed.size()+1].size() && int(index)-1 < 0)
        u = floor(uniform()*2);
      else
        u = 0;

      if(u == -1)
        modell_neu[var] = lambdavec[var-names_fixed.size()+1][index-1];
      else if(u == 1)
        modell_neu[var] = lambdavec[var-names_fixed.size()+1][index+1];
      }

    if(u == 0)  // Fall, da� Modell gar nicht ver�ndert wird!
      {
      kriterien_alle.push_back(kriterium_alt);
      priori.push_back(priori[priori.size()-1]);
      save_alle_betas(modell_alt);
      ST::string text = ausgabe[ausgabe.size()-1];
      maketext("  Previous Model:",modell_alt,kriterium_alt,text,false,trace,false);
      ausgabe.push_back(text);
      }
    else
      {
      fix_ganz_komplett(modell_neu);
      fullcond_komplett(modell_neu);
      fullcondp[0]->posteriormode_const();
      newmodel(kriterien_alle,modelle,ausgabe);

      kriterium_neu = kriterien_alle[kriterien_alle.size()-1];
      //prior = pow(fullcond_alle[0]->compute_df(),0.25);
      //for(j=1;j<fullcondp.size();j++)
      //  prior *= 1/((fullcondp[j]->compute_df()+1)*(fullcondp[j]->compute_df()+1));
      prior = 1;
      quotient = priori[priori.size()-1]/(exp(-0.5*(kriterium_alt-kriterium_neu))*prior);
      if(quotient > 1)
        quotient = 1;

      alpha = uniform();
      if(alpha > quotient)   // neues Modell wird nicht akzeptiert --> alten Zustand wieder herstellen
        {
        if(trace == "trace_minim" || trace == "trace_on")
          genoptions_mult[0]->out("  Acceptance rate = " + ST::doubletostring(quotient) + "; Model was not accepted!");

        modell_neu = modell_alt;
        kriterium_neu = kriterium_alt;
        kriterien_alle.erase(kriterien_alle.end()-1,kriterien_alle.end());
        modelle.erase(modelle.end()-1,modelle.end());
        ausgabe.erase(ausgabe.end()-1,ausgabe.end());

        fix_ganz_komplett(modell_alt);
        fullcond_komplett(modell_alt);
        fullcondp[0]->posteriormode_const();
        newmodel(kriterien_alle,modelle,ausgabe);
        prior = priori[priori.size()-1];
        }
      else
        {
        if(trace == "trace_minim" || trace == "trace_on")
          genoptions_mult[0]->out("  Acceptance rate = " + ST::doubletostring(quotient) + "; Model was accepted!");

        if(modell_neu[var] != 0)
          vorgekommen[var] = 1;
        }
        
      priori.push_back(prior);
      save_alle_betas(modell_neu);
      modell_alt = modell_neu;
      kriterium_alt = kriterium_neu;
      }

    }   // Ende: for(i=0;i<iterations;i++)

  outweight << "no  " << criterion << "  model" << endl << endl;
  for(i=0;i<kriterien_alle.size();i++)
     {
     //if(i==0)
       outweight << ST::inttostring(i) << "  " << ST::doubletostring(kriterien_alle[i],6)
                 << "  " << endl << endl;
     /*outweight << ST::inttostring(i) << "  " << ST::doubletostring(kriterien_alle[i],6)
                 << "  " << ausgabe[i]
                 << endl << endl; */
     kriterien_alle[i] = 1/(double(iterations)+1);
     }
  }



}





