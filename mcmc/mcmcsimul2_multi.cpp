 // DATE: 01.02.99

//---------------------------------------------------------------------------
#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include<StatwinFrame.h>
#endif

#include<mcmcsimul2_multi.h>
#include<time.h>
#include<clstring.h>
#include <stdlib.h>
#include<math.h>


namespace MCMC
{


STEPMULTIrun::STEPMULTIrun(MCMCoptions * go,DISTRIBUTION * dp,
vector<FULLCOND*> & fc)
: MCMCsimulate(go,dp,fc)
  {
  }


STEPMULTIrun::STEPMULTIrun(const STEPMULTIrun & st)
  : MCMCsimulate(MCMCsimulate(st))
  {
  }


const STEPMULTIrun & STEPMULTIrun::operator=(
const STEPMULTIrun & st)
  {
  if (this==&st)
    return *this;
  MCMCsimulate::operator=(MCMCsimulate(st));
  return *this;
  }


// -----------------------------------------------------------------------------
// ------- die grundlegenden Funktionen ----------------------------------------
// -----------------------------------------------------------------------------

bool STEPMULTIrun::stepwise(const ST::string & procedure, const ST::string & minimum,
        const ST::string & crit, const int & stp, const ST::string & trac,
        const int & number, const ST::string & stam, const int & inc, const bool & finet,
        const bool & fineloc, const bool & maveraging, int & fenster,
        const datamatrix & Da, const vector<ST::string> & modelvar,
        const ST::string & name, vector<FULLCOND*> & fullcond_z, ST::string & path,
        const bool & CI, bool & hier)
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
  finelocal = fineloc;
  trace = trac;
  bool modelaveraging = maveraging;
  window = fenster;
  smoothing = "global";
  hierarchical = hier;

  fullcond_alle = fullcondp;

  kategorien = likep_mult[0]->get_responsedim();
  katje = 0;
  anz_fullcond = unsigned(fullcond_alle.size()/kategorien);

  modell_alt.erase(modell_alt.begin(),modell_alt.end());
  vector<double> modell_final;
  double kriterium_final;
  ST::string text_final;

  ST::string tr_akt = "trace_on";
  posttitle.push_back("");

  lambdavec.erase(lambdavec.begin(),lambdavec.end());
  names_fixed.erase(names_fixed.begin(),names_fixed.end());
  names_nonp.erase(names_nonp.begin(),names_nonp.end());

  set_center(likep_mult[0],fullcond_alle,begin[0],end[0]);  // sorgt dafür, daß Funktionen zentriert werden!

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
  //outtex.open(path_tex.strtochar());
  //make_graphics(name,startindex);

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

     if(minim=="adaptiv")
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
     genoptions_mult[0]->out("  NOTE: This option is not available with this response distribution!");
  if(finelocal == true)
     genoptions_mult[0]->out("  NOTE: This option is not available with this response distribution!");

  if(modelaveraging == true)
     genoptions_mult[0]->out("  NOTE: This option is not available with this response distribution!");

  if(CI == true)
     genoptions_mult[0]->out("  NOTE: This option is not available with this response distribution!");

  fullcond_z = fullcondp;
  for(i=0;i<fullcond_z.size();i++)
    fullcond_z[i]->set_fcnumber(i);

  posteriormode(posttitle,false); // Problem: linearer Prädiktor bei "true" standardisiert! Hier wird zurückgerechnet!
                                    // danach nicht mehr compute_criterion() aufrufen!!!

  //make_tex_end(path,modell_final);

  // Files müssen wieder geschlossen werden!!!
  //outtex.close();
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


bool STEPMULTIrun::single_stepwise(const vector<unsigned> & start,
                            const vector<double> & startfix, const bool & tex)
  {
  modell_neu.erase(modell_neu.begin(),modell_neu.end());
  modellematrix.erase(modellematrix.begin(),modellematrix.end());
  steps_aktuell = 0;
  ST::string tr_akt = "trace_on";
  vector<vector<double> > startiteration;

  unsigned i;
  for(katje=0;katje<kategorien;katje++)
    {
    for(i=0;i<names_fixed.size()-1;i++)
       modell_neu.push_back(startfix[i]);

    for(i=0;i<anz_fullcond-1;i++)
       {
       double lambda = lambdavec[i][start[i]];
       modell_neu.push_back(lambda);
       }
    }
  katje=0;

  modell_alt = modell_neu;
  startiteration.push_back(modell_alt);
  modellematrix.push_back(startiteration);
  fix_komplett(modell_alt);
  fullcond_komplett(modell_alt);

  /*if(hierarchical == true)         // Abfrage, ob Startmodell hierarchisch ist! 
    {
    for(i=anz_fullcond-1;i>=1;i--)   // Hier müßte die Abfrage für die erste Kategorie ausreichen
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
       if(modell_alt[names_fixed.size()-2+i] == 0 && possible == "spfix")
         falsch = false;

       if(falsch == false)
         {
         genoptions_mult[0]->out("  NOTE: The startmodel is no hierarchical model! Please choose another one!);
         return true;
         }
       }
    } */

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
    {
    abbruch = koordabstieg();

    /*if(minim=="adaptiv")
      {
      while(modell_neu != modellematrix[0][0])
        {
        fertig = false;
        modell_neu = modell_alt;
        startiteration.erase(startiteration.begin(),startiteration.end());
        startiteration.push_back(modell_alt);
        modellematrix.erase(modellematrix.begin(),modellematrix.end());
        modellematrix.push_back(startiteration);
        posteriormode(posttitle,true);
        kriterium_alt = compute_criterion();
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



// -----------------------------------------------------------------------------
// ------------------ Funktionen für Stepmin -----------------------------------
// -----------------------------------------------------------------------------


bool STEPMULTIrun::stepfunctions(void)
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

       for(katje=0;katje<kategorien;katje++)
         {
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
              unsigned ind_fullc;
              if(katje==0)
                ind_fullc = (kategorien-1)*anz_fullcond;
              else
                ind_fullc = (katje-1)*anz_fullcond;
              fullcond_alle[ind_fullc]->posteriormode_const();
              posteriormode(posttitle,true);
              step_minfix(kriteriumiteration2,modeliteration,textiteration);
              unsigned z = step_minfactor(kriteriumiteration2,modeliteration,textiteration);
              stepmin_nonp(kriteriumiteration2,modeliteration,textiteration,z);
              }
           }    // End: else if(algorithm == "stepmin")

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

           if(kriterium_alt >= kriterium_neu)
             {
             fix_komplett(modell_neu);
             fullcond_komplett(modell_neu);
             if(katje < kategorien-1)
               {
               header = "\n\nStart Model for the next category:";
               neutext = false;
               maketext(header,modell_neu,kriterium_neu,text_neu,neutext,
                        trace,true);                    // Modell rausschreiben!!
               }
             modell_alt = modell_neu;
             kriterium_alt = kriterium_neu;
             text_alt = text_neu;
             }
           else  // if(kriterium_alt < kriterium_neu)
             {
             fix_komplett(modell_alt);
             fullcond_komplett(modell_alt);
             if(katje < kategorien-1)
               {
               header = "\n\nStart Model for the next category:";
               neutext = false;
               maketext(header,modell_alt,kriterium_alt,text_alt,neutext,
                        trace,true);                    // Modell rausschreiben!!
               }
             kriterium_neu = kriterium_alt;
             modell_neu = modell_alt;
             text_neu = text_alt;
             }
           }
         else
           {
           if(trace == "trace_on" || trace == "trace_minim")
             {
             genoptions_mult[0]->out("\n\n");
             genoptions_mult[0]->out("  There are no new models for this category! \n");
             }
           }
         }

       if(fertig == false)
         {
         header = "\n\nBest Model of this iteration:";
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


unsigned STEPMULTIrun::stepwise_fixfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration)
  {
  unsigned i;
  unsigned unten = katje*(anz_fullcond + names_fixed.size()-2);
  unsigned oben = katje*(anz_fullcond + names_fixed.size()-2) + names_fixed.size()-1;
  unsigned ind_name;

  for(i=unten;i<oben;i++)
    {
    ind_name = i + 1 - katje*(anz_fullcond + names_fixed.size()-2);;
    modell_neu = modell_alt;
    if(modell_alt[i]==-1)
      modell_neu[i]= 0;
    else if(modell_alt[i]==0)
      modell_neu[i] = -1;
    if(modelcomparison(modell_neu,modellematrix)==false)
      newmodel_fix(modell_neu[i],kriteriumiteration2,modeliteration,
                            textiteration,names_fixed[ind_name]);
    }

  unsigned z = 1 + katje*anz_fullcond;
  ind_name = z - katje*anz_fullcond - 1;
  unsigned ind_mod = z + (katje+1)*(names_fixed.size()-2);
  while(z<fullcond_alle.size() && fullcond_alle[z]->get_fctype()==factor)
     {
     modell_neu = modell_alt;
     if(modell_alt[ind_mod]==-1 && fullcond_alle[z]->get_forced()==false)
       modell_neu[ind_mod]= 0;
     else if(modell_alt[ind_mod]==0)
       modell_neu[ind_mod] = -1;
     if(modelcomparison(modell_neu,modellematrix)==false)
       newmodel_factor(modell_neu[ind_mod],z,kriteriumiteration2,
                 modeliteration,textiteration,names_nonp[ind_name]);
     z = z + 1;
     }
  return z;
  }


void STEPMULTIrun::stepwise_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  unsigned i;
  for(i=z;i<((katje+1)*anz_fullcond);i++)
    {
    unsigned ind_mod = i + (katje+1)*(names_fixed.size()-2);
    unsigned ind_lamb = i - (katje+1)*1;

    ST::string possible = "alles";
    if(hierarchical == true)
      fullcond_alle[i]->hierarchical(possible);

    unsigned sch;
    for(sch=1;sch<=unsigned(increment);sch++)
       {
       modell_neu = modell_alt;
       bool lambda_exist;    // zum Überprüfen, ob neues Lambda im Vektor enthalten ist
       unsigned index = search_lambdaindex(modell_alt[ind_mod],
                                lambdavec[ind_lamb],lambda_exist);
       lambda_exist = false;
       if(index < lambdavec[ind_lamb].size()-sch)
          lambda_exist = true;
       if(lambda_exist==true && hierarchical == true)
         {
         if(lambdavec[ind_lamb][index+sch] == 0 && (possible == "spline" || possible == "spfix"))
           lambda_exist = false;
         if(lambdavec[ind_lamb][index+sch] == -1 && (possible == "spline" || possible == "raus"))
           lambda_exist = false;
         if(lambdavec[ind_lamb][index+sch] > 0 && (possible == "rfix" || possible == "raus"))
           lambda_exist = false;
         }
       if(lambda_exist==true)
          {
          modell_neu[ind_mod] = lambdavec[ind_lamb][index+sch];
          if(modelcomparison(modell_neu,modellematrix)==false)
            newmodel_nonp(i,kriteriumiteration2,modeliteration,textiteration);
          }

       lambda_exist = false;
       modell_neu = modell_alt;
       if(index >= sch)
          lambda_exist = true;
       if(lambda_exist==true && hierarchical == true)
         {
         if(lambdavec[ind_lamb][index-sch] == 0 && (possible == "spline" || possible == "spfix"))
           lambda_exist = false;
         if(lambdavec[ind_lamb][index-sch] == -1 && (possible == "spline" || possible == "raus"))
           lambda_exist = false;
         if(lambdavec[ind_lamb][index-sch] > 0 && (possible == "rfix" || possible == "raus"))
           lambda_exist = false;
         }
       if(lambda_exist==true)
          {
          modell_neu[ind_mod] = lambdavec[ind_lamb][index-sch];
          if(modelcomparison(modell_neu,modellematrix)==false)
            newmodel_nonp(i,kriteriumiteration2,modeliteration,textiteration);
          }
       }
    }
  }


void STEPMULTIrun::stepmin_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  unsigned i;
  for(i=z;i<((katje+1)*anz_fullcond);i++)
    {
    unsigned ind_mod = i + (katje+1)*(names_fixed.size()-2);
    unsigned ind_fullc = katje*anz_fullcond;
    unsigned ind_lamb = i - (katje+1)*1;

    modell_neu = modell_alt;
    unsigned lambda_ind;

    unsigned y;
    for(y=1;y<fullcondp.size();y++)
      {
      if(fullcondp[y] == fullcond_alle[i])
        fullcondp[y]->remove_centering();
      }
    for(y=1;y<fullcond_alle.size();y++)
      fullcond_alle[y]->set_center(false);    // sorgt dafür, daß Funktionen nicht zentriert werden!

    vector<double> krit_fkt;
    if(modell_alt[ind_mod]==0)
      stepmin_nonp_leer(i,krit_fkt,kriterium_alt);
    else if(modell_alt[ind_mod]==-1)
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

    for(y=1;y<fullcond_alle.size();y++)
     {
     if(fullcond_alle[y]->is_identifiable() == false)
       fullcond_alle[y]->set_center(true);    // sorgt dafür, daß Funktionen zentriert werden!
     }

    modell_neu[ind_mod] = lambdavec[ind_lamb][lambda_ind];
    if(modell_neu[ind_mod] != modell_alt[ind_mod])
      {
      if(modelcomparison(modell_neu,modellematrix)==false)
        {
        newmodel_nonp(i,kriteriumiteration2,modeliteration,textiteration);
        // Stellt linearen Prädiktor wieder her. Besser wäre, den lin. Prädiktor zu speichern!!!
        fullcond_alle[ind_fullc]->posteriormode_const();
        posteriormode(posttitle,true);
        }
      }
    }
  }


void STEPMULTIrun::minexact_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  unsigned i;
  unsigned j;
  for(i=z;i<((katje+1)*anz_fullcond);i++)
    {
    unsigned ind_mod = i + (katje+1)*(names_fixed.size()-2);
    unsigned ind_fullc = katje*anz_fullcond;
    unsigned ind_lamb = i - (katje+1)*1;
    unsigned ind_name = i - katje*anz_fullcond - 1;

    modell_neu = modell_alt;
    unsigned lambda_ind;

    vector<double> krit_fkt;
    double df = 0;
    if(modell_alt[ind_mod]==0)
      {
      for(j=0;j<fullcondp.size();j++)
        df = df + fullcondp[j]->compute_df();
      fullcondp.push_back(fullcond_alle[i]);
      minexact_nonp_leer(i,krit_fkt,kriterium_alt,df);
      }
    else if(modell_alt[ind_mod]==-1)
      {
      for(j=0;j<fullcondp.size();j++)
        df = df + fullcondp[j]->compute_df();
      df = df - 1;
      reset_fix(names_nonp[ind_name][0]);
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

    modell_neu[ind_mod] = lambdavec[ind_lamb][lambda_ind];
    if(modell_neu[ind_mod] != modell_alt[ind_mod])
      {
      if(modelcomparison(modell_neu,modellematrix)==false)
        {
        newmodel_nonp(i,kriteriumiteration2,modeliteration,textiteration);
        // Stellt linearen Prädiktor wieder her.
        fullcond_alle[ind_fullc]->posteriormode_const();
        posteriormode(posttitle,true);
        }
      }
    }
  }


// -----------------------------------------------------------------------------
// ------------------ Funktionen für Stepmin -----------------------------------
// -----------------------------------------------------------------------------

void STEPMULTIrun::step_minfix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration)
  {
  unsigned i;
  unsigned unten = katje*(anz_fullcond + names_fixed.size()-2);
  unsigned oben = katje*(anz_fullcond + names_fixed.size()-2) + names_fixed.size()-1;
  unsigned ind;

  for(i=unten;i<oben;i++)
    {
    ind = i+1;
    if(modell_alt[i]==-1)
      stepmin_fix_leer(kriteriumiteration2,modeliteration,textiteration,ind);
    else if(modell_alt[i]==0)
      stepmin_leer_fix(kriteriumiteration2,modeliteration,textiteration,ind);
    }
  }


void STEPMULTIrun::stepmin_fix_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration, unsigned & i)
  {
  double df = 0;
  unsigned y;
  for(y=0;y<fullcondp.size();y++)
    df = df + fullcondp[y]->compute_df();
  df = df - 1;

  unsigned ind_name = i - katje*(anz_fullcond + names_fixed.size()-2);
  unsigned ind_fullc = katje*anz_fullcond;

  fullcond_alle[ind_fullc]->safe_const();
  reset_fix(names_fixed[ind_name]);
  fullcond_alle[ind_fullc]->posteriormode_const();
  kriterium_neu = criterion_min(df);
  fullcond_alle[ind_fullc]->set_const_old();
  if(minim == "approx_control")
    {
    posteriormode(posttitle,true);
    double kriterium_test = criterion_min(df);
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_fixed[ind_name] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    include_fix(names_fixed[ind_name]);
    fullcond_alle[ind_fullc]->posteriormode_const();
    posteriormode(posttitle,true);
    reset_fix(names_fixed[ind_name]);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_fixed[ind_name] + "\n");
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
      include_fix(names_fixed[ind_name]);
      fullcondp[0]->posteriormode_const();
      posteriormode(posttitle,true);
      }
    else
      {
      int c = column_for_fix(names_fixed[ind_name]);
      vector<ST::string> name_help;
      name_help.push_back(names_fixed[ind_name]);
      fullcond_alle[ind_fullc]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
      }
    modell_neu[i-1] = -1;
    }
  else
    {
    int c = column_for_fix(names_fixed[ind_name]);
    vector<ST::string> name_help;
    name_help.push_back(names_fixed[ind_name]);
    fullcond_alle[ind_fullc]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
    }
  }


void STEPMULTIrun::stepmin_leer_fix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration, unsigned & i)
  {
  double df = 0;
  unsigned y;
  for(y=0;y<fullcondp.size();y++)
    df = df + fullcondp[y]->compute_df();
  df = df + 1;

  unsigned ind_name = i - katje*(anz_fullcond + names_fixed.size()-2);
  unsigned ind_fullc = katje*anz_fullcond;

  fullcond_alle[ind_fullc]->safe_const();
  int c = column_for_fix(names_fixed[ind_name]);
  vector<ST::string> name_help;
  name_help.push_back(names_fixed[ind_name]);
  fullcond_alle[ind_fullc]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
  kriterium_neu = criterion_min(df);
  fullcond_alle[ind_fullc]->set_const_old();
  if(minim == "approx_control")
    {
    posteriormode(posttitle,true);
    double kriterium_test = criterion_min(df);
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_fixed[ind_name] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    reset_fix(names_fixed[ind_name]);
    fullcond_alle[ind_fullc]->posteriormode_const();
    posteriormode(posttitle,true);
    fullcond_alle[ind_fullc]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_fixed[ind_name] + "\n");
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
       reset_fix(names_fixed[ind_name]);
       fullcond_alle[ind_fullc]->posteriormode_const();
       posteriormode(posttitle,true);
       }
     else
       reset_fix(names_fixed[ind_name]);
     modell_neu[i-1] = 0;
     }
  else
     reset_fix(names_fixed[ind_name]);
  }


unsigned STEPMULTIrun::step_minfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration)
  {

  unsigned z = 1 + katje*anz_fullcond;
  unsigned ind_mod = z + (katje+1)*(names_fixed.size()-2);

  while(z<((katje+1)*anz_fullcond) && fullcond_alle[z]->get_fctype()==factor)
     {
     if(modell_alt[ind_mod]==-1 && fullcond_alle[z]->get_forced()==false)
       stepmin_factor_leer(kriteriumiteration2,modeliteration,textiteration,z);
     else if(modell_alt[ind_mod]==0)
       stepmin_leer_factor(kriteriumiteration2,modeliteration,textiteration,z);
     z = z + 1;
     }
  return z;
  }


void STEPMULTIrun::stepmin_factor_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  double df = 0;
  unsigned i;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();
  df = df - fullcond_alle[z]->get_data_forfixedeffects().cols();

  unsigned ind_name = z - katje*anz_fullcond - 1;
  unsigned ind_fullc = katje*anz_fullcond;
  unsigned ind_mod = z + (katje+1)*(names_fixed.size()-2);

  fullcond_alle[ind_fullc]->safe_const();
  for(i=0;i<names_nonp[ind_name].size();i++)
    reset_fix(names_nonp[ind_name][i]);
  fullcond_alle[ind_fullc]->posteriormode_const();
  kriterium_neu = criterion_min(df);
  fullcond_alle[ind_fullc]->set_const_old();
  if(minim == "approx_control")
    {
    posteriormode(posttitle,true);
    double kriterium_test = criterion_min(df);
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    fullcond_alle[ind_fullc]->include_effect(names_nonp[ind_name],fullcond_alle[z]->get_data_forfixedeffects());
    fullcond_alle[ind_fullc]->posteriormode_const();
    posteriormode(posttitle,true);
    for(i=0;i<names_nonp[ind_name].size();i++)
      reset_fix(names_nonp[ind_name][i]);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_alt,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(kriterium_neu < kriterium_alt)
     {
     modell_neu[ind_mod] = 0;
     if(modelcomparison(modell_neu,modellematrix)==false)
       {
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       fullcond_alle[ind_fullc]->include_effect(names_nonp[ind_name],fullcond_alle[z]->get_data_forfixedeffects());
       fullcond_alle[ind_fullc]->posteriormode_const();
       posteriormode(posttitle,true);
       }
     else
       fullcond_alle[ind_fullc]->posteriormode_single(names_nonp[ind_name],
                                  fullcond_alle[z]->get_data_forfixedeffects(),true);
     modell_neu[ind_mod] = -1;
     }
  else
     fullcond_alle[ind_fullc]->posteriormode_single(names_nonp[ind_name],
                                 fullcond_alle[z]->get_data_forfixedeffects(),true);
  }


void STEPMULTIrun::stepmin_leer_factor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  double df = 0;
  unsigned i;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();
  df = df + fullcond_alle[z]->get_data_forfixedeffects().cols();

  unsigned ind_name = z - katje*anz_fullcond - 1;
  unsigned ind_fullc = katje*anz_fullcond;
  unsigned ind_mod = z + (katje+1)*(names_fixed.size()-2);

  fullcond_alle[ind_fullc]->safe_const();
  fullcond_alle[ind_fullc]->posteriormode_single(names_nonp[ind_name],
                                 fullcond_alle[z]->get_data_forfixedeffects(),true);
  kriterium_neu = criterion_min(df);
  fullcond_alle[ind_fullc]->set_const_old();
  if(minim == "approx_control")
    {
    posteriormode(posttitle,true);
    double kriterium_test = criterion_min(df);
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    for(i=0;i<names_nonp[ind_name].size();i++)
      reset_fix(names_nonp[ind_name][i]);
    fullcond_alle[ind_fullc]->posteriormode_const();
    posteriormode(posttitle,true);
    fullcond_alle[ind_fullc]->posteriormode_single(names_nonp[ind_name],
                                 fullcond_alle[z]->get_data_forfixedeffects(),true);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_alt,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(kriterium_neu < kriterium_alt)
     {
     modell_neu[ind_mod] = -1;
     if(modelcomparison(modell_neu,modellematrix)==false)
       {
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       for(i=0;i<names_nonp[ind_name].size();i++)
         reset_fix(names_nonp[ind_name][i]);
       fullcond_alle[ind_fullc]->posteriormode_const();
       posteriormode(posttitle,true);
       }
     else
       {
       for(i=0;i<names_nonp[ind_name].size();i++)
         reset_fix(names_nonp[ind_name][i]);
       }
     modell_neu[ind_mod] = 0;
     }
  else
     {
     for(i=0;i<names_nonp[ind_name].size();i++)
       reset_fix(names_nonp[ind_name][i]);
     }
  }


void STEPMULTIrun::stepmin_nonp_nonp(unsigned & z, vector<double> & krit_fkt,double & kriterium)
  {
  if(smoothing == "local")
    fullcond_alle[z]->set_lambda_nr();

  unsigned i;
  double df = 0;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();

  unsigned ind_name = z - katje*anz_fullcond - 1;
  unsigned ind_fullc = katje*anz_fullcond;
  unsigned ind_mod = z + (katje+1)*(names_fixed.size()-2);
  unsigned ind_lamb = z - (katje+1)*1;

  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);

  if(minim == "adaptiv" || minim == "adap_exact")
    {
    fullcond_alle[z]->update_stepwise(modell_alt[ind_mod]);
    fullcond_alle[z]->posteriormode();
    fullcond_alle[z]->wiederholen(fullcond_alle[z],false);
    kriterium = criterion_min(df);
    }

  df = df - fullcond_alle[z]->compute_df();
  fullcond_alle[ind_fullc]->safe_const();
  for(i=0;i<lambdavec[ind_lamb].size();i++)
    {
    if(lambdavec[ind_lamb][i]!=modell_alt[ind_mod])
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[ind_lamb][i]!=-1 && lambdavec[ind_lamb][i]!=0)
        {
        if(possible == "alles" || possible == "spline" || possible == "spfix")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[ind_lamb][i]);
          fullcond_alle[z]->posteriormode();
          fullcond_alle[z]->wiederholen(fullcond_alle[z],false);
          kriterium_versuch = criterion_min(df + fullcond_alle[z]->compute_df());
          fullcond_alle[ind_fullc]->set_const_old();
          }
        }
      else if(lambdavec[ind_lamb][i]==-1)
        {
        if(possible == "alles" || possible == "spfix")
          {
          fullcond_alle[z]->set_inthemodel(-1);
          fullcond_alle[z]->reset_effect(0);
          fullcond_alle[ind_fullc]->posteriormode_single(names_nonp[ind_name],
                                  fullcond_alle[z]->get_data_forfixedeffects(),true);
          fullcond_alle[z]->wiederholen(fullcond_alle[z],false);
          kriterium_versuch = criterion_min(df + 1);
          fullcond_alle[ind_fullc]->set_const_old();
          reset_fix(names_nonp[ind_name][0]);
          }
        }
      else
        {
        if(possible == "alles")
          {
          fullcond_alle[z]->set_inthemodel(0);
          fullcond_alle[z]->reset_effect(0);
          fullcond_alle[ind_fullc]->posteriormode_const();
          fullcond_alle[z]->wiederholen(fullcond_alle[z],false);
          kriterium_versuch = criterion_min(df);
          fullcond_alle[ind_fullc]->set_const_old();
          }
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[ind_lamb][i],6).helpfill(8)
                               + "   " + ST::doubletostring(krit_fkt[i],12) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(minim == "approx_control")
    {
    vector<double> kriterium_control;
    genoptions_mult[0]->out("\n\n");
    fullcond_alle[z]->set_inthemodel(1);
    fullcond_alle[z]->update_stepwise(modell_alt[ind_mod]);
    unsigned y;
    for(y=1;y<fullcond_alle.size();y++)
      {
      if(fullcond_alle[y]->is_identifiable() == false)
        fullcond_alle[y]->set_center(true);
      }
    fullcond_alle[z]->posteriormode();
    fullcond_alle[z]->wiederholen(fullcond_alle[z],true);     // für Haupteffekt
    fullcond_alle[z]->remove_centering_fix();                 // für Interaktion mit fixem Haupteffekt
    fullcond_alle[ind_fullc]->update_linold();

    minexact_nonp_nonp(z,kriterium_control,kriterium,df);
    genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + " Testvalues: approx    exact \n");
    for(i=0;i<kriterium_control.size();i++)
      genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],12) + "   "
                        + ST::doubletostring(kriterium_control[i],6) + "\n");
    }
  else
    {
    fullcond_alle[z]->set_inthemodel(1);
    fullcond_alle[z]->update_stepwise(modell_alt[ind_mod]);
    if(fullcond_alle[z]->is_identifiable() == false)
      fullcond_alle[z]->set_center(true);
    fullcond_alle[z]->posteriormode();
    fullcond_alle[z]->wiederholen(fullcond_alle[z],true);          // für Haupteffekt
    fullcond_alle[z]->remove_centering_fix();                 // für Interaktion mit fixem Haupteffekt
    fullcond_alle[ind_fullc]->update_linold();
    }
  }


void STEPMULTIrun::stepmin_nonp_fix(unsigned & z, vector<double> & krit_fkt, double & kriterium)
  {
  unsigned i;
  double df = 0;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();

  unsigned ind_name = z - katje*anz_fullcond - 1;
  unsigned ind_fullc = katje*anz_fullcond;
  unsigned ind_lamb = z - (katje+1)*1;

  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);

  fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],-1,true);
  if(minim == "adaptiv" || minim == "adap_exact")
    {
    //reset_fix(names_nonp[ind_name][0]);
    fullcond_alle[ind_fullc]->posteriormode_single(names_nonp[ind_name],
                                  fullcond_alle[z]->get_data_forfixedeffects(),false);
    fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],1,false);
    kriterium = criterion_min(df);
    }
  else
    fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],1,false);

  df = df - 1;
  fullcond_alle[ind_fullc]->safe_const();
  reset_fix(names_nonp[ind_name][0]);
  fullcondp.push_back(fullcond_alle[z]);
  fullcond_alle[z]->set_inthemodel(1);

  for(i=0;i<lambdavec[ind_lamb].size();i++)
    {
    if(lambdavec[ind_lamb][i]!=-1)
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[ind_lamb][i]!=0)
        {
        if(possible == "alles" || possible == "spfix")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[ind_lamb][i]);
          fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],-1,false);
          fullcond_alle[z]->posteriormode();
          fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],1,false);
          kriterium_versuch = criterion_min(df + fullcond_alle[z]->compute_df());
          fullcond_alle[ind_fullc]->set_const_old();
          }
        }
      else
        {
        if(possible == "alles" || possible == "rfix")
          {
          fullcond_alle[z]->set_inthemodel(0);
          fullcond_alle[z]->reset_effect(0);
          fullcond_alle[ind_fullc]->posteriormode_const();
          kriterium_versuch = criterion_min(df);
          fullcond_alle[ind_fullc]->set_const_old();
          }
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[ind_lamb][i],6).helpfill(8)
                              + "   " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(minim == "approx_control")
    {
    vector<double> kriterium_control;
    fullcond_alle[z]->set_inthemodel(-1);
    fullcond_alle[z]->reset_effect(0);
    fullcondp.erase(fullcondp.end()-1,fullcondp.end());
    fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],-1,false);
    fullcond_alle[ind_fullc]->posteriormode_single(names_nonp[ind_name],
                                fullcond_alle[z]->get_data_forfixedeffects(),true);
    fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],1,true);
    fullcond_alle[ind_fullc]->update_linold();
    reset_fix(names_nonp[ind_name][0]);

    unsigned y;
    for(y=1;y<fullcond_alle.size();y++)
      {
      if(fullcond_alle[y]->is_identifiable() == false)
        fullcond_alle[y]->set_center(true);
      }
    fullcondp.push_back(fullcond_alle[z]);

    minexact_nonp_fix(z,kriterium_control,kriterium,df);
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + " Testvalues: approx    exact \n");
    for(i=0;i<kriterium_control.size();i++)
      genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],6) + "   "
             + ST::doubletostring(kriterium_control[i],6) + "\n");
    }
  else
    {
    fullcond_alle[z]->set_inthemodel(-1);
    fullcond_alle[z]->reset_effect(0);
    fullcondp.erase(fullcondp.end()-1,fullcondp.end());
    fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],-1,false);
    fullcond_alle[ind_fullc]->posteriormode_single(names_nonp[ind_name],
                                fullcond_alle[z]->get_data_forfixedeffects(),true);
    fullcond_alle[z]->wiederholen_fix(fullcond_alle[z],1,true);
    fullcond_alle[ind_fullc]->update_linold();
    }
  }


void STEPMULTIrun::stepmin_nonp_leer(unsigned & z, vector<double> & krit_fkt, double & kriterium)
  {
  unsigned i;

  double df = 0;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();

  unsigned ind_name = z - katje*anz_fullcond - 1;
  unsigned ind_fullc = katje*anz_fullcond;
  unsigned ind_lamb = z - (katje+1)*1;

  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);

  // if(minim == "adaptiv" || minim == "adap_exact")  hier nicht nötig (siehe koordmin_leer_fix)

  fullcond_alle[ind_fullc]->safe_const();
  fullcondp.push_back(fullcond_alle[z]);
  fullcond_alle[z]->set_inthemodel(1);

  for(i=0;i<lambdavec[ind_lamb].size();i++)
    {
    if(lambdavec[ind_lamb][i]!=0)
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[ind_lamb][i]!=-1)
        {
        if(possible == "alles")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[ind_lamb][i]);
          fullcond_alle[z]->posteriormode();
          kriterium_versuch = criterion_min(df + fullcond_alle[z]->compute_df());
          fullcond_alle[ind_fullc]->set_const_old();
          }
        }
      else
        {
        if(possible == "rfix" || possible == "alles")
          {
          fullcond_alle[z]->set_inthemodel(-1);
          fullcond_alle[z]->reset_effect(0);
          fullcond_alle[ind_lamb]->posteriormode_single(names_nonp[ind_name],
                                   fullcond_alle[z]->get_data_forfixedeffects(),true);
          kriterium_versuch = criterion_min(df + 1);
          reset_fix(names_nonp[ind_name][0]);
          fullcond_alle[ind_fullc]->set_const_old();
          }
        }
      krit_fkt.push_back(kriterium_versuch); // Länge des Vektors muß zu Anzahl Möglichkeiten passen!
      }
    else
      krit_fkt.push_back(kriterium);
    }

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[ind_lamb][i],6).helpfill(8)
                              + "   " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(minim == "approx_control")
    {
    vector<double> kriterium_control;
    fullcond_alle[z]->set_inthemodel(0);

    unsigned y;
    for(y=1;y<fullcond_alle.size();y++)
      {
      if(fullcond_alle[y]->is_identifiable() == false)
        fullcond_alle[y]->set_center(true);
      }

    minexact_nonp_leer(z,kriterium_control,kriterium,df);
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + " Testvalues: approx    exact \n");
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

void STEPMULTIrun::minexact_nonp_nonp(unsigned & z, vector<double> & krit_fkt,
                  double & kriterium, double & df)
  {
  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);

  unsigned ind_name = z - katje*anz_fullcond - 1;
  unsigned ind_fullc = katje*anz_fullcond;
  unsigned ind_mod = z + (katje+1)*(names_fixed.size()-2);
  unsigned ind_lamb = z - (katje+1)*1;

  unsigned i;
  for(i=0;i<lambdavec[ind_lamb].size();i++)
    {
    if(lambdavec[ind_lamb][i]!=modell_alt[ind_mod])
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[ind_lamb][i]!=-1 && lambdavec[ind_lamb][i]!=0)
        {
        if(possible == "alles" || possible == "spline" || possible == "spfix")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[ind_lamb][i]);
          fullcond_alle[ind_fullc]->posteriormode_const();
          posteriormode(posttitle,true);
          kriterium_versuch = criterion_min(df + fullcond_alle[z]->compute_df());
          }
        }
      else if(lambdavec[ind_lamb][i]==-1)
        {
        if(possible == "alles" || possible == "spfix")
          {
          fullcond_alle[z]->set_inthemodel(-1);
          vector<FULLCOND*> fullcond_start = fullcondp;
          vector<double> modell1 = modell_alt;
          modell1[ind_mod] = -1;
          fullcond_einzeln(modell1,modell_alt,z);  // hier muß der Fullcond-Vekor angepaßt werden!!!
          fullcond_alle[ind_fullc]->posteriormode_const();
          posteriormode(posttitle,true);
          kriterium_versuch = criterion_min(df + 1);
          fullcondp = fullcond_start;
          end[0] = fullcondp.size()-1;
          reset_fix(names_nonp[ind_name][0]);
          }
        }
      else
        {
        if(possible == "alles")
          {
          fullcond_alle[z]->set_inthemodel(0);
          vector<FULLCOND*> fullcond_start = fullcondp;
          vector<double> modell1 = modell_alt;
          modell1[ind_mod] = 0;
          fullcond_einzeln(modell1,modell_alt,z);  // hier muß der Fullcond-Vekor angepaßt werden!!!
          fullcond_alle[ind_fullc]->posteriormode_const();
          posteriormode(posttitle,true);
          kriterium_versuch = criterion_min(df);
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
  fullcond_alle[z]->update_stepwise(modell_alt[ind_mod]);
  fullcond_alle[ind_fullc]->posteriormode_const();
  posteriormode(posttitle,true);

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (exact): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[ind_lamb][i],6).helpfill(8) + "   " + 
ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }
  }


void STEPMULTIrun::minexact_nonp_fix(unsigned & z, vector<double> & krit_fkt,
          double & kriterium, double & df)
  {
  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);

  unsigned ind_name = z - katje*anz_fullcond - 1;
  unsigned ind_fullc = katje*anz_fullcond;
  unsigned ind_lamb = z - (katje+1)*1;

  unsigned i;
  end[0] = fullcondp.size()-1;
  fullcond_alle[z]->set_inthemodel(1);
  for(i=0;i<lambdavec[ind_lamb].size();i++)
    {
    if(lambdavec[ind_lamb][i]!=-1)
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[ind_lamb][i]!=0)
        {
        if(possible == "alles" || possible == "spfix")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[ind_lamb][i]);
          fullcond_alle[ind_fullc]->posteriormode_const();
          posteriormode(posttitle,true);
          kriterium_versuch = criterion_min(df + fullcond_alle[z]->compute_df());
          }
        }
      else
        {
        if(possible == "alles" || possible == "rfix")
          {
          fullcond_alle[z]->set_inthemodel(0);
          vector<FULLCOND*> fullcond_start = fullcondp;
          fullcondp.erase(fullcondp.end()-1,fullcondp.end());
          end[0] = fullcondp.size()-1;
          fullcond_alle[z]->reset_effect(0);
          fullcond_alle[ind_fullc]->posteriormode_const();
          posteriormode(posttitle,true);
          kriterium_versuch = criterion_min(df);
          fullcondp = fullcond_start;
          end[0] = fullcondp.size()-1;
          }
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }
  fullcond_alle[z]->set_inthemodel(-1);
  fullcond_alle[z]->reset_effect(0);
  fullcond_alle[ind_fullc]->include_effect(names_nonp[ind_name],
                                fullcond_alle[z]->get_data_forfixedeffects());
  fullcondp.erase(fullcondp.end()-1,fullcondp.end());
  end[0] = fullcondp.size()-1;
  fullcond_alle[ind_fullc]->posteriormode_const();
  posteriormode(posttitle,true);

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (exact): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[ind_lamb][i],6).helpfill(8) + "   " + 
ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }
  }


void STEPMULTIrun::minexact_nonp_leer(unsigned & z, vector<double> & krit_fkt,
                  double & kriterium, double & df)
  {
  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);

  unsigned ind_name = z - katje*anz_fullcond - 1;
  unsigned ind_fullc = katje*anz_fullcond;
  unsigned ind_lamb = z - (katje+1)*1;

  unsigned i;
  end[0] = fullcondp.size()-1;
  fullcond_alle[z]->set_inthemodel(1);
  for(i=0;i<lambdavec[ind_lamb].size();i++)
    {
    if(lambdavec[ind_lamb][i]!=0)
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[ind_lamb][i]!=-1)
        {
        if(possible == "alles")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[ind_lamb][i]);
          fullcond_alle[ind_fullc]->posteriormode_const();
          posteriormode(posttitle,true);
          kriterium_versuch = criterion_min(df + fullcond_alle[z]->compute_df());
          }
        }
      else
        {
        if(possible == "rfix" || possible == "alles")
          {
          fullcond_alle[z]->set_inthemodel(-1);
          vector<FULLCOND*> fullcond_start = fullcondp;
          fullcondp.erase(fullcondp.end()-1,fullcondp.end());
          end[0] = fullcondp.size()-1;
          fullcond_alle[z]->reset_effect(0);
          fullcond_alle[ind_fullc]->include_effect(names_nonp[ind_name],
                                fullcond_alle[z]->get_data_forfixedeffects());
          fullcond_alle[ind_fullc]->posteriormode_const();
          posteriormode(posttitle,true);
          kriterium_versuch = criterion_min(df + 1);
          fullcondp = fullcond_start;
          end[0] = fullcondp.size()-1;
          reset_fix(names_nonp[ind_name][0]);
          }
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }
  fullcond_alle[z]->set_inthemodel(0);
  fullcond_alle[z]->reset_effect(0);
  fullcondp.erase(fullcondp.end()-1,fullcondp.end());
  fullcond_alle[ind_fullc]->posteriormode_const();
  end[0] = fullcondp.size()-1;
  posteriormode(posttitle,true);

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (exact): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[ind_lamb][i],6).helpfill(8) + "   " + 
ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }
  }

//-----------------------------------------------------------------------------------

double STEPMULTIrun::criterion_min(const double & df)
  {
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
  //else //if(criterion=="AUC")
  //  kriterium = -1 * likep_mult[0]->compute_auc();

  return kriterium;
  }

double STEPMULTIrun::criterion_min(const double & df, const ST::string & auswahl)
  {
  double kriterium;

  if(auswahl=="GCV")
    kriterium = likep_mult[0]->compute_gcv(df);
  else if(auswahl=="AIC")
    kriterium = likep_mult[0]->compute_aic(df);
  else if(auswahl=="BIC")
    kriterium = likep_mult[0]->compute_bic(df);
  else if(auswahl=="AIC_imp")
    kriterium = likep_mult[0]->compute_improvedaic(df);
  else if(auswahl=="MSEP")
    kriterium = likep_mult[0]->compute_msep();
  //else //if(auswahl=="AUC")
  //  kriterium = -1 * likep_mult[0]->compute_auc();

  return kriterium;
  }

// -----------------------------------------------------------------------------
// ------------------ Funktionen für Koordinatenmethode ------------------------
// -----------------------------------------------------------------------------

bool STEPMULTIrun::koordabstieg(void)
  {
      // Schleife für Minimierung
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

       if(minim == "exact")
         {
         for(katje=0;katje<kategorien;katje++)
           {
           unsigned z = koordexact_fixfactor(kriteriumiteration2,modeliteration,
                              textiteration,kriterium_aktuell);
           koordexact_nonp(kriteriumiteration2,modeliteration,textiteration,z,kriterium_aktuell);
           }
         katje=0;
         }
       else
         {
         for(katje=0;katje<kategorien;katje++)
           {
           koord_minfix(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell);
           unsigned z = koord_minfactor(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell);
           koord_minnonp(kriteriumiteration2,modeliteration,textiteration,z,kriterium_aktuell);

           if((minim == "adaptiv" || minim == "adap_exact")
                     && likep_mult[0]->get_family() != "Gaussian")
             likep_mult[0]->compute_iwls();
           }
         katje=0;

         if((minim == "adaptiv" || minim == "adap_exact")
                                     && modellematrix.size()>=3)
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


void STEPMULTIrun::koord_minfix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell)
  {
  unsigned i;
  unsigned unten = katje*(anz_fullcond + names_fixed.size()-2);
  unsigned oben = katje*(anz_fullcond + names_fixed.size()-2) + names_fixed.size()-1;

  for(i=unten;i<oben;i++)
    {
    unsigned j = i+1;
    if(modell_alt[i]==-1)
      koord_fix_leer(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell,j);
    else if(modell_alt[i]==0)
      koord_leer_fix(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell,j);
    modell_alt = modell_neu;
    }
  }


void STEPMULTIrun::koord_fix_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & i)
  {
  double df = 0;
  unsigned y;
  for(y=0;y<fullcondp.size();y++)
    df = df + fullcondp[y]->compute_df();

  unsigned ind_name = i - katje*(anz_fullcond + names_fixed.size()-2);
  unsigned ind_fullc = katje*anz_fullcond;

  double kriterium_adaptiv = kriterium_aktuell;
  if(minim == "adaptiv" || minim == "adap_exact")
    {
    //reset_fix(names_fixed[ind_name]);
    int c = column_for_fix(names_fixed[ind_name]);
    vector<ST::string> name_help;
    name_help.push_back(names_fixed[ind_name]);
    fullcond_alle[ind_fullc]->posteriormode_single(name_help,datamatrix(D.getCol(c)),false);
    kriterium_aktuell = criterion_min(df);
    }

  df = df - 1;
  modell_neu[i-1] = 0;
  fullcond_alle[ind_fullc]->safe_const();
  reset_fix(names_fixed[ind_name]);
  fullcond_alle[ind_fullc]->posteriormode_const();
  kriterium_neu = criterion_min(df);
  fullcond_alle[ind_fullc]->set_const_old();
  if(minim == "approx_control")
    {
    posteriormode(posttitle,true);
    double kriterium_test = criterion_min(df);
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_fixed[ind_name] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    include_fix(names_fixed[ind_name]);
    posteriormode(posttitle,true);
    reset_fix(names_fixed[ind_name]);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_fixed[ind_name] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_aktuell,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
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
      int c = column_for_fix(names_fixed[ind_name]);
      vector<ST::string> name_help;
      name_help.push_back(names_fixed[ind_name]);
      fullcond_alle[ind_fullc]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
      modell_neu[i-1] = -1;
      if(kriterium_aktuell < kriterium_neu) // verhindert, daß "approx" schlechter wird!
        posteriormode(posttitle,true);
      }
    }

 if(kriterium_aktuell >= kriterium_neu)
   kriterium_aktuell = kriterium_neu;
 else //if(kriterium_neu > kriterium_aktuell)
    {
    int c = column_for_fix(names_fixed[ind_name]);
    vector<ST::string> name_help;
    name_help.push_back(names_fixed[ind_name]);
    fullcond_alle[ind_fullc]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
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

void STEPMULTIrun::koord_leer_fix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & i)
  {
  double df = 0;
  unsigned y;
  for(y=0;y<fullcondp.size();y++)
    df = df + fullcondp[y]->compute_df();

  unsigned ind_name = i - katje*(anz_fullcond + names_fixed.size()-2);
  unsigned ind_fullc = katje*anz_fullcond;

  double kriterium_adaptiv = kriterium_aktuell;
  // if(minim == "adaptiv" || minim == "adap_exact")   hier nicht nötig, weil Intercept immer erneuert wird!

  df = df + 1;
  modell_neu[i-1] = -1;
  fullcond_alle[ind_fullc]->safe_const();
  int c = column_for_fix(names_fixed[ind_name]);
  vector<ST::string> name_help;
  name_help.push_back(names_fixed[ind_name]);
  fullcond_alle[ind_fullc]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
  kriterium_neu = criterion_min(df);
  fullcond_alle[ind_fullc]->set_const_old();
  if(minim == "approx_control")
    {
    posteriormode(posttitle,true);
    double kriterium_test = criterion_min(df);
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_fixed[ind_name] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    reset_fix(names_fixed[ind_name]);
    fullcond_alle[ind_fullc]->posteriormode_const();
    posteriormode(posttitle,true);
    fullcond_alle[ind_fullc]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_fixed[ind_name] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_aktuell,6) + "\n");
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
       reset_fix(names_fixed[ind_name]);
       modell_neu[i-1] = 0;
       if(kriterium_aktuell < kriterium_neu)
         posteriormode(posttitle,true);
       }
     }

  if(kriterium_aktuell >= kriterium_neu)
    kriterium_aktuell = kriterium_neu;
  else //if(kriterium_neu >= kriterium_aktuell)
    {
    reset_fix(names_fixed[ind_name]);
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

unsigned STEPMULTIrun::koord_minfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell)
  {
  unsigned z = 1 + katje*anz_fullcond;
  unsigned ind_mod = z + (katje+1)*(names_fixed.size()-2);

  while(z<((katje+1)*anz_fullcond) && fullcond_alle[z]->get_fctype()==factor)
     {
     if(modell_alt[ind_mod]==-1 && fullcond_alle[z]->get_forced()==false)
       koord_factor_leer(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell,z);
     else if(modell_alt[ind_mod]==0)
       koord_leer_factor(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell,z);
     modell_alt = modell_neu;
     z = z + 1;
     }
  return z;
  }

void STEPMULTIrun::koord_factor_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & z)
  {
  double df = 0;
  unsigned i;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();

  unsigned ind_name = z - katje*anz_fullcond - 1;
  unsigned ind_fullc = katje*anz_fullcond;
  unsigned ind_mod = z + (katje+1)*(names_fixed.size()-2);

  double kriterium_adaptiv = kriterium_aktuell;
  if(minim == "adaptiv" || minim == "adap_exact")
    {
    //for(i=0;i<names_nonp[ind_name].size();i++)    // "posteriormode_single" ändern: ohne include_effect! (dann ist raus- + reinnehmen überflüssig!)    //  reset_fix(names_nonp[ind_name][i]);
    fullcond_alle[ind_fullc]->posteriormode_single(names_nonp[ind_name],
                                        fullcond_alle[z]->get_data_forfixedeffects(),false);
    kriterium_aktuell = criterion_min(df);
    }

  df = df - fullcond_alle[z]->get_data_forfixedeffects().cols();
  modell_neu[ind_mod] = 0;
  fullcond_alle[ind_fullc]->safe_const();
  for(i=0;i<names_nonp[ind_name].size();i++)
    reset_fix(names_nonp[ind_name][i]);
  fullcond_alle[ind_fullc]->posteriormode_const();
  kriterium_neu = criterion_min(df);
  fullcond_alle[ind_fullc]->set_const_old();
  if(minim == "approx_control")
    {
    posteriormode(posttitle,true);
    double kriterium_test = criterion_min(df);
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    fullcond_alle[0]->include_effect(names_nonp[ind_name],fullcond_alle[z]->get_data_forfixedeffects());
    fullcond_alle[ind_fullc]->posteriormode_const();
    posteriormode(posttitle,true);
    for(i=0;i<names_nonp[ind_name].size();i++)
      reset_fix(names_nonp[ind_name][i]);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_aktuell,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
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
       fullcond_alle[ind_fullc]->posteriormode_single(names_nonp[ind_name],
                                        fullcond_alle[z]->get_data_forfixedeffects(),true);
       modell_neu[ind_mod] = -1;
       if(kriterium_aktuell < kriterium_neu)
         posteriormode(posttitle,true);
       }
     }

  if(kriterium_aktuell >= kriterium_neu)
     kriterium_aktuell = kriterium_neu;
  else //if(kriterium_neu >= kriterium_aktuell)
    {
    fullcond_alle[ind_fullc]->posteriormode_single(names_nonp[ind_name],
                                     fullcond_alle[z]->get_data_forfixedeffects(),true);
    modell_neu[ind_mod] = -1;
    }

  if(minim == "adaptiv" || minim == "adap_exact")
    {
    if(fabs((kriterium_adaptiv - kriterium_aktuell)/kriterium_adaptiv) >= pow10(-6))
      fertig = false;
    if(modell_alt[ind_mod] != modell_neu[ind_mod]
                    && (trace == "trace_on" || trace == "trace_minim"))
      {
      ST::string text;
      maketext("  Trial:",modell_neu,kriterium_aktuell,text,true,trace,false);
      }
    kriterium_alt = kriterium_aktuell;
    modell_alt[ind_mod] = modell_neu[ind_mod];
    modeliteration.push_back(modell_alt);
    }
  }

void STEPMULTIrun::koord_leer_factor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & z)
  {
  double df = 0;
  unsigned i;
  for(i=0;i<fullcondp.size();i++)
    df = df + fullcondp[i]->compute_df();

  unsigned ind_name = z - katje*anz_fullcond - 1;
  unsigned ind_fullc = katje*anz_fullcond;
  unsigned ind_mod = z + (katje+1)*(names_fixed.size()-2);

  double kriterium_adaptiv = kriterium_aktuell;
  // if(minim == "adaptiv" || minim == "adap_exact")   hier überflüssig (siehe oben)!

  df = df + fullcond_alle[z]->get_data_forfixedeffects().cols();
  modell_neu[ind_mod] = -1;
  fullcond_alle[ind_fullc]->safe_const();
  fullcond_alle[ind_fullc]->posteriormode_single(names_nonp[ind_name],
                                 fullcond_alle[z]->get_data_forfixedeffects(),true);
  kriterium_neu = criterion_min(df);
  fullcond_alle[ind_fullc]->set_const_old();
  if(minim == "approx_control")
    {
    posteriormode(posttitle,true);
    double kriterium_test = criterion_min(df);
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    for(i=0;i<names_nonp[ind_name].size();i++)
      reset_fix(names_nonp[ind_name][i]);
    fullcond_alle[ind_fullc]->posteriormode_const();
    posteriormode(posttitle,true);
    fullcond_alle[ind_fullc]->posteriormode_single(names_nonp[ind_name],
                                  fullcond_alle[z]->get_data_forfixedeffects(),true);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[ind_name][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_aktuell,6) + "\n");
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
       for(i=0;i<names_nonp[ind_name].size();i++)
         reset_fix(names_nonp[ind_name][i]);
       modell_neu[ind_mod] = 0;
       if(kriterium_aktuell < kriterium_neu)
         posteriormode(posttitle,true);
       }
     }

  if(kriterium_aktuell >= kriterium_neu)
    kriterium_aktuell = kriterium_neu;
  else if(kriterium_neu >= kriterium_aktuell)
    {
    for(i=0;i<names_nonp[ind_name].size();i++)
      reset_fix(names_nonp[ind_name][i]);
    modell_neu[ind_mod] = 0;
    }

  if(minim == "adaptiv" || minim == "adap_exact")
    {
    if(fabs((kriterium_adaptiv - kriterium_aktuell)/kriterium_adaptiv) >= pow10(-6))
      fertig = false;
    if(modell_alt[ind_mod] != modell_neu[ind_mod]
                    && (trace == "trace_on" || trace == "trace_minim"))
      {
      ST::string text;
      maketext("  Trial:",modell_neu,kriterium_aktuell,text,true,trace,false);
      }
    kriterium_alt = kriterium_aktuell;
    modell_alt[ind_mod] = modell_neu[ind_mod];
    modeliteration.push_back(modell_alt);
    }
  }


void STEPMULTIrun::koord_minnonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z, double & kriterium_aktuell)
  {
  unsigned i;
  for(i=z;i<((katje+1)*anz_fullcond);i++)
    {
    unsigned lambda_ind;
    double kriterium_min;
    double kriterium_test = kriterium_aktuell;

    unsigned ind_mod = i + (katje+1)*(names_fixed.size()-2);
    unsigned ind_fullc = katje*anz_fullcond;
    unsigned ind_lamb = i - (katje+1)*1;
    unsigned ind_name = i - katje*anz_fullcond - 1;

    unsigned y;
    for(y=1;y<fullcondp.size();y++)
      {
      if(fullcondp[y] == fullcond_alle[i])
        fullcondp[y]->remove_centering();
      }
    for(y=1;y<fullcond_alle.size();y++)
      fullcond_alle[y]->set_center(false);    // sorgt dafür, daß Funktionen nicht zentriert werden!

    vector<double> krit_fkt;
    if(modell_alt[ind_mod]==0)
      stepmin_nonp_leer(i,krit_fkt,kriterium_aktuell);
    else if(modell_alt[ind_mod]==-1)
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

    for(y=1;y<fullcond_alle.size();y++)
     {
     if(fullcond_alle[y]->is_identifiable() == false)
       fullcond_alle[y]->set_center(true);    // sorgt dafür, daß Funktionen zentriert werden!
     }

    modell_neu[ind_mod] = lambdavec[ind_lamb][lambda_ind];
    if(minim != "adaptiv" && minim != "adap_exact")
      {
      if(modell_neu[ind_mod] != modell_alt[ind_mod])
        {
        bool neu = modelcomparison(modell_neu,modellematrix);
        fullcond_einzeln(modell_neu,modell_alt,i);
        if(neu==false)
          {
          fullcond_alle[ind_fullc]->posteriormode_const();
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
      if(modell_alt[ind_mod] != modell_neu[ind_mod])
        {
        fullcond_einzeln(modell_neu,modell_alt,i);
        if(modell_neu[ind_mod] == 0)      // noch mal überprüfen!!!
          {
          //fullcond_alle[i]->reset_effect(0);   // Nicht nötig, wegen "fullcond_einzeln"-> fullcond_alle[i] nicht in fullcondp!!!
          if(modell_alt[ind_mod] > 0)
            fullcond_alle[i]->wiederholen_fix(fullcond_alle[i],1,true);
          fullcond_alle[ind_fullc]->posteriormode_const();
          }
        else if(modell_neu[ind_mod] == -1)
          {
          if(modell_alt[ind_mod] == 0)
            fullcond_alle[i]->wiederholen_fix(fullcond_alle[i],-1,true);

          fullcond_alle[i]->reset_effect(0);
          //reset_fix(names_nonp[ind_name][0]);
          fullcond_alle[ind_fullc]->posteriormode_single(names_nonp[ind_name],
                                 fullcond_alle[i]->get_data_forfixedeffects(),false);

          //if(modell_alt[ind_mod] > 0)     // Warum???
          fullcond_alle[i]->wiederholen_fix(fullcond_alle[i],1,true);
          }
        else
          {
          if(modell_alt[ind_mod] > 0)
            fullcond_alle[i]->remove_centering();
          else
            fullcond_alle[i]->wiederholen_fix(fullcond_alle[i],-1,true);

          fullcond_alle[i]->update_stepwise(modell_neu[ind_mod]);
          fullcond_alle[i]->posteriormode();
          fullcond_alle[ind_fullc]->update_linold();
          fullcond_alle[i]->wiederholen(fullcond_alle[i],true);
          fullcond_alle[i]->remove_centering_fix();
          }
        if(trace == "trace_on" || trace == "trace_minim")
          {
          ST::string text;
          maketext("  Trial:",modell_neu,kriterium_min,text,true,trace,false);
          }
        }

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

unsigned STEPMULTIrun::koordexact_fixfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell)
  {
  unsigned i;
  unsigned unten = katje*(anz_fullcond + names_fixed.size()-2);
  unsigned oben = katje*(anz_fullcond + names_fixed.size()-2) + names_fixed.size()-1;

  unsigned ind_name;
  unsigned ind_fullc = katje*anz_fullcond;

  for(i=unten;i<oben;i++)
     {
     ind_name = i + 1 - katje*(anz_fullcond + names_fixed.size()-2);

     if(modell_alt[i]==-1)
       modell_neu[i]= 0;
     else if(modell_alt[i]==0)
       modell_neu[i] = -1;
     if(modelcomparison(modell_neu,modellematrix)==false)
       {
       if(modell_neu[i]==0)
         reset_fix(names_fixed[ind_name]);
       else
         include_fix(names_fixed[ind_name]);
       fullcond_alle[ind_fullc]->posteriormode_const();
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       if(kriteriumiteration2[kriteriumiteration2.size()-1] > kriterium_aktuell)
          {
          if(modell_neu[i]==0)
            include_fix(names_fixed[ind_name]);
          else
            reset_fix(names_fixed[ind_name]);
          modell_neu = modell_alt;
          }
       else
          {
          fertig = false;
          modell_alt = modell_neu;
          kriterium_aktuell = kriteriumiteration2[kriteriumiteration2.size()-1];
          }
       }
     else
       modell_neu = modell_alt;
     }


  unsigned z = 1 + katje*anz_fullcond;
  unsigned ind_mod;

  while(z<((katje+1)*anz_fullcond) && fullcond_alle[z]->get_fctype()==factor)
     {
     ind_name = z - katje*anz_fullcond - 1;
     ind_mod = z + (katje+1)*(names_fixed.size()-2);

     if(modell_alt[ind_mod]==-1 && fullcond_alle[z]->get_forced()==false)
       modell_neu[ind_mod]= 0;
     else if(modell_alt[ind_mod]==0)
       modell_neu[ind_mod] = -1;
     if(modelcomparison(modell_neu,modellematrix)==false)
        {
        if(modell_neu[ind_mod]==0)
          {
          for(i=0;i<names_nonp[ind_name].size();i++)
            reset_fix(names_nonp[ind_name][i]);
          }
        else
          fullcond_alle[ind_fullc]->include_effect(names_nonp[ind_name],fullcond_alle[z]->get_data_forfixedeffects());
        fullcond_alle[ind_fullc]->posteriormode_const();
        newmodel(kriteriumiteration2,modeliteration,textiteration);
        if(kriteriumiteration2[kriteriumiteration2.size()-1] > kriterium_aktuell)
           {
           if(modell_neu[ind_mod]==0)
             fullcond_alle[ind_fullc]->include_effect(names_nonp[ind_name],fullcond_alle[z]->get_data_forfixedeffects());
           else
             {
             for(i=0;i<names_nonp[ind_name].size();i++)
               reset_fix(names_nonp[ind_name][i]);
             }
           modell_neu = modell_alt;
           }
        else
           {
           fertig = false;
           modell_alt = modell_neu;
           kriterium_aktuell = kriteriumiteration2[kriteriumiteration2.size()-1];
           }
        }
     else
        modell_neu = modell_alt;
     z = z + 1;
     }
  return z;
  }


void STEPMULTIrun::koordexact_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z, double & kriterium_aktuell)
  {
  unsigned i;
  unsigned j;
  for(i=z;i<((katje+1)*anz_fullcond);i++)
    {

    unsigned ind_name = i - katje*anz_fullcond - 1;
    unsigned ind_fullc = katje*anz_fullcond;
    unsigned ind_mod = i + (katje+1)*(names_fixed.size()-2);
    unsigned ind_lamb = i - (katje+1)*1;

    modell_neu = modell_alt;
    unsigned lambda_ind;
    vector<double> krit_fkt;
    double df = 0;
    if(modell_alt[ind_mod]==0)
      {
      for(j=0;j<fullcondp.size();j++)
        df = df + fullcondp[j]->compute_df();
      fullcondp.push_back(fullcond_alle[i]);
      minexact_nonp_leer(i,krit_fkt,kriterium_aktuell,df);
      }
    else if(modell_alt[ind_mod]==-1)
      {
      for(j=0;j<fullcondp.size();j++)
        df = df + fullcondp[j]->compute_df();
      df = df - 1;
      reset_fix(names_nonp[ind_name][0]);
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

    modell_neu[ind_mod] = lambdavec[ind_lamb][lambda_ind];
    if(modell_neu[ind_mod] != modell_alt[ind_mod])
      {
      if(modelcomparison(modell_neu,modellematrix)==false)
        {
        fullcond_einzeln(modell_neu,modell_alt,i);
        fullcond_alle[ind_fullc]->posteriormode_const();
        newmodel(kriteriumiteration2,modeliteration,textiteration);
        kriterium_aktuell = kriteriumiteration2[kriteriumiteration2.size()-1];
        }
      }
    modell_alt = modell_neu;
    }
  }



// -----------------------------------------------------------------------------
// ------- Funktionen für die Erstellung des Startmodels -----------------------
// -----------------------------------------------------------------------------

bool STEPMULTIrun::vcm_doppelt(void)
  {
  unsigned i = 1;
  unsigned j;
  unsigned k;
  bool fehler_randomslope = false;
  while(i<anz_fullcond && fehler_randomslope == false)
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


void STEPMULTIrun::initialise_lambdas(vector<vector<ST::string> > & namen_nonp,
       vector<ST::string> & namen_fix, vector<vector<double> > & lambdavector,
       const int & number, const bool & gewichte)
  {
  namen_fix = fullcond_alle[0]->get_datanames();
  unsigned i;

       // berechnet ein Modell, um Gewichte für Augabe usw. zu erhalten!!!
  if(gewichte == true)
    {
    for(katje=0;katje<kategorien;katje++)
      {
      vector<double> modell_init;
      unsigned unten = katje*(anz_fullcond + namen_fix.size()-2);
      unsigned oben = katje*(anz_fullcond + namen_fix.size()-2) + namen_fix.size()-1;

      for(i=unten;i<oben;i++)
         modell_init.push_back(-1);

      unten = (katje+1)*anz_fullcond;
      oben = 1 + katje*anz_fullcond;
      for(i=oben;i<unten;i++)
        {
        if(fullcond_alle[i]->get_fctype() != MCMC::factor)
           fullcond_alle[i]->update_stepwise(100);
        else if(fullcond_alle[i]->get_fctype() == MCMC::factor)
           {
           fullcondp.erase(fullcondp.begin() + katje*anz_fullcond + 1 - katje);   //löscht das Element an Pos.1 aus fullcondp-Vektor           fullcond_alle[katje*anz_fullcond]->include_effect(fullcond_alle[i]->get_datanames(),
                                     fullcond_alle[i]->get_data_forfixedeffects());
           }
        }
      }
    katje=0;

//likep_mult[0]->compute_iwls();
    end[0] = fullcondp.size()-1;
    posteriormode(posttitle,true);
    }

  for(i=1;i<fullcond_alle.size();i++)
     {
     if(i%anz_fullcond != 0)   // 'i modulo anz_fullcond' -> Abfrage, ob fixe Effekte!
       {
       int nummer = number;
       if(fullcond_alle[i]->get_data_forfixedeffects().cols() > 1)  //bei Faktor-Variablen
          {
          if(i < anz_fullcond)
            namen_nonp.push_back(fullcond_alle[i]->get_datanames());
          }
       else
          {
          if(i < anz_fullcond)
            {
            vector<ST::string> names_help;
            names_help.push_back(fullcond_alle[i]->get_datanames()[fullcond_alle[i]->get_datanames().size()-1]);
            namen_nonp.push_back(names_help);
            }
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
  }


unsigned STEPMULTIrun::search_lambdaindex(const double & m,
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


unsigned STEPMULTIrun::search_lambdastartindex(const double & start,
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


void STEPMULTIrun::startwerte(const ST::string & startmodel,
          vector<vector<unsigned> > & startindex,
          vector<vector<double> > & startfix)
  {
  unsigned i;

  //for(katje=0;katje<kategorien;katje++)         // Für alle Kategorien dasselbe Startmodell!
  //  {
    if(startmodel == "empty" || startmodel == "both" || startmodel == "emplin")
      {
      vector<unsigned> indexhelp;
      vector<double> fixhelp;

      for(i=1;i<names_fixed.size();i++)
         fixhelp.push_back(0);

      for(i=1;i<anz_fullcond;i++)
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
      for(i=1;i<anz_fullcond;i++)
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
      for(i=1;i<anz_fullcond;i++)
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
      for(i=1;i<anz_fullcond;i++)
         {
         double start = fullcond_alle[i]->get_lambdastart();
         unsigned index = search_lambdastartindex(start,lambdavec[i-1]);
         indexhelp.push_back(index);
         }

      startindex.push_back(indexhelp);
      startfix.push_back(fixhelp);
      }
    //}
  }


// -----------------------------------------------------------------------------
// ------- Funktionen für die Berechnung neuer Modelle -------------------------
// -----------------------------------------------------------------------------

double STEPMULTIrun::compute_criterion(void)
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


void STEPMULTIrun::newmodel(vector<double> & krit,
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


void STEPMULTIrun::newmodel_fix(const double & mo, vector<double> & krit,
   vector<vector<double> > & mi, vector<ST::string> & textit,
   const ST::string & name)
  {
  if(mo==0)
    reset_fix(name);
  else
    include_fix(name);
  fullcond_alle[katje*anz_fullcond]->posteriormode_const();
  newmodel(krit,mi,textit);
  if(mo==0)
    include_fix(name);
  else
    reset_fix(name);
  }


void STEPMULTIrun::newmodel_factor(const double & mo, const unsigned & index,
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
    fullcond_alle[katje*anz_fullcond]->include_effect(name,fullcond_alle[index]->get_data_forfixedeffects());
  fullcond_alle[katje*anz_fullcond]->posteriormode_const();
  newmodel(krit,mi,textit);
  if(mo==0)
    fullcond_alle[katje*anz_fullcond]->include_effect(name,fullcond_alle[index]->get_data_forfixedeffects());
  else
    {
    for(i=0;i<name.size();i++)
      reset_fix(name[i]);
    }
  }


void STEPMULTIrun::newmodel_nonp(const unsigned & index,
    vector<double> & krit, vector<vector<double> > & mi,
    vector<ST::string> & textit)
  {
  fullcond_einzeln(modell_neu,modell_alt,index);
  fullcond_alle[katje*anz_fullcond]->posteriormode_const();
  newmodel(krit,mi,textit);
  fullcond_einzeln(modell_alt,modell_neu,index);
  }


bool STEPMULTIrun::modelcomparison(const vector<double> & m,
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

void STEPMULTIrun::fullcond_einzeln(const vector<double> & modell1,
         const vector<double> & modell2, const unsigned & index)
  {
  if(smoothing == "global")    // Versuch!!!
    {
    unsigned kat_alt = katje;
    vector<FULLCOND*> fullcond_neu;
    unsigned i,k;
    for(k=0;k<kategorien;k++)
      {
      katje = k;
      fullcond_neu.push_back(fullcond_alle[k*anz_fullcond]);

      unsigned z = 1 + k*anz_fullcond;
      while(z<((k+1)*anz_fullcond) && fullcond_alle[z]->get_fctype()==factor)
        {
        z = z + 1;
        }

      unsigned ind_mod, ind_name;

      for(i=z;i<((k+1)*anz_fullcond);i++)
         {
         ind_mod = i + (k+1)*(names_fixed.size()-2);
         ind_name = i - k*anz_fullcond - 1;

         fullcond_alle[i]->set_inthemodel(modell1[ind_mod]);
         if(modell2[ind_mod]==-1 && index==i)
           reset_fix(names_nonp[ind_name][0]);
         if(modell1[ind_mod]>0)
           {
           fullcond_neu.push_back(fullcond_alle[i]);
           if(i == index)
             fullcond_alle[i]->update_stepwise(modell1[ind_mod]);
           }
         else if(modell1[ind_mod]==0)    // && i == index)
           fullcond_alle[i]->reset_effect(0);
         else if(modell1[ind_mod] == -1) // && i == index)
           {
           fullcond_alle[i]->reset_effect(0);
           if(i == index)
             fullcond_alle[k*anz_fullcond]->include_effect(names_nonp[ind_name],
                                  fullcond_alle[i]->get_data_forfixedeffects());
           }
         }
      }

    fullcondp = fullcond_neu;
    end[0] = fullcondp.size()-1;
    katje = kat_alt;
    }
  /*else //if(smoothing == "local")        // noch für eine Kategorei!
    fullcond_alle[index]->update_stepwise(modell1[names_fixed.size()-2+index]);*/

  }


void STEPMULTIrun::fullcond_komplett(const vector<double> & m)
  {

  if(smoothing == "global")
    {
    unsigned kat_alt = katje;
    vector<FULLCOND*> fullcond_neu;
    unsigned i,k;
    for(k=0;k<kategorien;k++)
      {
      katje = k;
      fullcond_neu.push_back(fullcond_alle[k*anz_fullcond]);

      unsigned ind_mod, ind_name;
      for(i=(1+k*anz_fullcond);i<((k+1)*anz_fullcond);i++)
        {
        ind_mod = i + (k+1)*(names_fixed.size()-2);
        ind_name = i - k*anz_fullcond - 1;

        fullcond_alle[i]->set_inthemodel(m[ind_mod]);
        if(m[ind_mod]>0)
          {
          fullcond_alle[i]->update_stepwise(m[ind_mod]);
          fullcond_neu.push_back(fullcond_alle[i]);
          }
        else if(m[ind_mod]==0)
          fullcond_alle[i]->reset_effect(0);
        else if(m[ind_mod] == -1)
          {
          fullcond_alle[i]->reset_effect(0);
          fullcond_alle[anz_fullcond*k]->include_effect(names_nonp[ind_name],
                                  fullcond_alle[i]->get_data_forfixedeffects());
          }
        }
      fullcond_alle[anz_fullcond*k]->posteriormode_const();
      }

    fullcondp = fullcond_neu;
    end[0] = fullcondp.size()-1;
    katje = kat_alt;
    }
  /*else //if(smoothing == "local")     // noch für eine Kategorie!
    {
    unsigned i;
    for(i=1;i<fullcond_alle.size();i++)
      {
      fullcond_alle[i]->set_lambda_nr();
      fullcond_alle[i]->update_stepwise(m[names_fixed.size()-2+i]);
      }
    }*/  

  }


void STEPMULTIrun::fix_komplett(const vector<double> &  modell)
  {

  if(smoothing == "global")
    {
    unsigned kat_alt = katje;
    unsigned z,k,unten,oben, ind_name;
    for(k=0;k<kategorien;k++)
      {
      katje = k;
      unten = k * (anz_fullcond + names_fixed.size()-2);
      oben = k * (anz_fullcond + names_fixed.size()-2) + names_fixed.size()-1;
      for(z=unten;z<oben;z++)
        {
        ind_name = z - k*(anz_fullcond + names_fixed.size()-2) + 1;

        if(modell[z]==0)
          reset_fix(names_fixed[ind_name]);
        else if(modell[z]==-1)
          {
          unsigned i = 1;
          bool rein = false;
          while(i<fullcond_alle[k*anz_fullcond]->get_datanames().size() && rein==false)
            {
            if(fullcond_alle[k*anz_fullcond]->get_datanames()[i]==names_fixed[ind_name])
              rein = true;
            i = i + 1;
            }
          if(rein==false)
            include_fix(names_fixed[ind_name]);
          }
        }

      unten = (k+1) * (anz_fullcond + names_fixed.size()-2);
      for(z=oben;z<unten;z++)
        {
        ind_name = z - k*(anz_fullcond + names_fixed.size()-2) - (names_fixed.size()-1);

        bool gefunden = false;
        unsigned i = 1;
        while(i<fullcond_alle[k*anz_fullcond]->get_datanames().size() && gefunden==false)
          {
          if(fullcond_alle[k*anz_fullcond]->get_datanames()[i]==names_nonp[ind_name][0])    // ersetzen durch reset_fix?
            {
            gefunden = true;
            fullcond_alle[k*anz_fullcond]->reset_effect(i);
            }
            i = i + 1;
          }
        if(gefunden==true && names_nonp[ind_name].size()>1)
          {
          unsigned j;
          for(j=1;j<names_nonp[ind_name].size();j++)
            reset_fix(names_nonp[ind_name][j]);
          }
        }
      }
    katje = kat_alt;    
    }  // END: if(smoothing == "global")
  }


void STEPMULTIrun::reset_fix(const ST::string & name)
  {
  bool raus = false;
  unsigned j = 1;
  while(j<fullcond_alle[katje*anz_fullcond]->get_datanames().size() && raus==false)
     {
     if(fullcond_alle[katje*anz_fullcond]->get_datanames()[j]==name)
        {
        raus = true;
        fullcond_alle[katje*anz_fullcond]->reset_effect(j);
        }
     j = j + 1;
     }
  }


void STEPMULTIrun::include_fix(const ST::string & name)
  {
  int i = column_for_fix(name);
  vector<ST::string> help_name;
  help_name.push_back(name);
  fullcond_alle[katje*anz_fullcond]->include_effect(help_name,datamatrix(D.getCol(i)));
  }

int STEPMULTIrun::column_for_fix(const ST::string & name)
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

bool STEPMULTIrun::make_pause(void)
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


void STEPMULTIrun::maketext(const ST::string & h, const vector<double> & m,
                          const double & a, ST::string & text,
                          const bool & neutext, const ST::string & tr,
                          const bool & datei)
  {
  if(tr == "trace_on" || trace == "trace_minim")
    {
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out(h);
    }
  ST::string modeltext = "  ";
  if(neutext==true)
    {
    modeltext = modeltext + likep_mult[0]->get_responsename() + "_1 = ";
    unsigned i;
    unsigned k = 1;
    modeltext = modeltext + fullcond_alle[0]->get_effect();
    for(i=1;i<fullcondp.size();i++)
      {
      if(fullcondp[i] == fullcond_alle[k*anz_fullcond])
        {
        k++;
        modeltext = modeltext + "\n                  " + likep_mult[0]->get_responsename() + "_"
                    + ST::inttostring(k) + " = " + fullcondp[i]->get_effect();
        }
      else
        modeltext = modeltext + " + " + fullcondp[i]->get_effect();
      }
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


void STEPMULTIrun::options_text(const int & number,
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
     
  for(i=1;i<anz_fullcond;i++)
     {
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  OPTIONS FOR NONPARAMETRIC TERM: "
          + names_nonp[i-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     if(fullcond_alle[i]->get_lambdamin()!=0 && fullcond_alle[i]->get_lambdamin()!=-1)
          {
          genoptions_mult[0]->out("  Minimal value for the smoothing parameter: "
             + ST::doubletostring(fullcond_alle[i]->get_lambdamin()) + "\n");
          fullcond_alle[i]->update_stepwise(fullcond_alle[i]->get_lambdamin());
          if(fullcond_alle[i]->get_df_equidist()==false && fullcond_alle[i]->get_lambdamin_opt()==false)
            genoptions_mult[0]->out("  This is equivalent to degrees of freedom: "
               + ST::doubletostring(fullcond_alle[i]->compute_df(),6) + "\n");
          else
            genoptions_mult[0]->out("  This is equivalent to degrees of freedom: approximately "
               + ST::doubletostring(fullcond_alle[i]->get_df_lambdamin()) + ", exact "
               + ST::doubletostring(fullcond_alle[i]->compute_df(),6) + "\n");
          }
     if(fullcond_alle[i]->get_lambdamax()!=0 && fullcond_alle[i]->get_lambdamax()!=-1)
          {
          genoptions_mult[0]->out("  Maximal value for the smoothing parameter: "
             + ST::doubletostring(fullcond_alle[i]->get_lambdamax(),6) + "\n");
          fullcond_alle[i]->update_stepwise(fullcond_alle[i]->get_lambdamax());
          if(fullcond_alle[i]->get_df_equidist()==false && fullcond_alle[i]->get_lambdamax_opt()==false)
            genoptions_mult[0]->out("  This is equivalent to degrees of freedom: "
               + ST::doubletostring(fullcond_alle[i]->compute_df(),6) + "\n");
          else
            genoptions_mult[0]->out("  This is equivalent to degrees of freedom: approximately "
               + ST::doubletostring(fullcond_alle[i]->get_df_lambdamax()) + ", exact "
               + ST::doubletostring(fullcond_alle[i]->compute_df(),6) + "\n");
          }
     if(fullcond_alle[i]->get_df_equidist()==true)
         genoptions_mult[0]->out("  Number of different smoothing parameters with equidistant degrees of freedom: "
            + ST::doubletostring(fullcond_alle[i]->get_number()) + "\n");
     else if(fullcond_alle[i]->get_fctype()!=MCMC::factor)
         genoptions_mult[0]->out("  Number of different smoothing parameters on a logarithmic scale: "
            + ST::doubletostring(fullcond_alle[i]->get_number()) + "\n");
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
             fullcond_alle[i]->update_stepwise(lambdavec[i-1][startindex[j][i-1]]);
             genoptions_mult[0]->out("  This is equivalent to degrees of freedom: "
                + ST::doubletostring(fullcond_alle[i]->compute_df(),6) + "\n");
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

void STEPMULTIrun::make_graphics(const ST::string & name,
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


void STEPMULTIrun::make_tex_end(ST::string & path, const vector<double> & modell)
  {
  ST::string path_batch = path + "_graphics.prg";
  ST::string path_splus = path +  "_splus.txt";
  //ST::string path_stata = path +  "_stata.do";
  ST::string path_tex = path + "_model_summary.tex";

  outtex << "\n\\noindent {\\bf \\large Final Predictor:}\\\\" << endl;
  make_predictor();

  unsigned j, ind_mod;
  unsigned k = 0;
  outtex << "\n\\noindent {\\bf \\large Final Properties:}\\\\ \n\\\\" << endl;
  for(j=1;j<fullcond_alle.size();j++)
     {
     ind_mod = j + (k+1)*(names_fixed.size()-2);
     if(j%anz_fullcond != 0)
       {
       if(modell[ind_mod]!=0 && modell[ind_mod]!=-1)
         {
         vector<ST::string> prior = fullcond_alle[j]->get_priorassumptions();
         outtex << prior[0] << "\\\\" << endl
                << "smoothing parameter: $\\lambda = "
                << ST::doubletostring(modell[ind_mod],6)
                << " \\,\\, \\hat{=} \\,\\, \\mbox{df} = "
                << ST::doubletostring(fullcond_alle[j]->compute_df(),6)
                << "$ \\\\ \n\\\\" << endl;
         }
       }
     else
       {
       k++;
       outtex << ST::inttostring(k+1) << ". Response Category: \\\\" << endl << endl;
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


void STEPMULTIrun::make_options(void)
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


void STEPMULTIrun::make_predictor(void)
  {
  unsigned i,k;
  unsigned j = 0;
  char hcharu = '_';
  ST::string hstringu = "\\_";
  for(k=0;k<kategorien;k++)
    {
    ST::string term = "$\\eta_" + ST::inttostring(k+1) + "$ & $=$ & $\\gamma_0";

    for(i=1;i<fullcond_alle[k*anz_fullcond]->get_datanames().size();i++)
      term = term + " + " + fullcond_alle[anz_fullcond*k]->get_datanames()[i].insert_string_char(hcharu,hstringu);
    j++;
    if(k<kategorien-1)
      {
      while(j<fullcondp.size() && fullcondp[j] != fullcond_alle[(k+1)*anz_fullcond])
        {
        term = term + " + " + fullcondp[j]->get_term_symbolic();
        j++;
        }
      }
    else // if(k==kategorien-1)
      {
      while(j<fullcondp.size())
        {
        term = term + " + " + fullcondp[j]->get_term_symbolic();
        j++;
        }
      }

    outtex << endl << "\n\\begin{tabular}{ccp{12cm}}\n" << term
           << "$\n\\end{tabular}\n\\\\ \n\\\\" << endl;
    }  outtex << criterion.insert_string_char(hcharu,hstringu) << " = "         << ST::doubletostring(kriterium_tex,6) << " \\\\ \n\\\\" << endl;  }

void STEPMULTIrun::make_model(void)
  {
  //Vert-Fam wird übergeben
  ST::string fam = likep_mult[0]->get_family();
  fam = fam.replaceallsigns('_', ' ');

  //Anz. Beob. wird übergeben
  unsigned obs = likep_mult[0]->get_nrobs();

  //Name der Resp.-Var. übergeben
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


void STEPMULTIrun::make_prior(vector<vector<unsigned> > & startindex)
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
  for(i=1;i<anz_fullcond;i++)
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


void STEPMULTIrun::make_fixed_table(void)
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
  unsigned r;

  unsigned k;
  for(k=0;k<kategorien;k++)
    {
    r = 2;
    if(k==0)
      outtex << "\n\\newpage \n" << endl << "\n\\noindent {\\bf \\large Fixed Effects:}\\\\"
             << endl << "\\\\" << endl;
    else
      outtex << "\n\\newpage \n" << endl << "\n\\noindent {\\bf \\large Fixed Effects (" << ST::inttostring(k+1)
             << ". Response Category):}\\\\"
             << endl << "\\\\" << endl;

    outtex << "\\begin{tabular}{|r|r|}" << endl << "\\hline" << endl
           << "Variable & Mean \\\\" << endl << "\\hline" << endl;
    h = fullcond_alle[k*anz_fullcond]->get_results_latex();
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
  }
  

void STEPMULTIrun::make_plots(ST::string & path_batch,
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

      // Plotstyle: noplot, plotnonp, drawmap, drawmapgraph
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

}
