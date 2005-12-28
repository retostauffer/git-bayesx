//---------------------------------------------------------------------------
#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include<StatwinFrame.h>
#include<statwin_haupt.h>

#endif

#include<bayesreg.h>
#include<bayesreg2.h>
#include<typeinfo.h>
#include<stddef.h>


bool bayesreg::create_geokriging(const unsigned & collinpred)
  {

  ST::string pathnonp;
  ST::string pathres;

  long h;
  double nu,maxdist,p,q,lambda,startlambda;
  double a1,b1;
  unsigned nrknots, maxsteps;
  bool full;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial_geokriging.checkvector(terms,i) == true )
      {

      MCMC::fieldtype type;
      if (terms[i].options[0] == "geokriging")
        type = MCMC::kriging;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      nrknots = unsigned(h);

      f = (terms[i].options[2]).strtodouble(nu);
      if(nu!=0.5 && nu!=1.5&& nu!=2.5 && nu!=3.5)
        {
        outerror("ERROR: Invalid value for nu\n");
        return true;
        }
      f = (terms[i].options[3]).strtodouble(maxdist);
      if(maxdist<=0) // wähle maxdist so, dass Korrelation für Punkte mitmaximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;//4.605170186;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;//6.638352068;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;//8.022007057;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;//9.158140446;
          }
        }

      if(terms[i].options[4] == "true")
        {
        full=true;
        }
      else
        {
        full=false;
        }

      f = (terms[i].options[6]).strtodouble(p);
      f = (terms[i].options[7]).strtodouble(q);
      f = (terms[i].options[8]).strtolong(h);
      maxsteps = unsigned(h);

      f = (terms[i].options[9]).strtodouble(lambda);
      f = (terms[i].options[10]).strtodouble(startlambda);

      f = (terms[i].options[12]).strtodouble(a1);
      f = (terms[i].options[13]).strtodouble(b1);

      if (f==1)
        return true;

      mapobject * mapp;                           // pointer to mapobject
      int objpos = findstatobject(*statobj,terms[i].options[11],"map");
      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          outerror("ERROR: map object " + terms[i].options[11] + " is not existing\n");
        else
          outerror("ERROR: " + terms[i].options[11] + " is not a map object\n");
        return true;
        }
      MAP::map m = mapp->getmap();
      if(!m.centroids_existing())
        {
        outerror("ERROR: map object doesn´t contain centroids\n");
        return true;
        }

      datamatrix knotdata;
      if(terms[i].options[5]!="" && !full)
        {
        dataobject * datap;                           // pointer to datasetobject
        int objpos = findstatobject(*statobj,terms[i].options[5],"dataset");
        if (objpos >= 0)
          {
          statobject * s = statobj->at(objpos);
          datap = dynamic_cast<dataobject*>(s);
          if (datap->obs()==0 || datap->getVarnames().size()==0)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " does not contain any data\n");
            return true;
            }
          else if (datap->getVarnames().size()>2)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " contains more than two variables\n");
            return true;
            }
          }
        else
          {
          outerror("ERROR: dataset object " + terms[i].options[5] + " is not existing\n");
          return true;
          }
        list<ST::string> knotnames = datap->getVarnames();
        ST::string expr = "";
        datap->makematrix(knotnames,knotdata,expr);
        }
      else
        {
        knotdata = datamatrix(1,1,0);
        }

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_geokriging.raw","_geokriging.res","_geokriging");

      if ( check_gaussian() )
        {

        fckriging.push_back(
        FULLCOND_kriging2(&generaloptions[generaloptions.size()-1],
                distr[distr.size()-1],fcconst_intercept,
                D.getCol(j),m,terms[i].options[11],
                knotdata,
                nrknots,nu,maxdist,p,q,maxsteps,full,type,
                title,
                pathnonp,
                pathres,
                lambda,
                startlambda,
                collinpred
                ));

        }
      else if( check_iwls(true) )
        {

        ST::string proposal = terms[i].options[14];

        bool iwlsmode;
        if(proposal == "iwlsmode")
          iwlsmode = true;
        else
          iwlsmode = false;

        unsigned updateW;
        f = (terms[i].options[15]).strtolong(h);
        updateW = unsigned(h);

        bool updatetau;
        if(terms[i].options[16] == "false" || constlambda.getvalue() == true)
          updatetau = false;
        else
          updatetau = true;

        double fstart;
          f = (terms[i].options[17]).strtodouble(fstart);

        fckriging.push_back(
        FULLCOND_kriging2(&generaloptions[generaloptions.size()-1],
                distr[distr.size()-1],fcconst_intercept,
                D.getCol(j),m,terms[i].options[11],
                knotdata,
                nrknots,nu,maxdist,p,q,maxsteps,full,type,
                title,
                pathnonp,
                pathres,
                lambda,
                startlambda,
                iwlsmode,
                updateW,
                updatetau,
                fstart,
                collinpred
                ));

        }

      if (constlambda.getvalue() == true)
        fckriging[fckriging.size()-1].set_lambdaconst(lambda);

      vector<ST::string> na;
      na.push_back(terms[i].varnames[0]);

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                     "_geokriging_var.raw","_geokriging_var.res","_geokriging_variance");

      fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                     &fckriging[fckriging.size()-1],
                                     distr[distr.size()-1],
                                     a1,
                                     b1,
                                     title,pathnonp,pathres,
                                     false,collinpred)
                         );

      fckriging[fckriging.size()-1].init_names(na);
      fckriging[fckriging.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fckriging[fckriging.size()-1]);

     if (constlambda.getvalue() == false)
        {
        if(terms[i].options[18]=="true")
          fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
        }

      }
    }

  return false;
  }


bool bayesreg::create_spatial(const unsigned & collinpred)
  {

  ST::string pathnonpv;
  ST::string pathresv;

  long h;
  double hd;
  unsigned min,max;
  int f;
  double lambda,a1,b1,alpha;
  unsigned i;
  int j1,j2;
  bool iwls;
  bool varcoeff=false;
  unsigned updateW;
  bool updatetau,Laplace;
  double ftune;
  ST::string proposal;
  vector<ST::string> na;

  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial.checkvector(terms,i) == true)
      {

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);  // interacting var
      if (terms[i].type == "varcoeffspatial" || terms[i].type == "tvarcoeffspatial")
        {
        varcoeff=true;
        j2 = terms[i].varnames[1].isinlist(modelvarnamesv);  // effect modifier
        }

      mapobject * mapp;                           // pointer to mapobject

      int objpos = findstatobject(*statobj,terms[i].options[1],"map");

      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          {
          if ((terms[i].options[1] == "") || (terms[i].options[1] == " "))
            outerror("ERROR: map object must be specified to estimate a spatial effect\n");
          else
            outerror("ERROR: map object " + terms[i].options[1] + " is not existing\n");
          }
        else
          outerror("ERROR: " + terms[i].options[1] + " is not a map object\n");
        return true;
        }

      MAP::map m = mapp->getmap();
      bool isconnected = m.isconnected();
      if (isconnected==false)
        {
        outerror("ERROR: map is disconnected, spatial effect cannot be estimated\n");
        return true;
        }


      f = (terms[i].options[2]).strtolong(h);
      min = unsigned(h);

      f = (terms[i].options[3]).strtolong(h);
      max = unsigned(h);

      f = (terms[i].options[4]).strtodouble(hd);
      lambda = hd;

      f = (terms[i].options[5]).strtodouble(a1);

      f = (terms[i].options[6]).strtodouble(b1);

      proposal = terms[i].options[7];

      f = (terms[i].options[8]).strtolong(h);
      updateW = unsigned(h);

      if (terms[i].options[9] == "true")
        updatetau=true;
      else
        updatetau=false;

      f = (terms[i].options[10]).strtodouble(ftune);

      if (terms[i].options[16] == "true")
        Laplace=true;
      else
        Laplace=false;

      f = (terms[i].options[18]).strtodouble(alpha);

      if (f==1)
        return true;

      ST::string titlev;

      if (varcoeff == true)
        {
        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                   terms[i].varnames[0],"_spatial.raw","_spatial.res","_spatial");

        make_paths(collinpred,pathnonpv,pathresv,titlev,terms[i].varnames[1],
                   terms[i].varnames[0],"_spatial_var.raw","_spatial_var.res",
                   "_spatial_variance");
        }
      else
        {
        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_spatial.raw","_spatial.res","_spatial");

        make_paths(collinpred,pathnonpv,pathresv,titlev,terms[i].varnames[0],"",
                 "_spatial_var.raw","_spatial_var.res","_spatial_variance");
        }

      if (proposal != "cp")
        iwls=true;
      else
        iwls=false;

      if ( (check_gaussian()) || (check_iwls(iwls)) )
        {

        if (varcoeff == true)
          {
          fcnonpgaussian.push_back(
          FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
          distr[distr.size()-1],fcconst_intercept,m,terms[i].options[1],
          D.getCol(j2),D.getCol(j1),title,pathnonp,pathres,collinpred,lambda));

          na.push_back(terms[i].varnames[1]);
          na.push_back(terms[i].varnames[0]);
          fcnonpgaussian[fcnonpgaussian.size()-1].init_names(na);
          }
        else
          {
          fcnonpgaussian.push_back(
          FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
          distr[distr.size()-1],D.getCol(j1),fcconst_intercept,m,
          terms[i].options[1],title,pathnonp,pathres,collinpred,lambda));

          fcnonpgaussian[fcnonpgaussian.size()-1].init_name(terms[i].varnames[0]);
          }

        if(Laplace)
          fcnonpgaussian[fcnonpgaussian.size()-1].set_Laplace();

        if (fcnonpgaussian[fcnonpgaussian.size()-1].get_errors().size() > 0)
          {
          unsigned i;
          for(i=0;i<fcnonpgaussian[fcnonpgaussian.size()-1].get_errors().size();i++)
            errormessages.push_back(
            fcnonpgaussian[fcnonpgaussian.size()-1].get_errors()[i]);
          return true;
          }


        if (constlambda.getvalue() == true)
          {
          if (check_nongaussian())
            fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS(updateW,true);
          fcnonpgaussian[fcnonpgaussian.size()-1].set_lambdaconst(lambda);
          }
        else
          {

          if ( (check_nongaussian()) && (proposal == "iwls")
              && (updatetau==false) )
            fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS(updateW);

          if ( (check_nongaussian()) && (proposal == "iwlsmode")
            && (updatetau==false) )
            fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS(updateW,true);

          if ( (check_nongaussian()) && (proposal == "iwls")
            && (updatetau==true) )
            fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS_hyperblock(
                                                                 updateW,a1,b1);

          if ( (check_nongaussian()) && (proposal == "iwlsmode")
            && (updatetau==true) )
            fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS_hyperblock(
                                                      updateW,a1,b1,true);
          }

        if (terms[i].options[17] == "true")
          fcnonpgaussian[fcnonpgaussian.size()-1].set_stationary(alpha);

        fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);

        if (constlambda.getvalue() == false)
          {

          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
          &fcnonpgaussian[fcnonpgaussian.size()-1],distr[distr.size()-1],a1,b1,
          titlev,pathnonpv,pathresv,false,collinpred));

          if(terms[i].options[14]=="true")
            fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

          if ( (check_nongaussian()) && (proposal == "iwls")
            && (updatetau==true) )
            fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();

          if ( (check_nongaussian()) && (proposal == "iwlsmode")
            && (updatetau==true) )
            fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();

          bool alphafix = false;
          if (terms[i].options[19] == "true")
            alphafix = true;
          if (terms[i].options[17] == "true")
            fcvarnonp[fcvarnonp.size()-1].set_stationary(alpha,alphafix);

          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

          if(Laplace)
            fcvarnonp[fcvarnonp.size()-1].set_Laplace();

          }

        if ((terms[i].options[0] == "tspatial") ||
            (terms[i].options[0] == "tvarcoeffspatial")
           )
          {

          if(varcoeff==true)
            make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],terms[i].varnames[0],
                   "_spatial_tvar.raw","_spatial_tvar.res","_spatial_tvariance");
          else
            make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_spatial_tvar.raw","_spatial_tvar.res","_spatial_tvariance");

          unsigned v = nu.getvalue();
          unsigned nrrows;
          f = (terms[i].options[15]).strtolong(h);
          nrrows = unsigned(h);

          fctvariance2dim.push_back(FULLCOND_tvariance2dim(&generaloptions[generaloptions.size()-1],
          &fcnonpgaussian[fcnonpgaussian.size()-1],v,title,
          pathnonp,pathres,nrrows,false));

          fctvariance2dim[fctvariance2dim.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fctvariance2dim[fctvariance2dim.size()-1]);

          if(Laplace)
            fctvariance2dim[fctvariance2dim.size()-1].set_Laplace();

          }

        }
      else
        {

        if (varcoeff == true)
          {
          Pmatrices.push_back(PenaltyMatrix(D.getCol(j2),terms[i].varnames[1],
          m,min,max));
          }
        else
          Pmatrices.push_back(PenaltyMatrix(D.getCol(j1),terms[i].varnames[0],
          m,min,max));

        if (Pmatrices[Pmatrices.size()-1].get_errormessages().size() > 0)
          {
          outerror(Pmatrices[Pmatrices.size()-1].get_errormessages());
          return true;
          }

        fcnonp.push_back( FULLCOND_nonp(
                                       &generaloptions[generaloptions.size()-1],
                                        distr[distr.size()-1],
                                        &Pmatrices[Pmatrices.size()-1],
                                        fcconst_intercept,
                                        lambda,
                                        pathnonp,
                                        pathres,title,terms[i].options[1],
                                        collinpred
                                        )
                           );

        fcnonp[fcnonp.size()-1].init_name(terms[i].varnames[0]);

        fcnonp[fcnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcnonp[fcnonp.size()-1]);

        fcvarnonp.push_back(FULLCOND_variance_nonp(
        &generaloptions[generaloptions.size()-1],
                                         &fcnonp[fcnonp.size()-1],
                                         distr[distr.size()-1],
                                         a1,
                                         b1,
                                         titlev,pathnonpv,pathresv,
                                         false,collinpred)
                             );


        if (constlambda.getvalue() == true)
          fcvarnonp[fcvarnonp.size()-1].set_constlambda();
        if(terms[i].options[14]=="true")
          fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

        }


      }   // end: if ( nonpspatial.checkvector(terms,i) == true )

    } //  end:  for(i=0;i<terms.size();i++)


  return false;
  }


bool bayesreg::create_spatialxy(const unsigned & collinpred)
  {

  ST::string pathmap;
  ST::string mapname;

//  long h;
//  unsigned min,max;
  double a1,b1;
  double lambda;
  double maxdist;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatialxy.checkvector(terms,i) == true )
      {

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);


//      f = (terms[i].options[1]).strtolong(h);
//      min = unsigned(h);

//      f = (terms[i].options[2]).strtolong(h);
//      max = unsigned(h);

      f = (terms[i].options[3]).strtodouble(lambda);

      f = (terms[i].options[4]).strtodouble(a1);

      f = (terms[i].options[5]).strtodouble(b1);

      f = (terms[i].options[6]).strtodouble(maxdist);

      if (f==1)
        return true;

      pathmap = outfile.getvalue()+ add_name + "_" + terms[i].varnames[0] + "_" +
                terms[i].varnames[1] + "_dist" + ST::doubletostring(maxdist)
                + ".bnd";

      mapname = terms[i].varnames[0] + "_" + terms[i].varnames[1] + "_dist"
                + ST::doubletostring(maxdist) + add_name;

      ST::string help = terms[i].varnames[0]+"_"+terms[i].varnames[1];

      make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_spatial.raw","_spatial.res","_spatial");

      if (check_gaussian())
        {

        fcnonpgaussian.push_back( FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
                                        distr[distr.size()-1],
                                        D.getCol(j1),D.getCol(j2),
                                        fcconst_intercept,
                                        lambda,maxdist,mapname,title,
                                        pathnonp, pathres,pathmap,collinpred
                                        )
                           );

        fcnonpgaussian[fcnonpgaussian.size()-1].init_name("regionnr");

        help = terms[i].varnames[0]+"_"+terms[i].varnames[1];

        make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_spatial_var.raw","_spatial_var.res","_spatial_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                       &fcnonpgaussian[fcnonpgaussian.size()-1],
                                       distr[distr.size()-1],a1,b1,title,
                                       pathnonp,pathres,false,collinpred)
                           );

        if (constlambda.getvalue() == true)
          fcvarnonp[fcvarnonp.size()-1].set_constlambda();



        fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

        }
      else
        {


        }


      }   // end: if ( nonpspatial.checkvector(terms,i) == true )

    } //  end:  for(i=0;i<terms.size();i++)


  return false;
  }


bool bayesreg::create_randomslope(const unsigned & collinpred)
  {

  ST::string pathnonp2;
  ST::string pathres2;
  ST::string pathresfixed;
  ST::string title2;

  bool iwlsmode=false;
  bool updatetau;

  unsigned i;
  int j1,j2;
  bool inclf;
  double a1,b1,lambda;
  int f;
  for(i=0;i<terms.size();i++)
    {
    if ( randomeffslope.checkvector(terms,i) == true )
      {

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);

      if (terms[i].options[1] == "true")
        inclf = false;
      else
        inclf = true;

      f = (terms[i].options[2]).strtodouble(lambda);

      f = (terms[i].options[3]).strtodouble(a1);

      f = (terms[i].options[4]).strtodouble(b1);

      if (terms[i].options[5] == "iwlsmode")
        iwlsmode = true;

      if (terms[i].options[6] == "true")
        updatetau = true;
      else
        updatetau = false;

      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_random.raw","_random.res","_random");

      pathresfixed = pathres.substr(0,pathres.length()-4) + "_fixed.res";           

      make_paths(collinpred,pathnonp2,pathres2,title2,terms[i].varnames[1],
                 terms[i].varnames[0],"_random_var.raw",
                 "_random_var.res","_random_variance");




                 

      if (check_gaussian())
        {

        fcrandomgaussian.push_back(
        FULLCOND_random_gaussian(&generaloptions[generaloptions.size()-1],
                                                        distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        D.getCol(j1),
                                                        D.getCol(j2),
                                                        title,
                                                        pathnonp,
                                                        pathres,pathresfixed,
                                                        lambda,
                                                        inclf,collinpred
                                                        )
                          );
        fcrandomgaussian[fcrandomgaussian.size()-1].init_name(terms[i].varnames[1]);
        fcrandomgaussian[fcrandomgaussian.size()-1].set_fcnumber(fullcond.size());

        if (constlambda.getvalue() == true || terms[i].options[8]=="true")
          {
          fcrandomgaussian[fcrandomgaussian.size()-1].set_lambdaconst(lambda);
          fullcond.push_back(&fcrandomgaussian[fcrandomgaussian.size()-1]);

          }
        else
          {
          fullcond.push_back(&fcrandomgaussian[fcrandomgaussian.size()-1]);
          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fcrandomgaussian[fcrandomgaussian.size()-1],
                            distr[distr.size()-1],a1,b1,title2,pathnonp2,
                            pathres2,false,collinpred));
          if(terms[i].options[7]=="true")
                fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
          }


        }

      else
        {

        fcrandom.push_back(
        FULLCOND_random_nongaussian(&generaloptions[generaloptions.size()-1],
                                                        distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        D.getCol(j1),
                                                        D.getCol(j2),
                                                        title,
                                                        pathnonp,
                                                        pathres,pathresfixed,
                                                        lambda,iwlsmode,inclf,
                                                        collinpred
                                                        )
                          );

        fcrandom[fcrandom.size()-1].init_name(terms[i].varnames[1]);
        fcrandom[fcrandom.size()-1].set_fcnumber(fullcond.size());

        if (constlambda.getvalue() == true || terms[i].options[8]=="true")
          {
          fcrandom[fcrandom.size()-1].set_lambdaconst(lambda);
          fullcond.push_back(&fcrandom[fcrandom.size()-1]);
          }
        else
          {

          fullcond.push_back(&fcrandom[fcrandom.size()-1]);

          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fcrandom[fcrandom.size()-1],
                            distr[distr.size()-1],a1,b1,
                            title2,pathnonp2,pathres2,false,collinpred)
                            );

          if(updatetau)
            fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();
          if(terms[i].options[7]=="true")
                fcvarnonp[fcvarnonp.size()-1].set_uniformprior();


          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

          }

        }

      }

    }

  return false;

  }

bool bayesreg::create_random(const unsigned & collinpred)
  {

  ST::string pathnonp2;
  ST::string pathres2;
  ST::string title2;
  double a1,b1,lambda;
  int f;
  long h;
  bool iwlsmode=false;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( randomeff.checkvector(terms,i) == true )
      {

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtodouble(lambda);

      f = (terms[i].options[2]).strtodouble(a1);

      f = (terms[i].options[3]).strtodouble(b1);

      if (terms[i].options[4]=="iwlsmode")
        iwlsmode=true;


      if (f==1)
        return true;


      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_random.raw","_random.res","_random");

      make_paths(collinpred,pathnonp2,pathres2,title2,terms[i].varnames[0],"",
                   "_random_var.raw","_random_var.res","_random_variance");


      if (check_gaussian())
        {

        FULLCOND_nonp_gaussian * structuredp;
        unsigned structured=0;

        unsigned j1;
        for (j1=0;j1<fcnonpgaussian.size();j1++)
          {
          if  ( ((fcnonpgaussian[j1].get_datanames()).size() == 1) &&
                (fcnonpgaussian[j1].get_datanames()[0]==terms[i].varnames[0]) &&
                (fcnonpgaussian[j1].get_col() == collinpred) &&
                (fcnonpgaussian[j1].get_type() == MCMC::mrf)
              )
            {
            structuredp = &fcnonpgaussian[j1];
            structured ++;
            }
          }

        fcrandomgaussian.push_back(
        FULLCOND_random_gaussian(&generaloptions[generaloptions.size()-1],
                                                        distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        D.getCol(j),
                                                        title,
                                                        pathnonp,pathres,
                                                        lambda,collinpred
                                                        )
                          );

        if (structured==1)
          {

          ST::string pathnonpt = defaultpath + "\\temp\\" + name + add_name +
                 terms[i].varnames[0] +
                 "_spatialtotal.raw";

          ST::string pathrest = outfile.getvalue() + add_name + "_" + terms[i].varnames[0] +
                                "_spatialtotal.res";

          fcrandomgaussian[fcrandomgaussian.size()-1].init_spatialtotal(
          structuredp,pathnonpt,pathrest);
          }
        else if (structured==0)
          {
          }
        else
          {
          outerror("ERROR: more than one spatial effect specified for variable "
                   + terms[i].varnames[0] + "\n");
          return true;
          }

        fcrandomgaussian[fcrandomgaussian.size()-1].init_name(terms[i].varnames[0]);
        fcrandomgaussian[fcrandomgaussian.size()-1].set_fcnumber(fullcond.size());

        if (constlambda.getvalue() == true || terms[i].options[7]=="true")
          {
          fcrandomgaussian[fcrandomgaussian.size()-1].set_lambdaconst(lambda);
          fullcond.push_back(&fcrandomgaussian[fcrandomgaussian.size()-1]);
          }
        else
          {
          fullcond.push_back(&fcrandomgaussian[fcrandomgaussian.size()-1]);
          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fcrandomgaussian[fcrandomgaussian.size()-1],
                            distr[distr.size()-1],a1,b1,title2,
                            pathnonp2,pathres2,false,collinpred)
                            );
          if(terms[i].options[6]=="true")
            fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

          }

        } // end: gaussian, probit, etc. ...

      else
        {

        FULLCOND_nonp * structuredp1;
        unsigned structured1=0;

        unsigned j1;
        for (j1=0;j1<fcnonp.size();j1++)
          {
          if  ( ((fcnonp[j1].get_datanames()).size() == 1) &&
                (fcnonp[j1].get_datanames()[0]==terms[i].varnames[0]) &&
                (fcnonp[j1].get_col() == collinpred) //&&
                // (fcnonp[j1].get_type() == MCMC::mrf)
              )
            {
            structuredp1 = &fcnonp[j1];
            structured1 ++;
            }
          }


        FULLCOND_nonp_gaussian * structuredp2;
        unsigned structured2=0;

        for (j1=0;j1<fcnonpgaussian.size();j1++)
          {
          if  ( ((fcnonpgaussian[j1].get_datanames()).size() == 1) &&
                (fcnonpgaussian[j1].get_datanames()[0]==terms[i].varnames[0]) &&
                (fcnonpgaussian[j1].get_col() == collinpred) &&
                (fcnonpgaussian[j1].get_type() == MCMC::mrf)
              )
            {
            structuredp2 = &fcnonpgaussian[j1];
            structured2 ++;
            }
          }



        fcrandom.push_back( FULLCOND_random_nongaussian(
        &generaloptions[generaloptions.size()-1], distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        D.getCol(j),
                                                        title,
                                                        pathnonp,
                                                        pathres,
                                                        lambda,
                                                        iwlsmode,
                                                        collinpred
                                                        )
                          );


       if  ( ( (structured1==1) && (structured2==0) ) ||
             ( (structured1==0) && (structured2==1) )
           )
          {

          ST::string pathnonpt = defaultpath + "\\temp\\" + name + add_name +
                 terms[i].varnames[0] +
                 "_spatialtotal.raw";

          ST::string pathrest = outfile.getvalue()
                                + add_name + "_" + terms[i].varnames[0] +
                                "_spatialtotal.res";

          if (structured1==1)
            fcrandom[fcrandom.size()-1].init_spatialtotal(
            structuredp1,pathnonpt,pathrest);
          else
            fcrandom[fcrandom.size()-1].init_spatialtotal(
            structuredp2,pathnonpt,pathrest);
          }
        else if (structured1==0 && structured2==0)
          {
          }
        else
          {
          outerror("ERROR: more than one spatial effect specified for variable "
                   + terms[i].varnames[0] + "\n");
          return true;
          }


        fcrandom[fcrandom.size()-1].init_name(terms[i].varnames[0]);
        fcrandom[fcrandom.size()-1].set_fcnumber(fullcond.size());

        if (constlambda.getvalue() == true || terms[i].options[7]=="true")
          {
          fcrandom[fcrandom.size()-1].set_lambdaconst(lambda);
          fullcond.push_back(&fcrandom[fcrandom.size()-1]);
          }
        else
          {

          fullcond.push_back(&fcrandom[fcrandom.size()-1]);

          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fcrandom[fcrandom.size()-1],
                            distr[distr.size()-1],a1,b1,title2,
                            pathnonp2,pathres2,false,collinpred)
                            );
          if(terms[i].options[6]=="true")
            fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

          }

        }

      }

    }

  return false;

  }


bool bayesreg::create_mixture(const unsigned & collinpred)
  {
  int f;
  unsigned nrcomp,aclag;
  double wprior,mpriorm,mpriorv,vpriora,vpriorb;
  long h;
  bool nosamples=false,vpriorbunif=false,vpriorbgamma=false;
  ST::string ordertype="n";

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( mixtureeff.checkvector(terms,i) == true )
      {

      j = terms[i].varnames[0].isinlist(modelvarnamesv);
      f = (terms[i].options[1]).strtolong(h);
      nrcomp = unsigned(h);
      f = (terms[i].options[2]).strtodouble(wprior);
      f = (terms[i].options[3]).strtodouble(mpriorm);
      f = (terms[i].options[4]).strtodouble(mpriorv);
      f = (terms[i].options[5]).strtodouble(vpriora);
      f = (terms[i].options[6]).strtodouble(vpriorb);
      if(terms[i].options[7] == "true")
        nosamples = true;
      f = (terms[i].options[8]).strtolong(h);
      aclag = unsigned(h);
      if(terms[i].options[9] == "w")
        ordertype = "w";
      if(terms[i].options[10] == "true")
        vpriorbunif = true;
      if(terms[i].options[11] == "true")
        vpriorbgamma = true;

      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_mixture.raw","_mixture.res","_mixture");

      if (check_gaussian())
        {


           fcmixture.push_back(
           FULLCOND_mixture(&generaloptions[generaloptions.size()-1],
                                                        distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        D.getCol(j),
                                                        title,
                                                        pathnonp,pathres,
                                                        nrcomp,wprior,
                                                        mpriorm,mpriorv,
                                                        vpriora,vpriorb,
                                                        nosamples,aclag,
                                                        ordertype,
                                                        vpriorbunif,vpriorbgamma,
                                                        collinpred
                                                        )
                          );

          fcmixture[fcmixture.size()-1].init_name(terms[i].varnames[0]);
          fcmixture[fcmixture.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcmixture[fcmixture.size()-1]);
        } // end: gaussian,???

      else
        {
        outerror("ERROR: only family=gaussian allowed for variable "
                   + terms[i].varnames[0] + "\n");
        }

      }

    }

  return false;

  }


bool bayesreg::create_interactionspspline(const unsigned & collinpred)
  {

  ST::string proposal;

  long h;
  double lambda;
  bool reduced, singleblock;
  unsigned min,max,nrknots,degree,blocksize;
  double a1,b1;
  int gridsize;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpinteractpspline.checkvector(terms,i) == true )
      {

      MCMC::fieldtype type;
      if ((terms[i].options[0] == "pspline2dimrw1")   ||
          (terms[i].options[0] == "tpspline2dimrw1")
         )
        type = MCMC::mrflinear;
      else if (terms[i].options[0] == "pspline2dimrw2")
        type = MCMC::mrfquadratic8;
      else if ((terms[i].options[0] == "pspline2dimband")   ||
          (terms[i].options[0] == "tpspline2dimband")
         )
        type = MCMC::mrflinearband;
      else if (terms[i].options[0] == "psplinekrrw1")
        type = MCMC::mrfkr1;
      else if (terms[i].options[0] == "psplinekrrw2")
        type = MCMC::mrfkr2;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      min = unsigned(h);

      f = (terms[i].options[2]).strtolong(h);
      max = unsigned(h);

      f = (terms[i].options[3]).strtolong(h);
      degree = unsigned(h);

      f = (terms[i].options[4]).strtolong(h);
      nrknots = unsigned(h);

      f = (terms[i].options[5]).strtodouble(lambda);

      if (terms[i].options[6] == "false")
        reduced = false;
      else
        reduced = true;

      f = (terms[i].options[7]).strtodouble(a1);

      f = (terms[i].options[8]).strtodouble(b1);

      if (terms[i].options[9] == "false")
        singleblock = false;
      else
        singleblock = true;

      f = (terms[i].options[10]).strtolong(h);
      gridsize = unsigned(h);

      proposal = terms[i].options[11];

      MCMC::knotpos po;

      if (knots.getvalue() == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      if (f==1)
        return true;

      ST::string help  = terms[i].varnames[0] + "_" + terms[i].varnames[1];

      make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_pspline.raw","_pspline.res","_pspline");


      if (check_gaussian())
        {

        FULLCOND_pspline_gaussian * mainp1;
        FULLCOND_pspline_gaussian * mainp2;
        unsigned main1=0;
        unsigned main2=0;

        unsigned j;
        for (j=0;j<fcpsplinegaussian.size();j++)
          {
          if  ( ((fcpsplinegaussian[j].get_datanames()).size() == 1) &&
               (fcpsplinegaussian[j].get_datanames()[0] == terms[i].varnames[0]) &&
                fcpsplinegaussian[j].get_col() == collinpred
              )
            {
            mainp1 = &fcpsplinegaussian[j];
            if(mainp1->get_gridsize() != gridsize)
              {
              outerror("ERROR: gridsize for interaction term and main effects must be the same\n");
              return true;
              }
            if(mainp1->get_nrknots() != nrknots)
              {
              outerror("ERROR: number of knots for interaction term and main effects must be the same\n");
              return true;
              }
            if(mainp1->get_degree() != degree)
              {
              outerror("ERROR: degree for interaction term and main effects must be the same\n");
              return true;
              }
            main1 ++;
            }


          if  ( ((fcpsplinegaussian[j].get_datanames()).size() == 1) &&
              (fcpsplinegaussian[j].get_datanames()[0] == terms[i].varnames[1]) &&
              fcpsplinegaussian[j].get_col() == collinpred
              )
            {
            mainp2 = &fcpsplinegaussian[j];
            if(mainp2->get_gridsize() != gridsize)
              {
              outerror("ERROR: gridsize for interaction term and main effects must be the same\n");
              return true;
              }
            if(mainp2->get_nrknots() != nrknots)
              {
              outerror("ERROR: number of knots for interaction term and main effects must be the same\n");
              return true;
              }
            if(mainp2->get_degree() != degree)
              {
              outerror("ERROR: degree for interaction term and main effects must be the same\n");
              return true;
              }
            main2 ++;
            }

          }

        fcpsplinesurfgaussian.push_back(
        FULLCOND_pspline_surf_gaussian(&generaloptions[generaloptions.size()-1],distr[distr.size()-1],
                                      fcconst_intercept,
                                      D.getCol(j1),
                                      D.getCol(j2),
                                      title,
                                      nrknots,degree,po,
                                      lambda,
                                      gridsize,
                                      type,
                                      pathnonp,
                                      pathres,
                                      outfile.getvalue(),
                                      singleblock,
                                      collinpred
                                      ));

        if (constlambda.getvalue() == true)
          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_lambdaconst(lambda);

        if ( (main1==1) && (main2==1) )
          {

          ST::string pathnonpt;
          ST::string pathrest;
          ST::string titlet;

          make_paths(collinpred,pathnonpt,pathrest,titlet,help,"",
                 "_pspline_total.raw","_pspline_total.res","_pspline_total");

          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].
          init_maineffects(mainp1,mainp2,pathnonpt,pathrest);
          mainp1->set_interaction();
          mainp2->set_interaction();
          }
        else if ( (main1==0) && (main2==0) )
          {
          }
        else
          {
          // FEHLERMELDUNG
          }

        vector<ST::string> na;
        na.push_back(terms[i].varnames[0]);
        na.push_back(terms[i].varnames[1]);

        fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].init_names(na);

        make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_pspline_var.raw","_pspline_var.res","_pspline_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                       &fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );

        bool rowwise = false;
        if(terms[i].options[0] == "tpspline2dimband")
          rowwise = true;

        if ((terms[i].options[0] == "tpspline2dimrw1") ||
            (terms[i].options[0] == "tpspline2dimband")
           )
          {

          make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_pspline_tvar.raw","_pspline_tvar.res","_pspline_tvariance");

          unsigned v = nu.getvalue();
          f = (terms[i].options[16]).strtolong(h);
          blocksize = unsigned(h);

          fctvariance2dim.push_back(FULLCOND_tvariance2dim(&generaloptions[generaloptions.size()-1],
          &fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1],v,title,
          pathnonp,pathres,blocksize,rowwise));

          }

        fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1]);

        if (constlambda.getvalue() == false)
          {
          if(terms[i].options[15]=="true")
              fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
          }

        if (terms[i].options[0] == "tpspline2dimrw1" ||
           (terms[i].options[0] == "tpspline2dimband")
           )
          {
          fctvariance2dim[fctvariance2dim.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fctvariance2dim[fctvariance2dim.size()-1]);
          }

        }
      else
        {
        // NONGAUSSIAN CASE
        if(proposal=="iwls" || proposal=="iwlsmode")
          {
          // IWLS

          IWLS_pspline * mainp1;
          IWLS_pspline * mainp2;
          unsigned main1=0;
          unsigned main2=0;

          bool iwlsmode;
          if(proposal=="iwlsmode")
            iwlsmode = true;
          else
            iwlsmode = false;
          f = (terms[i].options[12]).strtolong(h);
          unsigned updateW;
          updateW = unsigned(h);
          bool updatetau;
          if(terms[i].options[13] == "false" || constlambda.getvalue() == true)
            updatetau = false;
          else
            updatetau = true;
          double fstart;
            f = (terms[i].options[14]).strtodouble(fstart);

          unsigned j;
          for (j=0;j<fciwlspspline.size();j++)
            {
            if  ( ((fciwlspspline[j].get_datanames()).size() == 1) &&
                   (fciwlspspline[j].get_datanames()[0] == terms[i].varnames[0]) &&
                   fciwlspspline[j].get_col() == collinpred
                )
              {
              mainp1 = &fciwlspspline[j];
              if(mainp1->get_gridsize() != gridsize)
                {
                outerror("ERROR: gridsize for interaction term and main effects must be the same\n");
                return true;
                }
              if(mainp1->get_nrknots() != nrknots)
                {
                outerror("ERROR: number of knots for interaction term and main effects must be the same\n");
                return true;
                }
              if(mainp1->get_degree() != degree)
                {
                outerror("ERROR: degree for interaction term and main effects must be the same\n");
                return true;
                }
              main1 ++;
              }


            if  ( ((fciwlspspline[j].get_datanames()).size() == 1) &&
                (fciwlspspline[j].get_datanames()[0] == terms[i].varnames[1]) &&
                fciwlspspline[j].get_col() == collinpred
                )
              {
              mainp2 = &fciwlspspline[j];
              if(mainp2->get_gridsize() != gridsize)
                {
                outerror("ERROR: gridsize for interaction term and main effects must be the same\n");
                return true;
                }
              if(mainp2->get_nrknots() != nrknots)
                {
                outerror("ERROR: number of knots for interaction term and main effects must be the same\n");
                return true;
                }
              if(mainp2->get_degree() != degree)
                {
                outerror("ERROR: degree for interaction term and main effects must be the same\n");
                return true;
                }
              main2 ++;
              }

            }



          fcpsplinesurfgaussian.push_back( FULLCOND_pspline_surf_gaussian(&generaloptions[generaloptions.size()-1],
                                         distr[distr.size()-1],
                                         fcconst_intercept,
                                         D.getCol(j1),
                                         D.getCol(j2),
                                         iwlsmode,
                                         title,
                                         nrknots,degree,
                                         po,
                                         lambda,
                                         updateW,
                                         updatetau,
                                         fstart,
                                         a1,b1,
                                         gridsize,
                                         type,
                                         pathnonp,
                                         pathres,
                                         outfile.getvalue(),
                                         true,
                                         singleblock,
                                         collinpred
                                        )
                        );

          if (constlambda.getvalue() == true)
            fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_lambdaconst(lambda);

          if ( (main1==1) && (main2==1) )
            {

            ST::string pathnonpt;
            ST::string pathrest;
            ST::string titlet;

            make_paths(collinpred,pathnonpt,pathrest,titlet,help,"",
                 "_pspline_total.raw","_pspline_total.res","_pspline_total");

            fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].
            init_maineffects(mainp1,mainp2,pathnonpt,pathrest);
            mainp1->set_interaction();
            mainp2->set_interaction();
            }
          else if ( (main1==0) && (main2==0) )
            {
            }
          else
            {
            // FEHLERMELDUNG
            }


          vector<ST::string> na;
          na.push_back(terms[i].varnames[0]);
          na.push_back(terms[i].varnames[1]);

          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].init_names(na);

          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_pspline_var.raw","_pspline_var.res","_pspline_variance");

          fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                       &fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );

          bool rowwise = false;
          if(terms[i].options[0] == "tpspline2dimband")
            rowwise = true;

          if ((terms[i].options[0] == "tpspline2dimrw1") ||
              (terms[i].options[0] == "tpspline2dimband")
             )
            {

            make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_pspline_tvar.raw","_pspline_tvar.res","_pspline_tvariance");

            unsigned v = nu.getvalue();
            f = (terms[i].options[16]).strtolong(h);
            blocksize = unsigned(h);

            fctvariance2dim.push_back(FULLCOND_tvariance2dim(&generaloptions[generaloptions.size()-1],
            &fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1],v,title,
            pathnonp,pathres,blocksize,rowwise));

            }

          if (constlambda.getvalue() == false)
            {
            if(terms[i].options[15]=="true")
              fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
            if(updatetau)
              fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();
            fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
            }

          if (terms[i].options[0] == "tpspline2dimrw1")
            {
            fctvariance2dim[fctvariance2dim.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fctvariance2dim[fctvariance2dim.size()-1]);
            }

          }
        else
          {
          // CONDITIONAL PRIOR

          if (terms[i].options[0] == "tpspline2dimrw1")
            {
            outerror("ERROR: '" + terms[i].options[0] + "' not available\n");
            return true;
            }

          FULLCOND_pspline * mainp1;
          FULLCOND_pspline * mainp2;
          unsigned main1=0;
          unsigned main2=0;

          unsigned j;
          for (j=0;j<fcpspline.size();j++)
            {
            if  ( ((fcpspline[j].get_datanames()).size() == 1) &&
                 (fcpspline[j].get_datanames()[0] == terms[i].varnames[0]) &&
                  fcpspline[j].get_col() == collinpred
                )
              {
              mainp1 = &fcpspline[j];
              if(mainp1->get_gridsize() != gridsize)
                {
                outerror("ERROR: gridsize for interaction term and main effects must be the same\n");
                return true;
                }
              if(mainp1->get_nrknots() != nrknots)
                {
                outerror("ERROR: number of knots for interaction term and main effects must be the same\n");
                return true;
                }
              if(mainp1->get_degree() != degree)
                {
                outerror("ERROR: degree for interaction term and main effects must be the same\n");
                return true;
                }
              main1 ++;
              }


            if  ( ((fcpspline[j].get_datanames()).size() == 1) &&
                (fcpspline[j].get_datanames()[0] == terms[i].varnames[1]) &&
                fcpspline[j].get_col() == collinpred
                )
              {
              mainp2 = &fcpspline[j];
              if(mainp2->get_gridsize() != gridsize)
                {
                outerror("ERROR: gridsize for interaction term and main effects must be the same\n");
                return true;
                }
              if(mainp2->get_nrknots() != nrknots)
                {
                outerror("ERROR: number of knots for interaction term and main effects must be the same\n");
                return true;
                }
              if(mainp2->get_degree() != degree)
                {
                outerror("ERROR: degree for interaction term and main effects must be the same\n");
                return true;
                }
              main2 ++;
              }

            }

          fcpsplinesurf.push_back( FULLCOND_pspline_surf(&generaloptions[generaloptions.size()-1],
                                         distr[distr.size()-1],
                                         fcconst_intercept,
                                         D.getCol(j1),
                                         D.getCol(j2),
                                         title,
                                         nrknots,degree,
                                         po,
                                         min,max,
                                         lambda,
                                         gridsize,
                                         type,
                                         pathnonp,
                                         pathres,
                                         outfile.getvalue(),
                                         collinpred
                                        )
                        );

          if (constlambda.getvalue() == true)
            fcpsplinesurf[fcpsplinesurf.size()-1].set_lambdaconst(lambda);

          if ( (main1==1) && (main2==1) )
            {

            ST::string pathnonpt;
            ST::string pathrest;
            ST::string titlet;

            make_paths(collinpred,pathnonpt,pathrest,titlet,help,"",
                   "_pspline_total.raw","_pspline_total.res","_pspline_total");

            fcpsplinesurf[fcpsplinesurf.size()-1].
            init_maineffects(mainp1,mainp2,pathnonpt,pathrest);
            mainp1->set_interaction();
            mainp2->set_interaction();
            }
          else if ( (main1==0) && (main2==0) )
            {
            }
          else
            {
            // FEHLERMELDUNG
            }

          vector<ST::string> na;
          na.push_back(terms[i].varnames[0]);
          na.push_back(terms[i].varnames[1]);

          fcpsplinesurf[fcpsplinesurf.size()-1].init_names(na);

          fcpsplinesurf[fcpsplinesurf.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcpsplinesurf[fcpsplinesurf.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_pspline_var.raw","_pspline_var.res","_variance");

          fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                       &fcpsplinesurf[fcpsplinesurf.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );


          if (constlambda.getvalue() == false)
            {
            if(terms[i].options[15]=="true")
              fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
            fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
            }

          }

        }

      }

    }


  return false;

  }


bool bayesreg::create_geospline(const unsigned & collinpred)
  {

  ST::string proposal;

  long h;
  double lambda;
  double a1,b1;
  bool reduced,singleblock;
  unsigned min,max,nrknots,degree;
  int gridsize;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpgeospline.checkvector(terms,i) == true )
      {

      MCMC::fieldtype type;
      if ((terms[i].options[0] == "geospline") || (terms[i].options[0] == "geosplinerw1"))
        type = MCMC::mrflinear;
      else if (terms[i].options[0] == "geosplinerw2")
        type = MCMC::mrfquadratic8;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      min = unsigned(h);

      f = (terms[i].options[2]).strtolong(h);
      max = unsigned(h);

      f = (terms[i].options[3]).strtolong(h);
      degree = unsigned(h);

      f = (terms[i].options[4]).strtolong(h);
      nrknots = unsigned(h);

      f = (terms[i].options[5]).strtodouble(lambda);

      if (terms[i].options[6] == "false")
        reduced = false;
      else
        reduced = true;

      mapobject * mapp;                           // pointer to mapobject

      int objpos = findstatobject(*statobj,terms[i].options[7],"map");

      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          outerror("ERROR: map object " + terms[i].options[7] + " is not existing\n");
        else
          outerror("ERROR: " + terms[i].options[7] + " is not a map object\n");
        return true;
        }

      MAP::map m = mapp->getmap();

      if(!m.centroids_existing())
        {
        outerror("ERROR: map object doesn´t contain centroids\n");
        return true;
        }

      if (terms[i].options[8] == "false")
        singleblock = false;
      else
        singleblock = true;

      f = (terms[i].options[9]).strtodouble(a1);
      f = (terms[i].options[10]).strtodouble(b1);

      proposal = terms[i].options[11];

      gridsize = -1;

      MCMC::knotpos po;

      if (knots.getvalue() == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_geospline.raw","_geospline.res","_geospline");

      if (check_gaussian())
        {

        fcpsplinesurfgaussian.push_back(
        FULLCOND_pspline_surf_gaussian(&generaloptions[generaloptions.size()-1],
                                      distr[distr.size()-1],fcconst_intercept,
                                      D.getCol(j),m,terms[i].options[7],
                                      title,
                                      nrknots,degree,po,
                                      lambda,
                                      gridsize,
                                      type,
                                      pathnonp,
                                      pathres,
                                      singleblock,
                                      collinpred
                                      ));

        if (constlambda.getvalue() == true)
          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_lambdaconst(lambda);

        vector<ST::string> na;
        na.push_back(terms[i].varnames[0]);

        fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].init_names(na);

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_geospline_var.raw","_geospline_var.res","_geospline_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                      &fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );

        fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1]);

        if (constlambda.getvalue() == false)
          {
          if(terms[i].options[15]=="true")
            fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
          }

        }
      else
        {
        // NONGAUSSIAN CASE
        if(proposal=="iwls" || proposal=="iwlsmode")
          {
          // IWLS
          bool iwlsmode;
          if(proposal=="iwlsmode")
            iwlsmode = true;
          else
            iwlsmode = false;
          f = (terms[i].options[12]).strtolong(h);
          unsigned updateW;
          updateW = unsigned(h);
          bool updatetau;
          if(terms[i].options[13] == "false" || constlambda.getvalue() == true)
            updatetau = false;
          else
            updatetau = true;
          double fstart;
            f = (terms[i].options[14]).strtodouble(fstart);

          fcpsplinesurfgaussian.push_back(
          FULLCOND_pspline_surf_gaussian(&generaloptions[generaloptions.size()-1],
                                      distr[distr.size()-1],fcconst_intercept,
                                      D.getCol(j),m,terms[i].options[7],
                                      iwlsmode,
                                      title,
                                      nrknots,degree,po,
                                      lambda,
                                      updateW,
                                      updatetau,
                                      fstart,
                                      a1,b1,
                                      gridsize,
                                      type,
                                      pathnonp,
                                      pathres,
                                      true,
                                      singleblock,
                                      collinpred
                                      ));

          if (constlambda.getvalue() == true)
            fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_lambdaconst(lambda);

          vector<ST::string> na;
          na.push_back(terms[i].varnames[0]);
          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].init_names(na);
          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_geospline_var.raw","_geospline_var.res","_geospline_variance");

          fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                        &fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );

          if (constlambda.getvalue() == false)
            {
            if(updatetau)
              fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();
            if(terms[i].options[15]=="true")
              fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
            fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
            }

          }
        else
          {
          // CONDITIONAL PRIOR
          fcpsplinesurf.push_back(
          FULLCOND_pspline_surf(&generaloptions[generaloptions.size()-1],distr[distr.size()-1],
                                      fcconst_intercept,
                                      D.getCol(j),m,terms[i].options[7],
                                      title,
                                      nrknots,degree,po,
                                      min,max,
                                      lambda,
                                      gridsize,
                                      type,
                                      pathnonp,
                                      pathres,
                                      collinpred
                                      ));

          if (constlambda.getvalue() == true)
            fcpsplinesurf[fcpsplinesurf.size()-1].set_lambdaconst(lambda);

          vector<ST::string> na;
          na.push_back(terms[i].varnames[0]);
          fcpsplinesurf[fcpsplinesurf.size()-1].init_names(na);
          fcpsplinesurf[fcpsplinesurf.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcpsplinesurf[fcpsplinesurf.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_geospline_var.raw","_geospline_var.res","_geospline_variance");

          fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                        &fcpsplinesurf[fcpsplinesurf.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );

          if (constlambda.getvalue() == false)
            {
            if(terms[i].options[15]=="true")
              fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
            fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
            }

          }

        }

      }

    }

  return false;

  }


bool bayesreg::create_varcoeff_geospline(const unsigned & collinpred)
  {

  ST::string proposal;

  long h;
  double lambda;
  double a1,b1;
  bool reduced,singleblock;
  unsigned min,max,nrknots,degree;
  int gridsize;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpvarcoeffgeospline.checkvector(terms,i) == true )
      {

      MCMC::fieldtype type = MCMC::mrflinear;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interaction var
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      min = unsigned(h);

      f = (terms[i].options[2]).strtolong(h);
      max = unsigned(h);

      f = (terms[i].options[3]).strtolong(h);
      degree = unsigned(h);

      f = (terms[i].options[4]).strtolong(h);
      nrknots = unsigned(h);

      f = (terms[i].options[5]).strtodouble(lambda);

      if (terms[i].options[6] == "false")
        reduced = false;
      else
        reduced = true;

      mapobject * mapp;                           // pointer to mapobject

      int objpos = findstatobject(*statobj,terms[i].options[7],"map");

      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          outerror("ERROR: map object " + terms[i].options[7] + " is not existing\n");
        else
          outerror("ERROR: " + terms[i].options[7] + " is not a map object\n");
        return true;
        }

      MAP::map m = mapp->getmap();

      if(!m.centroids_existing())
        {
        outerror("ERROR: map object doesn´t contain centroids\n");
        return true;
        }

      if (terms[i].options[8] == "false")
        singleblock = false;
      else
        singleblock = true;

      f = (terms[i].options[9]).strtodouble(a1);
      f = (terms[i].options[10]).strtodouble(b1);

      proposal = terms[i].options[11];

      gridsize = -1;

      MCMC::knotpos po;

      if (knots.getvalue() == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      if(f==1)
        return false;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],terms[i].varnames[0],
                 "_geospline.raw","_geospline.res","_geospline");

      if (check_gaussian())
        {

        fcpsplinesurfgaussian.push_back(
        FULLCOND_pspline_surf_gaussian(&generaloptions[generaloptions.size()-1],distr[distr.size()-1],
                                      fcconst_intercept,
                                      D.getCol(j1),D.getCol(j2),m,terms[i].options[7],
                                      title,
                                      nrknots,degree,po,
                                      lambda,
                                      gridsize,
                                      type,
                                      pathnonp,
                                      pathres,
                                      singleblock,
                                      collinpred
                                      ));

        if(constlambda.getvalue() == true)
          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_lambdaconst(lambda);

        vector<ST::string> na;
        na.push_back(terms[i].varnames[1]);
        na.push_back(terms[i].varnames[0]);

        fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].init_names(na);

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],terms[i].varnames[1],
                   "_geospline_var.raw","_geospline_var.res","_geospline_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                      &fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );

        fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1]);

        if (constlambda.getvalue() == false)
          {
          if(terms[i].options[15]=="true")
            fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
          }

        }
      else
        {
        // NONGAUSSIAN CASE
        if(proposal=="iwls" || proposal=="iwlsmode")
          {
          // IWLS
          bool iwlsmode;
          if(proposal=="iwlsmode")
            iwlsmode = true;
          else
            iwlsmode = false;
          f = (terms[i].options[12]).strtolong(h);
          unsigned updateW;
          updateW = unsigned(h);
          bool updatetau;
          if(terms[i].options[13] == "false")
            updatetau = false;
          else
            updatetau = true;
          double fstart;
            f = (terms[i].options[14]).strtodouble(fstart);

          fcpsplinesurfgaussian.push_back(
          FULLCOND_pspline_surf_gaussian(&generaloptions[generaloptions.size()-1],distr[distr.size()-1],
                                      fcconst_intercept,
                                      D.getCol(j1),D.getCol(j2),m,terms[i].options[7],
                                      title,
                                      nrknots,degree,po,
                                      lambda,
                                      updateW,
                                      updatetau,
                                      fstart,
                                      a1,b1,
                                      gridsize,
                                      type,
                                      pathnonp,
                                      pathres,
                                      iwlsmode,
                                      singleblock,
                                      collinpred
                                      ));

          if (constlambda.getvalue() == true)
            fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_lambdaconst(lambda);

          vector<ST::string> na;
          na.push_back(terms[i].varnames[1]);
          na.push_back(terms[i].varnames[0]);
          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].init_names(na);
          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_geospline_var.raw","_geospline_var.res","_geospline_variance");

          fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                        &fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );

          if (constlambda.getvalue() == false)
            {
            if(updatetau)
              fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();
            if(terms[i].options[15]=="true")
              fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
            fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
            }

          }
        else
          {
          // CONDITIONAL PRIOR
          fcpsplinesurf.push_back(
          FULLCOND_pspline_surf(&generaloptions[generaloptions.size()-1],distr[distr.size()-1],
                                      fcconst_intercept,
                                      D.getCol(j1),D.getCol(j2),m,terms[i].options[7],
                                      title,
                                      nrknots,degree,po,
                                      min,max,
                                      lambda,
                                      gridsize,
                                      type,
                                      pathnonp,
                                      pathres,
                                      collinpred
                                      ));

          if (constlambda.getvalue() == true)
            fcpsplinesurf[fcpsplinesurf.size()-1].set_lambdaconst(lambda);

          vector<ST::string> na;
          na.push_back(terms[i].varnames[1]);
          na.push_back(terms[i].varnames[0]);
          fcpsplinesurf[fcpsplinesurf.size()-1].init_names(na);
          fcpsplinesurf[fcpsplinesurf.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcpsplinesurf[fcpsplinesurf.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_geospline_var.raw","_geospline_var.res","_geospline_variance");

          fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                        &fcpsplinesurf[fcpsplinesurf.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );


          if (constlambda.getvalue() == false)
            {
            if(terms[i].options[15]=="true")
              fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
            fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
            }

          }

        }

      }

    }

  return false;

  }


bool bayesreg::create_baseline(const unsigned & collinpred)
  {

  ST::string proposal;

  long h;
  unsigned min,max,degree,nrknots;
  double lambda,a1,b1;
  bool ub, wb;
  int gridsize;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( baseline.checkvector(terms,i) == true )
      {

      // --------------- reading options, term information ---------------------

      MCMC::fieldtype type;
      type = MCMC::RW2;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      min = unsigned(h);

      f = (terms[i].options[2]).strtolong(h);
      max = unsigned(h);

      f = (terms[i].options[3]).strtolong(h);
      degree = unsigned(h);

      f = (terms[i].options[4]).strtolong(h);
      nrknots = unsigned(h);

      f = (terms[i].options[5]).strtodouble(lambda);

      f = (terms[i].options[6]).strtodouble(a1);

      f = (terms[i].options[7]).strtodouble(b1);

      if (terms[i].options[8] == "false")
        ub = false;
      else
        ub = true;

      f = (terms[i].options[9]).strtolong(h);
      gridsize = unsigned(h);

      MCMC::knotpos po;

      if (knots.getvalue() == "equidistant" && terms[i].options[14] == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      proposal = terms[i].options[11];

      if (terms[i].options[12] == "false")
        wb = false;
      else
        wb = true;

      if (f==1)
        return true;

      datamatrix beg;
      if (begin.getvalue() == "")
        beg = datamatrix(1,1);
      else
        beg = D.getCol(begpos);

      // -------------end: reading options, term information -------------------

      //--------- creating path for samples and and results, creating title ----

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_logbaseline.raw","_logbaseline.res","_logbaseline");

      //----- end: creating path for samples and and results, creating title ---


      if(proposal == "cp")
        {
        fcbaseline.push_back( pspline_baseline(&generaloptions[generaloptions.size()-1],
                                                distr[distr.size()-1],
                                                fcconst_intercept,
                                                D.getCol(j),
                                                a1,
                                                b1,
                                                nrknots,
                                                degree,
                                                po,
                                                lambda,
                                                min,
                                                max,
                                                type,
                                                title,
                                                pathnonp,
                                                pathres,
                                                gridsize,
                                                collinpred,
                                                beg,
                                                wb
                                               )
                             );

        if (constlambda.getvalue() == true)
          fcbaseline[fcbaseline.size()-1].set_lambdaconst(lambda);

        fcbaseline[fcbaseline.size()-1].init_name(terms[i].varnames[0]);
        fcbaseline[fcbaseline.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcbaseline[fcbaseline.size()-1]);

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_logbaseline_var.raw","_logbaseline_var.res","_logbaseline_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                &fcbaseline[fcbaseline.size()-1],
                                distr[distr.size()-1],a1,b1,
                                title,pathnonp,pathres,ub,collinpred)
                                );

        }
      else if(proposal == "iwls")
        {

/*        fcbaselineiwls.push_back( IWLS_baseline(&generaloptions[generaloptions.size()-1],
                                                distr[distr.size()-1],
                                                fcconst_intercept,
                                                D.getCol(j),
                                                false,
                                                nrknots,
                                                degree,
                                                po,
                                                lambda,
                                                type,
                                                "unrestricted",
                                                1,
                                                false,
                                                2,
                                                a1,
                                                b1,
                                                title,
                                                pathnonp,
                                                pathres,
                                                false,
                                                gridsize,
                                                false,
                                                collinpred,
                                                beg
                                               )
                             );

        if (constlambda.getvalue() == true)
          fcbaselineiwls[fcbaselineiwls.size()-1].set_lambdaconst(lambda);

        fcbaselineiwls[fcbaselineiwls.size()-1].init_name(terms[i].varnames[0]);
        fcbaselineiwls[fcbaselineiwls.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcbaselineiwls[fcbaselineiwls.size()-1]);

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_logbaseline_var.raw","_logbaseline_var.res","_logbaseline_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                &fcbaselineiwls[fcbaselineiwls.size()-1],
                                distr[distr.size()-1],a1,b1,
                                title,pathnonp,pathres,ub,collinpred)
                                );   */
        }

      if (constlambda.getvalue() == false)
        {
        if(terms[i].options[10]=="true")
          fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
        }
      }

    }

  return false;
  }

