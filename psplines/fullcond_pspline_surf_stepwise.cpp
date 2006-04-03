
#include "first.h"

#include "fullcond_pspline_surf_stepwise.h"

namespace MCMC
{


void FULLCOND_pspline_surf_stepwise::init_maineffects(
                        FULLCOND_pspline_stepwise * mp1,FULLCOND_pspline_stepwise * mp2,
                        const ST::string & pnt,const ST::string & prt)
  {

  mainpoi1 = mp1;
  mainpoi2 = mp2;

  interaction = true;

  centertotal = false;
  maineffectsexisting = 11;

  fctotalrespath = prt;

  datamatrix h(1,1,0);
  if(gridsize < 0)
    fctotal = FULLCOND(optionsp,h,title+"total",nrdiffobs,1,pnt);
  else
    fctotal = FULLCOND(optionsp,h,title+"total",gridsize,1,pnt);
  fctotal.setflags(MCMC::norelchange | MCMC::nooutput);
  fctotal.set_transform(transform);

  beta1 = datamatrix(nrpar1dim,1,0);
  beta2 = datamatrix(nrpar1dim,1,0);
  he1 = datamatrix(xv.size(),1,0);
  he2 = datamatrix(yv.size(),1,0);

  }


void FULLCOND_pspline_surf_stepwise::init_maineffect(FULLCOND_pspline_stepwise * mp1,
         const ST::string & pnt,const ST::string & prt, const unsigned & number)
  {
  interaction = true;

  centertotal = false;
  centerboth = number;  // neu: normalerweise = 0, hier 1 oder 2, gibt damit an,
                        // daß nur ein Haupteffekt vorhanden ist und welcher!

  fctotalrespath = prt;

  datamatrix h(1,1,0);
  if(gridsize < 0)
    fctotal = FULLCOND(optionsp,h,title+"total",nrdiffobs,1,pnt);
  else
    fctotal = FULLCOND(optionsp,h,title+"total",gridsize,1,pnt);
  fctotal.setflags(MCMC::norelchange | MCMC::nooutput);
  fctotal.set_transform(transform);

  if(number == 1)
    {
    maineffectsexisting = 1;
    mainpoi1 = mp1;
    beta1 = datamatrix(nrpar1dim,1,0);
    he1 = datamatrix(xv.size(),1,0);
    }
  else
    {
    maineffectsexisting = 10;
    mainpoi2 = mp1;
    beta2 = datamatrix(nrpar1dim,1,0);
    he2 = datamatrix(yv.size(),1,0);
    }
  }


void FULLCOND_pspline_surf_stepwise::search_maineffects(void)
  {
  centerboth = 0;
  centerfix = 0;
  centertotal = false;

  bool h1 = false;                    
  bool f1 = false;
  bool h2 = false;
  bool f2 = false;
  if(maineffectsexisting == 1 || maineffectsexisting == 11)
    mainpoi1->get_inthemodel(h1,f1);
  if(maineffectsexisting == 10 || maineffectsexisting == 11)
    mainpoi2->get_inthemodel(h2,f2);
  if(h1 == false && h2 == false)
    {
    if(f1 == false && f2 == false)
      centertotal = true;
    else if(f1 == true && f2 == true)
      centerfix = 11;
    else if(f1 == false && f2 == true)
      centerfix = 10;
    else
      centerfix = 1;
    }
  else if(h1 == true && h2 == false)
    {
    centerboth = 1;
    if(f2 == true)
      centerfix = 10;
    }
  else if(h1 == false && h2 == true)
    {
    centerboth = 10;
    if(f1 == true)
      centerfix = 1;
    }
  else if(h1 == true && h2 == true)
    centerboth = 11;
  }


void FULLCOND_pspline_surf_stepwise::create(const datamatrix & v1, const datamatrix & v2, const datamatrix & intact)
  {

  unsigned i;
  data_forfixed = datamatrix(v1.rows(),1,0);

  // neu:
  for(i=0;i<v2.rows();i++)
    {
    data_forfixed(i,0) = v1(i,0) * v2(i,0);
    }

  }

  // CONSTRUCTOR

FULLCOND_pspline_surf_stepwise::FULLCOND_pspline_surf_stepwise(MCMCoptions * o,
                      DISTRIBUTION * dp, FULLCOND_const * fcc,
                      const datamatrix & v1, const datamatrix & v2,
                      const ST::string & ti,
                      const unsigned & nrk, const unsigned & degr,
                      const knotpos & kp, const double & l, const int & gs,
                      const fieldtype & ft, const ST::string & fp,
                      const ST::string & pres, const ST::string & of, const bool & gauss, const unsigned & c)
  : FULLCOND_pspline_surf_gaussian(o,dp,fcc,v1,v2,ti,nrk,degr,kp,l,gs,ft,fp,pres,of,true,c)
  {

  if(gauss==false)
    utype = iwls;

  create(v1,v2);

  if(type == mrflinear)
    grenzfall = 0;
  else if(type == mrfquadratic8)
    grenzfall = 1;

  interhaupt = 0.0;
  interhaupt1 = 0.0;
  interhaupt2 = 0.0;

  centerboth = 0;
  maineffectsexisting = 0;
  centerfix = 0;

  }


  // CONSTRUCTOR 3: geosplines

FULLCOND_pspline_surf_stepwise::FULLCOND_pspline_surf_stepwise(
                      MCMCoptions * o, DISTRIBUTION * dp,
                      FULLCOND_const * fcc,
                      const datamatrix & region, const MAP::map & mp, const ST::string & mn,
                      const ST::string & ti, const unsigned & nrk,
                      const unsigned & degr, const knotpos & kp, const double & l,
                      const int & gs, const fieldtype & ft, const ST::string & fp,
                      const ST::string & pres, const bool & gauss, const unsigned & c)
  : FULLCOND_pspline_surf_gaussian(o,dp,fcc,region,mp,mn,ti,nrk,degr,kp,l,gs,ft,fp,pres,true,c)
  {

  if(gauss==false)
    utype = iwls;

  datamatrix v1 = datamatrix(likep->get_nrobs(),1,0.0);
  datamatrix v2 = datamatrix(likep->get_nrobs(),1,0.0);

  ST::string regname;
  for(unsigned i=0;i<likep->get_nrobs();i++)
    {
    regname = ST::doubletostring(region(i,0));
    regionnames.push_back(regname);
    v1(i,0) = m.get_region(m.getnr(regname)).get_xcenter();
    v2(i,0) = m.get_region(m.getnr(regname)).get_ycenter();
    }

  create(v1,v2);

  if(type == mrflinear)
    grenzfall = 0;
  else if(type == mrfquadratic8)
    grenzfall = 1;

  }


  // CONSTRUCTOR 5: varying coefficients

FULLCOND_pspline_surf_stepwise::FULLCOND_pspline_surf_stepwise(MCMCoptions * o,
                      DISTRIBUTION * dp, FULLCOND_const * fcc,
                      const datamatrix & intact,
                      const datamatrix & v1, const datamatrix & v2,
                      const ST::string & ti, const unsigned & nrk, const unsigned & degr,
                      const knotpos & kp, const double & l, const int & gs,
                      const fieldtype & ft, const ST::string & fp,
                      const ST::string & pres, const ST::string & of, const bool & gauss, const unsigned & c)
  : FULLCOND_pspline_surf_gaussian(o,dp,fcc,intact,v1,v2,ti,nrk,degr,kp,l,gs,ft,fp,pres,of,true,c)
  {

  if(gauss==false)
    utype = iwls;

  create(v1,v2,intact);

  if(type == mrflinear)
    grenzfall = 1;
  else if(type == mrfquadratic8)
    grenzfall = 2;

  }


  // CONSTRUCTOR 7: geosplines varying coefficients

FULLCOND_pspline_surf_stepwise::FULLCOND_pspline_surf_stepwise(
                      MCMCoptions * o, DISTRIBUTION * dp, FULLCOND_const * fcc,
                      const datamatrix & intact,
                      const datamatrix & region, const MAP::map & mp, const ST::string & mn,
                      const ST::string & ti, const unsigned & nrk,
                      const unsigned & degr, const knotpos & kp, const double & l,
                      const int & gs, const fieldtype & ft, const ST::string & fp,
                      const ST::string & pres, const bool & gauss, const unsigned & c)
  : FULLCOND_pspline_surf_gaussian(o,dp,fcc,intact,region,mp,mn,ti,nrk,degr,kp,l,gs,ft,fp,pres,true,c)
  {

  if(gauss==false)
    utype = iwls;

  datamatrix v1 = datamatrix(likep->get_nrobs(),1,0.0);
  datamatrix v2 = datamatrix(likep->get_nrobs(),1,0.0);

  ST::string regname;
  for(unsigned i=0;i<likep->get_nrobs();i++)
    {
    regname = ST::doubletostring(region(i,0));
    regionnames.push_back(regname);
    v1(i,0) = m.get_region(m.getnr(regname)).get_xcenter();
    v2(i,0) = m.get_region(m.getnr(regname)).get_ycenter();
    }

  create(v1,v2,intact);

  if(type == mrflinear)
    grenzfall = 1;
  else if(type == mrfquadratic8)
    grenzfall = 2;

  }


FULLCOND_pspline_surf_stepwise::FULLCOND_pspline_surf_stepwise(const FULLCOND_pspline_surf_stepwise & fc)
  : FULLCOND_pspline_surf_gaussian(FULLCOND_pspline_surf_gaussian(fc))
  {
  mainpoi1 = fc.mainpoi1;
  mainpoi2 = fc.mainpoi2;

  interhaupt = fc.interhaupt;
  interhaupt1 = fc.interhaupt1;
  interhaupt2 = fc.interhaupt2;

  centerboth = fc.centerboth;
  maineffectsexisting = fc.maineffectsexisting;
  centerfix = fc.centerfix;
  }


const FULLCOND_pspline_surf_stepwise & FULLCOND_pspline_surf_stepwise::operator=(
                                            const FULLCOND_pspline_surf_stepwise & fc)
  {
  if (this == &fc)
    return *this;
  FULLCOND_pspline_surf_gaussian::operator=(FULLCOND_pspline_surf_gaussian(fc));
  
  mainpoi1 = fc.mainpoi1;
  mainpoi2 = fc.mainpoi2;

  interhaupt = fc.interhaupt;
  interhaupt1 = fc.interhaupt1;
  interhaupt2 = fc.interhaupt2;
  
  centerboth = fc.centerboth;
  maineffectsexisting = fc.maineffectsexisting;
  centerfix = fc.centerfix;

  return *this;
  }


bool FULLCOND_pspline_surf_stepwise::posteriormode(void)
  {
  if(maineffectsexisting != 0)
    search_maineffects();      //schlecht, weil Fkt. in jeder Iteration eines Backfitting aufgerufen wird!!!

  bool converged = false;
  bool converged1 = true;
  bool converged2 = true;
  if(!centertotal && (centerboth == 11 || centerboth == 1))
    converged1 = false;
  if(!centertotal && (centerboth == 11 || centerboth == 10))
    converged2 = false;

  transform = likep->get_trmult(column);
  fchelp.set_transform(transform);
  fctotal.set_transform(transform);

  unsigned i;

  if(utype == gaussian)
    likep->substr_linearpred_m(spline,column,true);

  compute_XWXenv(likep->get_weightiwls(),column);

  prec_env.addto(XX_env,Kenv,1.0,lambda);

  if(utype != gaussian)
    likep->substr_linearpred_m(spline,column,true);

  likep->compute_workingresiduals(column);

  compute_XWtildey(likep->get_weightiwls(),likep->get_workingresiduals(),1.0,column);

  prec_env.solve(muy,beta);

  add_linearpred_multBS(beta);

  if(center)
    {
    if(centertotal)
      {
      compute_intercept();
      for(i=0;i<nrpar;i++)
        beta(i,0) -= intercept;
      for(i=0;i<likep->get_nrobs();i++)
        spline(i,0) -= intercept;
      fcconst->posteriormode_intercept(intercept);
      intercept = 0.0;
      }
    else
      {
      beta_uncentered.assign(beta);

      intercept = 0.0;
      if(centerboth==11 || (centerboth==1 && centerfix==10) || (centerboth==10 && centerfix==1) || centerfix==11)
        compute_intercept();
      compute_main();
      //compute_beta();

      //if(centerboth==11)      // alt!!!
      //  fcconst->posteriormode_intercept(intercept);
      //else if(centerfix==11 || (centerboth==1 && centerfix==10) || (centerboth==10 && centerfix==1))
      if(centerfix==11 || (centerboth==1 && centerfix==10) || (centerboth==10 && centerfix==1) || centerboth==11)
        fcconst->update_intercept(intercept);

      if(centerboth == 11)       // neu!!!
        {
        interhaupt1 = mainpoi1->get_intercept();
        interhaupt2 = mainpoi2->get_intercept();
        }
      else
        {
        interhaupt1 = 0.0;
        interhaupt2 = 0.0;
        }                       // bis hier!!!

      if(centerboth==1 && centerfix==10)
        {
        double test = mainpoi1->get_intercept();
        if(test != 0)
          {
          interhaupt = test;
          //fcconst->update_interceptold(interhaupt);
          }
        else
          fcconst->posteriormode_intercept(interhaupt);
        }
      else if(centerboth==10 && centerfix==1)
        {
        double test = mainpoi2->get_intercept();
        if(test != 0)
          {
          interhaupt = test;
          //fcconst->update_interceptold(interhaupt);
          }
        else
          fcconst->posteriormode_intercept(interhaupt);
        }
      else
        interhaupt = 0.0;

      if(centerboth == 11 || centerboth == 1)
         converged1 = mainpoi1->changeposterior(he1,intercept);
      if(centerboth == 11 ||centerboth == 10)
         converged2 = mainpoi2->changeposterior(he2,intercept);

// Gesamteffekt in fctotal schreiben

      double * fctotalbetap = fctotal.getbetapointer();

      if(gridsize < 0)
        {
        vector<int>::iterator freqwork = freq.begin();
        int * workindex = index.getV();
        for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
          {
          if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
            {
            if(centerboth == 11)
              *fctotalbetap = spline(*workindex,0)
                            + mainpoi1->get_spline()(*workindex,0)
                            + mainpoi2->get_spline()(*workindex,0);
            else if(centerboth == 1)
              *fctotalbetap = spline(*workindex,0)
                            + mainpoi1->get_spline()(*workindex,0);
            else if(centerboth == 10)
              *fctotalbetap = spline(*workindex,0)
                            + mainpoi2->get_spline()(*workindex,0);
            fctotalbetap++;
            }
          }
        }
      else
        {
        multDG(splinehelp,beta);
        unsigned k,l;
        for(k=0;k<gridsizex;k++)
          for(l=0;l<gridsizey;l++,fctotalbetap++)
            {
            if(centerboth == 11)
              *fctotalbetap = splinehelp(k*gridsizey + l,0)
                            + mainpoi1->get_splinehelp()(k,0)
                            + mainpoi2->get_splinehelp()(l,0);
            else if(centerboth == 1)
              *fctotalbetap = splinehelp(k*gridsizey + l,0)
                            + mainpoi1->get_splinehelp()(k,0);
            else if(centerboth == 10)
              *fctotalbetap = splinehelp(k*gridsizey + l,0)
                            + mainpoi2->get_splinehelp()(l,0);
            }
        }

      fctotal.posteriormode();

      } // END: interaction
    } // END: if(center)

// Interaktionseffekt in fchelp schreiben
    double * fchelpbetap = fchelp.getbetapointer();

    if(gridsize < 0)
      {
      if(varcoeff)
        multBout(splinehelp,beta);

      vector<int>::iterator freqwork = freq.begin();
      int * workindex = index.getV();
      for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
        {
        if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
          {
          if(varcoeff)
            *fchelpbetap = splinehelp(i,0);
          else
            *fchelpbetap = spline(*workindex,0);
          fchelpbetap++;
          }
        }
      }
    else
      {
      multDG(splinehelp,beta);  // NICHT zeilen- und spaltenweise zentriert!!!
                                // dazu: if(center && !centertotal){multDG(splinehelp,beta);}
                                // und: 'splinehelp' ändern 'in compute_maineffects'
      for(i=0;i<gridsize;i++,fchelpbetap++)
        *fchelpbetap = splinehelp(i,0);
      }

  fchelp.posteriormode();
  converged = FULLCOND_nonp_basis::posteriormode();

  if(converged && converged1 && converged2)
    return true;
  else
    return false;

  }


//------------------------ FÜR INTERAKTION -------------------------------------

void FULLCOND_pspline_surf_stepwise::compute_main(void)
  {

  unsigned i;

  double *workspline;
  int *workindex = index.getV();
  vector<int>::iterator freqwork;

  datamatrix hilf = datamatrix(spline.rows(),1,0);

  // Haupteffekte berechnen
  if(maineffectsexisting == 11 || maineffectsexisting == 1)
    he1.mult(betaweightx,beta);
  if(maineffectsexisting == 11 || maineffectsexisting == 10)
    he2.mult(betaweighty,beta);

  if(centerboth == 11 || centerboth == 1 || centerfix == 11 || centerfix == 1)
    {
    // Haupteffekt berechnen
    //he1.mult(betaweightx,beta);
    // 'spline' ändern
    freqwork = mainpoi1->get_freqit();
    workindex = mainpoi1->get_indexp();
    for(i=0;i<spline.rows();i++,freqwork++,workindex++)
      spline(*workindex,0) -= he1(*freqwork,0);

    if(centerfix == 11 || centerfix == 1)   // wichtig!
      {
      freqwork = mainpoi1->get_freqit();
      workindex = mainpoi1->get_indexp();
      for(i=0;i<spline.rows();i++,freqwork++,workindex++)
         hilf(*workindex,0) = -he1(*freqwork,0) + intercept;
      likep->add_linearpred(hilf);
      }
    }
  if(centerboth == 11 || centerboth == 10 || centerfix == 11 || centerfix == 10)
    {
    // Haupteffekt berechnen
    //he2.mult(betaweighty,beta);
    // 'spline' ändern
    freqwork = mainpoi2->get_freqit();
    workindex = mainpoi2->get_indexp();
    for(i=0;i<spline.rows();i++,freqwork++,workindex++)
      spline(*workindex,0) -= he2(*freqwork,0);

    if(centerfix == 11 || centerfix == 10)
       {
       freqwork = mainpoi2->get_freqit();
       workindex = mainpoi2->get_indexp();
       for(i=0;i<spline.rows();i++,freqwork++,workindex++)
          hilf(*workindex,0) = -he2(*freqwork,0) + intercept;
       likep->add_linearpred(hilf);
       }
    }

  if(centerboth==11 || centerfix==11 || (centerboth==1 && centerfix==10) || (centerboth==10 && centerfix==1))   
    {
    workspline = spline.getV();
    for(i=0;i<spline.rows();i++,workspline++)
      *workspline += intercept;
    }
  }


void FULLCOND_pspline_surf_stepwise::reset_effect(const unsigned & pos)
  {
  remove_centering();

  likep->substr_linearpred_m(spline,column,true);

  unsigned i;
  double * work;
  work = spline.getV();
  for(i=0;i<spline.rows();i++,work++)
    *work = 0.0;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;
  intercept = 0.0;
  }


void FULLCOND_pspline_surf_stepwise::remove_centering(void)
  {
  if(center && !centertotal)
    {
    unsigned i;
    int *workindex = index.getV();
    vector<int>::iterator freqwork;

    if(centerboth == 11 || centerboth == 1)
      {
      mainpoi1->changeposterior(-he1,-intercept);
      freqwork = mainpoi1->get_freqit();
      workindex = mainpoi1->get_indexp();
      for(i=0;i<spline.rows();i++,freqwork++,workindex++)
        {
        spline(*workindex,0) += he1(*freqwork,0);
        if(*freqwork>0 && *freqwork!=*(freqwork-1))
          he1(*freqwork-1,0) = 0;
        }
      he1(*(freqwork-1),0) = 0;
      }
    if(centerboth == 11 ||centerboth == 10)
      {
      mainpoi2->changeposterior(-he2,-intercept);
      freqwork = mainpoi2->get_freqit();
      workindex = mainpoi2->get_indexp();
      for(i=0;i<spline.rows();i++,freqwork++,workindex++)
        {
        spline(*workindex,0) += he2(*freqwork,0);
        if(*freqwork>0 && *freqwork!=*(freqwork-1))
          he2(*freqwork-1,0) = 0;
        }
      he2(*(freqwork-1),0) = 0;
      }
    if(centerboth==11 || (centerboth==1 && centerfix==10) || (centerboth==10 && centerfix==1))
      {
      double addieren = -intercept;
      if(centerboth==11)
        {
        double * workspline = spline.getV();
        for(i=0;i<spline.rows();i++,workspline++)
          *workspline += addieren;
        }
      if((centerboth==1 && centerfix==10) || (centerboth==10 && centerfix==1))
        addieren -= interhaupt;
      fcconst->posteriormode_intercept(addieren);
      fcconst->update_linold();
      }
    intercept = 0.0;  
    //interhaupt = 0.0;
    }
  }


void FULLCOND_pspline_surf_stepwise::remove_centering_fix(void)
  {

  if(center && !centertotal && centerboth!=11 && centerfix!=0)
    {

    datamatrix hilf = datamatrix(spline.rows(),1,0);
    unsigned i;

    int *workindex = index.getV();
    vector<int>::iterator freqwork;

    if(centerfix == 1 || centerfix == 11)
      {
      freqwork = mainpoi1->get_freqit();
      workindex = mainpoi1->get_indexp();
      for(i=0;i<spline.rows();i++,freqwork++,workindex++)
        {
        spline(*workindex,0) += he1(*freqwork,0);
        hilf(*workindex,0) = he1(*freqwork,0) - intercept;
        }
      likep->add_linearpred(hilf);
      }
    if(centerfix == 10 || centerfix == 11)
      {
      freqwork = mainpoi2->get_freqit();
      workindex = mainpoi2->get_indexp();
      for(i=0;i<spline.rows();i++,freqwork++,workindex++)
        {
        spline(*workindex,0) += he2(*freqwork,0);
        hilf(*workindex,0) = he2(*freqwork,0) - intercept;
        }
      likep->add_linearpred(hilf);
      }

    fcconst->update_linold();

    if((centerboth==1 && centerfix==10) || (centerboth==10 && centerfix==1) || centerfix==11)
      {
      double * workspline = spline.getV();
      for(i=0;i<spline.rows();i++,workspline++)
        *workspline -= intercept;
      fcconst->posteriormode_const();
      }
    else
      intercept = 0.0;

    // Interaktionseffekt in fchelp schreiben
    double * fchelpbetap = fchelp.getbetapointer();
    if(gridsize < 0)
      {
      if(varcoeff)
        multBout(splinehelp,beta);

      vector<int>::iterator freqwork = freq.begin();
      int * workindex = index.getV();
      for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
        {
        if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
          {
          if(varcoeff)
            *fchelpbetap = splinehelp(i,0);
          else
            *fchelpbetap = spline(*workindex,0);
          fchelpbetap++;
          }
        }
      }
    else
      {
      multDG(splinehelp,beta);  // NICHT zeilen- und spaltenweise zentriert!!!
                                // dazu: if(center && !centertotal){multDG(splinehelp,beta);}
                                // und: 'splinehelp' ändern 'in compute_maineffects'
      for(i=0;i<gridsize;i++,fchelpbetap++)
        *fchelpbetap = splinehelp(i,0) - intercept;
      }
    }
  fchelp.posteriormode();
  }


void FULLCOND_pspline_surf_stepwise::get_zentrierung(FULLCOND * haupt, bool & konst)
  {
  double m2 = -interhaupt;
  if(centerboth == 11)
    {
    if(haupt == mainpoi1)
      m2 = -interhaupt1;
      //m2 = he2.mean(0);    // entspricht dem Zentrierungsfaktor des Haupteffekts bei nur einer Interaktion!
    else
      m2 = -interhaupt2;
      //m2 = he1.mean(0);
    }
  double m = -intercept+m2;

  if(haupt == mainpoi1)
    m2 = -1 * mainpoi1->get_intercept();
  else
    m2 = -1 * mainpoi2->get_intercept();  

  if(haupt == mainpoi1)
    mainpoi1->changeposterior2(he1,m);
  else if(haupt == mainpoi2)
    mainpoi2->changeposterior2(he2,m);

  if(konst == true)
    {
    fcconst->posteriormode_intercept(m2);
    }
  fcconst->update_linold();
  }


void FULLCOND_pspline_surf_stepwise::set_zentrierung(FULLCOND * haupt, int & vorzeichen, bool & inter)
  {
  unsigned i;
  double m = interhaupt;

  if(haupt == mainpoi1)
    {
    if(centerfix == 11)
      m = -he1.mean(0);
    else if(centerfix == 1)
      m = he1.mean(0);

    int *workindex = index.getV();
    vector<int>::iterator freqwork;
    freqwork = mainpoi1->get_freqit();
    workindex = mainpoi1->get_indexp();
    datamatrix hilf = datamatrix(spline.rows(),1,0);
    for(i=0;i<spline.rows();i++,freqwork++,workindex++)
      {
      if(inter == true)
        spline(*workindex,0) += vorzeichen * (he1(*freqwork,0) + m);
      hilf(*workindex,0) = vorzeichen * (he1(*freqwork,0) + m);
      }
    likep->add_linearpred(hilf);
    }
  else if(haupt == mainpoi2)
    {
    if(centerfix == 11)
      m = -he2.mean(0);
    else if(centerfix == 10)
      m = he2.mean(0);

    int *workindex = index.getV();
    vector<int>::iterator freqwork;
    freqwork = mainpoi2->get_freqit();
    workindex = mainpoi2->get_indexp();
    datamatrix hilf = datamatrix(spline.rows(),1,0);
    for(i=0;i<spline.rows();i++,freqwork++,workindex++)
      {
      if(inter == true)
        spline(*workindex,0) += vorzeichen * he2(*freqwork,0 + m);
      hilf(*workindex,0) = vorzeichen * (he2(*freqwork,0) + m);
      }
    likep->add_linearpred(hilf);
    }

  fcconst->posteriormode_const();
  }


void FULLCOND_pspline_surf_stepwise::hierarchical(ST::string & possible)
  {
  bool raus = false;
  bool fix = false;
  bool spline1, fix1;

  if(maineffectsexisting == 1 || maineffectsexisting == 11)
    {
    mainpoi1->get_inthemodel(spline1,fix1);
    if(spline1 == false && fix1 == false)
      raus = true;
    else if(fix1 == true)
      fix = true;
    }

  if(maineffectsexisting == 10 || maineffectsexisting == 11)
    {
    mainpoi2->get_inthemodel(spline1,fix1);
    if(spline1 == false && fix1 == false)
      raus = true;
    else if(fix1 == true)
      fix = true;
    }
    
  if(raus == true)
    possible = "raus";
  else if(fix == false && raus == false)
    possible = "alles";
  else
    possible = "rfix";
  }


void FULLCOND_pspline_surf_stepwise::hierarchie_rw1(vector<double> & untervector)
  {

  unsigned number = untervector.size()-1;

  update_stepwise(untervector[0]);
  double df_max = compute_df();

  update_stepwise(untervector[number]);
  double df_min = compute_df();

  if(df_max > 1 && df_min < 1)
     {
     bool geordnet = false;
     unsigned stelle_oben = number;
     unsigned stelle_unten = 0;
     while(geordnet==false)
        {
        unsigned stelle = stelle_oben + stelle_unten;
        update_stepwise(untervector[stelle/2]);
        double df_mitteunten = compute_df();
        update_stepwise(untervector[stelle/2 + 1]);
        double df_mitteoben = compute_df();

        if(df_mitteunten > 1 && df_mitteoben > 1)
          stelle_unten = stelle/2;
        else if(df_mitteunten < 1 && df_mitteoben < 1)
          stelle_oben = stelle/2 + 1;
        else
          {
          geordnet = true;
          vector<double> hilf;
          unsigned i;
          stelle_unten = stelle/2;
          stelle_oben = stelle/2 + 1;
          for(i=0;i<=stelle_unten;i++)
             hilf.push_back(untervector[i]);
          hilf.push_back(-1);
          for(i=stelle_oben;i<untervector.size();i++)
            hilf.push_back(untervector[i]);
          untervector = hilf;
          }
        }
     }
  else if(df_min >= 1)
     {
     untervector.push_back(-1);
     }
  else
     {
     vector<double> hilf;
     hilf.push_back(-1);
     unsigned i;
     for(i=0;i<untervector.size();i++)
        hilf.push_back(untervector[i]);
     untervector = hilf;
     }
  }


void FULLCOND_pspline_surf_stepwise::compute_lambdavec(
vector<double> & lvec, int & number)
  {
  if(number>0)
    {
    if (get_df_equidist()==true && number>1)
       FULLCOND::compute_lambdavec_equi(lvec,number);
    else
       FULLCOND::compute_lambdavec(lvec,number);
    }

  if((type==RW1 || type==mrflinear) && number>0)
    hierarchie_rw1(lvec);
  else  // if(varcoeff || type==RW2 || (type==RW1 && number==-1))
    lvec.push_back(-1);

  get_forced();
  if(forced_into==false)
     lvec.push_back(0);
  }


double FULLCOND_pspline_surf_stepwise::compute_df(void)
  {
  double df = FULLCOND_pspline_surf_gaussian::compute_df();

  if(maineffectsexisting != 0)
    {
    df += 1;
    double p = sqrt(df);
    df = p*p - 1;

    search_maineffects();

    if(centerboth == 1 || centerboth == 11)        // Versuch: Spur von he1 berechnen und von df abziehen!
      {
      df -= p-1;
      //datamatrix einheit = datamatrix(xv.size(),1,0);
      //double * e = einheit.getV();
      //datamatrix hilf1 = datamatrix(nrpar,1,0);
      //datamatrix hilf2 = datamatrix(xv.size(),1,0);
      //double * h2;
      //datamatrix diagonal = datamatrix(xv.size(),1,0);
      //double * d = diagonal.getV();
      //unsigned i;
      //for(i=0;i<xv.size();i++,d++,e++)
      //  {
      //  *e = 1;
      //  compute_XWtildey(likep->get_weightiwls(),einheit,1.0,column);
      //  prec_env.solve(muy,hilf1);
      //  hilf2.mult(betaweightx,hilf1);
      //  h2 = hilf2.getV()+i;
      //  *d = *h2;
      //  *e = 0;
      //  }
      //double spur = 0;
      //d = diagonal.getV();
      //for(i=0;i<xv.size();i++,d++)
      //  spur += *d;
      //df -= spur;
      }
    if(centerboth == 10 || centerboth == 11)
      {
      df -= p-1;
      //datamatrix einheit = datamatrix(yv.size(),1,0);
      //double * e = einheit.getV();
      //datamatrix hilf1 = datamatrix(nrpar,1,0);
      //datamatrix hilf2 = datamatrix(yv.size(),1,0);
      //double * h2;
      //datamatrix diagonal = datamatrix(yv.size(),1,0);
      //double * d = diagonal.getV();
      //unsigned i;
      //for(i=0;i<yv.size();i++,d++,e++)
      //  {
      //  *e = 1;
      //  compute_XWtildey(likep->get_weightiwls(),einheit,1.0,column);
      //  prec_env.solve(muy,hilf1);
      //  hilf2.mult(betaweighty,hilf1);
      //  h2 = hilf2.getV()+i;
      //  *d = *h2;
      //  *e = 0;
      //  }
      //double spur = 0;
      //d = diagonal.getV();
      //for(i=0;i<yv.size();i++,d++)
      //  spur += *d;
      //df -= spur;
      }

    if(centerfix == 1 || centerfix == 11)
      df -= 1;
    if(centerfix == 10 || centerfix == 11)
      df -= 1;
    }
    
  return df;
  }


const datamatrix & FULLCOND_pspline_surf_stepwise::get_data_forfixedeffects(void)
  {
  return data_forfixed;
  }


ST::string FULLCOND_pspline_surf_stepwise::get_effect(void)
  {
  ST::string h;

  if(varcoeff)
    {
    h = datanames[1] + "*" + datanames[0] + "(geospline";
    }
  else
    {
    if(type==mrflinear)
      h = datanames[0] + "(pspline2dimrw1";
    else
      h = datanames[0] + "(pspline2dimrw2";
    }

  h = h + ",df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";

  return h;
  }


ST::string FULLCOND_pspline_surf_stepwise::get_befehl(void)
  {
  ST::string h;

  if(varcoeff)
    {
    h = datanames[1] + "*" + datanames[0] + "(geospline";
    }
  else
    {
    if(type==mrflinear)
      h = datanames[0] + "(pspline2dimrw1";
    else
      h = datanames[0] + "(pspline2dimrw2";
    }

  h = h + ", lambda=" + ST::doubletostring(lambda,6);
  if(degree!=3)
    h = h + ", degree=" + ST::inttostring(degree);
  if(nrknots!=20)
    h = h + ",nrknots=" + ST::inttostring(nrknots);
  h = h + ")";

  return h;
  }

} // end: namespace MCMC
