
#include "baseline.h"


namespace MCMC
{


//------------------------------------------------------------------------------
//---------- CLASS: pspline_baseline (implementation of member functions) ------
//------------------------------------------------------------------------------







pspline_baseline::pspline_baseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & time, const datamatrix & z,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l,const unsigned & minb,const unsigned & maxb,
                    const fieldtype & ft,const ST::string & ti,
                    const ST::string & fp, const ST::string & pres,
                    const int & gs,const unsigned & c,const datamatrix & anfang)
 : FULLCOND_pspline(o,dp,fcc,ft,ti,nrk,degr,kp,fp,pres,false,gs,c)

{


  vc_dummy1=false;
  baselinep = vector<pspline_baseline*>(0);

  lambda = l;
  sigma2 = 1.0/l;

  zi=time;
  z_vc=z;


  if(anfang.rows()==1)
    {
     begin0 = true;
     beg_i = datamatrix(time.rows(),1,0);
    }
  else
    {
     begin0 = false;
     beg_i = anfang;
    }




//--------------bei linkstrunkierten Daten oder zeitlich var. Kovariablen Beginn der Beob.zeit berücksichtigen---------------------
   if(begin0==false)
  {

  zi_ges=datamatrix(2*zi.rows(),1,0);

  for(unsigned i=0;i<zi.rows();i++)
  {
  zi_ges(i,0)=zi(i,0);
  zi_ges(zi.rows()+i,0)=beg_i(i,0);
  }

  ges_index = statmatrix<int>(zi_ges.rows(),1);
  ges_index.indexinit();
  zi_ges.indexsort(ges_index,0,zi_ges.rows()-1,0,0);

  testmat = MCMC::bsplinemat(zi_ges,nrk,degr,kp,true);
  beg_i=anfang;
  }
//------------------------------------*/

  oldacceptance = 0;
  oldnrtrials = 0;

  min = minb;
  max = maxb;
  mintoobig = false;
  maxtoobig = false;

  unsigned i;


  varcoeff = true;
  setbeta(nrknots-1+degree,1,0);

  make_index(time,z);
  make_Bspline(time,true);
  make_BS(z);

  betaold = datamatrix(nrpar,1,0);

// xvalues und fchelp initialisieren
  ST::string pnt = fp.substr(0,fp.length()-4)+"_fchelp.raw";
  vector<int>::iterator freqwork = freqoutput.begin();
  int * workindex = index.getV();
  if(gridsize < 0)
    {
    xvalues = datamatrix(nrdiffobs,1,0);
    for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        xvalues(*freqwork,0) = time(*workindex,0);
    fchelp = FULLCOND(optionsp,datamatrix(1,1,0),title+"fchelp",nrdiffobs,1,pnt);
    splinehelp = datamatrix(likep->get_nrobs(),1,0);
    }
  else
    {
    double xmin = time.min(0);
    double xmax = time.max(0);
//    xmin=0.0;
//    xmax=2.0;
    xvalues = datamatrix(gridsize,1);
    for(i=0;i<gridsize;i++)
      xvalues(i,0) = xmin + i*(xmax-xmin)/double(xvalues.rows()-1);
    fchelp = FULLCOND(optionsp,datamatrix(1,1,0),title+"fchelp",gridsize,1,pnt);
    splinehelp = datamatrix(gridsize,1,0);
    make_DG();
    }
  fchelp.setflags(MCMC::norelchange | MCMC::nooutput);

  compute_Kweights();

  if (type == RW1)
    {
    K = Krw1(weight);
    rankK = K.get_rows()-1;
    }
  else if (type == RW2)
    {
    K = Krw2(weight);
    rankK = K.get_rows()-2;
    }

  if(minb == 0 && maxb == 0)
    {
    automatic = true;

    min = 1;
    max = rankK;
    minauto = nrpar/5;
    maxauto = nrpar/3;
    if(minauto < 1)
      minauto = 1;
    }
  else
    {
    automatic = false;

    if(max > rankK || max == 0)
      {
      maxtoobig = true;
      max = rankK;
      }
    if(min > max || min == 0)
      {
      mintoobig = true;
      min = 1;
      }
    }

  for(i=0;i<max;i++)
    {
    fc_random.push_back(datamatrix(i+1,1,0));
    randnorm.push_back(datamatrix(i+1,1,0));
    }

  make_Kab_list();

  identifiable = true;

  compute_betaweight();

//--------------------

 double knot_min = 0.0;
 double knot_max = zi.max(0);


    int_knots=datamatrix (50,1,0);

    int j;
    for(j=0;j<int_knots.rows();j++)
      {int_knots(j,0) = knot_min + j*(knot_max-knot_min)/double(int_knots.rows()-1); }



   int_D = datamatrix(int_knots.rows(),nrpar,0.0);

    datamatrix bsp;



    for(int i=0;i<int_knots.rows();i++)
      {
      bsp = bspline(int_knots(i,0));

      for(int j=0;j<nrpar;j++)
        {
        int_D(i,j) = bsp(j,0);
        }

      }



 spline_ges = datamatrix(2*likep->get_nrobs(),1,0);
 int_ti_help = datamatrix(2*likep->get_nrobs(),1,0);

 spline_zi = datamatrix(likep->get_nrobs(),1,0);
// int_ti_help = datamatrix(likep->get_nrobs(),1,0);


//-----------------------


  }



pspline_baseline::pspline_baseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & d,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l,const unsigned & minb,const unsigned & maxb,
                    const fieldtype & ft,const ST::string & ti,
                    const ST::string & fp, const ST::string & pres,
                    const int & gs, const unsigned & c,const datamatrix & anfang)
  : FULLCOND_pspline(o,dp,fcc,ft,ti,nrk,degr,kp,fp,pres,false,gs,c)
  {
  vc_dummy1=false;
  baseline = true;
  baselinep = vector<pspline_baseline*>(0);

  lambda = l;
  sigma2 = 1.0/l;

  zi=d;

  if(anfang.rows()==1)
    {
    begin0 = true;
    beg_i = datamatrix(d.rows(),1,0);
    }
  else
    {
    begin0 = false;
    beg_i = anfang;
    }




  if(begin0==false)
   {
   zi_ges=datamatrix(2*zi.rows(),1,0);

   for(unsigned i=0;i<zi.rows();i++)
    {
    zi_ges(i,0)=zi(i,0);
    zi_ges(zi.rows()+i,0)=beg_i(i,0);
    }

   ges_index = statmatrix<int>(zi_ges.rows(),1);
   ges_index.indexinit();
   zi_ges.indexsort(ges_index,0,zi_ges.rows()-1,0,0);

   testmat = MCMC::bsplinemat(zi_ges,nrk,degr,kp,true);
   beg_i=anfang;
   }




  oldacceptance = 0;
  oldnrtrials = 0;

  min = minb;
  max = maxb;
  mintoobig = false;
  maxtoobig = false;

  unsigned i;

  varcoeff = false;
  setbeta(nrknots-1+degree,1,0);

  betaold = datamatrix(nrpar,1,0);

  make_index(d);
  make_Bspline(d,true);


// xvalues und fchelp initialisieren
  ST::string pnt = fp.substr(0,fp.length()-4)+"_fchelp.raw";
  vector<int>::iterator freqwork = freq.begin();
  int * workindex = index.getV();
  if(gridsize < 0)
    {
    xvalues = datamatrix(nrdiffobs,1,0);
    for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
      if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
        xvalues(*freqwork,0) = d(*workindex,0);
    fchelp = FULLCOND(optionsp,datamatrix(1,1,0),title+"fchelp",nrdiffobs,1,pnt);
    splinehelp = datamatrix(likep->get_nrobs(),1,0);
    }
  else
    {
    double xmin = d.min(0);
    double xmax = d.max(0);
//    xmin = 0.0;
//    xmax = 2.0;
    xvalues = datamatrix(gridsize,1);
    for(i=0;i<gridsize;i++)
      xvalues(i,0) = xmin + i*(xmax-xmin)/double(xvalues.rows()-1);
    fchelp = FULLCOND(optionsp,datamatrix(1,1,0),title+"fchelp",gridsize,1,pnt);
    splinehelp = datamatrix(gridsize,1,0);
    make_DG();
    }
  fchelp.setflags(MCMC::norelchange | MCMC::nooutput);

  compute_Kweights();

  if (type == RW1)
    {
    K = Krw1(weight);
    rankK = K.get_rows()-1;
    }
  else if (type == RW2)
    {
    K = Krw2(weight);
    rankK = K.get_rows()-2;
    }

  if(minb == 0 && maxb == 0)
    {
    automatic = true;

    min = 1;
    max = rankK;
    minauto = nrpar/5;
    maxauto = nrpar/3;
    if(minauto < 1)
      minauto = 1;
    }
  else
    {
    automatic = false;

    if(max > rankK || max == 0)
      {
      maxtoobig = true;
      max = rankK;
      }
    if(min > max || min == 0)
      {
      mintoobig = true;
      min = 1;
      }
    }

  for(i=0;i<max;i++)
    {
    fc_random.push_back(datamatrix(i+1,1,0));
    randnorm.push_back(datamatrix(i+1,1,0));
    }

  make_Kab_list();

  identifiable = false;

  compute_betaweight();

//--------------------

 double knot_min = 0.0;
 double knot_max = zi.max(0);


    int_knots=datamatrix (50,1,0);

    int j;
    for(j=0;j<int_knots.rows();j++)
      {int_knots(j,0) = knot_min + j*(knot_max-knot_min)/double(int_knots.rows()-1); }



   int_D = datamatrix(int_knots.rows(),nrpar,0.0);

    datamatrix bsp;



    for(int i=0;i<int_knots.rows();i++)
      {
      bsp = bspline(int_knots(i,0));

      for(int j=0;j<nrpar;j++)
        {
        int_D(i,j) = bsp(j,0);
        }

      }



 spline_ges = datamatrix(2*likep->get_nrobs(),1,0);
 int_ti_help = datamatrix(2*likep->get_nrobs(),1,0);

 spline_zi = datamatrix(likep->get_nrobs(),1,0);
// int_ti_help = datamatrix(likep->get_nrobs(),1,0);

//------------------------------------------------------------------------


  }


pspline_baseline::pspline_baseline(const pspline_baseline & fc)
  :FULLCOND_pspline(FULLCOND_pspline(fc))
  {

  begin0 = fc.begin0;
  int_knots = fc.int_knots;
  int_D = fc.int_D;
  testmat = fc.testmat;
  zi=fc.zi;
  vc_dummy1=fc.vc_dummy1;
  beg_i=fc.beg_i;
  zi_ges=fc.zi_ges;
  z_vc=fc.z_vc;
  spline_ges=fc.spline_ges;
  spline_beg=fc.spline_beg;
  spline_zi=fc.spline_zi;
  ges_index=fc.ges_index;
  beg_index=fc.beg_index;
  int_ti_help=fc.int_ti_help;
  baselinep=fc.baselinep;

  }


const pspline_baseline & pspline_baseline::operator=(const pspline_baseline & fc)
  {
  if (this == &fc)
    return *this;
  FULLCOND_pspline::operator=(FULLCOND_pspline(fc));

  begin0 = fc.begin0;
  int_knots = fc.int_knots;
  int_D = fc.int_D;
  testmat = fc.testmat;
  zi=fc.zi;
  beg_i=fc.beg_i;
  zi_ges=fc.zi_ges;
  z_vc=fc.z_vc;
  vc_dummy1=fc.vc_dummy1;
  spline_ges=fc.spline_ges;
  spline_beg=fc.spline_beg;
  spline_zi=fc.spline_zi;
  ges_index=fc.ges_index;
  beg_index=fc.beg_index;
  int_ti_help=fc.int_ti_help;
  baselinep=fc.baselinep;

  return *this;
  }


void pspline_baseline::outoptions(void)
  {

  if(varcoeff)
    optionsp->out("  OPTIONS FOR P-SPLINE TERM: " + title + "\n",true);
  else
    optionsp->out("  OPTIONS FOR P-SPLINE TERM: " + title + " (log(baseline))\n",true);

  if(maxtoobig || mintoobig)
    optionsp->out("\n");

  if(maxtoobig)
    optionsp->out("NOTE:  Maximum blocksize is missing or too big, "
                    + ST::inttostring(max) + " has been used\n");
  if(mintoobig)
    optionsp->out("NOTE:  Minimum blocksize is missing or too big, "
                      + ST::inttostring(min) + " has been used\n");

  spline_basis::outoptions();

  if(automatic)
    {
    optionsp->out("  Initial minimum blocksize for automatic tuning: " + ST::inttostring(minauto) + "\n");
    optionsp->out("  Initial maximum blocksize for automatic tuning: " + ST::inttostring(maxauto) + "\n");
    }
  else
    {
    optionsp->out("  Minimum blocksize: " + ST::inttostring(min) + "\n");
    optionsp->out("  Maximum blocksize: " + ST::inttostring(max) + "\n");
    }

  optionsp->out("\n");

  }


/*void pspline_baseline::outresults(void)
  {

  spline_basis::outresults();

  unsigned i,j;

  ST::string pnt = pathresult.substr(0,pathresult.length()-15)+"baseline.res";
  ofstream outres(pnt.strtochar());

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  outres << "intnr" << "   ";
  if(varcoeff)
    outres << datanames[1] << "   ";
  else
    outres << datanames[0] << "   ";
  outres << "pmean   ";
  outres << "pqu"  << l1  << "   ";
  outres << "pqu"  << l2  << "   ";
  outres << "pmed   ";
  outres << "pqu"  << u1  << "   ";
  outres << "pqu"  << u2  << "   ";
  outres << "pcat" << level1 << "   ";
  outres << "pcat" << level2 << "   ";

  outres << endl;

  datamatrix sample(optionsp->get_samplesize(),1);

  double* workxvalues = xvalues.getV();
  double* workmean = fchelp.get_betameanp();
  double* wqu1l = fchelp.get_beta_lower1_p();
  double* wqu2l = fchelp.get_beta_lower2_p();
  double* wqu50 = fchelp.get_betaqu50p();
  double* wqu1u = fchelp.get_beta_upper1_p();
  double* wqu2u = fchelp.get_beta_upper2_p();

  for(i=0;i<xvalues.rows();i++,workxvalues++,workmean++,wqu1l++,wqu2l++,wqu50++,wqu1u++,wqu2u++)
    {
    fchelp.readsample(sample,i);

    for(j=0;j<sample.rows();j++)
      sample(j,0) = exp(sample(j,0));

    outres << (i+1) << "   ";
    outres << *workxvalues << "   ";
    outres << exp(*workmean) << "   ";

    outres << sample.quantile(lower1,0) << "   ";
    outres << sample.quantile(lower2,0) << "   ";
    outres << sample.quantile(50,0) << "   ";
    outres << sample.quantile(upper1,0) << "   ";
    outres << sample.quantile(upper2,0) << "   ";
// wegen plotnonp!
    outres << 1 << "   ";
    outres << 1 << "   ";

    outres << endl;

    }

  }     */







void pspline_baseline::update(void)
  {

/*/ ---------------- linearer baseline Effekt --------------------------------//
  unsigned i;
  double m;
  double var;
  double logold;
  double logprop;

  betaold(0,0) = beta(0,0);

//---------Integral berechnen----------

  m = betaold(0,0);

  for(i=0;i<zi.rows();i++)
  {
  spline(i,0)=m*zi(i,0);
  }

  compute_int_ti_linear(m);

//-------------------------------------

  logold = likep->loglikelihood(true);

  likep->substr_linearpred_m(spline,column,true);

  var = 0.0025;
  beta(0,0) = betaold(0,0) + sqrt(var)*rand_normal();

 //---------Integral berechnen----------

  m = beta(0,0);

  for(i=0;i<zi.rows();i++)
  {
  spline(i,0)=m*zi(i,0);
  }
  compute_int_ti_linear(m);

  likep->add_linearpred_m(spline,column,true);

 //-------------------------------------

  logprop = likep->loglikelihood(true);

  double u = log(uniform());

  if (u <= (logprop-logold))
    {
    acceptance++;
    }
  else
    {
    beta(0,0) = betaold(0,0);
    for(i=0;i<zi.rows();i++)
    {
    spline(i,0)=(betaold(0,0)-m)*zi(i,0);
    }

    likep->add_linearpred_m(spline,column,true);

 //---------Integral berechnen----------

      m = beta(0,0);

      for(i=0;i<zi.rows();i++)
      {
      spline(i,0)=m*zi(i,0);
      }
      compute_int_ti_linear(m);

 //-------------------------------------

    }

  if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
      (optionsp->get_nriter() % (optionsp->get_step()) == 0) )
    {

    double * fchelpbetap = fchelp.getbetapointer();

    vector<int>::iterator freqwork = freqoutput.begin();
    int * workindex = index.getV();
    for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
      {
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        {
        *fchelpbetap = spline(*workindex,0);
        fchelpbetap++;
        }
      }
    }

// ---------------- ENDE: linearer baseline Effekt ---------------------------*/


//----------------------------------------------------------------------------//

  if(baselinep.size()==2 && begin0==true)
   {
    int i=0;
    double * z_vc_help=baselinep[1]->get_z_vc();
    while((*z_vc_help==0.0||*z_vc_help==1.0)&&i<zi.rows())
     {
      i++;
      z_vc_help++;
     }
    if(i==zi.rows()) vc_dummy1=true;
   }

  unsigned blocksize;

  if(automatic)
    {
     if(optionsp->get_nriter()%100==0 && optionsp->get_nriter()<optionsp->get_burnin())
      adjust_blocksize(30,70);
      blocksize = minauto + random(maxauto-minauto+1);
    }
  else
    {
      blocksize = min + random(max-min+1);
    }

  double u;
  unsigned an = 1;
  unsigned en = blocksize;

  unsigned beg;

  double logold;
  double logprop;
  double * workbeta;
  unsigned i,j,k;



  for(j=0;j<matquant[blocksize-min];j++)
    {

    nrtrials++;

    compute_fc(beta,blocksize,an,en,sqrt(sigma2));

    logold = 0;
    logprop = 0;

    beg = firstnonzero[an-1];


//    logold += likep->loglikelihood( );
    logold += likep->loglikelihood(beg,likep->get_nrobs()-1,index);

    likep->assign(false);

    add_linearpred_multBS_Block(an-1,en-1,fc_random[en-an]);

    double * workbetaold = betaold.getV()+an-1;
    workbeta = beta.getV()+an-1;          // change beta
    for(k=an-1;k<en;k++,workbeta++,workbetaold++)
      {
      *workbetaold = *workbeta;
      *workbeta = fc_random[en-an](k-an+1,0);
      }

    //---------Integral berechnen--------

    if(baselinep.size()>0)
      {
      if(vc_dummy1==true)
      {
        vector <double *> splinevec;
        vector <double *> betavec;
        for(i=0;i<(baselinep.size());i++)
          splinevec.push_back(baselinep[i]->get_spline_zi());
        for(i=0;i<(baselinep.size());i++)
          betavec.push_back(baselinep[i]->getbetapointer());
//        multBS(spline,beta);
//        multBS(spline_zi,beta);
        compute_int_ti_vc_d0(splinevec,betavec);
        compute_int_ti_vc_d1(splinevec,betavec);
      }
      else
      {
//       testmat.mult_index(spline_ges,beta);
//        compute_int_ti_vc(beg);
       vector<double *> splinevec;
       vector<double *> betavec;
       for(i=0;i<(baselinep.size());i++)
          splinevec.push_back(baselinep[i]->get_spline_ges());
       for(i=0;i<(baselinep.size());i++)
          betavec.push_back(baselinep[i]->getbetapointer());

       compute_int_ti_vc(beg,splinevec,betavec);
      }
      }
    else
      {
      if(begin0==false)
       {
        testmat.mult(spline_ges,beta);
        compute_int_ti(beta);
       }
      else
       {
        multBS(spline,beta);
        compute_int_ti(beta);
       }
      }
    //------------------------------------

//    logprop += likep->loglikelihood(false);
    logprop += likep->loglikelihood(beg,likep->get_nrobs()-1,index,false);

    u = log(uniform());

    if (u <= (logprop-logold))
      {
	  acceptance++;
      likep->swap_linearpred();

      }
    else
      {
      workbetaold = betaold.getV()+an-1;
      workbeta = beta.getV()+an-1;
      for(k=an-1;k<en;k++,workbeta++,workbetaold++)
        *workbeta = *workbetaold;

     //---------Integral berechnen--------

    if(baselinep.size()>0)
      {
      if(vc_dummy1==true)
      {
        vector <double *> splinevec;
        vector <double *> betavec;
        for(i=0;i<(baselinep.size());i++)
          splinevec.push_back(baselinep[i]->get_spline_zi());
        for(i=0;i<(baselinep.size());i++)
          betavec.push_back(baselinep[i]->getbetapointer());
//        multBS(spline,beta);
//        multBS(spline_zi,beta);
        compute_int_ti_vc_d0(splinevec,betavec);
        compute_int_ti_vc_d1(splinevec,betavec);
      }
      else
      {
//        testmat.mult_index(spline_ges,beta);
//        compute_int_ti_vc(beg);
       vector<double *> splinevec;
       vector<double *> betavec;
       for(i=0;i<(baselinep.size());i++)
          splinevec.push_back(baselinep[i]->get_spline_ges());
       for(i=0;i<(baselinep.size());i++)
          betavec.push_back(baselinep[i]->getbetapointer());

       compute_int_ti_vc(beg,splinevec,betavec);
      }
      }
    else
      {
      if(begin0==false)
       {
        testmat.mult(spline_ges,beta);
        compute_int_ti(beta);
       }
      else
       {
        multBS(spline,beta);
        compute_int_ti(beta);
       }
      }

    //------------------------------------/

      }

    an+=blocksize;
    if (j == matquant[blocksize-min]-2)
      en = nrpar;
    else
      en+=blocksize;

    } // end: for(j=0;j<matquant[blocksize-min];j++)

  if (center)
    {
    compute_intercept();

    for(i=0;i<nrpar;i++)
      beta(i,0) -= intercept;

//-------------------Spline zentrieren------------------------

    if(baselinep.size()>0)
      {
      if(vc_dummy1==true)
       {
       for(i=0;i<likep->get_nrobs();i++)
        {
        spline(i,0) -= intercept;
        spline_zi(i,0) -= intercept;
        }
       }
      else
       {
       for(i=0;i<2.0*likep->get_nrobs();i++)
       spline_ges(i,0) -= intercept;
       }
      }
    else
      {
      if(begin0==false)
       {
       for(i=0;i<2.0*likep->get_nrobs();i++)
       spline_ges(i,0) -= intercept;

       }
      else
       {
       for(i=0;i<likep->get_nrobs();i++)
       spline(i,0) -= intercept;
       }
      }
    //-------------------------------------*/

//    likep->add_linearpred_m(-intercept,column);
    fcconst->update_intercept(intercept);
    }

  if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
      (optionsp->get_nriter() % (optionsp->get_step()) == 0) )
    {

    double * fchelpbetap = fchelp.getbetapointer();

    if(gridsize < 0)
      {
      multBS(splinehelp,beta);
      vector<int>::iterator freqwork = freqoutput.begin();
      int * workindex = index.getV();
      for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
        {
        if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
          {
          *fchelpbetap = splinehelp(i,0);
          fchelpbetap++;
          }
        }
      }
    else
      {
      multDG(splinehelp,beta);
      for(i=0;i<gridsize;i++,fchelpbetap++)
        *fchelpbetap = splinehelp(i,0);
      }
    }

  fchelp.update();
  FULLCOND::update();

  } // end function update


void pspline_baseline::compute_int_ti_mean(void)
  {
      if(baselinep.size()>0)
      {
      if(vc_dummy1==true)
      {
        vector <double *> splinevec;
        vector <double *> betavec;
        for(int i=0;i<(baselinep.size());i++)
          splinevec.push_back(baselinep[i]->get_spline_zi_mean());
        for(int i=0;i<(baselinep.size());i++)
          betavec.push_back(baselinep[i]->get_betamean());
//        multBS(spline,beta);
//        multBS(spline_zi,beta);
        compute_int_ti_vc_d0(splinevec,betavec);
        compute_int_ti_vc_d1(splinevec,betavec);
      }
      else
      {
//        testmat.mult_index(spline_ges,beta);
//        compute_int_ti_vc(beg);
       vector<double *> splinevec;
       vector<double *> betavec;
       for(int i=0;i<(baselinep.size());i++)
          splinevec.push_back(baselinep[i]->get_spline_ges_mean());
       for(int i=0;i<(baselinep.size());i++)
          betavec.push_back(baselinep[i]->get_betamean());

       compute_int_ti_vc(splinevec,betavec);
      }
      }
  else
  {
  if(begin0==false)
    {
      testmat.mult(spline_ges,betamean);
      compute_int_ti(betamean);
    }

  else
    {
    multBS(spline,betamean);
    compute_int_ti(betamean);
//    compute_int_ti_linear(betamean(0,0));  //für linearen baseline
     }
  }
  }


void pspline_baseline::compute_int_ti_linear(const double & b)
{
double * int_ti_p=likep->get_integral_ti();
for(unsigned i=0;i<zi.rows();i++,int_ti_p++)
{
if(b==0)
*int_ti_p = zi(i,0)/(exp(b*zi(i,0)));
else
*int_ti_p =(1/b*(exp(b*zi(i,0))-1.0))/(exp(b*zi(i,0)));
}

}










  void pspline_baseline::compute_int_ti(const datamatrix & b)
{

//------------------------left truncation----------------------------
 if(begin0==false)
    {

    double * int_D_help;
    double * betap;
    double dist_knots=int_knots(1,0)-int_knots(0,0);
    unsigned i,j,k;
    k=1;
    double erg,spline_u,spline_o;
    erg = 0.0;
    double * int_ti_p=likep->get_integral_ti();
    double * int_ti_help_p=int_ti_help.getV();
    double * int_ti_help_p2=int_ti_help.getV();




    spline_o=0.0;
    spline_u=0.0;
    int_D_help =int_D.getV();
    betap=b.getV();

    for(j=0;j<nrpar;j++,int_D_help++,betap++)
      spline_o += *betap* *int_D_help;
    spline_u=spline_o;

  //------------------erster Integralwert------------------------

    while(k<int_knots.rows() && int_knots(k,0)<=zi_ges(ges_index(0,0),0) )
      {

      spline_u=spline_o;
      spline_o=0.0;
      betap = b.getV();

      for(j=0;j<nrpar;j++,int_D_help++,betap++)
        spline_o += *betap* *int_D_help;

      erg=erg+(exp(spline_u)+exp(spline_o));

      k=k+1;
      }


    erg=erg*dist_knots;

    erg=erg+(exp(spline_ges(0,0))+exp(spline_o))*(zi_ges(ges_index(0,0),0)-int_knots(k-1,0));

    int_ti_p=likep->get_integral_ti()+ges_index(0,0);
    *int_ti_p =erg*0.5/(exp(spline_ges(0,0)));

    int_ti_help_p=int_ti_help.getV()+ges_index(0,0);
    *int_ti_help_p =erg*0.5;

  //------------------------------------------------------------


    for(i=1;i<zi_ges.rows();i++)
      {

      if(k==int_knots.rows())
         k=int_knots.rows()-1;

      if(k<int_knots.rows() && zi_ges(ges_index(i,0),0)<=int_knots(k,0))
         erg=erg+(zi_ges(ges_index(i,0),0)-zi_ges(ges_index(i-1,0),0))*(exp(spline_ges(i-1,0))+exp(spline_ges(i,0))) ;

      else
         {
         spline_u=spline_o;
         spline_o=0.0;
         betap = b.getV();
         for(j=0;j<nrpar;j++,int_D_help++,betap++)
            spline_o += *betap* *int_D_help;

         erg=erg+(int_knots(k,0)-zi_ges(ges_index(i-1,0),0))*(exp(spline_ges(i-1,0))+exp(spline_o)) ;

         k++;


         while(k<int_knots.rows() && int_knots(k,0)<=zi_ges(ges_index(i,0),0) )
            {
            spline_u=spline_o;
            spline_o=0.0;

            betap = b.getV();
            for(j=0;j<nrpar;j++,int_D_help++,betap++)
                spline_o += *betap* *int_D_help;

            erg=erg+dist_knots*(exp(spline_u)+exp(spline_o));
            k++;
            }

         erg=erg+(exp(spline_ges(i,0))+exp(spline_o))*(zi_ges(ges_index(i,0),0)-int_knots(k-1,0));
         }

      int_ti_p=likep->get_integral_ti()+ges_index(i,0);
      *int_ti_p =erg*0.5/(exp(spline_ges(i,0)));

      int_ti_help_p=int_ti_help.getV()+ges_index(i,0);
      *int_ti_help_p =erg*0.5;


      }

//-------------------------------------------------------------

    i=0;
    for(i=likep->get_nrobs();i<2*likep->get_nrobs();i++)
       {
       if(zi_ges(i,0)!=0)
         {
         int_ti_p=likep->get_integral_ti()+i-likep->get_nrobs();
         int_ti_help_p=int_ti_help.getV()+i-likep->get_nrobs();
         int_ti_help_p2=int_ti_help.getV()+i;

         *int_ti_p=(*int_ti_help_p-*int_ti_help_p2)/exp(spline_ges(ges_index(i-likep->get_nrobs(),0),0));
         assert(*int_ti_p>=0.0);

         }
       }
    }//left_trunc


//--------------------Beginn=0 ---------------------
 else
   {

   double * int_D_help;
   double * betap;
   double dist_knots=int_knots(1,0)-int_knots(0,0);
   unsigned i,j,k;
   k=1;
   double erg,spline_u,spline_o;
   erg = 0.0;
   double * int_ti_p=likep->get_integral_ti();
   double * int_ti_help_p=int_ti_help.getV();
//   double * int_ti_help_p2=int_ti_help.getV();

   spline_o=0.0;
   spline_u=0.0;

   int_D_help =int_D.getV();
   betap=b.getV();

   for(j=0;j<nrpar;j++,int_D_help++,betap++)
      spline_o += *betap* *int_D_help;
   spline_u=spline_o;

  //------------------erster Integralwert------------------------

   while(k<int_knots.rows() && int_knots(k,0)<=zi(index(0,0),0) )
     {

     spline_u=spline_o;
     spline_o=0.0;

     betap = b.getV();
     for(j=0;j<nrpar;j++,int_D_help++,betap++)
         spline_o += *betap* *int_D_help;

     erg=erg+(exp(spline_u)+exp(spline_o));
     k=k+1;
     }


   erg=erg*dist_knots;

   erg=erg+(exp(spline(0,0))+exp(spline_o))*(zi(index(0,0),0)-int_knots(k-1,0));

   int_ti_p=likep->get_integral_ti()+index(0,0);
   *int_ti_p =erg*0.5/(exp(spline(0,0)));

   int_ti_help_p=int_ti_help.getV()+index(0,0);
   *int_ti_help_p =erg*0.5;

  //------------------------------------------------------------


   for(i=1;i<zi.rows();i++)
      {
      if(k==int_knots.rows())
         k=int_knots.rows()-1;

      if(k<int_knots.rows() && zi(index(i,0),0)<=int_knots(k,0))
         erg=erg+(zi(index(i,0),0)-zi(index(i-1,0),0))*(exp(spline(i-1,0))+exp(spline(i,0))) ;

      else
        {
        spline_u=spline_o;
        spline_o=0.0;
        betap = b.getV();
        for(j=0;j<nrpar;j++,int_D_help++,betap++)
            spline_o += *betap* *int_D_help;
//if(k<int_knots.rows())
        erg=erg+(int_knots(k,0)-zi(index(i-1,0),0))*(exp(spline(i-1,0))+exp(spline_o)) ;
//if(k<int_knots.rows()-1)
        k++;
//cout<<k<<" "<<endl;


        while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
          {

          spline_u=spline_o;
          spline_o=0.0;

          betap = b.getV();
          for(j=0;j<nrpar;j++,int_D_help++,betap++)
              spline_o += *betap* *int_D_help;

          erg=erg+dist_knots*(exp(spline_u)+exp(spline_o));

          k++;

          }

        erg=erg+(exp(spline(i,0))+exp(spline_o))*(zi(index(i,0),0)-int_knots(k-1,0));
        }

      int_ti_p=likep->get_integral_ti()+index(i,0);
      *int_ti_p =erg*0.5/(exp(spline(i,0)));

      int_ti_help_p=int_ti_help.getV()+index(i,0);
      *int_ti_help_p =erg*0.5;
      }
  }//else, i.e. not begin0==false

} //compute_ti



void pspline_baseline::compute_int_ti_vc(vector<double *>splinevector, vector<double *>betavector)
{

//------------------------left truncation----------------------------



 double * betap;
 double dist_knots=int_knots(1,0)-int_knots(0,0);
 unsigned i,j,k,obs;
 double erg,spline_u,spline_o,spline_bi,spline_ti;

 double * int_ti_p=likep->get_integral_ti();
 statmatrix<double*> int_D_help;
 int_D_help = statmatrix<double*>(baselinep.size(),1);
 statmatrix<double*> z_vc_help;
 z_vc_help = statmatrix<double*>(baselinep.size()-1,1);
 statmatrix<double*> spline_ges_help;
 spline_ges_help = statmatrix<double*>(baselinep.size(),1);



 for(i=0;i<(baselinep.size()-1);i++)
  {
  z_vc_help(i,0) = baselinep[i+1]->get_z_vc();
  }

 for(i=0;i<(baselinep.size());i++)
  {
//  spline_ges_help(i,0) = baselinep[i]->get_spline_ges();
   spline_ges_help(i,0) = splinevector[i];

  }

 double help;


//--------------------------------------------------------------------------

 for(obs=0;obs<zi.rows();obs++,int_ti_p++)
  {


  k=0;
  erg=0.0;
  spline_u=0.0;
  spline_o=0.0;
  spline_bi=0.0;
  spline_ti=0.0;

 //-------------------Splines und Kov. der nächsten Beob. -- Anfang der Designmatrix--------
  for(i=0;i<baselinep.size();i++)
     {
     int_D_help(i,0) = baselinep[i]->get_int_D();

     if(obs>0)  spline_ges_help(i,0)++;

     if(obs>0 && i>0) z_vc_help(i-1,0)++;
     }

//-------------Spline bei ti und bi berechnen-----------------------
  for(i=0;i<baselinep.size();i++)
     {
     help= * (spline_ges_help(i,0)+zi.rows());
     if(i>0)  help = help* *z_vc_help(i-1,0);
     spline_bi +=  help;
     }

  for(i=0;i<baselinep.size();i++)
     {
     help= * (spline_ges_help(i,0));
     if(i>0)  help = help* *z_vc_help(i-1,0);
     spline_ti +=  help;
     }

//-----------ersten Knoten finden, der größer ist als bi-----------------------------

  while(k<int_knots.rows() && int_knots(k,0)<=beg_i(obs,0))
     {
     k++;
     for(i=0;i<baselinep.size();i++)
        {
        for(j=0;j<baselinep[i]->get_nrpar();j++)
           {
           int_D_help(i,0)++;
           }
        }
     }

//---------------falls kein Knoten zwischen bi und ti liegt----------------
  if(zi(obs,0)<=int_knots(k,0))
     {

     *int_ti_p=(exp(spline_bi)+exp(spline_ti))*(zi(obs,0)-beg_i(obs,0))*0.5/(exp(spline_ti));
     }
//----------------------Knoten zwischen bi und ti-------------------------
  else
     {
 //--------------Spline am ersten Knoten, der größer als bi ist (und kleiner als ti) berechenen und erstes Integralstück berechnen---------
     for(i=0;i<baselinep.size();i++)
        {
//        betap = baselinep[i]->getbetapointer();
        betap = betavector[i];
        help = 0.0;
        for(j=0;j<baselinep[i]->get_nrpar();j++,int_D_help(i,0)++,betap++)
              help += *betap* *int_D_help(i,0);
        if(i>0)  help = help* *z_vc_help(i-1,0);
          spline_o +=  help;
        }
     erg += (exp(spline_bi)+exp(spline_o))*(int_knots(k,0)-beg_i(obs,0));
     k++;

//------------------Integralstücke zwischen folgenden Knoten, die kleiner sind als ti-------------
     while(k<int_knots.rows() && int_knots(k,0)<=zi(obs,0))
        {
        spline_u=spline_o;
        spline_o=0.0;
        for(i=0;i<baselinep.size();i++)
           {
 //          betap = baselinep[i]->getbetapointer();
           betap = betavector[i];
           help = 0.0;
           for(j=0;j<baselinep[i]->get_nrpar();j++,int_D_help(i,0)++,betap++)
                   help += *betap* *int_D_help(i,0);
           if(i>0)  help = help* *z_vc_help(i-1,0);
             spline_o +=  help;
           }

        erg += (exp(spline_o)+exp(spline_u))*dist_knots;
        k++;

        }

//----------------letztes Integralstück zwischen ti und dem Knoten davor--------------
       erg += (exp(spline_ti)+exp(spline_o))*(zi(obs,0)-int_knots(k-1,0));
       *int_ti_p=erg*0.5/exp(spline_ti);

     }//else
  }//for obs
}

//-----------------------Version mit beg----------------
void pspline_baseline::compute_int_ti_vc(unsigned beg,const vector<double *>splinevector,const vector<double *>betavector)
{
double * betap;
 double dist_knots=int_knots(1,0)-int_knots(0,0);
 unsigned i,j,k,obs;
 double erg,spline_u,spline_o,spline_bi,spline_ti;

 double * int_ti_p=likep->get_integral_ti();
 statmatrix<double*> int_D_help;
 int_D_help = statmatrix<double*>(baselinep.size(),1);
 statmatrix<double*> z_vc_help;
 z_vc_help = statmatrix<double*>(baselinep.size()-1,1);
 statmatrix<double*> spline_ges_help;
 spline_ges_help = statmatrix<double*>(baselinep.size(),1);

 double beg_help = zi(index(beg,0),0);


 for(i=0;i<(baselinep.size()-1);i++)
  {
  z_vc_help(i,0) = baselinep[i+1]->get_z_vc();
  }

 for(i=0;i<(baselinep.size());i++)
  {
 // spline_ges_help(i,0) = baselinep[i]->get_spline_ges();
  spline_ges_help(i,0) = splinevector[i];
  }

 double help;


//--------------------------------------------------------------------------

 for(obs=0;obs<zi.rows();obs++,int_ti_p++)
  {


  k=0;
  erg=0.0;
  spline_u=0.0;
  spline_o=0.0;
  spline_bi=0.0;
  spline_ti=0.0;

 //-------------------Splines und Kov. der nächsten Beob. -- Anfang der Designmatrix--------
  for(i=0;i<baselinep.size();i++)
     {
     int_D_help(i,0) = baselinep[i]->get_int_D();

     if(obs>0)  spline_ges_help(i,0)++;

     if(obs>0 && i>0) z_vc_help(i-1,0)++;
     }

//--------------------muss int_ti neu berechnet werden?-----------------------------------------------------
  if(zi(obs,0)>=beg_help)
  {
//-------------Spline bei ti und bi berechnen-----------------------
  for(i=0;i<baselinep.size();i++)
     {
     help= * (spline_ges_help(i,0)+zi.rows());
     if(i>0)  help = help* *z_vc_help(i-1,0);
     spline_bi +=  help;
     }

  for(i=0;i<baselinep.size();i++)
     {
     help= * (spline_ges_help(i,0));
     if(i>0)  help = help* *z_vc_help(i-1,0);
     spline_ti +=  help;
     }

//-----------ersten Knoten finden, der größer ist als bi-----------------------------

  while(k<int_knots.rows() && int_knots(k,0)<=beg_i(obs,0))
     {
     k++;
     for(i=0;i<baselinep.size();i++)
        {
        for(j=0;j<baselinep[i]->get_nrpar();j++)
           {
           int_D_help(i,0)++;
           }
        }
     }

//---------------falls kein Knoten zwischen bi und ti liegt----------------
  if(zi(obs,0)<=int_knots(k,0))
     {

     *int_ti_p=(exp(spline_bi)+exp(spline_ti))*(zi(obs,0)-beg_i(obs,0))*0.5/(exp(spline_ti));
     }
//----------------------Knoten zwischen bi und ti-------------------------
  else
     {
 //--------------Spline am ersten Knoten, der größer als bi ist (und kleiner als ti) berechenen und erstes Integralstück berechnen---------
     for(i=0;i<baselinep.size();i++)
        {
 //       betap = baselinep[i]->getbetapointer();
        betap = betavector[i];
        help = 0.0;
        for(j=0;j<baselinep[i]->get_nrpar();j++,int_D_help(i,0)++,betap++)
              help += *betap* *int_D_help(i,0);
        if(i>0)  help = help* *z_vc_help(i-1,0);
          spline_o +=  help;
        }
     erg += (exp(spline_bi)+exp(spline_o))*(int_knots(k,0)-beg_i(obs,0));
     k++;

//------------------Integralstücke zwischen folgenden Knoten, die kleiner sind als ti-------------
     while(k<int_knots.rows() && int_knots(k,0)<=zi(obs,0))
        {
        spline_u=spline_o;
        spline_o=0.0;
        for(i=0;i<baselinep.size();i++)
           {
//           betap = baselinep[i]->getbetapointer();
           betap = betavector[i];
           help = 0.0;
           for(j=0;j<baselinep[i]->get_nrpar();j++,int_D_help(i,0)++,betap++)
                   help += *betap* *int_D_help(i,0);
           if(i>0)  help = help* *z_vc_help(i-1,0);
             spline_o +=  help;
           }

        erg += (exp(spline_o)+exp(spline_u))*dist_knots;
        k++;

        }

//----------------letztes Integralstück zwischen ti und dem Knoten davor--------------
       erg += (exp(spline_ti)+exp(spline_o))*(zi(obs,0)-int_knots(k-1,0));
       *int_ti_p=erg*0.5/exp(spline_ti);

     }//else
  }//for obs
  }//if
}



void pspline_baseline::compute_int_ti_vc_d0(const vector<double *> splinevector,const vector<double *> betavector)
{

 double * int_D_help;
 double * betap;
 double dist_knots=int_knots(1,0)-int_knots(0,0);
 unsigned i,j,k,i_help;
 i=0;
 i_help=0;
 double erg,spline_u,spline_o,spline_ti_help,spline_ti;
 double * int_ti_p=likep->get_integral_ti();
 double * int_ti_help_p=int_ti_help.getV();
 double * z_vc_help;
 double * spline_help;
 z_vc_help = baselinep[1]->get_z_vc();
 k=1;
 erg=0.0;
 spline_o=0.0;
 spline_u=0.0;
// spline_help=baselinep[0]->get_spline_zi();
 spline_help = splinevector[0];

 int_D_help =baselinep[0]->get_int_D();
// betap=baselinep[0]->getbetapointer();
 betap = betavector[0];

 for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help++,betap++)
    spline_o += *betap* *int_D_help;
 spline_u=spline_o;
  //------i0 berechnen---------
 while(*(z_vc_help+index(i,0))!=0.0)
   {
   i++;
   spline_help++;
   }

 i_help=i;

  //------------------erster Integralwert------------------------

 while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
   {
   spline_u=spline_o;
   spline_o=0.0;

//   betap = baselinep[0]->getbetapointer();
   betap = betavector[0];
   for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help++,betap++)
      spline_o += *betap* *int_D_help;

   erg=erg+(exp(spline_u)+exp(spline_o));
   k=k+1;
   }


 erg=erg*dist_knots;
 spline_ti=*spline_help;
 spline_ti_help=spline_ti;
 erg=erg+(exp(spline_ti)+exp(spline_o))*(zi(index(i,0),0)-int_knots(k-1,0));

 int_ti_p=likep->get_integral_ti()+index(i,0);
 *int_ti_p =erg*0.5/(exp(spline_ti));

 int_ti_help_p=int_ti_help.getV()+index(i,0);
 *int_ti_help_p =erg*0.5;

  //------------------------------------------------------------

 i++;
 spline_help++;
 while(i<zi.rows())
   {
   if(*(z_vc_help+index(i,0))!=0)
     {
     i++;
     spline_help++;
     }
   else
     {
     if(k== int_knots.rows())
         k=int_knots.rows()-1;
     spline_ti=*spline_help;
     if(k<int_knots.rows() && zi(index(i,0),0)<=int_knots(k,0))
         erg=erg+(zi(index(i,0),0)-zi(index(i_help,0),0))*(exp(spline_ti_help)+exp(spline_ti)) ;
     else
        {
        spline_u=spline_o;
        spline_o=0.0;
//        betap = baselinep[0]->getbetapointer();
        betap = betavector[0];
        for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help++,betap++)
              spline_o += *betap* *int_D_help;
        erg=erg+(int_knots(k,0)-zi(index(i_help,0),0))*(exp(spline_ti_help)+exp(spline_o)) ;
        k++;

        while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
           {
           spline_u=spline_o;
           spline_o=0.0;

//           betap = baselinep[0]->getbetapointer();
           betap = betavector[0];
           for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help++,betap++)
                  spline_o += *betap* *int_D_help;
           erg=erg+dist_knots*(exp(spline_u)+exp(spline_o));
           k++;
           }

        erg=erg+(exp(spline_ti)+exp(spline_o))*(zi(index(i,0),0)-int_knots(k-1,0));
        }

     int_ti_p=likep->get_integral_ti()+index(i,0);
     *int_ti_p =erg*0.5/(exp(spline_ti));

     int_ti_help_p=int_ti_help.getV()+index(i,0);
     *int_ti_help_p =erg*0.5;

     i_help=i;
     i++;
     spline_help++;
     spline_ti_help=spline_ti;
     }//--------------- else: d.h. z_vc(index(i,0),0)==0 -----------------------
 } //while
}



void pspline_baseline::compute_int_ti_vc_d1(const vector<double *> splinevector,const vector<double *> betavector)
{

 double * betap;
 double dist_knots=int_knots(1,0)-int_knots(0,0);
 unsigned i,j,k,i_help;
 i=0;
 i_help=0;
 double erg,spline_u,spline_o,spline_ti_help,spline_ti;
 double * int_ti_p=likep->get_integral_ti();
 double * int_ti_help_p=int_ti_help.getV();
 double * z_vc_help;


 spline_ti=0.0;
 spline_ti_help=0.0;

 statmatrix<double*> int_D_help_1;
 int_D_help_1=statmatrix<double*>(baselinep.size(),1);

 statmatrix<double*> spline_zi_help;
 spline_zi_help=statmatrix<double*>(baselinep.size(),1);



 double help;
 int i_vc;
 k=1;
 erg=0.0;
 spline_o=0.0;
 spline_u=0.0;
 z_vc_help = baselinep[1]->get_z_vc();

 for(i_vc=0;i_vc<(baselinep.size());i_vc++)
   {
 //  spline_zi_help(i_vc,0)=baselinep[i_vc]->get_spline_zi();
   spline_zi_help(i_vc,0) = splinevector[i_vc];
   int_D_help_1(i_vc,0)= baselinep[i_vc]->get_int_D();
   }

 for(i_vc=0;i_vc<baselinep.size();i_vc++)
   {
//   betap = baselinep[i_vc]->getbetapointer();
   betap = betavector[i_vc];
   help=0.0;
   for(j=0;j<baselinep[i_vc]->get_nrpar();j++,int_D_help_1(i_vc,0)++,betap++)
         help+=*betap* *int_D_help_1(i_vc,0);
   spline_o +=help;
   }

 spline_u=spline_o;
  //------kleinstes ti mit z_vc==1 und splines an der Stelle suchen---------
 while(*(z_vc_help+index(i,0))==0.0)
   {
   i++;
   for(i_vc=0;i_vc<(baselinep.size());i_vc++)
      {
      spline_zi_help(i_vc,0)++;
      }
   }

 i_help=i;
//------------------erster Integralwert------------------------

 while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
   {
    spline_u=spline_o;
    spline_o=0.0;

    for(i_vc=0;i_vc<baselinep.size();i_vc++)
      {
//      betap = baselinep[i_vc]->getbetapointer();
      betap = betavector[i_vc];
      help=0.0;
      for(j=0;j<baselinep[i_vc]->get_nrpar();j++,int_D_help_1(i_vc,0)++,betap++)
          help+=*betap* *int_D_help_1(i_vc,0);
           //   if(i_vc>0) help= help * *z_vc_help;
      spline_o +=help;
      }


    erg=erg+(exp(spline_u)+exp(spline_o));
    k=k+1;
   }


 erg=erg*dist_knots;

 //---------------Spline an der Stelle ti ---------------

 for(i_vc=0;i_vc<baselinep.size();i_vc++)
   {
   help= *spline_zi_help(i_vc,0);
    //    if(i_vc>0) help= help * *z_vc(index(i1,0),0);
   spline_ti +=help;
   }

 spline_ti_help=spline_ti;
 erg=erg+(exp(spline_ti)+exp(spline_o))*(zi(index(i,0),0)-int_knots(k-1,0));

 int_ti_p=likep->get_integral_ti()+index(i,0);
 *int_ti_p =erg*0.5/(exp(spline_ti));

 int_ti_help_p=int_ti_help.getV()+index(i,0);
 *int_ti_help_p =erg*0.5;

//------------------------------------------------------------

 i++;

 for(i_vc=0;i_vc<(baselinep.size());i_vc++)
   {
   spline_zi_help(i_vc,0)++;
   }

 while(i<zi.rows())
   {
//-----------falls z_vc==1: zu nächster Beobachtung---------------
   if(*(z_vc_help+index(i,0))==0)
      {
      i++;
      for(i_vc=0;i_vc<(baselinep.size());i_vc++)
         {
         spline_zi_help(i_vc,0)++;
         }
      }
//-------------falls z_vc==0: Integral berechenen---------------------
   else
      {
      if(k== int_knots.rows())
          k=int_knots.rows()-1;

         //Spline an der Stelle ti
      spline_ti=0.0;
      for(i_vc=0;i_vc<baselinep.size();i_vc++)
          {
          help= *spline_zi_help(i_vc,0);
             //    if(i_vc>0) help= help * z_vc(index(i,0),0);
          spline_ti +=help;
          }

      if(k<int_knots.rows() && zi(index(i,0),0)<=int_knots(k,0))
      erg=erg+(zi(index(i,0),0)-zi(index(i_help,0),0))*(exp(spline_ti_help)+exp(spline_ti)) ;

      else
          {
          spline_u=spline_o;
          spline_o=0.0;

          for(i_vc=0;i_vc<baselinep.size();i_vc++)
             {
//             betap = baselinep[i_vc]->getbetapointer();
             betap = betavector[i_vc];
             help=0.0;
             for(j=0;j<baselinep[i_vc]->get_nrpar();j++,int_D_help_1(i_vc,0)++,betap++)
                    help+=*betap* *int_D_help_1(i_vc,0);
               //    if(i_vc>0) help= help * z_vc(index(i,0),0);

             spline_o +=help;
             }

          erg=erg+(int_knots(k,0)-zi(index(i_help,0),0))*(exp(spline_ti_help)+exp(spline_o)) ;
          k++;

          while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
             {
             spline_u=spline_o;
             spline_o=0.0;

             for(i_vc=0;i_vc<baselinep.size();i_vc++)
                 {
//                 betap = baselinep[i_vc]->getbetapointer();
                 betap = betavector[i_vc];
                 help=0.0;
                 for(j=0;j<baselinep[i_vc]->get_nrpar();j++,int_D_help_1(i_vc,0)++,betap++)
                           help+=*betap* *int_D_help_1(i_vc,0);
                      //    if(i_vc>0) help= help * z_vc(index(i,0),0);

                 spline_o +=help;
                 }


             erg=erg+dist_knots*(exp(spline_u)+exp(spline_o));
             k++;
             }

          erg=erg+(exp(spline_ti)+exp(spline_o))*(zi(index(i,0),0)-int_knots(k-1,0));
          }

      int_ti_p=likep->get_integral_ti()+index(i,0);
      *int_ti_p =erg*0.5/(exp(spline_ti));

      int_ti_help_p=int_ti_help.getV()+index(i,0);
      *int_ti_help_p =erg*0.5;

      i_help=i;
      spline_ti_help=spline_ti;
      i++;
         //  z_vc_help++;
      for(i_vc=0;i_vc<(baselinep.size());i_vc++)
          {
          spline_zi_help(i_vc,0)++;
          }


      }//--------------- else: d.h. z_vc(index(i,0),0)==1 -----------------------
  } //while
}//--------------compute_int_ti_vc_d1---------------------




void pspline_baseline::compute_int_ti(unsigned beg)
{
 double * int_D_help;
 double * betap;
 double dist_knots=int_knots(1,0)-int_knots(0,0);
 unsigned i,j,k;
 k=1;
 double erg,spline_u,spline_o;
 erg = 0.0;
 double * int_ti_p=likep->get_integral_ti();
 double * int_ti_help_p=int_ti_help.getV();
 spline_o=0.0;
 spline_u=0.0;
 int_D_help =int_D.getV();
 betap=beta.getV();
 if(beg==0)
  {
  for(j=0;j<nrpar;j++,int_D_help++,betap++)
         spline_o += *betap* *int_D_help;
  spline_u=spline_o;

 //------------------erster Integralwert------------------------

 while(k<int_knots.rows() && int_knots(k,0)<=zi(index(0,0),0) )
    {
    spline_u=spline_o;
    spline_o=0.0;

    betap = beta.getV();
    for(j=0;j<nrpar;j++,int_D_help++,betap++)
              spline_o += *betap* *int_D_help;
    erg=erg+(exp(spline_u)+exp(spline_o));
    k=k+1;
    }


 erg=erg*dist_knots;

 erg=erg+(exp(spline(0,0))+exp(spline_o))*(zi(index(0,0),0)-int_knots(k-1,0));

 int_ti_p=likep->get_integral_ti()+index(0,0);
 *int_ti_p =erg*0.5/(exp(spline(0,0)));

 int_ti_help_p=int_ti_help.getV()+index(0,0);
 *int_ti_help_p =erg*0.5;

 //------------------------------------------------------------


 for(i=1;i<zi.rows();i++)
     {
     if(k==int_knots.rows())
        k=int_knots.rows()-1;
    

     if(k<int_knots.rows() && zi(index(i,0),0)<=int_knots(k,0))
            erg=erg+(zi(index(i,0),0)-zi(index(i-1,0),0))*(exp(spline(i-1,0))+exp(spline(i,0))) ;
     else
        {
        spline_u=spline_o;
        spline_o=0.0;
        betap = beta.getV();
        for(j=0;j<nrpar;j++,int_D_help++,betap++)
                spline_o += *betap* *int_D_help;

        erg=erg+(int_knots(k,0)-zi(index(i-1,0),0))*(exp(spline(i-1,0))+exp(spline_o)) ;
        k++;

        while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
           {
           spline_u=spline_o;
           spline_o=0.0;

           betap = beta.getV();
           for(j=0;j<nrpar;j++,int_D_help++,betap++)
                 spline_o += *betap* *int_D_help;

           erg=erg+dist_knots*(exp(spline_u)+exp(spline_o));
           k++;
           }

        erg=erg+(exp(spline(i,0))+exp(spline_o))*(zi(index(i,0),0)-int_knots(k-1,0));
        }

     int_ti_p=likep->get_integral_ti()+index(i,0);
     *int_ti_p =erg*0.5/(exp(spline(i,0)));

     int_ti_help_p=int_ti_help.getV()+index(i,0);
     *int_ti_help_p =erg*0.5;
     }
   }//----------------if beg==0-------------------------

    //-----------------if beg!=0-----------------------------
 else
   {
//---------------------erg=Integral bei "beg-1"------------------
   int_ti_help_p=int_ti_help.getV()+index(beg-1,0);
   erg=*int_ti_help_p*2.0;
     //ersten Knoten finden, der größer ist als "beg-1"
   while(k<int_knots.rows() && int_knots(k,0)<=zi(index(beg-1,0),0) )
      {
      for(j=0;j<nrpar;j++)
           int_D_help++;

      k++;
      }
//--------Wert des Splines an diesem Knoten ausrechnen----------
   for(j=0;j<nrpar;j++,int_D_help++,betap++)
         spline_o += *betap* *int_D_help;
   spline_u=spline_o;


//-----------------------------------------------------------------------

   for(i=beg;i<zi.rows();i++)
      {
      if(k==int_knots.rows())
              k=int_knots.rows()-1;
      if(k<int_knots.rows() && zi(index(i,0),0)<=int_knots(k,0))
              erg=erg+(zi(index(i,0),0)-zi(index(i-1,0),0))*(exp(spline(i-1,0))+exp(spline(i,0))) ;
      else
         {
         spline_u=spline_o;
         spline_o=0.0;
         betap = beta.getV();
         for(j=0;j<nrpar;j++,int_D_help++,betap++)
                spline_o += *betap* *int_D_help;

         erg=erg+(int_knots(k,0)-zi(index(i-1,0),0))*(exp(spline(i-1,0))+exp(spline_o)) ;
         k++;

         while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
              {
              spline_u=spline_o;
              spline_o=0.0;

              betap = beta.getV();
              for(j=0;j<nrpar;j++,int_D_help++,betap++)
                     spline_o += *betap* *int_D_help;

              erg=erg+dist_knots*(exp(spline_u)+exp(spline_o));
              k++;
              }

         erg=erg+(exp(spline(i,0))+exp(spline_o))*(zi(index(i,0),0)-int_knots(k-1,0));
         }

      int_ti_p=likep->get_integral_ti()+index(i,0);
      *int_ti_p =erg*0.5/(exp(spline(i,0)));

      int_ti_help_p=int_ti_help.getV()+index(i,0);
      *int_ti_help_p =erg*0.5;
      }
  }// ----------else, d.h. beg!=0
} //compute_ti



void pspline_baseline::compute_int_ti(void)
{

//------------------------left truncation----------------------------
 if(begin0==false)
    {

    double * int_D_help;
    double * betap;
    double dist_knots=int_knots(1,0)-int_knots(0,0);
    unsigned i,j,k;
    k=1;
    double erg,spline_u,spline_o;
    erg = 0.0;
    double * int_ti_p=likep->get_integral_ti();
    double * int_ti_help_p=int_ti_help.getV();
    double * int_ti_help_p2=int_ti_help.getV();




    spline_o=0.0;
    spline_u=0.0;
    int_D_help =int_D.getV();
    betap=beta.getV();

    for(j=0;j<nrpar;j++,int_D_help++,betap++)
      spline_o += *betap* *int_D_help;
    spline_u=spline_o;

  //------------------erster Integralwert------------------------

    while(k<int_knots.rows() && int_knots(k,0)<=zi_ges(ges_index(0,0),0) )
      {

      spline_u=spline_o;
      spline_o=0.0;
      betap = beta.getV();

      for(j=0;j<nrpar;j++,int_D_help++,betap++)
        spline_o += *betap* *int_D_help;

      erg=erg+(exp(spline_u)+exp(spline_o));

      k=k+1;
      }


    erg=erg*dist_knots;

    erg=erg+(exp(spline_ges(0,0))+exp(spline_o))*(zi_ges(ges_index(0,0),0)-int_knots(k-1,0));

    int_ti_p=likep->get_integral_ti()+ges_index(0,0);
    *int_ti_p =erg*0.5/(exp(spline_ges(0,0)));

    int_ti_help_p=int_ti_help.getV()+ges_index(0,0);
    *int_ti_help_p =erg*0.5;

  //------------------------------------------------------------


    for(i=1;i<zi_ges.rows();i++)
      {

      if(k==int_knots.rows())
         k=int_knots.rows()-1;

      if(k<int_knots.rows() && zi_ges(ges_index(i,0),0)<=int_knots(k,0))
         erg=erg+(zi_ges(ges_index(i,0),0)-zi_ges(ges_index(i-1,0),0))*(exp(spline_ges(i-1,0))+exp(spline_ges(i,0))) ;

      else
         {
         spline_u=spline_o;
         spline_o=0.0;
         betap = beta.getV();
         for(j=0;j<nrpar;j++,int_D_help++,betap++)
            spline_o += *betap* *int_D_help;

         erg=erg+(int_knots(k,0)-zi_ges(ges_index(i-1,0),0))*(exp(spline_ges(i-1,0))+exp(spline_o)) ;

         k++;


         while(k<int_knots.rows() && int_knots(k,0)<=zi_ges(ges_index(i,0),0) )
            {
            spline_u=spline_o;
            spline_o=0.0;

            betap = beta.getV();
            for(j=0;j<nrpar;j++,int_D_help++,betap++)
                spline_o += *betap* *int_D_help;

            erg=erg+dist_knots*(exp(spline_u)+exp(spline_o));
            k++;
            }

         erg=erg+(exp(spline_ges(i,0))+exp(spline_o))*(zi_ges(ges_index(i,0),0)-int_knots(k-1,0));
         }

      int_ti_p=likep->get_integral_ti()+ges_index(i,0);
      *int_ti_p =erg*0.5/(exp(spline_ges(i,0)));

      int_ti_help_p=int_ti_help.getV()+ges_index(i,0);
      *int_ti_help_p =erg*0.5;


      }

//-------------------------------------------------------------

    i=0;
    for(i=likep->get_nrobs();i<2*likep->get_nrobs();i++)
       {
       if(zi_ges(i,0)!=0)
         {
         int_ti_p=likep->get_integral_ti()+i-likep->get_nrobs();
         int_ti_help_p=int_ti_help.getV()+i-likep->get_nrobs();
         int_ti_help_p2=int_ti_help.getV()+i;

         *int_ti_p=(*int_ti_help_p-*int_ti_help_p2)/exp(spline_ges(ges_index(i-likep->get_nrobs(),0),0));
         assert(*int_ti_p>=0.0);

         }
       }
    }//left_trunc


//--------------------Beginn=0 ---------------------
 else
   {

   double * int_D_help;
   double * betap;
   double dist_knots=int_knots(1,0)-int_knots(0,0);
   unsigned i,j,k;
   k=1;
   double erg,spline_u,spline_o;
   erg = 0.0;
   double * int_ti_p=likep->get_integral_ti();
   double * int_ti_help_p=int_ti_help.getV();
 //  double * int_ti_help_p2=int_ti_help.getV();

   spline_o=0.0;
   spline_u=0.0;

   int_D_help =int_D.getV();
   betap=beta.getV();

   for(j=0;j<nrpar;j++,int_D_help++,betap++)
      spline_o += *betap* *int_D_help;
   spline_u=spline_o;

  //------------------erster Integralwert------------------------

   while(k<int_knots.rows() && int_knots(k,0)<=zi(index(0,0),0) )
     {

     spline_u=spline_o;
     spline_o=0.0;

     betap = beta.getV();
     for(j=0;j<nrpar;j++,int_D_help++,betap++)
         spline_o += *betap* *int_D_help;

     erg=erg+(exp(spline_u)+exp(spline_o));
     k=k+1;
     }


   erg=erg*dist_knots;

   erg=erg+(exp(spline(0,0))+exp(spline_o))*(zi(index(0,0),0)-int_knots(k-1,0));

   int_ti_p=likep->get_integral_ti()+index(0,0);
   *int_ti_p =erg*0.5/(exp(spline(0,0)));

   int_ti_help_p=int_ti_help.getV()+index(0,0);
   *int_ti_help_p =erg*0.5;

  //------------------------------------------------------------


   for(i=1;i<zi.rows();i++)
      {
      if(k==int_knots.rows())
         k=int_knots.rows()-1;

      if(k<int_knots.rows() && zi(index(i,0),0)<=int_knots(k,0))
         erg=erg+(zi(index(i,0),0)-zi(index(i-1,0),0))*(exp(spline(i-1,0))+exp(spline(i,0))) ;

      else
        {
        spline_u=spline_o;
        spline_o=0.0;
        betap = beta.getV();
        for(j=0;j<nrpar;j++,int_D_help++,betap++)
            spline_o += *betap* *int_D_help;
//if(k<int_knots.rows())
        erg=erg+(int_knots(k,0)-zi(index(i-1,0),0))*(exp(spline(i-1,0))+exp(spline_o)) ;
//if(k<int_knots.rows()-1)
        k++;
//cout<<k<<" "<<endl;


        while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
          {

          spline_u=spline_o;
          spline_o=0.0;

          betap = beta.getV();
          for(j=0;j<nrpar;j++,int_D_help++,betap++)
              spline_o += *betap* *int_D_help;

          erg=erg+dist_knots*(exp(spline_u)+exp(spline_o));

          k++;

          }

        erg=erg+(exp(spline(i,0))+exp(spline_o))*(zi(index(i,0),0)-int_knots(k-1,0));
        }

      int_ti_p=likep->get_integral_ti()+index(i,0);
      *int_ti_p =erg*0.5/(exp(spline(i,0)));

      int_ti_help_p=int_ti_help.getV()+index(i,0);
      *int_ti_help_p =erg*0.5;
      }
  }//else, i.e. not begin0==false

} //compute_ti







} // end: namespace MCMC



