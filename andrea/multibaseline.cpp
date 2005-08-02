
#include "multibaseline.h"


namespace MCMC
{


//------------------------------------------------------------------------------
//-----------------------Multi-Baseline-----------------------------------------
//------------------------------------------------------------------------------
pspline_multibaseline::pspline_multibaseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & d, const double & a, const double & b,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l,const unsigned & minb,const unsigned & maxb,
                    const fieldtype & ft,const ST::string & ti,
                    const ST::string & fp, const ST::string & pres,
                    const int & gs, const unsigned & c,const datamatrix & zustand, const datamatrix & anfang, const bool & wb)
  : FULLCOND_pspline(o,dp,fcc,ft,ti,nrk,degr,kp,fp,pres,false,gs,c)
  {
  unsigned i,j;

  baseline = true;
  baselinep = vector<pspline_multibaseline*>(0);

  lambda = l;
  sigma2 = 1.0/l;

  zi = d;
  col =c;

//  Weibull = wb;

  state_i = zustand;

  if(anfang.rows()==1)
    {
    begin0 = true;
    beg_i = datamatrix(zi.rows(),1,0);
    }
  else
    {
    begin0 = false;
    beg_i = anfang;
    }

  zi_ges = datamatrix(2*zi.rows(),1,0);

  for(i=0;i<zi.rows();i++)
    {
    zi_ges(i,0) = zi(i,0);
    zi_ges(zi.rows()+i,0) = beg_i(i,0);
    }
  ges_index = statmatrix<int>(zi_ges.rows(),1);
  ges_index.indexinit();
  zi_ges.indexsort(ges_index,0,zi_ges.rows()-1,0,0);

  testmat = MCMC::bsplinemat(zi_ges,nrk,degr,kp,true);



  double maxzi=0.0;
  for(i=0;i<zi.rows();i++)
    if (zi(i,0)>maxzi) maxzi=zi(i,0);

//-----------------------------------------------------------------------------



  oldacceptance = 0;
  oldnrtrials = 0;

  min = minb;
  max = maxb;
  mintoobig = false;
  maxtoobig = false;

  varcoeff = false;
  setbeta(nrknots-1+degree,1,0);
  betaold = datamatrix(nrpar,1,0);

  make_index(d);
  make_index2();
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

//------------------Designmatrix int_D für P-Spline an Knoten-------------------

  double knot_min = 0.0;
  double knot_max = zi.max(0);
  int_knots=datamatrix (50,1,0);
  for(j=0;j<int_knots.rows();j++)
    int_knots(j,0) = knot_min + j*(knot_max-knot_min)/double(int_knots.rows()-1);

  int_D = datamatrix(int_knots.rows(),nrpar,0.0);
  datamatrix bsp;
  for(i=0;i<int_knots.rows();i++)
    {
    bsp = bspline(int_knots(i,0));
    for(j=0;j<nrpar;j++)
      {
      int_D(i,j) = bsp(j,0);
      }
    }
//------------------------------------------------------------------------------

  spline_ges = datamatrix(2*likep->get_nrobs(),1,0);
  spline_ges2 = datamatrix(2*likep->get_nrobs(),1,0);
  int_ti_help = datamatrix(2*likep->get_nrobs(),1,0);
  spline_zi = datamatrix(likep->get_nrobs(),1,0);
//------------------------------------------------------------------------

  double sum_logti = 0.0;
  for(i=0;i<zi.rows();i++)
    {
    sum_logti = sum_logti + log(zi(i,0))*likep->get_response(i,0);
    }

  }


pspline_multibaseline::pspline_multibaseline(const pspline_multibaseline & fc)
  :FULLCOND_pspline(FULLCOND_pspline(fc))
  {

  begin0 = fc.begin0;
  int_knots = fc.int_knots;
  int_D = fc.int_D;
  testmat = fc.testmat;
  zi = fc.zi;
  vc_dummy1 = fc.vc_dummy1;
  beg_i = fc.beg_i;
  state_i = fc.state_i;
  zi_ges = fc.zi_ges;
  z_vc = fc.z_vc;
  spline_ges = fc.spline_ges;
  spline_ges2 = fc.spline_ges2;
  spline_zi = fc.spline_zi;
  ges_index = fc.ges_index;
  int_ti_help = fc.int_ti_help;
  baselinep = fc.baselinep;
  }


const pspline_multibaseline & pspline_multibaseline::operator=(const pspline_multibaseline & fc)
  {
  if (this == &fc)
    return *this;
  FULLCOND_pspline::operator=(FULLCOND_pspline(fc));

  begin0 = fc.begin0;
  int_knots = fc.int_knots;
  int_D = fc.int_D;
  testmat = fc.testmat;
  zi=fc.zi;
  beg_i = fc.beg_i;
  state_i = fc.state_i;
  zi_ges = fc.zi_ges;
  z_vc = fc.z_vc;
  vc_dummy1 = fc.vc_dummy1;
  spline_ges = fc.spline_ges;
  spline_ges2 = fc.spline_ges2;
  spline_zi = fc.spline_zi;
  ges_index = fc.ges_index;
  int_ti_help = fc.int_ti_help;
  baselinep = fc.baselinep;

  return *this;
  }


void pspline_multibaseline::outoptions(void)
  {

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
  }







void pspline_multibaseline::update(void)
  {



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
//----------------------------------------------------------
    update_multibaseline();
//----------------------------------------------------------

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

//---------Integral berechnen für vorgeschlagenes beta--------
      update_multibaseline();
//------------------------------------/

      }

    an+=blocksize;
    if (j == matquant[blocksize-min]-2)
      en = nrpar;
    else
      en+=blocksize;

    } // end: for(j=0;j<matquant[blocksize-min];j++)  */

  if (center)
    {
    compute_intercept();

    for(i=0;i<nrpar;i++)
      beta(i,0) -= intercept;

//-------------------Spline zentrieren------------------------------------------



        for(i=0;i<2.0*likep->get_nrobs();i++)
          {
          spline_ges(i,0) -= intercept;
          spline_ges2(i,0) -= intercept;
          }


//------------------------------------------------------------------------------

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

//--------Für's DIC-------------------------------------------------------------
void pspline_multibaseline::compute_int_ti_mean(void)
  {
  testmat.mult(spline_ges,betamean);
  testmat.mult_index(spline_ges2,betamean);
  compute_int_ti(betamean);
  }





//--------berechnet int_0^{t_i}/exp(logbaseline(t_i))---------------------------
//--------mit Hilfe von geordneten t_i und Knoten-------------------------------
void pspline_multibaseline::compute_int_ti(const datamatrix & b)
{

//------------------------left truncation----------------------------
if(begin0==false)
  {
  double * int_D_help;
  double * betap;
  double dist_knots = int_knots(1,0)-int_knots(0,0);
  unsigned i,j,k;
  k=1;
  double erg,spline_u,spline_o;
  erg = 0.0;
  double * int_ti_p = likep->get_integral_ti()+col;
  double * int_ti_help_p = int_ti_help.getV();
  double * int_ti_help_p2 = int_ti_help.getV();

  spline_o=0.0;
  spline_u=0.0;
  int_D_help = int_D.getV();
  betap=b.getV();

  for(j=0;j<nrpar;j++,int_D_help++,betap++)
    spline_o += *betap* *int_D_help;
  spline_u=spline_o;

//--------------------------------erster Integralwert---------------------------

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

  int_ti_p=likep->get_integral_ti()+(ges_index(0,0)*likep->get_responsedim()+col);
  *int_ti_p =erg*0.5/(exp(spline_ges(0,0)));

  int_ti_help_p=int_ti_help.getV()+ges_index(0,0);
     *int_ti_help_p = erg*0.5;

//------------------------------------------------------------

  for(i=1;i<zi_ges.rows();i++)
    {
    if(k==int_knots.rows())
      k=int_knots.rows()-1;
    if(k<int_knots.rows() && zi_ges(ges_index(i,0),0)<=int_knots(k,0))
      erg=erg+(zi_ges(ges_index(i,0),0)-zi_ges(ges_index(i-1,0),0))*(exp(spline_ges(i-1,0))+exp(spline_ges(i,0)));
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

    int_ti_p=likep->get_integral_ti()+ges_index(i,0)*likep->get_responsedim()+col;
    *int_ti_p =erg*0.5/(exp(spline_ges(i,0)));

    int_ti_help_p=int_ti_help.getV()+ges_index(i,0);
    *int_ti_help_p =erg*0.5;
    }
//------------------------------------------------------------------------------

  i=0;
  for(i=likep->get_nrobs();i<2*likep->get_nrobs();i++)
    {
    if(zi_ges(i,0)!=0)
      {
      int_ti_p=likep->get_integral_ti()+(i-likep->get_nrobs())*likep->get_responsedim()+col;
      int_ti_help_p=int_ti_help.getV()+i-likep->get_nrobs();
      int_ti_help_p2=int_ti_help.getV()+i;
      *int_ti_p = (*int_ti_help_p-*int_ti_help_p2)/exp(spline_ges2(i-likep->get_nrobs(),0));
      assert(*int_ti_p>=0.0);
      }
    }
  }//left_trunc



} //compute_ti











void pspline_multibaseline::update_multibaseline()
{
//---------Integral berechnen---------------------------------

    testmat.mult(spline_ges,beta);
    testmat.mult_index(spline_ges2,beta);
    compute_int_ti(beta);
}
//------------------------------------------------------------------------------









} // end: namespace MCMC



