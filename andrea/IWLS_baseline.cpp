//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include "IWLS_baseline.h"



namespace MCMC
{


//------------------------------------------------------------------------------
//---------- CLASS: IWLS_baseline (implementation of member functions) ------
//------------------------------------------------------------------------------


IWLS_baseline::IWLS_baseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & d,const bool & mode,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l,const fieldtype & ft,const ST::string & monotone,
                    const unsigned & upW, const bool & updatetau, const double & fstart,
                    const double & a, const double & b, const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const bool & diag, const unsigned & c,const datamatrix & anfang)
  : IWLS_pspline(o,dp,fcc,d,mode,nrk,degr,kp,l,ft,monotone,upW,updatetau,fstart,ti,fp,pres,deriv,gs,diag,c)
  {

  a_invgamma = a;
  b_invgamma = b;
//  prec_env = envmatdouble(0,nrpar,2);

  unsigned i,j;
  gauss_n = 9;
  vc_dummy1 = false;
  baseline = true;
  baselinep = vector<IWLS_baseline*>(0);

  response_help = datamatrix(d.rows(),1,0);
  Xdelta = datamatrix(beta.rows(),1,0);
  Adelta = datamatrix(beta.rows(),1,0);
  Eins = datamatrix(d.rows(),1,1.0);
  for(int i=0;i<d.rows();i++)
    response_help(i,0)=likep->get_response(index(i,0),0);

  compute_XWtildey(Eins,response_help,1.0);
  Xdelta = muy;
  ofstream test3("d:\\temp\\Xdelta.txt");
  Xdelta.prettyPrint(test3);
  test3.close();
  zi = d;

//-----------------linkstrunkiert oder zeitl. variierende Kovariablen?--------
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

  vector<datamatrix> gaussy(gauss_n);
  for(i=0;i<gauss_n;i++)
    {
    gaussy[i] = datamatrix(zi.rows()+1,1,0);
    }

  for(i=0;i<zi.rows();i++)
    {
    zi_ges(i,0) = zi(i,0);
    zi_ges(zi.rows()+i,0) = beg_i(i,0);
    }
  ges_index = statmatrix<int>(zi_ges.rows(),1);
  ges_index.indexinit();
  zi_ges.indexsort(ges_index,0,zi_ges.rows()-1,0,0);

  testmat = MCMC::bsplinemat(zi_ges,nrk,degr,kp,true);

  coeff = datamatrix(gauss_n,1,0);
  coeff(0,0) = 0.330239355001260;
  coeff(1,0) = 0.312347077040003;
  coeff(2,0) = 0.312347077040003;
  coeff(3,0) = 0.260610696402935;
  coeff(4,0) = 0.260610696402935;
  coeff(5,0) = 0.180648160694857;
  coeff(6,0) = 0.180648160694857;
  coeff(7,0) = 0.081274388361574;
  coeff(8,0) = 0.081274388361574;

  double help1,help2;
  for(i=0;i<zi.rows();i++)
    {
    help1 = (zi(i,0)-beg_i(i,0))*0.5;
    help2 = (zi(i,0)+beg_i(i,0))*0.5;
//    if(gauss_n==3)
//      {
//      gaussy[0](i,0)= help2;
//      gaussy[1](i,0)= help2 + 0.774596669241483*help1;
//      gaussy[2](i,0)= help2 - 0.774596669241483*help1;
//      }
//    if(gauss_n==5)
//      {
//      gaussy[0](i,0) = help2;
//      gaussy[1](i,0) = help2 + 0.538469310105683*help1;
//      gaussy[2](i,0) = help2 - 0.538469310105683*help1;
//      gaussy[3](i,0) = help2 + 0.906179845938664*help1;
//      gaussy[4](i,0) = help2 - 0.906179845938664*help1;
//      }
    gaussy[0](i,0) = help2;
    gaussy[1](i,0) = help2 + 0.324253423403809*help1;
    gaussy[2](i,0) = help2 - 0.324253423403809*help1;
    gaussy[3](i,0) = help2 + 0.613371432700590*help1;
    gaussy[4](i,0) = help2 - 0.613371432700590*help1;
    gaussy[5](i,0) = help2 + 0.836031107326636*help1;
    gaussy[6](i,0) = help2 - 0.836031107326636*help1;
    gaussy[7](i,0) = help2 + 0.968160239507626*help1;
    gaussy[8](i,0) = help2 - 0.968160239507626*help1;
    }

  double maxzi=0.0;
  for(i=0;i<zi.rows();i++)
    if (zi(i,0)>maxzi) maxzi=zi(i,0);

  gaussmat = vector<MCMC::bsplinemat>(gauss_n);
  for(i=0;i<gauss_n;i++)
    {
    gaussy[i](zi.rows(),0)=maxzi;
    gaussmat[i] = MCMC::bsplinemat(gaussy[i],nrk,degr,kp,true);
    }
//------------------Designmatrix int_D für P-Spline an Knoten-------------------

  double knot_min = 0.0;
  double knot_max = zi.max(0);
  int_knots=datamatrix (50,1,0);
  for(int j=0;j<int_knots.rows();j++)
    int_knots(j,0) = knot_min + j*(knot_max-knot_min)/double(int_knots.rows()-1);

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
//------------------------------------------------------------------------------

  spline_ges = datamatrix(2*likep->get_nrobs(),1,0);
  int_ti_help = datamatrix(2*likep->get_nrobs(),1,0);
  gaussspline = datamatrix(zi.rows()+1,gauss_n,0);
  spline_zi = datamatrix(likep->get_nrobs(),1,0);
//------------------------------------------------------------------------

  A = datamatrix(beta.rows()-2,beta.rows(),0);
  for(int i=0;i<beta.rows()-2;i++)
    {
    A(i,i)=1.0/6.0;
    A(i,i+1)=2.0/3.0;
    A(i,i+2)=1.0/6.0;
    }
  An = datamatrix(d.rows(),beta.rows(),0);
  score = datamatrix(beta.rows(),1,0);
  AWA = datamatrix(beta.rows(),beta.rows(),0);
  Wbase = datamatrix(beta.rows()-2,1,0);
//  Wbase = datamatrix(d.rows(),1,0);
  distance = datamatrix(beta.rows()-2,1,0);
  distance(0,0)=(knot[4]-knot[3])/2.0;
  for(int i=1;i<beta.rows()-3;i++)
    {
    distance(i,0)=(knot[i+4]-knot[i+2])/2.0;
    }
  distance(beta.rows()-3,0)=(knot[beta.rows()]-knot[beta.rows()-1])/2.0;

  ofstream test("d:\\temp\\dist.txt");
  distance.prettyPrint(test);
  test.close();

  interval = datamatrix(zi.rows(),1,0);
  for(int i=0;i<zi.rows();i++)
    {
    for(int j=0;j<beta.rows()-2;j++)
      {
      double testknot=knot[j+3];
      if(zi(index(i,0),0)<=((knot[j+4]+knot[j+3])/2.0) && zi(index(i,0),0) > ((knot[j+3]+knot[j+2])/2.0))
        interval(i,0) = j;
      }
    }
  ofstream test2("d:\\temp\\int.txt");
  interval.prettyPrint(test2);
  test2.close();

  for(int i=0;i<d.rows();i++)
    {
    for(int j=0;j<beta.rows()-2;j++)
      {
      if(interval(i,0)==j)
        {
        An(i,j)=1.0/6.0;
        An(i,j+1)=2.0/3.0;
        An(i,j+2)=1.0/6.0;
        }
      }
    }

  Adelta.mult(An.transposed(),response_help);
  ofstream Adeltaof("d:\\temp\\Adelta.txt");
  Adelta.prettyPrint(Adeltaof);
  Adeltaof.close();

  deltaexact = datamatrix(zi.rows(),beta.rows()-2,0);
  for(int i=0;i<zi.rows();i++)
    {
    for(int j=0;j<beta.rows()-2;j++)
      {
      if(interval(i,0)>j)
        deltaexact(i,j)=distance(j,0);
      else
        {
        if(interval(i,0)==j)
          {
          if(j==0)
            deltaexact(i,j)=zi(index(i,0),0);
          else
            deltaexact(i,j)=zi(index(i,0),0) - ((knot[j+3]+knot[j+2])/2.0);
          }
        else
          deltaexact(i,j)=0.0;
        }
      }
    }
  ofstream exactof("d:\\temp\\deltaexact.txt");
  deltaexact.prettyPrint(exactof);
  exactof.close();

  DeltaN = datamatrix(d.rows(),d.rows(),0);

  for(int i=0;i<d.rows();i++)
    DeltaN(i,0)=zi(index(0,0),0);

  for(int i=1;i<d.rows();i++)
    {
    for(int j=1;j<=i;j++)
      {
       DeltaN(i,j)=zi(index(j,0),0)-zi(index(j-1,0),0);
      }
    }

  ofstream DeltaNof("d:\\temp\\deltaN.txt");
  DeltaN.prettyPrint(DeltaNof);
  DeltaNof.close();


  }


IWLS_baseline::IWLS_baseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & effmod,const datamatrix & intact,const bool & mode,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l, const fieldtype & ft, const ST::string & monotone,
                    const unsigned & upW, const bool & updatetau, const double & fstart,
                    const double & a, const double & b, const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const bool & diag, const unsigned & c,const datamatrix & anfang)
  : IWLS_pspline(o,dp,fcc,effmod,intact,mode,nrk,degr,kp,l,ft,monotone,upW,updatetau,fstart,a,b,ti,fp,pres,deriv,gs,diag,c)
  {

  }


IWLS_baseline::IWLS_baseline(const IWLS_baseline & fc)
  :IWLS_pspline(IWLS_pspline(fc))
  {

  begin0 = fc.begin0;
  int_knots = fc.int_knots;
  int_D = fc.int_D;
  testmat = fc.testmat;
  gaussmat = fc.gaussmat;
  coeff = fc.coeff;
  gauss_n = fc.gauss_n;
  zi=fc.zi;
  beg_i = fc.beg_i;
  zi_ges = fc.zi_ges;
  z_vc = fc.z_vc;
  vc_dummy1 = fc.vc_dummy1;
  spline_ges = fc.spline_ges;
  spline_zi = fc.spline_zi;
  gaussspline = fc.gaussspline;
  ges_index = fc.ges_index;
  int_ti_help = fc.int_ti_help;
  baselinep = fc.baselinep;
  A = fc.A;
  distance = fc.distance;
  interval = fc.interval;
  AWA = fc.AWA;
  Wbase = fc.Wbase;
  response_help = fc.response_help;
  Xdelta = fc.Xdelta;
  Eins = fc.Eins;
  score = fc.score;
  deltaexact = fc.deltaexact;
  An = fc.An;
  Adelta = fc.Adelta;
  DeltaN = fc.DeltaN;
  }


const IWLS_baseline & IWLS_baseline::operator=(const IWLS_baseline & fc)
  {
  if (this == &fc)
    return *this;
  IWLS_pspline::operator=(IWLS_pspline(fc));

  begin0 = fc.begin0;
  int_knots = fc.int_knots;
  int_D = fc.int_D;
  testmat = fc.testmat;
  gaussmat = fc.gaussmat;
  coeff = fc.coeff;
  gauss_n = fc.gauss_n;
  zi=fc.zi;
  beg_i = fc.beg_i;
  zi_ges = fc.zi_ges;
  z_vc = fc.z_vc;
  vc_dummy1 = fc.vc_dummy1;
  spline_ges = fc.spline_ges;
  spline_zi = fc.spline_zi;
  gaussspline = fc.gaussspline;
  ges_index = fc.ges_index;
  int_ti_help = fc.int_ti_help;
  baselinep = fc.baselinep;
  A = fc.A;
  distance = fc.distance;
  interval = fc.interval;
  AWA = fc.AWA;
  Wbase = fc.Wbase;
  response_help = fc.response_help;
  Xdelta = fc.Xdelta;
  Eins = fc.Eins;
  score = fc.score;
  deltaexact = fc.deltaexact;
  An = fc.An;
  Adelta = fc.Adelta;
  DeltaN = fc.DeltaN;
  return *this;
  }


bool IWLS_baseline::posteriormode_converged(const unsigned & itnr)
  {
//  if(increasing || decreasing || diagtransform)
  if(diagtransform)
    return true;
  else
    return likep->posteriormode_converged_fc(beta,beta_mode,itnr);
  }


bool IWLS_baseline::posteriormode(void)
  {

//  compute_XWXenv(likep->get_weightiwls(),column);
  prec_env.addto(XX_env,Kenv,1.0,lambda);
  lambda_prec = lambda;

  likep->substr_linearpred_m(spline,column,true);
  likep->compute_workingresiduals(column);

//  compute_XWtildey(likep->get_weightiwls(),likep->get_workingresiduals(),1.0,column);

  prec_env.solve(muy,beta);
  add_linearpred_multBS(beta,true);

  if(center)
    {
    compute_intercept();
    for(unsigned i=0;i<nrpar;i++)
      beta(i,0) -= intercept;
    fcconst->posteriormode_intercept(intercept);
    for(unsigned i=0;i<likep->get_nrobs();i++)
      spline(i,0) -= intercept;
    intercept = 0.0;
    }

  write_spline();
  write_derivative();

  if(derivative)
    fcderivative.posteriormode();

  fchelp.posteriormode();
  return FULLCOND::posteriormode();

  }


void IWLS_baseline::update(void)
  {

  if(lambdaconst == true)
    sigma2 = likep->get_scale(column)/lambda;

  if(optionsp->get_nriter()==1)       // posterior mode Schätzung übernehmen
    betaold.assign(beta);


  if(utype == iwls)
    {
    update_IWLS();
    }
  else if(utype == iwlsmode)
    {
    update_IWLS_mode();
    }


  if(predictright || predictleft)
    update_prediction();

//  spline in fchelp schreiben
  if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
      ((optionsp->get_nriter()-optionsp->get_burnin()-1) % optionsp->get_step() == 0) )
    {
    if(diagtransform)
      write_spline(G*beta);
    else
      write_spline();

    write_derivative();
    }

  if(derivative)
    fcderivative.update();

  fchelp.update();
  FULLCOND::update();

  } // end: function update


void IWLS_baseline::update_IWLS_mode(void)
  {

  unsigned i;

  double logold = likep->loglikelihood(true)
                - 0.5*Kenv.compute_quadformblock(betaold,0,nrparpredictleft,nrpar-nrparpredictright-1)/sigma2;

  add_linearpred_multBS(beta_mode,betaold,true);
//  multBS_index(spline,beta_mode);

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
//    likep->compute_weight(W,column,true);
//    likep->tilde_y(mu,spline,column,true,W);
    likep->compute_IWLS_weight_tildey(W,mu,column,true);
    mu.plus(spline,mu);
//    compute_XWXenv(W);
//    compute_XWtildey(W,1.0);
    compute_XWXenv_XWtildey(W,1.0);
    }
  else
    {
    likep->tilde_y(mu,spline,column,true,W);
    compute_XWtildey(W,1.0);
    }

  prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);

  prec_env.solve(muy,betahelp);

  double * work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  prec_env.solveU(beta,betahelp);

  add_linearpred_multBS(beta,beta_mode,true);
  beta_mode.assign(betahelp);

  betahelp.minus(beta,beta_mode);
  double qold = - 0.5*prec_env.compute_quadformblock(betahelp,0,nrparpredictleft,nrpar-nrparpredictright-1);

  betahelp.minus(betaold,beta_mode);
  double lognew = likep->loglikelihood(true)
                - 0.5*Kenv.compute_quadformblock(beta,0,nrparpredictleft,nrpar-nrparpredictright-1)/sigma2;
  double qnew = - 0.5*prec_env.compute_quadformblock(betahelp,0,nrparpredictleft,nrpar-nrparpredictright-1);

  double alpha = lognew + qnew - logold - qold;
  double u = log(uniform());

  if(center)
    {
    compute_intercept();
    for(i=0;i<nrpar;i++)
      beta_mode(i,0) -= intercept;
    intercept = 0.0;
    }

  if(u<=alpha)
    {
    acceptance++;
    if(center)
      {
      compute_intercept();
      for(i=0;i<nrpar;i++)
        beta(i,0) -= intercept;
      fcconst->update_intercept(intercept);
      for(i=0;i<likep->get_nrobs();i++)
        spline(i,0) -= intercept;
      intercept = 0.0;
      }
    betaold.assign(beta);
    }
  else
    {
    add_linearpred_multBS(betaold,beta,true);
    beta.assign(betaold);
    }

  }


void IWLS_baseline::update_IWLS(void)
  {


  unsigned i;

/*  if(optionsp->get_nriter()==1)
    {
    beta(0,0) = -6.481843;
    beta(1,0) = -5.261254;
    beta(2,0) = -4.035323;
    beta(3,0) = -2.855564;
    beta(4,0) = -1.813996;
    beta(5,0) = -1.116795;
    beta(6,0) = -0.668608;
    beta(7,0) = -0.3481318;
    beta(8,0) = -0.0590283;
    beta(9,0) =  0.2067058;
    beta(10,0) = 0.3909008;
    beta(11,0) = 0.5491;
    beta(12,0) = 0.7070876;
    beta(13,0) = 0.8955011;
    beta(14,0) = 1.092252;
    beta(15,0) = 1.306334;
    beta(16,0) = 1.51953;
    beta(17,0) = 1.702514;
    beta(18,0) = 1.890162;
    beta(19,0) = 2.137843;
    beta(20,0) = 2.416735;
    beta(21,0) = 2.70698;
    add_linearpred_multBS(beta);
    betaold.assign(beta);
    }*/

  double logold = - 0.5*Kenv.compute_quadformblock(betaold,0,nrparpredictleft,nrpar-nrparpredictright-1)/sigma2;

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
    multBS(spline,beta);
    compute_Wbase();
    likep->tilde_y(mu,spline,column,true,W);
    update_baseline();
    logold += likep->loglikelihood();
//    mu.plus(spline,mu);

//    compute_XWXenv(W);
//    compute_XWtildey(W,1.0);
    compute_XWXenv_XWtildey(W,1.0);
    }
  else
    {
//    likep->tilde_y(mu,spline,column,true,W);
    logold += likep->compute_IWLS(W,mu,false,column,true);
    mu.plus(mu,spline);

    compute_XWtildey(W,1.0);
    }

  prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);

  prec_env.solve(muy,betahelp);

  double * work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  prec_env.solveU(beta,betahelp);

  ofstream allesout("d:\\temp\\alles.txt");
  allesout<<"prec_env"<<endl;
  prec_env.print4(allesout);
  allesout<<"betahelp"<<endl;
  betahelp.prettyPrint(allesout);
  allesout<<"W"<<endl;
  (W.transposed()).prettyPrint(allesout);
  allesout<<"linpred"<<endl;
  (likep->get_linearpred().transposed()).prettyPrint(allesout);
  allesout<<"beta"<<endl;
  beta.prettyPrint(allesout);


  add_linearpred_multBS(beta,betaold,true);
  betahelp.minus(beta,betahelp);

  double qold = - 0.5*prec_env.compute_quadformblock(betahelp,0,nrparpredictleft,nrpar-nrparpredictright-1);

  double lognew = - 0.5*Kenv.compute_quadformblock(beta,0,nrparpredictleft,nrpar-nrparpredictright-1)/sigma2;

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
    qold += 0.5*prec_env.getLogDet();

    multBS(spline,beta);
    allesout<<"spline"<<endl;
    spline.prettyPrint(allesout);
    allesout.close();
    compute_Wbase();
    likep->tilde_y(mu,spline,column,true,W);
    update_baseline();
    lognew += likep->compute_IWLS(W,mu,true,column,true);
//    mu.plus(spline,mu);

//    compute_XWXenv(W);
//    compute_XWtildey(W,1.0);
    compute_XWXenv_XWtildey(W,1.0);
    prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);
    }
  else
    {
//    likep->tilde_y(mu,spline,column,true,W);
    lognew += likep->compute_IWLS(W,mu,false,column,true);
    mu.plus(mu,spline);

    compute_XWtildey(W,1.0);
    }

  prec_env.solve(muy,betahelp);

  betahelp.minus(betaold,betahelp);
  double qnew = - 0.5*prec_env.compute_quadformblock(betahelp,0,nrparpredictleft,nrpar-nrparpredictright-1);

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
    qnew += 0.5*prec_env.getLogDet();
    }

  double alpha = lognew + qnew - logold - qold;
  double u = log(uniform());

  if(u<=alpha)
    {
    acceptance++;
    if(center)
      {
      compute_intercept();
      for(i=0;i<nrpar;i++)
        beta(i,0) -= intercept;
      fcconst->update_intercept(intercept);
      for(i=0;i<likep->get_nrobs();i++)
        spline(i,0) -= intercept;
      intercept = 0.0;
      }
    betaold.assign(beta);
    }
  else
    {
    add_linearpred_multBS(betaold,beta,true);
    beta.assign(betaold);
    }
   update_baseline();



/*  unsigned i;

  double logold = - 0.5*Kenv.compute_quadformblock(betaold,0,nrparpredictleft,nrpar-nrparpredictright-1)/sigma2;

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
//    likep->compute_weight(W,column,true);
//    likep->tilde_y(mu,spline,column,true,W);
    update_baseline();
    logold += likep->loglikelihood(true);

//    mu.plus(A'beta,mu);

    multBS(spline,beta);
    compute_score();
//    compute_Wbase();
    compute_AWA();

    muy = Xdelta - A.transposed()*Wbase + AWA*beta;
//    muy = Adelta - A.transposed()*Wbase + AWA*beta;

    ofstream test5("d:\\temp\\muy.txt");
    muy.prettyPrint(test5);
    test5.close();

    XX_env = envmatrix<double>(AWA);
    }
  else
    {
//    likep->tilde_y(mu,spline,column,true,W);
    update_baseline();
    logold += likep->compute_IWLS(W,mu,false,column,true);
//    mu.plus(mu,spline);

    compute_XWtildey(W,1.0);
    }

  prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);

  ofstream test7("d:\\temp\\prec.txt");
    test7<<"prec_env"<<endl;
    prec_env.print4(test7);
    test7<<endl<<"Kenv"<<endl;
    Kenv.print4(test7);
    test7<<endl<<"XX_env"<<endl;
    XX_env.print4(test7);
    test7<<"Wbase"<<endl;
    Wbase.prettyPrint(test7);
  test7.close();

  prec_env.solve(muy,betahelp);

  double * work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  prec_env.solveU(beta,betahelp);

  add_linearpred_multBS(beta,betaold,true);
  betahelp.minus(beta,betahelp);

  double qold = - 0.5*prec_env.compute_quadformblock(betahelp,0,nrparpredictleft,nrpar-nrparpredictright-1);

  double lognew = - 0.5*Kenv.compute_quadformblock(beta,0,nrparpredictleft,nrpar-nrparpredictright-1)/sigma2;

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
    qold += 0.5*prec_env.getLogDet();

    update_baseline();
    lognew += likep->loglikelihood(true);

    multBS(spline,beta);
    compute_score();
//    compute_Wbase();
    compute_AWA();

    muy = Xdelta - A.transposed()*Wbase + AWA*beta;
//    muy = Adelta - A.transposed()*Wbase + AWA*beta;
    XX_env = envmatrix<double>(AWA);

    prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);

    ofstream test8("d:\\temp\\prec2.txt");
      test8<<"prec_env"<<endl;
      prec_env.print4(test8);
      test8<<endl<<"muy"<<endl;
      muy.prettyPrint(test8);
      test8<<endl<<"XX_env"<<endl;
      XX_env.print4(test8);
      test8<<"Wbase"<<endl;
      Wbase.prettyPrint(test8);
    test8.close();


    }
  else
    {
//    likep->tilde_y(mu,spline,column,true,W);
    update_baseline();
    lognew += likep->compute_IWLS(W,mu,false,column,true);
//    mu.plus(mu,spline);

    compute_XWtildey(W,1.0);
    }

  prec_env.solve(muy,betahelp);

  betahelp.minus(betaold,betahelp);
  double qnew = - 0.5*prec_env.compute_quadformblock(betahelp,0,nrparpredictleft,nrpar-nrparpredictright-1);

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
    qnew += 0.5*prec_env.getLogDet();
    }

  double alpha = lognew + qnew - logold - qold;
  double u = log(uniform());

  if(u<=alpha)
    {
     acceptance++;
    if(center)
      {
      compute_intercept();
      for(i=0;i<nrpar;i++)
        beta(i,0) -= intercept;
      fcconst->update_intercept(intercept);
      for(i=0;i<likep->get_nrobs();i++)
        spline(i,0) -= intercept;
      intercept = 0.0;
      }
    betaold.assign(beta);
    }
  else
    {
    add_linearpred_multBS(betaold,beta,true);
    beta.assign(betaold);
    }
  update_baseline();     */
  }


void IWLS_baseline::predict(const datamatrix & newX, datamatrix & linpred)
  {

  unsigned i,j;
  datamatrix betac(beta.rows(),beta.cols());
  datamatrix bspline(1,nrpar,0);
  double * worklin = linpred.getV();

  for(i=0;i<nrpar;i++)
    bspline(0,i) = bspline_rek(degree,i,newX);

  if(varcoeff)
    {
    double help;
    for(i=0;i<optionsp->get_samplesize();i++,worklin++)
      {
      help = 0.0;
      readsample2(betac,i);
      for(j=0;j<nrpar;j++)
        {
        help += betac(j,0) * bspline(0,j);
        }
      *worklin += help * newX(0,1);
      }
    }
  else
    {
    for(i=0;i<optionsp->get_samplesize();i++,worklin++)
      {
      readsample2(betac,i);
      for(j=0;j<nrpar;j++)
        {
        *worklin += betac(j,0) * bspline(0,j);
        }
      }
    }

  }


double IWLS_baseline::compute_quadform(void)
  {
  if(predictright || predictleft)
    {
/*
    double quadform = Kenv.compute_quadform(beta,0);
    if(predictright)
      for(unsigned i=0;i<nrparpredictright;i++)
        quadform -= beta(nrpar-1-i,0)*beta(nrpar-1-i,0);
    if(predictleft)
      for(unsigned i=0;i<nrparpredictleft;i++)
        quadform -= beta(i,0)*beta(i,0);
    return quadform;
*/
    return Kenv.compute_quadformblock(beta,0,nrparpredictleft,nrpar-nrparpredictright-1);
    }
  else
    return Kenv.compute_quadform(beta,0);
  }


//--------- für baseline -------------------------------------------------------

//--------Für's DIC-------------------------------------------------------------
void IWLS_baseline::compute_int_ti_mean(void)
  {
  unsigned i;
  if(baselinep.size()>1)
    {
    if(vc_dummy1==true)
      {
      vector <double *> splinevec;
      vector <double *> betavec;
      for(i=0;i<(baselinep.size());i++)
        splinevec.push_back(baselinep[i]->get_spline_zi_mean());
      for(i=0;i<(baselinep.size());i++)
        betavec.push_back(baselinep[i]->get_betamean());
      compute_int_ti_vc_di0(splinevec,betavec);
      for(i=1;i<baselinep.size();i++)
        {
        compute_int_ti_vc_di(i,splinevec,betavec);
        }
      }
    else
      {
      compute_int_gauss_DIC();
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


void IWLS_baseline::compute_int_ti_linear(const double & b)
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


//--------berechnet int_0^{t_i}/exp(logbaseline(t_i))---------------------------
//--------mit Hilfe von geordneten t_i und Knoten-------------------------------
void IWLS_baseline::compute_int_ti(const datamatrix & b)
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
  double * int_ti_p = likep->get_integral_ti();
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

    int_ti_p=likep->get_integral_ti()+ges_index(i,0);
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
      int_ti_p=likep->get_integral_ti()+i-likep->get_nrobs();
      int_ti_help_p=int_ti_help.getV()+i-likep->get_nrobs();
      int_ti_help_p2=int_ti_help.getV()+i;

      double help = (*int_ti_help_p-*int_ti_help_p2)/exp(spline_ges(ges_index(i-likep->get_nrobs(),0),0));
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

//  ofstream testof("d:\\temp\\beta.txt");
//  for(i=0;i<22;i++)
//    testof<<b(i,0)<<endl;
//  testof.close();

  double erg,spline_u,spline_o;
  erg = 0.0;
  double * int_ti_p=likep->get_integral_ti();
  double * int_ti_help_p=int_ti_help.getV();

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
      erg=erg+(int_knots(k,0)-zi(index(i-1,0),0))*(exp(spline(i-1,0))+exp(spline_o)) ;

      k++;

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




void IWLS_baseline::compute_int_ti(unsigned beg)
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

    k++;
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




void IWLS_baseline::compute_int_ti_vc_di0(const vector<double *> splinevector,const vector<double *> betavector)
{
//ofstream oftest_spline ("f:\\baseline\\spline_di0.txt");
double * int_D_help;
double * betap;
double dist_knots=int_knots(1,0)-int_knots(0,0);
unsigned i,j,k,i_help,i_vc;
i=0;
i_help=0;
double erg,spline_u,spline_o,spline_ti_help,spline_ti;
double * int_ti_p=likep->get_integral_ti();
double * int_ti_help_p=int_ti_help.getV();
double * spline_help;
statmatrix<double*>z_vc_help;
z_vc_help = statmatrix<double*>(baselinep.size()-1,1);
for(i_vc=1;i_vc<baselinep.size();i_vc++)
  z_vc_help(i_vc-1,0) = baselinep[i_vc]->get_z_vc();
k=1;
erg=0.0;
spline_o=0.0;
spline_u=0.0;
spline_help = splinevector[0];

int_D_help =baselinep[0]->get_int_D();
betap = betavector[0];

for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help++,betap++)
  spline_o += *betap* *int_D_help;
spline_u=spline_o;
//------i0 berechnen---------
double z_vc_sum =0.0;
for(i_vc=1;i_vc<baselinep.size();i_vc++)
  z_vc_sum = z_vc_sum + *(z_vc_help(i_vc-1,0) + index(i,0));
while(z_vc_sum!=0.0 && i < index.rows()-1)
  {
  i++;
  spline_help++;
  z_vc_sum =0.0;
  for(i_vc=1;i_vc<baselinep.size();i_vc++)
    {
    double test1= *(z_vc_help(i_vc-1,0) + index(i,0));
    double test2 = index(i,0);
    z_vc_sum = z_vc_sum + *(z_vc_help(i_vc-1,0) + index(i,0));
    }
  }

i_help=i;
//------------------erster Integralwert------------------------

while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
  {
  spline_u=spline_o;
  spline_o=0.0;
  betap = betavector[0];
  for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help++,betap++)
    spline_o += *betap* *int_D_help;
  erg=erg+(exp(spline_u)+exp(spline_o));

  k++;
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
  z_vc_sum =0.0;
  for(i_vc=1;i_vc<baselinep.size();i_vc++)
    z_vc_sum = z_vc_sum + *(z_vc_help(i_vc-1,0) + index(i,0));
  if(z_vc_sum!=0.0)
    {
    i++;
    spline_help++;
    }
  else
    {
    if(k == int_knots.rows())
      k=int_knots.rows()-1;
    spline_ti=*spline_help;
    if(k<int_knots.rows() && zi(index(i,0),0)<=int_knots(k,0))
      erg=erg+(zi(index(i,0),0)-zi(index(i_help,0),0))*(exp(spline_ti_help)+exp(spline_ti));
    else
      {
      spline_u=spline_o;
      spline_o=0.0;
      betap = betavector[0];
      for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help++,betap++)
        spline_o += *betap* *int_D_help;
      erg=erg+(int_knots(k,0)-zi(index(i_help,0),0))*(exp(spline_ti_help)+exp(spline_o)) ;

      k++;

      while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
        {
        spline_u=spline_o;
        spline_o=0.0;
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
//oftest_spline.close();  
}



void IWLS_baseline::compute_int_ti_vc_di(const int dummy, const vector<double *> splinevector,const vector<double *> betavector)
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
int_D_help_1=statmatrix<double*>(2,1);
statmatrix<double*> spline_zi_help;
spline_zi_help=statmatrix<double*>(2,1);
double help;
int i_vc;
k=1;
erg=0.0;
spline_o=0.0;
spline_u=0.0;
z_vc_help = baselinep[dummy]->get_z_vc();
spline_zi_help(0,0) = splinevector[0];
int_D_help_1(0,0)= baselinep[0]->get_int_D();
spline_zi_help(1,0) = splinevector[dummy];
int_D_help_1(1,0)= baselinep[dummy]->get_int_D();

betap = betavector[0];
help=0.0;
for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help_1(0,0)++,betap++)
  help+=*betap* *int_D_help_1(0,0);
spline_o +=help;

betap = betavector[dummy];
help=0.0;
for(j=0;j<baselinep[dummy]->get_nrpar();j++,int_D_help_1(1,0)++,betap++)
  help+=*betap* *int_D_help_1(1,0);
spline_o +=help;

spline_u=spline_o;

//------kleinstes ti mit z_vc==1 und splines an der Stelle suchen---------
while(*(z_vc_help+index(i,0))==0.0)
  {
  i++;
  for(i_vc=0;i_vc<2;i_vc++)
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

  betap = betavector[0];
  help=0.0;
  for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help_1(0,0)++,betap++)
    help+=*betap* *int_D_help_1(0,0);
  spline_o +=help;

  betap = betavector[dummy];
  help=0.0;
  for(j=0;j<baselinep[dummy]->get_nrpar();j++,int_D_help_1(1,0)++,betap++)
    help+=*betap* *int_D_help_1(1,0);
  spline_o +=help;

  erg=erg+(exp(spline_u)+exp(spline_o));

  k++;
  }
erg=erg*dist_knots;

//---------------Spline an der Stelle ti ---------------
for(i_vc=0;i_vc<2;i_vc++)
  {
  help= *spline_zi_help(i_vc,0);
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
for(i_vc=0;i_vc<2;i_vc++)
  spline_zi_help(i_vc,0)++;
while(i<zi.rows())
  {
//-----------falls z_vc==0: zu nächster Beobachtung---------------
  if(*(z_vc_help+index(i,0))==0)
    {
    i++;
    for(i_vc=0;i_vc<2;i_vc++)
      spline_zi_help(i_vc,0)++;
    }
//-------------falls z_vc==1: Integral berechenen---------------------
  else
    {
    if(k == int_knots.rows())
      k=int_knots.rows()-1;

//Spline an der Stelle ti
    spline_ti=0.0;
    for(i_vc=0;i_vc<2;i_vc++)
      {
      help= *spline_zi_help(i_vc,0);
      spline_ti +=help;
      }
    if(k<int_knots.rows() && zi(index(i,0),0)<=int_knots(k,0))
      erg=erg+(zi(index(i,0),0)-zi(index(i_help,0),0))*(exp(spline_ti_help)+exp(spline_ti)) ;
    else
      {
      spline_u=spline_o;
      spline_o = 0.0;

      betap = betavector[0];
      help = 0.0;
      for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help_1(0,0)++,betap++)
        help+=*betap* *int_D_help_1(0,0);
      spline_o +=help;

      betap = betavector[dummy];
      help=0.0;
      for(j=0;j<baselinep[dummy]->get_nrpar();j++,int_D_help_1(1,0)++,betap++)
        help+=*betap* *int_D_help_1(1,0);
      spline_o +=help;

      erg=erg+(int_knots(k,0)-zi(index(i_help,0),0))*(exp(spline_ti_help)+exp(spline_o)) ;

      k++;

      while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
        {
        spline_u=spline_o;
        spline_o=0.0;

        betap = betavector[0];
        help=0.0;
        for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help_1(0,0)++,betap++)
          help+=*betap* *int_D_help_1(0,0);
        spline_o +=help;

        betap = betavector[dummy];
        help=0.0;
        for(j=0;j<baselinep[dummy]->get_nrpar();j++,int_D_help_1(1,0)++,betap++)
          help+=*betap* *int_D_help_1(1,0);
        spline_o +=help;

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
    for(i_vc=0;i_vc<2;i_vc++)
      spline_zi_help(i_vc,0)++;
    }//--------------- else: d.h. z_vc(index(i,0),0)==1 -----------------------
  } //while
}//--------------compute_int_ti_vc_di---------------------



void IWLS_baseline::compute_int_gauss(void)
{

vector<double *> betavec;
unsigned i,j,k;
vector<double *> splinevec;
vector<double *> gausssplinevec;
vector<datamatrix > z_vc_help;

for(i=0;i<(baselinep.size());i++)
  splinevec.push_back(baselinep[i]->get_spline_ges());
for(i=0;i<(baselinep.size());i++)
  gausssplinevec.push_back(baselinep[i]->get_gaussspline());
for(i=0;i<(baselinep.size()-1);i++)
  z_vc_help.push_back(baselinep[i+1]->get_z_vc_np());

double * int_ti_p=likep->get_integral_ti();

double help,splinehelp1,splinehelp2;

for(i=0;i<zi.rows();i++,int_ti_p++)
  {
  help = 0.0;
  splinehelp2 = 0.0;
  for(j=0;j<gauss_n;j++)
    {
    splinehelp2=0.0;
    for(k=0;k<baselinep.size();k++)
      {
      splinehelp1 = *(gausssplinevec[k]);
      if(k>0) splinehelp1 = splinehelp1* (z_vc_help[k-1](i,0));
      splinehelp2 = splinehelp2 + splinehelp1;
      gausssplinevec[k]++;
      }
    help = help + coeff(j,0)*exp(splinehelp2);
    }
  splinehelp2 = 0.0;
  for(k=0;k<baselinep.size();splinevec[k]++,k++)
    {
    splinehelp1 = *(splinevec[k]);
    if(k>0)
      {
      splinehelp1 = splinehelp1* (z_vc_help[k-1](i,0));
      }
    splinehelp2 = splinehelp2 + splinehelp1;
    }
  *int_ti_p = ((zi(i,0)-beg_i(i,0))*0.5)*help/(exp(splinehelp2));
  }
}


void IWLS_baseline::compute_int_gauss_DIC(void)
{

unsigned i,j,k;
vector<double *> splinevec;
vector<double *> gausssplinevec;
vector<double *> z_vc_help;

for(i=0;i<(baselinep.size());i++)
  splinevec.push_back(baselinep[i]->get_spline_ges_mean());
for(i=0;i<(baselinep.size());i++)
  gausssplinevec.push_back(baselinep[i]->get_gaussspline_mean());
for(i=0;i<(baselinep.size()-1);i++)
  z_vc_help.push_back(baselinep[i+1]->get_z_vc());

double * int_ti_p=likep->get_integral_ti();

double help,splinehelp1,splinehelp2;

for(i=0;i<zi.rows();i++,int_ti_p++)
  {
  help = 0.0;
  splinehelp2 = 0.0;
  for(j=0;j<gauss_n;j++)
    {
    splinehelp2 = 0.0;
    for(k=0;k<baselinep.size();gausssplinevec[k]++,k++)
      {
      splinehelp1 = *(gausssplinevec[k]);
      if(k>0) splinehelp1 = splinehelp1* *(z_vc_help[k-1]);
      splinehelp2 = splinehelp2 + splinehelp1;
      }

    help = help + coeff(j,0)*exp(splinehelp2);
    }
  splinehelp2 = 0.0;
  for(k=0;k<baselinep.size();splinevec[k]++,k++)
    {
    splinehelp1 = *(splinevec[k]);
    if(k>0)
      {
      splinehelp1 = splinehelp1* *(z_vc_help[k-1]);
      z_vc_help[k]++;
      }
    splinehelp2 = splinehelp2 + splinehelp1;
    }
  *int_ti_p = ((zi(i,0)-beg_i(i,0))*0.5)*help/(exp(splinehelp2));
  }
}


void IWLS_baseline::update_baseline()
{
//---------Integral berechnen---------------------------------
int i;
if(baselinep.size()>=1)
  {
  if(vc_dummy1==true)   //keine Linkstrunkierung, zeitl. var. Effekt für dummykod. Variable
    {
    vector <double *> splinevec;
    vector <double *> betavec;
    for(i=0;i<(baselinep.size());i++)
      splinevec.push_back(baselinep[i]->get_spline_zi());
    for(i=0;i<(baselinep.size());i++)
      betavec.push_back(baselinep[i]->getbetapointer());
    compute_int_ti_vc_di0(splinevec,betavec);
    for(i=1;i<baselinep.size();i++)
      {
      compute_int_ti_vc_di(i,splinevec,betavec);
      }
    }
  else    //zeitl. var. Effekt für beliebige Kovariablen, Linkstrunkierung
    {
    compute_int_gauss();
    }
  }
else    //kein zeitl. var. Effekt
  {
  if(begin0==false)     //Linkstrunkierung
    {
    testmat.mult(spline_ges,beta);
    compute_int_ti(beta);
    }
  else    //keine Linkstrunkierung
    {
    multBS(spline,beta);
    compute_int_ti(beta);
    }
  }
}
//------------------------------------------------------------------------------


double * IWLS_baseline::get_gaussspline()
  {
  datamatrix egon (gaussspline.rows(),1,0);
  for(unsigned i=0;i<gauss_n;i++)
    {
    egon = gaussspline.getCol(i);
    gaussmat[i].mult_index(egon,beta);
    gaussspline.putCol(i,egon);
    }
  return gaussspline.getV();
  }

double * IWLS_baseline::get_gaussspline_mean()
  {
  datamatrix egon (gaussspline.rows(),1,0);
  for(unsigned i=0;i<gauss_n;i++)
    {
    egon = gaussspline.getCol(i);
    gaussmat[i].mult_index(egon,betamean);
    gaussspline.putCol(i,egon);
    }
  return gaussspline.getV();
  }


void IWLS_baseline::compute_Wbase(void)
  {
  int i;

  ofstream splineof ("d:\\temp\\spline.txt");
  spline.prettyPrint(splineof);
  splineof.close();


  datamatrix W2(zi.rows(),1,0.0);
  datamatrix etaminus;
  etaminus = datamatrix(zi.rows(),1,0.0);
  for(i=0;i<zi.rows();i++)
    etaminus(i,0) = exp(likep->get_linearpred(index(i,0),0)-spline(i,0));

  W2.mult(DeltaN.transposed(),etaminus);
  for(i=0;i<zi.rows();i++)
    W(index(i,0),0) = W2(i,0)* exp(spline(i,0));

/*  datamatrix betatilde;
//  datamatrix Wbase_help;
//  double Wbase_sum = 0.0;
  int i,j;
  double sumlinpred_help = 0.0;

  betatilde = datamatrix(beta.rows()-2,1,0);
//  Wbase_help = datamatrix(betatilde.rows(),1,0);
  betatilde.mult(A,beta);
  ofstream betaof ("d:\\temp\\betatilde.txt");
  betatilde.prettyPrint(betaof);
  betaof.close();
  ofstream splineof ("d:\\temp\\spline.txt");
  spline.prettyPrint(splineof);
  splineof.close();

  ofstream lpof ("d:\\temp\\sumlinpred.txt");

  for(i=betatilde.rows()-1;i>-1;i--)
    {
    for(j=0;j<zi.rows();j++)
      {
      double linpredtest= exp(likep->get_linearpred(index(j,0),0));
      double splinetest = exp(spline(j,0));
      double etaminustest = exp(likep->get_linearpred(index(j,0),0)-spline(j,0));
      double betatildetest= exp(betatilde(i,0));
      double distancetest= distance(i,0);
      if(interval(j,0)==i)
        sumlinpred_help=sumlinpred_help + exp(likep->get_linearpred(index(j,0),0)-spline(j,0));
      }
    Wbase(i,0) = exp(betatilde(i,0))*distance(i,0)*sumlinpred_help;
    lpof <<sumlinpred_help<<endl;
    }
  lpof.close();  */
  }

void IWLS_baseline::compute_AWA(void)
  {
  int i,j;
  datamatrix Wbase_help;
  datamatrix AWA_help;
  Wbase_help = datamatrix((beta.rows()-2),(beta.rows()-2),0);
  AWA_help = datamatrix(beta.rows(),beta.rows()-2,0);
  for(i=0;i<beta.rows()-2;i++)
    Wbase_help(i,i) = Wbase(i,0);
  AWA_help.mult(A.transposed(),Wbase_help);
  AWA.mult(AWA_help,A);
  }

void IWLS_baseline::compute_score(void)
  {
  int i,j;
  datamatrix betatilde;
  datamatrix AWA_test;
  AWA_test = datamatrix(beta.rows(),beta.rows(),0);
  double etaiminus=0.0;
  betatilde = datamatrix(beta.rows()-2,1,0);
  betatilde.mult(A,beta);
  ofstream betatof ("d:\\temp\\betatilde.txt");
  betatilde.prettyPrint(betatof);
  betatof.close();
  ofstream splinetof ("d:\\temp\\spline.txt");
  spline.prettyPrint(splinetof);
  splinetof.close();
  datamatrix sumlinpred;
  sumlinpred = datamatrix(betatilde.rows(),1,0);
  datamatrix Wtest;
  Wtest = datamatrix(betatilde.rows(),1,0);
  for(j=0;j<betatilde.rows();j++)
    {
    for(i=0;i<zi.rows();i++)
      {
      etaiminus = exp(likep->get_linearpred(index(i,0),0)-spline(i,0));
      sumlinpred(j,0) = sumlinpred(j,0)+ etaiminus*deltaexact(i,j);
      }
    }



  score(0,0) = Xdelta(0,0)-( 1.0/6.0*exp(betatilde(0,0))*sumlinpred(0,0) );
  score(1,0) = Xdelta(1,0)-( (2.0/3.0*exp(betatilde(0,0))*sumlinpred(0,0)) + (1.0/6.0*exp(betatilde(1,0))*sumlinpred(1,0)) );
  score(2,0) = Xdelta(2,0)-( (1.0/6.0*exp(betatilde(0,0))*sumlinpred(0,0)) + (2.0/3.0*exp(betatilde(1,0))*sumlinpred(1,0)) + (1.0/6.0*exp(betatilde(2,0))*sumlinpred(2,0)) );
  score(3,0) = Xdelta(3,0)-( (1.0/6.0*exp(betatilde(1,0))*sumlinpred(1,0)) + (2.0/3.0*exp(betatilde(2,0))*sumlinpred(2,0)) + (1.0/6.0*exp(betatilde(3,0))*sumlinpred(3,0)) );
  score(4,0) = Xdelta(4,0)-( (1.0/6.0*exp(betatilde(2,0))*sumlinpred(2,0)) + (2.0/3.0*exp(betatilde(3,0))*sumlinpred(3,0)) + (1.0/6.0*exp(betatilde(4,0))*sumlinpred(4,0)) );
  score(5,0) = Xdelta(5,0)-( (1.0/6.0*exp(betatilde(3,0))*sumlinpred(3,0)) + (2.0/3.0*exp(betatilde(4,0))*sumlinpred(4,0)) + (1.0/6.0*exp(betatilde(5,0))*sumlinpred(5,0)) );
  score(6,0) = Xdelta(6,0)-( (1.0/6.0*exp(betatilde(4,0))*sumlinpred(4,0)) + (2.0/3.0*exp(betatilde(5,0))*sumlinpred(5,0)) + (1.0/6.0*exp(betatilde(6,0))*sumlinpred(6,0)) );
  score(7,0) = Xdelta(7,0)-( (1.0/6.0*exp(betatilde(5,0))*sumlinpred(5,0)) + (2.0/3.0*exp(betatilde(6,0))*sumlinpred(6,0)) + (1.0/6.0*exp(betatilde(7,0))*sumlinpred(7,0)) );
  score(8,0) = Xdelta(8,0)-( (1.0/6.0*exp(betatilde(6,0))*sumlinpred(6,0)) + (2.0/3.0*exp(betatilde(7,0))*sumlinpred(7,0)) + (1.0/6.0*exp(betatilde(8,0))*sumlinpred(8,0)) );
  score(9,0) = Xdelta(9,0)-( (1.0/6.0*exp(betatilde(7,0))*sumlinpred(7,0)) + (2.0/3.0*exp(betatilde(8,0))*sumlinpred(8,0)) + (1.0/6.0*exp(betatilde(9,0))*sumlinpred(9,0)) );
  score(10,0) = Xdelta(10,0)-( (1.0/6.0*exp(betatilde(8,0))*sumlinpred(8,0)) + (2.0/3.0*exp(betatilde(9,0))*sumlinpred(9,0)) + (1.0/6.0*exp(betatilde(10,0))*sumlinpred(10,0)) );
  score(11,0) = Xdelta(11,0)-( (1.0/6.0*exp(betatilde(9,0))*sumlinpred(9,0)) + (2.0/3.0*exp(betatilde(10,0))*sumlinpred(10,0)) + (1.0/6.0*exp(betatilde(11,0))*sumlinpred(11,0)) );
  score(12,0) = Xdelta(12,0)-( (1.0/6.0*exp(betatilde(10,0))*sumlinpred(10,0)) + (2.0/3.0*exp(betatilde(11,0))*sumlinpred(11,0)) + (1.0/6.0*exp(betatilde(12,0))*sumlinpred(12,0)) );
  score(13,0) = Xdelta(13,0)-( (1.0/6.0*exp(betatilde(11,0))*sumlinpred(11,0)) + (2.0/3.0*exp(betatilde(12,0))*sumlinpred(12,0)) + (1.0/6.0*exp(betatilde(13,0))*sumlinpred(13,0)) );
  score(14,0) = Xdelta(14,0)-( (1.0/6.0*exp(betatilde(12,0))*sumlinpred(12,0)) + (2.0/3.0*exp(betatilde(13,0))*sumlinpred(13,0)) + (1.0/6.0*exp(betatilde(14,0))*sumlinpred(14,0)) );
  score(15,0) = Xdelta(15,0)-( (1.0/6.0*exp(betatilde(13,0))*sumlinpred(13,0)) + (2.0/3.0*exp(betatilde(14,0))*sumlinpred(14,0)) + (1.0/6.0*exp(betatilde(15,0))*sumlinpred(15,0)) );
  score(16,0) = Xdelta(16,0)-( (1.0/6.0*exp(betatilde(14,0))*sumlinpred(14,0)) + (2.0/3.0*exp(betatilde(15,0))*sumlinpred(15,0)) + (1.0/6.0*exp(betatilde(16,0))*sumlinpred(16,0)) );
  score(17,0) = Xdelta(17,0)-( (1.0/6.0*exp(betatilde(15,0))*sumlinpred(15,0)) + (2.0/3.0*exp(betatilde(16,0))*sumlinpred(16,0)) + (1.0/6.0*exp(betatilde(17,0))*sumlinpred(17,0)) );
  score(18,0) = Xdelta(18,0)-( (1.0/6.0*exp(betatilde(16,0))*sumlinpred(16,0)) + (2.0/3.0*exp(betatilde(17,0))*sumlinpred(17,0)) + (1.0/6.0*exp(betatilde(18,0))*sumlinpred(18,0)) );
  score(19,0) = Xdelta(19,0)-( (1.0/6.0*exp(betatilde(17,0))*sumlinpred(17,0)) + (2.0/3.0*exp(betatilde(18,0))*sumlinpred(18,0)) + (1.0/6.0*exp(betatilde(19,0))*sumlinpred(19,0)) );
  score(20,0) = Xdelta(20,0)-( (1.0/6.0*exp(betatilde(18,0))*sumlinpred(18,0)) + (2.0/3.0*exp(betatilde(19,0))*sumlinpred(19,0)) );
  score(21,0) = Xdelta(21,0)-( 1.0/6.0*exp(betatilde(19,0))*sumlinpred(19,0));

  for(j=0;j<betatilde.rows();j++)
    {
    Wbase(j,0) = sumlinpred(j,0)*exp(betatilde(j,0));
     Wtest(j,0) = sumlinpred(j,0)*exp(betatilde(j,0));
    }

  AWA_test(0,0) = (score(0,0)-Xdelta(0,0))/6.0;
  AWA_test(1,0) = (score(0,0)-Xdelta(0,0))*2.0/3.0;
  AWA_test(2,0) = AWA_test(0,0);
  AWA_test(0,1) = AWA_test(1,0);
  AWA_test(1,1) = 4.0/9.0*exp(betatilde(0,0))*sumlinpred(0,0) + 1.0/36.0*exp(betatilde(1,0))*sumlinpred(1,0);
  AWA_test(1,2) = 2.0/18.0*exp(betatilde(0,0))*sumlinpred(0,0) + 2.0/18.0*exp(betatilde(1,0))*sumlinpred(1,0);
  AWA_test(2,1) = AWA_test(1,2);
  AWA_test(3,1) = 1.0/36.0*exp(betatilde(1,0))*sumlinpred(1,0);
  AWA_test(2,2) = 1.0/36.0*exp(betatilde(0,0))*sumlinpred(0,0) + 4.0/9.0*exp(betatilde(1,0))*sumlinpred(1,0) + 1.0/36.0*exp(betatilde(2,0))*sumlinpred(2,0);


  ofstream sco("d:\\temp\\score.txt");
  sco<<"score bzw. muy"<<endl;
  score.prettyPrint(sco);
  sco<<"sumlinpred exact"<<endl;
  sco<<endl;
  sumlinpred.prettyPrint(sco);
  sco<<"Wbase Test"<<endl;
  Wtest.prettyPrint(sco);
  sco<<"AWA_test"<<endl;
  AWA_test.prettyPrint(sco);
  sco.close();


  }

//------------ ENDE: für baseline ----------------------------------------------

} // end: namespace MCMC


//---------------------------------------------------------------------------
#pragma package(smart_init)


