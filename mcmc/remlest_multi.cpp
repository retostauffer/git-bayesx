#include <remlest_multi.h>

#if defined(BORLAND_OUTPUT_WINDOW)
#include "StatResults.h"
#include "statwinframe.h"

#endif

//------------------------------------------------------------------------------
//------------------------ CLASS: remlest_multinomial --------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//----------------------------- Constructor ------------------------------------
//------------------------------------------------------------------------------

remlest_multinomial::remlest_multinomial(
#if defined(JAVA_OUTPUT_WINDOW)
administrator_basic * adb,
#endif
vector<MCMC::FULLCOND*> & fc,datamatrix & re,
                const ST::string & family, const ST::string & ofile,
                const int & maxiter, const double & lowerlimit,
                const double & epsi, const datamatrix & categories, ostream * lo)
  {

  nrcat2=categories.rows();
  nrcat=nrcat2+1;
  cats=categories;

  nrobs=re.rows();

  #if defined(JAVA_OUTPUT_WINDOW)
  adminb_p = adb;
  #endif

  logout = lo;
  respfamily=family;
  outfile=ofile;

  maxit=maxiter;
  lowerlim=lowerlimit;
  eps=epsi;

  fullcond = fc;
  unsigned i,j;

  xcut.push_back(0);
  xcutbeta.push_back(0);
  zcut.push_back(0);
  zcutbeta.push_back(0);

  for(i=0;i<fullcond.size();i++)
    {
    xcut.push_back(xcut[i]+fullcond[i]->get_dimX());
    if (i>0)
      {
      zcut.push_back(zcut[i-1]+fullcond[i]->get_dimZ());
      }
    }

  for(j=0; j<nrcat2; j++)
    {
    for(i=0;i<fullcond.size();i++)
      {
      xcutbeta.push_back(xcutbeta[xcutbeta.size()-1]+fullcond[i]->get_dimX());
      if (i>0)
        {
        zcutbeta.push_back(zcutbeta[zcutbeta.size()-1]+fullcond[i]->get_dimZ());
        }
      }
    }


  X = datamatrix(re.rows(),xcut[xcut.size()-1],0);
  Z = datamatrix(re.rows(),zcut[zcut.size()-1],0);

  fullcond[0]->createreml(X,Z,xcut[0],0);

  for(i=1;i<fullcond.size();i++)
    {
    fullcond[i]->createreml(X,Z,xcut[i],zcut[i-1]);
    }

  partialnrpar=X.cols()+Z.cols();
  partialnrfixed=X.cols();
  totalnrfixed=X.cols()*nrcat2;
  partialnrrandom=Z.cols();
  totalnrpar=partialnrpar*nrcat2;
  partialvar=zcut.size()-1;

  beta=statmatrix<double>(partialnrpar*(nrcat2),1,0);

  theta=statmatrix<double>(zcutbeta.size()-1,1,0);

  for(j=0; j<nrcat2; j++)
    {
    for(i=1; i<fullcond.size(); i++)
      {
      theta(j*(fullcond.size()-1)+i-1,0) = fullcond[i]->get_startlambda();
      }
    }
  }

//------------------------------------------------------------------------------
//----------------------------- Estimation -------------------------------------
//------------------------------------------------------------------------------

bool remlest_multinomial::estimate(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight)
  {
  unsigned i, j, k, l;

  outoptions();
  out("\n");

  for(i=0;i<fullcond.size();i++)
    fullcond[i]->outoptionsreml();

  out("\n");
  out("REML ESTIMATION STARTED\n",true);
  out("\n");

  bool stop = check_pause();
  if (stop)
    return true;

  double help;
  datamatrix helpmat(nrcat2,1,0);
  datamatrix helpmat2(resp.rows(),1,0);

  // Matrix to store old versions of beta and theta
  statmatrix<double>betaold(beta.rows(),1,0);
  statmatrix<double>thetaold(theta.rows(),1,0);

  // Score-function and expected Fisher information for theta
  statmatrix<double>score(theta.rows(),1,0);
  statmatrix<double>Fisher(theta.rows(),theta.rows(),0);

  // Number of iterations
  unsigned it=1;

  // Criteria to detemine convergence
  double crit1=1;                //relative changes in regression parameters
  double crit2=1;                //relative changes in variance parameters
  bool test=true;

  vector<double>stopcrit(theta.rows(),10);
  vector<int>its(theta.rows(),0);
  vector<int>signs(theta.rows(),1);

  // Linear predictor and indicator response
  statmatrix<double>respind((nrcat2)*resp.rows(),1,0);
  statmatrix<double>eta(respind.rows(),1,0);
  compute_respind(resp,respind);

  // Working observations and weights
  statmatrix<double>worky(respind.rows(),1,0);
  statmatrix<double>workweight(respind.rows(),nrcat2,0);
  statmatrix<double>mu(respind.rows(),1,0);

  // Matrix containing the inverse covariance matrix of the random effects
  statmatrix<double>Qinv(Z.cols()*nrcat2,1,0);

  // Matrices for Fisher scoring (regression parameters)
  statmatrix<double>H(beta.rows(),beta.rows(),0);
  statmatrix<double>H1(beta.rows(),1,0);

  // Matrices for Fisher scoring (variance parameters)
  statmatrix<double>Hinv(beta.rows(),beta.rows(),0);
  statmatrix<double>wresid(nrcat2,resp.rows(),0);                  //matrix containing row vectors !!

  // Transform smoothing paramater starting values to variances
  for(i=0; i<theta.rows(); i++)
    {
    theta(i,0)=1/theta(i,0);
    }

  while(test==true)
    {

    // store current values in betaold and thetaold and compute Qinv
    betaold=beta;
    thetaold=theta;
    for(i=0, l=0; i<theta.rows(); i++)
      {
      for(k=zcutbeta[i]; k<zcutbeta[i+1]; k++, l++)
        {
        Qinv(l,0)=1/theta(i,0);
        }
      }

    compute_eta2(eta);

    compute_weights(mu,workweight,worky,eta,respind);
    compute_sscp2(H,workweight);

    stop = check_pause();
    if (stop)
      return true;

    compute_sscp_resp2(H1,workweight,worky);
    H.addtodiag(Qinv,totalnrfixed,totalnrpar);

    // Fisher-Scoring for beta
    beta=H.solve(H1);

    // update linear predictor and compute residuals
    compute_eta2(eta);
    worky = worky - eta;
    for(i=0; i<resp.rows(); i++)
      {
      for(j=0; j<nrcat2; j++)
        {
        wresid(j,i)=(workweight.getRow(i*nrcat2+j)*worky.getRowBlock(i*nrcat2,(i+1)*nrcat2))(0,0);
        }
      }

    // transform theta
    for(i=0; i<theta.rows(); i++)
      {
      thetaold(i,0)=signs[i]*sqrt(thetaold(i,0));
      }

    Hinv=H.inverse();
    H.subfromdiag(Qinv,totalnrfixed,totalnrpar);

    stop = check_pause();
    if (stop)
      return true;

    // compute score-function and expected fisher information

    for(i=0; i<theta.rows(); i++)
      {
      score(i,0)=-0.5*(H.getBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i+1],totalnrfixed+zcutbeta[i+1])*thetaold(i,0)*2).trace()+
                 0.5*((H.getRowBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i+1]))*Hinv*(H.getColBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i+1]))*thetaold(i,0)*2).trace();
      for(k=0; k<theta.rows(); k++)
        {
        Fisher(i,k)=0.5*(H.getBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[k],totalnrfixed+zcutbeta[i+1],totalnrfixed+zcutbeta[k+1])*H.getBlock(totalnrfixed+zcutbeta[k],totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[k+1],totalnrfixed+zcutbeta[i+1])*thetaold(i,0)*4*thetaold(k,0)).trace()-
                    (H.getRowBlock(totalnrfixed+zcutbeta[k],totalnrfixed+zcutbeta[k+1])*Hinv*H.getColBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i+1])*H.getBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[k],totalnrfixed+zcutbeta[i+1],totalnrfixed+zcutbeta[k+1])*thetaold(i,0)*4*thetaold(k,0)).trace()+
                    0.5*(H.getRowBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i+1])*Hinv*H.getColBlock(totalnrfixed+zcutbeta[k],totalnrfixed+zcutbeta[k+1])*H.getRowBlock(totalnrfixed+zcutbeta[k],totalnrfixed+zcutbeta[k+1])*Hinv*H.getColBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i+1])*thetaold(i,0)*4*thetaold(k,0)).trace();
        Fisher(k,i)=Fisher(i,k);
        }
      }

    for(j=0; j<nrcat2; j++)
      {
      for(i=0; i<partialvar; i++)
        {
        for(l=zcut[i]; l<zcut[i+1]; l++)
          {
          help = (wresid.getRow(j)*(Z.getCol(l)))(0,0);
          score(j*partialvar+i,0) += 0.5*help*help*thetaold(j*partialvar+i,0)*2;
          }
        }
      }

    // fisher scoring for theta
    theta = thetaold + Fisher.solve(score);

    // transform theta back to original parameterisation

    for(i=0; i<theta.rows(); i++)
      {
      signs[i] = -1*(theta(i,0)<0)+1*(theta(i,0)>=0);
      theta(i,0) *= theta(i,0);
      thetaold(i,0) *= thetaold(i,0);
      }

    // test whether to stop estimation of theta[i]

    //compute norm of eta for the different catetgories
    helpmat=datamatrix(nrcat2,1,0);
    for(i=0; i<resp.rows(); i++)
      {
      for(j=0; j<nrcat2; j++)
        {
        helpmat(j,0) += eta(i*nrcat2+j,0)*eta(i*nrcat2+j,0);
        }
      }
    for(j=0; j<nrcat2; j++)
      {
      helpmat(j,0) = sqrt(helpmat(j,0));
      }
    // compute norm of the random parts
    for(j=0; j<nrcat2; j++)
      {
      for(i=0; i<partialvar; i++)
        {
        helpmat2 = Z.getColBlock(zcut[i],zcut[i+1])*beta.getRowBlock(totalnrfixed+zcutbeta[j*partialvar+i],totalnrfixed+zcutbeta[j*partialvar+i+1]);
        stopcrit[j*partialvar+i]=helpmat2.norm(0)/helpmat(j,0);
        if(stopcrit[j*partialvar+i]<lowerlim)
          {
          theta(j*partialvar+i,0)=thetaold(j*partialvar+i,0);
          }
        else
          {
          its[j*partialvar+i]=it;
          }
        }
      }

    // compute convergence criteria
    help=betaold.norm(0);
    if(help==0)
      {
      help=0.00001;
      }
    betaold.minus(betaold,beta);
    crit1 = betaold.norm(0)/help;

    help=thetaold.norm(0);
    if(help==0)
      {
      help=0.00001;
      }
    thetaold.minus(thetaold,theta);
    crit2 = thetaold.norm(0)/help;

    stop = check_pause();
    if (stop)
      return true;

    out("  iteration "+ST::inttostring(it)+"\n");
    out("  relative changes in the regression coefficients: "+
         ST::doubletostring(crit1,6)+"\n");
    out("  relative changes in the variance parameters:     "+
         ST::doubletostring(crit2,6)+"\n");
    out("\n");

    // test criterion
    test=((crit1>eps) || (crit2>eps)) && (it<(unsigned)maxit);

    // count iteration
    it=it+1;
    }

  if(it<(unsigned)maxit)
    {
    out("\n");
    out("REML ESTIMATION CONVERGED\n",true);
    out("\n");
    }
  else
    {
    out("\n");
    outerror("WARNING: Number of iterations reached " + ST::inttostring(maxit) + "\n");
    outerror("         REML ESTIMATION DID NOT CONVERGE\n");
    out("\n");
    }
  out("ESTIMATION RESULTS:\n",true);
  out("\n");

  ofstream outit((outfile+"_it.raw").strtochar());
  outit << it-1;
  outit.close();

  datamatrix thetareml(theta.rows(),3,0);
  thetareml.putCol(0,theta);
  for(i=0; i<theta.rows(); i++)
    {
    if(stopcrit[i]<lowerlim)
      {
      thetareml(i,1)=1;
      }
    thetareml(i,2)=its[i];
    }

  for(j=0; j<nrcat2; j++)
    {
    out("\n");
    out("RESULTS FOR CATEGORY "+ST::doubletostring(cats(j,0),6)+":\n",true);
    out("\n");
    for(i=1; i<fullcond.size(); i++)
      {
      beta(j*partialnrfixed,0) += fullcond[i]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[i],zcut[i-1],j*partialvar+i-1,false,xcutbeta[j*fullcond.size()+i],totalnrfixed+zcutbeta[j*partialvar+i-1],cats(j,0),true,j*fullcond.size()+i);
      }
    beta(j*partialnrfixed,0) += fullcond[0]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[0],0,0,false,xcutbeta[j*fullcond.size()],0,cats(j,0),true,0);
    }

  double loglike=0;
  double aic=0;
  double bic=0;
  double gcv=0;
  double df=(H*Hinv).trace();
  double refprob;

  for(i=0; i<resp.rows(); i++)
    {
    k=0;
    refprob=0;
    for(j=0; j<nrcat2; j++)
      {
      if(respind(i*nrcat2+j,0)==1)
        {
        loglike += log(mu(i*nrcat2+j,0));
        k=1;
        }
      else
        {
        refprob += mu(i*nrcat2+j,0);
        }
      }
    if(k==0)
      {
      loglike += log(1-refprob);
      }
    }
  loglike *= -2;
  gcv = loglike/(double)resp.rows()*(1-(double)df/(double)eta.rows())*(1-(double)df/(double)eta.rows());
  aic = loglike + 2*df;
  bic = loglike + log(resp.rows())*df;

  out("\n");
  out("  Model Fit\n",true);
  out("\n");
  out("\n");
  out("  -2*log-likelihood:                 " + ST::doubletostring(loglike,6) + "\n");
  out("  Degrees of freedom:                " + ST::doubletostring(df,6) + "\n");
  out("  (conditional) AIC:                 " + ST::doubletostring(aic,6) + "\n");
  out("  (conditional) BIC:                 " + ST::doubletostring(bic,6) + "\n");
  out("  GCV (based on deviance residuals): " + ST::doubletostring(gcv,6) + "\n");
  out("\n");

  out("\n");
  out("  Additive predictors and expectations\n",true);
  out("\n");
  out("\n");
  out("  Additive predictors and expectations for each observation\n");
  out("  and category are stored in file\n");
  out("  "+outfile+"_predict.raw\n");
  out("\n");
  out("\n");

  ofstream outpredict((outfile+"_predict.raw").strtochar());
  for(j=0; j< nrcat2; j++)
    {
    outpredict << "eta" << cats(j,0) << " ";
    outpredict << "mu" << cats(j,0) << " ";
    }
  outpredict << endl;
  for(i=0; i<nrobs; i++)
    {
    for(j=0; j<nrcat2; j++)
      {
      outpredict << eta(i*nrcat2+j,0) << " " << mu(i*nrcat2+j,0) << " ";
      }
    outpredict << endl;
    }
  outpredict.close();

  return false;
  }

bool remlest_multinomial::estimate_glm(const datamatrix resp,
                  const datamatrix & offset, const datamatrix & weight)
  {
  unsigned i,j;

  outoptions();
  out("\n");

  for(i=0;i<fullcond.size();i++)
    fullcond[i]->outoptionsreml();

  out("\n");
  out("REML ESTIMATION STARTED\n",true);
  out("\n");

  bool stop = check_pause();
  if (stop)
    return true;

  double help;

  // Matrix to store old version of beta
  statmatrix<double>betaold(beta.rows(),1,0);

  // Number of iterations
  unsigned it=1;

  // Criteria to detemine convergence
  double crit1=1;                //relative changes in regression parameters
  bool test=true;

  // Linear predictor and indicator response
  statmatrix<double>respind((nrcat2)*resp.rows(),1,0);
  statmatrix<double>eta(respind.rows(),1,0);
  compute_respind(resp,respind);

  // Working observations and weights
  statmatrix<double>worky(respind.rows(),1,0);
  statmatrix<double>workweight(respind.rows(),nrcat2,0);
  statmatrix<double>mu(respind.rows(),1,0);

  // Matrices for Fisher scoring (regression parameters)
  statmatrix<double>H(beta.rows(),beta.rows(),0);
  statmatrix<double>H1(beta.rows(),1,0);

  // Estimation loop
  while(test==true)
    {
    // store current values in betaold
    betaold=beta;

    compute_eta(eta);

    compute_weights(mu,workweight,worky,eta,respind);
    compute_sscp(H,workweight);
    compute_sscp_resp(H1,workweight,worky);

    stop = check_pause();
    if (stop)
      return true;

    // Fisher-Scoring for beta
    beta=H.solve(H1);

    // compute convergence criteria
    help=betaold.norm(0);
    if(help==0)
      {
      help=0.00001;
      }
    betaold.minus(betaold,beta);
    crit1 = betaold.norm(0)/help;

    stop = check_pause();
    if (stop)
      return true;

    out("  iteration "+ST::inttostring(it)+"\n");
    out("  relative changes in the regression coefficients: "+
         ST::doubletostring(crit1,6)+"\n");
    out("\n");

    // test criterion
    test=(crit1>eps) && (it<(unsigned)maxit);

    // count iteration
    it=it+1;

    }

  if(it<(unsigned)maxit)
    {
    out("\n");
    out("REML ESTIMATION CONVERGED\n",true);
    out("\n");
    }
  else
    {
    out("\n");
    outerror("WARNING: Number of iterations reached " + ST::inttostring(maxit) + "\n");
    outerror("         REML ESTIMATION DID NOT CONVERGE\n");
    out("\n");
    }
  out("ESTIMATION RESULTS:\n",true);
  out("\n");

  H=H.inverse();

  for(j=0; j<nrcat2; j++)
    {
    out("\n");
    out("RESULTS FOR CATEGORY "+ST::doubletostring(cats(j,0),6)+":\n",true);
    out("\n");
    for(i=1; i<fullcond.size(); i++)
      {
      beta(j*partialnrfixed,0) += fullcond[i]->outresultsreml(X,Z,beta,H,datamatrix(1,1,0),xcut[i],zcut[i-1],j*partialvar+i-1,false,xcutbeta[j*fullcond.size()+i],totalnrfixed+zcutbeta[j*partialvar+i-1],cats(j,0),true,j*fullcond.size()+i);
      }
    beta(j*partialnrfixed,0) += fullcond[0]->outresultsreml(X,Z,beta,H,datamatrix(1,1,0),xcut[0],0,0,false,xcutbeta[j*fullcond.size()],0,cats(j,0),true,0);
    }

  double loglike=0;
  double aic=0;
  double bic=0;
  double gcv=0;
  double df=beta.rows();
  double refprob;
  unsigned k;

  for(i=0; i<resp.rows(); i++)
    {
    k=0;
    refprob=0;
    for(j=0; j<nrcat2; j++)
      {
      if(respind(i*nrcat2+j,0)==1)
        {
        loglike += log(mu(i*nrcat2+j,0));
        k=1;
        }
      else
        {
        refprob += mu(i*nrcat2+j,0);
        }
      }
    if(k==0)
      {
      loglike += log(1-refprob);
      }
    }
  loglike *= -2;
  gcv = loglike/(double)resp.rows()*(1-(double)df/(double)eta.rows())*(1-(double)df/(double)eta.rows());
  aic = loglike + 2*df;
  bic = loglike + log(resp.rows())*df;

  out("\n");
  out("  Model Fit\n",true);
  out("\n");
  out("\n");
  out("  -2*log-likelihood:                 " + ST::doubletostring(loglike,6) + "\n");
  out("  Degrees of freedom:                " + ST::doubletostring(df,6) + "\n");
  out("  (conditional) AIC:                 " + ST::doubletostring(aic,6) + "\n");
  out("  (conditional) BIC:                 " + ST::doubletostring(bic,6) + "\n");
  out("  GCV (based on deviance residuals): " + ST::doubletostring(gcv,6) + "\n");
  out("\n");


  out("\n");
  out("  Linear predictors and expectations\n",true);
  out("\n");
  out("\n");
  out("  Linear predictors and expectations for each observation\n");
  out("  and category are stored in file\n");
  out("  "+outfile+"_predict.raw\n");
  out("\n");
  out("\n");

  ofstream outpredict((outfile+"_predict.raw").strtochar());
  for(j=0; j< nrcat2; j++)
    {
    outpredict << "eta" << cats(j,0) << " ";
    outpredict << "mu" << cats(j,0) << " ";
    }
  outpredict << endl;
  for(i=0; i<nrobs; i++)
    {
    for(j=0; j<nrcat2; j++)
      {
      outpredict << eta(i*nrcat2+j,0) << " " << mu(i*nrcat2+j,0) << " ";
      }
    outpredict << endl;
    }
  outpredict.close();

  return false;
  }

//------------------------------------------------------------------------------
//------------- Weights, expectation, linear predictor, etc --------------------
//------------------------------------------------------------------------------

void remlest_multinomial::compute_respind(const datamatrix & re, datamatrix & respind)
  {
  unsigned i,j;
  for(i=0; i<nrobs; i++)
    {
    for(j=0; j<nrcat2; j++)
      {
      if(re(i,0)==cats(j,0))
        {
        respind(i*nrcat2+j,0)=1;
        }
      }
    }
  }

void remlest_multinomial::compute_eta(datamatrix & eta)
  {
  unsigned i,j;
  for(i=0; i<nrobs; i++)
    {
    for(j=0; j<nrcat2; j++)
      {
      eta(i*nrcat2+j,0)=((X.getRow(i))*(beta.getRowBlock(j*partialnrpar,(j+1)*partialnrpar)))(0,0);
      }
    }
  }

void remlest_multinomial::compute_eta2(datamatrix & eta)
  {
  unsigned i,j;
  for(i=0; i<nrobs; i++)
    {
    for(j=0; j<nrcat2; j++)
      {
      eta(i*nrcat2+j,0)=((X.getRow(i))*(beta.getRowBlock(j*partialnrfixed,(j+1)*partialnrfixed)))(0,0)+
                        ((Z.getRow(i))*(beta.getRowBlock(totalnrfixed+j*partialnrrandom,totalnrfixed+(j+1)*partialnrrandom)))(0,0);
      }
    }
  }

void remlest_multinomial::compute_weights(datamatrix & mu, datamatrix & weights,
                  datamatrix & worky, datamatrix & eta, datamatrix & respind)
  {
  unsigned i,j,k1,k2,l;

// Compute mu
  datamatrix expos(nrcat2,1,0);
  double exposum;
  for(i=0; i<nrobs; i++)
    {
    exposum=0;
    for(j=0; j<nrcat2; j++)
      {
      expos(j,0)=exp(eta(i*nrcat2+j,0));
      exposum+=expos(j,0);
      }
    exposum+=1;
    for(j=0; j<nrcat2; j++)
      {
      mu(i*nrcat2+j,0)=expos(j,0)/exposum;
      }
    }

// Compute weights

  for(i=0; i<nrobs; i++)
    {
    for(j=i*nrcat2,l=0; l<nrcat2; j++, l++)
      {
      weights(j,l)=mu(j,0)*(1-mu(j,0));
      for(k1=j+1,k2=l+1; k2<nrcat2; k1++, k2++)
        {
        weights(j,k2)=-mu(j,0)*mu(k1,0);
        weights(k1,l)=weights(j,k2);
        }
      }
    }

// Compute worky;

  for(i=0; i<nrobs; i++)
    {
    worky.putRowBlock(i*nrcat2,(i+1)*nrcat2,eta.getRowBlock(i*nrcat2,(i+1)*nrcat2)+
                      weights.getRowBlock(i*nrcat2,(i+1)*nrcat2).inverse()*
                      (respind.getRowBlock(i*nrcat2,(i+1)*nrcat2)-mu.getRowBlock(i*nrcat2,(i+1)*nrcat2)));
    }

  }

void remlest_multinomial::compute_sscp(datamatrix & H, datamatrix & workweight)
  {
  unsigned i;
  H=datamatrix(H.rows(),H.cols(),0);
  datamatrix Htemp=datamatrix(H.rows(),H.cols(),0);
  for(i=0; i<nrobs; i++)
    {
    Htemp=kronecker(workweight.getRowBlock(i*nrcat2,(i+1)*nrcat2),(X.getRow(i).transposed())*(X.getRow(i)));
    H.plus(Htemp);
    }
  }

void remlest_multinomial::compute_sscp2(datamatrix & H, datamatrix & workweight)
  {
  unsigned i;
  H=datamatrix(H.rows(),H.cols(),0);
  datamatrix Htemp=datamatrix(H.rows(),H.cols(),0);
  datamatrix weighttemp=datamatrix(nrcat2,nrcat2,0);
  for(i=0; i<nrobs; i++)
    {
    weighttemp=workweight.getRowBlock(i*nrcat2,(i+1)*nrcat2);
    Htemp.putBlock(kronecker(weighttemp,(X.getRow(i).transposed())*(X.getRow(i))),0,0,totalnrfixed,totalnrfixed);
    Htemp.putBlock(kronecker(weighttemp,(Z.getRow(i).transposed())*(Z.getRow(i))),totalnrfixed,totalnrfixed,totalnrpar,totalnrpar);
    Htemp.putBlock(kronecker(weighttemp,(X.getRow(i).transposed())*(Z.getRow(i))),0,totalnrfixed,totalnrfixed,totalnrpar);
    Htemp.putBlock(Htemp.getBlock(0,totalnrfixed,totalnrfixed,totalnrpar).transposed(),totalnrfixed,0,totalnrpar,totalnrfixed);
    H.plus(Htemp);
    }
  }

void remlest_multinomial::compute_sscp_resp(datamatrix & H1, datamatrix & workweight, datamatrix & worky)
  {
  unsigned i;
  H1=datamatrix(H1.rows(),1,0);
  datamatrix H1temp = datamatrix(nrcat2,1,0);
  datamatrix xtemp;
  for(i=0; i<nrobs;i++)
    {
    H1temp = workweight.getRowBlock(i*nrcat2,(i+1)*nrcat2)*worky.getRowBlock(i*nrcat2,(i+1)*nrcat2);
    xtemp = X.getRow(i).transposed();
    H1.plus(kronecker(H1temp,xtemp));
    }
  }

void remlest_multinomial::compute_sscp_resp2(datamatrix & H1, datamatrix & workweight, datamatrix & worky)
  {
  unsigned i;
  H1=datamatrix(H1.rows(),1,0);
  datamatrix H1temp = datamatrix(H1.rows(),1,0);
  datamatrix weightytemp=datamatrix(nrcat2,1,0);
  for(i=0; i<nrobs;i++)
    {
    weightytemp=workweight.getRowBlock(i*nrcat2,(i+1)*nrcat2)*worky.getRowBlock(i*nrcat2,(i+1)*nrcat2);
    H1temp.putRowBlock(0,totalnrfixed,kronecker(weightytemp,X.getRow(i).transposed()));
    H1temp.putRowBlock(totalnrfixed,totalnrpar,kronecker(weightytemp,Z.getRow(i).transposed()));
    H1.plus(H1temp);
    }
  }

//------------------------------------------------------------------------------
//----------------------------- Object description -----------------------------
//------------------------------------------------------------------------------

  void remlest_multinomial::outoptions()
    {
    out("\n");
    out("GENERAL OPTIONS:\n",true);
    out("\n");
    out("  Maxmimum number of iterations:          "+ST::inttostring(maxit)+"\n");
    out("  Termination criterion:                  "+ST::doubletostring(eps,7)+"\n");
    out("  Stopping criterion for small variances: "+ST::doubletostring(lowerlim,6)+"\n");
    out("\n");
    out("RESPONSE DISTRIBUTION:\n",true);
    out("\n");
    ST::string familyname = "multinomial logit";
    out("  Family:                 "+familyname+"\n");
    out("  Number of observations: "+ST::inttostring(X.rows())+"\n");
    }

//------------------------------------------------------------------------------
//----------------------------- Writing results --------------------------------
//------------------------------------------------------------------------------

void remlest_multinomial::make_plots(ofstream & outtex,ST::string path_batch,
                         ST::string path_splus)
  {

  char hcharu = '_';
  ST::string hstringu = "\\_";

  unsigned i,j;
  ST::string pathresult;
  bool stil = false;

// Schleife überprüft, ob es ein fullcond-Object
// gibt, bei dem Effekt gezeichnet werden kann
  MCMC::plotstyles plst;
  for(j=0;j<fullcond.size();j++)
    {
    plst = fullcond[j]->get_plotstyle();
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

#if defined(JAVA_OUTPUT_WINDOW)
    out("  --------------------------------------------------------------------------- \n");
    out("\n");
    out("  Batch file for visualizing effects of nonlinear functions is stored in file \n");
    out("  " + path_batch + "\n");
    out("\n");
#endif

    bool stil2 = true;
    for(j=0;j<fullcond.size();j++)  //Schleife überprüft, ob es map-Objekt gibt
      {
      plst = fullcond[j]->get_plotstyle();
      if(plst == MCMC::drawmap)
        stil2 = false;
      }

    if(stil2 == true)
      {
      out("  --------------------------------------------------------------------------- \n");
      out("\n");
      out("  Batch file for visualizing effects of nonlinear functions ");
      out("  in S-Plus is stored in file \n");
      out("  " + path_splus + "\n");
      out("\n");
      }

    if(stil2 == false)
      {
      out("  NOTE: 'input filename' must be substituted by the filename of the boundary-file \n");
      out("\n");
      out("  --------------------------------------------------------------------------- \n");
      out("\n");
      out("  Batch file for visualizing effects of nonlinear functions ");
      out("  in S-Plus is stored in file \n");
      out("  " + path_splus + "\n");
      out("\n");
      out("  NOTE: 'input filename' must be substituted by the filename of the boundary-file \n");
      out("\n");
      }


    outbatch << "% usefile " << path_batch << endl;

    // falls andere Quantile gewünscht werden
    double u = fullcond[0]->get_level1();
    double o = fullcond[0]->get_level2();
    ST::string u_str = ST::doubletostring(u,0);
    ST::string o_str = ST::doubletostring(o,0);

    for(i=0; i<nrcat2; i++)
      {
    // durchlaufen der Fullconditionals
    for(j=0;j<fullcond.size();j++)
      {

      // Pfad der Regr.-Ergebnisse
      pathresult = fullcond[j]->get_pathresult();
      pathresult = pathresult.insert_after_string(ST::doubletostring(cats(i,0),6)+"_","_f_");

      // Plotstyle: noplot, plotnonp, drawmap
      plst = fullcond[j]->get_plotstyle();

      if (plst != MCMC::noplot)
        {

        // Pfade für ps-, tex-, SPlus-files
        ST::string pathps = pathresult.substr(0, pathresult.length()-4);
        ST::string pathgr = pathps.replaceallsigns('\\', '/');

        char hchar = '\\';
        ST::string hstring = "\\\\";

        ST::string pathps_spl = pathps.insert_string_char(hchar,hstring);
        ST::string pathres_spl = pathresult.insert_string_char(hchar,hstring);

        if (plst == MCMC::plotnonp)
          {
          outbatch << "\n";                // Befehle f. d. batch-file
          outbatch << "dataset _dat" << endl;
          outbatch << "_dat.infile using " << pathresult << endl;
          outbatch << "graph _g" << endl;
          vector<ST::string> varnames = fullcond[j]->get_datanames();
          ST::string xvar = varnames[0];
          outbatch << "_g.plot " << xvar
                   << " pmode ci" << u_str << "lower ci"
                   << o_str.replaceallsigns('.','p') << "lower ci"
                   << o_str.replaceallsigns('.','p') << "upper ci"
                   << u_str.replaceallsigns('.','p') << "upper, "
                   << "title = \"Effect of " << xvar << "\" xlab = " << xvar
                   << " ylab = \" \" " << "outfile = " << pathps
                   << ".ps replace using _dat" << endl;
          outbatch << "drop _dat" << endl;
          outbatch << "drop _g" << endl;
          // Plot-Befehle f. d. SPlus-file
          outsplus << "plotnonp(\"" << pathres_spl << "\", psname = \""
                   << pathps_spl << ".ps\")" << endl;
          // Plot-Befehle f. d. tex-file
          ST::string effect = xvar;
          if(varnames.size()>1)
            {
            effect = varnames[1] + "*" + effect;
            }
          outtex << "\n\\begin{figure}[h!]" << endl
                  << "\\centering" << endl
                  << "\\includegraphics[scale=0.6]{" << pathgr << ".ps}" << endl
                  << "\\caption{Non--linear Effect of '" <<
                  effect.insert_string_char(hcharu,hstringu) << "'";
          outtex << " (Category " << cats(i,0) << ")." << endl << "Shown are the posterior modes together with "
                 << u_str << "\\% and " << o_str
                 << "\\% pointwise credible intervals.}" << endl
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
          vector<ST::string> varnames = fullcond[j]->get_datanames();
          ST::string regionvar = varnames[0];
          outbatch << "_g.drawmap " << "pmode" << " " << regionvar
                   << ", map = _map color outfile = " << pathps
                   << "_pmode.ps replace using _dat" << endl;
          outbatch << "_g.drawmap " << "pcat" << u_str << " " << regionvar
                     << ", map = _map nolegend pcat outfile = " << pathps
                     << "_pcat" << u_str << ".ps replace using _dat" << endl;
          outbatch << "_g.drawmap " << "pcat" << o_str << " " << regionvar
                     << ", map = _map nolegend pcat outfile = " << pathps
                     << "_pcat" << o_str << ".ps replace using _dat" << endl;
          outbatch << "drop _dat" << endl;
          outbatch << "drop _g" << endl;
          outbatch << "drop _map" << endl;
          // Plot-Befehle f. d. SPlus-file
          outsplus << "# NOTE: 'input_filename' must be substituted by the "
                   << "filename of the boundary-file \n"
                   << "# NOTE: choose a 'name' for the map \n" << endl
                   << "readbndfile(\"'input_filename'\", \"'name'\")" << endl
                   << "drawmap(map = 'name', outfile = \"" << pathps_spl
                   <<"_pmode.ps\", dfile = \"" << pathres_spl
                   << "\" ,plotvar = \"pmode\", regionvar = \""
                   << regionvar << "\", color = T)" << endl;
          outsplus << "drawmap(map = 'name', outfile = \"" << pathps_spl
                    <<"_pcat" << u_str << ".ps\", dfile = \"" << pathres_spl
                    << "\" ,plotvar = \"pcat" << u_str << "\", regionvar = \""
                  << regionvar << "\", legend = F, pcat = T)" << endl;
          outsplus << "drawmap(map = 'name', outfile = \"" << pathps_spl
                    <<"_pcat" << o_str << ".ps\", dfile = \"" << pathres_spl
                    << "\",plotvar = \"pcat" << o_str << "\", regionvar = \""
                    << regionvar << "\", legend = F, pcat = T)" << endl;
            // Plot-Befehle f. d. tex-file
          ST::string effect = regionvar;
          if(varnames.size()>1)
            {
            effect = varnames[1] + "*" + effect;
            }
          outtex << "\n\\begin{figure}[h!]" << endl
                 << "\\centering" << endl
                 << "\\includegraphics[scale=0.6]{" << pathgr << "_pmode.ps}"
                 << endl
                 << "\\caption{Non--linear Effect of '" <<
                 effect.insert_string_char(hcharu,hstringu) << "'";
          outtex << " (Category " << cats(i,0) << "). Shown are the posterior modes.}" << endl
                 << "\\end{figure}" << endl;
          outtex << "\n\\begin{figure}[htb]" << endl
                 << "\\centering" << endl
                 << "\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                 << u_str << ".ps}" << endl
                 << "\\caption{Non--linear Effect of '" << effect << "'";
          outtex << " (Category " << cats(i,0) << "). Posterior probabilities for a nominal level of "
                 << u_str << "\\%." << endl
                 << "Black denotes regions with strictly negative credible intervals,"
                 << endl
                 << "white denotes regions with strictly positive credible intervals.}"
                 << endl << "\\end{figure}" << endl;
          outtex << "\n\\begin{figure}[htb]" << endl
                 << "\\centering" << endl
                 << "\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                 << o_str << ".ps}" << endl
                 << "\\caption{Non--linear Effect of '" << effect << "'";
          outtex << " (Category " << cats(i,0) << "). Posterior probabilities for a nominal level of "
                 << o_str << "\\%." << endl
                 << "Black denotes regions with strictly negative credible intervals,"
                 << endl
                 << "white denotes regions with strictly positive credible intervals.}"
                 << endl << "\\end{figure}" << endl;
          } // end: else if
        } // end: if
      } // end: for
      }
    }
  }

void remlest_multinomial::make_model(ofstream & outtex, const ST::string & rname)
  {
  ST::string familyname;
  if(respfamily=="multinomial")
    {
    familyname="multinomial logit";
    }

  //Anz. Beob. wird übergeben
  unsigned obs = X.rows();

  char charh = '_';
  ST::string stringh = "\\_";
  ST::string helprname = rname.insert_string_char(charh,stringh);

  //schreibt das Modell und die Priori-annahmen ins Tex-File
  outtex << "\n\\noindent {\\bf \\large Response:}" << endl
         << "\\begin{tabbing}\n"
         << "Number of observations: \\= " << obs << "\\\\" << endl
         << "Response Variable: \\> " << helprname << "\\\\" << endl
         << "Family: \\> " << familyname << "\\\\" << endl
         << "\\end{tabbing}" << endl
         << "\n\\noindent {\\bf \\large Predictor:}\\\\" << endl;
  }

void remlest_multinomial::make_predictor(ofstream & outtex)
  {

  unsigned j;

  ST::string term2 = fullcond[0]->get_term_symbolic();
  ST::string term = "$\\eta$ & $=$ & $" + term2;    //linearer Prädiktor wird erweitert
  for(j=1;j<fullcond.size();j++)
    {
    out(fullcond[j]->get_results_type());
    term2 = fullcond[j]->get_term_symbolic();
    term = term + " + " + term2;    //linearer Prädiktor wird erweitert
    }

  outtex << endl << "\n\\begin{tabular}{ccp{12cm}}\n";
  for(j=0; j<nrcat2; j++)
    {
    term2 = term.insert_after_string("^{("+ST::doubletostring(cats(j,0),6)+")}","\\eta");
    term2 = term2.insert_after_all_string("^{("+ST::doubletostring(cats(j,0),6)+")}","\\gamma");
    term2 = term2.insert_after_all_string("^{("+ST::doubletostring(cats(j,0),6)+")}","+ f");

    outtex << term2 << "$\\\\\n";
    }
  outtex << "\\end{tabular}\n\\\\ \n\\\\" << endl;
  }

void remlest_multinomial::make_prior(ofstream & outtex)
  {
  unsigned i,j;
  outtex << "\n\\noindent {\\bf \\large Priors:}\\\\" << endl << "\\\\" << endl;
  for(j=0;j<fullcond.size();j++)
    {
    vector<ST::string> prior = fullcond[j]->get_priorassumptions();
    if(prior.size() != 0)// nur wenn Priors da sind (d.h. Vektor hat Elemente)
      {
      if(fullcond[j]->get_results_type()!="fixed")
        {
        prior[0] = prior[0].insert_after_string("^{(j)}","f");
        }
      for(i=0;i<prior.size();i++)
        {
        if( j!=0 || i<prior.size()-1)
          {
          outtex << prior[i] << "\\\\" << endl;
          }
        }
      outtex << "\\\\" <<endl;
      }
    }
  }

void remlest_multinomial::make_options(ofstream & outtex)
  {
  double l1 = fullcond[0]->get_level1();
  double l2 = fullcond[0]->get_level2();

  //schreibt REML options ins Tex-File
  outtex << "\n\\noindent {\\bf \\large General Options:}" << endl
         << "\\begin{tabbing}" << endl
         << "Levels for credible intervals: \\hspace{2cm}\\= \\\\" << endl
         << "Level 1: \\> " << l1 << "\\\\" << endl
         << "Level 2: \\> " << l2 << "\\\\" << endl
         << "Maxmimum number of iterations: \\> " << maxit << "\\\\" << endl
         << "Termination criterion: \\> " << eps << "\\\\" << endl
         << "Stopping criterion for small variances: \\> " << lowerlim << endl
         << "\\end{tabbing}\n"  << "\\vspace{0.5cm}" <<  endl;
  }

void remlest_multinomial::make_fixed_table(ofstream & outtex)
  {

  // falls andere Quantile gewünscht werden
  double u = fullcond[0]->get_level1();
  ST::string u_str = ST::doubletostring(u,0);

  vector<ST::string> h;
  h = fullcond[0]->get_results_latex();

  unsigned i,j, r;
  unsigned partialhsize=h.size()/nrcat2;

  for(i=0; i<nrcat2; i++)
    {
    r=2;
    // Tabelle im Tex-File mit fixen Effekten
    outtex << "\n\\newpage \n" << endl << "\n\\noindent {\\bf \\large Fixed Effects (Category "
         << ST::doubletostring(cats(i,0),6) << "):}\\\\"
         << endl << "\\\\" << endl;

    outtex << "\\begin{tabular}{|r|rrrrr|}" << endl << "\\hline" << endl
           << "Variable & Post. Mode & Std. Dev. & p-value & \\multicolumn{2}{r|}{" << u << "\\% confidence interval}\\\\"
           << endl << "\\hline" << endl;

    for (j=i*partialhsize;j<(i+1)*partialhsize;j++)
      {
      r++;
      if (r < 39)
        {
        outtex << h[j] << endl;
        }
      else
        {
        r=1;
        outtex << "\\hline \n\\end{tabular}" << endl;

        outtex << "\n\\newpage \n" << endl
               << "\n\\noindent {\\bf \\large Fixed Effects (Category "
               << ST::doubletostring(cats(i,0),6) << "continued):}\\\\"
               << endl << "\\\\" << endl;

        outtex << "\\begin{tabular}{|r|rrrrr|}" << endl << "\\hline" << endl
               << "Variable & Post. Mode & Std. Dev. & p-value & \\multicolumn{2}{r|}{" << u << "\\% confidence interval}\\\\"
               << endl << "\\hline" << endl;

        outtex << h[j] << endl;
        }
      }
    outtex << "\\hline \n\\end{tabular}" << endl;
    }
  }

void remlest_multinomial::make_graphics(const ST::string & title,
                     const ST::string & path_batch,
                     const ST::string & path_tex,
                     const ST::string & path_splus,
                     const ST::string & rname)
  {
  ST::string pathresult;                 //Pfad des Ergebnis-Files

  vector<ST::string> distr_results;

 // erzeugt Tex-File
  ofstream outtex(path_tex.strtochar());

  //erzeugt den Kopf des Tex-Files
  outtex << "\\documentclass[a4paper, 12pt]{article}" << endl
         << "\n" << "\\usepackage{graphicx}" << endl
         << "\\parindent0em" << endl
         << "\n\\begin{document}" << endl
         << "\\begin{center}" << endl
         << "\\LARGE{\\bf " << title << "}"
         << endl << "\\end{center} \n\\vspace{1cm}" << endl;

  make_model(outtex,rname);

  make_predictor(outtex);

  make_prior(outtex);

  make_options(outtex);

  make_fixed_table(outtex);

  // Pfade der Files
  //werden im BayesX-Output angegeben
  out("  Files of model summary: \n" , true);
  out("\n");

  make_plots(outtex,path_batch,path_splus);

  out("  --------------------------------------------------------------------------- \n");
  out("\n");
  out("  Latex file of model summaries is stored in file \n");
  out("  " + path_tex + "\n");
  out("\n");
  out("  --------------------------------------------------------------------------- \n");
  out("\n");


  outtex << "\\end{document}" << endl;

  }

bool remlest_multinomial::check_pause()
  {
#if defined(BORLAND_OUTPUT_WINDOW)
  Application->ProcessMessages();
  if (Frame->stop)
    {
    return true;
    }

  if (Frame->pause)
    {
    out("\n");
    out("ESTIMATION PAUSED\n");
    out("Click CONTINUE to proceed\n");
    out("\n");

    while (Frame->pause)
      {
      Application->ProcessMessages();
      }

    out("ESTIMATION CONTINUED\n");
    out("\n");
    }
  return false;
#elif defined(JAVA_OUTPUT_WINDOW)
  return adminb_p->breakcommand();
#endif
  }

void remlest_multinomial::out(const ST::string & s,bool thick,bool italic,
                      unsigned size,int r,int g, int b)
  {
#if defined(BORLAND_OUTPUT_WINDOW)
  ST::string sh = s;
  sh = sh.replaceallsigns('\n',' ');
  if (!Frame->suppoutput)
    Results->ResultsRichEdit->Lines->Append(sh.strtochar());
 if (!(logout->fail()))
    (*logout) << s << flush;
#elif defined(JAVA_OUTPUT_WINDOW)
  ST::string sh = s;
  sh = sh.replaceallsigns('\n',' ');
  sh = sh+"\n";
  if (!adminb_p->get_suppressoutput())
    adminb_p->Java->CallVoidMethod(adminb_p->BayesX_obj, adminb_p->javaoutput,
    adminb_p->Java->NewStringUTF(sh.strtochar()),thick,italic,size,r,g,b);
  if (!(logout->fail()))
    (*logout) << s << flush;
#else
  (*logout) << s << flush;
#endif
  }


void remlest_multinomial::outerror(const ST::string & s)
  {
  out(s,true,true,12,255,0,0);
  }


