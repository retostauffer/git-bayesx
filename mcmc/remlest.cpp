#include <remlest.h>

#if defined(BORLAND_OUTPUT_WINDOW)
#include "StatResults.h"
#include "statwinframe.h"

#endif

//------------------------------------------------------------------------------
//----------------------------- Constructor ------------------------------------
//------------------------------------------------------------------------------

  remlest::remlest(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adb,
  #endif
  vector<MCMC::FULLCOND*> & fc,datamatrix & re,bool dispers,
                   const ST::string & family, const ST::string & ofile,
                  const int & maxiter, const double & lowerlimit,
                  const double & epsi, ostream * lo)
    {

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
    unsigned i;

    xcut.push_back(0);
    zcut.push_back(0);

    for(i=0;i<fullcond.size();i++)
      {
      xcut.push_back(xcut[i]+fullcond[i]->get_dimX());
      if (i>0)
        {
        zcut.push_back(zcut[i-1]+fullcond[i]->get_dimZ());
        }
      }

    X = datamatrix(re.rows(),xcut[xcut.size()-1],0);

    Z = datamatrix(re.rows(),zcut[zcut.size()-1],0);

    fullcond[0]->createreml(X,Z,xcut[0],0);

    for(i=1;i<fullcond.size();i++)
      {
      fullcond[i]->createreml(X,Z,xcut[i],zcut[i-1]);
      }

    beta=statmatrix<double>(X.cols()+Z.cols(),1,0);
    if(dispers==true)
      {
      theta=statmatrix<double>(zcut.size(),1,0);
      }
    else
      {
      theta=statmatrix<double>(zcut.size()-1,1,0);
      }
    for(i=1; i<fullcond.size(); i++)
      {
      theta(i-1,0) = fullcond[i]->get_startlambda();
      }
    }

//------------------------------------------------------------------------------
//----------------------------- REML estimation --------------------------------
//------------------------------------------------------------------------------


  bool remlest::estimate(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight)
  {
  // indices
  unsigned i,k,l;

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

  // some doubles
  double help;

  // Matrices to store old versions of beta and theta
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

  // Linear predictor
  statmatrix<double>eta(resp.rows(),1,0);

  // Working observations and weights
  statmatrix<double>worky(resp.rows(),1,0);
  statmatrix<double>workweight(resp.rows(),1,0);
  statmatrix<double>dinv(resp.rows(),1,0);
  statmatrix<double>mu(resp.rows(),1,0);

  // Matrix containing the inverse covariance matrix of the random effects
  statmatrix<double>Qinv(Z.cols(),1,0);

  // Matrices for Fisher scoring (regression parameters)
  statmatrix<double>H(beta.rows(),beta.rows(),0);
  statmatrix<double>H1(beta.rows(),1,0);

  // Matrices for Fisher scoring (variance parameters)
  statmatrix<double>Hinv(beta.rows(),beta.rows(),0);
  statmatrix<double>wresid(1,resp.rows(),0);                     //row vector !!

  // Transform smoothing paramater starting values to variances
  for(i=0; i<theta.rows(); i++)
    {
    theta(i,0)=1/theta(i,0);
    }

  // Estimation loop
  while(test==true)
    {
    // store current values in betaold and thetaold and compute Qinv
    betaold=beta;
    thetaold=theta;
    for(i=0, l=0; i<theta.rows(); i++)
      {
      for(k=zcut[i]; k<zcut[i+1]; k++, l++)
        {
        Qinv(l,0)=1/theta(i,0);
        }
      }

    // update linear predictor, working observations, weights, etc.
    eta=offset+X*beta.getRowBlock(0,X.cols())+Z*beta.getRowBlock(X.cols(),beta.rows());

    if(respfamily=="poisson")
      {
      for(i=0; i<eta.rows(); i++)
        {
        mu(i,0) = exp(eta(i,0));
        workweight(i,0) = mu(i,0)*weight(i,0);
        dinv(i,0) = 1/mu(i,0);
        }
      }
    else if(respfamily=="binomial")
      {
      for(i=0; i<eta.rows(); i++)
        {
        mu(i,0) = exp(eta(i,0))/(1+exp(eta(i,0)));
        dinv(i,0) = 1/(mu(i,0)*(1-mu(i,0)));
        workweight(i,0) = weight(i,0)/dinv(i,0);
        }
      }
    else if(respfamily=="binomialprobit")
      {
      for(i=0; i<eta.rows(); i++)
        {
        mu(i,0) = randnumbers::Phi2(eta(i,0));
        workweight(i,0) = weight(i,0)*(randnumbers::phi(eta(i,0)))*(randnumbers::phi(eta(i,0)))/(mu(i,0)*(1-mu(i,0)));
        dinv(i,0) = 1/randnumbers::phi(eta(i,0));
        }
      }
    dinv.elemmult(resp-mu);
    worky = eta - offset + dinv;

    // compute H and H1
    H.weightedsscp2(X,Z,workweight);
    H.addtodiag(Qinv,X.cols(),beta.rows());

    stop = check_pause();
    if (stop)
      return true;

    H1.weightedsscp_resp2(X,Z,worky,workweight);

    // Fisher-Scoring for beta
    beta=H.solve(H1);

    // update linear predictor and compute residuals
    eta=X*beta.getRowBlock(0,X.cols())+Z*beta.getRowBlock(X.cols(),beta.rows());
    for(i=0; i<eta.rows(); i++)
      {
      wresid(0,i)=(worky(i,0)-eta(i,0))*workweight(i,0);
      eta(i,0) += offset(i,0);
      }

    // transform theta
    for(i=0; i<theta.rows(); i++)
      {
      thetaold(i,0)=signs[i]*sqrt(thetaold(i,0));
      }

    Hinv=H.inverse();
    H.subfromdiag(Qinv,X.cols(),beta.rows());

    stop = check_pause();
    if (stop)
      return true;

    // compute score-function and expected fisher information

    for(i=0; i<theta.rows(); i++)
      {
      score(i,0)=-0.5*(H.getBlock(X.cols()+zcut[i],X.cols()+zcut[i],X.cols()+zcut[i+1],X.cols()+zcut[i+1])*thetaold(i,0)*2).trace()+
                 0.5*((H.getRowBlock(X.cols()+zcut[i],X.cols()+zcut[i+1]))*Hinv*(H.getColBlock(X.cols()+zcut[i],X.cols()+zcut[i+1]))*thetaold(i,0)*2).trace();
      for(l=zcut[i]; l<zcut[i+1]; l++)
        {
        help = (wresid*(Z.getCol(l)))(0,0);
        score(i,0) += 0.5*help*help*thetaold(i,0)*2;
        }

      for(k=0; k<theta.rows(); k++)
        {
        Fisher(i,k)=0.5*(H.getBlock(X.cols()+zcut[i],X.cols()+zcut[k],X.cols()+zcut[i+1],X.cols()+zcut[k+1])*H.getBlock(X.cols()+zcut[k],X.cols()+zcut[i],X.cols()+zcut[k+1],X.cols()+zcut[i+1])*thetaold(i,0)*4*thetaold(k,0)).trace()-
                    (H.getRowBlock(X.cols()+zcut[k],X.cols()+zcut[k+1])*Hinv*H.getColBlock(X.cols()+zcut[i],X.cols()+zcut[i+1])*H.getBlock(X.cols()+zcut[i],X.cols()+zcut[k],X.cols()+zcut[i+1],X.cols()+zcut[k+1])*thetaold(i,0)*4*thetaold(k,0)).trace()+
                    0.5*(H.getRowBlock(X.cols()+zcut[i],X.cols()+zcut[i+1])*Hinv*H.getColBlock(X.cols()+zcut[k],X.cols()+zcut[k+1])*H.getRowBlock(X.cols()+zcut[k],X.cols()+zcut[k+1])*Hinv*H.getColBlock(X.cols()+zcut[i],X.cols()+zcut[i+1])*thetaold(i,0)*4*thetaold(k,0)).trace();
        Fisher(k,i)=Fisher(i,k);
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
   help=eta.norm(0);
   for(i=0; i<theta.rows(); i++)
     {
     dinv=Z.getColBlock(zcut[i],zcut[i+1])*beta.getRowBlock(X.cols()+zcut[i],X.cols()+zcut[i+1]);
     stopcrit[i]=dinv.norm(0)/help;
     if(stopcrit[i]<lowerlim)
       {
       theta(i,0)=thetaold(i,0);
       }
     else
       {
       its[i]=it;
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

  // update workweight and H
  eta=offset+X*beta.getRowBlock(0,X.cols())+Z*beta.getRowBlock(X.cols(),beta.rows());
  if(respfamily=="poisson")
    {
    for(i=0; i<eta.rows(); i++)
      {
      mu(i,0) = exp(eta(i,0));
      workweight(i,0) = mu(i,0)*weight(i,0);
      }
    }
  else if(respfamily=="binomial")
    {
    for(i=0; i<eta.rows(); i++)
      {
      mu(i,0) = exp(eta(i,0))/(1+exp(eta(i,0)));
      dinv(i,0) = 1/(mu(i,0)*(1-mu(i,0)));
      workweight(i,0) = weight(i,0)/dinv(i,0);
      }
    }
  else if(respfamily=="binomialprobit")
    {
    for(i=0; i<eta.rows(); i++)
      {
      mu(i,0) = randnumbers::Phi2(eta(i,0));
      workweight(i,0) = weight(i,0)*(randnumbers::phi(eta(i,0)))*(randnumbers::phi(eta(i,0)))/(mu(i,0)*(1-mu(i,0)));
      }
    }
  H.weightedsscp2(X,Z,workweight);
  for(i=0, l=0; i<theta.rows(); i++)
    {
    for(k=zcut[i]; k<zcut[i+1]; k++, l++)
      {
      Qinv(l,0)=1/theta(i,0);
      }
    }
  H.addtodiag(Qinv,X.cols(),beta.rows());
  Hinv=H.inverse();

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

  for(i=1;i<fullcond.size();i++)
    {
    beta(0,0) += fullcond[i]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[i],zcut[i-1],i-1,false,xcut[i],X.cols()+zcut[i-1],0,false,i);
    }
  beta(0,0) += fullcond[0]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[0],0,0,false,xcut[0],0,0,false,0);

  out("\n");
  out("  Additive predictor and expectations\n",true);
  out("\n");
  out("\n");
  out("  Additive predictor and expectation for each observation are stored in file\n");
  out("  "+outfile+"_predict.raw\n");
  out("\n");

  ofstream outpredict((outfile+"_predict.raw").strtochar());
  outpredict << "eta mu" << endl;
  for(i=0; i<eta.rows(); i++)
    {
    outpredict << eta(i,0) << " " << mu(i,0) << endl;
    }
  outpredict.close();

  return false;
  }

bool remlest::estimate_glm(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight)
  {
  // indices
  unsigned i;//,k,l;

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

  // some doubles
  double help;

  // Matrix to store old version of beta
  statmatrix<double>betaold(beta.rows(),1,0);

  // Number of iterations
  unsigned it=1;

  // Criteria to deterine convergence
  double crit1=1;                //relative changes in regression parameters
  bool test=true;

  // Linear predictor
  statmatrix<double>eta(resp.rows(),1,0);

  // Working observations and weights
  statmatrix<double>worky(resp.rows(),1,0);
  statmatrix<double>workweight(resp.rows(),1,0);
  statmatrix<double>dinv(resp.rows(),1,0);
  statmatrix<double>mu(resp.rows(),1,0);

  // Matrices for Fisher scoring (regression parameters)
  statmatrix<double>H(beta.rows(),beta.rows(),0);
  statmatrix<double>H1(beta.rows(),1,0);

  // Estimation loop
  while(test==true)
    {
    // store current values in betaold
    betaold=beta;

    // update linear predictor, working observations, weights, etc.
    eta=offset+X*beta;

    if(respfamily=="poisson")
      {
      for(i=0; i<eta.rows(); i++)
        {
        mu(i,0) = exp(eta(i,0));
        workweight(i,0) = mu(i,0)*weight(i,0);
        dinv(i,0) = 1/mu(i,0);
        }
      }
    else if(respfamily=="binomial")
      {
      for(i=0; i<eta.rows(); i++)
        {
        mu(i,0) = exp(eta(i,0))/(1+exp(eta(i,0)));
        dinv(i,0) = 1/(mu(i,0)*(1-mu(i,0)));
        workweight(i,0) = weight(i,0)/dinv(i,0);
        }
      }
    else if(respfamily=="binomialprobit")
      {
      for(i=0; i<eta.rows(); i++)
        {
        mu(i,0) = randnumbers::Phi2(eta(i,0));
        workweight(i,0) = weight(i,0)*(randnumbers::phi(eta(i,0)))*(randnumbers::phi(eta(i,0)))/(mu(i,0)*(1-mu(i,0)));
        dinv(i,0) = 1/randnumbers::phi(eta(i,0));
        }
      }
    dinv.elemmult(resp-mu);
    worky = eta - offset + dinv;

    // compute H and H1
    H.weightedsscp(X,workweight);
    H1.weightedsscp_resp(X,worky,workweight);

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

  // update workweight and H
  eta=offset+X*beta;
  if(respfamily=="poisson")
    {
    for(i=0; i<eta.rows(); i++)
      {
      mu(i,0) = exp(eta(i,0));
      workweight(i,0) = mu(i,0)*weight(i,0);
      }
    }
  else if(respfamily=="binomial")
    {
    for(i=0; i<eta.rows(); i++)
      {
      mu(i,0) = exp(eta(i,0))/(1+exp(eta(i,0)));
      dinv(i,0) = 1/(mu(i,0)*(1-mu(i,0)));
      workweight(i,0) = weight(i,0)/dinv(i,0);
      }
    }
  else if(respfamily=="binomialprobit")
    {
    for(i=0; i<eta.rows(); i++)
      {
      mu(i,0) = randnumbers::Phi2(eta(i,0));
      workweight(i,0) = weight(i,0)*(randnumbers::phi(eta(i,0)))*(randnumbers::phi(eta(i,0)))/(mu(i,0)*(1-mu(i,0)));
      }
    }
  H.weightedsscp(X,workweight);
  H=H.inverse();

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


  for(i=1;i<fullcond.size();i++)
    {
    beta(0,0) += fullcond[i]->outresultsreml(X,Z,beta,H,datamatrix(1,1,0),xcut[i],zcut[i-1],i-1,false,xcut[i],X.cols()+zcut[i-1],0,false,i);
    }
  beta(0,0) += fullcond[0]->outresultsreml(X,Z,beta,H,datamatrix(1,1,0),xcut[0],0,0,false,xcut[0],0,0,false,0);

  out("\n");
  out("  Linear predictor and expectations\n",true);
  out("\n");
  out("\n");
  out("  Linear predictor and expectation for each observation are stored in file\n");
  out("  "+outfile+"_predict.raw\n");
  out("\n");

  ofstream outpredict((outfile+"_predict.raw").strtochar());
  outpredict << "eta mu" << endl;
  for(i=0; i<eta.rows(); i++)
    {
    outpredict << eta(i,0) << " " << mu(i,0) << endl;
    }
  outpredict.close();


  return false;
  }


//------------------------------------------------------------------------------
//-------------- REML estimation with dispersion parameter ---------------------
//------------------------------------------------------------------------------

bool remlest::estimate_dispers(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight)
  {
  // indices
  unsigned i,k,l;

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

  // some doubles
  double help;

  // Matrices to store old versions of beta and theta
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
  vector<unsigned>its(theta.rows(),0);
  vector<int>signs(theta.rows(),1);

  // Linear predictor
  statmatrix<double>eta(resp.rows(),1,0);

  // Working observations and weights
  statmatrix<double>worky(resp.rows(),1,0);
  statmatrix<double>workweight(resp.rows(),1,0);
  statmatrix<double>dinv(resp.rows(),1,0);
  statmatrix<double>mu(resp.rows(),1,0);

  // Matrix containing the inverse covariance matrix of the random effects
  statmatrix<double>Qinv(Z.cols(),1,0);

  // Matrices for Fisher scoring (regression parameters)
  statmatrix<double>H(beta.rows(),beta.rows(),0);
  statmatrix<double>H1(beta.rows(),1,0);

  // Matrices for Fisher scoring (variance parameters)
  statmatrix<double>Hinv(beta.rows(),beta.rows(),0);
  statmatrix<double>wresid(1,resp.rows(),0);                     //row vector !!
  statmatrix<double>w1resid(resp.rows(),1,0);

  // Transform smoothing paramater starting values to variances
  theta(theta.rows()-1,0)=resp.var(0);
  for(i=0; i<theta.rows()-1; i++)
    {
    theta(i,0)=theta(theta.rows()-1,0)/theta(i,0);
    }

  // Substract offset from response for gaussian response.

  if(respfamily=="gaussian")
    {
    for(i=0; i<resp.rows(); i++)
      {
      worky(i,0) = resp(i,0)-offset(i,0);
      }
    }

  // Estimation loop
  while(test==true)
    {
    // store current values in betaold and thetaold and compute Qinv
    betaold=beta;
    thetaold=theta;

    for(i=0, l=0; i<theta.rows()-1; i++)
      {
      for(k=zcut[i]; k<zcut[i+1]; k++, l++)
        {
        Qinv(l,0)=1/theta(i,0);
        }
      }

    // update linear predictor, working observations, weights, etc.
    eta=offset+X*beta.getRowBlock(0,X.cols())+Z*beta.getRowBlock(X.cols(),beta.rows());

    if(respfamily=="gaussian")
      {
      workweight=1/theta(theta.rows()-1,0)*weight;
      }
    else if(respfamily=="poissondispers")
      {
      for(i=0; i<eta.rows(); i++)
        {
        mu(i,0) = exp(eta(i,0));
        workweight(i,0) = 1/theta(theta.rows()-1,0)*mu(i,0)*weight(i,0);
        dinv(i,0) = 1/mu(i,0);
        worky(i,0)=eta(i,0)-offset(i,0)+dinv(i,0)*(resp(i,0)-mu(i,0));
        }
      }
    else if(respfamily=="binomialdispers")
      {
      for(i=0; i<eta.rows(); i++)
        {
        mu(i,0) = exp(eta(i,0))/(1+exp(eta(i,0)));
        dinv(i,0) = 1/(mu(i,0)*(1-mu(i,0)));
        workweight(i,0) = 1/theta(theta.rows()-1,0)*weight(i,0)/dinv(i,0);
        worky(i,0)=eta(i,0)-offset(i,0)+dinv(i,0)*(resp(i,0)-mu(i,0));
        }
      }
    else if(respfamily=="binomialprobitdispers")
      {
      for(i=0; i<eta.rows(); i++)
        {
        mu(i,0) = randnumbers::Phi2(eta(i,0));
        workweight(i,0) = 1/theta(theta.rows()-1,0)*weight(i,0)*(randnumbers::phi(eta(i,0)))*(randnumbers::phi(eta(i,0)))/(mu(i,0)*(1-mu(i,0)));
        dinv(i,0) = 1/randnumbers::phi(eta(i,0));
        worky(i,0)=eta(i,0)-offset(i,0)+dinv(i,0)*(resp(i,0)-mu(i,0));
        }
      }
    else if(respfamily=="gamma")
      {
      for(i=0; i<eta.rows(); i++)
        {
        mu(i,0) = exp(eta(i,0));
        workweight(i,0) = 1/theta(theta.rows()-1,0)*weight(i,0);
        dinv(i,0) = 1/mu(i,0);
//        worky(i,0)=eta(i,0)-offset(i,0)+dinv(i,0)*(resp(i,0)-mu(i,0));
        worky(i,0)=eta(i,0)-offset(i,0)+dinv(i,0)*resp(i,0)-1;
        }
      }

    // compute H and H1
    H.weightedsscp2(X,Z,workweight);
    H.addtodiag(Qinv,X.cols(),beta.rows());

    stop = check_pause();
    if (stop)
      return true;

    H1.weightedsscp_resp2(X,Z,worky,workweight);

    // Fisher-Scoring for beta
    beta=H.solve(H1);

    // update linear predictor and compute weighted residuals
    eta=X*beta.getRowBlock(0,X.cols())+Z*beta.getRowBlock(X.cols(),beta.rows());
    for(i=0; i<eta.rows(); i++)
      {
      wresid(0,i)=(worky(i,0)-eta(i,0))*workweight(i,0);
      w1resid(i,0)=(worky(i,0)-eta(i,0))*sqrt(workweight(i,0));
      eta(i,0) += offset(i,0);
      }

    // transform theta
    for(i=0; i<theta.rows(); i++)
      {
      thetaold(i,0)=signs[i]*sqrt(thetaold(i,0));
      }

    Hinv=H.inverse();
    H.subfromdiag(Qinv,X.cols(),beta.rows());

    stop = check_pause();
    if (stop)
      return true;

    // compute score-function and expected fisher information

    for(i=0; i<theta.rows()-1; i++)
      {
      score(i,0)=-0.5*(H.getBlock(X.cols()+zcut[i],X.cols()+zcut[i],X.cols()+zcut[i+1],X.cols()+zcut[i+1])*thetaold(i,0)*2).trace()+
                 0.5*((H.getRowBlock(X.cols()+zcut[i],X.cols()+zcut[i+1]))*Hinv*(H.getColBlock(X.cols()+zcut[i],X.cols()+zcut[i+1]))*thetaold(i,0)*2).trace();
      for(l=zcut[i]; l<zcut[i+1]; l++)
        {
        help = (wresid*(Z.getCol(l)))(0,0);
        score(i,0) += 0.5*help*help*thetaold(i,0)*2;
        }
      Fisher(theta.rows()-1,i)=0.5*(H.getBlock(X.cols()+zcut[i],X.cols()+zcut[i],X.cols()+zcut[i+1],X.cols()+zcut[i+1])*4*thetaold(i,0)/thetaold(theta.rows()-1,0)).trace()-
                  (H.getColBlock(X.cols()+zcut[i],X.cols()+zcut[i+1])*H.getRowBlock(X.cols()+zcut[i],X.cols()+zcut[i+1])*Hinv*4*thetaold(i,0)/thetaold(theta.rows()-1,0)).trace()+
                  0.5*(H*Hinv*H.getColBlock(X.cols()+zcut[i],X.cols()+zcut[i+1])*H.getRowBlock(X.cols()+zcut[i],X.cols()+zcut[i+1])*Hinv*4*thetaold(i,0)/thetaold(theta.rows()-1,0)).trace();
      Fisher(i,theta.rows()-1)=Fisher(theta.rows()-1,i);
      for(k=0; k<theta.rows()-1; k++)
        {
        Fisher(i,k)=0.5*(H.getBlock(X.cols()+zcut[i],X.cols()+zcut[k],X.cols()+zcut[i+1],X.cols()+zcut[k+1])*H.getBlock(X.cols()+zcut[k],X.cols()+zcut[i],X.cols()+zcut[k+1],X.cols()+zcut[i+1])*thetaold(i,0)*4*thetaold(k,0)).trace()-
                    (H.getRowBlock(X.cols()+zcut[k],X.cols()+zcut[k+1])*Hinv*H.getColBlock(X.cols()+zcut[i],X.cols()+zcut[i+1])*H.getBlock(X.cols()+zcut[i],X.cols()+zcut[k],X.cols()+zcut[i+1],X.cols()+zcut[k+1])*thetaold(i,0)*4*thetaold(k,0)).trace()+
                    0.5*(H.getRowBlock(X.cols()+zcut[i],X.cols()+zcut[i+1])*Hinv*H.getColBlock(X.cols()+zcut[k],X.cols()+zcut[k+1])*H.getRowBlock(X.cols()+zcut[k],X.cols()+zcut[k+1])*Hinv*H.getColBlock(X.cols()+zcut[i],X.cols()+zcut[i+1])*thetaold(i,0)*4*thetaold(k,0)).trace();
        Fisher(k,i)=Fisher(i,k);
        }
      }
    score(theta.rows()-1,0)=-(double)resp.rows()/thetaold(theta.rows()-1,0)+
          0.5*(2/thetaold(theta.rows()-1,0)*H*Hinv).trace()+
          w1resid.sum2(0)/thetaold(theta.rows()-1,0);

    Fisher(theta.rows()-1,theta.rows()-1)=2*(double)resp.rows()/(thetaold(theta.rows()-1,0)*thetaold(theta.rows()-1,0))-
           (H*Hinv*4/(thetaold(theta.rows()-1,0)*thetaold(theta.rows()-1,0))).trace()+
           0.5*(H*Hinv*H*Hinv*4/(thetaold(theta.rows()-1,0)*thetaold(theta.rows()-1,0))).trace();

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
   help=eta.norm(0);
   for(i=0; i<theta.rows()-1; i++)
     {
     dinv=Z.getColBlock(zcut[i],zcut[i+1])*beta.getRowBlock(X.cols()+zcut[i],X.cols()+zcut[i+1]);
     stopcrit[i]=dinv.norm(0)/help;
     if(stopcrit[i]<lowerlim)
       {
       theta(i,0)=thetaold(i,0);
       }
     else
       {
       its[i]=it;
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

  // update workweight and H
  eta=offset+X*beta.getRowBlock(0,X.cols())+Z*beta.getRowBlock(X.cols(),beta.rows());
  if(respfamily=="gaussian")
    {
    workweight=1/theta(theta.rows()-1,0)*weight;
    }
  else if(respfamily=="poissondispers")
    {
    for(i=0; i<eta.rows(); i++)
      {
      mu(i,0) = exp(eta(i,0));
      workweight(i,0) = 1/theta(theta.rows()-1,0)*mu(i,0)*weight(i,0);
      }
    }
  else if(respfamily=="binomialdispers")
    {
    for(i=0; i<eta.rows(); i++)
      {
      mu(i,0) = exp(eta(i,0))/(1+exp(eta(i,0)));
      dinv(i,0) = 1/(mu(i,0)*(1-mu(i,0)));
      workweight(i,0) = 1/theta(theta.rows()-1,0)*weight(i,0)/dinv(i,0);
      }
    }
  else if(respfamily=="binomialprobitdispers")
    {
    for(i=0; i<eta.rows(); i++)
      {
      mu(i,0) = randnumbers::Phi2(eta(i,0));
      workweight(i,0) = 1/theta(theta.rows()-1,0)*weight(i,0)*(randnumbers::phi(eta(i,0)))*(randnumbers::phi(eta(i,0)))/(mu(i,0)*(1-mu(i,0)));
      }
    }
  else if(respfamily=="gamma")
    {
    for(i=0; i<eta.rows(); i++)
      {
      mu(i,0) = exp(eta(i,0));
      workweight(i,0) = 1/theta(theta.rows()-1,0)*weight(i,0);
      }
    }
  H.weightedsscp2(X,Z,workweight);
  for(i=0, l=0; i<theta.rows()-1; i++)
    {
    for(k=zcut[i]; k<zcut[i+1]; k++, l++)
      {
      Qinv(l,0)=1/theta(i,0);
      }
    }
  H.addtodiag(Qinv,X.cols(),beta.rows());
  Hinv=H.inverse();

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

  datamatrix thetareml(theta.rows(),3,0);
  thetareml.putCol(0,theta);
  for(i=0; i<theta.rows()-1; i++)
    {
    if(stopcrit[i]<lowerlim)
      {
      thetareml(i,1)=1;
      }
    thetareml(i,2)=its[i];
    }

  out("\n");
  out("  Estimated scale parameter: " + ST::doubletostring(theta(theta.rows()-1,0),6) + "\n");
  out("\n");
  out("  Scale parameter is also stored in file\n");
  out("  " + outfile + "_scale.res");
  out("\n");

  ofstream outscale((outfile + "_scale.res").strtochar());
  outscale << "scale" << endl;
  outscale << theta(theta.rows()-1,0) << endl;
  outscale.close();

  for(i=1;i<fullcond.size();i++)
    {
    beta(0,0) += fullcond[i]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[i],zcut[i-1],i-1,true,xcut[i],X.cols()+zcut[i-1],0,false,i);
    }
  beta(0,0) += fullcond[0]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[0],0,0,true,xcut[0],0,0,false,0);

  out("\n");
  out("  Additive predictor and expectations\n",true);
  out("\n");
  out("\n");
  out("  Additive predictor and expectation for each observation are stored in file\n");
  out("  "+outfile+"_predict.raw\n");
  out("\n");

  ofstream outpredict((outfile+"_predict.raw").strtochar());
  outpredict << "eta mu" << endl;
  if(respfamily=="gaussian")
    {
    for(i=0; i<eta.rows(); i++)
      {
      outpredict << eta(i,0) << " " << eta(i,0) << endl;
      }
    }
  else
    {
    for(i=0; i<eta.rows(); i++)
      {
      outpredict << eta(i,0) << " " << mu(i,0) << endl;
      }
    }
  outpredict.close();

  return false;
  }

bool remlest::estimate_glm_dispers(const datamatrix resp,
                                   const datamatrix & offset,
                                   const datamatrix & weight)
  {

  theta(0,0)=resp.var(0);

  // indices
  unsigned i;

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

  // some doubles
  double help;

  // Matrix to store old version of beta
  statmatrix<double>betaold(beta.rows(),1,0);
  statmatrix<double>thetaold(1,1,0);

  // Score-function and expected Fisher information for theta
  statmatrix<double>score(1,1,0);
  statmatrix<double>Fisher(1,1,0);

  // Number of iterations
  unsigned it=1;

  // Criteria to deterine convergence
  double crit1=1;                //relative changes in regression parameters
  double crit2=1;                //relative changes in variance parameters
  bool test=true;

  // Linear predictor
  statmatrix<double>eta(resp.rows(),1,0);

  // Working observations and weights
  statmatrix<double>worky(resp.rows(),1,0);
  statmatrix<double>workweight(resp.rows(),1,0);
  statmatrix<double>dinv(resp.rows(),1,0);
  statmatrix<double>mu(resp.rows(),1,0);

  // Matrices for Fisher scoring (variance parameters)
  statmatrix<double>w1resid(resp.rows(),1,0);

  // Matrices for Fisher scoring (regression parameters)
  statmatrix<double>H(beta.rows(),beta.rows(),0);
  statmatrix<double>H1(beta.rows(),1,0);

  // Substract offset from response for gaussian response.

  if(respfamily=="gaussian")
    {
    for(i=0; i<resp.rows(); i++)
      {
      worky(i,0) = resp(i,0)-offset(i,0);
      }
    }

  // Estimation loop
  while(test==true)
    {
    // store current values in betaold and thetaold
    betaold=beta;
    thetaold=theta;

    // update linear predictor, working observations, weights, etc.
    eta=offset+X*beta;

    if(respfamily=="gaussian")
      {
      workweight=1/theta(0,0)*weight;
      }
    else if(respfamily=="poissondispers")
      {
      for(i=0; i<eta.rows(); i++)
        {
        mu(i,0) = exp(eta(i,0));
        workweight(i,0) = 1/theta(theta.rows()-1,0)*mu(i,0)*weight(i,0);
        dinv(i,0) = 1/mu(i,0);
        worky(i,0)=eta(i,0)-offset(i,0)+dinv(i,0)*(resp(i,0)-mu(i,0));
        }
      }
    else if(respfamily=="binomialdispers")
      {
      for(i=0; i<eta.rows(); i++)
        {
        mu(i,0) = exp(eta(i,0))/(1+exp(eta(i,0)));
        dinv(i,0) = 1/(mu(i,0)*(1-mu(i,0)));
        workweight(i,0) = 1/theta(theta.rows()-1,0)*weight(i,0)/dinv(i,0);
        worky(i,0)=eta(i,0)-offset(i,0)+dinv(i,0)*(resp(i,0)-mu(i,0));
        }
      }
    else if(respfamily=="binomialprobitdispers")
      {
      for(i=0; i<eta.rows(); i++)
        {
        mu(i,0) = randnumbers::Phi2(eta(i,0));
        workweight(i,0) = 1/theta(theta.rows()-1,0)*weight(i,0)*(randnumbers::phi(eta(i,0)))*(randnumbers::phi(eta(i,0)))/(mu(i,0)*(1-mu(i,0)));
        dinv(i,0) = 1/randnumbers::phi(eta(i,0));
        worky(i,0)=eta(i,0)-offset(i,0)+dinv(i,0)*(resp(i,0)-mu(i,0));
        }
      }
    else if(respfamily=="gamma")
      {
      for(i=0; i<eta.rows(); i++)
        {
        mu(i,0) = exp(eta(i,0));
        workweight(i,0) = 1/theta(theta.rows()-1,0)*weight(i,0);
        dinv(i,0) = 1/mu(i,0);
        worky(i,0)=eta(i,0)-offset(i,0)+dinv(i,0)*(resp(i,0)-mu(i,0));
        }
      }

    // compute H and H1
    H.weightedsscp(X,workweight);
    H1.weightedsscp_resp(X,worky,workweight);

    stop = check_pause();
    if (stop)
      return true;

    // Fisher-Scoring for beta
    beta=H.solve(H1);

    // update linear predictor and compute weighted residuals
    eta=X*beta;
    for(i=0; i<eta.rows(); i++)
      {
      w1resid(i,0)=(worky(i,0)-eta(i,0))*sqrt(workweight(i,0));
      eta(i,0) += offset(i,0);
      }

    // compute score-function and expected fisher information

    score(0,0)=-0.5*(resp.rows()-beta.rows())+0.5*w1resid.sum2(0);
    Fisher(0,0)=(resp.rows()-beta.rows())/(2*thetaold(0,0));

    // fisher scoring for theta
    theta(0,0) = thetaold(0,0) + score(0,0)/Fisher(0,0);

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
    test=((crit1>eps) || (crit2>eps) && (it<(unsigned)maxit));

    // count iteration
    it=it+1;

    }

  // update workweight and H
  eta=X*beta;
  if(respfamily=="gaussian")
    {
    workweight=1/theta(0,0)*weight;
    }
  else if(respfamily=="poissondispers")
    {
    for(i=0; i<eta.rows(); i++)
      {
      mu(i,0) = exp(eta(i,0));
      workweight(i,0) = 1/theta(theta.rows()-1,0)*mu(i,0)*weight(i,0);
      }
    }
  else if(respfamily=="binomialdispers")
    {
    for(i=0; i<eta.rows(); i++)
      {
      mu(i,0) = exp(eta(i,0))/(1+exp(eta(i,0)));
      dinv(i,0) = 1/(mu(i,0)*(1-mu(i,0)));
      workweight(i,0) = 1/theta(theta.rows()-1,0)*weight(i,0)/dinv(i,0);
      }
    }
  else if(respfamily=="binomialprobitdispers")
    {
    for(i=0; i<eta.rows(); i++)
      {
      mu(i,0) = randnumbers::Phi2(eta(i,0));
      workweight(i,0) = 1/theta(theta.rows()-1,0)*weight(i,0)*(randnumbers::phi(eta(i,0)))*(randnumbers::phi(eta(i,0)))/(mu(i,0)*(1-mu(i,0)));
      }
    }
  else if(respfamily=="gamma")
    {
    for(i=0; i<eta.rows(); i++)
      {
      mu(i,0) = exp(eta(i,0));
      workweight(i,0) = 1/theta(theta.rows()-1,0)*weight(i,0);
      }
    }
  H.weightedsscp(X,workweight);
  H=H.inverse();

  //write results

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

  out("\n");
  out("  Estimated scale parameter: " + ST::doubletostring(theta(theta.rows()-1,0),6) + "\n");
  out("\n");
  out("  Scale parameter is also stored in file\n");
  out("  " + outfile + "_scale.res");
  out("\n");

  ofstream outscale((outfile + "_scale.res").strtochar());
  outscale << "scale" << endl;
  outscale << theta(theta.rows()-1,0) << endl;
  outscale.close();

  for(i=1;i<fullcond.size();i++)
    {
    beta(0,0) += fullcond[i]->outresultsreml(X,Z,beta,H,theta,xcut[i],zcut[i-1],i-1,true,xcut[i],X.cols()+zcut[i-1],0,false,i);
    }
  beta(0,0) += fullcond[0]->outresultsreml(X,Z,beta,H,theta,xcut[0],0,0,true,xcut[0],0,0,false,0);

  out("\n");
  out("  Linear predictor and expectations\n",true);
  out("\n");
  out("\n");
  out("  Linear predictor and expectation for each observation are stored in file\n");
  out("  "+outfile+"_predict.raw\n");
  out("\n");

  ofstream outpredict((outfile+"_predict.raw").strtochar());
  outpredict << "eta mu" << endl;
  if(respfamily=="gaussian")
    {
    for(i=0; i<eta.rows(); i++)
      {
      outpredict << eta(i,0) << " " << eta(i,0) << endl;
      }
    }
  else
    {
    for(i=0; i<eta.rows(); i++)
      {
      outpredict << eta(i,0) << " " << mu(i,0) << endl;
      }
    }
  outpredict.close();

  return false;
  }

//------------------------------------------------------------------------------
//------------------------------- Survival data --------------------------------
//------------------------------------------------------------------------------

bool remlest::estimate_survival(const datamatrix resp,
                const datamatrix & offset, const datamatrix & weight)
  {

  unsigned i, j, k, l;
  double help;

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

  // Matrix to store old versions of beta and theta
  statmatrix<double>betaold(beta.rows(),1,0);
  statmatrix<double>thetaold(theta.rows(),1,0);

  // Score-function and expected Fisher information for theta
  statmatrix<double>score(theta.rows(),1,0);
  statmatrix<double>Fisher(theta.rows(),theta.rows(),0);

  // Matrices for Fisher scoring (regression parameters)
  statmatrix<double>H(beta.rows(),beta.rows(),0);
  statmatrix<double>Hinv(beta.rows(),beta.rows(),0);
  statmatrix<double>H1(beta.rows(),1,0);

  // Number of iterations, nr of observations
  unsigned it=1;
  unsigned nrobs=Z.rows();
  unsigned xcols = X.cols();
  unsigned zcols = Z.cols();

  // Criteria to detemine convergence
  double crit1=1;                //relative changes in regression parameters
  double crit2=1;                //relative changes in variance parameters
  bool test=true;

  vector<double>stopcrit(theta.rows(),10);
  vector<int>its(theta.rows(),0);
  vector<int>signs(theta.rows(),1);

  // Matrix containing the inverse covariance matrix of the random effects
  statmatrix<double>Qinv(zcols,1,0);

  // Inzidenzmatrix, die für jeden Wert in fullcond bzw. beta angibt, ob er zur Baseline-HR beiträgt
  vector<int>isbaseline(fullcond.size(),0);
  int nrbaseline=0;
  for(i=0; i<fullcond.size(); i++)
    {
    if(fullcond[i]->is_baseline()==true)
      {
      isbaseline[i]=1;
      nrbaseline++;
      }
    }

  vector<int>isbaselinebeta(beta.rows(),0);
  vector<int>fc_pos(beta.rows(),0);
  vector<int>dm_pos(beta.rows(),0);
  l=0;
  for(i=0; i<fullcond.size(); i++)
    {
    if(isbaseline[i]==1)
      {
      k=0;
      for(j=zcut[i-1]; j<zcut[i]; j++, k++)
        {
        isbaselinebeta[xcols+j]=1;
        fc_pos[xcols+j]=l;
        dm_pos[xcols+j]=k;
        }
      k=0;
      for(j=xcut[i]; j<xcut[i+1]; j++, k++)
        {
        isbaselinebeta[j]=1;
        fc_pos[j]=l;
        dm_pos[j]=k;
        if(xcut[i+1]==xcut[i]+1)
          {
          dm_pos[j]=1;
          }
        }
      l++;
      }
    }

  bool check;
  if(nrbaseline>1)
    {
    check=false;
    }
  else
    {
    check=true;
    }

  // Matrices and variables for baseline effects
//  double tstep;
  statmatrix<double> tsteps;
  datamatrix t_X;
  datamatrix t_Z;
  vector<unsigned> tstart;
  vector<unsigned> tend;
  datamatrix interactvar(nrobs,nrbaseline,0);
  j=0;
  for(i=0; i<fullcond.size(); i++)
    {
    if(isbaseline[i]==1)
      {
      fullcond[i]->initialize_baseline(j,t_X,t_Z,tstart,tend,interactvar,tsteps);
      j++;
      }
    }
  datamatrix basef(t_X.rows(),nrbaseline,0);
  statmatrix<double>baseline;
  // Baseline-HR, linear predictor
  if(check)
    {
    baseline = statmatrix<double>(t_X.rows(),1,0);
    }
  else
    {
    baseline = statmatrix<double>(nrobs,t_X.rows(),0);
    }
  statmatrix<double>cumbaseline(nrobs,1,0);
  statmatrix<double>cumhazard(nrobs,1,0);
  statmatrix<double>eta(nrobs,1,0);
  statmatrix<double>baseline_eta(nrobs,1,0);
  statmatrix<double>mult_eta(nrobs,1,0);
  statmatrix<double>mult_hazard(nrobs,1,0);
  statmatrix<double>helpmat(nrobs,1,0);

  // Transform smoothing paramater starting values to variances
  for(i=0; i<theta.rows(); i++)
    {
    theta(i,0)=1/theta(i,0);
    }

  // Startwert für beta0 ist der ML-Schätzer bei konstanter Rate + Poisson-Verteilung
  beta(0,0) = log(resp.sum(0)/t_X(t_X.rows()-1,1));

  while(test==true)
    {

    // store current values in betaold and thetaold
    betaold=beta;
    thetaold=theta;

    // compute Qinv
    for(i=0, j=0; i<theta.rows(); i++)
      {
      for(k=zcut[i]; k<zcut[i+1]; k++, j++)
        {
        Qinv(j,0)=1/theta(i,0);
        }
      }

    // compute basef

    j=0;
    for(i=0; i<fullcond.size(); i++)
      {
      if(isbaseline[i]==1)
        {
        if(xcut[i+1]==xcut[i]+1)
          {
          basef.putCol(j,beta(xcut[i],0)*t_X.getCol(1)+t_Z*beta.getRowBlock(X.cols()+zcut[i-1],X.cols()+zcut[i]));
          }
        else
          {
          basef.putCol(j,t_X*beta.getRowBlock(xcut[i],xcut[i+1])+t_Z*beta.getRowBlock(X.cols()+zcut[i-1],X.cols()+zcut[i]));
          }
        j++;
        }
      }

    // CHECK
    if(check)
      {
      for(i=0; i< t_X.rows();i++)
        {
        baseline(i,0)=exp(basef(i,0));
        }
      }
    else
      {
      baseline = interactvar*basef.transposed();
      for(i=0; i<baseline.rows(); i++)
        {
        for(j=tstart[i]; j<tend[i];j++)
          {
          baseline(i,j) = exp(baseline(i,j));
          }
        }
      }

    // compute eta, baseline-hazard, hazard, cum. baseline, cum. hazard

    baseline_eta=datamatrix(nrobs,1,0);
    for(i=0; i<fullcond.size(); i++)
      {
      if(isbaseline[i]==1)
        {
        baseline_eta=baseline_eta+Z.getColBlock(zcut[i-1],zcut[i])*beta.getRowBlock(X.cols()+zcut[i-1],X.cols()+zcut[i]);
        if(xcut[i]<xcut[i+1])
          {
          baseline_eta=baseline_eta+X.getColBlock(xcut[i],xcut[i+1])*beta.getRowBlock(xcut[i],xcut[i+1]);
          }
        }
      }
    eta=X*beta.getRowBlock(0,X.cols())+Z*beta.getRowBlock(X.cols(),beta.rows());
    mult_eta=eta-baseline_eta;
    cumbaseline=datamatrix(nrobs,1,0);

    //CHECK
    if(check)
      {
      for(i=0; i<nrobs; i++)
        {
        for(k=tstart[i]; k<tend[i]-1; k++)
          {
          cumbaseline(i,0) += 0.5*tsteps(k,0)*(baseline(k,0)+baseline(k+1,0));
          }
        }
      }
    else
      {
      for(i=0; i<nrobs; i++)
        {
        for(k=tstart[i]; k<tend[i]-1; k++)
          {
//          cumbaseline(i,0) += 0.5*tstep*(exp((interactvar.getRow(i).transposed()*basef.getRow(k))(0,0))
//                          +exp((interactvar.getRow(i).transposed()*basef.getRow(k+1))(0,0)));
          cumbaseline(i,0) += 0.5*tsteps(k,0)*(baseline(i,k)+baseline(i,k+1));
          }
        }
      }

    for(i=0; i<nrobs; i++)
      {
      mult_hazard(i,0)=exp(mult_eta(i,0));
      cumhazard(i,0)=cumbaseline(i,0)*mult_hazard(i,0);
      }

    // Score-Funktion für beta

    for(j=0; j<X.cols(); j++)
      {
      if(isbaselinebeta[j]==1)
        {
        helpmat = datamatrix(nrobs,1,0);
        //CHECK
        if(check)
          {
          for(i=0; i<nrobs; i++)
            {
            for(k=tstart[i]; k<tend[i]-1; k++)
              {
              helpmat(i,0) += 0.5*tsteps(k,0)*(t_X(k,dm_pos[j])*baseline(k,0)+t_X(k+1,dm_pos[j])*baseline(k+1,0));
              }
            }
          }
        else
          {
          for(i=0; i<nrobs; i++)
            {
            for(k=tstart[i]; k<tend[i]-1; k++)
              {
//              helpmat(i,0) += 0.5*tstep*interactvar(i,fc_pos[j])*
//                               (t_X(k,dm_pos[j])*exp((interactvar.getRow(i).transposed()*basef.getRow(k))(0,0))
//                               +t_X(k+1,dm_pos[j])*exp((interactvar.getRow(i).transposed()*basef.getRow(k+1))(0,0)));
              helpmat(i,0) += 0.5*tsteps(k,0)*interactvar(i,fc_pos[j])*
                               (t_X(k,dm_pos[j])*baseline(i,k)+t_X(k+1,dm_pos[j])*baseline(i,k+1));
              }
            }
          }
        for(i=0; i<nrobs; i++)
          {
          helpmat(i,0) *= mult_hazard(i,0);
          }
        H1(j,0)=(resp.transposed()*X.getCol(j))(0,0)-helpmat.sum(0);
        }
      else
        {
        H1(j,0)=(resp.transposed()*X.getCol(j))(0,0)-(cumhazard.transposed()*X.getCol(j))(0,0);
        }
      }

    for(j=0; j<zcols; j++)
      {
      if(isbaselinebeta[xcols+j]==1)
        {
        helpmat = datamatrix(nrobs,1,0);
        //CHECK
        if(check)
          {
          for(i=0; i<nrobs; i++)
            {
            for(k=tstart[i]; k<tend[i]-1; k++)
              {
              helpmat(i,0) += 0.5*tsteps(k,0)*(t_Z(k,dm_pos[xcols+j])*baseline(k,0)
                               +t_Z(k+1,dm_pos[xcols+j])*baseline(k+1,0));
              }
            }
          }
        else
          {
          for(i=0; i<nrobs; i++)
            {
            for(k=tstart[i]; k<tend[i]-1; k++)
              {
//              helpmat(i,0) += 0.5*tstep*interactvar(i,fc_pos[xcols+j])*
//                              (t_Z(k,dm_pos[xcols+j])*exp((interactvar.getRow(i).transposed()*basef.getRow(k))(0,0))
//                               +t_Z(k+1,dm_pos[xcols+j])*exp((interactvar.getRow(i).transposed()*basef.getRow(k+1))(0,0)));
              helpmat(i,0) += 0.5*tsteps(k,0)*interactvar(i,fc_pos[xcols+j])*
                              (t_Z(k,dm_pos[xcols+j])*baseline(i,k)
                               +t_Z(k+1,dm_pos[xcols+j])*baseline(i,k+1));
              }
            }
          }
        for(i=0; i<nrobs; i++)
          {
          helpmat(i,0) *= mult_hazard(i,0);
          }
        H1(X.cols()+j,0)=(resp.transposed()*Z.getCol(j))(0,0)-helpmat.sum(0)-Qinv(j,0)*beta(X.cols()+j,0);
        }
      else
        {
        H1(X.cols()+j,0)=(resp.transposed()*Z.getCol(j))(0,0)-(cumhazard.transposed()*Z.getCol(j))(0,0)-Qinv(j,0)*beta(X.cols()+j,0);
        }
      }

  // Fisher-Info für beta

  // X & X
    for(j=0; j<X.cols(); j++)
      {
      for(k=j; k<X.cols(); k++)
        {
        if(isbaselinebeta[j]==1 && isbaselinebeta[k]==1)
          {
          helpmat = datamatrix(nrobs,1,0);
          //CHECK
          if(check)
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
                helpmat(i,0) += 0.5*tsteps(l,0)*(t_X(l,dm_pos[j])*t_X(l,dm_pos[k])*baseline(l,0)
                                +t_X(l+1,dm_pos[j])*t_X(l+1,dm_pos[k])*baseline(l+1,0));
                }
              }
            }
          else
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
//                helpmat(i,0) += 0.5*tstep*interactvar(i,fc_pos[j])*interactvar(i,fc_pos[k])*
//                                (t_X(l,dm_pos[j])*t_X(l,dm_pos[k])*exp((interactvar.getRow(i).transposed()*basef.getRow(l))(0,0))
//                                 +t_X(l+1,dm_pos[j])*t_X(l+1,dm_pos[k])*exp((interactvar.getRow(i).transposed()*basef.getRow(l+1))(0,0)));
                helpmat(i,0) += 0.5*tsteps(l,0)*interactvar(i,fc_pos[j])*interactvar(i,fc_pos[k])*
                                (t_X(l,dm_pos[j])*t_X(l,dm_pos[k])*baseline(i,l)
                                 +t_X(l+1,dm_pos[j])*t_X(l+1,dm_pos[k])*baseline(i,l+1));
                }
              }
            }
          for(i=0; i<nrobs; i++)
            {
            helpmat(i,0) *= mult_hazard(i,0);
            }
          H(j,k)=helpmat.sum(0);
          }
        else if(isbaselinebeta[j]==0 && isbaselinebeta[k]==1)
          {
          helpmat = datamatrix(nrobs,1,0);
          //CHECK
          if(check)
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
                helpmat(i,0) += 0.5*tsteps(l,0)*(t_X(l,dm_pos[k])*baseline(l,0)
                                 +t_X(l+1,dm_pos[k])*baseline(l+1,0));
                }
              }
            }
          else
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
//                helpmat(i,0) += 0.5*tstep*interactvar(i,fc_pos[k])*
//                                (t_X(l,dm_pos[k])*exp((interactvar.getRow(i).transposed()*basef.getRow(l))(0,0))
//                                 +t_X(l+1,dm_pos[k])*exp((interactvar.getRow(i).transposed()*basef.getRow(l+1))(0,0)));
                helpmat(i,0) += 0.5*tsteps(l,0)*interactvar(i,fc_pos[k])*
                                (t_X(l,dm_pos[k])*baseline(i,l)
                                 +t_X(l+1,dm_pos[k])*baseline(i,l+1));
                }
              }
            }
          for(i=0; i<nrobs; i++)
            {
            helpmat(i,0) *= mult_hazard(i,0);
            }
          H(j,k)=(helpmat.transposed()*X.getCol(j))(0,0);
          }
        else if(isbaselinebeta[j]==1 && isbaselinebeta[k]==0)
          {
          helpmat = datamatrix(nrobs,1,0);
          //CHECK
          if(check)
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
                helpmat(i,0) += 0.5*tsteps(l,0)*(t_X(l,dm_pos[j])*baseline(l,0)
                                 +t_X(l+1,dm_pos[j])*baseline(l+1,0));
                }
              }
            }
          else
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
//                helpmat(i,0) += 0.5*tstep*interactvar(i,fc_pos[j])*
//                                (t_X(l,dm_pos[j])*exp((interactvar.getRow(i).transposed()*basef.getRow(l))(0,0))
//                                 +t_X(l+1,dm_pos[j])*exp((interactvar.getRow(i).transposed()*basef.getRow(l+1))(0,0)));
                helpmat(i,0) += 0.5*tsteps(l,0)*interactvar(i,fc_pos[j])*
                                (t_X(l,dm_pos[j])*baseline(i,l)
                                 +t_X(l+1,dm_pos[j])*baseline(i,l+1));
                }
              }
            }
          for(i=0; i<nrobs; i++)
            {
            helpmat(i,0) *= mult_hazard(i,0);
            }
          H(j,k)=(helpmat.transposed()*X.getCol(k))(0,0);
          }
        else
          {
          for(i=0; i<nrobs; i++)
            {
            helpmat(i,0)=X(i,j)*X(i,k);
            }
          H(j,k)=(helpmat.transposed()*cumhazard)(0,0);
          }
        H(k,j)=H(j,k);
        }
      }

  // Z & Z
    for(j=0; j<zcols; j++)
      {
      for(k=j; k<zcols; k++)
        {
        if(isbaselinebeta[xcols+j]==1 && isbaselinebeta[xcols+k]==1)
          {
          helpmat = datamatrix(nrobs,1,0);
          //CHECK
          if(check)
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
                helpmat(i,0) += 0.5*tsteps(l,0)*(t_Z(l,dm_pos[xcols+j])*t_Z(l,dm_pos[xcols+k])*baseline(l,0)
                                 +t_Z(l+1,dm_pos[xcols+j])*t_Z(l+1,dm_pos[xcols+k])*baseline(l+1,0));
                }
              }
            }
          else
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
//                helpmat(i,0) += 0.5*tstep*interactvar(i,fc_pos[xcols+j])*interactvar(i,fc_pos[xcols+k])*
//                                (t_Z(l,dm_pos[xcols+j])*t_Z(l,dm_pos[xcols+k])*exp((interactvar.getRow(i).transposed()*basef.getRow(l))(0,0))
//                                 +t_Z(l+1,dm_pos[xcols+j])*t_Z(l+1,dm_pos[xcols+k])*exp((interactvar.getRow(i).transposed()*basef.getRow(l+1))(0,0)));
                helpmat(i,0) += 0.5*tsteps(l,0)*interactvar(i,fc_pos[xcols+j])*interactvar(i,fc_pos[xcols+k])*
                                (t_Z(l,dm_pos[xcols+j])*t_Z(l,dm_pos[xcols+k])*baseline(i,l)
                                 +t_Z(l+1,dm_pos[xcols+j])*t_Z(l+1,dm_pos[xcols+k])*baseline(i,l+1));
                }
              }
            }
          for(i=0; i<nrobs; i++)
            {
            helpmat(i,0) *= mult_hazard(i,0);
            }
          H(xcols+j,xcols+k)=helpmat.sum(0);
          }
        else if(isbaselinebeta[xcols+j]==0 && isbaselinebeta[xcols+k]==1)
          {
          helpmat = datamatrix(nrobs,1,0);
          //CHECK
          if(check)
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
                helpmat(i,0) += 0.5*tsteps(l,0)*(t_Z(l,dm_pos[xcols+k])*baseline(l,0)
                                 +t_Z(l+1,dm_pos[xcols+k])*baseline(l+1,0));
                }
              }
            }
          else
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
//                helpmat(i,0) += 0.5*tstep*interactvar(i,fc_pos[xcols+k])*
//                                (t_Z(l,dm_pos[xcols+k])*exp((interactvar.getRow(i).transposed()*basef.getRow(l))(0,0))
//                                 +t_Z(l+1,dm_pos[xcols+k])*exp((interactvar.getRow(i).transposed()*basef.getRow(l+1))(0,0)));
                helpmat(i,0) += 0.5*tsteps(l,0)*interactvar(i,fc_pos[xcols+k])*
                                (t_Z(l,dm_pos[xcols+k])*baseline(i,l)
                                 +t_Z(l+1,dm_pos[xcols+k])*baseline(i,l+1));
                }
              }
            }
          for(i=0; i<nrobs; i++)
            {
            helpmat(i,0) *= mult_hazard(i,0);
            }
          H(xcols+j,xcols+k)=(helpmat.transposed()*Z.getCol(j))(0,0);
          }
        else if(isbaselinebeta[xcols+j]==1 && isbaselinebeta[xcols+k]==0)
          {
          helpmat = datamatrix(nrobs,1,0);
          //CHECK
          if(check)
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
                helpmat(i,0) += 0.5*tsteps(l,0)*(t_Z(l,dm_pos[xcols+j])*baseline(l,0)
                                 +t_Z(l+1,dm_pos[xcols+j])*baseline(l+1,0));
                }
              }
            }
          else
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
//                helpmat(i,0) += 0.5*tstep*interactvar(i,fc_pos[xcols+j])*
//                                (t_Z(l,dm_pos[xcols+j])*exp((interactvar.getRow(i).transposed()*basef.getRow(l))(0,0))
//                                 +t_Z(l+1,dm_pos[xcols+j])*exp((interactvar.getRow(i).transposed()*basef.getRow(l+1))(0,0)));
                helpmat(i,0) += 0.5*tsteps(l,0)*interactvar(i,fc_pos[xcols+j])*
                                (t_Z(l,dm_pos[xcols+j])*baseline(i,l)
                                 +t_Z(l+1,dm_pos[xcols+j])*baseline(i,l+1));
                }
              }
            }
          for(i=0; i<nrobs; i++)
            {
            helpmat(i,0) *= mult_hazard(i,0);
            }
          H(xcols+j,xcols+k)=(helpmat.transposed()*Z.getCol(k))(0,0);
          }
        else
          {
          for(i=0; i<nrobs; i++)
            {
            helpmat(i,0)=Z(i,j)*Z(i,k);
            }
          H(xcols+j,xcols+k)=(helpmat.transposed()*cumhazard)(0,0);
          }
        H(xcols+k,xcols+j)=H(xcols+j,xcols+k);
        }
      }

  // X & Z
    for(j=0; j<xcols; j++)
      {
      for(k=0; k<zcols; k++)
        {
        if(isbaselinebeta[j]==1 && isbaselinebeta[xcols+k]==1)
          {
          helpmat = datamatrix(nrobs,1,0);
          //CHECK
          if(check)
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
                helpmat(i,0) += 0.5*tsteps(l,0)*(t_X(l,dm_pos[j])*t_Z(l,dm_pos[xcols+k])*baseline(l,0)
                                 +t_X(l+1,dm_pos[j])*t_Z(l+1,dm_pos[xcols+k])*baseline(l+1,0));
                }
              }
            }
          else
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
//                helpmat(i,0) += 0.5*tstep*interactvar(i,fc_pos[j])*interactvar(i,fc_pos[xcols+k])*
//                                (t_X(l,dm_pos[j])*t_Z(l,dm_pos[xcols+k])*exp((interactvar.getRow(i).transposed()*basef.getRow(l))(0,0))
//                                 +t_X(l+1,dm_pos[j])*t_Z(l+1,dm_pos[xcols+k])*exp((interactvar.getRow(i).transposed()*basef.getRow(l+1))(0,0)));
                helpmat(i,0) += 0.5*tsteps(l,0)*interactvar(i,fc_pos[j])*interactvar(i,fc_pos[xcols+k])*
                                (t_X(l,dm_pos[j])*t_Z(l,dm_pos[xcols+k])*baseline(i,l)
                                 +t_X(l+1,dm_pos[j])*t_Z(l+1,dm_pos[xcols+k])*baseline(i,l+1));
                }
              }
            }
          for(i=0; i<nrobs; i++)
            {
            helpmat(i,0) *= mult_hazard(i,0);
            }
          H(j,xcols+k)=helpmat.sum(0);
          }
        else if(isbaselinebeta[j]==0 && isbaselinebeta[xcols+k]==1)
          {
          helpmat = datamatrix(nrobs,1,0);
          //CHECK
          if(check)
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
                helpmat(i,0) += 0.5*tsteps(l,0)*(t_Z(l,dm_pos[xcols+k])*baseline(l,0)
                                 +t_Z(l+1,dm_pos[xcols+k])*baseline(l+1,0));
                }
              }
            }
          else
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
//                helpmat(i,0) += 0.5*tstep*interactvar(i,fc_pos[xcols+k])*
//                                (t_Z(l,dm_pos[xcols+k])*exp((interactvar.getRow(i).transposed()*basef.getRow(l))(0,0))
//                                 +t_Z(l+1,dm_pos[xcols+k])*exp((interactvar.getRow(i).transposed()*basef.getRow(l+1))(0,0)));
                helpmat(i,0) += 0.5*tsteps(l,0)*interactvar(i,fc_pos[xcols+k])*
                                (t_Z(l,dm_pos[xcols+k])*baseline(i,l)
                                 +t_Z(l+1,dm_pos[xcols+k])*baseline(i,l+1));
                }
              }
            }
          for(i=0; i<nrobs; i++)
            {
            helpmat(i,0) *= mult_hazard(i,0);
            }
          H(j,xcols+k)=(helpmat.transposed()*X.getCol(j))(0,0);
          }
        else if(isbaselinebeta[j]==1 && isbaselinebeta[xcols+k]==0)
          {
          helpmat = datamatrix(nrobs,1,0);
          //CHECK
          if(check)
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
                helpmat(i,0) += 0.5*tsteps(l,0)*(t_X(l,dm_pos[j])*baseline(l,0)
                                 +t_X(l+1,dm_pos[j])*baseline(l+1,0));
                }
              }
            }
          else
            {
            for(i=0; i<nrobs; i++)
              {
              for(l=tstart[i]; l<tend[i]-1; l++)
                {
//                helpmat(i,0) += 0.5*tstep*interactvar(i,fc_pos[j])*
//                                (t_X(l,dm_pos[j])*exp((interactvar.getRow(i).transposed()*basef.getRow(l))(0,0))
//                                 +t_X(l+1,dm_pos[j])*exp((interactvar.getRow(i).transposed()*basef.getRow(l+1))(0,0)));
                helpmat(i,0) += 0.5*tsteps(l,0)*interactvar(i,fc_pos[j])*
                                (t_X(l,dm_pos[j])*baseline(i,l)
                                 +t_X(l+1,dm_pos[j])*baseline(i,l+1));
                }
              }
            }
          for(i=0; i<nrobs; i++)
            {
            helpmat(i,0) *= mult_hazard(i,0);
            }
          H(j,xcols+k)=(helpmat.transposed()*Z.getCol(k))(0,0);
          }
        else
          {
          for(i=0; i<nrobs; i++)
            {
            helpmat(i,0)=X(i,j)*Z(i,k);
            }
          H(j,xcols+k)=(helpmat.transposed()*cumhazard)(0,0);
          }
        H(xcols+k,j)=H(j,xcols+k);
        }
      }

    H.addtodiag(Qinv,xcols,beta.rows());

    // Fisher-scoring für beta
    beta = betaold + H.solve(H1);

    stop = check_pause();
    if (stop)
      return true;

  //////////////////////////////////////////////
  // Marginale Likelihood optimieren          //
  //////////////////////////////////////////////

    Hinv=H.inverse();

    // transform theta
    for(i=0; i<theta.rows(); i++)
      {
      thetaold(i,0)=signs[i]*sqrt(thetaold(i,0));
      theta(i,0)=signs[i]*sqrt(theta(i,0));
      }

    // Score-Funktion für theta

   for(j=0; j<theta.rows(); j++)
      {
      score(j,0)=-1*((zcut[j+1]-zcut[j])/theta(j,0)-
                       (Hinv.getBlock(X.cols()+zcut[j],X.cols()+zcut[j],X.cols()+zcut[j+1],X.cols()+zcut[j+1])).trace()/(theta(j,0)*theta(j,0)*theta(j,0))-
                       (beta.getRowBlock(X.cols()+zcut[j],X.cols()+zcut[j+1]).transposed()*beta.getRowBlock(X.cols()+zcut[j],X.cols()+zcut[j+1]))(0,0)/(theta(j,0)*theta(j,0)*theta(j,0)));
      }

    // Fisher-Info für theta

    for(j=0; j<theta.rows(); j++)
      {
      for(k=j; k< theta.rows(); k++)
        {
        Fisher(j,k) = 2*((Hinv.getBlock(X.cols()+zcut[j],X.cols()+zcut[k],X.cols()+zcut[j+1],X.cols()+zcut[k+1])*Hinv.getBlock(X.cols()+zcut[k],X.cols()+zcut[j],X.cols()+zcut[k+1],X.cols()+zcut[j+1])).trace())/(theta(j,0)*theta(j,0)*theta(j,0)*theta(k,0)*theta(k,0)*theta(k,0));
        Fisher(k,j) = Fisher(j,k);
        }
      }

    //Fisher-scoring für theta

    theta = thetaold + Fisher.solve(score);

    // transform theta back to original parameterisation

    for(i=0; i<theta.rows(); i++)
      {
      signs[i] = -1*(theta(i,0)<0)+1*(theta(i,0)>=0);
      theta(i,0) *= theta(i,0);
      thetaold(i,0) *= thetaold(i,0);
      }

    // update linear predictor
    eta=X*beta.getRowBlock(0,xcols)+Z*beta.getRowBlock(xcols,beta.rows());

    // test whether to stop estimation of theta[i]
   help=eta.norm(0);
   for(i=0; i<theta.rows(); i++)
     {
     helpmat=Z.getColBlock(zcut[i],zcut[i+1])*beta.getRowBlock(X.cols()+zcut[i],X.cols()+zcut[i+1]);
     stopcrit[i]=helpmat.norm(0)/help;
     if(stopcrit[i]<lowerlim)
       {
       theta(i,0)=thetaold(i,0);
       }
     else
       {
       its[i]=it;
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

    // test criterion
    test=((crit1>eps) || (crit2>eps)) && (it<(unsigned)maxit);

    out("  iteration "+ST::inttostring(it)+"\n");
    out("  relative changes in the regression coefficients: "+
         ST::doubletostring(crit1,6)+"\n");
    out("  relative changes in the variance parameters:     "+
         ST::doubletostring(crit2,6)+"\n");
    out("\n");

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

  for(i=1;i<fullcond.size();i++)
    {
    beta(0,0) += fullcond[i]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[i],zcut[i-1],i-1,false,xcut[i],X.cols()+zcut[i-1],0,false,i);
    }
  beta(0,0) += fullcond[0]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[0],0,0,false,xcut[0],0,0,false,0);

  return false;
  }


//------------------------------------------------------------------------------
//----------------------------- Object description -----------------------------
//------------------------------------------------------------------------------

  void remlest::outoptions()
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
    ST::string familyname;
    ST::string respname;
    if(respfamily=="gaussian")
      {
      familyname="gaussian";
      respname="identity";
      }
    else if(respfamily=="gamma")
      {
      familyname="gamma";
      respname="exponential";
      }
    else if(respfamily=="poisson")
      {
      familyname="poisson";
      respname="exponential";
      }
    else if(respfamily=="poissondispers")
      {
      familyname="poisson (with overdispersion)";
      respname="exponential";
      }
    else if(respfamily=="binomial")
      {
      familyname="binomial";
      respname="logistic distribution function (logit link)";
      }
    else if(respfamily=="binomialdispers")
      {
      familyname="binomial (with overdispersion)";
      respname="logistic distribution function (logit link)";
      }
    else if(respfamily=="binomialprobit")
      {
      familyname="binomial";
      respname="standard normal (probit link)";
      }
    else if(respfamily=="binomialprobitdispers")
      {
      familyname="binomial (with overdispersion)";
      respname="standard normal (probit link)";
      }

    if(respfamily=="cox")
      {
      out("  Family:                 cox\n");
      out("  Number of observations: "+ST::inttostring(X.rows())+"\n");
      }
    else
      {
      out("  Family:                 "+familyname+"\n");
      out("  Response function:      "+respname+"\n");
      out("  Number of observations: "+ST::inttostring(X.rows())+"\n");
      }
    }

//------------------------------------------------------------------------------
//----------------------------- Writing results --------------------------------
//------------------------------------------------------------------------------

void remlest::make_plots(ofstream & outtex,ST::string path_batch,
                         ST::string path_splus)
  {

  char hcharu = '_';
  ST::string hstringu = "\\_";

  unsigned j;
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

    // durchlaufen der Fullconditionals
    for(j=0;j<fullcond.size();j++)
      {

      // Pfad der Regr.-Ergebnisse
      pathresult = fullcond[j]->get_pathresult();

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
          outtex << "\n\\begin{figure}[h!]" << endl
                  << "\\centering" << endl
                  << "\\includegraphics[scale=0.6]{" << pathgr << ".ps}" << endl
                  << "\\caption{Non--linear Effect of '" <<
                  xvar.insert_string_char(hcharu,hstringu) << "'";
          outtex << "." << endl << "Shown are the posterior modes together with "
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
          outtex << "\n\\begin{figure}[h!]" << endl
                 << "\\centering" << endl
                 << "\\includegraphics[scale=0.6]{" << pathgr << "_pmode.ps}"
                 << endl
                 << "\\caption{Non--linear Effect of '" <<
                 regionvar.insert_string_char(hcharu,hstringu) << "'";
          outtex << ". Shown are the posterior modes.}" << endl
                 << "\\end{figure}" << endl;
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
          } // end: else if
        } // end: if
      } // end: for
    }
  }

void remlest::make_model(ofstream & outtex, const ST::string & rname)
  {
    ST::string familyname;
    ST::string respname;
    if(respfamily=="gaussian")
      {
      familyname="gaussian";
      respname="identity";
      }
    else if(respfamily=="gamma")
      {
      familyname="gamma";
      respname="exponential";
      }
    else if(respfamily=="poisson")
      {
      familyname="poisson";
      respname="exponential";
      }
    else if(respfamily=="poissondispers")
      {
      familyname="poisson (with overdispersion)";
      respname="exponential";
      }
    else if(respfamily=="binomial")
      {
      familyname="binomial";
      respname="logistic distribution function (logit link)";
      }
    else if(respfamily=="binomialdispers")
      {
      familyname="binomial (with overdispersion)";
      respname="logistic distribution function (logit link)";
      }
    else if(respfamily=="binomialprobit")
      {
      familyname="binomial";
      respname="standard normal (probit link)";
      }
    else if(respfamily=="binomialprobitdispers")
      {
      familyname="binomial (with overdispersion)";
      respname="standard normal (probit link)";
      }
    else if(respfamily=="cox")
      {
      familyname="cox";
      }

  //Anz. Beob. wird übergeben
  unsigned obs = X.rows();

  //schreibt das Modell und die Priori-annahmen ins Tex-File
  outtex << "\n\\noindent {\\bf \\large Response:}" << endl
         << "\\begin{tabbing}\n"
         << "Number of observations: \\= " << obs << "\\\\" << endl
         << "Response Variable: \\> " << rname << "\\\\" << endl
         << "Family: \\> " << familyname << "\\\\" << endl;

  if(respfamily!="cox")
    outtex << "Response function: \\> " << respname << "\\\\" << endl;

  outtex << "\\end{tabbing}" << endl
         << "\n\\noindent {\\bf \\large Predictor:}\\\\" << endl;

  }

void remlest::make_predictor(ofstream & outtex)
  {

  unsigned j;

  ST::string term2 = fullcond[0]->get_term_symbolic();
  ST::string term = "$\\eta$ & $=$ & $" + term2;    //linearer Prädiktor wird erweitert
  for(j=1;j<fullcond.size();j++)
    {
    term2 = fullcond[j]->get_term_symbolic();
    term = term + " + " + term2;    //linearer Prädiktor wird erweitert
    }

  outtex << endl << "\n\\begin{tabular}{ccp{12cm}}\n" << term
         << "$\n\\end{tabular}\n\\\\ \n\\\\" << endl;
  }

void remlest::make_prior(ofstream & outtex)
  {
  unsigned i,j;
  outtex << "\n\\noindent {\\bf \\large Priors:}\\\\" << endl << "\\\\" << endl;
  for(j=0;j<fullcond.size();j++)
    {
    vector<ST::string> prior = fullcond[j]->get_priorassumptions();
    if(prior.size() != 0)// nur wenn Priors da sind (d.h. Vektor hat Elemente)
      {
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

void remlest::make_options(ofstream & outtex)
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

void remlest::make_fixed_table(ofstream & outtex)
  {

  // falls andere Quantile gewünscht werden
  double u = fullcond[0]->get_level1();
  ST::string u_str = ST::doubletostring(u,0);

  vector<ST::string> h;

  unsigned j;
  unsigned r;

  r = 2;

  // Tabelle im Tex-File mit fixen Effekten
  outtex << "\n\\newpage \n" << endl << "\n\\noindent {\\bf \\large Fixed Effects:}\\\\"
       << endl << "\\\\" << endl;

  outtex << "\\begin{tabular}{|r|rrrrr|}" << endl << "\\hline" << endl
         << "Variable & Post. Mode & Std. Dev. & p-value & \\multicolumn{2}{r|}{" << u << "\\% confidence interval}\\\\"
         << endl << "\\hline" << endl;

  h = fullcond[0]->get_results_latex();
  for (j=0;j<h.size();j++)
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
             << "\n\\noindent {\\bf \\large Fixed Effects (continued):}\\\\"
             << endl << "\\\\" << endl;

      outtex << "\\begin{tabular}{|r|rrrrr|}" << endl << "\\hline" << endl
             << "Variable & Post. Mode & Std. Dev. & p-value & \\multicolumn{2}{r|}{" << u << "\\% confidence interval}\\\\"
             << endl << "\\hline" << endl;

      outtex << h[j] << endl;
      }
    }
  outtex << "\\hline \n\\end{tabular}" << endl;
  }

void remlest::make_graphics(const ST::string & title,
                     const ST::string & path_batch,
                     const ST::string & path_tex,
                     const ST::string & path_splus,
                     const ST::string & rname,
                     const bool & dispers)
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

  if(dispers==true)
    {
    outtex<<"{\\bf \\large Estimated scale parameter:} " << theta(theta.rows()-1,0) <<"\\\\"<<endl;
    }

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

bool remlest::check_pause()
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

void remlest::out(const ST::string & s,bool thick,bool italic,
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


void remlest::outerror(const ST::string & s)
  {
  out(s,true,true,12,255,0,0);
  }


