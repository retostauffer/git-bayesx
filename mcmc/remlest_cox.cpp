#if defined(BORLAND_OUTPUT_WINDOW)
#include "StatResults.h"
#include "statwinframe.h"
#endif

#include<remlest.h>
#include<remlest_cox.h>

//------------------------------------------------------------------------------
//---------- Survival data with interval censoring & lef truncation-------------
//------------------------------------------------------------------------------

bool remlest::estimate_survival_interval2(datamatrix resp,
                const datamatrix & offset, const datamatrix & weight)
  {

  unsigned i, j, k, l;
  double help, former, helpint;

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
  unsigned nrobs = Z.rows();
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

  // Inzidenzmatrix, die für jeden Eintrag in fullcond bzw. beta angibt, ob er zur Baseline-HR beiträgt
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
  vector<int>dmat_pos(beta.rows(),0);
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
  l=0;
  for(i=0; i<beta.rows(); i++)
    {
    if(isbaselinebeta[i]==1)
      {
      dmat_pos[i]=l;
      l++;
      }
    }

  bool timevarying;
  if(nrbaseline>1)
    {
    timevarying=true;
    }
  else
    {
    timevarying=false;
    }

// Matrices and variables for baseline effects
  statmatrix<double> tsteps;
  datamatrix t_X;
  datamatrix t_Z;
  vector<unsigned> tleft;
  vector<unsigned> tright;
  vector<unsigned> ttrunc(nrobs,0);
  datamatrix interactvar(nrobs,nrbaseline,0);
  statmatrix<int> index;
  j=0;
  for(i=0; i<fullcond.size(); i++)
    {
    if(isbaseline[i]==1)
      {
      fullcond[i]->initialize_baseline(j,t_X,t_Z,tleft,tright,interactvar,tsteps,index);
      j++;
      }
    }
// first derivative of the cumulated baseline hazard
  datamatrix Dmat(t_X.rows(),t_X.cols()+t_Z.cols(),0);

// time-varying effects. the first row corresponds to the log-baseline
  datamatrix basef(t_X.rows(),nrbaseline,0);
  statmatrix<double>baseline = datamatrix(t_X.rows(),1,0);

  statmatrix<double>cumbaseline(t_X.rows(),1,0);
  statmatrix<double>negcumbaseline(t_X.rows(),1,0);
  statmatrix<double>cumhazard(nrobs,1,0);
  statmatrix<double>negcumhazard(nrobs,1,0);
  statmatrix<double>eta(nrobs,1,0);
  statmatrix<double>baseline_eta(nrobs,1,0);
  statmatrix<double>mult_eta(nrobs,1,0);
  statmatrix<double>mult_hazard(nrobs,1,0);
  statmatrix<double>helpmat(nrobs,1,0);
// indicator for interval or left censoring
  vector<bool>interval(nrobs,false);
  for(i=0; i<nrobs; i++)
    {
    if(resp(i,0)==-1)
      {
      interval[i] = true;
      resp(i,0)=0;
      }
    }

  // Transform smoothing paramater starting values to variances
  for(i=0; i<theta.rows(); i++)
    {
    theta(i,0)=1/theta(i,0);
    }

  beta(0,0) = log(10/t_X(t_X.rows()-1,1));

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

    // compute baseline

    for(i=0; i<t_X.rows(); i++)
      {
      baseline(i,0) = exp(basef(i,0));
      }

    // compute cumulated baseline

    former=0;
    cumbaseline = datamatrix(cumbaseline.rows(),1,0);
    for(i=0; i<t_X.rows()-1; i++)
      {
      cumbaseline(i,0) = former + 0.5*tsteps(i,0)*(baseline(i,0)+baseline(i+1,0));
      former = cumbaseline(i,0);
      }
    cumbaseline(t_X.rows()-1,0) = former;
    negcumbaseline = -cumbaseline;

    // compute mult_hazard = exp(x'beta) without time-varying covariates x(t)

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
    for(i=0; i<nrobs; i++)
      {
      mult_hazard(i,0)=exp(mult_eta(i,0));
      }

    // compute cumulated hazard

    for(i=0; i<nrobs; i++)
      {
      cumhazard(i,0)=cumbaseline(tright[i]-1,0)*mult_hazard(i,0);
      }
    negcumhazard = - cumhazard;

    // compute derivative matrix D

    Dmat = datamatrix(Dmat.rows(),Dmat.cols(),0);
    for(j=0; j<xcols; j++)
      {
      if(isbaselinebeta[j]==1)
        {
        former=0;
        for(i=0; i<t_X.rows()-1; i++)
          {
          Dmat(i,dmat_pos[j]) = former - 0.5*tsteps(i,0)*(t_X(i,dm_pos[j])*baseline(i,0)+t_X(i+1,dm_pos[j])*baseline(i+1,0));
          former = Dmat(i,dmat_pos[j]);
          }
        Dmat(t_X.rows()-1,dmat_pos[j]) = former;
        }
      }
    for(j=0; j<zcols; j++)
      {
      if(isbaselinebeta[xcols+j]==1)
        {
        former=0;
        for(i=0; i<t_Z.rows()-1; i++)
          {
          Dmat(i,dmat_pos[xcols+j]) = former - 0.5*tsteps(i,0)*(t_Z(i,dm_pos[xcols+j])*baseline(i,0)+t_Z(i+1,dm_pos[xcols+j])*baseline(i+1,0));
          former = Dmat(i,dmat_pos[xcols+j]);
          }
        Dmat(t_Z.rows()-1,dmat_pos[xcols+j]) = former;
        }
      }

    // to be deleted
    datamatrix Survivor(nrobs,2,0);
    for(i=0; i<nrobs; i++)
      {
      if(interval[i])
        {
        Survivor(i,0) = exp(-cumbaseline(tleft[i],0)*mult_hazard(i,0));
        Survivor(i,1) = exp(-cumbaseline(tright[i]-1,0)*mult_hazard(i,0));
        }
      }

    // Score-Funktion für beta

    // X

    for(j=0; j<xcols; j++)
      {
      H1(j,0)=(resp.transposed()*X.getCol(j))(0,0);

      // x_j gehört zu Baseline
      if(isbaselinebeta[j]==1)
        {
        for(i=0; i<nrobs; i++)
          {
          if(interval[i])
            {
            helpint = exp(-(cumbaseline(tright[i]-1,0) - cumbaseline(tleft[i],0)) * mult_hazard(i,0));
            H1(j,0) += Dmat(tleft[i],dmat_pos[j])*mult_hazard(i,0) -
                       helpint*(Dmat(tright[i]-1,dmat_pos[j]) - Dmat(tleft[i],dmat_pos[j]))*mult_hazard(i,0)/(1-helpint);
            }
          else
            {
            H1(j,0) += Dmat(tright[i]-1,dmat_pos[j])*mult_hazard(i,0);
            }
          if(ttrunc[i] > 0)
            {
            H1(j,0) -= Dmat(ttrunc[i],dmat_pos[j])*mult_hazard(i,0);
            }
          }
        }
      // x_j gehört nicht zur Baseline
      else
        {
        for(i=0; i<nrobs; i++)
          {
          if(interval[i])
            {
            helpint = (cumbaseline(tright[i]-1,0) - cumbaseline(tleft[i],0)) * mult_hazard(i,0);
            H1(j,0) += negcumbaseline(tleft[i],0)*X(i,j)*mult_hazard(i,0) +
                       exp(-helpint)*X(i,j)*helpint/(1-exp(-helpint));
            }
          else
            {
            H1(j,0) += negcumhazard(i,0)*X(i,j);
            }
          if(ttrunc[i] > 0)
            {
            H1(j,0) -= cumbaseline(ttrunc[i],0)*X(i,j)*mult_hazard(i,0);
            }
          }
        }
      }

    // Z

    for(j=0; j<zcols; j++)
      {
      H1(xcols+j,0)=(resp.transposed()*Z.getCol(j))(0,0);

      // z_j gehört zu Baseline
      if(isbaselinebeta[xcols+j]==1)
        {
        for(i=0; i<nrobs; i++)
          {
          if(interval[i])
            {
            helpint = exp(-(cumbaseline(tright[i]-1,0) - cumbaseline(tleft[i],0)) * mult_hazard(i,0));
            H1(xcols+j,0) += Dmat(tleft[i],dmat_pos[xcols+j])*mult_hazard(i,0) -
                       helpint*(Dmat(tright[i]-1,dmat_pos[xcols+j]) - Dmat(tleft[i],dmat_pos[xcols+j]))*mult_hazard(i,0)/(1-helpint);
            }
          else
            {
            H1(xcols + j,0) += Dmat(tright[i]-1,dmat_pos[xcols+j])*mult_hazard(i,0);
            }
          if(ttrunc[i] > 0)
            {
            H1(xcols+j,0) -= Dmat(ttrunc[i],dmat_pos[xcols+j])*mult_hazard(i,0);
            }
          }
        }
      // z_j gehört nicht zur Baseline
      else
        {
        for(i=0; i<nrobs; i++)
          {
          if(interval[i])
            {
            helpint = (cumbaseline(tright[i]-1,0) - cumbaseline(tleft[i],0)) * mult_hazard(i,0);
            H1(xcols+j,0) += negcumbaseline(tleft[i],0)*Z(i,j)*mult_hazard(i,0) +
                       exp(-helpint)*Z(i,j)*helpint/(1-exp(-helpint));
            }
          else
            {
            H1(xcols + j,0) += negcumhazard(i,0)*Z(i,j);
            }
          if(ttrunc[i] > 0)
            {
            H1(xcols+j,0) -= cumbaseline(ttrunc[i],0)*Z(i,j)*mult_hazard(i,0);
            }
          }
        }
      }

    for(j=0; j<zcols;j++)
      {
      H1(xcols+j,0) -= Qinv(j,0)*beta(xcols+j,0);
      }

    // Fisher-Information for beta

    //clear H-matrix
    H = datamatrix(H.rows(),H.cols(),0);

    // X & X
    for(j=0; j<xcols; j++)
      {
      for(k=j; k<xcols; k++)
        {
        // helpmat contains the second derivatives of the cumulated baseline
        // with respect to beta_j and beta_k
        helpmat = datamatrix(t_X.rows(),1,0);

        if(isbaselinebeta[j]==1 && isbaselinebeta[k]==1)
          {
          former=0;
          for(i=0; i<t_X.rows()-1; i++)
            {
            helpmat(i,0) = former - 0.5*tsteps(i,0)*
                                    (t_X(i,dm_pos[j])*t_X(i,dm_pos[k])*baseline(i,0)
                                    + t_X(i+1,dm_pos[j])*t_X(i+1,dm_pos[k])*baseline(i+1,0));
            former=helpmat(i,0);
            }
          helpmat(t_X.rows()-1,0) = former;
          for(i=0; i<nrobs; i++)
            {
            if(interval[i])
              {
              H(j,k) += (
                         (
                          helpmat(tleft[i],0) + Dmat(tleft[i],dmat_pos[j])*Dmat(tleft[i],dmat_pos[k])*mult_hazard(i,0)
                         )*Survivor(i,0)
                        -
                         (
                          helpmat(tright[i]-1,0) + Dmat(tright[i]-1,dmat_pos[j])*Dmat(tright[i]-1,dmat_pos[k])*mult_hazard(i,0)
                         )*Survivor(i,1)
                        ) * mult_hazard(i,0) / (Survivor(i,0)-Survivor(i,1))
                        -
                        (
                         (
                          Dmat(tleft[i],dmat_pos[j])*Survivor(i,0)-Dmat(tright[i]-1,dmat_pos[j])*Survivor(i,1)
                         ) *mult_hazard(i,0)
                        *
                         (
                         Dmat(tleft[i],dmat_pos[k])*Survivor(i,0)-Dmat(tright[i]-1,dmat_pos[k])*Survivor(i,1)
                         ) *mult_hazard(i,0)
                        )
                        /
                        (
                         (Survivor(i,0)-Survivor(i,1))
                         *
                         (Survivor(i,0)-Survivor(i,1))
                        );
              }
            else
              {
              H(j,k) += helpmat(tright[i]-1,0)*mult_hazard(i,0);
              }
            }
          }

        else if(isbaselinebeta[j]==0 && isbaselinebeta[k]==1)
          {
          former=0;
          for(i=0; i<t_X.rows()-1; i++)
            {
            helpmat(i,0) = former - 0.5*tsteps(i,0)*(t_X(i,dm_pos[k])*baseline(i,0)
                                +t_X(i+1,dm_pos[k])*baseline(i+1,0));
            former = helpmat(i,0);
            }
          helpmat(t_X.rows()-1,0) = former;
          for(i=0; i<nrobs; i++)
            {
            if(interval[i])
              {
              H(j,k) += (
                         (
                          helpmat(tleft[i],0) + negcumbaseline(tleft[i],0)*Dmat(tleft[i],dmat_pos[k])*mult_hazard(i,0)
                         ) * Survivor(i,0)
                         -
                         (
                          helpmat(tright[i]-1,0) + negcumbaseline(tright[i]-1,0)*Dmat(tright[i]-1,dmat_pos[k])*mult_hazard(i,0)
                         ) * Survivor(i,1)
                        ) *mult_hazard(i,0) * X(i,j) / (Survivor(i,0)-Survivor(i,1))
                        -
                        (
                         (
                          negcumbaseline(tleft[i],0)*Survivor(i,0) - negcumbaseline(tright[i]-1,0)*Survivor(i,1)
                         ) * mult_hazard(i,0) * X(i,j)
                         *
                         (
                          Dmat(tleft[i],dmat_pos[k])*Survivor(i,0) - Dmat(tright[i]-1,dmat_pos[k])*Survivor(i,1)
                         ) * mult_hazard(i,0)
                        )
                        /
                        (
                         (Survivor(i,0)-Survivor(i,1))
                         *
                         (Survivor(i,0)-Survivor(i,1))
                        );
              }
            else
              {
              H(j,k) += helpmat(tright[i]-1,0)*X(i,j)*mult_hazard(i,0);
              }
            }
          }

        else if(isbaselinebeta[j]==1 && isbaselinebeta[k]==0)
          {
          former=0;
          for(i=0; i<t_X.rows()-1; i++)
            {
            helpmat(i,0) = former - 0.5*tsteps(i,0)*(t_X(i,dm_pos[j])*baseline(i,0)
                                +t_X(i+1,dm_pos[j])*baseline(i+1,0));
            former = helpmat(i,0);
            }
          helpmat(t_X.rows()-1,0) = former;
          for(i=0; i<nrobs; i++)
            {
            if(interval[i])
              {
              H(j,k) += (
                         (
                          helpmat(tleft[i],0) + Dmat(tleft[i],dmat_pos[j])*negcumbaseline(tleft[i],0)*mult_hazard(i,0)
                         ) * Survivor(i,0)
                         -
                         (
                          helpmat(tright[i]-1,0) + Dmat(tright[i]-1,dmat_pos[j])*negcumbaseline(tright[i]-1,0)*mult_hazard(i,0)
                         ) * Survivor(i,1)
                        ) * mult_hazard(i,0) * X(i,k) / (Survivor(i,0)-Survivor(i,1))
                        -
                        (
                         (
                          Dmat(tleft[i],dmat_pos[j])*Survivor(i,0)-Dmat(tright[i]-1,dmat_pos[j])*Survivor(i,1)
                         ) * mult_hazard(i,0)
                         *
                         (
                          negcumbaseline(tleft[i],0)*Survivor(i,0)-negcumbaseline(tright[i]-1,0)*Survivor(i,1)
                         ) * mult_hazard(i,0) * X(i,k)
                        )
                        /
                        (
                         (Survivor(i,0)-Survivor(i,1)) * (Survivor(i,0)-Survivor(i,1))
                        );
              }
            else
              {
              H(j,k) += helpmat(tright[i]-1,0)*X(i,k)*mult_hazard(i,0);
              }
            }
          }

        else
          {
          for(i=0; i<nrobs; i++)
            {
            if(interval[i])
              {
              H(j,k) += (
                         (
                          negcumbaseline(tleft[i],0) + negcumbaseline(tleft[i],0)*negcumbaseline(tleft[i],0)*mult_hazard(i,0)
                         ) * Survivor(i,0)
                         -
                         (
                          negcumbaseline(tright[i]-1,0) + negcumbaseline(tright[i]-1,0)*negcumbaseline(tright[i]-1,0)*mult_hazard(i,0)
                         ) * Survivor(i,1)
                        ) * mult_hazard(i,0) * X(i,j) * X(i,k) / (Survivor(i,0)-Survivor(i,1))
                        -
                        (
                         (
                          negcumbaseline(tleft[i],0)*Survivor(i,0) - negcumbaseline(tright[i]-1,0)*Survivor(i,1)
                         )*mult_hazard(i,0)
                         *
                         (
                          negcumbaseline(tleft[i],0)*Survivor(i,0) - negcumbaseline(tright[i]-1,0)*Survivor(i,1)
                         ) * mult_hazard(i,0) * X(i,j) * X(i,k)
                        )
                        /
                        (
                         (Survivor(i,0)-Survivor(i,1)) * (Survivor(i,0)-Survivor(i,1))
                        );
              }
            else
              {
              H(j,k) += negcumbaseline(tright[i]-1,0)*X(i,j)*X(i,k)*mult_hazard(i,0);
              }
            }
          }
        H(k,j)=H(j,k);
        }
      }

    // X & Z
    for(j=0; j<xcols; j++)
      {
      for(k=0; k<zcols; k++)
        {
        // helpmat contains the second derivatives of the cumulated baseline
        // with respect to beta_j and beta_k
        helpmat = datamatrix(t_X.rows(),1,0);

        if(isbaselinebeta[j]==1 && isbaselinebeta[xcols+k]==1)
          {
          former=0;
          for(i=0; i<t_X.rows()-1; i++)
            {
            helpmat(i,0) = former - 0.5*tsteps(i,0)*(t_X(i,dm_pos[j])*t_Z(i,dm_pos[xcols+k])*baseline(i,0)
                                +t_X(i+1,dm_pos[j])*t_Z(i+1,dm_pos[xcols+k])*baseline(i+1,0));
            former = helpmat(i,0);
            }
          helpmat(t_X.rows()-1,0) = former;
          for(i=0; i<nrobs; i++)
            {
            if(interval[i])
              {
              H(j,xcols+k) += (
                               (
                                helpmat(tleft[i],0) + Dmat(tleft[i],dmat_pos[j])*Dmat(tleft[i],dmat_pos[xcols+k])*mult_hazard(i,0)
                               ) * Survivor(i,0)
                               -
                              (
                               helpmat(tright[i]-1,0) + Dmat(tright[i]-1,dmat_pos[j])*Dmat(tright[i]-1,dmat_pos[xcols+k])*mult_hazard(i,0)
                              ) * Survivor(i,1)
                             ) * mult_hazard(i,0) / (Survivor(i,0)-Survivor(i,1))
                             -
                             (
                              (
                               Dmat(tleft[i],dmat_pos[j])*Survivor(i,0) - Dmat(tright[i]-1,dmat_pos[j])*Survivor(i,1)
                              ) * mult_hazard(i,0)
                              *
                             (
                              Dmat(tleft[i],dmat_pos[xcols+k])*Survivor(i,0) - Dmat(tright[i]-1,dmat_pos[xcols+k])*Survivor(i,1)

                             ) * mult_hazard(i,0)
                            )
                            /
                            (
                             (Survivor(i,0)-Survivor(i,1)) * (Survivor(i,0)-Survivor(i,1))
                            );
              }
            else
              {
              H(j,xcols+k) += helpmat(tright[i]-1,0)*mult_hazard(i,0);
              }
            }
          }

        else if(isbaselinebeta[j]==0 && isbaselinebeta[xcols+k]==1)
          {
          former=0;
          for(i=0; i<t_X.rows()-1; i++)
            {
            helpmat(i,0) = former - 0.5*tsteps(i,0)*(t_Z(i,dm_pos[xcols+k])*baseline(i,0)
                                +t_Z(i+1,dm_pos[xcols+k])*baseline(i+1,0));
            former = helpmat(i,0);
            }
          helpmat(t_X.rows()-1,0) = former;
          for(i=0; i<nrobs; i++)
            {
            if(interval[i])
              {
              H(j,xcols+k) += (
                               (
                                helpmat(tleft[i],0) + negcumbaseline(tleft[i],0)*Dmat(tleft[i],dmat_pos[xcols+k])*mult_hazard(i,0)
                               ) * Survivor(i,0)
                               -
                               (
                                helpmat(tright[i]-1,0) + negcumbaseline(tright[i]-1,0)*Dmat(tright[i]-1,dmat_pos[xcols+k])*mult_hazard(i,0)
                               ) * Survivor(i,1)
                              ) * mult_hazard(i,0) * X(i,j) / (Survivor(i,0)-Survivor(i,1))
                              -
                             (
                              (
                               negcumbaseline(tleft[i],0)*Survivor(i,0) - negcumbaseline(tright[i]-1,0)*Survivor(i,1)
                              ) * mult_hazard(i,0) * X(i,j)
                              *
                              (
                               Dmat(tleft[i],dmat_pos[xcols+k])*Survivor(i,0) - Dmat(tright[i]-1,dmat_pos[xcols+k])*Survivor(i,1)
                              ) * mult_hazard(i,0)
                             )
                             /
                             (
                              (Survivor(i,0)-Survivor(i,1)) * (Survivor(i,0)-Survivor(i,1))
                             );
              }
            else
              {
              H(j,xcols+k) += helpmat(tright[i]-1,0)*X(i,j)*mult_hazard(i,0);
              }
            }
          }

        else if(isbaselinebeta[j]==1 && isbaselinebeta[xcols+k]==0)
          {
          former=0;
          for(i=0; i<t_X.rows()-1; i++)
            {
            helpmat(i,0) = former - 0.5*tsteps(i,0)*(t_X(i,dm_pos[j])*baseline(i,0)
                                +t_X(i+1,dm_pos[j])*baseline(i+1,0));
            former=helpmat(i,0);
            }
          helpmat(t_X.rows()-1,0) = former;
          for(i=0; i<nrobs; i++)
            {
            if(interval[i])
              {
              H(j,xcols+k) += (
                               (
                                helpmat(tleft[i],0) + Dmat(tleft[i],dmat_pos[j])*negcumbaseline(tleft[i],0)*mult_hazard(i,0)
                               ) * Survivor(i,0)
                               -
                               (
                                helpmat(tright[i]-1,0) + Dmat(tright[i]-1,dmat_pos[j])*negcumbaseline(tright[i]-1,0)*mult_hazard(i,0)
                               ) * Survivor(i,1)
                              ) * mult_hazard(i,0) * Z(i,k) / (Survivor(i,0)-Survivor(i,1))
                              -
                              (
                               (
                                Dmat(tleft[i],dmat_pos[j])*Survivor(i,0)-Dmat(tright[i]-1,dmat_pos[j])*Survivor(i,1)
                               ) * mult_hazard(i,0)
                               *
                              (
                               negcumbaseline(tleft[i],0)*Survivor(i,0)-negcumbaseline(tright[i]-1,0)*Survivor(i,1)
                              ) * mult_hazard(i,0) * Z(i,k)
                             )
                             /
                             (
                              (Survivor(i,0)-Survivor(i,1))*(Survivor(i,0)-Survivor(i,1))
                             );
              }
            else
              {
              H(j,xcols+k) += helpmat(tright[i]-1,0)*Z(i,k)*mult_hazard(i,0);
              }
            }
          }

        else
          {
          for(i=0; i<nrobs; i++)
            {
            if(interval[i])
              {
              H(j,xcols+k) += (
                               (
                                negcumbaseline(tleft[i],0) + negcumbaseline(tleft[i],0)*negcumbaseline(tleft[i],0)*mult_hazard(i,0)
                               ) * Survivor(i,0)
                               -
                               (
                                negcumbaseline(tright[i]-1,0) + negcumbaseline(tright[i]-1,0)*negcumbaseline(tright[i]-1,0)*mult_hazard(i,0)
                               ) * Survivor(i,1)
                              ) * mult_hazard(i,0) * X(i,j) * Z(i,k) / (Survivor(i,0)-Survivor(i,1))
                              -
                              (
                               (
                                negcumbaseline(tleft[i],0)*Survivor(i,0)-negcumbaseline(tright[i]-1,0)*Survivor(i,1)
                               ) * mult_hazard(i,0)
                               *
                              (
                               negcumbaseline(tleft[i],0)*Survivor(i,0)-negcumbaseline(tright[i]-1,0)*Survivor(i,1)
                              ) * mult_hazard(i,0) * X(i,j) * Z(i,k)
                             )
                             /
                             (
                              (Survivor(i,0)-Survivor(i,1))*(Survivor(i,0)-Survivor(i,1))
                             );
              }
            else
              {
              H(j,xcols+k) += negcumbaseline(tright[i]-1,0)*X(i,j)*Z(i,k)*mult_hazard(i,0);
              }
            }
          }
        H(xcols+k,j)=H(j,xcols+k);
        }
      }

    // Z & Z
    for(j=0; j<zcols; j++)
      {
      for(k=j; k<zcols; k++)
        {
        // helpmat contains the second derivatives of the cumulated baseline
        // with respect to beta_j and beta_k
        helpmat = datamatrix(t_X.rows(),1,0);

        if(isbaselinebeta[xcols+j]==1 && isbaselinebeta[xcols+k]==1)
          {
          former=0;
          for(i=0; i<t_X.rows()-1; i++)
            {
            helpmat(i,0) = former - 0.5*tsteps(i,0)*(t_Z(i,dm_pos[xcols+j])*t_Z(i,dm_pos[xcols+k])*baseline(i,0)
                                +t_Z(i+1,dm_pos[xcols+j])*t_Z(i+1,dm_pos[xcols+k])*baseline(i+1,0));
            former = helpmat(i,0);
            }
          helpmat(t_X.rows()-1,0) = former;
          for(i=0; i<nrobs; i++)
            {
            if(interval[i])
              {
              H(xcols+j,xcols+k) += (
                                     (
                                      helpmat(tleft[i],0) + Dmat(tleft[i],dmat_pos[xcols+j])*Dmat(tleft[i],dmat_pos[xcols+k])*mult_hazard(i,0)
                                     ) * Survivor(i,0)
                                     -
                                     (
                                      helpmat(tright[i]-1,0) + Dmat(tright[i]-1,dmat_pos[xcols+j])*Dmat(tright[i]-1,dmat_pos[xcols+k])*mult_hazard(i,0)
                                     ) * Survivor(i,1)
                                    ) * mult_hazard(i,0) / (Survivor(i,0)-Survivor(i,1))
                                    -
                                    (
                                     (
                                      Dmat(tleft[i],dmat_pos[xcols+j])*Survivor(i,0) - Dmat(tright[i]-1,dmat_pos[xcols+j])*Survivor(i,1)
                                     ) * mult_hazard(i,0)
                                     *
                                     (
                                      Dmat(tleft[i],dmat_pos[xcols+k])*Survivor(i,0) - Dmat(tright[i]-1,dmat_pos[xcols+k])*Survivor(i,1)
                                     ) * mult_hazard(i,0)
                                    )
                                    /
                                    (
                                     (Survivor(i,0)-Survivor(i,1))*(Survivor(i,0)-Survivor(i,1))
                                    );
              }
            else
              {
              H(xcols+j,xcols+k) += helpmat(tright[i]-1,0)*mult_hazard(i,0);
              }
            }
          }

        else if(isbaselinebeta[xcols+j]==0 && isbaselinebeta[xcols+k]==1)
          {
          former=0;
          for(i=0; i<t_X.rows()-1; i++)
            {
            helpmat(i,0) = former - 0.5*tsteps(i,0)*(t_Z(i,dm_pos[xcols+k])*baseline(i,0)
                                +t_Z(i+1,dm_pos[xcols+k])*baseline(i+1,0));
            former = helpmat(i,0);
            }
          helpmat(t_X.rows()-1,0) = former;
          for(i=0; i<nrobs; i++)
            {
            if(interval[i])
              {
              H(xcols+j,xcols+k) += (
                                     (
                                      helpmat(tleft[i],0) + negcumbaseline(tleft[i],0)*Dmat(tleft[i],dmat_pos[xcols+k])*mult_hazard(i,0)
                                     ) * Survivor(i,0)
                                     -
                                     (
                                      helpmat(tright[i]-1,0) + negcumbaseline(tright[i]-1,0)*Dmat(tright[i]-1,dmat_pos[xcols+k])*mult_hazard(i,0)
                                     ) * Survivor(i,1)
                                    ) * mult_hazard(i,0) * Z(i,j) / (Survivor(i,0)-Survivor(i,1))
                                    -
                                    (
                                     (
                                      negcumbaseline(tleft[i],0)*Survivor(i,0) - negcumbaseline(tright[i]-1,0)*Survivor(i,1)
                                     ) * mult_hazard(i,0) * Z(i,j)
                                     *
                                     (
                                      Dmat(tleft[i],dmat_pos[xcols+k])*Survivor(i,0) - Dmat(tright[i]-1,dmat_pos[xcols+k])*Survivor(i,1)
                                     ) * mult_hazard(i,0)
                                    )
                                    /
                                    (
                                     (Survivor(i,0)-Survivor(i,1)) * (Survivor(i,0)-Survivor(i,1))
                                    );
              }
            else
              {
              H(xcols+j,xcols+k) += helpmat(tright[i]-1,0)*Z(i,j)*mult_hazard(i,0);
              }
            }
          }

        else if(isbaselinebeta[xcols+j]==1 && isbaselinebeta[xcols+k]==0)
          {
          former=0;
          for(i=0; i<t_X.rows()-1; i++)
            {
            helpmat(i,0) = former - 0.5*tsteps(i,0)*(t_Z(i,dm_pos[xcols+j])*baseline(i,0)
                                +t_Z(i+1,dm_pos[xcols+j])*baseline(i+1,0));
            former = helpmat(i,0);
            }
          helpmat(t_X.rows()-1,0) = former;
          for(i=0; i<nrobs; i++)
            {
            if(interval[i])
              {
              H(xcols+j,xcols+k) += (
                                     (
                                      helpmat(tleft[i],0) + Dmat(tleft[i],dmat_pos[xcols+j])*negcumbaseline(tleft[i],0)*mult_hazard(i,0)
                                     ) * Survivor(i,0)
                                     -
                                     (
                                      helpmat(tright[i]-1,0) + Dmat(tright[i]-1,dmat_pos[xcols+j])*negcumbaseline(tright[i]-1,0)*mult_hazard(i,0)
                                     ) * Survivor(i,1)
                                    ) * mult_hazard(i,0) * Z(i,k) / (Survivor(i,0)-Survivor(i,1))
                                    -
                                    (
                                     (
                                      Dmat(tleft[i],dmat_pos[xcols+j])*Survivor(i,0) - Dmat(tright[i]-1,dmat_pos[xcols+j])*Survivor(i,1)
                                      ) * mult_hazard(i,0)
                                     *
                                     (
                                      negcumbaseline(tleft[i],0)*Survivor(i,0)-negcumbaseline(tright[i]-1,0)*Survivor(i,1)
                                     ) * mult_hazard(i,0) * Z(i,k)
                                    )
                                    /
                                    (
                                     (Survivor(i,0)-Survivor(i,1))*(Survivor(i,0)-Survivor(i,1))
                                    );
              }
            else
              {
              H(xcols+j,xcols+k) += helpmat(tright[i]-1,0)*Z(i,k)*mult_hazard(i,0);
              }
            }
          }

        else
          {
          for(i=0; i<nrobs; i++)
            {
            if(interval[i])
              {
              H(xcols+j,xcols+k) += (
                                     (
                                      negcumbaseline(tleft[i],0) + negcumbaseline(tleft[i],0)*negcumbaseline(tleft[i],0)*mult_hazard(i,0)
                                     ) * Survivor(i,0)
                                     -
                                     (
                                      negcumbaseline(tright[i]-1,0) + negcumbaseline(tright[i]-1,0)*negcumbaseline(tright[i]-1,0)*mult_hazard(i,0)
                                     ) * Survivor(i,1)
                                    ) * mult_hazard(i,0) * Z(i,j) * Z(i,k) / (Survivor(i,0)-Survivor(i,1))
                                    -
                                    (
                                     (
                                      negcumbaseline(tleft[i],0)*Survivor(i,0) - negcumbaseline(tright[i]-1,0)*Survivor(i,1)
                                     ) * mult_hazard(i,0)
                                     *
                                     (
                                      negcumbaseline(tleft[i],0)*Survivor(i,0) - negcumbaseline(tright[i]-1,0)*Survivor(i,1)
                                     ) * mult_hazard(i,0) * Z(i,j) * Z(i,k)
                                    )
                                    /
                                    (
                                     (Survivor(i,0)-Survivor(i,1))*(Survivor(i,0)-Survivor(i,1))
                                    );
              }
            else
              {
              H(xcols+j,xcols+k) += negcumbaseline(tright[i]-1,0)*Z(i,j)*Z(i,k)*mult_hazard(i,0);
              }
            }
          }
        H(xcols+k,xcols+j)=H(xcols+j,xcols+k);
        }
      }
    H = -H;

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

