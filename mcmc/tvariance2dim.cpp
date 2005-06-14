
#include<tvariance2dim.h>


namespace MCMC
{

FULLCOND_tvariance2dim::FULLCOND_tvariance2dim(MCMCoptions * o,
                   FULLCOND_pspline_surf_gaussian * p,
                     unsigned & v,const ST::string & ti, const ST::string & fp,
                     const ST::string & pres,const unsigned & bs,const bool & rw)
           : FULLCOND(o,datamatrix(1,1),ti,1,1,fp)
  {
  spatial = false;

  rowwise = rw;
  Kp = p;
  pathresults = pres;
  nu = v;

  if (!rowwise)
    {
    m = sqrt(p->get_nrpar());
    nrpar = (m-1)*(m-1)*2+2*(m-1);

    SparseMatrix Ksp = Kmrflinear(m,m);
    datamatrix de(m*m-1,1,0.0);
    datamatrix ud = datamatrix(m*m-1,m,0.0);
    for(unsigned i=0;i<m*m-2;i++)
      {
      de(i,0) = Ksp(i,i);
      for(unsigned j=0;j<ud.cols();j++)
        {
        if (i+j+1 < m*m)
          ud(i,j) = Ksp(i,i+j+1);
        }
      }
    de(m*m-2,0) = Ksp(m*m-2,m*m-2);
    K11 = envmatdouble(bandmatdouble(de,ud));
    detalt = K11.getLogDet();
    detneu = 0.0;
    nrrows = unsigned(bs/2);
    betakvec = vector<double>(0);
    rowvec = vector<double>(0);
    colvec = vector<double>(0);
    deltapropvec = vector<double>(0);
    }
  else
    {
    m = sqrt(p->get_nrpar());
    nrpar = m*m;
    u = datamatrix (nrpar,1,0);
    }

  setbeta(nrpar,1,1);

  }


FULLCOND_tvariance2dim::FULLCOND_tvariance2dim(MCMCoptions * o,
                   FULLCOND_nonp_gaussian * p,
                     unsigned & v,const ST::string & ti, const ST::string & fp,
                     const ST::string & pres,const unsigned & bs,const bool & rw)
           : FULLCOND(o,datamatrix(1,1),ti,1,1,fp)
  {
  spatial = true;

  rowwise = rw;
  Kp_spat = p;
  pathresults = pres;
  nu = v;

    envmatdouble Kenv = Kp_spat->get_K();

    nrpar = 0;
    for(unsigned i=0;i<Kenv.getDim();i++)
      nrpar += Kenv.getDiag(i);
    nrpar = nrpar/2;

    datamatrix help(Kenv.getDim()-1,Kenv.getDim()-1,0);

    for(unsigned i=0;i<help.rows();i++)
      {
      for(unsigned j=0;j<help.cols();j++)
        {
        help(i,j) = Kenv(i,j);
        }
      }

    unsigned k = 0;
    indexmat = statmatrix<int>(nrpar,3,0);

    for(unsigned i=0;i<help.rows()+1;i++)
      {
      for(unsigned j=0;j<help.cols()+1;j++)
        {
        if(Kenv(i,j) != 0 && j>i)
          {
          indexmat(k,0) = i;
          indexmat(k,1) = j;
          k++;
          }
        }
      }

    K11 = envmatdouble(help);
    detalt = K11.getLogDet();
    detneu = 0.0;
    nrrows = bs;
    betakvec = vector<double>(0);
    rowvec = vector<double>(0);
    colvec = vector<double>(0);
    deltapropvec = vector<double>(0);

  setbeta(nrpar,1,1);

  }


FULLCOND_tvariance2dim::FULLCOND_tvariance2dim(const FULLCOND_tvariance2dim & t)
  : FULLCOND(FULLCOND(t))
  {
  rowwise = t.rowwise;
  Kp = t.Kp;
  Kp_spat = t.Kp_spat;
  spatial = t.spatial;
  indexmat = t.indexmat;
  pathresults = t.pathresults;
  nu = t.nu;
  m = t.m;
  u = t.u;

  K11 = t.K11;
  detalt = t.detalt;
  detneu = t.detneu;
  betakvec = t.betakvec;
  rowvec = t.rowvec;
  colvec = t.colvec;
  deltapropvec = t.deltapropvec;
  nrrows = t.nrrows;
  }


const FULLCOND_tvariance2dim &
FULLCOND_tvariance2dim::operator=(const FULLCOND_tvariance2dim & t)
  {
  if (this == &t)
    return *this;
  FULLCOND::operator=(FULLCOND(t));
  rowwise = t.rowwise;
  Kp = t.Kp;
  Kp_spat = t.Kp_spat;
  spatial = t.spatial;
  indexmat = t.indexmat;  
  pathresults = t.pathresults;
  nu = t.nu;
  m = t.m;
  u = t.u;

  K11 = t.K11;
  detalt = t.detalt;
  detneu = t.detneu;
  betakvec = t.betakvec;
  rowvec = t.rowvec;
  colvec = t.colvec;
  deltapropvec = t.deltapropvec;
  nrrows = t.nrrows;
  return *this;
  }


void FULLCOND_tvariance2dim::update(void)
  {
  if(spatial)
    update_spat();
  else
    update_2dim();
  }

void FULLCOND_tvariance2dim::update_2dim(void)
  {

  if (!rowwise)
    {

    unsigned i,j,l;
    int k = 0;
    double aneu = double(nu)/2;
    double bneu;

    double alpha,u,betak;
    double deltaprop;
    unsigned row,col;

    unsigned dim = m*m;

    for (row=0;row<dim;row++)
      {

      i = row/m;
      j = row%m;

      col=row+1;
      if(j < m-1)
        {
        betak = beta(k,0);

        bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i,j+1,m);
        deltaprop = randnumbers::rand_gamma(aneu,bneu);

        deltapropvec.push_back(deltaprop);
        rowvec.push_back(row);
        colvec.push_back(col);
        betakvec.push_back(betak);

        K11.setDiag(row,K11(row,row) + deltaprop - beta(k,0));
        if(col<K11.getDim())
          {
          K11.set(row,col,-deltaprop);
          K11.setDiag(col,K11(col,col) + deltaprop - beta(k,0));
          }

        k++;
        }

      col=row+m;
      if(i < m-1)
        {
        betak = beta(k,0);

        bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i+1,j,m);
        deltaprop = randnumbers::rand_gamma(aneu,bneu);

        deltapropvec.push_back(deltaprop);
        rowvec.push_back(row);
        colvec.push_back(col);
        betakvec.push_back(betak);

        K11.setDiag(row,K11(row,row) + deltaprop - beta(k,0));
        if(col<K11.getDim())
          {
          K11.set(row,col,-deltaprop);
          K11.setDiag(col,K11(col,col) + deltaprop - beta(k,0));
          }

        k++;
        }

      if((row+1)%nrrows == 0 || (row+1)==dim)
        {

        if(detalt==detneu)
          K11.decomp2(row-nrrows+1);
        detneu = K11.getLogDet();

        alpha = 0.5*(detneu - detalt);
        u = log(uniform());

        nrtrials++;

        if(u <= alpha)
          {
          for(l=0;l<deltapropvec.size();l++)
	        {
    	    beta(k-deltapropvec.size()+l,0) = deltapropvec[l];
            Kp->setK(rowvec[l],colvec[l],-deltapropvec[l]);
	        }
          detalt = detneu;
          acceptance++;
          }
        else
          {
          for(l=0;l<deltapropvec.size();l++)
	        {
            K11.setDiag(rowvec[l],K11(rowvec[l],rowvec[l]) - deltapropvec[l] + betakvec[l]);
            if(colvec[l]<K11.getDim())
              {
              K11.set(rowvec[l],colvec[l],-betakvec[l]);
              K11.setDiag(colvec[l],K11(colvec[l],colvec[l]) - deltapropvec[l] + betakvec[l]);
              }
    	    }
          }

        deltapropvec = vector<double>(0);
	    rowvec = vector<double>(0);
        colvec = vector<double>(0);
        betakvec = vector<double>(0);

        } // END:       if(row%nrrows == 0)
      } // END:    for (row=0;row<dim;row++)

/*
    unsigned i,j;
    int k = 0;
//    double aneu = double(nu)/2;
    double aneu = double(nu)/2+0.5;
    double bneu;

    double deltaprop,deltaprop2,alpha,u;
    unsigned row,col;

    for (i=0;i<m;i++)
      {
      for(j=0;j<m;j++)
        {

        if ( (i < m-1) && (j < m-1) )
          {

          row = i*m+j;
          col = i*m+j+1;

          bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i,j+1,m);
          deltaprop = randnumbers::rand_gamma(aneu,bneu);

          K11.set(row,col,-deltaprop);
          K11.setDiag(row,K11(row,row) + deltaprop - beta(k,0));
          K11.setDiag(col,K11(col,col) + deltaprop - beta(k,0));

          k++;

          row = i*m+j;
          col = (i+1)*m+j;

          bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i+1,j,m);
          deltaprop2 = randnumbers::rand_gamma(aneu,bneu);

          K11.set(row,col,-deltaprop2);
          K11.setDiag(row,K11(row,row) + deltaprop2 - beta(k,0));
          K11.setDiag(col,K11(col,col) + deltaprop2 - beta(k,0));

          if(detalt==detneu)
            K11.decomp2(row);
          detneu = K11.getLogDet();

          alpha = 0.5*(detneu - detalt);
          alpha += 0.5*( log(beta(k-1,0)) + log(beta(k,0)) - log(deltaprop) - log(deltaprop2) );
          u = log(uniform());
          nrtrials++;

          if(u <= alpha)
            {
            beta(k-1,0) = deltaprop;
            beta(k,0) = deltaprop2;

            Kp->setK(row,col-m+1,-deltaprop);
            Kp->setK(row,col,-deltaprop2);

            detalt = detneu;
            acceptance++;
            }
          else
            {
            col = i*m+j+1;
            K11.set(row,col,-beta(k-1,0));
            K11.setDiag(row,K11(row,row) - deltaprop + beta(k-1,0));
            K11.setDiag(col,K11(col,col) - deltaprop + beta(k-1,0));

            col = (i+1)*m+j;
            K11.set(row,col,-beta(k,0));
            K11.setDiag(row,K11(row,row) - deltaprop2 + beta(k,0));
            K11.setDiag(col,K11(col,col) - deltaprop2 + beta(k,0));
            }

          k++;

          }  // end: if ( (i < m-1) && (j < m-1) )

        if ( (i< m-1) && (j == m-1) )
          {

          row = i*m+j;
          col = (i+1)*m+j;

          bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i+1,j,m);
          deltaprop = randnumbers::rand_gamma(aneu,bneu);

          K11.setDiag(row,K11(row,row) + deltaprop - beta(k,0));
          if(col<K11.getDim())
            {
            K11.set(row,col,-deltaprop);
            K11.setDiag(col,K11(col,col) + deltaprop - beta(k,0));
            }

          if(detalt==detneu)
            K11.decomp2(row);
          detneu = K11.getLogDet();

          alpha = 0.5*(detneu - detalt);
          alpha += 0.5*( log(beta(k,0)) - log(deltaprop) );
          u = log(uniform());
          nrtrials++;

          if(u <= alpha)
            {
            beta(k,0) = deltaprop;

            Kp->setK(row,col,-deltaprop);

            detalt = detneu;
            acceptance++;
            }
          else
            {
            K11.setDiag(row,K11(row,row) - deltaprop + beta(k,0));
            if(col<K11.getDim())
              {
              K11.set(row,col,-beta(k,0));
              K11.setDiag(col,K11(col,col) - deltaprop + beta(k,0));
              }
            }

          k++;

          } // end: if ( (i< m-1) && (j == m-1) )

        if ( (i == m-1) && (j < m-1) )
          {

          row = i*m+j;
          col = i*m+j+1;

          bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i,j+1,m);
          deltaprop = randnumbers::rand_gamma(aneu,bneu);

          K11.setDiag(row,K11(row,row) + deltaprop - beta(k,0));
          if(col<K11.getDim())
            {
            K11.set(row,col,-deltaprop);
            K11.setDiag(col,K11(col,col) + deltaprop - beta(k,0));
            }

          if(detalt==detneu)
            K11.decomp2(row);
          detneu = K11.getLogDet();

          alpha = 0.5*(detneu - detalt);
          alpha += 0.5*( log(beta(k,0)) - log(deltaprop) );
          u = log(uniform());
          nrtrials++;

          if(u <= alpha)
            {
            beta(k,0) = deltaprop;

            Kp->setK(row,col,-deltaprop);

            detalt = detneu;
            acceptance++;
            }
          else
            {
            K11.setDiag(row,K11(row,row) - deltaprop + beta(k,0));
            if(col<K11.getDim())
              {
              K11.set(row,col,-beta(k,0));
              K11.setDiag(col,K11(col,col) - deltaprop + beta(k,0));
              }
            }

          k++;

          } // end: if ( (i == m-1) && (j < m-1) )

        }  // end: for(j=0;j<m;j++)

      } // end: for (i=0;i<m;i++)
*/
/*
    for (i=0;i<m;i++)
      {
      for(j=0;j<m;j++)
        {

        if ( (i < m-1) && (j < m-1) )
          {

          bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i,j+1,m);
          beta(k,0) = randnumbers::rand_gamma(aneu,bneu);
          *Kmatupper = -beta(k,0);
          k++;

          bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i+1,j,m);
          beta(k,0) = randnumbers::rand_gamma(aneu,bneu);
          Kmatupper+= m-1;
          *Kmatupper = -beta(k,0);
          *Kmatupper++;
          k++;

          }  // end: if ( (i < m-1) && (j < m-1) )

        if ( (i< m-1) && (j == m-1) )
          {

          bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i+1,j,m);
          beta(k,0) = randnumbers::rand_gamma(aneu,bneu);
          Kmatupper+= m-1;
          *Kmatupper = -beta(k,0);
          Kmatupper++;
          k++;

          } // end: if ( (i< m-1) && (j == m-1) )

        if ( (i == m-1) && (j < m-1) )
          {
          bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i,j+1,m);
          beta(k,0) = randnumbers::rand_gamma(aneu,bneu);
          *Kmatupper = -beta(k,0);
          Kmatupper+= m;
          k++;

          } // end: if ( (i == m-1) && (j < m-1) )

        }  // end: for(j=0;j<m;j++)

      } // end: for (i=0;i<m;i++)
*/

    unsigned size = m*m;
    Kmatdiag = Kp->getdiagpointer();
    Kmatupper = Kp->getupperpointer();
    double * phelp;

    for (i=0;i<size;i++,Kmatdiag++)
      {

      *Kmatdiag = 0;

      if (i < size-1)
        {
        *Kmatdiag += - *Kmatupper;
        Kmatupper += m-1;
        *Kmatdiag += - *Kmatupper;
        Kmatupper ++;
        }

      if (i > 0)
        {

        phelp =  Kp->getupperpointer()+(i-1)*m;

        *Kmatdiag += - *phelp;

        if (i >= m)
          {

          phelp =  Kp->getupperpointer()+(i-m)*m + m-1;
          *Kmatdiag += - *phelp;

          }

        } // end: if (i > 0)

      }

    }         // end: !rowwise
  else        // rowwise
    {

    unsigned i,j;

    Kp->compute_squareddiff(u);

    Kmatdiag = Kp->getdiagpointer();
    Kmatupper = Kp->getupperpointer();
    double * workbeta = beta.getV()+1;
    double * worku = u.getV()+1;
    double wold,wnew;
    double v = nu/2.0;
    double b,t;



    for (i=0;i<m;i++,workbeta++,worku++,Kmatdiag++,Kmatupper++)
      {
      *workbeta = rand_invgamma(v+0.5,v+0.5* *worku);
      wold=1.0/(*workbeta);
      *Kmatdiag = wold;
      *Kmatupper = -wold;
      Kmatdiag++;
      Kmatupper++;
      workbeta++;
      worku++;

      for (j=1;j<m-1;j++,Kmatdiag++,Kmatupper++,workbeta++,worku++)
        {
        *workbeta = rand_invgamma(v+0.5,v+0.5* *worku);
        wnew = 1.0/(*workbeta);
        *Kmatdiag = wold+wnew;
        *Kmatupper = -wnew;
        wold = wnew;
        }

      *Kmatdiag = wold;
      *Kmatupper = 0;


      }

    acceptance++;

    }

  FULLCOND::update();
  }


void FULLCOND_tvariance2dim::update_spat(void)
  {

    unsigned i,j,l;
    int k = 0;
    double aneu = double(nu)/2;
    double bneu;

    double alpha,u,betak;
    double deltaprop;
    unsigned row,col;

    while(k<nrpar)
      {

      betak = beta(k,0);
      row = indexmat(k,0);
      col = indexmat(k,1);

      bneu = nu*0.5 + 0.5*Kp_spat->compute_squareddiff(row,col);
      deltaprop = randnumbers::rand_gamma(aneu,bneu);

      deltapropvec.push_back(deltaprop);
      rowvec.push_back(row);
      colvec.push_back(col);
      betakvec.push_back(betak);

      K11.setDiag(row,K11(row,row) + deltaprop - beta(k,0));
      if(col<K11.getDim())
        {
        K11.set(row,col,-deltaprop);
        K11.setDiag(col,K11(col,col) + deltaprop - beta(k,0));
        }

      k++;

      if((row+1)%nrrows == 0 || k==nrpar)
        {

        if(detalt==detneu)
          K11.decomp2(row-nrrows+1);
        detneu = K11.getLogDet();

        alpha = 0.5*(detneu - detalt);
        u = log(uniform());

        nrtrials++;

        if(u <= alpha)
          {
          for(l=0;l<deltapropvec.size();l++)
	        {
    	    beta(k-deltapropvec.size()+l,0) = deltapropvec[l];
            Kp_spat->setK(rowvec[l],colvec[l],-deltapropvec[l]);

            Kp_spat->setK(colvec[l],colvec[l],Kp_spat->get(colvec[l],colvec[l])+deltapropvec[l]-betakvec[l]);
            Kp_spat->setK(rowvec[l],rowvec[l],Kp_spat->get(rowvec[l],rowvec[l])+deltapropvec[l]-betakvec[l]);

	        }
          detalt = detneu;
          acceptance++;
          }
        else
          {
          for(l=0;l<deltapropvec.size();l++)
	        {
            K11.setDiag(rowvec[l],K11(rowvec[l],rowvec[l]) - deltapropvec[l] + betakvec[l]);
            if(colvec[l]<K11.getDim())
              {
              K11.set(rowvec[l],colvec[l],-betakvec[l]);
              K11.setDiag(colvec[l],K11(colvec[l],colvec[l]) - deltapropvec[l] + betakvec[l]);
              }
    	    }
          }

        deltapropvec = vector<double>(0);
	    rowvec = vector<double>(0);
        colvec = vector<double>(0);
        betakvec = vector<double>(0);

        } // END:       if(row%nrrows == 0)
      } // END:    for (row=0;row<dim;row++)
/*
ofstream out("c:\\bayesx\\K11.raw");
K11.print4(out);
out.close();

ofstream out2("c:\\bayesx\\Kenv.raw");
(Kp_spat->get_K()).print4(out2);
out2.close();
*/

/*/ Diagonalelemente ausrechnen

  unsigned size = (Kp_spat->get_K()).getDim();

  double help;
  k=0;
  for(unsigned i=0;i<size;i++)
    {
    help = 0.0;
    row = indexmat(i,0);
    while(indexmat(i,0) == row)
      {
      help += beta(k,0);
      i++;
      }
    Kp_spat->setK(row,row,-help);
    k++;
    }
*/
//  acceptance++;

  FULLCOND::update();
  }


void FULLCOND_tvariance2dim::outresults(void)
  {

  FULLCOND::outresults();

  optionsp->out("  Results are stored in file " + pathresults + "\n");
  optionsp->out("\n");

  unsigned i;

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  ofstream outres(pathresults.strtochar());

  ST::string name = title;

  outres << "intnr" << "   ";
  outres << name << "   ";
  outres << "pmean   ";
  outres << "pqu"  << l1  << "   ";
  outres << "pqu"  << l2  << "   ";
  outres << "pmed   ";
  outres << "pqu"  << u1  << "   ";
  outres << "pqu"  << u2  << "   ";
  outres << endl;

  double * workmean = betamean.getV();
  double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
  double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
  double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
  double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
  double * workbetaqu50 = betaqu50.getV();

  for(i=0;i<nrpar;i++,workmean++,workbetaqu_l1_lower_p++,
                           workbetaqu_l2_lower_p++,workbetaqu50++,
                           workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++)
    {
    outres << (i+1) << "   ";
    outres << (i+1) << "   ";
    outres << *workmean << "   ";
    outres << *workbetaqu_l1_lower_p << "   ";
    outres << *workbetaqu_l2_lower_p << "   ";
    outres << *workbetaqu50 << "   ";
    outres << *workbetaqu_l2_upper_p << "   ";
    outres << *workbetaqu_l1_upper_p << "   ";
    outres << endl;
    }

  }                                         


void FULLCOND_tvariance2dim::outoptions(void)
  {
  optionsp->out("  OPTIONS FOR NONPARAMETRIC TERM: " + title +
                " (weights)\n",true);
  optionsp->out("\n");

  optionsp->out("  Hyperprior nu for variance parameter: " +
                ST::inttostring(nu) + "\n" );
  optionsp->out("  Blocksize for updating variances: ");
  if(spatial)
    optionsp->out(ST::inttostring(nrrows) + " rows of penalty matrix\n" );
  else
    optionsp->out(ST::inttostring(nrrows*2) + "\n" );  
  optionsp->out("\n");

  }




} // end: namespace MCMC



