
#include<tvariance2dim.h>


namespace MCMC
{

FULLCOND_tvariance2dim::FULLCOND_tvariance2dim(MCMCoptions * o,
                   FULLCOND_pspline_surf_gaussian * p,
                     unsigned & v,const ST::string & ti, const ST::string & fp,
                     const ST::string & pres,const bool & rw)
           : FULLCOND(o,datamatrix(1,1),ti,180,1,fp)
  {
  rowwise = rw;
  Kp = p;
  pathresults = pres;
  nu = v;

  if (!rowwise)
    {
    m = sqrt(p->get_nrpar());
    nrpar = (m-1)*(m-1)*2+2*(m-1);
    }
  else
    {
    m = sqrt(p->get_nrpar());
    nrpar = m*m;
    u = datamatrix (nrpar,1,0);
    }

  setbeta(nrpar,1,1);

  }


FULLCOND_tvariance2dim::FULLCOND_tvariance2dim(const FULLCOND_tvariance2dim & t)
  : FULLCOND(FULLCOND(t))
  {
  rowwise = t.rowwise;
  Kp = t.Kp;
  pathresults = t.pathresults;
  nu = t.nu;
  m = t.m;
  u = t.u;
  }


const FULLCOND_tvariance2dim &
FULLCOND_tvariance2dim::operator=(const FULLCOND_tvariance2dim & t)
  {
  if (this == &t)
    return *this;
  FULLCOND::operator=(FULLCOND(t));
  rowwise = t.rowwise;
  Kp = t.Kp;
  pathresults = t.pathresults;
  nu = t.nu;
  m = t.m;
  u = t.u;
  return *this;
  }



void FULLCOND_tvariance2dim::update(void)
  {

  Kmatdiag = Kp->getdiagpointer();
  Kmatupper = Kp->getupperpointer();

  if (!rowwise)
    {
    unsigned i,j;
    int k = 0;
    double aneu = double(nu)/2+0.5;
    double bneu;

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

    unsigned size = m*m;
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

    }


  acceptance++;

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
  optionsp->out("\n");

  }




} // end: namespace MCMC



