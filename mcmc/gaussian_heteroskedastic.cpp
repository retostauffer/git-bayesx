
#include "first.h"

#include "gaussian_heteroskedastic.h"

namespace MCMC
{


void DISTRIBUTION_gaussianh::standardise(void)
  {

  double s = sqrt(response.var(0,weight));

  trmult = datamatrix(1,1,s);

  unsigned i;
  double * workresp = response.getV();
  double * worklin = (*linpred_current).getV();
  for (i=0;i<nrobs;i++,workresp+=2,worklin+=2)
   {
   *workresp = *workresp/trmult(0,0);
   *worklin = *worklin/trmult(0,0);
   }


  datamatrix tr(1,1,trmult(0,0)*trmult(0,0));
  Scalesave.set_transformmult(tr);

  }


DISTRIBUTION_gaussianh::DISTRIBUTION_gaussianh(const double & a,
                   const datamatrix & b, MCMCoptions * o, const datamatrix & r,
                   const datamatrix & w)
  : DISTRIBUTION(o,r,w)
  {

  nrcat = response.cols(); //neu

  family = "Gaussian with heteroscedastic errors";

  standardise();

  scale(0,0) = 1;
  scaleexisting = false;

  constant_iwlsweights=true;


  }


DISTRIBUTION_gaussianh::DISTRIBUTION_gaussianh(
const DISTRIBUTION_gaussianh & nd)
: DISTRIBUTION(DISTRIBUTION(nd))
    {
        nrcat = nd.nrcat;//neu
    }


const DISTRIBUTION_gaussianh & DISTRIBUTION_gaussianh::operator=(
const DISTRIBUTION_gaussianh & nd)
  {

  if (this==&nd)
    return *this;
  DISTRIBUTION::operator=(DISTRIBUTION(nd));

  nrcat = nd.nrcat;//neu

  return *this;
  }



double DISTRIBUTION_gaussianh::loglikelihood(double * response,
                      double * linpred,
                      double * weight,const int & i) const//für eine Beob.
  {
        double * worklin = linpred;
        double eta1 = (*worklin); //erster Prediktor, also der für mu
        worklin++;
        double eta2 = exp((*worklin)); //zweiter Prediktor, also der für die
                                        //Varianz/ist die Anwendung von exp auf
                                        //linearen Prediktor hier notwendig,
                                        //dies hängt vom Schätzverfahren ab
        double help = (*response) - eta1;
        return -0.5 * (*linpred) - 0.5*(help*help)/eta2;

        // DISTRIBUTION_gaussian abschauen
        // DISTRIBUTION_multinomial abschauen

  return 0;
  }



void DISTRIBUTION_gaussianh::compute_mu(const double * linpred,double * mu)
                                           const
  {

  // DISTRIBUTION_gaussian abschauen
  // DISTRIBUTION_multinomial abschauen
  // zweite Spalte: exp(linpred für varianz)


  }


void DISTRIBUTION_gaussianh::compute_mu_notransform(const double * linpred,
double * mu) const
  {

  // DISTRIBUTION_gaussian abschauen
  // DISTRIBUTION_multinomial abschauen
  // zweite Spalte: exp(linpred für varianz)

  }


void DISTRIBUTION_gaussianh::compute_deviance(const double * response,
                             const double * weight,const double * mu,
                             double * deviance,double * deviancesat,
                             const datamatrix & scale,const int & i) const
  {

  // DISTRIBUTION_gaussian abschauen
  // DISTRIBUTION_multinomial abschauen
  // zweite Spalte: exp(linpred für varianz)

  }


double DISTRIBUTION_gaussianh::compute_weight(double * linpred, double * weight,
                        const int & i, const unsigned & col) const
  {

  // vgl. distribution h datei

  }


double DISTRIBUTION_gaussianh::compute_IWLS(double * response,double * linpred,
                                           double * weight,
                      const int & i,double * weightiwls,double * tildey,
                      bool weightyes, const unsigned & col)
  {

  // vgl. distribution h datei

  }

void DISTRIBUTION_gaussianh::compute_IWLS_weight_tildey(double * response,
                              double * linpred,
                              double * weight,const int & i,
                              double * weightiwls,double * tildey,
                              const unsigned & col)
  {

  // vgl. distribution h datei

  }

double DISTRIBUTION_gaussianh::compute_gmu(double * linpred,
const unsigned & col) const
  {

  // vgl. distribution h datei

  }


void DISTRIBUTION_gaussianh::outoptions(void)
  {
  DISTRIBUTION::outoptions();

  optionsp->out("\n");

  }


void DISTRIBUTION_gaussianh::update(void)
  {

//  DISTRIBUTION::update();

  }


void DISTRIBUTION_gaussianh::update_predict(void)
  {
  DISTRIBUTION::update_predict();
  }


void DISTRIBUTION_gaussianh::outresults(void)
  {

  }

bool DISTRIBUTION_gaussianh::posteriormode(void)
  {

  return true;

  }


bool DISTRIBUTION_gaussianh::posteriormode_converged_fc(const datamatrix & beta,
                                  const datamatrix & beta_mode,
                                  const unsigned & itnr)
  {
  return true;
  }


void DISTRIBUTION_gaussianh::compute_iwls(void)
  {

  // vgl. distribution h datei
  tildey.assign(response);
  }



} // end: namespace MCMC

