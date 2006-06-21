
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
        double workmu = (*worklin); //erster Prediktor, also der für mu
        worklin++;
        double s = exp((*worklin)); //zweiter Prediktor, also der für die
                                        //Varianz/ist die Anwendung von exp auf
                                        //linearen Prediktor hier notwendig,
                                        //dies hängt vom Schätzverfahren ab
        double help = (*response) - workmu;
        return -0.5 * (*worklin) - 0.5*(help*help)/s;

        // DISTRIBUTION_gaussian abschauen
        // DISTRIBUTION_multinomial abschauen
        // Die Variablen weight ist hier überflüssig und wird nicht genutzt

  return 0;
  }



void DISTRIBUTION_gaussianh::compute_mu(const double * linpred,double * mu)
                                           const
  {
        const double * worklin = linpred;
        double * workmu = mu;

        //Zunächst den Wert für den Erwartungswertschätzer zuweisen
        (*mu) = trmult(0,0)* (*linpred);

        //Zeiger auf die Einträge des Varianzschätzers richten
        worklin++;
        workmu++;

        //Wert für den Varianzschätzer zuweisen, Transformation nicht notwendig,
        //da in Funktion standardise nicht vorgenommen
        (*workmu) = exp((*worklin));

  // DISTRIBUTION_gaussian abschauen
  // DISTRIBUTION_multinomial abschauen
  // zweite Spalte: exp(linpred für varianz)

  }


void DISTRIBUTION_gaussianh::compute_mu_notransform(const double * linpred,
double * mu) const
  {
        const double * worklin = linpred;
        double * workmu = mu;

        //Zunächst den Wert für den Erwartungswertschätzer zuweisen
        (*mu) = (*linpred);

        //Zeiger auf die Einträge des Varianzschätzers richten
        worklin++;
        workmu++;

        //Wert für den Varianzschätzer zuweisen
        (*workmu) = exp((*worklin));


  // DISTRIBUTION_gaussian abschauen
  // DISTRIBUTION_multinomial abschauen
  // zweite Spalte: exp(linpred für varianz)

  }


void DISTRIBUTION_gaussianh::compute_deviance(const double * response,
                             const double * weight,const double * mu,
                             double * deviance,double * deviancesat,
                             const datamatrix & scale,const int & i) const
  {
      const double  * workmu = mu;
      double r = (*response)* trmult(0,0) - (*mu);
      workmu++;
      double s = (*workmu)*pow(trmult(0,0),2);
      *deviance =  ((double)1.0/s)*r*r+log(2*M_PI*s);
      *deviancesat = ((double)1.0/s)*r*r;

  // DISTRIBUTION_gaussian abschauen
  // DISTRIBUTION_multinomial abschauen
  // zweite Spalte: exp(linpred für varianz)
  // die Variable weight als Argument der Funktion ist wahrscheinlich
  // überflüssig

  }


double DISTRIBUTION_gaussianh::compute_weight(double * linpred, double * weight,
                        const int & i, const unsigned & col) const
  {

  // vgl. distribution h datei
  // diese Funktion dürfte eigentlich überflüssig sein, da eine gewichtete
  // Regression, da die Heteroskedastizität der Varianz bereits durch
  // sigma_{i}^{2} modelliert wird.

  }


double DISTRIBUTION_gaussianh::compute_IWLS(double * response,double * linpred,
                                           double * weight,
                      const int & i,double * weightiwls,double * tildey,
                      bool weightyes, const unsigned & col)
  {

    double * workweightiwls = weightiwls;
    double * worktildey = tildey;
    double * worklinpred = linpred;
    double workmu = (*worklinpred);
    worklinpred++; //zeigt jetzt auf eta
    double s = exp((*worklinpred));
    double help = (*response)-workmu;


    if(col == 0) //Berechnung für den Prädiktor des Mittelwertes
    {
        (*workweightiwls) = (double)1.0/s;//1/exp(eta)

        (*worktildey) =  (*response); //Für den Erwartungswertschätzer stimmen
                                   //working-obs mit obs überein
    }
    if(col == 1) //Berechnung für den Prädiktor der Varianz
    {
        workweightiwls++;
        worktildey++;

        (*workweightiwls) = (double)0.5;

        (*worktildey) = (*worklinpred) + ((help*help)/s) - 1;
    }

    return - 0.5 * (*worklinpred) - 0.5 * (help*help)/s;


  // vgl. distribution h datei
  // Berechnet die IWLS-Gewichte und speichert diese in weightiwls. Als iwls
  // Gewichte werden die negativen Erwartungswerte der jeweiligen Gewichts-
  // matrizen verwendet. Dies entspricht Fisher-Scoring. Inwiefern dies im
  // Rahmen der Programmlogik zulässig ist, ist mir unklar.
  // Berechnet die working-observations ~y und speichert diese in tildey
  // Der i-te Summand der Loglikelihood wird zurückgegeben
  // Die Berechnung erfolgt für die i-te Beobachtung
  // Der Übergabewert weight wird nicht beachtet, da Heteroskedastizität mittels
  // sigma_{i}^{2} modelliert wird.
  // Der Übergabewert col kennzeichnet, ob Berechnungen für den Prädiktor
  // des Mittelwertes (col=0) oder für den der Varianz (col=1) durchgeführt
  // werden.
  // weightiwls gibt an, ob die IWLS-Gewichte berechnet werden sollen, d.h.
  // weightiwls = true, falls dies erfolgen sollen, weightiwls = false, sonst

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
  // wir benötigt, um die working observations zu berechnen

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

