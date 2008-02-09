// Modul zur Erzeugung von Zufallszahlen

#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined( RANDOM_INCLUDED )

#define RANDOM_INCLUDED


#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "tmatrix.h"
#include<vector>

namespace randnumbers
{

#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#define ABS(a) ((a) >= 0 ? (a) : -(a))
#define SQR(x)    ((x)*(x))
#define PI 3.141592654
#define PI2 9.869604401
#define sqrt_pi 2.506628275


// Erzeugen von auf (0,1) gleichverteilten Zufallszahlen

double __EXPORT_TYPE uniform(void);


inline double __EXPORT_TYPE phi(const double & x)
  {
  return 0.39894228*exp(-0.5*x*x);
  }


// Berechnen der Verteilungsfunktion der Standardnormalverteilung an der
// Stelle x

double __EXPORT_TYPE Phi2(const double & x);


double __EXPORT_TYPE invPhi (const double & p);

double __EXPORT_TYPE invPhi2 (const double & p);

double __EXPORT_TYPE Phi(const double & x);

// FUNCTION: ksdist
// TASK: Approximate Kolmogorow-Smirnov distribution function at value lambda
//       by a finite sum with 2*kmax +1 terms from -kmax to +kmax

double __EXPORT_TYPE ksdist(int kmax, double lambda);

static double t = 0.75;
static double Ft = ksdist(1000, t);
static double ts = SQR(PI)/(8*SQR(t));

// FUNCTION: kssample
// TASK: generates random numbers from a kolmogorov smirnov distribution

double __EXPORT_TYPE kssample(void);

// Erzeugen von standardnormalverteilten Zufallszahlen

double __EXPORT_TYPE rand_normal(void);

// Erzeugen von Zufallszahlen gemäß einer truncated normal distribution

double __EXPORT_TYPE trunc_normal(const double & a,const double & b,const double & mu,
                    const double & s = 1);

double __EXPORT_TYPE trunc_normal2(const double & a,const double & b,const double & mu,
                    const double & s = 1);

double __EXPORT_TYPE trunc_normal3(const double & a,const double & b,const double & mu,
                    const double & s = 1);

double __EXPORT_TYPE truncnormal(const double & a,const double & b);

// Erzeugen von exponentialverteilten Zufallszahlen mit Parameter lambda

inline double __EXPORT_TYPE rand_expo(double lambda)
  {
  return (- 1/lambda)*log(uniform());
  }


// Erzeugen einer gammaverteilten Zufallszahl mit Parametern a und b
// Für a > 1 Best's Rejection Algorithmus (vgl. Devroye (1986) S.410)
// Für a = 1 Exponentialverteilung
// Für a < 1 Stuart's Theorem (vgl. Devroye (1986) S.182)
// Dichte der Gammaverteilung:
// f(x) = b^a * Gamma(a)^-1 * x^a-1 * exp(-bx)
// E(X) = a/b     Var(X) = a/b^2

double __EXPORT_TYPE rand_gamma(double a,double b);


// erzeugen einer invers gammaverteilten Zufallszahl mit Parametern a,b
// E(X) = b/(a-1) für a > 1
// Var(X) = b^2/((a-1)^2 * (a-2))

inline double __EXPORT_TYPE rand_invgamma(double a,double b)
  {
  return 1/rand_gamma(a,b);
  }

// Erzeugung von Chi-Quadrat verteilten Zufallszahlen
// mit n Freiheitsgraden

inline double __EXPORT_TYPE rand_chisquare (unsigned n)
  {
  return rand_gamma ((double(n)/2.0),0.5);
  }

// Erzeugen eines standardnormalverteilten Zufallsvektors (Spaltenvektor)
// mit Dimension dim !

Matrix<double> __EXPORT_TYPE rand_normvek(unsigned dim);

// Erzeugung einer Wishart verteilten Zufallsmatrix mit n Freiheitsgraden
// und Skalenparameter Sigma und Dimension q x q

void __EXPORT_TYPE rand_wishart(Matrix<double> & Sigma,const unsigned & n,Matrix<double> & res);


//Erzeugen einer Inverse Gaussian Zufallszahl mit Parametern mu, lambda
//E(X)=mu
//Var(X)=mu^3/lambda
//Devroye 1986 (S. 148)

double __EXPORT_TYPE rand_inv_gaussian(const double mu, const double lambda);

// Erzeugen einer Zufallszahl x ~ 1+1/x auf dem Intervall [1/f,f]

double __EXPORT_TYPE rand_variance(const double f);


// distribution function of standard logistic distribution

double __EXPORT_TYPE invlogit(double x);

// inverse distribution function of standard logistic distribution

double __EXPORT_TYPE logit(double x);

// generates a sample from the logistic distribution
// with mean mean truncated to values only left from zero

double __EXPORT_TYPE trunc_logistic_left(double mean);

// left = 1: sampling left from zero
// left = 0: sampling right from zero

double __EXPORT_TYPE trunc_logistic(double mean, int left);

// Generator for inverse Gaussian distribution
// taken from Devroye, page 149

double __EXPORT_TYPE IG(double mu, double lambda);


double __EXPORT_TYPE GIG(double lambda, double psi, double chi);


double __EXPORT_TYPE f1old(double x, int j);


double __EXPORT_TYPE f2old(double x, int j);


double __EXPORT_TYPE lambda_fc(double chi);

// Erzeugen einer betaverteilten Zufallszahl mit Parametern a und b
// (vgl. Devroye (1986) S.430)
// Dichte der Betaverteilung:
// f(x) = Gamma(a+b)/Gamma(a)*Gamma(b) * x^a-1 * (1-x)^b-1
// E(X) = a/a+b     Var(X) = a*b/(b+b)^2*(a+b+1)

double __EXPORT_TYPE rand_beta(double a, double b);

// Erzeugen einer dirichletverteilten Zufallszahl mit Parameternvektor
// alpha=(a_1,....,a_nrpar)
// (vgl. Devroye (1986) S.594)
// Dichte der Dirichletverteilung:
// f(x) = Gamma(a_1+...+a_nrpar)/Gamma(a_1)*...*Gamma(a_nrpar) * x_1^a_1-1 * ... * x_nrpar^a_nrpar-1
// mit x_nrpar = 1-x_1-...-x_nrpar-1

vector<double> __EXPORT_TYPE rand_dirichlet(double nrpar, vector<double> alpha);


}



#endif




