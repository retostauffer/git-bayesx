// Modul zur Erzeugung von Zufallszahlen

#include "random.h"

namespace randnumbers
{


double Phi2(const double & x)
  {

  double xt;
  if (x < 0)
    xt = -x;
  else
    xt = x;
  double t=xt*xt;
  double pol = 1+0.196854*xt+0.115194*t;
  t*=xt;
  pol+=0.000344*t;
  t*=xt;
  pol+=0.019527*t;

  pol = 1.0/pow(pol,4);

  if (x < 0)
    return 0.5*pol;
  else
    return 1.0-0.5*pol;
  }


double __EXPORT_TYPE invPhi (const double & p);


double invPhi2 (const double & p)
  {

    double pt;
    if (p < 0.5)
      pt = p;
    else
      pt = 1-p;

    if(pt < 1.0e-100)
      pt = 1.0e-100;

    double t = sqrt(-2.0*log(pt));
    double t2 = t*t;
    double t3 = t2*t;
    double r1 = 2.515517+0.802853*t+0.010328*t2;
    double r2 = 1+1.432788*t+0.189269*t2+0.001308*t3;

    if (p<0.5)
      return -t+r1/r2;
    else
      return t-r1/r2;

  }


double uniform(void)
  {


/*
  randomize();
  static double  m = 2147483648;
  static double a = 843314861;
  static double c = 453816693;


  static double x = int(rand()+5);

  x = fmod(a*x+c,m);
  return x/m;
*/

  int zufall = 0;
  while ((zufall == 0) || (zufall == RAND_MAX))
	 {
	 zufall = rand();
	 }
  return double(zufall)/double(RAND_MAX);

  }



double rand_gamma(double a,double b)
  {
  if (a > 1)
	 {
	 double h1 = a-1;         // h1 entspricht b in Devroye (1986)
	 double h2 = 3*a-0.75;    // h2 entspricht c in Devroye (1986)
	 double U,V,W,Y,X,Z;
	 int accept = 0;
	 do
		{
		U = uniform();
		V = uniform();
		W = U*(1-U);
		Y = sqrt(h2/W)*(U-0.5);
		X = h1 + Y;
		if (X > 0)
			{
			Z = 64*W*W*W*V*V;
			if ( Z <= (1 - (2*Y*Y)/X) )
			  accept = 1;
			else
			  {
			  if ( ((X/h1) > 0) &&  ( log(Z) <= ( 2*(h1*log(X/h1) - Y) ) ) )
				 accept = 1;
			  }
			}
		}
	 while (accept == 0);
	 return X/b;
	 }
  else
	 {
	 if (a == 1)
		return rand_expo(b);
	 else
		{
		double X = rand_gamma(a+1,1)*pow(uniform(),1/a);
		return X/b;
		}
	 }
  }


/*
double Phi(const double & x)
  {
  if (x==0)
    return 0.5;
  else
    {
    double a,b;
    if (x > 0)
      {
      a = 0;
      b = x;
      }
    else if (x<0)
      {
      a=x;
      b=0;
      }

    double h = (b-a)/50.0;
    double x2 = -0.5*a*a;
    double sum= exp(x2);
    double xhelp = a+h;
    unsigned i;
    for (i=1;i<=25;i++)
      {
      x2 = -0.5*xhelp*xhelp;
      sum+=4*exp(x2);
      xhelp+=h;

      x2 = -0.5*xhelp*xhelp;
      if (i == 25)
        sum+=exp(x2);
      else
        sum+=2*exp(x2);
      xhelp+=h;

      }

    double c = 0.13298076*h;

    sum *= c;

    if (x > 0)
      return sum+0.5;
    else
      return 0.5-sum;
    }

  }

*/

double ksdist(int kmax, double lambda)
  {
  double result = 0.0;
  int k;
  for(k=-kmax; k<=kmax; k++)
    if((k%2)==0)
      result += exp(-2*SQR(lambda)*SQR(k));
    else
      result -= exp(-2*SQR(lambda)*SQR(k));
  return(result);
}

double kssample(void)
  {
  double U = uniform();
  int accept = 0;
  double E0, E1, n, G, X, W, Z, P, Q, U2, E, U3, dummy;
  int j;
  if(U < Ft)  // generator for the leftmost interval
    {
    do{

      do{
        E0 = -log(uniform());        // exponential
        E1 = -log(uniform());       // exponential
        E0 = E0/(1-1/(2*ts));
        E1 = 2*E1;
        G = ts+E0;
        accept = (SQR(E0) <= (ts*E1*(G+ts)));

        if(!accept)
          accept = ((G/ts-1-log(G/ts)) <= (E1));
        }
      while (!accept);

      X = PI/sqrt(8*G);
      W = 0.0;
      Z = 1/(2*G);
      P = exp(-G);
      n = 1.0;
      Q = 1.0;
      U2 = uniform();

      do{
        W = W + Z*Q;
        if(U2 >= W)
          return(X);
        n = n+2;
        Q = P;
      for(j=2; j<=(SQR(n)-1); j++)
         Q *= P;
         W = W - SQR(n)*Q;
         }
      while(U2 >= W);

      }
    while(1);
    }
  else  // generator for the rightmost interval
    {

    do
      {
      E = -log(uniform()); // exponential
      U3 = uniform();
      X = sqrt(SQR(t) + E/2.0);
      W = 0.0;
      n = 1.0;
      Z = exp(-2*SQR(X));
      do
        {
        n++;
        dummy = Z;
        for(j=2; j<=(SQR(n)-1); j++)
          dummy *= Z;
        W = W + SQR(n)*dummy;
        if (U3 >= W)
          return(X);
          n++;
        dummy = Z;
        for(j=2; j<=(SQR(n)-1); j++)
          dummy *= Z;
        W = W - SQR(n)*dummy;

        }
      while(U3 > W);
      }
    while(1);
    }

  }


double trunc_normal(const double & a,const double & b,const double & mu,
                   const double & s)
  {
  bool accept = false;
  double rand;
  while (!accept)
    {
    rand = mu+s*rand_normal();
    if ((rand <= b) && (rand >= a))
      accept = true;
    }

  return rand;
  }

double trunc_normal2(const double & a,const double & b,const double & mu,
                    const double & s)
  {
  double at = Phi2((a-mu)/s);
  double bt = Phi2((b-mu)/s);
  double u = at+(bt-at)*uniform();
  double r = mu+s*invPhi2(u);
  if (r < a)
    r = a+0.00000001;
  if (r > b)
    r = b-0.00000001;

  return r;
  }

double trunc_normal3(const double & a,const double & b,const double & mu,
                    const double & s)
  {
  double z;
  double at = (a-mu)/s;
  double bt = (b-mu)/s;
  double u = 1;
  double r = 0;
  while(u > r)
    {
    z = (bt-at)*uniform()+at;

    if(0 < at)
      r = exp((at*at-z*z)/2);
    else if (0 > bt)
      r = exp((bt*bt-z*z)/2);
    else
      r = exp((-z*z)/2);

    u = uniform();
    }

  return mu+s*z;
  }

double truncnormal(const double & a,const double & b)
  {
  bool accept = false;
  double rand;
  if (a > 2.5)
    {
    double u;
    while (!accept)
      {
      u = uniform();
      rand = a+(b-a)*uniform();
      if ( u <= (Phi(rand)/Phi(a)) )
        accept=true;
      }
    }
  else if (b < -2.5)
    {
    double u;
    while (!accept)
      {
      u = uniform();
      rand = a+(b-a)*uniform();
      if ( u <= (Phi(rand)/Phi(b)) )
        accept=true;
      }
    }
  else
    {

    while (!accept)
      {
      rand = rand_normal();
      if ((rand <= b) && (rand >= a))
        accept = true;
      }

    }

  return rand;
  }



// Erzeugen eines standardnormalverteilten Zufallsvektors (Spaltenvektor)
// mit Dimension dim !

Matrix<double> rand_normvek(unsigned dim)
  {
  Matrix<double> stnorm(dim,1);
  for (unsigned i=0; i < dim; i++)
	 stnorm(i,0) = rand_normal();
  return stnorm;
  }

// Erzeugung einer Wishart verteilten Zufallsmatrix mit n Freiheitsgraden
// und Skalenparameter Sigma und Dimension q x q
// w gibt an, ob es sich bei Sigma schon um die Choleskyzerlegung einer
// Kovarianzmatrix handelt, oder ob Sigma noch zerlegt werden muß
// w = 1 (default) entspricht schon zerlegt ansonsten unzerlegt

void rand_wishart(Matrix<double> & Sigma,const unsigned & n,Matrix<double> & res)
  {
  unsigned p = Sigma.rows();
  Matrix<double> V(p,p);
  Matrix<double> E(p,p);

  unsigned i,j,k;

  double zeta;
  double sum;

  for(i=0;i<p;i++)
    for(j=0;j<p;j++)
      E(i,j) = rand_normal();

  for(i=0;i<p;i++)
    {

    zeta = rand_chisquare(n-i);

    sum=0;
    for(k=0;k<i;k++)
      sum += E(k,i)*E(k,i);

    V(i,i) = zeta + sum;

    for(j=i+1;j<p;j++)
      {
      sum =  0;
      for(k=0;k<i;k++)
        sum += E(k,i)*E(k,j);
      V(i,j) = E(i,j)*sqrt(zeta) + sum;
      V(j,i) = V(i,j);
      }

    }

  Matrix<double> R = Sigma.root();

  res = R*V*R.transposed();

  }


  //Inverse Gaussian random numbers: Devroye 1986
double rand_inv_gaussian(const double mu, const double lambda)
{
    double N;
    double Y;
    double X1;
    N = rand_normal();
    Y = N*N;
    X1 = mu + mu*mu*Y/(2*lambda)-mu*sqrt(4*mu*lambda*Y+mu*mu*Y*Y)/(2*lambda);
    if(uniform()<= mu/(mu+X1))
    {
        return X1;
    }
    else
    {
        return mu*mu/X1;
    }
}


double rand_variance(const double f)
  {
  double length = f - 1/f;
  if (f == 1.0)
    return 1.0;
  if (uniform() < length/(length+2*log(f)))
    return (1/f + length*uniform());
  else
    return pow(f, 2.0*uniform()-1.0);
  }


double invlogit(double x)
  {
  return 1/(1+exp(-x));
  }


double logit(double x)
  {
  return log(x/(1-x));
  }


double trunc_logistic_left(double mean)
  {

  double Fa, arg, result;

  Fa = invlogit(mean);
  arg = Fa + uniform()*(1.0-Fa);
  result = mean - logit(arg);
  return result;
  }


double trunc_logistic(double mean, int left)
  {
  assert(left == 1 || left == 0);

  if(left == 1)
    return trunc_logistic_left(mean);
  else 
    return -trunc_logistic_left(-mean);
  }


double IG(double mu, double lambda)
  {

  double N, Y, X1, mu2;

  mu2 = SQR(mu);
  N = rand_normal();
  Y = SQR(N);
  X1 = mu + mu2*Y/(2*lambda)-mu/(2*lambda)*
       sqrt(4*mu*lambda*Y+mu2*SQR(Y));
  if(uniform() <= mu/(mu+X1))
    return X1;
  else
    return mu2/X1;
  }


double GIG(double chi)
  {
  // returns a sample from a generalized inverse gaussian
  // distribution with parameters lambda=0.5, chi, psi=1
  // Devroye, page 478 parametrization
  // assume lambda > 0

  double r = sqrt(chi);
  return r/IG(1., r);

  }


double f1old(double x, int j)
  {
  // first series approximation, guranteed to be monotone for x > 1.333
  double a;

  a = SQR(j+1)*exp(-0.5*x*(SQR(j+1)-1));
  return a;
  }


double f2old(double x, int j)
  {

  // second series approximation, guranteed to be monotone for x < 1.333
  double a;

  if ((j%2)==1)
    // odd
   a = (x/PI2)*exp(-((SQR(j)-1.)*PI2)/(2*x));
  else
    // even
   a =  SQR(j+1)*exp(-((SQR(j+1)-1.)*PI2)/(2*x));

  return a;
  }


double lambda_fc(double chi)
  {

  double fn_val,u;
  int ok_out;

  ok_out = 0;


  int  j, t=0;
  double aa, upper, lower, factor;

   while (!ok_out){

     fn_val = GIG(chi);

     // u is the random variable used to test for acceptance prob, 4 is the
     // supremum under the generalised-inverse-gaussian density
     u = uniform()*4.0;

     // using Devroye (and transformation of RVs) we can find an alternating
     // series approximation. In particular we can find a series of monotone
     // "squeezing" functions that bound the true ratio
     // (true_function(x)/sampling_function(x)) from above and below ever
     // more closely. However as in Devroye, in order to do this we need to
     // split up the x space into two regions within which we find a guaranteed
     // monotone series.

     // so in the first region.....

     if (fn_val > 1.334){

       j=1;
       factor = 4.0;
       upper = 1.0;
       // apply squeezing
       while (1){
         // first adjust the lower bound using the alternating series
         // f1() - given below
         lower = upper - f1old(fn_val,j);
         if (u < factor*lower){
           // if the draw of the acceptance prob is below the lower bound
           // then you are definatly a draw from the density - ACCEPT
           ok_out = 1;
           break;
         }
         // now adjust the upper bound
         upper = lower + f1old(fn_val,j+1);
         if (u > factor*upper){
           // if the draw of the acceptance prob is above the upper bound
           // then you are definatly NOT a draw from the density - REJECT
           ok_out=0;
           break;
         }
         // else u lies somewhere inbetween (lower, upper) so we're not sure
         // and we must continue
         j+=2;
       }
     }
     else{ // you are in the other region (the above series f1() is not
           // guarenteed to be monotone hence we must find other monotone series
       j=1;

       // aa is simply the supremum under the rejection sampler which lies at
       // the boundary
       aa = 0.5*log(2*PI)+2*log(PI)+log(4.0)-2.5*log(fn_val)-
           (SQR(PI)/(2*fn_val))+0.5*fn_val;
       factor = exp(aa);
       upper = 1;
       while (1){
         // this bit is the same as above but we use the series f2()
         lower = upper - f2old(fn_val,j);
         if (u < factor*lower){
           ok_out = 1;
           break;
         }
         upper = lower + f2old(fn_val,j+1);
         if (u > factor*upper){
           ok_out=0;
           break;
         }
         j=j+2;
       }
     }

     // t counts the number of reject steps in the sampler
     t += 1;
   }
   return fn_val;
}



} // end: namespace randnumbers
