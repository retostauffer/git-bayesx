
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (REALVAR_INCLUDED)

#define REALVAR_INCLUDED


#include<list>
#include"vectorn.h"
#include"realobs.h"

//------------------------------------------------------------------------------
//-------------------------- CLASS realvar -------------------------------------
//------------------------------------------------------------------------------

namespace realob
{

class __EXPORT_TYPE realvar : public vectorrealobs
  {

  public:


  //------------------------- PUBLIC FUNCTIONS ---------------------------------

  // DEFAULT CONSTRUCTOR

  realvar(void) : vectornum<realobs>()
	 {
	 }

  // CONSTRUCTOR

  realvar(int nr) : vectornum<realobs>(nr)
	 {
	 }

  // CONSTRUCTOR

  realvar(int nr,realobs o) : vectornum<realobs>(nr,o)
	 {
	 }

  // COPY CONSTRUCTOR

  realvar(const realvar & v) : vectornum<realobs>(vectornum<realobs>(v))
	 {
	 }

  realvar(const vectornum<realobs> & v) : vectornum<realobs>(v)
	 {
	 }

  realvar applied(realobs (*func)(realobs &)) const
    {
    realvar h(size());
    unsigned i;
    realobs help;
    for (i=0;i<size();i++)
      {
      help = operator[](i);
      h[i] = help;
      }
//    for (i=0;i<size();i++)
//	  h[i] = func(operator[](i));
    return h;
    }

  realobs min();

  realobs sum();


  };


//------------------------------------------------------------------------------
//---------------- FUNCTIONS THAT CAN BE APPLIED TO REALVAR --------------------
//------------------------------------------------------------------------------


//----------------------- MATHEMATICAL FUNCTIONS -------------------------------


//  | POSSIBLE FUNCTIONS                ABBREVIATION
//  | ----------------------------------------------
//  | square root                      | sqrt
//  | absolute value                   | abs
//  | expontential                     | exp
//	| cosinus                          | sin
//	| sinus                            | cos
//	| logarithm                        | log
//	| logarithm with basis 10          | log10
//  | floor function                   | floor
//  | power function                   | pow
//  | lag function                     | lag


realvar __EXPORT_TYPE sqrt(const realvar & v)
  {
  return v.applied(realob::sqrt);
  }


realvar __EXPORT_TYPE abs(const realvar & v)
  {
  return v.applied(realob::abs);
  }


realvar __EXPORT_TYPE exp(const realvar & v)
  {
  return v.applied(realob::exp);
  }


realvar __EXPORT_TYPE cos(const realvar & v)
  {
  return v.applied(realob::cos);
  }


realvar __EXPORT_TYPE sin(const realvar & v)
  {
  return v.applied(realob::sin);
  }


realvar __EXPORT_TYPE log(const realvar & v)
  {
  return v.applied(realob::log);
  }


realvar __EXPORT_TYPE log10(const realvar & v)
  {
  return v.applied(realob::log10);
  }


realvar __EXPORT_TYPE floor(const realvar & v)
  {
  return v.applied(realob::floor);
  }


// FUNCTION: lag
// TASK: computes the lagged vector of 'v' and returns the result

realvar __EXPORT_TYPE lag(realvar & v);

// FUNCTION: lag
// TASK: computes the lagged vector of 'v' and returns the result
// ADDITIONAL INFORMATION: assumes, that the 'true' order of the elements
//                         of vector 'v' is stored in 'index'

realvar __EXPORT_TYPE lagrealvar(realvar & v,vector<int> & index);

realvar __EXPORT_TYPE power(const realvar & v1,const realvar & v2);

//------------------------ STATISTICAL FUNCTIONS -------------------------------

realvar __EXPORT_TYPE cumulnorm(realvar & v);

// FUNCTION: cumul
// TASK: returns the empirical cumulative distribution function of 'v'

realvar __EXPORT_TYPE cumul(realvar & v,vector<int> & index);

// FUNCTION: uniform
// TASK: returns a realvar that contains 'obs' in (0,1) uniformly distributed
//       random numbers

realvar __EXPORT_TYPE uniform(unsigned obs);

// FUNCTION: normal
// TASK: returns a realvar that contains 'obs' standard normal distributed
//       random numbers

realvar __EXPORT_TYPE normal(unsigned obs);

// FUNCTION: exponential
// TASK: returns a realvar that contains 'obs' exponentially distributed
//       random numbers (with parameter lambda)

realvar __EXPORT_TYPE exponential(unsigned obs,realobs lambda);

realvar __EXPORT_TYPE exponential(realvar & lambda);

// FUNCTION: bernoulli
// TASK: returns a realvar that contains 'obs' bernoulli distributed random
//       numbers (with parameter p)

realvar __EXPORT_TYPE bernoulli(realvar & p);

realvar __EXPORT_TYPE binomial(realvar & n,realvar & p);

realvar __EXPORT_TYPE gamma(realvar & mu,realvar & nu);


//--------------------------- TYPE DEFINITIONS ---------------------------------

typedef std::list<realvar>::iterator __EXPORT_TYPE variterator;

}   // end: namespace realob


#endif

