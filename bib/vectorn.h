// Date: 17.12.97


#if !defined (VECTORN_INCLUDED)

#define VECTORN_INCLUDED


#include<math.h>
#include<vector>
#include<realobs.h>


//------------------------------------------------------------------------------
//----------------------------- CLASS vectornum --------------------------------
//------------------------------------------------------------------------------


template <class T>
class vectornum : public vector<T>
  {


  public:


  //---------------------------- PUBLIC FUNCTIONS ------------------------------

  // DEFAULT CONSTRUCTOR
  // TASK: creates a vector with no elements

  vectornum(void) : vector<T>()
	 {
	 }

  // CONSTRUCTOR
  // TASK: creates a vector with dimension nr

  vectornum(int nr) : vector<T> (nr)
	 {
	 }

  // CONSTRUCTOR
  // TASK: creates a nr dimensional vector with all elements equal to v

  vectornum(int nr,T v) : vector<T> (nr,v)
	 {
	 }

  // COPY CONSTRUCTOR (copies a vectornum)

  vectornum(const vectornum & v) : vector<T> ( vector<T>(v) )
	 {
	 }

  // COPY CONSTRUCTOR (copies a vector)

  vectornum(const vector<T> & v) : vector<T> (v)
	 {
	 }

  // OVERLOADED ASSIGNMENT OPERATORS

  // OVERLOADED ASSIGNMENT OPERATOR (assigns a vectornum)

  const vectornum<T> & operator=(const vectornum<T> & v) 
	 {
	 vector<T>::operator=( vector<T>(v) );
	 return *this;
	 }

  //  OVERLOADED ASSIGNMENT OPERATOR (assigns a vector to a vectornum)

  const vectornum<T> & operator=(const vector<T> & v)
	 {
	 vector<T>::operator=(v);
	 return *this;
	 }

  // OVERLOADED ASSIGNMENT OPERATOR
  // TASK: returns a vector with all elements equal to v

  const vectornum<T> & operator=(const T v);

  // BINARY ADDITION OPERATORS

  // BINARY ADDITION OPERATOR
  // TASK: adds to vectornum objects (element by element)

  friend vectornum<T> operator+(const vectornum<T> & v1,const vectornum<T> & v2);

  // BINARY ADDITION OPERATOR
  // TASK: creates a vector with the same length as v2 and all elements
  //       equal to v1 and finaly adds the two vectors (element by element)

  friend vectornum<T> operator+(const T & v1,const vectornum<T> & v2);

  // BINARY ADDITION OPERATOR
  // TASK: creates a vector with the same length as the calling vector
  //       and all elements equal to v1
  //       and finaly adds the two vectors (element by element)

  vectornum<T> operator+(const T v2)
	 {
	 return v2+(*this);
	 }

  // UNARY + OPERATOR

  vectornum<T> operator+()
	 {
	 return *this;
	 }

  // BINARY MINUS OPERATORS
  // (see BIARY ADDITION OPERATORS for description)

  friend vectornum<T> operator-(const vectornum<T> & v1,const vectornum<T> & v2);

  friend vectornum<T> operator-(const T & v1,const vectornum<T> & v2);

  friend vectornum<T> operator-(const vectornum<T> & v1,const T & v2);

  // UNARY - OPERATOR

  vectornum<T> operator-()
	 {
	 return T(-1)*(*this);
	 }

  // BINARY MULTIPLIKATION OPERATORS
  // (see BIARY ADDITION OPERATORS for description)

  friend vectornum<T> operator*(const vectornum<T> & v1,const vectornum<T> & v2);

  friend vectornum<T> operator*(const T & v1,const vectornum<T> & v2);

  friend vectornum<T> operator*(vectornum<T> & v1, const T & v2)
	 {
	 return v2*v1;
	 }

  // BINARY DIVISION OPERATORS
  // (see BIARY ADDITION OPERATORS for description)

  friend vectornum<T> operator/(const vectornum<T> & v1,const vectornum<T> & v2);

  friend vectornum<T> operator/(const T & v1,const vectornum<T> & v2);

  vectornum<T> operator/(const T & v2) const
	 {
	 return (*this) * (T(1)/v2);
	 }

  // OVERLOADED COMPARISON OPERATORS

  // > OPERATORS

  // > OPERATOR
  // TASK: compares the two vectors element by element
  //        returns 1 in the i. row of the resulting vector, if
  //        v1[i] > v2[i], else 0

  friend vectornum<T> operator>(const vectornum<T> & v1,const vectornum<T> & v2);

  friend vectornum<T> operator>(const vectornum<T> & v1,const T & v2);

  friend vectornum<T> operator>(const T & v1, const vectornum<T> & v2);

  // < OPERATORS
  // (see > OPERATORS for description)

  friend vectornum<T> operator<(const vectornum<T> & v1,const vectornum<T> & v2);

  friend vectornum<T> operator<(const vectornum<T> & v1,const T & v2);

  friend vectornum<T> operator<(const T & v2,const vectornum<T> & v1);

  // == OPERATORS
  // (see > OPERATORS for description)

  vectornum<T> isequal(vectornum<T> & v2);

  friend vectornum<T> isequal(vectornum<T> & v1,const T v2);

  friend vectornum<T> isequal(T v2, vectornum<T> & v1);

  // != OPERATORS
  // (see > OPERATORS for description)

  vectornum<T> isnotequal(vectornum<T> & v2);

  friend vectornum<T> isnotequal(vectornum<T> & v1,const T v2);

  friend vectornum<T> isnotequal(T v2, vectornum<T> & v1);

  // >= OPERATORS
  // (see > OPERATORS for description)

  friend vectornum<T> operator>=(const vectornum<T> & v1,const vectornum<T> & v2);

  friend vectornum<T> operator>=(vectornum<T> & v1,const T v2);

  friend vectornum<T> operator>=(const T v2, const vectornum<T> & v1);

  // <= OPERATORS
  // (see > OPERATORS for description)

  friend vectornum<T> operator<=(const vectornum<T> & v1,const vectornum<T> & v2);

  friend vectornum<T> operator<=(vectornum<T> & v1,const T v2);

  friend vectornum<T> operator<=(const T v2, vectornum<T> & v1);

  // OVERLOADED LOGICAL OPERATORS
  // (see > OPERATORS for description)

  friend vectornum<T> operator||(const vectornum<T> & v1,const vectornum<T> & v2);

  friend vectornum<T> operator&&(const vectornum<T> & v1,const vectornum<T> & v2);

  // OVERLOADED << OPERATOR

  friend ostream & operator<<(ostream & out,vectornum<T> & v);

  // FUNCTION: applied
  // TASK: applies function func (which must be implemented for type T)
  //       to all elements of the vector and returns the resulting vector

  vectornum<T> applied(T (*func)(const T &));

  // FUNCTION: applied
  // TASK: applies function func (which must be implemented for type T)
  //       to all elements of the vector and returns the resulting vector

  vectornum<T> applied(T (*func)(T));

  // FUNCTION: sort
  // TASK: assumes that the 'true' order of the elements of the vector is
  //       stored in 'index', performs a indexsort for rows left to right

  void sort(vector<int> & index,int left,int right);

  // FUNCTION: clear
  // TASK: deletes all elements (length = 0)

  void clear()
	 {
	 if (size() > 0)
		erase(begin(),end());
	 }


  };

/*

//------------------------------------------------------------------------------
//---------------- FUNCTIONS THAT CAN BE APPLIED TO VECTORNUM ------------------
//------------------------------------------------------------------------------


// | POSSIBLE FUNCTIONS                ABBREVIATION
// | ----------------------------------------------
// | square root                      | sqrt
// | absolute value                   | abs
// | expontential                     | exp
//	| cosinus                          | sin
//	| sinus                            | cos
//	| logarithm                        | log
//	| logarithm with basis 10          | log10
// |
// | floor function                   | floor


// FUNCTION: sqrt
// TASK: returns sqrt(o)

template <class T>
vectornum<T> sqrt(const vectornum<T> & o)
  {
  return o.applied(sqrt);
  }


  // FUNCTION: abs
// TASK: returns abs(o)

template <class T>
vectornum<T> abs(const vectornum<T> & o)
  {
  return o.applied(fabs);
  }


// FUNCTION: exp
// TASK: returns exp(o)

template <class T>
vectornum<T> exp(const vectornum<T> & o)
  {
  return o.applied(exp);
  }


// FUNCTION: cos
// TASK: returns cos(o)

template <class T>
vectornum<T> cos(const vectornum<T> & o)
  {
  return o.applied(cos);
  }


// FUNCTION: sin
// TASK: returns sin

template <class T>
vectornum<T> sin(const vectornum<T> & o)
  {
  return o.applied(sin);
  }


// FUNCTION: log
// TASK: returns log(o)

template <class T>
vectornum<T> log(const vectornum<T> & o)
  {
  return o.applied(log);
  }


// FUNCTION: log10
// TASK: returns log10(o)

template <class T>
vectornum<T> log10(const vectornum<T> & o)
  {
  return o.applied(log10);
  }


// FUNCTION: pow
// TASK: returns v1^v2

template<class T>
vectornum<T> pow(vectornum<T> & v1,vectornum<T> & v2)
  {
  vectornum<T>::iterator pos1,pos2,pos3;
  vectornum<T> h(v1.size());
  for (pos1=v1.begin(),pos2=v2.begin(),pos3=h.begin();
		 pos1 != v1.end();++pos1,++pos2,++pos3)
	 *pos3 = pow(*pos1,*pos2);
  return h;
  }


// FUNCTION: floor
// TASK: returns floor(o)
//       floor(x) is the greatest integer value, not greater than x

template <class T>
vectornum<T> floor(const vectornum<T> & o)
  {
  return o.applied(floor);
  }

*/

// --------------------  SOME TYPEDEFINITIONS ----------------------------------


typedef vectornum<double> vectord;
typedef vectornum<float> vectorf;
typedef vectornum<int> vectori;
typedef vectornum<realob::realobs> vectorrealobs; 


#include<vectorn.cpp>


#endif

