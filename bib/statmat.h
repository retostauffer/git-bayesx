// DATE: 02.01.99


#if !defined(STATMATRIX_INCLUDED)

#define STATMATRIX_INCLUDED

#include <tmatrix.h>
#include <vector>
# include <clstring.h>

using std::vector;

//------------------------------------------------------------------------------
//---------------------------- CLASS: statmatrix -------------------------------
//------------------------------------------------------------------------------

class SparseMatrix;

// AENDERUNG (Eva)
class adja;

template <class T>
class statmatrix : public Matrix<T>

  {

  public:

  // ---------------------------- CONSTRUCTORS ---------------------------------

  // DEFAULT CONSTRUCTOR

  statmatrix(void) : Matrix<T> () {}

  // CONSTRUCTOR 1

  statmatrix(unsigned rows, unsigned cols = 1) : Matrix<T> (rows,cols) {}

  // CONSTRUCTOR 2

  statmatrix(unsigned rows, unsigned cols, const T & init)
				 : Matrix<T> (rows, cols, init) {}

  // CONSTRUCTOR 3

  statmatrix(const SparseMatrix & m);

  // CONSTRUCTOR 4

  statmatrix(const vector<T> & v);

  // COPY CONSTRUCTORS

  statmatrix(const Matrix<T> & m) : Matrix<T> (m) {}
  statmatrix(const statmatrix & s) : Matrix<T> (s) {}

  // OVERLOADED ASSIGNMENT OPERATORS

  const statmatrix & operator=(const SparseMatrix & m);

  // AENDERUNG (Eva)
  // OVERLOADED ASSIGNMENT OPERATOR
  statmatrix   operator= (const adja & a);


  // ------------------------- PUBLIC FUNCTIONS --------------------------------

  // FUNCTION: assign
  // TASK: assigns the elements of A to the elements of the calling matrix
  //        faster than B = A or B.putBlock(A,0,0,B.rows,B.cols())

  void assign(const statmatrix & A);

  // FUNCTION: plus
  // TASK: assigns A+B to the calling matrix
  //       faster than C = A+B

  void plus(const statmatrix & A,const statmatrix & B);

  // FUNCTION: plus
  // TASK: adds A to the calling matrix

  void plus(const statmatrix & A);

  // FUNCTION: minus
  // TASK: assigns A-B to the calling matrix
  //       faster than C = A-B

  void minus(const statmatrix & A,const statmatrix & B);

  // FUNCTION: minus
  // TASK: substracts the colA th column of A from the colB th column B and
  //       assigns the result to the calling matrix (must be column vector)

  void minus(const statmatrix & A,const statmatrix & B,const unsigned & colA,
             const unsigned & colB);

  // FUNCTION: mult
  // TASK: assigns A*B to the calling matrix

  void mult(const statmatrix & A,const statmatrix & B);

  // FUNCTION: addmult
  // TASK: computes A*B and adds the result to the calling matrix

  void addmult(const statmatrix & A,const statmatrix & B);

  // FUNCTION: addmult
  // TASK: computes A*B and adds the result to the calling matrix
  //       A is assumed to be symmetric with elements stored in the lower
  //       triangular

  void addmultsym(const statmatrix & A,const statmatrix & B);

  // FUNCTION: inverse
  // TASK: computes the inverse of the calling matrix

  statmatrix<T> inverse(void);

  // FUNCTION: multdiagback
  // TASK: multiplies the calling matrix with the diagonal matrix D (from the
  //       right). The diagonal elements of D are stored in d

  void multdiagback(const statmatrix & d);

  // FUNCTION: multdiagfront
  // TASK: multiplies the calling matrix with the diagonal matrix D (from the
  //       left). The diagonal elements of D are stored in d

  void multdiagfront(const statmatrix & d);

  // FUNCTION: multdiagfront
  // TASK: multiplies A with the diagonal matrix D (from the left) and assigns
  //       the result to the calling matrix. The diagonal elements of D are
  //       stored in d

  void multdiagfront(const statmatrix & A, const statmatrix & d);

  // FUNCTION: addtodiag
  // TASK: Adds the elements in d to the diagonal elements of the calling matrix
  //       The first element of d is added to A(first,first), the last element
  //       of d is added to A(last-1,last-1)
  //       The calling matrix is assumed to have rows==cols

  void addtodiag(const statmatrix & d, unsigned first, unsigned last);

  // FUNCTION: subfromdiag
  // TASK: Subcracts the elements in d to the diagonal elements of the calling
  //       matrix. The first element of d is subtracted from A(first,first),
  //       the last element of d is subtracted from A(last-1,last-1)
  //       The calling matrix is assumed to have rows==cols

  void subfromdiag(const statmatrix & d, unsigned first, unsigned last);

  // FUNCTION: elemmult
  // TASK: computes the elementwise product of the calling matrix and A

  void elemmult(const statmatrix & A);

  // FUNCTION: elemquot
  // TASK: computes the elementwise ratio B/A and assigns it to the calling
  //       matrix B

  void elemquot(const statmatrix & A);

  // FUNCTION: weightedsscp
  // TASK: Computes the weighted sums of squares and crossproducts
  //       matrix X'WX

  void weightedsscp(const statmatrix & X, const statmatrix & w);

  // FUNCTION: weightedsscp2
  // TASK: Computes the weighted sums of squares and crossproducts
  //       matrix (X Z)'W(X Z)

  void weightedsscp2(const statmatrix & X, const statmatrix & Z,
                     const statmatrix & w);

  // FUNCTION: weightedsscp_resp
  // TASK: Computes X'W y

  void weightedsscp_resp(const statmatrix & X, const statmatrix & y,
                          const statmatrix & w);

  // FUNCTION: weightedsscp_resp2
  // TASK: Computes (X Z)'W y

  void weightedsscp_resp2(const statmatrix & X, const statmatrix & Z,
                     const statmatrix & y, const statmatrix & w);

  // ------------------ functions for sorting a column -------------------------

  // FUNCTION: sort
  // TASK: sorts the statmatrix between the 'start' th and the 'ende' th row
  //       according to its 'col' th column

  void sort (int start ,int ende,int col);

  // FUNCTION: sort
  // TASK: sorts the 'col' th column of the statmatrix between the 'start' th
  //       and the 'ende' th row

  void sortcol (int start ,int ende,int col);

  // FUNCTION: indexinit
  // TASK: initializes a column vector with elements (i,0) = i

  void indexinit (void);

  // FUNCTION: indexsort
  // TASK: index sort of the 'col' th column within the 'start' th and 'ende'
  //       th row. after sorting 'indexcol' of index contains the rank vector
  //       of the 'col' th column of the calling matrix

  void indexsort (statmatrix<int> & index,int start,int ende,
                  int col,int indexcol) const;


//------------------------- end: sorting a column ------------------------------


  // FUNCTION: sum
  // TASK: returns the sum of the elements of the col. column of the calling
  //       matrix

  T  sum (const unsigned & col) const;

  // FUNCTION: sum2
  // TASK: returns the sum of the squared elements of the col th column of the
  //       calling matrix

  T  sum2 (const unsigned & col) const;

  T  sum2 (const unsigned & col,const statmatrix<T> & weight) const;

  // FUNCTION: sum2
  // TASK: returns a column vector whose elements are the squared column sums of
  //       the calling matrix

  statmatrix<T> sum2();

  // FUNCTION: sumcomplete
  // TASK: returns the sum of the elements of the calling
  //       matrix

  T  sumcomplete(void) const;

  // FUNCTION: sum
  // TASK: returns a column vector whose elements are the column sums of the
  //       calling matrix

  statmatrix<T> sum() const;

  // FUNCTION: mean
  // TASK: computes the mean of column 'col' of the calling matrix

  T mean(const unsigned & col) const
    {
    return sum(col)/T(rows());
    }

  T mean(const unsigned & col,const statmatrix<T> & weight) const;

  // FUNCTION: meancomplete
  // TASK:  returns the mean of the elements of the calling matrix

  T meancomplete(void) const
    {
    return sumcomplete()/T(rows());
    }

  // FUNCTION: norm
  // TASK: computes the euclidean norm of column col

  T norm(unsigned col);

  // FUNCTION: norm
  // TASK: returns a column vector of the norms of the columns of the calling
  //       matrix

  statmatrix<T> norm();

  // FUNCTION: mean
  // TASK: computes the mean of the columns of the calling matrix
  //       returns a column vector of means

  statmatrix<T> mean() const;

  // FUNCTION: var
  // TASK: computes the variance of column 'col' of the calling matrix

  T var(const unsigned & col) const;

  T var(const unsigned & col,const statmatrix<double> & weight) const;

  // FUNCTION: min
  // TASK: returns the minimum of column 'col'

  T min(const unsigned & col) const;

  // FUNCTION: max
  // TASK: returns the maximum of column 'col'

  T max(const unsigned & col) const;


  // FUNCTION: quantile
  // TASK: returns the 'percent' (0 < percent < 100) percent quantile
  //       of the 'col' th column

  T quantile  (const T & percent,const unsigned & col) const;

  // FUNCTION: quantile
  // TASK: computes the 'percent' percent quantile of the columns of the
  //       calling matrix, returns a column vektor with the quantiles

  statmatrix<T> quantile  (T percent);

  // FUNCTION: autocorr
  // TASK: computes the autocorrelation function for the 'col'. column and
  //			  lag 'lag'

  T             autocorr  (const unsigned & lag,const unsigned & col) const;

  // FUNCTION: autocorr
  // TASK: computes the autocorrelation function for lags 1 to 'lag' for
  //       all columns
  //       returns a  lag x column matrix of autocorrelations

  statmatrix<T> autocorr  (const unsigned & lag) const;

  // FUNCTION: autocorr
  // TASK: computes the autocorrelation function for lags 'beginlag' to 'endlag'
  //       for column 'col'
  //       returns a column vector of autocorrelations

  statmatrix<T> autocorr (const unsigned & beginlag,const unsigned & endlag,
                          const unsigned & col) const;


  // FUNCTION: cov
  // TASK: computes the covariance matrix of the calling matrix

  statmatrix<T> cov();

  statmatrix<T> corr();

  T compute_quadform(const statmatrix<T> & x,const unsigned & c=0);

  	//	Datenzeiger zugreifbar

  T *getV() const { return m_v; }


   statmatrix<T> strike (unsigned int k);

   statmatrix<T> get_cov_iX (int i, int j);

   statmatrix<T> partial_var(void);

  };


typedef statmatrix<double> datamatrix;

// FUNCTION: multdiagback
// TASK: Computes X*D, where D is a diagonal matrix whose elements are stored
//       in d

statmatrix<double> multdiagback(datamatrix X, const datamatrix & d);

// FUNCTION: multdiagback
// TASK: Computes D*X, where D is a diagonal matrix whose elements are stored
//       in d

statmatrix<double> multdiagfront(datamatrix X, const datamatrix & d);

// FUNCTION: diffmat
// TASK: Computes a difference matrix with dimension (d-k x d), where k is the
//       difference order (k=1,2)

statmatrix<double> diffmat(const int k, const int d);

// FUNCTION: diffmat
// TASK: Computes a difference matrix with dimension (d-k x d), where k is the
//       difference order (k=0,1,2,...,d-1)

statmatrix<double> diffmat_k(const int k, const int d);

// FUNCTION: weighteddiffmat
// TASK: Computes a weighted difference matrix with dimension (d-k x d), where
//       k is the difference order and d is given by weight.size()

statmatrix<double> weighteddiffmat(const int k, const vector<double> & weight);

// FUNCTION: seasonalfactor
// TASK:     computes the factor of the penalty matrix for a seasonal effect
//           with period per and s time periods.

statmatrix<double> seasonalfactor(const unsigned & per, const unsigned & s);

// FUNCTION: seasonalX
// TASK:     computes the deterministic part of a seasonal effect
//           with period per and s time periods.

statmatrix<double> seasonalX(const unsigned & per, const unsigned & s);

// FUNCTION: rotate
// TASK:     help-function for the computation of eigenvalues and -vectors

void rotate(statmatrix<double> & a, const double & s, const double & tau,
            const int & i, const int & j, const int & k, const int & l);

// FUNCTION: tridiag
// TASK:     reduces the (symmetric) matrix a to tridiagonal form. On output,
//           the diagonal elements are stored in d and the subdiagonal elements
//           are stored in e. d and e are assumed to be n x 1 matrices.
//           The first element of e is zero on output.
// NOTE:     a is replaced by the orthogonal matrix effecting its transformation !!

void tridiag(statmatrix<double> & a, statmatrix<double> & d,
             statmatrix<double> & e);

// FUNCTION: eigentridiag
// TASK:     computes the eigenvalues and -vectors of a tridiagonal matrix. d
//           and e are assumed to be n x 1 matrices containing the diagonal
//           elements and the subdiagonal elements (with e(0,0) arbitrary)
//           respectively. If the eigenvectors of an originally tridiagonal
//           matrix are desired, z is input as the identity matrix, otherwise
//           z contains the output a from tridiag. On output d contains the
//           eigenvalues and z contains the eigenvectors. e ist destroyed on
//           output.

bool eigentridiag(statmatrix<double> & d, statmatrix<double> & e,
                  statmatrix<double> & z);

// FUNCTION: pythag
// TASK:     computes a^2+b^2 without destructive underflow or overflow

double pythag(const double & a, const double & b);

double sqr(const double & a);

double SIGN(const double & a, const double & b);


// FUNCTION: eigen
// TASK:     computes eigenvalues and eigenvectors of a and
//           stores them in values or vectors respectively.
//           The return value gives the number of iterations, that where needed
//           in the computation. If the return value equals 50, no convergence
//           could be achieved.
// NOTE:     a is modified in the computation of the eigenvalues and -vectors!!

int eigen(statmatrix<double> & a, statmatrix<double> & values,
           statmatrix<double> & vectors);

// FUNCTION: eigen2
// TASK:     computes eigenvalues and eigenvectors of a. On output, eigenvectors
//           are stored in a and eigenvalues are stored in d. If the return
//           value is false, no convergence could be achieved

bool eigen2(statmatrix<double> & a, statmatrix<double> & d);

// FUNCTION: eigensort
// TASK:     sorts the eigenvalues in values into descending order and
//           rearranges the columns of vectors correspondingly

void eigensort(statmatrix<double> & values, statmatrix<double> & vectors);

// FUNCTION: kronecker
// TASK:     computes the kronecker product of the matrices Aand B

statmatrix<double> kronecker(const statmatrix<double> & A, const statmatrix<double> & B);

void compare(const datamatrix & ref, const datamatrix & neu, double limit, unsigned col, const ST::string & colname, vector<ST::string> & out);
void compare_nonp(const ST::string & ref, const ST::string & neu, double limit, vector<ST::string> & out);


#include<statmat.cpp>


#endif
