// tpremat.cc 3.3 97/08/07 17:01:19 
//
// Muster fuer die outline-Funktionen der template-Klasse 
// PreMatrix
//
// SCHALTER :
//
// SAVE_ALGORITHMS
//
// Mit dem Schalter kann zwischen einer "theorienahen" und einer
// Speicherdarstellungsnahen Variante der Algorithmen gewaehlt
// werden. Die an der Speicherdarstellung orientierten Algorithmen
// sind fuer Fliesskommadatentypen den "theorienahen" Algorithmen
// nicht deutlich ueberlegen (Sie sparen einige Prozent Rechenzeit)
// Interessanter sind diese Algorithmen bei Festkommaarithmetik,
// wenn der Overhead der Algorithmen gegenueber der Rechenzeit
// ins Gewicht faellt.

#include <tlinklst.h>
#include <tarray.h>
#include <math.h>
#include <limits.h>

#if defined(_MSC_VER2)
#include <strstrea.h>
#else
#include <strstream.h>
#endif
#include <string.h>

// PreMatrix operator+( const PreMatrix &m ) const
//
// (Elementweise) NonRealMatrixaddition
//
// Parameter :
//
// m - zweiter Operand der NonRealMatrixaddition
//
// Ergebnis :
//
// (*this) + m
//
// Voraussetzungen :
//
// (*this) und m stimmen in ihrer Zeilen und ihrer Spaltendimension
// ueberein.

template <class T>
PreMatrix<T> PreMatrix<T>::operator+( const PreMatrix<T> &m ) const
{
	assert( !operator!( ) );
	assert( m );
	assert( m.rows( ) == rows( ) );
	assert( m.cols( ) == cols( ) );

	if ( operator!( ) || !m || m.rows( ) != rows( ) ||
		m.cols( ) != cols( )  )
		return PreMatrix<T>( 0 );

	PreMatrix<T> result( rows( ), cols( ) );
	assert( result );
	if ( !result )
		return result;

#if defined( SAVE_ALGORITHMS )

	// A - PreMatrix an die die Botschaft gerichtet ist (*this) (m x n)
	// B - Parameter der Botschaft (m)  (m x n)
	// R - Ergebnis von (*this) + m   (m x n)
	//
	// r_{ij} = a_{ij} + b_{ij}

	register unsigned i, j;
	for ( i = 0; i < rows( ); i++ )
		for( j = 0; j < cols( ); j++ )
			result( i, j ) = get( i, j ) + m.get( i, j );

#else

	// An der Speicherdarstellung orientierte Variante:
	// Je ein Zeiger laeuft durch jede der drei Matrizen,
	// die ja von gleicher Gestalt sind.

	unsigned size = rows( ) * cols( );

	// Laufvariablen (Zeiger und Z�aehler)

	T *workA = getV( );
	T *workB = m.getV( );
	T *workR = result.getV( );
	register unsigned i;

	// NonRealMatrixaddition

	for ( i = 0;
			i < size;
			i++, workA++, workB++, workR++ )
		*workR = *workA + *workB;

#endif

	return result;
}


// unaerer Operator +
//
template <class T>
PreMatrix<T> PreMatrix<T>::operator+()
{
	if ( !(*this) ) return PreMatrix<T>(0);
	return PreMatrix<T>( *this );
}

 
// PreMatrix operator-( const PreMatrix &m ) const
//
// (Elementweise) NonRealMatrixsubtraktion
//
// Parameter :
//
// m - zweiter Operand der NonRealMatrixsubtraktion
//
// Ergebnis :
//
// (*this) - m
//
// Voraussetzungen :
//
// (*this) und m stimmen in ihrer Zeilen und ihrer Spaltendimension
// ueberein.

template <class T>
PreMatrix<T> PreMatrix<T>::operator-( const PreMatrix<T> &m ) const
{
	assert( !operator!( ) );
	assert( m );
	assert( m.rows( ) == rows( ) );
	assert( m.cols( ) == cols( ) );

	if ( operator!( ) || !m || m.rows( ) != rows( ) ||
		m.cols( ) != cols( )  )
		return PreMatrix<T>( 0 );

	PreMatrix<T> result( rows( ), cols( ) );
	assert( result );
	if ( !result )
		return result;

#if defined( SAVE_ALGORITHMS )

	// A - PreMatrix an die die Botschaft gerichtet ist (*this) (m x n)
	// B - Parameter der Botschaft (m)  (m x n)
	// R - Ergebnis von (*this) - m  (m x n)
	//
	// r_{ij} = a_{ij} - b_{ij}

	register unsigned i, j;
	for ( i = 0; i < rows( ); i++ )
		for( j = 0; j < cols( ); j++ )
			result( i, j ) = get( i, j ) - m.get( i, j );
#else

	// An der Speicherdarstellung orientierte Variante:
	// Je ein Zeiger laeuft durch jede der drei Matrizen,
	// die ja von gleicher Gestalt sind.

	unsigned size = rows( ) * cols( );

	// Laufvariablen

	T *workA = getV( );
	T *workB = m.getV( );
	T *workR = result.getV( );
	register unsigned i;

	// NonRealMatrixsubtraktion

	for ( i = 0;
			i < size;
			i++, workA++, workB++, workR++ )
		*workR = *workA - *workB;

#endif

	return result;
}


// unaerer Operator -
//
template <class T>
PreMatrix<T> PreMatrix<T>::operator-()
{
	if ( !(*this) ) return PreMatrix<T>(0);
	return PreMatrix<T>( (*this) * T(-1) );
}


// PreMatrix operator*( const PreMatrix &m ) const
//
// NonRealMatrixmultiplikation
//
// Parameter :
//
// m - zweiter Operand der NonRealMatrixmultiplikation
//
// Ergebnis :
//
// (*this) * m
//
// Voraussetzungen :
//
// Die Spaltendimension des ersten Operanden stimmt mit der Zeilen-
// dimension des zweiten Operanden ueberein, d.h. die Matrizen sind
// multiplizierbar.

template <class T>
PreMatrix<T> PreMatrix<T>::operator*( const PreMatrix<T> &m ) const
{
	assert( !operator!( ) );
	assert( m );
	assert( m.rows( ) == cols( ) );

	if ( operator!( ) || !m || m.rows( ) != cols( ) )
		return PreMatrix<T>( 0 );

	PreMatrix<T> result( rows( ), m.cols( ) );
	assert( result );
	if ( !result )
		return PreMatrix<T>( 0 );

	T sum;

#if defined( SAVE_ALGORITHMS )

	// A - PreMatrix an die die Botschaft gerichtet ist (*this) (m x l)
	// B - Parameter der Botschaft (m) (l x n)
	// R - Ergebnis von (*this) - m  (m x n)
	//
	// r_{ij} = \sum_{k=1}^l a_{ik} * b_{kj}

	register unsigned i, j, k;
	for ( i = 0; i < result.rows( ); i++ )
	{
		for( j = 0; j < result.cols( ); j++ )
		{
			sum = T( 0 );
			for ( k = 0; k < cols( ); k++ )
			{
				T aik = get( i, k );
				T bkj = m.get( k, j );
				if ( aik != T( 0 ) && bkj != T( 0 ) )
					sum += aik * bkj;
			}
			result( i, j ) = sum;
		}
	}

#else

	// Schnelle NonRealMatrixmultiplikation

	// Ergebnis-Zeiger

	T *workR = result.getV( );

	// "Dritte" Dimension

	unsigned n = result.cols( );

	// Elemente der Ergebnismatrix

	unsigned size = result.rows( ) * n;
	register unsigned i, k;

	// Ergebnismatrix elementweise besetzen

	for (  i = 0; i < size; i++ )
	{
		sum = T( 0 );

		// workA zeigt auf den Anfang der Zeile i/n

		T *workA = getV( ) +  (i / n) * cols( );

		// workB zeigt auf das entsprechende Element in der ersten Zeile

		T *workB = m.getV( ) +  i % n;

		// Beim Durchlaufen der beiden Operandenmatrizen wird der Zeiger
		// workA in jedem Schritt um ein Element, der Zeiger workB in
		// jedem Schritt um eine Zeile weitergezaehlt.

		for (k = 0; k < cols(); ++k, ++workA, workB += n )
			if (!(*workA == T(0) || *workB == T(0)))
				sum += *workA * *workB;

		// Eintrag in die Ergebnismatrix

		*workR++ = sum;
	}

#endif

	return result;
}

// PreMatrix operator*( const T v ) const
//
// Multiplikation einer PreMatrix mit einem Skalar
//
// Parameter :
//
// v - Skalar, mit dem die PreMatrix multipliziert werden soll
//
// Ergebnis :
//
// Mit dem Skalar multiplizierte PreMatrix

template <class T>
PreMatrix<T> PreMatrix<T>::operator*( const T v ) const
{
	assert( !operator!( ) );
	if ( operator!( ) )
		return PreMatrix<T>( 0 );

	if ( v == T( 1 ) )
		return PreMatrix<T>( *this );
   else if ( v == T( 0 ) )
		return PreMatrix<T>( rows( ), cols( ), T( 0 ) );

	//	PreMatrix uninitialisiert anlgen und mit der PreMatrix, deren 
	//	Methode aufgerufen wurde und dem Argument der Methode 
	//	besetzten (spart die initialisierung mit Werten, die 
	//	im naechsten Schritt ueberschrieben werden.

	PreMatrix<T> result( rows( ), cols( ) );
	assert( result );
	if ( !result )
		return PreMatrix<T>( 0 );

#if defined( SAVE_ALGORITHMS )

	// A - PreMatrix an die die Botschaft gerichtet ist (*this) (m x n)
	// v - Parameter der Botschaft
	// R - Ergebnis von (*this) * v  (m x n)
	//
	// r_{ij} = a_{ik} * v

	register unsigned i, j;
	for ( i = 0; i < rows( ); i++ )
		for( j = 0; j < cols( ); j++ )
			result( i, j ) = get( i, j ) * v;

#else

	unsigned size = rows( ) * cols( );
	register unsigned i;
	T *work  = result.getV( );
	T *workA = getV( );

	for ( i = 0; i < size; i++, work++, workA++ )
		*work = *workA * v;

#endif

	return result;
}

// PreMatrix operator/( const T v ) const
//
// Multiplikation einer PreMatrix mit dem Kehrwert eines Skalars
//
// Parameter :
//
// v - Skalar, mit dessen Kehrwert die PreMatrix multipliziert werden soll
//
// Ergebnis :
//
// Mit dem Kehrwert des Skalars multiplizierte PreMatrix

template <class T>
PreMatrix<T> PreMatrix<T>::operator/( const T v ) const
{
	assert(!operator!());	
	assert(!(v == T(0)));

	if ( operator!( ) || v == T( 0 ) )
		return PreMatrix<T>( 0 );

	if ( v == T( 1 ) )
		return PreMatrix<T>( *this );

	PreMatrix<T> result( rows( ), cols( ) );
	assert( result );
	if ( !result )
		return result;

#if defined( SAVE_ALGORITHMS )

	// A - PreMatrix an die die Botschaft gerichtet ist (*this) (m x n)
	// v - Parameter der Botschaft 
	// R - Ergebnis von (*this) * v^-1  (m x n)
	//
	// r_{ij} = a_{ik} * v

	register unsigned i, j;
	for ( i = 0; i < rows( ); i++ )
		for( j = 0; j < cols( ); j++ )
			result( i, j ) = get( i, j ) / v;

#else

	unsigned size = rows( ) * cols( );
	register unsigned i;
	T *work  = result.getV( );
	T *workA = getV( );
	for ( i = 0; i < size; i++, work++, workA++ )
		*work = *workA / v;

#endif

	return result;
}

template <class T>
const PreMatrix<T> &PreMatrix<T>::operator+=( const PreMatrix<T> &m )
{
	assert( !operator!( ) );
	assert( m );
	assert( m.rows( ) == rows( ) );
	assert( m.cols( ) == cols( ) );

	if ( operator!( ) )
	{
		return *this;
	}
	else if ( !m || m.rows( ) != rows( ) ||  m.cols( ) != cols( )  )
	{
		*this = PreMatrix<T>( 0 );
		return *this;
	}

#if defined( SAVE_ALGORITHMS )

	register unsigned i, j;
	for ( i = 0; i < rows( ); i++ )
		for ( j = 0; j < cols( ); j++ )
			operator()( i, j ) += m.get( i, j );

#else
	
	unsigned size = rows( ) * cols( );
	register unsigned i;
	T * work = getV( );
	T * workB = m.getV( );
	for ( i = 0; i < size; i++ )
		*work++ += *workB++;

#endif

	return *this;
}

template <class T>
const PreMatrix<T> &PreMatrix<T>::operator-=( const PreMatrix<T> &m )
{
	assert( !operator!( ) );
	assert( m );
	assert( m.rows( ) == rows( ) );
	assert( m.cols( ) == cols( ) );

	if ( operator!( ) )
	{
		return *this;
	}
	else if ( !m || m.rows( ) != rows( ) || m.cols( ) != cols( )  )
	{
		*this = PreMatrix<T>( 0 );
		return *this;
	}

#if defined( SAVE_ALGORITHMS )

	register unsigned i, j;
	for ( i = 0; i < rows( ); i++ )
		for ( j = 0; j < cols( ); j++ )
			operator()( i, j ) += m.get( i, j );
#else

	unsigned size = rows( ) * cols( );
	register unsigned i;
	T * work  = getV( );
	T * workB = m.getV( );
	for ( i = 0; i < size; i++ )
		*work++ -= *workB++;

#endif

	return *this;
}

template <class T>
const PreMatrix<T> &PreMatrix<T>::operator*=( const T v )
{
	assert( !operator!( ) );

	if ( operator!( ) )
		return *this;

#if defined( SAVE_ALGORITHMS )

	register unsigned i, j;
	for ( i = 0; i < rows( ); i++ )
		for ( j = 0; j < cols( ); j++ )
			operator()( i, j ) *= v;

#else

	unsigned size = rows( ) * cols( );
	register unsigned i;
	T * work  = getV( );
	for ( i = 0; i < size; i++ )
		*work++ *= v;

#endif

	return *this;
}

template <class T>
const PreMatrix<T> &PreMatrix<T>::operator/=( const T v )
{
	assert( !operator!( ) );

	if ( operator!( ) )
		return *this;

	assert( v != T( 0 ) );

	if ( v == T( 0 ) )
	{
		*this = PreMatrix<T>( 0 );
		return *this;
	}

#if defined( SAVE_ALGORITHMS )

	register unsigned i, j;
	for ( i = 0; i < rows( ); i++ )
		for ( j = 0; j < cols( ); j++ )
			operator()( i, j ) /= v;

#else

	unsigned size = rows( ) * cols( );
	register unsigned i;
	T * work  = getV( );
	for ( i = 0; i < size; i++ )
		*work++ /= v;

#endif

	return *this;
}


template <class T>
PreMatrix<T> PreMatrix<T>::transposed( void ) const
{
	if ( operator!( ) )
		return PreMatrix<T>( 0 );

	PreMatrix<T> result( cols( ), rows( ) );
	if ( !result )
		return result;

	register unsigned i, j;
	for ( i = 0; i < rows( ); i++ )
		for ( j = 0; j < cols( ); j++ )
			result( j, i ) = get( i, j );
	return result;
}

template <class T>
PreMatrix<T> PreMatrix<T>::sscp( void ) const
{
	assert( !operator!( ) );

	if ( operator!( ) )
		return PreMatrix<T>( 0 );

	PreMatrix<T> result( cols( ), cols( ) );
	assert( result );
	if ( !result )
		return result;

#if defined( SAVE_ALGORITHMS )

	register unsigned i, j, k;
	for ( i = 0; i < cols( ); i++ )
	{
		for ( j = i; j < cols( ); j++ )
		{
			T sum = T( 0 );
			for ( k = 0; k < rows( ); k++ )
				sum += get( k, i ) * get( k, j );
			result( i, j ) = sum;
			if ( i != j )
				result( j, i ) = sum;
		}
	}

#else

	unsigned n = cols( );
	register unsigned i, j, k;
	for ( i = 0; i < n; i++ )
	{
		for ( j = i; j < n; j++ )
		{
			T sum = T( 0 );
			T *workA = getV( ) +  i;
			T *workB = getV( ) +  j;
			for ( k = 0; k < rows( ); k++, workA += n, workB += n )
				sum += (*workA) * (*workB);
			result( i, j ) = sum;
			if ( i != j )
				result( j, i ) = sum;
		}
	}

#endif

	return result;
}

template <class T>
void PreMatrix<T>::prettyPrint( ostream &out)
{
	register unsigned int i, j;

	int *w = new int[ cols( ) ];
	for ( j = 0; j < cols( ); j++ )
	{
		w[ j ] = 0;
		for ( i = 0; i < rows( ); i++ )
		{
			char buffer[128];
			ostrstream item(buffer, sizeof buffer);

			item << get( i, j ) << ends;

			int currW = strlen( item.str( ) );
			if ( currW > w[ j ] )
				w[ j ] = currW;
		}
	}
	for ( i = 0; i < rows( ); i++ )
	{
		for ( j = 0; j < cols( ); j++ )
		{
			out.width( w[ j ] + 2 );
			out << get( i, j );
		}
		out << endl;
	}
	delete [] w;
}

//	int prettyScan(istream &in)
//
// Eine Matrix aus einem formatierten Eingabestrom auslesen.
//
//	Parameter:
//	in - Eingabestrom
//
//	Ergebnis:
// War die Operation erfolgreich?
//
//	Bemerkung:
//	Die PreMatrix muss zu diesem Zweck nicht initialisiert
// sein, die Dimension wird beim Lesen ermittelt, d.h.
//
// PreMatrix<number> X;
//
// X.prettyScan(cin);
//
// ist zulaessig und richtig. 


template <class T>
int 
PreMatrix<T>::
prettyScan( istream &in )
{
	in >> ws;
	
	// Die gelesenen Daten werden zeilenweise, innerhalb der
	// gelesenen Zeilen spaltenweise abgekellert, die so
	// zwischengespeicherte Matrix wird dann Zeilenweise von
	// unten nach oben, spaltenweise von hinten nach vorn in die
	// zu diesem Zweck bereitgestellte PreMatrix eingetragen.

	// Die Eingabe wird zuerst in einem Stack von Zeilen gesammelt

	Stack< Array<T> > inputBuffer;

	//	Zeilenpuffer zum Einlesen einer Datenzeile

	static const unsigned int buffSize = 8192;
	char *lineBuffer = new char[ buffSize ];

	unsigned cols = 0;
	do
	{
		//	Zeile einlesen; abbrechen, wenn das Ende der Eingabedatei
		//	erreicht ist.

		in.getline( lineBuffer, buffSize );
		if ( in.bad( ) || in.eof( ) )
			break;

		//	Die Datenelemente in der Zeile werden ebenfalls erst in
		//	einem Stack von Elementen (dataRow) gesammelt  

		istrstream line( lineBuffer, buffSize );

		Stack<T> dataRow;

		int done = 0;
		do
		{
			//	Datenelemente einzeln lesen

			T curr;

			//	Weisse Leerzeichen ueberspringen

			line >> ws;
			int c = line.peek( );

			//	Gelesene Elemente in die Liste dataRow einfuegen;
			// abbrechen am Ende der Zeile

			if ( c && c != EOF )
			{
				line >> curr;
				if ( !line.bad( ) && !line.eof( ) )
					dataRow.insert( curr );
			}
			else
			{
				done = 1;
			}
		} while ( !done && !line.eof( ) && !line.bad( ) );

		//	Wenn noch keine Zeile mit Mehr als 0 Spalten gelesen wurde,
		//	legt die erste Zeile mit mehr als 0 Spalten die Spaltendimension 
		//	der eingelesenen PreMatrix fest. (Die PreMatrix wird als beeendet
		//	angesehen, wenn eine Zeile mit einer von dieser Zahl abweichenden 
		//	Zahl von Spalten gelesen wird.)

		if ( !cols )
			cols = dataRow.len( );
		else if ( cols != dataRow.len( ) )
			break;

		if ( cols )
		{
			Array<T> row( cols );
			unsigned at = cols - 1;

			while( !dataRow.empty() )
			{
				row( at ) = dataRow.top();
				dataRow.remove();
				--at;
			}
			inputBuffer.insert( row );
		}
	} while( cols && !in.eof( ) && !in.bad( ) );

	//	Gelesene Zeilen in eine PreMatrix eintragen

	unsigned rows = inputBuffer.len( );
	if ( rows && cols )
	{
		discard( );
		m_rows = rows;
		m_cols = cols;
		create( );

		unsigned at = rows - 1;

		while(!inputBuffer.empty())
		{
			const Array<T> &curr = inputBuffer.top();
			register unsigned j;
			for ( j = 0; j < cols; j++ )
				put( at, j, curr( j ) );
			inputBuffer.remove();
			--at;
		}
		return 1;
	}
	return 0;
}


// Ausgabe mit einem bestimmten Delimiter
template <class T>
void PreMatrix<T>::print( ostream& out, char delimiter ) const
{  
  for ( register unsigned i=0; i<rows(); i++ ) 
  {
    for ( register unsigned j=0; j<cols()-1; j++ )
      out << get(i,j) << delimiter;
    out << get(i, cols()-1) << endl;	
  }
}

template <class T>
void PreMatrix<T>::print( ostream& out, char* delimiter ) const
{  
  for ( register unsigned i=0; i<rows(); i++ ) 
  {
    for ( register unsigned j=0; j<cols()-1; j++ )
      out << get(i,j) << delimiter;
    out << get(i, cols()-1) << endl;	
  }
}

template <class T>
bool
PreMatrix<T>::
zero(const T epsilon) const
{
	register unsigned int i, j;
	
	assert(!operator!());
   assert(epsilon >= T(0));
	
	for (i = 0; i < rows(); ++i)
		for (j = 0; j < cols(); ++j)
		   {
			T curr = get(i, j);
			if (curr > epsilon || curr < -epsilon)
				return false;
		   }
   return true;
}

template <class T>
bool 
PreMatrix<T>::
symmetric(const T epsilon) const
{
  register unsigned i, j;

  assert(!operator!());
  assert(rows() == cols());
  assert(epsilon >= T(0));
  
  for (i = 0; i < rows(); i++)
     {
     for (j = i + 1; j < rows(); j++)
        {
	T diff = get(i, j) - get(j, i);
	if (diff < T(0))
	   diff = -diff;
	if (epsilon < diff)
	   return false;
	}
     }
  return true;
}

template 
<class T>
PreMatrix<T> 
PreMatrix<T>::
diag() const
{
	assert(!operator!());
	assert(rows() == cols());

	PreMatrix<T> res(rows(), 1);

	for (register unsigned int i = 0; i < rows(); ++i)
		res.put(i, 0, get(i, i));
	
	return res;
}

template <class T>	
PreMatrix<T> 
PreMatrix<T>::
diag(unsigned int dim, const T v)
{
	assert(dim > 0);
	
	PreMatrix<T> res(dim, dim, T(0));
		
	if (res)
	{
		for (register unsigned i = 0; i < dim; ++i)
			res.put(i, i, v);
	}
	return res;
}


template <class T>
PreMatrix<T> 
PreMatrix<T>::
diag(const PreMatrix<T> &v)
{
	assert(!(!v));
	assert(v.rows() == 1 || v.cols() == 1);

	unsigned int dim = v.rows() > v.cols() ? v.rows() : v.cols();

	PreMatrix<T> res(dim, dim, T(0));
	if (v.cols() == v.rows())
	{
		T value = v.get(0,0);
		for (register unsigned int i = 0; i < dim; ++i)
			res.put(i, i, value);
	}
	else if (v.cols() == 1)
	{
		assert(v.rows() == dim);
		for (register unsigned int i = 0; i < v.rows(); ++i)
			res.put(i, i, v.get(i, 0));
	}
	else if (v.rows() == 1)
	{
		assert(v.cols() == dim);
		for (register unsigned int j = 0; j < v.cols(); ++j)
			res.put(j, j, v.get(0, j));
	}
	return res;
}


template <class T>
PreMatrix<T> 
PreMatrix<T>::
tridiag(const PreMatrix<T> &m, const PreMatrix<T> &lm, const PreMatrix<T> &um)
{
	assert(!m.operator!());
	assert(!lm.operator!());
	assert(!um.operator!());

	assert(m.cols() == 1);
	assert(lm.cols() == 1);
	assert(um.cols() == 1);

	assert(m.rows() >= 2);
	assert(lm.rows() == m.rows() - 1);
	assert(um.rows() == m.rows() - 1);

	PreMatrix<T> res(m.rows(), m.rows(), T(0));

	for (register unsigned int i = 0; i < m.rows(); ++i)
	{
		res.put(i, i, m.get(i, 0));
		if (i < m.rows() - 1)
		{
			res.put(i, i + 1, um.get(i, 0));
			res.put(i + 1, i, lm.get(i, 0));
		}
	}

	return res;
}

// operator&
//
// Untereinanderhaengen in Operatorschreibweise

template <class T>
PreMatrix<T>
PreMatrix<T>:: 
operator&(const PreMatrix<T> &m) const
{
  PreMatrix<T> res;

  Array2D<T>::vcat(m).purge(res);
  return res;
}

// operator|
//
// Nebeneinanderhaengen in Operatorschreibweise
	
template <class T>
PreMatrix<T>
PreMatrix<T>:: 
operator|(const PreMatrix<T> &m) const
{
  PreMatrix<T> res;

  Array2D<T>::hcat(m).purge(res);
  return res;
}

//
//
//
//
// Parameter:
//
// Ergebnis:
//
//

template <class T>
const PreMatrix<T> &
PreMatrix<T>:: 
operator&=(const PreMatrix<T> &m)
{
  PreMatrix<T> res;

  Array2D<T>::vcat(m).purge(res);
  return operator=(res);
}
	
//
//
//
//
// Parameter:
//
// Ergebnis:
//
//

template <class T>
const PreMatrix<T> &
PreMatrix<T>:: 
operator|=(const PreMatrix<T> &m)
{
  PreMatrix<T> res;

  Array2D<T>::hcat(m).purge(res);
  return operator=(res);
}

//
//
//
//
// Parameter:
//
// Ergebnis:
//
//

template <class T>
const PreMatrix<T> &
PreMatrix<T>::
operator=(const PreMatrix<T> &from)
{
  Array2D<T>::operator=(from);
  return *this;
}

//
//
//
//
// Parameter:
//
// Ergebnis:
//
//

template <class T>
PreMatrix<T> 
PreMatrix<T>:: 
blockdiag(const PreMatrix<T> &m)
{
  PreMatrix<T> res(rows() + m.rows(), cols() + m.cols(), T(0));

  res.putBlock(*this, 0, 0);
  res.putBlock(m, rows(), cols());
  return res;
}

//
//
//
//
// Parameter:
//
// Ergebnis:
//
//

template <class T>
PreMatrix<T> 
PreMatrix<T>::
decompLU( int *Index, int *IsEven, int unique )  const
{
  int i, imax, j, k, n;
  T Big, Summe, Temp;

  assert( !(operator!( ) ) );
  assert( rows( ) == cols( ) );

  n = rows( );
  PreMatrix<T> Scalings( n, 1, T( 1 ) );
  assert( Scalings );
  if ( !Scalings )
    return PreMatrix<T>(0);

  if ( IsEven )
    *IsEven = 1;

  PreMatrix<T> Dest(*this);
  assert(Dest);
  if (!Dest)
    return PreMatrix<T>(0);

  for ( i = 0; i < n; i++ )
    {
    Big = T( 0 );

    for ( j = 0; j < n; j++ )
      {
      Temp = Dest( i, j );
      if ( Temp < T( 0 ) )
	Temp = -Temp;
      if (Big < Temp)
	Big = Temp;
      }

    T invBig =  T( 1 ) / Big;
    if ( Big == T( 0 ) || invBig == T( 0 ) )
      return PreMatrix<T>(0);

    Scalings( i, 0 ) = invBig;
    }

  for ( j = 0; j < n; j++ )
    {
    for ( i = 0; i < j; i++ )
      {
      Summe = Dest( i, j );
      for ( k = 0; k < i; k++ )
	Summe -= Dest( i, k ) * Dest( k, j );
      Dest( i, j ) = Summe;
      }

    Big = T( 0 );
    imax = -1;
    for ( i = j; i < n; i++ )
      {
      Summe = Dest( i, j );
      for ( k = 0; k < j; k++ )
	Summe -= Dest( i, k ) * Dest( k, j );
      Dest( i, j ) = Summe;
      if ( Summe < T( 0 ) )
	Summe = -Summe;

      if ( Big < ( Temp = Scalings( i, 0 ) * Summe ) )
	{
	Big  = Temp;
	imax = i;
	}
      }

    assert( imax != -1 );

    if ( j != imax )
      {
      for( k = 0; k < n; k++ )
	{
	Temp = Dest( imax, k );
	Dest( imax, k ) = Dest( j, k );
	Dest( j, k ) = Temp;
	}

      if ( IsEven )
	*IsEven = - (*IsEven);
      Scalings( imax, 0 ) = Scalings( j, 0 );
      }

    if (Index)
      Index[j] = imax;

    if (Dest(j, j) == T(0))
      {
      if (unique)
	return PreMatrix<T>(0);
      else
	Dest(j, j) = T(1) / T(32000);
      }
    if ( j != n - 1 )              
      {
      Temp = T( 1 ) / Dest( j, j );
      for ( i = j + 1; i < n; i++ )
			Dest( i, j ) *= Temp;
      }
    }
  return Dest;
}

//
//
//
//
// Parameter:
//
// Ergebnis:
//
//

template <class T>
PreMatrix<T> 
PreMatrix<T>::
backsubstLU(const PreMatrix<T> &LU, const PreMatrix<T> &bIn, int *Index)
{
  if ( !LU || !bIn )
    return PreMatrix<T>( 0 );

  int i, ii = -1, Pivot, j, n;
  T Summe;

  n = LU.rows( );

  assert( LU.rows( ) == LU.cols( ) );

  PreMatrix<T> b = bIn;
  if ( !b )
    return PreMatrix<T>( 0 );

  for ( i = 0; i < n; i++ )
    {
    Pivot = Index[ i ];
    Summe = b( Pivot, 0 );
    b( Pivot, 0 ) = b( i, 0 );

    if( ii >= 0 )
      {
      for ( j = ii; j < i; j++ )
	Summe -= LU.get( i, j ) * b( j, 0 );
      }
    else if (!(Summe == T(0)))
      ii = i;

    b( i, 0 ) = Summe;
    }

  for( i = n - 1; i >= 0; i-- )
    {
    Summe = b( i, 0 );
    for( j = i + 1; j < n; j++ )
      Summe -= LU.get( i, j ) * b( j , 0 );
    b( i, 0 ) = Summe / LU.get( i, i );
    }

  return b;
}

//
//
//
//
// Parameter:
//
// Ergebnis:
//
//

template <class T>
PreMatrix<T> 
PreMatrix<T>::
luinverse() const
{
  assert( !( operator!( ) ) );
  assert( rows( ) == cols( ) );
  assert( rows( ) );

  if ( rows( ) == 1 )
    {
    T v = get( 0, 0 );

    if ( v == T( 0 ) )
      return PreMatrix<T>( 0 );

    return PreMatrix<T>( 1, 1, T( T(1) / v ) );
    }

  PreMatrix<T> Inverse( rows( ), cols( ) );
  if ( !Inverse )
    return PreMatrix<T>( 0 );

  int *index = new int[ rows( ) ];
  if ( !index )
    return PreMatrix<T>( 0 );

  PreMatrix<T> LU = decompLU( index );
  if ( !LU )
    {
    delete [] index;
    return PreMatrix<T>( 0 );
    }
  for ( register unsigned j = 0; j < cols( ); j++ )
    {
    PreMatrix<T> select( rows( ), 1, T( 0 ) );
    if ( !select )
      {
      delete [] index;
      return PreMatrix<T>(0);
      }

    select( j, 0 ) = T( 1 );

    PreMatrix<T> colInverse = backsubstLU( LU, select, index );
    if ( !colInverse )
      {
      delete [] index;
      return PreMatrix<T>( 0 );
      }

    for ( register unsigned i = 0; i < rows( ); i++ )
      Inverse( i, j ) = colInverse( i, 0 );
    }
  delete [] index;
  return Inverse;
}

//
//
//
//
// Parameter:
//
// Ergebnis:
//
//

template <class T>
PreMatrix<T> 
PreMatrix<T>::
solve ( const PreMatrix<T> &bIn ) const
{
  assert( !( operator!( ) ) );
  assert( rows( ) == cols( ) );
  assert( bIn.rows() == rows() && bIn.cols() == 1 );

  if ( rows( ) == 1 )
    {
    if (get( 0, 0 ) == T(0))
      return PreMatrix<T>( 0 );

    return PreMatrix<T>( 1, 1, T( bIn.get(0,0) / get( 0, 0 ) ) );
    }

  int *index = new int[ rows( ) ];
  if ( !index )
    return PreMatrix<T>( 0 );

  PreMatrix<T> LU = decompLU( index );
  if ( !LU )
    {
    delete [] index;
    return PreMatrix<T>( 0 );
    }
  PreMatrix<T> Solution = backsubstLU( LU, bIn, index );
  if ( !Solution )
    {
    delete [] index;
    return PreMatrix<T>( 0 );
    }
  delete [] index;
  return Solution;
}

//
//
//
//
// Parameter:
//
// Ergebnis:
//
//

template <class T>
PreMatrix<T>
PreMatrix<T>::
kronecker(const PreMatrix<T> &m) const
{
   PreMatrix<T> res(rows() * m.rows(), cols() * m.cols());
   unsigned int i, j;

   for (i = 0; i < rows(); ++i)
      for (j = 0; j < cols(); ++j)
         res.putBlock(get(i, j) * m, i * m.rows(), j * m.cols());
   return res;
}

// PreMatrix vec()
//
// Matrixelelemente in einen Vektor zusammenfassen
//
// Ergebnis:
//
// Spaltenvektor, der die Elemente der spaltenweise durchlaufenen Matrix enth�lt.

template <class T>
PreMatrix<T>
PreMatrix<T>::
vec() const
{
   PreMatrix<T> res(rows() * cols(), 1);
   unsigned int j;

   for (j = 0; j < cols(); ++j)
      res.putBlock(getCol(j), j * rows(), 0);
   return res;
}

// PreMatrix vech()
//
// Horizontal durch eine Matrix gehen und einen Vektor aus den Elementen bilden
//
// Ergebnis:
//
// Zeilenvektor, der die Elemente der zeilenweise durchlaufenen Matrix enth�lt.

template <class T>
PreMatrix<T>
PreMatrix<T>::
vech() const
{
   PreMatrix<T> res(1, cols() * rows());
   unsigned int i;

   for (i = 0; i < rows(); ++i)
      res.putBlock(getRow(i), 0, i * cols());
   return res;
}

// T det()
//
// Determinante berechnen
//
// Ergebnis:
//
// Determinante einer Matrix

template <class T>
T
PreMatrix<T>::
det() const
{
	assert(!operator!());
	assert(rows() == cols());

        int isEven;
	PreMatrix<T> LU = decompLU(0, &isEven);
	if (!LU)
		return PreMatrix<T>(0);

	T res = T(1);
	for (unsigned int i = 0; i < rows(); ++i)
		res *= LU.get(i, i);
	return res * T(isEven);
}


// T trace()
//
// Berechnen der Spur
//
// Ergebnis:
//
// Spur einer quadratischen Matrix

template <class T>
T
PreMatrix<T>::trace() const
{
	assert(!operator!());
	assert(rows() == cols());

	T res = 0.0;
	for (unsigned int i = 0; i < rows(); ++i)
		res += get(i, i);
	return res;
}


	




