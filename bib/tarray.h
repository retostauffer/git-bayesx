//	/home/kurt/src/templib/SCCS/s.tarray.h 3.1 97/08/05 01:55:49
//
//	Deklaration der template-Klasse Array

#if !defined(TARRAY_H_INCLUDED)

#define TARRAY_H_INCLUDED

#include <iostream.h>
#include <assert.h>

//	Array - Feld fester Groesse

template <class T> class Array
{
   public :

	//	1. Konstruktoren und Destruktor

	//	Default-Konstruktor

	Array()
		{ m_size = 1; create(); }

	//	Konstruktor fuer uninitalisiertes Array

	Array(unsigned int size)
		{ m_size = size; create(); }

	//	Konstruktor fuer mit festem Wert initialisiertes Array

	Array(unsigned int size, const T init);

	//	Konstruktor fuer mit mehreren Werten initialisiertes Array.
	//	Der Vektor init muss ein Zeiger auf size Objekte der Klasse
	//	T sein.
	//
	//	Dieser Konstruktor erm"oglicht die Umwandlung eines built-in
	//	Arrays in ein Array der Bibliothek.

	Array(unsigned int size, const T *init);

	//	Kopierkonstruktor

	Array(const Array &init);

	//	Destruktor

	virtual ~Array()
		{ discard(); }

	//	2. "Offentliche Elementfunktionen

	//	Groesse - Anzahl der Elemente im Array

	const unsigned int size() const
		{ return m_size; }

	//	Lesender Zugriff auf ein Element des Arrays. Dieser ist
	//	auch dann zugelassen, wenn der Aufrufer nur ein const Array
	//	hat.

	const T &get(unsigned int at) const
		{ assert(at < m_size); return m_v[ at ]; }

	//   schreibender Zugriff auf ein Array-Element

	void put(unsigned int at, const T &v)
		{ assert(at < m_size); m_v[ at ] = v; }

	//	3. Operatoren

	//	Zugriff auf eine Array-Element als lvalue

	T &operator()(unsigned int at)
		{ assert(at < m_size); return m_v[at]; }

	const T &operator()(unsigned int at) const
		{ assert(at < m_size); return m_v[at]; }

	//	Zuweisungsoperator fuer Arrays

	const Array &operator=(const Array &from);

	//	negative Validitaetspruefung

	int operator!() const
		{ return m_v ? 0 : 1; }

	//	 positive Validitaetspruefung

	operator int() const
		{ return m_v ? 1 : 0; }

	//	4. Befreundete Klassen und Funktionen

	friend ostream &operator << (ostream &out, const Array<T> &fld)
		{ fld.writeOn(out); return out; }

	friend istream &operator >> (istream &in, Array<T> &fld)
		{ fld.readFrom(in); return in; }

	int operator==(const Array<T> &other) const;

  private :

	//	5. Instanzvariablen

	//	Daten - Zeiger auf die Elemente des Arrays

	T *m_v;

	//	Groesse - Elemente des Arrays

	unsigned int m_size;

	//	6. Hilfsfunktionen

	//	Klartextausgabe

	void writeOn(ostream &out) const;

	//	Klartexteingabe

	void readFrom(istream &in);

	//	Speicher freigeben

	void discard()
		{ if (m_v) delete [] m_v; m_v = 0; }

	//	Speicher bereitstellen

	void create();


	//	Inhalte kopieren

	void copyContents(const Array &from);
};

#if defined(TEMPL_INCL_DEF)
#	if defined(CC_SOURCE)
#		include <tarray.cc>
#	else
#		include <tarray.cpp>
#	endif
#endif

#endif
