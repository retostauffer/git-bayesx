#if !defined(INCL_ARITHM_H)
#define INCL_ARITHM_H

class arithm
{
public:
   arithm() : m_v(0.0) {};
   arithm(int v) : m_v(v) {}
   arithm(const arithm &from) { m_v = from.m_v; }

   arithm operator+(const arithm &ob) const
      { return arithm(m_v + ob.m_v); }
   arithm operator-(const arithm &ob) const
      { return arithm(m_v - ob.m_v); }
   arithm operator*(const arithm &ob) const
      { return arithm(m_v * ob.m_v); }
   arithm operator/(const arithm &ob) const
      { return arithm(m_v / ob.m_v); }

   arithm operator-()
     { return arithm(-m_v); }

   int operator==(const arithm &other) const
      { return m_v == other.m_v; } 
   int operator<(const arithm &other)
     { return m_v < other.m_v; } 
   int operator<=(const arithm &other)
     { return m_v <= other.m_v; } 
   arithm operator=(const arithm &from)
      { m_v = from.m_v; return *this; }
   arithm operator+=(const arithm &from)
      { m_v += from.m_v; return *this; }
   arithm operator-=(const arithm &from)
      { m_v -= from.m_v; return *this; }
   arithm operator*=(const arithm &from)
      { m_v *= from.m_v; return *this; }
   friend ostream &operator << (ostream &os, const arithm &ob)
      { os.precision(3); os << ob.m_v; return os; } 
   friend istream &operator >> (istream &is, arithm &ob)
      { is >> ob.m_v; return is; } 

 private:
   arithm(double v) : m_v(v) {}
   double m_v;
};

#endif
