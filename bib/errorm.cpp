// DATE: 14.11.97


#include<errorm.h>

namespace errorm
{

const messages & messages::operator=(const messages & e)
  {
  if (this == &e)
	 return *this;
  vector<ST::string>::operator=(vector<ST::string>(e));
  return *this;
  }


ostream & operator<<(ostream & c, const messages & em)
  {
  if (em.empty())
	 return c;
  else
	 {
	 int i;
	 for (i=0;i<em.size();i++)
		c << em[i];
	 return c;
	 }
  }

}

