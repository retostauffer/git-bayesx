// DATE: 14.11.97


#if !defined (ERRORM_INCLUDED)

#define ERRORM_INCLUDED


#include<iostream.h>
#include<clstring.h>
#include<vector>


namespace errorm
{

	class messages : public vector<ST::string>
  {


  public:


  //--------------------------- PUBLIC FUNCTIONS -------------------------------

  // DEFAULT CONSTRUCTOR

	  messages(void) : vector<ST::string>() {}

  // COPY CONSTRUCTOR

	  messages(const messages & e) : vector<ST::string>(vector<ST::string>(e)) {}

  // OVERLOADED ASSIGNMENT OPERATORS

  const messages & operator=(const messages & e);

  // FUNKTION: insert_back

  void insert_back(const messages & e)
	 {
	 if (! e.empty())
		insert(end(),e.begin(),e.end());
	 }

  // OVERLOADED << OPERATOR

  friend ostream & operator<<(ostream & c, const messages & em);

  // FUNCTIONS: clear
  // TASK: deletes all errormessages

  void clear(void)
	 {
	 if (! empty())
		erase(begin(),end());
	 }


  };

}

#endif
