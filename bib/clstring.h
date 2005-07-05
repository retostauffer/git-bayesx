// DATE: 15.12.97

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __declspec(dllimport)
#endif




#if !defined (CLSTRING_INCLUDED)

#define CLSTRING_INCLUDED


#include<string.h>
#include<stdio.h>
#include<iostream.h>
#include<fstream.h>
#include<assert.h>
#include<sys\stat.h>
#include<list.h>
#include<vector.h>




using std::vector;
using std::list;

//------------------------------------------------------------------------------
//----------------------------- CLASS string -----------------------------------
//------------------------------------------------------------------------------

namespace ST
 {

class __EXPORT_TYPE string
  {

  private:


  //------------------------- PRIVATE VARIABLES --------------------------------

  // stores the string



  // length of the string




  //------------------------- PRIVATE FUNCTIONS --------------------------------

  // FUNCTION: checkindex
  // TASK: checks, if access to charakter 'i' (of the calling string) is //       possible

  void checkindex(unsigned i) const
	 {
	 assert (len > 0);
	 assert(i < len);
	 }


  public:

  char * str;

  unsigned len;
  //------------------------- PUBLIC FUNCTIONS ---------------------------------

  // DEFAULT CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - length = 0

  //string();
  string()
    {
    len = 0;
    str = new char[1];
    strcpy(str,"");
    }


  // CONSTRUCTOR
  // TASK: creates a string that consists of 'len' signs 'sign';

  //string(const char & sign,const unsigned & len);
  string(const char & sign,const unsigned & l)
    {
    len = l;
    str = new char[len+1];
    for (unsigned i=0;i<len;i++)
  	 str[i] = sign;
    str[len] = '\0';
    }


  // COPY CONSTRUCTORS

  //string(const char * s);
  string(const char * s)
    {
    len = strlen(s);
    str = new char[len+1];
    strcpy(str,s);
     }

  //string(const string & s);
  string(const string & s)
    {
    len = s.len;
    str = new char[len+1];
    strcpy(str,s.str);
    }

  //string(const std::string & s);
  string(const std::string & s)
    {
    len = s.length();
    str = new char[len+1];
    strcpy(str,s.c_str());
   }

  // DESTRUCTOR

  ~string()
	 {
	 delete [] str;
	 }

    // FUNCTION: strtochar
    // TASK: converts the calling string to char* and returns the result



    //char * strtochar() const;
    char * strtochar() const
	  {
	  char * h = new char[len+1];
	  strcpy(h,str);
	  return h;
      }

   // FUNCTION: to_bstr
   // TASK: converts the string to Std::string an returns the result

   //std::string to_bstr(void) const;
   std::string to_bstr(void) const
     {
     return std::string(str);
     }

  // OVERLOADED ASSIGNMENT OPERATORS

  //const string & operator=(const string & s);
  const string & operator=(const string & s)
    {
    if (this == &s)
  	 return *this;
    delete [] str;
    len = s.len;
    str = new char[len+1];
    strcpy(str,s.str);
    return *this;
  }

  //const string & operator=(string & s);
  const string & operator=(string & s)
    {
    if (this == &s)
  	 return *this;
    delete [] str;
    len = s.len;
    str = new char[len+1];
    strcpy(str,s.str);
    return *this;
  }

  //const string & operator=(const char * s);
  const string & operator=(const char * s)
    {
    delete [] str;
    len = strlen(s);
    str = new char[len+1];
    strcpy(str,s);
    return *this;
  }

  //const string & operator=(const std::string & s);
  const string & operator=(const std::string & s)
    {
    delete[] str;
    len = s.length();
    str = new char[len+1];
    strcpy(str,s.c_str());
    return *this;
  }

  // OVERLOADED + OPERATORS

  friend string __EXPORT_TYPE operator+(const string & st1,const string & st2);
  friend string __EXPORT_TYPE operator+(const char * s,const string & st);

  // OVERLOADED OUTPUT OPERATOR

  friend ostream & __EXPORT_TYPE operator<<(ostream & c, const string & s);


  // OVERLOADED INPUT OPERATOR

  friend istream & __EXPORT_TYPE operator>>(istream & i, string & s);

  // FRIEND FUNCTION: getline
  // TASK: takes signs from the input stream 'i' until the delimeter
  //       is reached or number of signs is 'maxlen'.
  //       the resulting string is stored in 's'

  friend istream & __EXPORT_TYPE getline(istream & i,unsigned int maxlen,
									string & s, char delim = '\n');

  // FRIEND FUNCTION: getline
  // TASK: see function above
  // ADDITIONAL INFORMATION:
  // - maxlen = 256 Bytes

  friend istream & __EXPORT_TYPE getline(
  istream & i, string & s, char delim = '\n');

  // FRIEND FUNCTION: open
  // TASk:

//friend void (std::ifstream & fin,string & s, int mode);

  friend void __EXPORT_TYPE open(
  std::ifstream & fin,string & s,int mode = ios::in);

//MICRO  friend void open(ofstream & out,string & s,ios_base::openmode mode = ios_base::out);
  friend void __EXPORT_TYPE open(
  std::ofstream & out,string & s,int mode = ios::out);


  // OVERLOADED [] OPERATOR
  // TASK: Zugriff auf ein einzelnes Zeichen
  // ADDITIONAL INFORMATION:  if i is invalid, the program terminates and an
  // errormessage occurs

  //char & operator[] (int i);
  char & operator[](int i)
    {
    checkindex(i);
    return str[i];
    }

  // OVERLOADED ASSIGNMENT OPERATORS

  friend int __EXPORT_TYPE operator==(const string & s1, const char * s2);

  friend int __EXPORT_TYPE operator==(const string & s1, const string & s2);

  friend int __EXPORT_TYPE operator!=(const string & s1, char * s2);

  friend int __EXPORT_TYPE operator!=(string & s1, string & s2);


  int operator<(const string & s2) const

  	  {
  	  if (len < s2.len)
  		 return 1;
  	  else if (len > s2.len)
  		 return 0;
  	  else
  		 {
         return strcmp(str,s2.str) < 0;
  		 }
  	  }



  friend int __EXPORT_TYPE operator>(string & s1,string & s2);



  // FUNCTION: length
  // TASK: returns the length of the string

  int length(void) const
	 {
	 return len;
	 }

  // Liefert die L‰nge des Teilsstrings der ausschlieﬂlich aus Zeichen besteht,
  // die in s enthalten sind

  int spn(const string & s) const
	 {
	 return strspn(str,s.str);
	 }

  int spn(char * s) const
	 {
	 return strspn(str,s);
	 }

	// FUNCTION: firstpos
	// TASK: returns the position in the calling string of the first appearence
	//       of 'sign', else -1

	//int firstpos (char sign) const;
	int firstpos (char sign) const
	  {
	  int pos = -1;
	  unsigned i = 0;
	  while ( (i < len) && (pos == -1) )
		 {
		 if (str[i] == sign)
			pos = i;
		 i++;
		 }
	  return pos;
  }

  // FUNCTION: helpfill
  // TASK: returns the string + the number of whitespaces that
  //       are necessary to complete the wanted width of the output column

  //string string::helpfill(unsigned n);
  string string::helpfill(unsigned n)
    {
    unsigned d;
    ST::string wert;
    if (n>=len)
      {
      d = n - len;
      wert = substr(0, len);
      }
    else
      {
      wert = substr(0, n-2) + "~";
      d = 1;
      }
    ST::string platz = string(' ',d);
    ST::string out = platz + wert;
    return out;
    }


  // FUNCTION: substr
  // TASK: returns a substring of the calling string
  //       copies 'nr' signs of the string starting at position 'pos'
  // POSSIBLE ERRORS:
  // - program terminates if  pos+nr > length of the string

  //string substr(unsigned pos,unsigned nr) const;
  string substr(unsigned pos, unsigned nr) const
    {
    assert(pos+nr <= len);
    assert(nr > 0);

    char * help = new char[nr+1];
    strncpy(help,str+pos,nr);
    help[nr] = '\0';
    string ret(help);
    delete [] help;
    return ret;
  }

  // FUNCTION: deletesign
  // TASK: deletetes the sign at position 'pos' of the calling string
  //       and returns the result
  // ADDITIONAL INFORMATION:
  // asserts if pos is a valid index (leads to an assertion failure,
  // abnormal program termination)

  //string deletesign(unsigned pos) const;
  string deletesign(unsigned pos) const
    {
    checkindex(pos);
    if (pos == 0)
      {
      if (len > 1)
        return substr(1,len-1);
      else
        return "";
      }
    else if (pos == len-1)
  	 return substr(0,len-1);
    else
      return substr(0,pos) + substr(pos+1,len-pos-1);
  }

  // FUNCTION: deleteallsigns
  // TASK: deletes all 'sign' signs of the calling string and returns the
  //       result
  // EXAMPLE:
  // string test = "today"
  // test.deleteallsigns('d') returns the string "toay"

  //string deleteallsigns(char sign) const;
  string deleteallsigns(char sign) const
    {
    string result = *this;
    int i=0;
    while ( i < result.length() )
  	 {
  	 if (result[i] == sign)
  		{
  		result = result.deletesign(i);
  		}
  	 else
  		i++;
  	 }
    return result;
  }

  // FUNCTION: replaceallsigns
  // TASK: replaces all 'oldsign' signs of the calling string with 'newsign' signs
  //       and returns the result
  // EXAMPLE:
  // string test = "today"
  // test.replaceallsigns('d','k') returns the string "tokay"

  //string replaceallsigns(char oldsign, char newsign) const;
  string replaceallsigns(char oldsign, char newsign) const
    {
    int i=0;
    string result = *this;
    while (i< result.length())
  	 {
  	 if (result[i] == oldsign)
  		result[i] = newsign;
  	 i++;
  	 }
    return result;
  }

  //FUNCTION: insert_string_num
  //TASK: inserts the string 'str' at the position 'pos' into a string
  //EXAMPLE:
  // string test = "today"
  // test.insert_string(2, "hi") returns the string "tohiday"

  //string insert_string_num(unsigned pos, string & str) const;
  string insert_string_num(unsigned pos, string & str) const
    {
    string s = *this;
    assert(pos<s.length());
    string s1 = s.substr(0, pos);
    string s2 = s.substr(pos, s.length()-pos);
    string result = s1 + str + s2;
    return result;
  }

  //FUNCTION: insert_string_char
  //TASK: inserts the string 'str' instead of each character 'p' into a string
  //EXAMPLE:
  // string test = "today"
  // test.insert_string('o', "hi") returns the string "thiday"

  //string insert_string_char(char p, string & str) const;
  string insert_string_char(char p, const string & str) const
    {
    string s = *this;
    string result = s;
    unsigned l = str.length();
    unsigned k = 0;
    for (unsigned i=0;i<s.length()-1;i++)
        {
        char z = s[i];
        if (z == p)
           {
           string s1 = result.substr(0, i+k*l-k);
           string s2 = result.substr(i+1+k*l-k, result.length()-(i+1+k*l-k));
           result = s1 + str + s2;
           k = k + 1;
           }
        }
    return result;
  }

  //FUNCTION: insert_after_string
  //TASK: inserts the string 's1' after the first appearance of the string 's2'
  //NOTE: returns the calling string if there is no appearance of 's2'
  //EXAMPLE:
  // string test = "today"
  // test.insert_after_string("hi", "od") returns the string "todhiay"

  //string insert_after_string(string s1, string s2) const;
  string insert_after_string(string s1, string s2) const
    {
    string s = *this;
    string shelp;
    unsigned n2 = s2.length();
    unsigned n = s.length()-n2;
    unsigned i;
    for(i=0; i<n+1; i++)
      {
      shelp=s.substr(i,n2);
      if(shelp==s2)
        {
        shelp = s.substr(0,i+n2)+s1;
        if(i<n-1)
          {
          shelp = shelp+s.substr(i+n2,n-i);
          }
        return shelp;
        }
      }
    return s;
  }

  //FUNCTION: insert_after_all_string
  //TASK: inserts the string 's1' after every appearance of the string 's2'
  //NOTE: returns the calling string if there is no appearance of 's2'
  //EXAMPLE:
  // string test = "today"
  // test.insert_after_string("hi", "od") returns the string "todhiay"

  //string insert_after_all_string(string s1, string s2) const;
  string insert_after_all_string(string s1, string s2) const
    {
    string s = *this;
    string shelp;
    string result = " ";

    unsigned n2 = s2.length();
    unsigned n = s.length()-n2;
    unsigned i,k;
    unsigned j=0;

    for(i=0; i<n+1; i++)
      {
      shelp=s.substr(i,n2);
      if(shelp==s2)
        {
        if(j==0)
          {
          result=s.substr(0,i+n2)+s1;
          j++;
          }
        else
          {
          result = result+s.substr(k+n2,i-k)+s1;
          }
        k=i;
        }
      }
    if(j==1)
      {
      if(k<n-1)
        {
        result = result + s.substr(k+n2,n-k);
        }
      return result;
      }
    return s;
    }


  // FUNCTION: eatwhitespace
  // TASK: returns a string, where (possible) leading and ... whitespace signs
  //       of the calling string are deleted
  // EXAMPLE:
  // string test = "  today ";
  // string result = test.eatwhitespace();
  // result evaluates to "today"

  //string eatwhitespace(void) const;
  string eatwhitespace(void) const
    {

    int beginpos;
    int endpos;

    int i=0;

    while ( (i<len) && (str[i] == ' ') )
  	 i++;
    if (i == len)
  	 return "";
    beginpos = i;

    i=len-1;
    while ( (i >= 0) && (str[i] == ' ') )
  	 i--;
    endpos=i;

    return substr(beginpos,endpos-beginpos+1);

  }

  // FUNCTION: eatallwhitespace
  // TASK: deletes all whitespacesigns of the calling string and returns the
  //       result
  // EXAMPLE:
  // string test = " hello world"
  // test.eatallwhitespace() returns the string "helloworld"

  string eatallwhitespace(void) const
	 {
	 return deleteallsigns(' ');
	 }

  // FUNCTION: closingbracketpos
  // TASK: returns the position of the closing bracket ) in the calling string,
  //       if 'bracketpos' is the position of an opening bracket
  //       returns -1 if the closing bracket is missing
  // ADDITIONAL INFORMATION:
  // - asserts 0 <= bracketpos < length  of the string
  // - asserts if position 'bracketpos' contains an opening bracket '('
  // EXAMPLE:
  // string test = "sin(x+y)+4"
  // test.closingbracketpos(3) returns 7
  // string test = "sin(4+x"
  // test.closingbracketpos(3) returns -1 (closing bracket is missing)

  //int closingbracketpos(const unsigned bracketpos) const;
  int closingbracketpos(const unsigned bracketpos) const
    {
    assert (bracketpos < len);
    assert (str[bracketpos] == '(');

    unsigned i=bracketpos+1;
    unsigned nropen = 1;
    while ( (i < len) && (nropen > 0) )
  	 {
  	 if (str[i] == '(')
  		nropen++;
  	 else if (str[i] == ')')
  		nropen--;
  	 i++;
  	 }
    if (nropen == 0)
  	 return i-1;
    else
  	 return -1;
  }

  // FUNCTION: closingbracketpos
  // TASK: returns the position of the closing bracket ] in the calling string,
  //       if 'bracketpos' is the position of an opening bracket
  //       returns -1 if the closing bracket is missing
  // ADDITIONAL INFORMATION:
  // - asserts 0 <= bracketpos < length  of the string
  // - asserts if position 'bracketpos' contains an opening bracket '('
  // EXAMPLE:
  // string test = "sin[x+y]+4"
  // test.closingbracketpos(3) returns 7
  // string test = "sin[4+x"
  // test.closingbracketpos(3) returns -1 (closing bracket is missing)

  //int closingbracketpos2(const unsigned bracketpos) const;
  int closingbracketpos2(const unsigned bracketpos) const
    {
    assert (bracketpos < len);
    assert (str[bracketpos] == '[');

    unsigned i=bracketpos+1;
    unsigned nropen = 1;
    while ( (i < len) && (nropen > 0) )
  	 {
  	 if (str[i] == '[')
  		nropen++;
  	 else if (str[i] == ']')
  		nropen--;
  	 i++;
  	 }
    if (nropen == 0)
  	 return i-1;
    else
  	 return -1;
  }

  // FUNCTION: lowestprecedencepos
  // TASK:

  //int lowestprecedencepos(string & sign) const;
  int lowestprecedencepos(string & sign) const
    {
    int i = 0;
    int pos = -1;
    int minp = 6;
    while (i < len)
  	 {
  	 if (str[i] == '(')
  		{
  		i = closingbracketpos(i);
  		if (i == -1)
  		  return -2;
  		i++;
  		}
       else if (str[i] == '[')
         {

  		i = closingbracketpos2(i);
  		if (i == -1)
  		  return -2;
  		i++;

         }
  	 else if ( (str[i] == '+') || (str[i] == '-') )
  		{
  		if (minp >= 3)
  		  {
  		  minp = 3;
  		  pos = i;
  		  if (str[i] == '+')
  			 sign = "+";
  		  else
  			 sign = "-";
  		  }
  		i++;
  		}
  	 else if ( (str[i] == '*') || (str[i] == '/') )
  		{
  		if (minp >= 4)
  		  {
  		  minp = 4;
  		  pos = i;
  		  if (str[i] == '*')
  			 sign = "*";
  		  else
  			 sign = "/";
  		  }
  		i++;
  		}
  	 else if ( (str[i] == '^') )
  		{
  		if (minp >= 5)
  		  {
  		  minp = 5;
  		  pos = i;
  		  sign = "^";
  		  }
  		i++;
  		}
  	 else if (str[i] == '=')
  		{
  		if (minp >= 2)
  		  {
  		  minp = 2;
  		  pos = i;
  		  sign = "=";
  		  }
  		i++;
  		}
  	 else if ( (str[i] == '>') && (i+1 < len) && (str[i+1] == '=') )
  		{
  		if (minp >= 2)
  		  {
  		  minp = 2;
  		  pos = i;
  		  sign = ">=";
  		  }
  		i=i+2;
  		}
  	 else if ( (str[i] == '<') && (i+1 < len) && (str[i+1] == '=') )
  		{
  		if (minp >= 2)
  		  {
  		  minp = 2;
  		  pos = i;
  		  sign = "<=";
  		  }
  		i=i+2;
  		}
  	 else if ( (str[i] == '!') && (i+1 < len) && (str[i+1] == '=') )
  		{
  		if (minp >= 2)
  		  {
  		  minp = 2;
  		  pos = i;
  		  sign = "!=";
  		  }
  		i=i+2;
  		}
  	 else if (str[i] == '>')
  		{
  		if (minp >= 2)
  		  {
  		  minp = 2;
  		  pos = i;
  		  sign = ">";
  		  }
  		i++;
  		}
  	 else if (str[i] == '<')
  		{
  		if (minp >= 2)
  		  {
  		  minp = 2;
  		  pos = i;
  		  sign = "<";
  		  }
  		i++;
  		}
  	 else if (str[i] == '&')
  		{
  		if (minp >= 1)
  		  {
  		  minp = 1;
  		  pos = i;
  		  sign = "&";
  		  }
  		i++;
  		}
  	 else if (str[i] == '|')
  		{
  		if (minp >= 1)
  		  {
  		  minp = 1;
  		  pos = i;
  		  sign = "|";
  		  }
  		i++;
  		}
  	 else
  		i++;
  	 }
    return pos;
  }


  // FUNCTION: isfunction
  // TASK: returns 1, if the calling string is a function, i.e. of type
  //       functionname(argument)
  //       in 'functionname' the name of the function will be stored,
  //       in 'argument' the argument of the function will be stored
  //       returns 0, if calling string is not a function
  //        returns -1 if closing bracket is missing
  // ADDITIONAL INFORMATION:
  // - whitespace signs should be deleted, before calling the function
  // EXAMPLE:
  // sin(x+z-3).isfunction(functionname,argument) will return 1 and
  // functionname = "sin" and argument = "x+z-3"

  //int isfunction(string & functionname,string & argument) const;
  int isfunction(string & functionname,string & argument) const
    {
    int startbr = firstpos('(');
    if (startbr > 0)
      {
      int endbr = closingbracketpos(startbr);
      if (endbr == -1)
        return -1;
      if (endbr == len-1)
        {
        functionname = substr(0,startbr);
        if (len -startbr-2 > 0)
          argument = substr(startbr+1,len-startbr-2);
        else
          argument = "";
        return 1;
        }
      else
        return 0;
      }
    else
  	 return 0;
  }

  //int issubscribing(string & varname, string & argument) const;
  int issubscribing(string & varname,string & argument) const
    {
    int startbr = firstpos('[');
    if ( (startbr > 0) && (startbr != -1) )
  	 {
  	 int endbr = closingbracketpos2(startbr);
  	 if (endbr == len-1)
  		{
  		varname = substr(0,startbr);
  		argument = substr(startbr+1,len-startbr-2);
  		return 1;
  		}
  	 else
  		return -1;
  	 }
    else
  	 return -1;
  }

  // FUNCTION: isexistingfile
  // TASK: checks, if calling string contains a path to an existing file,
  //       that can be opened for reading
  //       returns 0 if path is valid
  //               1 if path is invalid (error)

  //int isexistingfile(void) const;
  int isexistingfile(void) const
    {
    ifstream fin(str,ios::in);
    if (fin.fail() != 0)
  	 return 1;
    else
  	 return 0;
    }


  // FUNCTION: isvalidfile
  // TASK: checks, if calling string contains a path to a file, that can be
  //       opened for writing
  //       returns: 1 if file can not be opened for writing
  //                0 if file can be opened for writing and is not existing
  //                -1 if file can be opened for writing but is already existing

  //int isvalidfile(void) const;
  int isvalidfile(void) const
    {
    struct stat statbuf;
    int existing = stat(str,&statbuf);
    if (existing == 0)
  	 {
  	 ofstream fout(str,ios::app);
  	 if (fout.fail())
  		return 1;
  	 else
  		return -1;
  	 }
    else
  	 {
  	 ofstream fout(str);
  	 if (fout.fail())
  		{
  		fout.close();
  		remove(str);
  		return 1;
  		}
  	 else
  		{
  		fout.close();
  		remove(str);
  		return 0;
  		}
  	 }
    }



  // FUNCTION: isvarname
  // TASK: checks, if calling string is a  valid variable name
  //       returns: 1 = error, invalid variablename
  //                0 = no error
  // VALID VARIABLE NAME:
  // - first sign is a literal
  // (abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_)
  // - following signs are either literals or numbers (0123456789)

  //int isvarname() const;
  int isvarname() const
    {
    if (len > 0)
  	 {
  	 string valid =
  	 "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789";
  	 unsigned anz = spn(valid);                       // number of valid signs
  	 if (anz < len)
  		return 1;
  	 else
  		{
  		string numbers = "0123456789";
  		string firstsign = substr(0,1);
  		if (firstsign.spn(numbers) != 0)
  		  return 1;
  		else
  		  return 0;
  		}
  	 }
    else
  	 return 1;
  }



  int removefile(void) const
	 {
	 return remove(str);
	 }


  //int isint(void) const;
  int isint(void) const
    {
    if ((len == 0))
  	 return 0;
    else
  	 {
  	 int h = 1;
  	 int i;
  	 if ((str[0] == '-')  || (str[0] == '+'))
  		i = 1;
  	 else
  		i = 0;
  	 string firstsign = "123456789";
  	 string nextsigns = "0123456789";
  	 if (firstsign.checksign(str[i]) == -1)
  		h = 0;
  	 i++;
  	 while ((i < len) && (h == 1))
  		{
  		if (nextsigns.checksign(str[i]) == -1)
  		  h = 0;
  		i++;
  		}
  	 return h;

  	 }
  }

  // FUNCTIONS: strtolong,strtochar,strtouchar,strtodouble
  // TASK: functions, that convert strings into numbers
  //       Converts the calling string into the long value
  //       (char value, unsigned char value, double value) "value"
  //       returns:  0 if no error occured
  //                 1 = if an error occured (value will not be changed)

  //int strtolong(long & value) const;
  int strtolong(long & value) const
    {
    if (len > 0)
  	 {
  	 char * sentinel;
  	 long h = strtol(str, & sentinel,10);
  	 if (sentinel != str+len)
  		return 1;
  	 else
  		{
  		value = h;
  		return 0;
  		}
  	 }
    else
  	 return 1;
  }


  //int strtochar(char & value) const;
  int strtochar(char & value) const
    {
    long h;
    if (strtolong(h) == 1)
  	 return 1;
    else if ((h < CHAR_MIN) || (h > CHAR_MAX))
  	 return 1;
    else
  	 {
  	 value = h;
  	 return 0;
  	 }
  }

  //int strtouchar(unsigned char & value) const;
  int strtouchar(unsigned char & value) const
    {
    long h;
    if (strtolong(h) == 1)
  	 return 1;
    else if ((h < 0) || (h > UCHAR_MAX))
  	 return 1;
    else
  	 {
  	 value = h;
  	 return 0;
  	 }
  }

  //int strtodouble(double & value) const;
  int strtodouble(double & value) const
    {
    if (len > 0)
  	 {
  	 char * sentinel;
  	 double h = strtod(str, & sentinel);
  	 if (sentinel != str+len)
  		return 1;
  	 else
  		{
  		value = h;
  		return 0;
  		}
  	 }
    else
  	 return 1;
  }

  // FUNCTION: checksign
  // TASK: checks, if the sign 'signs' is a member of the calling string
  //       returns the position of the first sign in the calling string,
  //       that is equal to 'signs'
  //       returns -1 if 'sign' is not a member of the calling string

  //int checksign(const char sign) const;
  int checksign(const char sign) const
    {
    unsigned i=0;
    int isfrom = -1;
    while ((i < len) && (isfrom == -1))
  	 {
  	 if (str[i] == sign)
  		isfrom = i;
  	 i++;
  	 }
    return isfrom;
    }


  // FUNCTION: isinlist
  // TASK: returns -1, if calling string is not a meber of 'stringlist'
  //                position of the string in stringlist otherwise

  //int isinlist(const vector<string> & stringlist) const;
  int isinlist(const vector<string> & stringlist) const
    {
    unsigned i = 0;
    int isinl = -1;
    while ((i < stringlist.size()) && (isinl == -1))
  	 {
  	 if (*this == stringlist[i])
  		isinl = i;
  	 i++;
  	 }
    return isinl;
   }

  // FUNCTION: getFirstToken
  // TASK: returns the first token of the calling string
  //       if length = 0 an empty string will be returned
  //       delimeters are stored in 'parsingsigns'

  //string getFirstToken(const string & parsingsigns) const;
  string getFirstToken(const string & parsingsigns) const
    {
    if (len > 0)
  	 {
  	 unsigned i;
  	 unsigned j=0;
  	 while (str[j] == ' ')
  		j++;
  	 i = j;
  	 while ( (i < len) && (parsingsigns.checksign(str[i]) == -1))
  		i++;
       if (i > j)
  	   return substr(j,i-j);
       else
         return string();
  	 }
    else
  	 return string();
  }

  // FUNCTION: strtoken
  // TASK: returns the tokens of the string as a vector of strings
  //       delimeters are stored in 'parsingsigns'

  //vector<string> strtoken(const string & parsingsigns,bool inclsigns = true) const;
  vector<string> strtoken(const string & parsingsigns,bool inclsigns = true) const
    {
    vector<string>  hilfe;


    if (len > 0)
  	 {
  	 int i=0;                             // looping variable
  	 while (i < len)
  		{
  		if (parsingsigns.checksign(str[i]) != -1)   // str[i] is a parsingsign
  		  if (str[i] == ' ')
  			 while ( (i < len) && (str[i] == ' ') )
  				i++;
  		  else
  			 {
               if (inclsigns)
  			   hilfe.push_back(substr(i,1));
  			 i++;
  			 }
  		else                                       // str[i] is not a parsingsign
  		  {
  		  int anf = i;
  		  while ( (i < len) && (parsingsigns.checksign(str[i]) == -1) )
  			 i++;
  		  hilfe.push_back(substr(anf,i-anf));
  		  }
  		}
  	 }
    return hilfe;
  }

  //int strtoken_quot(vector<string> & hilfe,const string & parsingsigns,
  //bool inclsigns = true) const;
  int strtoken_quot(vector<string> & hilfe, const string & parsingsigns,
                                       bool inclsigns = true) const
    {

    bool ok = true;

    if (len > 0)
  	 {
  	 int i=0;                             // looping variable
  	 while (i < len)
  		{
  		if (parsingsigns.checksign(str[i]) != -1)   // str[i] is a parsingsign
  		  if (str[i] == ' ')
  			 while ( (i < len) && (str[i] == ' ') )
  				i++;
  		  else
  			 {
               if (inclsigns)
  			   hilfe.push_back(substr(i,1));
  			 i++;
  			 }
  		else                                       // str[i] is not a parsingsign
  		  {
  		  int anf = i;

            if (str[i] != '"')
              {
  		    while ( (i < len) && (parsingsigns.checksign(str[i]) == -1) )
                i++;
  		    hilfe.push_back(substr(anf,i-anf));
              }
            else
              {
              i++;
              ok = false;
              while ( ( i< len) && (str[i] != '"') )
                i++;

              if ( ( i < len) && (str[i] == '"') )
                ok = true;

              i++;
              if (i-anf-2 > 0)
                hilfe.push_back(substr(anf+1,i-anf-2));
              else
                hilfe.push_back("");
              }

  		  }
  		}
  	 }

    return ok;

  }


  // FUNCTION: strtoken2
  // TASK: returns the tokens of the string as a vector of strings
  //       delimeters are stored in 'parsingsigns'
  // ADDITIONAL INFORMATION:
  // Difference to strtoken: signs in brackets are ignored

  //vector<string> strtoken2(const string & parsingsigns,bool & bracketmiss) const;
  vector<string> strtoken2(const string & parsingsigns,bool & bracketmiss) const
    {
    vector<string>  hilfe;

    bracketmiss = false;

    if (len > 0)
  	 {
  	 int i=0;                             // looping variable
  	 int anf = 0;
  	 while (i < len)
  		{
  		if (str[i] == '(')
  		  {
  		  i = closingbracketpos(i);
  		  if (i == -1)
              {
              bracketmiss = true;
              return vector<string>();
              }
  		  else
  			 i++;
  		  }
  		else if (parsingsigns.checksign(str[i]) != -1)
  		  {
            if (i-anf > 0)
  		    hilfe.push_back(substr(anf,i-anf));
  		  i++;
  		  anf = i;
  		  }
  		else
  		  i++;
  		}  // end: while (i <len)
  	 if (anf < len)
  		hilfe.push_back(substr(anf,len-anf));
  	 }

    return hilfe;

    }


  //vector<string> strtoken2_quot(const string & parsingsigns,
  //                                bool & bracketmiss, bool & quotmiss) const;
  vector<string> strtoken2_quot(const string & parsingsigns,
                                        bool & bracketmiss,bool & quotmiss) const
    {
    vector<string>  hilfe;

    bracketmiss = false;
    quotmiss = false;

    if (len > 0)
  	 {
  	 int i=0;                             // looping variable
  	 int anf = 0;
  	 while (i < len)
  		{
  		if (str[i] == '(')
  		  {
  		  i = closingbracketpos(i);
  		  if (i == -1)
              {
              bracketmiss = true;
              return vector<string>();
              }
  		  else
  			 i++;
  		  }
          else if (str[i] == '"')
            {
            i++;
            while ( ( i< len) && (str[i] != '"') )
              i++;

            if ( ( i < len) && (str[i] == '"') )
              {
              }
            else
              quotmiss = true;

            i++;
            }
  		else if (parsingsigns.checksign(str[i]) != -1)
  		  {
            if (i-anf > 0)
  		    hilfe.push_back(substr(anf,i-anf));
  		  i++;
  		  anf = i;
  		  }
  		else
  		  i++;
  		}  // end: while (i <len)
  	 if (anf < len)
  		hilfe.push_back(substr(anf,len-anf));
  	 }

    return hilfe;

    }


  // FUNCTION: strtoken
  // TASK: returns the tokens of the string as a vector of strings
  //       delimeters are stored in 'parsingtokens'

  //vector<string> strtoken(const vector<string> & parsingtokens) const;
  vector<string> strtoken(const vector<string> & parsingtoken) const
    {

    vector<string> hilfe;

    vector<string> token = strtoken(" ");

    unsigned i = 0;
    while (i < token.size())
  	 {
  	 string h = token[i];
  	 i++;
  	 if (h.isinlist(parsingtoken) == -1)
  		while ((i < token.size()) && (token[i].isinlist(parsingtoken) == -1))
  		  {
  		  h = h + " " + token[i];
  		  i++;
  		  }
  	 hilfe.push_back(h);
  	 }
    return hilfe;
  }

  // FUNCTION: strtoken2
  // TASK: returns the tokens of the string as a vector of strings
  //       delimeters are stored in 'parsingtokens'
  // ADDITIONAL INFORMATION:
  // Difference to strtoken: signs in brackets are ignored
  // If closing brackets are missing, bracketmiss = true, else false

  //vector<string> strtoken2(const vector<string> & parsingtokens,
  //                         bool & bracketmiss) const;
  vector<string> strtoken2(const vector<string> & parsingtoken,
                                   bool & bracketmiss) const
    {
    vector<string> hilfe;

    vector<string> token = strtoken2(" ",bracketmiss);

    if (!bracketmiss)
    {
    unsigned i = 0;
    while (i < token.size())
  	 {
  	 string h = token[i];
  	 i++;
  	 if (h.isinlist(parsingtoken) == -1)
  		while ((i < token.size()) && (token[i].isinlist(parsingtoken) == -1))
  		  {
  		  h = h + " " + token[i];
  		  i++;
  		  }
  	 hilfe.push_back(h);
  	 }
     }

    return hilfe;

  }


  //vector<string> strtoken2_quot(const vector<string> & parsingtokens,
  //                                    bool & bracketmiss,bool & quotmiss) const;
  vector<string> strtoken2_quot(const vector<string> & parsingtoken,
                                        bool & bracketmiss, bool & quotmiss) const
    {
    vector<string> hilfe;

    vector<string> token = strtoken2_quot(" ",bracketmiss,quotmiss);

    if ((!bracketmiss) && (!quotmiss))
     {
     unsigned i = 0;
     while (i < token.size())
       {
  	 string h = token[i];
  	 i++;
  	 if (h.isinlist(parsingtoken) == -1)
  		while ((i < token.size()) && (token[i].isinlist(parsingtoken) == -1))
  		  {
  		  h = h + " " + token[i];
  		  i++;
  		  }
  	 hilfe.push_back(h);
  	 }
     }

    return hilfe;

  }


  //list<string> strtokenlist(const string & parsingsigns,bool  inclsigns = true) const;
  list<string> strtokenlist(const string & parsingsigns,bool inclsigns = true) const
    {
    list<string> help;                     // will be returned at the
  													  // of the function
    if (len > 0)
  	 {
  	 int i=0;                           // looping variable
       while (i < len)
  		{
  		if (parsingsigns.checksign(str[i]) != -1)
  		  if (str[i] == ' ')
  			 while ( (i < len) && (str[i] == ' ') )
  				i++;
  		  else
  			 {
               if (inclsigns)
  			   help.push_back(substr(i,1));
  			 i++;
  			 }
  		else
  		  {
  		  int anf = i;
  		  while ( (i < len) && (parsingsigns.checksign(str[i]) == -1) )
              i++;

  		  help.push_back(substr(anf,i-anf));
  		  }
  		}  // end: while (i < text.length())
  	 }  // end: if (text.length() > 0)
    return help;
  }

  // FUNCTION: endswith
  // TASK: checks if the calling string ends with 'c'

  //bool endswith(const char * c) const;
  bool endswith(const char * c) const
      {
      int lenc = strlen(c);
      bool endwith = true;
      for(int i=0;i<lenc;i++)
        if(c[lenc-1-i] != str[len-1-i])
          return false;
      return true;
    }


  };  // end: class string



string __EXPORT_TYPE outresults(const unsigned & l,const string & name,
                  const double & mean,
                  const double & std, const double & qu10, const double & qu50,
                  const double & qu90);

string __EXPORT_TYPE make_latextable(vector<string> & v);



//------------ functions, that convert a number into a string ------------------

// FUNCTION: inttostring
// TASK: converts the integer number 'value' into a string

string __EXPORT_TYPE inttostring(int value);

// FUNCTION: doubletostring
// TASK: converts the double number 'value' into a string

string __EXPORT_TYPE doubletostring(double value,int dec=15);


//------------------------ type definitons -------------------------------------

typedef vector<string> __EXPORT_TYPE stringvec;
typedef list<string> __EXPORT_TYPE stringlist;

/*
template<class T>
void pr(const vector<T> & v)
  {
  if (! v.empty())
	 for (int i = 0;i < v.size();i++)
		cout << v[i] << "\n";
  cout << "\n";
  }

template<class T>
void pr(list<T> & v)
  {
  if (!v.empty())
	 {
	 list<T>::iterator i;
	 for(i=v.begin();i != v.end();i++)
		cout << (*i) << endl;
	 }
  }
*/


}  // end: namespace ST

#endif


