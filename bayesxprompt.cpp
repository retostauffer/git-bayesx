// g++ -I/usr/local/include bayesxprompt.cpp -c -o bayesxprompt.o
// g++ -I/usr/local/include bayesxprompt.cpp clstring.so -o bayesxprompt

#if !defined (__BUILDING_GNU)
#define __BUILDING_GNU
#endif

#if !defined (TEMPL_INCL_DEF)
#define TEMPL_INCL_DEF
#endif

#if !defined (_MSC_VER2)
#define _MSC_VER2
#endif

#if !defined (NO_TEMPLATE_FRIENDS)
#define NO_TEMPLATE_FRIENDS
#endif


#include "clstring.h"
#include "random.h"
#include "graph.h"
#include "map.h"
#include "statmat.h"
#include "model.h"
#include "envmatrix.h"
#include "envmatrix_penalty.h"
#include "bandmat.h"
#include "adminparse_gnu.h"
#include <iostream>
#include <string>

int main()
  {
  // terminating commands
  ST::string* stop1 = new ST::string("quit") ;
  ST::string* stop2 = new ST::string("exit") ;

  bool run=false;
  admin_gnu a;

  while(!run)
    {
    std::cout <<"BayesX>";

    // read from the command line
    char array[256];
    std::cin.getline(array, sizeof(array), '\n');
    const char* p=array;
    ST::string* s=new ST::string(p) ;


/*    int e= s->firstpos('e');
    std::cout << "  Erstes e in '" << array << "' an Position " << e << "\n";
    std::cout << "  pi = " << ST::doubletostring(3.14) << "\n";

    std::cout << "  phi(0)=" << randnumbers::phi(0.0) << "\n";

    std::cout << "  Sample from IG(0.01,0.01): " << randnumbers::rand_gamma(0.01,0.01) << "\n";

    graph g = graph();
    std::cout << "  g.getlinenr(): " << g.getlinenr() << "\n";

    MAP::map m = MAP::map();
    std::cout << "  m.get_maxn(): " << m.get_maxn() << "\n";
    
    datamatrix d = datamatrix(2,2,1);
    datamatrix d2 = datamatrix(2,1,3);

    datamatrix d3 = multdiagback(d,d2);
    std::cout << "  d3: " << d3 << "\n";

    modelStandard mod = modelStandard();
    std::cout << "  mod.getModelText(): " << mod.getModelText() << "\n";

    symbandmatrix<double> sym = symbandmatrix<double>();

    envmatrix<double> env = Kseasonenv(2,5);
    std::cout << "  env.getBandwidth(): " << env.getBandwidth() << "\n";*/

    run = a.parse(*s);

    // check for terminating condition
//    if(*s==*stop1 || *s==*stop2) 
//      run=-1;
    } //end while run 
  }

