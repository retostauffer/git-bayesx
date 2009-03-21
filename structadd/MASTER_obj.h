
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (MASTEROBJ)

#define MASTEROBJ

#include"distr.h"

#if defined(JAVA_OUTPUT_WINDOW)
#include"adminparse_basic.h"
#endif


namespace MCMC
{


//------------------------------------------------------------------------------
//--------------------------- CLASS: MASTER_OBJ --------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE MASTER_OBJ
  {


  protected:


  public:

  DISTR * level1_likep;

  // DEFAULT CONSTRUCTOR

  MASTER_OBJ(void);

  // COPY CONSTRUCTOR

  MASTER_OBJ(const MASTER_OBJ & o);

  // OVERLOADED ASSIGNMENT OPERATOR

  const MASTER_OBJ & operator=(const MASTER_OBJ & o);


  // DESTRUCTOR

  ~MASTER_OBJ() {}


  };

//------------------------------------------------------------------------------
//------------------------ End: CLASS MASTER_OBJ -------------------------------
//------------------------------------------------------------------------------

} // end: namespace MCMC

#endif