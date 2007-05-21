/******************************************************************************/
/* Hier werden die Bedingungen für die GNU Version definiert. (Bei der        */
/* C++-Version ist BORLAND_OUTPUT_WINDOW definiert, bei der JAVA-Version ist  */
/* __BUILDING_THE_DLL definiert (in Borland), andernfalls werden die          */
/* Bedingungen _BUILDING_GNU, JAVA_OUTPUT_WINDOW, TEMPL_INCL_DEF, _MSC_VER2   */
/* und NO_TEMPLATE_FRIENDS hier definiert.)                                   */
/******************************************************************************/

#if !defined (__BUILDING_THE_DLL) && !defined (BORLAND_OUTPUT_WINDOW)

#if !defined (__BUILDING_GNU)
#define __BUILDING_GNU
#endif

//#if !defined (JAVA_OUTPUT_WINDOW)
//#define JAVA_OUTPUT_WINDOW
//#endif

#if !defined (TEMPL_INCL_DEF)
#define TEMPL_INCL_DEF
#endif

#if !defined (_MSC_VER2)
#define _MSC_VER2
#endif

#if !defined (NO_TEMPLATE_FRIENDS)
#define NO_TEMPLATE_FRIENDS
#endif

#if !defined (INCLUDE_REML)
#define INCLUDE_REML
#endif

#if !defined (INCLUDE_MCMC)
#define INCLUDE_MCMC
#endif

#endif


