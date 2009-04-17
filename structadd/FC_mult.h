
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE  __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FCmultINCLUDED)

#define FCmultINCLUDED

#include"../values.h"
#include<fstream>
#include"clstring.h"
#include"FC_nonp.h"
#include"design.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_mult -----------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE FC_mult  : public FC
  {

  protected:

  MASTER_OBJ * masterp;

  bool multexp;

  FC FCmulteffect;
  FC FCmulteffect_mean;

  bool samplemult;
  bool compmeaneffect;

  DESIGN * dp1;
  DESIGN * dp2;
  FC_nonp * FCnp;
  FC_nonp * FCnp2;

  bool RE_update;

  datamatrix effect;

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_mult(void);

  // CONSTRUCTOR

  FC_mult(bool reu,bool mexp=false);

  // COPY CONSTRUCTOR

  FC_mult(const FC_mult & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_mult & operator=(const FC_mult & m);

  // DESTRUCTOR

  ~FC_mult()
    {
    }

  void update(void);

  bool posteriormode(void);

  void outresults(ofstream & out_stata, ofstream & out_R,
                  const ST::string & pathresults);

  void reset(void);

  void set_effectp(DESIGN * d,FC_nonp * fp);

  void set_intp(DESIGN * d,FC_nonp * fp);

  void set_multeffects(MASTER_OBJ * mp,GENERAL_OPTIONS * o,const ST::string & t,
           const ST::string & fp,bool sm,bool sstore,bool meane);

  void compute_autocorr_all(const ST::string & path, unsigned lag,
                                    ofstream & outg) const;

  void update_multeffect(void);

  };


} // end: namespace MCMC

#endif


