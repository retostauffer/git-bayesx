//---------------------------------------------------------------------------
#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (kriging_INCLUDED)
#define kriging_INCLUDED

#include<fullcond.h>
#include<mcmc_nonpbasis.h>

namespace MCMC
{

//------------------------------------------------------------------------------
//------------------------------- class: kriging -------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE FULLCOND_kriging : public FULLCOND_nonp_basis
  {

  protected:

  unsigned nrknots;
  double nu;
  double rho;
  double maxdist;
  bool full;
  bool spacefill;

  MAP::map m;                         // Variablen für geokriging
  bool mapexisting;
  ST::string mapname;
  vector<ST::string> regionnames;

  double p;                           // p und q für Space-Filling-Algorithmus
  double q;
  unsigned maxsteps;

  vector<double> xknots;              // x-und y-Koordinaten der Knoten
  vector<double> yknots;

  vector<double> xvalues;             // unterschiedliche Werte der Kovariablen
  vector<double> yvalues;

  datamatrix xorig;                   // Original x- und y-Variable.
  datamatrix yorig;

  vector<int> index2;

  vector<int>freq;
  vector<int>freqoutput;
  unsigned nrdiffobs;

  void make_index(const datamatrix & var1,const datamatrix & var2);
  void make_xy_values(const datamatrix & var1,const datamatrix & var2);
  void compute_knots(const vector<double> & xvals,const vector<double> & yvals);

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_kriging(void) : FULLCOND_nonp_basis()
    {
    }

  // COPY CONSTRUCTOR

  FULLCOND_kriging(const FULLCOND_kriging & kr);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_kriging & operator=(const FULLCOND_kriging & kr);

  // Constructor1

  FULLCOND_kriging(MCMCoptions * o, const datamatrix & v1,
               const datamatrix & v2, const datamatrix & knotdata,
               const unsigned & nrk, const double & n, const double & maxd,
               const double & pval, const double & qval, const unsigned & maxst,
               const bool & fu, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl);

  // Constructor2: geokriging

  FULLCOND_kriging(MCMCoptions * o, const datamatrix & region,
               const MAP::map & mp, const ST::string & mn, const datamatrix & knotdata,
               const unsigned & nrk, const double & n, const double & maxd,
               const double & pval, const double & qval, const unsigned & maxst,
               const bool & fu, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl);
  // DESTRUCTOR

  ~FULLCOND_kriging(){}

  void createreml(datamatrix & X,datamatrix & Z,const unsigned & Xpos,
                  const unsigned & Zpos);
                  
  double outresultsreml(datamatrix & X,datamatrix & Z,
                                     datamatrix & betareml,
                                     datamatrix & betacov,
                                     datamatrix & thetareml,
                                     const unsigned & Xpos,
                                     const unsigned & Zpos,
                                     const unsigned & thetapos,
                                     const bool & dispers,
                                     const unsigned & betaXpos,
                                     const unsigned & betaZpos,
                                     const double & category,
                                     const bool & ismultinomial,
                                     const unsigned plotpos);

  void outoptionsreml();

  void init_names(const vector<ST::string> & na);

  ST::string getinfo(void);

  };


} // end: namespace MCMC

//---------------------------------------------------------------------------
#endif
