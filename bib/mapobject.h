
#if !defined (MAPOBJECT_INCLUDED)

#define MAPOBJECT_INCLUDED

#include"../export_type.h"
#if defined(JAVA_OUTPUT_WINDOW)
#include"adminparse_pointers.h"
#endif

#include"map.h"
#include"statobj.h"
#include"dataobj.h"
#include"graph.h"


//------------------------------------------------------------------------------
//------------------------- CLASS: mapobject -----------------------------------
//------------------------------------------------------------------------------


// HINZUFÜGEN EINER NEUEN METHODE

// 1. private - Teil ergänzen
//   - optionlist für die neue Methode definieren
//   - optionen definieren
//   - Modell Klasse definieren
//   - use Klasse definieren
//   - run Funktion definieren
// 2. create - Teil ergänzen (in mapobject.cpp)
// 3. run funktion schreiben


//using MAP::map;


class __EXPORT_TYPE mapobject : public statobject
  {

  private:

  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_pointer * adminp_p;
  #endif

  // vector of pointers to current statobjects

  vector<statobject*> * statobj;

  // mapexisting is true, if a map is currently existing

  bool mapexisting;

  // map object

  MAP::map mapinfo;

  // pointer to functions

  typedef void ( *runpointer ) (mapobject & m);

  runpointer functions[5];

  //----------------------------------------------------------------------------
  //------------------------- for method 'infile' ------------------------------
  //----------------------------------------------------------------------------

  modelStandard mod;
  usePathRead uread;
  vector<ST::string> weightdefs;
  stroption weightdef;
  simpleoption neighbors;
  simpleoption graf;
  simpleoption centroids;
  optionlist infileoptions;

  friend void infilerun(mapobject & m);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  //------------------------ for method createmap ------------------------------
  //----------------------------------------------------------------------------

  modelStandard modcreate;
  use udata;
  vector<ST::string> weightdefscr;
  stroption weightdefcr;
  doubleoption maxdif;
  optionlist createmapoptions;

  friend void createmaprun(mapobject & m);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  //---------------------- for method 'computeneighbors' -----------------------
  //----------------------------------------------------------------------------

  modelStandard modn;
  usePathWrite uwrite;
  simpleoption replace;
  optionlist computeneighborsoptions;

  friend void computeneighborsrun(mapobject & m);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  //--------------------------- for method 'outfile' ---------------------------
  //----------------------------------------------------------------------------

  modelStandard modo;
  usePathWrite uwriteo;
  simpleoption grafo;
  simpleoption centroidso;
  simpleoption replaceo;
  simpleoption includeweights;
  optionlist outfileoptions;

  friend void outfilerun(mapobject & m);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //-------------------------- for method 'reorder' ----------------------------
  //----------------------------------------------------------------------------

  modelStandard modre;
  use udatare;

  optionlist reorderoptions;

  friend void reorderrun(mapobject & m);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // FUNCTION create
  // TASK: intitializes methods

  void create(void);


  public:

  // DEFAULT CONSTRUCTOR:

  mapobject(void) : statobject()
	 {
	 type = "map";
	 mapexisting = false;
	 }

  // CONSTRUCTOR

  #if defined(JAVA_OUTPUT_WINDOW)
  mapobject(administrator_basic * adb, administrator_pointer * adp,
            const ST::string & n,ofstream * lo,istream * i,ST::string p,
            vector<statobject*> * st);
  #else
  mapobject(const ST::string & n,ofstream * lo,istream * i,ST::string p,
            vector<statobject*> * st);
  #endif

  // COPY CONSTRUCTOR

  mapobject(const mapobject & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const mapobject & operator=(const mapobject & m);

  int parse(const ST::string & c);

  MAP::map getmap(void) const
    {
    return mapinfo;
    }


  void describe(const optionlist & globaloptions = optionlist());

  bool getexisting(void) const
    {
    return mapexisting;
    }


  };

#if defined (__BUILDING_GNU)
// ----------------------- forward friends decls ---------------------

void infilerun(mapobject & m);
void createmaprun(mapobject & m);
void computeneighborsrun(mapobject & m);
void outfilerun(mapobject & m);
void reorderrun(mapobject & m);
#endif

#endif
