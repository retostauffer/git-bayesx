
//---------------------------------------------------------------------------
#ifndef statwin_hauptH
#define statwin_hauptH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include<fstream.h>
#include<errorm.h>
#include<data.h>
#include<statobj.h>
#include<dataobj.h>

#if defined(INCLUDE_MCMC)
#include<bayesreg.h>
#endif

#if defined(INCLUDE_REML)
#include<remlreg.h>
#endif

#if defined(INCLUDE_STEP)
#include<stepwisereg.h>
#endif

#include<superbayesreg.h>

#include<mapobject.h>
#if defined(INCLUDE_DAG)
#include<dagobject.h>
#endif
//#include<diseaseobj.h>
#include <ComCtrls.hpp>
//---------------------------------------------------------------------------
class Thauptformular : public TForm
{
__published:    // Komponenten, die von der IDE verwaltet werden
    TMemo *commandedit;
    void __fastcall commandeditKeyPress(TObject *Sender, char &Key);



    void __fastcall FormCreate(TObject *Sender);

    void __fastcall commandeditKeyUp(TObject *Sender, WORD &Key,
          TShiftState Shift);
private:        // Benutzerdeklarationen

  //------------------------ PRIVATE VARIABLES ---------------------------------

  ST::string cmd;          // command

  char delim;              // sign that indicates the end of a command

  ofstream logout;

  istream * input;

  ST::string defaultpath;

  bool logfileopen;

  ST::string logfilepath;

  // contains (valid) objecttyps
  // valid types:
  // - dataset
  // - bayesreg
  // - mcmcreg
  // - map
  // - dag
  // - diseasemap

  vector<ST::string> objecttyps;

  // contains pointers to current (stat-)objects

  vector<statobject*> objects;

  // contains current errmormessages

//  errorm::messages errormessages;
  vector<ST::string> errormessages;

  // 'dataobjects' contains current dataobjects

  vector<dataobject> dataobjects;

  // 'bayesregobjects' contains current bayesreg objects

#if defined(INCLUDE_MCMC)
  vector<bayesreg> bayesregobjects;
#endif

  // 'mcmcregobjects' contains current bayesreg objects

  vector<superbayesreg> mcmcregobjects;

  // 'remlregobjects' contains current remlreg objects

#if defined(INCLUDE_REML)
  vector<remlreg> remlregobjects;
#endif

  // 'stepwiseregobjects' contains current stepwisereg objects

#if defined(INCLUDE_STEP)
  vector<stepwisereg> stepwiseregobjects;
#endif

  // 'mapobjects' contains current map objects

  vector<mapobject> mapobjects;

  // 'dagobjects' contains current dag objects

#if defined(INCLUDE_DAG)
  vector<dagobject> dagobjects;
#endif

  // 'diseasemapobjects' contains current diseasemap objects

//  vector<diseaseobj> diseasemapobjects;


  //------------------------ PRIVATE FUNCTIONS ---------------------------------

  void out(const ST::string & c);

  void out(const vector<ST::string> & m);

  // FUNCTION: alreadyexisting
  // TASK: returns 'true', if object with objectname 'name' is already existing

  bool alreadyexisting(const ST::string & name);

  // FUNCTION: create

  ST::string create(const ST::string & in);

  void adjustobjects(void);

  // FUNCTION: drop

  void dropobjects(ST::string name,ST::string type);


  // FUNCTION: parseexisting
  // TASK: parses command 'com' for object with name 'objectname'
  //       objectname should be an object, that is still existing
  // POSSIBLE ERRORS:
  // - object with name 'objectname' is not existing
  // - command 'com' is invalid (i.e. contains errors) (depending on the
  //   special structure of object with name 'objectname')

  void parseexisting(const ST::string & objectname,const ST::string & com);

  void parsespecial(const ST::string & com);

  bool parse(ST::string & in);

public:         // Benutzerdeklarationen

  void parsecommand(ST::string & in);

  int CheckOutputSaved(void);

  vector<ST::string> & get_objecttype(void)
    {
    return objecttyps;
    }

  vector<statobject*> & get_objects(void)
    {
    return objects;
    }

  vector<mapobject> & get_mapobjects(void)
    {
    return mapobjects;
    }

  vector<dataobject> & get_dataobjects(void)
    {
    return dataobjects;
    }

  // FUNCTION: breakcommand
  // returns true, if current process should be interrupted

  bool breakcommand(void);

    __fastcall Thauptformular(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE Thauptformular *hauptformular;
//---------------------------------------------------------------------------
#endif

