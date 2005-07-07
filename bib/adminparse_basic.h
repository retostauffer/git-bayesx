

#ifdef __BUILDING_THE_DLL
#define __EXPORT_TYPE __export
#else
#define __EXPORT_TYPE __import
#endif


//---------------------------------------------------------------------------
#ifndef ADMINPARSEBASIC
#define ADMINPARSEBASIC
//---------------------------------------------------------------------------

#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include "StatReview.h"
#include<StatwinFrame.h>
#include<statwin_haupt.h>
#elif defined(JAVA_OUTPUT_WINDOW)
#include<jni.h>
#endif


//------------------------------------------------------------------------------


class __EXPORT_TYPE administrator_basic
  {

  private:

  //------------------------ PRIVATE VARIABLES ---------------------------------


  public:

  bool pause;
  bool stop;
  bool processrunning;
  bool suppressoutput;


#if defined(JAVA_OUTPUT_WINDOW)
  JNIEnv* Java;
  jclass BayesX_cls;
  jobject BayesX_obj;
  jmethodID javaoutput;
#endif

  // CONSTRUCTOR

  administrator_basic(void);

  // DESTRUCTOR

  ~administrator_basic() {}

  bool breakcommand(void);

  bool get_pause(void)
    {
    return pause;
    }

  bool get_stop(void)
    {
    return stop;
    }

  bool get_processrunning(void)
    {
    return processrunning;
    }

  bool get_suppressoutput(void)
    {
    return suppressoutput;
    }

  void set_pause(bool b)
    {
    pause = b;
    }

  void set_stop(bool b)
    {
    stop = b;
    }

  void set_processrunning(bool b)
    {
    processrunning = b;
    }

  void set_suppressoutput(bool b)
    {
    suppressoutput = b;
    }


 
  };


#endif
