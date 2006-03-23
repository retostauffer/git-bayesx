
#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __declspec(dllexport)
#else
#define __EXPORT_TYPE __import
#endif

#ifndef ADMINPARSEPOINTER
#define ADMINPARSEPOINTER

#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include "StatReview.h"
#include<StatwinFrame.h>
#include<statwin_haupt.h>
#endif

#include<fstream.h>
#include<string.h>
#include"data.h"
#include"map.h"
//#include<dir.h>
//------------------------------------------------------------------------------


class __EXPORT_TYPE administrator_pointer
  {

  private:

  //------------------------ PRIVATE VARIABLES ---------------------------------

  dataset * datap;
  MAP::map * mapinfo;
  datamatrix * Dp;
  vector<ST::string> * varnamesp;


  public:


  // CONSTRUCTOR

  administrator_pointer(void);

  // DESTRUCTOR

  ~administrator_pointer() {}

  dataset * get_datap(void)
    {
    return datap;
    }

  void set_datap(dataset *d)
    {
    datap = d;
    }

  MAP::map * get_mapinfo(void)
    {
    return mapinfo;
    }

  void set_mapinfo(MAP::map *m)
    {
    mapinfo = m;
    }

  datamatrix * get_Dp(void)
    {
    return Dp;
    }

  void set_Dp(datamatrix *D)
    {
    Dp = D;
    }

  vector<ST::string> * get_varnamesp(void)
    {
    return varnamesp;
    }

  void set_varnamesp(vector<ST::string> *varn)
    {
    varnamesp = varn;
    }


  };


#endif
