/* Replace "dll.h" with the name of your header */
#include "dll.h"
#include <windows.h>

#include"../export_type.h"

#include <jni.h>
#include "BayesX.h"

#include <stdio.h>
#include "adminparse.h"

DllClass::DllClass()
{

}


DllClass::~DllClass ()
{

}


// Vorschlag:
administrator admin;// = administrator();

//--------------- Funktionen zum parsen --------------------------------------

JNIEXPORT jboolean JNICALL Java_BayesX_parse
  (JNIEnv *env, jobject obj, jstring str)

  {


  admin.adminb.Java = env;

  admin.adminb.BayesX_obj = obj;
  admin.adminb.BayesX_cls = admin.adminb.Java->GetObjectClass(obj);
//  admin.javaoutput = admin.Java->GetMethodID(admin.BayesX_cls, "JavaOutput", "(Ljava/lang/String;)V");
  admin.adminb.javaoutput = admin.adminb.Java->GetMethodID(admin.adminb.BayesX_cls, "JavaOutput", "(Ljava/lang/String;ZZSIII)V");

  if (admin.adminb.javaoutput == 0) {
     return false;
  }

  const char *helpstr = env->GetStringUTFChars(str, 0);
  ST::string inp = helpstr;
  env->ReleaseStringUTFChars(str, helpstr);

  inp = inp.deleteallsigns(admin.get_delim());
  inp = inp.replaceallsigns('\n',' ');
  inp = inp.replaceallsigns('\r',' ');
  inp = inp.eatwhitespace();

  jboolean stop = false;

  stop = admin.parse(inp);

  return stop;

  }



JNIEXPORT void JNICALL Java_BayesX_parsecommand2
  (JNIEnv *env, jobject obj, jobject oA, jstring str)
      {

        jclass cls = env->GetObjectClass(oA);
        jmethodID mid = env->GetMethodID(cls, "append", "(Ljava/lang/String;)V");
        if (mid == 0) {
            return;
        }

        env->CallVoidMethod(oA, mid, str);
        env->CallVoidMethod(oA, mid, env->NewStringUTF("\nENDE\n"));

      return;

      }


JNIEXPORT void JNICALL Java_BayesX_parsecommand
      (JNIEnv *env, jobject obj, jstring str)
      {

      admin.adminb.Java = env;
      admin.adminb.BayesX_obj = obj;
      admin.adminb.BayesX_cls = admin.adminb.Java->GetObjectClass(obj);
//      admin.javaoutput = admin.Java->GetMethodID(admin.BayesX_cls, "JavaOutput", "(Ljava/lang/String;)V");
//      admin.javaoutput = admin.Java->GetMethodID(admin.BayesX_cls, "JavaOutput", "(Ljava/lang/String;ZZS)V");
      admin.adminb.javaoutput = admin.adminb.Java->GetMethodID(
      admin.adminb.BayesX_cls, "JavaOutput", "(Ljava/lang/String;ZZSIII)V");
      if (admin.adminb.javaoutput == 0) {
         return;
      }

      bool stop = false;

      const char *helpstr = env->GetStringUTFChars(str, 0);
      ST::string inp = helpstr;
      env->ReleaseStringUTFChars(str, helpstr);

      inp = inp.deleteallsigns(admin.get_delim());
      inp = inp.replaceallsigns('\n',' ');
      inp = inp.replaceallsigns('\r',' ');
      inp = inp.eatwhitespace();

      stop = admin.parse(inp);


      return;

      }

// --------------------- Zugang zu Objekten und Objekttypen -----------------

JNIEXPORT void JNICALL Java_BayesX_setObjectList
  (JNIEnv *env, jobject obj, jobject v, jstring type)

  {


  jclass cls = env->GetObjectClass(obj);

  jmethodID addvector = env->GetMethodID(cls,"addtoVector","(Ljava/util/Vector;Ljava/lang/String;)V");

  if (addvector == 0) {

     return;

  }

  const char* t = env->GetStringUTFChars(type,0);

  for(unsigned i=0;i<admin.get_objects().size();i++)

    {

    if(strcmp(admin.get_objects()[i]->gettype().strtochar(),t) == 0)

      {

      env->CallVoidMethod(obj,addvector,v,env->NewStringUTF(admin.get_objects()[i]->getname().strtochar()));

      }

    }


  return;

  }

JNIEXPORT void JNICALL Java_BayesX_setObjectTypeList
  (JNIEnv *env, jobject obj, jobject v)

  {


  jclass cls = env->GetObjectClass(obj);

  jmethodID addvector = env->GetMethodID(cls,"addtoVector","(Ljava/util/Vector;Ljava/lang/String;)V");

  if (addvector == 0) {

     return;

  }

  for(unsigned i=0;i<admin.get_objecttype().size();i++)

    {

    env->CallVoidMethod(obj,addvector,v,env->NewStringUTF(admin.get_objecttype()[i].strtochar()));

    }


  return;

  }


// -------------------- Zugang zu 'dataset' ---------------------------------

JNIEXPORT void JNICALL Java_BayesX_setVarnames
  (JNIEnv *env, jobject obj, jobject v)

  {


  jclass cls = env->GetObjectClass(obj);

  jmethodID addvector = env->GetMethodID(cls,"addtoVector","(Ljava/util/Vector;Ljava/lang/String;)V");

  if (addvector == 0) {

     return;

  }

  list<ST::string> liste = admin.adminp.get_datap()->getVarnames();

  list<ST::string>::iterator it = liste.begin();
  for(unsigned i=0;i<admin.adminp.get_datap()->getVarnames().size();i++,it++)

    {

    env->CallVoidMethod(obj,addvector,v,env->NewStringUTF((*it).strtochar()));

    }


  return;

  }

JNIEXPORT jstring JNICALL Java_BayesX_getValue
  (JNIEnv *env, jobject obj, jint i, jint j)
  {
  jstring value;
  if(i<admin.adminp.get_datap()->obs() && j<admin.adminp.get_datap()->getVarnames().size())
    {
    admin.adminp.get_datap()->set_iterator(j+1);
    double v = admin.adminp.get_datap()->getvalue(i);
    if(v==MAXDOUBLE)
      value = env->NewStringUTF(".");
    else
      value = env->NewStringUTF(ST::doubletostring(v,8).strtochar());
    }
  else
    {
    value = env->NewStringUTF("");
    }
  return value;
  }

JNIEXPORT jdouble JNICALL Java_BayesX_getDoubleValue
  (JNIEnv *env, jobject obj, jint i, jint j)

  {

  jdouble value;

  value = (*(admin.adminp).get_Dp())(i,j);
  return value;
  }


JNIEXPORT jint JNICALL Java_BayesX_getRows
  (JNIEnv *env, jobject obj)
  {
  jint rows = admin.adminp.get_datap()->obs();
  return rows;
  }

JNIEXPORT jint JNICALL Java_BayesX_getDRows
  (JNIEnv *, jobject)
  {
  jint rows = (*(admin.adminp).get_Dp()).rows();
  return rows;
  }

JNIEXPORT jint JNICALL Java_BayesX_getDCols
  (JNIEnv *, jobject)
  {
  jint cols = (*(admin.adminp).get_Dp()).cols();
  return cols;
  }

JNIEXPORT jdouble JNICALL Java_BayesX_getMax
  (JNIEnv *, jobject, jint col)

  {

  jdouble max = (*(admin.adminp).get_Dp()).max(col);

  return max;

  }


JNIEXPORT jdouble JNICALL Java_BayesX_getMin
  (JNIEnv *, jobject, jint col)

  {

  jdouble min = (*(admin.adminp).get_Dp()).min(col);

  return min;

  }


// ------------------------------------------------------------------

JNIEXPORT jstring JNICALL Java_BayesX_getVarname
  (JNIEnv *env, jobject obj, jint i)
  {
  jstring varname;
  if(i<(*(admin.adminp).get_varnamesp()).size())
    varname = env->NewStringUTF((*(admin.adminp).get_varnamesp())[i].strtochar());
  else
    varname = env->NewStringUTF("");
  return varname;
  }

// ---------------------- Stop, Pause etc. --------------------------

JNIEXPORT void JNICALL Java_BayesX_setStop
  (JNIEnv *env, jobject obj, jboolean stop)

  {

  admin.adminb.set_stop(stop);

  }


JNIEXPORT void JNICALL Java_BayesX_setPause

  (JNIEnv *env, jobject obj, jboolean pause)

  {

  admin.adminb.set_pause(pause);

  }


JNIEXPORT void JNICALL Java_BayesX_setProcessrunning

  (JNIEnv *env, jobject obj, jboolean processrunning)

  {

  admin.adminb.set_processrunning(processrunning);

  }


JNIEXPORT void JNICALL Java_BayesX_setSuppressoutput

  (JNIEnv *env, jobject obj, jboolean suppressoutput)

  {

  admin.adminb.set_suppressoutput(suppressoutput);

  }



// ---------------------- Zugang zu 'map' -----------------------------------

JNIEXPORT void JNICALL Java_BayesX_getline
  (JNIEnv *env, jobject obj, jdoubleArray line, jint i, jint j, jint k)

  {

  if(env->GetArrayLength(line)!=4)

    return;

  jdouble *help = env->GetDoubleArrayElements(line, 0);

  help[0] = admin.adminp.get_mapinfo()->get_region(i).get_polygone(j).get_line(k).x1;

  help[1] = admin.adminp.get_mapinfo()->get_region(i).get_polygone(j).get_line(k).x2;

  help[2] = admin.adminp.get_mapinfo()->get_region(i).get_polygone(j).get_line(k).y1;

  help[3] = admin.adminp.get_mapinfo()->get_region(i).get_polygone(j).get_line(k).y2;

  env->ReleaseDoubleArrayElements(line, help, 0);

  return;

  }


JNIEXPORT void JNICALL Java_BayesX_getboundaries
  (JNIEnv *env, jobject obj, jdoubleArray bnd)

  {

  if(env->GetArrayLength(bnd)!=4)

    return;

  jdouble *help = env->GetDoubleArrayElements(bnd, 0);

  help[0] = admin.adminp.get_mapinfo()->get_minX();

  help[1] = admin.adminp.get_mapinfo()->get_minY();

  help[2] = admin.adminp.get_mapinfo()->get_maxX();

  help[3] = admin.adminp.get_mapinfo()->get_maxY();

  env->ReleaseDoubleArrayElements(bnd, help, 0);

  return;

  }

JNIEXPORT jint JNICALL Java_BayesX_getnrregions
  (JNIEnv *env, jobject obj)

  {

  return admin.adminp.get_mapinfo()->get_nrregions();

  }


JNIEXPORT jint JNICALL Java_BayesX_getnrpoly
  (JNIEnv *env, jobject obj, jint i)

  {

  if(i<admin.adminp.get_mapinfo()->get_nrregions())

    return admin.adminp.get_mapinfo()->get_region(i).get_nrpoly();

  else

    return 0;

  }

JNIEXPORT jint JNICALL Java_BayesX_getnrlines
  (JNIEnv *env, jobject obj, jint i, jint j)

  {

  if(i<admin.adminp.get_mapinfo()->get_nrregions() &&

     j<admin.adminp.get_mapinfo()->get_region(i).get_nrpoly())

    return admin.adminp.get_mapinfo()->get_region(i).get_polygone(j).get_nrlines();

  else

    return 0;

  }

JNIEXPORT jboolean JNICALL Java_BayesX_isin
  (JNIEnv *, jobject, jint i)

  {

  assert(i<admin.adminp.get_mapinfo()->get_nrregions());

  if(admin.adminp.get_mapinfo()->get_region(i).get_isin() == "")
     return false;
  else
     return true;
  }

JNIEXPORT jstring JNICALL Java_BayesX_getregionname
  (JNIEnv *env, jobject, jint i)

  {

  assert(i<admin.adminp.get_mapinfo()->get_nrregions());

  jstring varname;

  varname = env->NewStringUTF(admin.adminp.get_mapinfo()->get_region(i).get_name().strtochar());
  return varname;
  }

JNIEXPORT void JNICALL Java_BayesX_getcentroid
  (JNIEnv *env, jobject, jdoubleArray centroid, jint i)

  {

  if(env->GetArrayLength(centroid)!=2)

    return;
  if(i<admin.adminp.get_mapinfo()->get_nrregions())
    {

    jdouble *help = env->GetDoubleArrayElements(centroid, 0);
    help[0] = admin.adminp.get_mapinfo()->get_region(i).get_xcenter();
    help[1] = admin.adminp.get_mapinfo()->get_region(i).get_ycenter();
    env->ReleaseDoubleArrayElements(centroid, help, 0);
    }
  return;
  }


JNIEXPORT jdouble JNICALL Java_BayesX_getname
  (JNIEnv *, jobject, jint i)
  {
  assert(i<admin.adminp.get_mapinfo()->get_nrregions());
  jdouble name;
  admin.adminp.get_mapinfo()->get_region(i).get_name().strtodouble(name);
  return name;
  }


BOOL APIENTRY DllMain (HINSTANCE hInst     /* Library instance handle. */ ,
                       DWORD reason        /* Reason this function is being called. */ ,
                       LPVOID reserved     /* Not used. */ )
{
    switch (reason)
    {
      case DLL_PROCESS_ATTACH:
        break;

      case DLL_PROCESS_DETACH:
        break;

      case DLL_THREAD_ATTACH:
        break;

      case DLL_THREAD_DETACH:
        break;
    }

    /* Returns TRUE on success, FALSE on failure */
    return TRUE;
}
