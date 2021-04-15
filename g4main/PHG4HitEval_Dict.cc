// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4HitEval_Dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "PHG4HitEval.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4HitEval(void *p = 0);
   static void *newArray_PHG4HitEval(Long_t size, void *p);
   static void delete_PHG4HitEval(void *p);
   static void deleteArray_PHG4HitEval(void *p);
   static void destruct_PHG4HitEval(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4HitEval*)
   {
      ::PHG4HitEval *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4HitEval >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4HitEval", ::PHG4HitEval::Class_Version(), "", 26,
                  typeid(::PHG4HitEval), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4HitEval::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4HitEval) );
      instance.SetNew(&new_PHG4HitEval);
      instance.SetNewArray(&newArray_PHG4HitEval);
      instance.SetDelete(&delete_PHG4HitEval);
      instance.SetDeleteArray(&deleteArray_PHG4HitEval);
      instance.SetDestructor(&destruct_PHG4HitEval);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4HitEval*)
   {
      return GenerateInitInstanceLocal((::PHG4HitEval*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4HitEval*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4HitEval::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4HitEval::Class_Name()
{
   return "PHG4HitEval";
}

//______________________________________________________________________________
const char *PHG4HitEval::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4HitEval*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4HitEval::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4HitEval*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4HitEval::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4HitEval*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4HitEval::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4HitEval*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4HitEval::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4HitEval.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4HitEval::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4HitEval::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4HitEval(void *p) {
      return  p ? new(p) ::PHG4HitEval : new ::PHG4HitEval;
   }
   static void *newArray_PHG4HitEval(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4HitEval[nElements] : new ::PHG4HitEval[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4HitEval(void *p) {
      delete ((::PHG4HitEval*)p);
   }
   static void deleteArray_PHG4HitEval(void *p) {
      delete [] ((::PHG4HitEval*)p);
   }
   static void destruct_PHG4HitEval(void *p) {
      typedef ::PHG4HitEval current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4HitEval

namespace {
  void TriggerDictionaryInitialization_PHG4HitEval_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4HitEval_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4HitEval;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4HitEval_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $Id: $

/*!
 * \file PHG4HitEval.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4MAIN_PHG4HITEVAL_H
#define G4MAIN_PHG4HITEVAL_H

#include "PHG4Hitv1.h"

class PHG4Hit;
class PHObject;

/*!
 * \brief PHG4HitEval for evaluating PHG4Hitv1 in CINT readable evaluation trees
 */
class PHG4HitEval : public PHG4Hitv1
{
 public:
  PHG4HitEval();

  PHG4HitEval(const PHG4Hit *g4hit);

  virtual ~PHG4HitEval() {}

  virtual void CopyFrom(const PHObject *phobj);

  float
  get_eion() const
  {
    return eion;
  }
  void
  set_eion(const float f)
  {
    eion = f;
  }

  int get_scint_id() const
  {
    return scint_id;
  }

  void
  set_scint_id(const int i)
  {
    scint_id = i;
  }

  float
  get_light_yield() const
  {
    return light_yield;
  }

  void
  set_light_yield(float lightYield)
  {
    light_yield = lightYield;
  }

  float
  get_path_length() const
  {
    return path_length;
  }

  void
  set_path_length(float pathLength)
  {
    path_length = pathLength;
  }

 protected:
  float eion;

  int scint_id;

  //! a number proportional to the scintillation light yield.
  float light_yield;

  //! path length of the track to the hit
  float path_length;

  ClassDef(PHG4HitEval, 1)
};

#endif /* PHG4HITEVAL_H_ */

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4HitEval", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4HitEval_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4HitEval_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4HitEval_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4HitEval_Dict() {
  TriggerDictionaryInitialization_PHG4HitEval_Dict_Impl();
}
