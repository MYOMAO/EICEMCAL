// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4Particlev2_Dict
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
#include "PHG4Particlev2.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4Particlev2(void *p = 0);
   static void *newArray_PHG4Particlev2(Long_t size, void *p);
   static void delete_PHG4Particlev2(void *p);
   static void deleteArray_PHG4Particlev2(void *p);
   static void destruct_PHG4Particlev2(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4Particlev2*)
   {
      ::PHG4Particlev2 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4Particlev2 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4Particlev2", ::PHG4Particlev2::Class_Version(), "", 15,
                  typeid(::PHG4Particlev2), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4Particlev2::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4Particlev2) );
      instance.SetNew(&new_PHG4Particlev2);
      instance.SetNewArray(&newArray_PHG4Particlev2);
      instance.SetDelete(&delete_PHG4Particlev2);
      instance.SetDeleteArray(&deleteArray_PHG4Particlev2);
      instance.SetDestructor(&destruct_PHG4Particlev2);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4Particlev2*)
   {
      return GenerateInitInstanceLocal((::PHG4Particlev2*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4Particlev2*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4Particlev2::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4Particlev2::Class_Name()
{
   return "PHG4Particlev2";
}

//______________________________________________________________________________
const char *PHG4Particlev2::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Particlev2*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4Particlev2::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Particlev2*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4Particlev2::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Particlev2*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4Particlev2::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Particlev2*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4Particlev2::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4Particlev2.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4Particlev2::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4Particlev2::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4Particlev2(void *p) {
      return  p ? new(p) ::PHG4Particlev2 : new ::PHG4Particlev2;
   }
   static void *newArray_PHG4Particlev2(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4Particlev2[nElements] : new ::PHG4Particlev2[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4Particlev2(void *p) {
      delete ((::PHG4Particlev2*)p);
   }
   static void deleteArray_PHG4Particlev2(void *p) {
      delete [] ((::PHG4Particlev2*)p);
   }
   static void destruct_PHG4Particlev2(void *p) {
      typedef ::PHG4Particlev2 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4Particlev2

namespace {
  void TriggerDictionaryInitialization_PHG4Particlev2_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4Particlev2_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4Particlev2;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4Particlev2_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PARTICLEV2_H
#define G4MAIN_PHG4PARTICLEV2_H

#include "PHG4Particlev1.h"

#include <iostream>
#include <string>

class PHG4Particle;

class PHG4Particlev2 : public PHG4Particlev1
{
 public:
  PHG4Particlev2();
  PHG4Particlev2(const std::string &name, const int pid, const double px, const double py, const double pz);
  PHG4Particlev2(const PHG4Particle *in);

  virtual ~PHG4Particlev2() {}

  int get_track_id() const { return trkid; }
  int get_vtx_id() const { return vtxid; }
  int get_parent_id() const { return parentid; }
  int get_primary_id() const { return primaryid; }
  double get_e() const { return fe; }

  void set_track_id(const int i) { trkid = i; }
  void set_vtx_id(const int i) { vtxid = i; }
  void set_parent_id(const int i) { parentid = i; }
  void set_primary_id(const int i) { primaryid = i; }
  void set_e(const double e) { fe = e; }

  void identify(std::ostream &os = std::cout) const;

 protected:
  int trkid;
  int vtxid;
  int parentid;
  int primaryid;
  double fe;

  ClassDef(PHG4Particlev2, 2)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4Particlev2", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4Particlev2_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4Particlev2_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4Particlev2_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4Particlev2_Dict() {
  TriggerDictionaryInitialization_PHG4Particlev2_Dict_Impl();
}
