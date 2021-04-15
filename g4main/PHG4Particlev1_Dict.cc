// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4Particlev1_Dict
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
#include "PHG4Particlev1.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4Particlev1(void *p = 0);
   static void *newArray_PHG4Particlev1(Long_t size, void *p);
   static void delete_PHG4Particlev1(void *p);
   static void deleteArray_PHG4Particlev1(void *p);
   static void destruct_PHG4Particlev1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4Particlev1*)
   {
      ::PHG4Particlev1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4Particlev1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4Particlev1", ::PHG4Particlev1::Class_Version(), "", 13,
                  typeid(::PHG4Particlev1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4Particlev1::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4Particlev1) );
      instance.SetNew(&new_PHG4Particlev1);
      instance.SetNewArray(&newArray_PHG4Particlev1);
      instance.SetDelete(&delete_PHG4Particlev1);
      instance.SetDeleteArray(&deleteArray_PHG4Particlev1);
      instance.SetDestructor(&destruct_PHG4Particlev1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4Particlev1*)
   {
      return GenerateInitInstanceLocal((::PHG4Particlev1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4Particlev1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4Particlev1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4Particlev1::Class_Name()
{
   return "PHG4Particlev1";
}

//______________________________________________________________________________
const char *PHG4Particlev1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Particlev1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4Particlev1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Particlev1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4Particlev1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Particlev1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4Particlev1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Particlev1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4Particlev1::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4Particlev1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4Particlev1::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4Particlev1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4Particlev1(void *p) {
      return  p ? new(p) ::PHG4Particlev1 : new ::PHG4Particlev1;
   }
   static void *newArray_PHG4Particlev1(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4Particlev1[nElements] : new ::PHG4Particlev1[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4Particlev1(void *p) {
      delete ((::PHG4Particlev1*)p);
   }
   static void deleteArray_PHG4Particlev1(void *p) {
      delete [] ((::PHG4Particlev1*)p);
   }
   static void destruct_PHG4Particlev1(void *p) {
      typedef ::PHG4Particlev1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4Particlev1

namespace {
  void TriggerDictionaryInitialization_PHG4Particlev1_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4Particlev1_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4Particlev1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4Particlev1_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PARTICLEV1_H
#define G4MAIN_PHG4PARTICLEV1_H

#include "PHG4Particle.h"

#include <iostream>
#include <string>

class PHG4Particlev1 : public PHG4Particle
{
 public:
  PHG4Particlev1();
  PHG4Particlev1(const std::string &name, const int pid, const double px, const double py, const double pz);
  PHG4Particlev1(const PHG4Particle *in);

  virtual ~PHG4Particlev1() {}

  int get_pid() const { return fpid; }
  std::string get_name() const { return fname; }
  double get_px() const { return fpx; }
  double get_py() const { return fpy; }
  double get_pz() const { return fpz; }

  int get_barcode() const { return barcode; }

  void set_name(const std::string &name) { fname = name; }
  void set_pid(const int i) { fpid = i; }
  void set_px(const double x) { fpx = x; }
  void set_py(const double x) { fpy = x; }
  void set_pz(const double x) { fpz = x; }

  void set_barcode(const int bcd) { barcode = bcd; }

  void identify(std::ostream &os = std::cout) const;

 protected:
  std::string fname;
  int fpid;
  double fpx, fpy, fpz;
  int barcode;

  ClassDef(PHG4Particlev1, 1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4Particlev1", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4Particlev1_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4Particlev1_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4Particlev1_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4Particlev1_Dict() {
  TriggerDictionaryInitialization_PHG4Particlev1_Dict_Impl();
}
