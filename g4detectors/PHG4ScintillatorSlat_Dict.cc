// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4ScintillatorSlat_Dict
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
#include "PHG4ScintillatorSlat.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_PHG4ScintillatorSlat(void *p);
   static void deleteArray_PHG4ScintillatorSlat(void *p);
   static void destruct_PHG4ScintillatorSlat(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4ScintillatorSlat*)
   {
      ::PHG4ScintillatorSlat *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4ScintillatorSlat >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4ScintillatorSlat", ::PHG4ScintillatorSlat::Class_Version(), "", 19,
                  typeid(::PHG4ScintillatorSlat), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4ScintillatorSlat::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4ScintillatorSlat) );
      instance.SetDelete(&delete_PHG4ScintillatorSlat);
      instance.SetDeleteArray(&deleteArray_PHG4ScintillatorSlat);
      instance.SetDestructor(&destruct_PHG4ScintillatorSlat);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4ScintillatorSlat*)
   {
      return GenerateInitInstanceLocal((::PHG4ScintillatorSlat*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4ScintillatorSlat*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4ScintillatorSlat::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4ScintillatorSlat::Class_Name()
{
   return "PHG4ScintillatorSlat";
}

//______________________________________________________________________________
const char *PHG4ScintillatorSlat::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4ScintillatorSlat*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4ScintillatorSlat::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4ScintillatorSlat*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4ScintillatorSlat::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4ScintillatorSlat*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4ScintillatorSlat::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4ScintillatorSlat*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4ScintillatorSlat::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4ScintillatorSlat.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4ScintillatorSlat::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4ScintillatorSlat::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PHG4ScintillatorSlat(void *p) {
      delete ((::PHG4ScintillatorSlat*)p);
   }
   static void deleteArray_PHG4ScintillatorSlat(void *p) {
      delete [] ((::PHG4ScintillatorSlat*)p);
   }
   static void destruct_PHG4ScintillatorSlat(void *p) {
      typedef ::PHG4ScintillatorSlat current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4ScintillatorSlat

namespace {
  void TriggerDictionaryInitialization_PHG4ScintillatorSlat_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4ScintillatorSlat_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4ScintillatorSlat;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4ScintillatorSlat_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4SCINTILLATORSLAT_H
#define G4DETECTORS_PHG4SCINTILLATORSLAT_H

#include "PHG4ScintillatorSlatDefs.h"

#include <g4main/PHG4HitDefs.h>
#include <phool/PHObject.h>
#include <cmath>
#include <map>
#include <set>

class PHG4ScintillatorSlat : public PHObject
{
 public:
  
  virtual ~PHG4ScintillatorSlat(){}

  virtual void identify(std::ostream& os = std::cout) const {
    os << "PHG4ScintillatorSlat base class" << std::endl;
  }
  
  virtual void add_edep(const double edep, const double e, const double light_yield) {return;}

  virtual void set_key(const PHG4ScintillatorSlatDefs::keytype) {return;}
  virtual void add_hit_key(PHG4HitDefs::keytype) {return;}

  virtual short get_column() const {return -1;}
  virtual short get_row() const {return -1;}
  virtual PHG4ScintillatorSlatDefs::keytype get_key() const {return 0xFFFFFFFF;}

  virtual double get_edep() const {return NAN;}
  virtual double get_eion() const {return NAN;}
  virtual double get_light_yield() const  {return NAN;}
  virtual std::pair<std::set<PHG4HitDefs::keytype>::const_iterator, std::set<PHG4HitDefs::keytype>::const_iterator> get_hit_ids() const = 0;

  
 protected:

  PHG4ScintillatorSlat() {}
  ClassDef(PHG4ScintillatorSlat,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4ScintillatorSlat", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4ScintillatorSlat_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4ScintillatorSlat_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4ScintillatorSlat_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4ScintillatorSlat_Dict() {
  TriggerDictionaryInitialization_PHG4ScintillatorSlat_Dict_Impl();
}
