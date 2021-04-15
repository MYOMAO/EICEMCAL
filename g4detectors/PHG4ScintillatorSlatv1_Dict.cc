// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4ScintillatorSlatv1_Dict
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
#include "PHG4ScintillatorSlatv1.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4ScintillatorSlatv1(void *p = 0);
   static void *newArray_PHG4ScintillatorSlatv1(Long_t size, void *p);
   static void delete_PHG4ScintillatorSlatv1(void *p);
   static void deleteArray_PHG4ScintillatorSlatv1(void *p);
   static void destruct_PHG4ScintillatorSlatv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4ScintillatorSlatv1*)
   {
      ::PHG4ScintillatorSlatv1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4ScintillatorSlatv1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4ScintillatorSlatv1", ::PHG4ScintillatorSlatv1::Class_Version(), "", 21,
                  typeid(::PHG4ScintillatorSlatv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4ScintillatorSlatv1::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4ScintillatorSlatv1) );
      instance.SetNew(&new_PHG4ScintillatorSlatv1);
      instance.SetNewArray(&newArray_PHG4ScintillatorSlatv1);
      instance.SetDelete(&delete_PHG4ScintillatorSlatv1);
      instance.SetDeleteArray(&deleteArray_PHG4ScintillatorSlatv1);
      instance.SetDestructor(&destruct_PHG4ScintillatorSlatv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4ScintillatorSlatv1*)
   {
      return GenerateInitInstanceLocal((::PHG4ScintillatorSlatv1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4ScintillatorSlatv1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4ScintillatorSlatv1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4ScintillatorSlatv1::Class_Name()
{
   return "PHG4ScintillatorSlatv1";
}

//______________________________________________________________________________
const char *PHG4ScintillatorSlatv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4ScintillatorSlatv1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4ScintillatorSlatv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4ScintillatorSlatv1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4ScintillatorSlatv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4ScintillatorSlatv1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4ScintillatorSlatv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4ScintillatorSlatv1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4ScintillatorSlatv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4ScintillatorSlatv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4ScintillatorSlatv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4ScintillatorSlatv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4ScintillatorSlatv1(void *p) {
      return  p ? new(p) ::PHG4ScintillatorSlatv1 : new ::PHG4ScintillatorSlatv1;
   }
   static void *newArray_PHG4ScintillatorSlatv1(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4ScintillatorSlatv1[nElements] : new ::PHG4ScintillatorSlatv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4ScintillatorSlatv1(void *p) {
      delete ((::PHG4ScintillatorSlatv1*)p);
   }
   static void deleteArray_PHG4ScintillatorSlatv1(void *p) {
      delete [] ((::PHG4ScintillatorSlatv1*)p);
   }
   static void destruct_PHG4ScintillatorSlatv1(void *p) {
      typedef ::PHG4ScintillatorSlatv1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4ScintillatorSlatv1

namespace ROOT {
   static TClass *setlEULong64_tgR_Dictionary();
   static void setlEULong64_tgR_TClassManip(TClass*);
   static void *new_setlEULong64_tgR(void *p = 0);
   static void *newArray_setlEULong64_tgR(Long_t size, void *p);
   static void delete_setlEULong64_tgR(void *p);
   static void deleteArray_setlEULong64_tgR(void *p);
   static void destruct_setlEULong64_tgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const set<ULong64_t>*)
   {
      set<ULong64_t> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(set<ULong64_t>));
      static ::ROOT::TGenericClassInfo 
         instance("set<ULong64_t>", -2, "set", 94,
                  typeid(set<ULong64_t>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &setlEULong64_tgR_Dictionary, isa_proxy, 0,
                  sizeof(set<ULong64_t>) );
      instance.SetNew(&new_setlEULong64_tgR);
      instance.SetNewArray(&newArray_setlEULong64_tgR);
      instance.SetDelete(&delete_setlEULong64_tgR);
      instance.SetDeleteArray(&deleteArray_setlEULong64_tgR);
      instance.SetDestructor(&destruct_setlEULong64_tgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Insert< set<ULong64_t> >()));

      ::ROOT::AddClassAlternate("set<ULong64_t>","std::set<unsigned long long, std::less<unsigned long long>, std::allocator<unsigned long long> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const set<ULong64_t>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *setlEULong64_tgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const set<ULong64_t>*)0x0)->GetClass();
      setlEULong64_tgR_TClassManip(theClass);
   return theClass;
   }

   static void setlEULong64_tgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_setlEULong64_tgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<ULong64_t> : new set<ULong64_t>;
   }
   static void *newArray_setlEULong64_tgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<ULong64_t>[nElements] : new set<ULong64_t>[nElements];
   }
   // Wrapper around operator delete
   static void delete_setlEULong64_tgR(void *p) {
      delete ((set<ULong64_t>*)p);
   }
   static void deleteArray_setlEULong64_tgR(void *p) {
      delete [] ((set<ULong64_t>*)p);
   }
   static void destruct_setlEULong64_tgR(void *p) {
      typedef set<ULong64_t> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class set<ULong64_t>

namespace {
  void TriggerDictionaryInitialization_PHG4ScintillatorSlatv1_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4ScintillatorSlatv1_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4ScintillatorSlatv1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4ScintillatorSlatv1_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4SCINTILLATORSLATV1_H
#define G4DETECTORS_PHG4SCINTILLATORSLATV1_H

#include "PHG4ScintillatorSlat.h"

#include "PHG4ScintillatorSlatDefs.h"  // for keytype

#include <g4main/PHG4HitDefs.h>

#include <iostream>                    // for cout, ostream
#include <set>
#include <utility>                     // for make_pair, pair

class PHG4ScintillatorSlatv1 : public PHG4ScintillatorSlat
{
 public:

  PHG4ScintillatorSlatv1();
  virtual ~PHG4ScintillatorSlatv1(){}

  void identify(std::ostream& os = std::cout) const;

  void add_edep(const double f, const double e, const double ly) {edep+=f; eion+= e; light_yield+=ly;}
  void add_hit_key(PHG4HitDefs::keytype key) {hit_id.insert(key);}
  
  void set_key(PHG4ScintillatorSlatDefs::keytype i) {key = i;}
  void set_light_yield(const double lightYield)  {light_yield = lightYield;}

  short get_row() const;
  short get_column() const;
  PHG4ScintillatorSlatDefs::keytype get_key() const {return key;}
  double get_edep() const {return edep;}
  double get_eion() const {return eion;}
  double get_light_yield() const {return light_yield;}
  std::pair<std::set<PHG4HitDefs::keytype>::const_iterator, std::set<PHG4HitDefs::keytype>::const_iterator> get_hit_ids() const {return std::make_pair(hit_id.begin(),hit_id.end());}


 protected:
  PHG4ScintillatorSlatDefs::keytype key;
  double edep;
  double eion;
  double light_yield;
  std::set<PHG4HitDefs::keytype> hit_id;

   
  ClassDef(PHG4ScintillatorSlatv1,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4ScintillatorSlatv1", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4ScintillatorSlatv1_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4ScintillatorSlatv1_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4ScintillatorSlatv1_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4ScintillatorSlatv1_Dict() {
  TriggerDictionaryInitialization_PHG4ScintillatorSlatv1_Dict_Impl();
}
