// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4EventHeaderv1_Dict
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
#include "PHG4EventHeaderv1.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4EventHeaderv1(void *p = 0);
   static void *newArray_PHG4EventHeaderv1(Long_t size, void *p);
   static void delete_PHG4EventHeaderv1(void *p);
   static void deleteArray_PHG4EventHeaderv1(void *p);
   static void destruct_PHG4EventHeaderv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4EventHeaderv1*)
   {
      ::PHG4EventHeaderv1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4EventHeaderv1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4EventHeaderv1", ::PHG4EventHeaderv1::Class_Version(), "", 13,
                  typeid(::PHG4EventHeaderv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4EventHeaderv1::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4EventHeaderv1) );
      instance.SetNew(&new_PHG4EventHeaderv1);
      instance.SetNewArray(&newArray_PHG4EventHeaderv1);
      instance.SetDelete(&delete_PHG4EventHeaderv1);
      instance.SetDeleteArray(&deleteArray_PHG4EventHeaderv1);
      instance.SetDestructor(&destruct_PHG4EventHeaderv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4EventHeaderv1*)
   {
      return GenerateInitInstanceLocal((::PHG4EventHeaderv1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4EventHeaderv1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4EventHeaderv1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4EventHeaderv1::Class_Name()
{
   return "PHG4EventHeaderv1";
}

//______________________________________________________________________________
const char *PHG4EventHeaderv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4EventHeaderv1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4EventHeaderv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4EventHeaderv1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4EventHeaderv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4EventHeaderv1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4EventHeaderv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4EventHeaderv1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4EventHeaderv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4EventHeaderv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4EventHeaderv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4EventHeaderv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4EventHeaderv1(void *p) {
      return  p ? new(p) ::PHG4EventHeaderv1 : new ::PHG4EventHeaderv1;
   }
   static void *newArray_PHG4EventHeaderv1(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4EventHeaderv1[nElements] : new ::PHG4EventHeaderv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4EventHeaderv1(void *p) {
      delete ((::PHG4EventHeaderv1*)p);
   }
   static void deleteArray_PHG4EventHeaderv1(void *p) {
      delete [] ((::PHG4EventHeaderv1*)p);
   }
   static void destruct_PHG4EventHeaderv1(void *p) {
      typedef ::PHG4EventHeaderv1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4EventHeaderv1

namespace {
  void TriggerDictionaryInitialization_PHG4EventHeaderv1_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4EventHeaderv1_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4EventHeaderv1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4EventHeaderv1_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4EVENTHEADERV1_H
#define G4MAIN_PHG4EVENTHEADERV1_H

#include "PHG4EventHeader.h"

#include <iostream>

///
class PHG4EventHeaderv1: public PHG4EventHeader
{
 public:

  PHG4EventHeaderv1();

  /// dtor
  virtual ~PHG4EventHeaderv1() {}

  /// Clear Event
  virtual void Reset();

  /** identify Function from PHObject
      @param os Output Stream 
   */
  void identify(std::ostream& os = std::cout) const;

  /// isValid returns non zero if object contains valid data
  int isValid() const;

  /// get Event Number
  int get_EvtSequence() const {return evtseq;}
  /// set Event Number
  void set_EvtSequence(const int ival) {evtseq = ival;}

  float get_ImpactParameter() const {return bimp;}
  void set_ImpactParameter(const float b) {bimp = b;}

  float get_EventPlaneAngle() const {return rplane;}
  void set_EventPlaneAngle(const float r) {rplane = r;}

 protected:
  int evtseq;
  float bimp;
  float rplane;

 private: // prevent doc++ from showing ClassDef
  ClassDef(PHG4EventHeaderv1,1)

};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4EventHeaderv1", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4EventHeaderv1_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4EventHeaderv1_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4EventHeaderv1_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4EventHeaderv1_Dict() {
  TriggerDictionaryInitialization_PHG4EventHeaderv1_Dict_Impl();
}
