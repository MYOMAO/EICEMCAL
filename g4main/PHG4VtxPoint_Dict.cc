// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4VtxPoint_Dict
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
#include "PHG4VtxPoint.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_PHG4VtxPoint(void *p);
   static void deleteArray_PHG4VtxPoint(void *p);
   static void destruct_PHG4VtxPoint(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4VtxPoint*)
   {
      ::PHG4VtxPoint *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4VtxPoint >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4VtxPoint", ::PHG4VtxPoint::Class_Version(), "", 14,
                  typeid(::PHG4VtxPoint), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4VtxPoint::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4VtxPoint) );
      instance.SetDelete(&delete_PHG4VtxPoint);
      instance.SetDeleteArray(&deleteArray_PHG4VtxPoint);
      instance.SetDestructor(&destruct_PHG4VtxPoint);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4VtxPoint*)
   {
      return GenerateInitInstanceLocal((::PHG4VtxPoint*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4VtxPoint*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4VtxPoint::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4VtxPoint::Class_Name()
{
   return "PHG4VtxPoint";
}

//______________________________________________________________________________
const char *PHG4VtxPoint::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4VtxPoint*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4VtxPoint::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4VtxPoint*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4VtxPoint::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4VtxPoint*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4VtxPoint::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4VtxPoint*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4VtxPoint::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4VtxPoint.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4VtxPoint::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4VtxPoint::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PHG4VtxPoint(void *p) {
      delete ((::PHG4VtxPoint*)p);
   }
   static void deleteArray_PHG4VtxPoint(void *p) {
      delete [] ((::PHG4VtxPoint*)p);
   }
   static void destruct_PHG4VtxPoint(void *p) {
      typedef ::PHG4VtxPoint current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4VtxPoint

namespace {
  void TriggerDictionaryInitialization_PHG4VtxPoint_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4VtxPoint_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4VtxPoint;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4VtxPoint_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4VTXPOINT_H
#define G4MAIN_PHG4VTXPOINT_H

#include <phool/PHObject.h>

#include <cmath>
#include <climits>
#include <iostream>

class PHG4VtxPoint: public PHObject
{
 public:
  virtual ~PHG4VtxPoint() {}

  virtual void set_x(const double r) {}
  virtual void set_y(const double r) {}
  virtual void set_z(const double r) {}
  virtual void set_t(const double r) {}
  virtual void set_id(const int i) {}

  virtual double get_x() const {return NAN;}
  virtual double get_y() const {return NAN;}
  virtual double get_z() const {return NAN;}
  virtual double get_t() const {return NAN;}
  virtual int get_id() const {return INT_MIN;}

  virtual void identify(std::ostream& os = std::cout) const;

  //! comparison of vertex value only, not on the id, per algorithm requirement in PHG4TruthInfoContainer::AddPrimaryVertex
  bool  operator== (const PHG4VtxPoint &) const ;


 protected:
  PHG4VtxPoint(){}
  ClassDef(PHG4VtxPoint,1)

};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4VtxPoint", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4VtxPoint_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4VtxPoint_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4VtxPoint_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4VtxPoint_Dict() {
  TriggerDictionaryInitialization_PHG4VtxPoint_Dict_Impl();
}
