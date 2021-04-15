// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4VtxPointv1_Dict
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
#include "PHG4VtxPointv1.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4VtxPointv1(void *p = 0);
   static void *newArray_PHG4VtxPointv1(Long_t size, void *p);
   static void delete_PHG4VtxPointv1(void *p);
   static void deleteArray_PHG4VtxPointv1(void *p);
   static void destruct_PHG4VtxPointv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4VtxPointv1*)
   {
      ::PHG4VtxPointv1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4VtxPointv1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4VtxPointv1", ::PHG4VtxPointv1::Class_Version(), "", 14,
                  typeid(::PHG4VtxPointv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4VtxPointv1::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4VtxPointv1) );
      instance.SetNew(&new_PHG4VtxPointv1);
      instance.SetNewArray(&newArray_PHG4VtxPointv1);
      instance.SetDelete(&delete_PHG4VtxPointv1);
      instance.SetDeleteArray(&deleteArray_PHG4VtxPointv1);
      instance.SetDestructor(&destruct_PHG4VtxPointv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4VtxPointv1*)
   {
      return GenerateInitInstanceLocal((::PHG4VtxPointv1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4VtxPointv1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4VtxPointv1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4VtxPointv1::Class_Name()
{
   return "PHG4VtxPointv1";
}

//______________________________________________________________________________
const char *PHG4VtxPointv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4VtxPointv1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4VtxPointv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4VtxPointv1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4VtxPointv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4VtxPointv1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4VtxPointv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4VtxPointv1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4VtxPointv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4VtxPointv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4VtxPointv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4VtxPointv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4VtxPointv1(void *p) {
      return  p ? new(p) ::PHG4VtxPointv1 : new ::PHG4VtxPointv1;
   }
   static void *newArray_PHG4VtxPointv1(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4VtxPointv1[nElements] : new ::PHG4VtxPointv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4VtxPointv1(void *p) {
      delete ((::PHG4VtxPointv1*)p);
   }
   static void deleteArray_PHG4VtxPointv1(void *p) {
      delete [] ((::PHG4VtxPointv1*)p);
   }
   static void destruct_PHG4VtxPointv1(void *p) {
      typedef ::PHG4VtxPointv1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4VtxPointv1

namespace {
  void TriggerDictionaryInitialization_PHG4VtxPointv1_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4VtxPointv1_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4VtxPointv1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4VtxPointv1_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4VTXPOINTV1_H
#define G4MAIN_PHG4VTXPOINTV1_H

#include "PHG4VtxPoint.h"

#include <climits> // for INT_MIN
#include <cmath> // def of NAN
#include <iostream>   // for cout, ostream

class PHG4VtxPointv1: public PHG4VtxPoint
{
 public:
  PHG4VtxPointv1():
    vx(NAN),
    vy(NAN),
    vz(NAN),
    t0(NAN),
    id(INT_MIN)
      {}

  PHG4VtxPointv1(const PHG4VtxPoint *vtx):
    vx(vtx->get_x()),
    vy(vtx->get_y()),
    vz(vtx->get_z()),
    t0(vtx->get_t()),
    id(vtx->get_id())
      {}

  PHG4VtxPointv1(const double x, const double y, const double z, const double t, const int id_value = INT_MIN):
    vx(x),
    vy(y),
    vz(z),
    t0(t),
    id(id_value)
      {}

  virtual ~PHG4VtxPointv1() {}

  void set_x(const double r) {vx = r;}
  void set_y(const double r) {vy = r;}
  void set_z(const double r) {vz = r;}
  void set_t(const double r) {t0 = r;}
  void set_id(const int i) {id = i;}

  double get_x() const {return vx;}
  double get_y() const {return vy;}
  double get_z() const {return vz;}
  double get_t() const {return t0;}
  int get_id() const {return id;}

  void identify(std::ostream& os = std::cout) const;

 protected:

  double vx;
  double vy;
  double vz;
  double t0;

  //! id tag for this vertex
  int id;

  ClassDef(PHG4VtxPointv1,2)
};




#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4VtxPointv1", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4VtxPointv1_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4VtxPointv1_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4VtxPointv1_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4VtxPointv1_Dict() {
  TriggerDictionaryInitialization_PHG4VtxPointv1_Dict_Impl();
}
