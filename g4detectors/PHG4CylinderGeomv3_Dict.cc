// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4CylinderGeomv3_Dict
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
#include "PHG4CylinderGeomv3.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4CylinderGeomv3(void *p = 0);
   static void *newArray_PHG4CylinderGeomv3(Long_t size, void *p);
   static void delete_PHG4CylinderGeomv3(void *p);
   static void deleteArray_PHG4CylinderGeomv3(void *p);
   static void destruct_PHG4CylinderGeomv3(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderGeomv3*)
   {
      ::PHG4CylinderGeomv3 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderGeomv3 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderGeomv3", ::PHG4CylinderGeomv3::Class_Version(), "", 15,
                  typeid(::PHG4CylinderGeomv3), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderGeomv3::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderGeomv3) );
      instance.SetNew(&new_PHG4CylinderGeomv3);
      instance.SetNewArray(&newArray_PHG4CylinderGeomv3);
      instance.SetDelete(&delete_PHG4CylinderGeomv3);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderGeomv3);
      instance.SetDestructor(&destruct_PHG4CylinderGeomv3);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderGeomv3*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderGeomv3*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderGeomv3*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderGeomv3::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderGeomv3::Class_Name()
{
   return "PHG4CylinderGeomv3";
}

//______________________________________________________________________________
const char *PHG4CylinderGeomv3::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeomv3*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderGeomv3::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeomv3*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderGeomv3::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeomv3*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderGeomv3::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeomv3*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4CylinderGeomv3::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderGeomv3.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderGeomv3::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderGeomv3::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4CylinderGeomv3(void *p) {
      return  p ? new(p) ::PHG4CylinderGeomv3 : new ::PHG4CylinderGeomv3;
   }
   static void *newArray_PHG4CylinderGeomv3(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4CylinderGeomv3[nElements] : new ::PHG4CylinderGeomv3[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4CylinderGeomv3(void *p) {
      delete ((::PHG4CylinderGeomv3*)p);
   }
   static void deleteArray_PHG4CylinderGeomv3(void *p) {
      delete [] ((::PHG4CylinderGeomv3*)p);
   }
   static void destruct_PHG4CylinderGeomv3(void *p) {
      typedef ::PHG4CylinderGeomv3 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderGeomv3

namespace {
  void TriggerDictionaryInitialization_PHG4CylinderGeomv3_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4CylinderGeomv3_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4CylinderGeomv3;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4CylinderGeomv3_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERGEOMV3_H
#define G4DETECTORS_PHG4CYLINDERGEOMV3_H

#include "PHG4CylinderGeomv2.h"

#include <iostream>              // for cout, ostream

class PHG4CylinderGeomv3: public PHG4CylinderGeomv2
{
 public:
  PHG4CylinderGeomv3();
  PHG4CylinderGeomv3(const double r, const double zmi, const double zma, const double thickn, const int n_scint,
                     const double tangl, const double phi_slat_null):
    PHG4CylinderGeomv2(r,zmi,zma,thickn,n_scint),
    tiltangle(tangl),
    phi_slat_zero(phi_slat_null)
      {}

  virtual ~PHG4CylinderGeomv3() {}

  void identify(std::ostream& os = std::cout) const;
  void set_tiltangle (const double phi) {tiltangle=phi;}
  void set_phi_slat_zero (const double phi) {phi_slat_zero=phi;}

  double get_phi_slat_zero() const {return phi_slat_zero;}
  double get_tiltangle() const {return tiltangle;}

 protected:
  double tiltangle;
  double phi_slat_zero;

  ClassDef(PHG4CylinderGeomv3,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4CylinderGeomv3", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4CylinderGeomv3_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4CylinderGeomv3_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4CylinderGeomv3_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4CylinderGeomv3_Dict() {
  TriggerDictionaryInitialization_PHG4CylinderGeomv3_Dict_Impl();
}
