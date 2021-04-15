// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4CylinderGeomv1_Dict
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
#include "PHG4CylinderGeomv1.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4CylinderGeomv1(void *p = 0);
   static void *newArray_PHG4CylinderGeomv1(Long_t size, void *p);
   static void delete_PHG4CylinderGeomv1(void *p);
   static void deleteArray_PHG4CylinderGeomv1(void *p);
   static void destruct_PHG4CylinderGeomv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderGeomv1*)
   {
      ::PHG4CylinderGeomv1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderGeomv1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderGeomv1", ::PHG4CylinderGeomv1::Class_Version(), "", 17,
                  typeid(::PHG4CylinderGeomv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderGeomv1::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderGeomv1) );
      instance.SetNew(&new_PHG4CylinderGeomv1);
      instance.SetNewArray(&newArray_PHG4CylinderGeomv1);
      instance.SetDelete(&delete_PHG4CylinderGeomv1);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderGeomv1);
      instance.SetDestructor(&destruct_PHG4CylinderGeomv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderGeomv1*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderGeomv1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderGeomv1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderGeomv1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderGeomv1::Class_Name()
{
   return "PHG4CylinderGeomv1";
}

//______________________________________________________________________________
const char *PHG4CylinderGeomv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeomv1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderGeomv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeomv1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderGeomv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeomv1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderGeomv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeomv1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4CylinderGeomv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderGeomv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderGeomv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderGeomv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4CylinderGeomv1(void *p) {
      return  p ? new(p) ::PHG4CylinderGeomv1 : new ::PHG4CylinderGeomv1;
   }
   static void *newArray_PHG4CylinderGeomv1(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4CylinderGeomv1[nElements] : new ::PHG4CylinderGeomv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4CylinderGeomv1(void *p) {
      delete ((::PHG4CylinderGeomv1*)p);
   }
   static void deleteArray_PHG4CylinderGeomv1(void *p) {
      delete [] ((::PHG4CylinderGeomv1*)p);
   }
   static void destruct_PHG4CylinderGeomv1(void *p) {
      typedef ::PHG4CylinderGeomv1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderGeomv1

namespace {
  void TriggerDictionaryInitialization_PHG4CylinderGeomv1_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4CylinderGeomv1_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4CylinderGeomv1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4CylinderGeomv1_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERGEOMV1_H
#define G4DETECTORS_PHG4CYLINDERGEOMV1_H

#include "PHG4CylinderGeom.h"

#include <iostream>  // for cout, ostream

class PHParameters;

class PHG4CylinderGeomv1 : public PHG4CylinderGeom
{
 public:
  PHG4CylinderGeomv1();
  PHG4CylinderGeomv1(const double r, const double zmi, const double zma, const double thickn)
    : layer(-1)
    , radius(r)
    , zmin(zmi)
    , zmax(zma)
    , thickness(thickn)
  {
  }

  virtual ~PHG4CylinderGeomv1() {}

  void identify(std::ostream& os = std::cout) const;
  int get_layer() const { return layer; }
  double get_radius() const { return radius; }
  double get_thickness() const { return thickness; }
  double get_zmin() const { return zmin; }
  double get_zmax() const { return zmax; }

  void set_layer(const int i) { layer = i; }
  void set_radius(const double r) { radius = r; }
  void set_thickness(const double t) { thickness = t; }
  void set_zmin(const double z) { zmin = z; }
  void set_zmax(const double z) { zmax = z; }

  //! load parameters from PHParameters, which interface to Database/XML/ROOT files
  virtual void ImportParameters(const PHParameters& param);

 protected:
  int layer;
  double radius;
  double zmin;
  double zmax;
  double thickness;

  ClassDef(PHG4CylinderGeomv1, 1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4CylinderGeomv1", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4CylinderGeomv1_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4CylinderGeomv1_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4CylinderGeomv1_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4CylinderGeomv1_Dict() {
  TriggerDictionaryInitialization_PHG4CylinderGeomv1_Dict_Impl();
}
