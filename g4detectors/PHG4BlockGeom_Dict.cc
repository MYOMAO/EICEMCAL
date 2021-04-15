// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4BlockGeom_Dict
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
#include "PHG4BlockGeom.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_PHG4BlockGeom(void *p);
   static void deleteArray_PHG4BlockGeom(void *p);
   static void destruct_PHG4BlockGeom(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4BlockGeom*)
   {
      ::PHG4BlockGeom *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4BlockGeom >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4BlockGeom", ::PHG4BlockGeom::Class_Version(), "", 18,
                  typeid(::PHG4BlockGeom), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4BlockGeom::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4BlockGeom) );
      instance.SetDelete(&delete_PHG4BlockGeom);
      instance.SetDeleteArray(&deleteArray_PHG4BlockGeom);
      instance.SetDestructor(&destruct_PHG4BlockGeom);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4BlockGeom*)
   {
      return GenerateInitInstanceLocal((::PHG4BlockGeom*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4BlockGeom*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4BlockGeom::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4BlockGeom::Class_Name()
{
   return "PHG4BlockGeom";
}

//______________________________________________________________________________
const char *PHG4BlockGeom::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockGeom*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4BlockGeom::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockGeom*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4BlockGeom::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockGeom*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4BlockGeom::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockGeom*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4BlockGeom::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4BlockGeom.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4BlockGeom::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4BlockGeom::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PHG4BlockGeom(void *p) {
      delete ((::PHG4BlockGeom*)p);
   }
   static void deleteArray_PHG4BlockGeom(void *p) {
      delete [] ((::PHG4BlockGeom*)p);
   }
   static void destruct_PHG4BlockGeom(void *p) {
      typedef ::PHG4BlockGeom current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4BlockGeom

namespace {
  void TriggerDictionaryInitialization_PHG4BlockGeom_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4BlockGeom_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4BlockGeom;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4BlockGeom_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BLOCKGEOM_H
#define G4DETECTORS_PHG4BLOCKGEOM_H

#include <phool/PHObject.h>

#include <phool/phool.h>

#include <iostream>          // for cout, ostream
#include <cmath>

class PHG4BlockGeom: public PHObject
{
 public:

  virtual ~PHG4BlockGeom() {}

  virtual void identify(std::ostream& os = std::cout) const;

  virtual int get_layer() const {PHOOL_VIRTUAL_WARN("get_layer()"); return -99999;}
  virtual double get_size_x() const {PHOOL_VIRTUAL_WARN("get_size_x()");return NAN;}
  virtual double get_size_y() const {PHOOL_VIRTUAL_WARN("get_size_y()");return NAN;}
  virtual double get_size_z() const {PHOOL_VIRTUAL_WARN("get_size_z()");return NAN;}
  virtual double get_center_x() const {PHOOL_VIRTUAL_WARN("get_place_x()");return NAN;}
  virtual double get_center_y() const {PHOOL_VIRTUAL_WARN("get_place_y()");return NAN;}
  virtual double get_center_z() const {PHOOL_VIRTUAL_WARN("get_place_z()");return NAN;}
  virtual double get_z_rot() const {PHOOL_VIRTUAL_WARN("get_z_rot()");return NAN;}

  virtual double get_width() const {PHOOL_VIRTUAL_WARN("get_width()");return NAN;}
  virtual double get_length() const {PHOOL_VIRTUAL_WARN("get_length()");return NAN;}
  virtual double get_thickness() const {PHOOL_VIRTUAL_WARN("get_thickness()");return NAN;}

  virtual double get_rot_matrix(const int, const int) const {PHOOL_VIRTUAL_WARN("get_rot_matrix(const int, const int)"); return NAN;}

  virtual void set_layer(const int i) {PHOOL_VIRTUAL_WARN("set_layer(const int)");}
  virtual void set_size(const double sizex, const double sizey, const double sizez)
    {PHOOL_VIRTUAL_WARN("set_size(const double, const double, const double)");}
  virtual void set_place(const double placex, const double placey, const double placez)
    {PHOOL_VIRTUAL_WARN("set_place(const double, const double, const double)");}
  virtual void set_z_rot(const double z) {PHOOL_VIRTUAL_WARN("set_z_rot(const double)");}

  virtual void convert_local_to_global(const double, const double, const double,
                                       double &, double &, double &) const
    {PHOOL_VIRTUAL_WARN("convert_local_to_global(const double, const double, const double, double &, double &, double &)");}
  virtual void convert_global_to_local(const double, const double, const double, 
                                         double &, double &, double &) const
    {PHOOL_VIRTUAL_WARN("convert_global_to_local(const double, const double, const double, double &, double &, double &)");}

 protected:
  PHG4BlockGeom() {}

  ClassDef(PHG4BlockGeom,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4BlockGeom", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4BlockGeom_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4BlockGeom_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4BlockGeom_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4BlockGeom_Dict() {
  TriggerDictionaryInitialization_PHG4BlockGeom_Dict_Impl();
}
