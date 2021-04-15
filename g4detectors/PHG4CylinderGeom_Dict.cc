// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4CylinderGeom_Dict
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
#include "PHG4CylinderGeom.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_PHG4CylinderGeom(void *p);
   static void deleteArray_PHG4CylinderGeom(void *p);
   static void destruct_PHG4CylinderGeom(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderGeom*)
   {
      ::PHG4CylinderGeom *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderGeom >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderGeom", ::PHG4CylinderGeom::Class_Version(), "", 20,
                  typeid(::PHG4CylinderGeom), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderGeom::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderGeom) );
      instance.SetDelete(&delete_PHG4CylinderGeom);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderGeom);
      instance.SetDestructor(&destruct_PHG4CylinderGeom);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderGeom*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderGeom*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderGeom*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderGeom::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderGeom::Class_Name()
{
   return "PHG4CylinderGeom";
}

//______________________________________________________________________________
const char *PHG4CylinderGeom::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderGeom::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderGeom::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderGeom::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4CylinderGeom::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderGeom.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderGeom::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderGeom::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PHG4CylinderGeom(void *p) {
      delete ((::PHG4CylinderGeom*)p);
   }
   static void deleteArray_PHG4CylinderGeom(void *p) {
      delete [] ((::PHG4CylinderGeom*)p);
   }
   static void destruct_PHG4CylinderGeom(void *p) {
      typedef ::PHG4CylinderGeom current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderGeom

namespace {
  void TriggerDictionaryInitialization_PHG4CylinderGeom_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4CylinderGeom_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4CylinderGeom;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4CylinderGeom_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERGEOM_H
#define G4DETECTORS_PHG4CYLINDERGEOM_H

#include <phool/PHObject.h>

#include <phool/phool.h>

#include <cmath>
#include <iostream>          // for cout, ostream

class PHParameters;

class PHG4CylinderGeom: public PHObject
{
 public:

  virtual ~PHG4CylinderGeom() {}

  virtual void identify(std::ostream& os = std::cout) const;
  virtual int get_layer() const {PHOOL_VIRTUAL_WARN("get_layer()"); return -99999;}
  virtual double get_radius() const {PHOOL_VIRTUAL_WARN("get_radius()");return NAN;}
  virtual double get_thickness() const {PHOOL_VIRTUAL_WARN("get_thickness()");return NAN;}
  virtual double get_zmin() const {PHOOL_VIRTUAL_WARN("get_zmin()");return NAN;}
  virtual double get_zmax() const {PHOOL_VIRTUAL_WARN("get_zmax()");return NAN;}
  virtual int get_nscint() const {PHOOL_VIRTUAL_WARN("get_nscint()"); return -99999;}
  virtual double get_tiltangle() const {PHOOL_VIRTUAL_WARN("get_tiltangle()");return NAN;}
  virtual double get_phi_slat_zero() const {PHOOL_VIRTUAL_WARN("get_phi_slat_zero()"); return NAN;}

  virtual void set_layer(const int i) {PHOOL_VIRTUAL_WARN("set_layer(const int)");}
  virtual void set_radius(const double r) {PHOOL_VIRTUAL_WARN("set_radius(const double)");}
  virtual void set_thickness(const double t) {PHOOL_VIRTUAL_WARN("set_thickness(const double)");}
  virtual void set_zmin(const double z) {PHOOL_VIRTUAL_WARN("set_zmin(const double)");}
  virtual void set_zmax(const double z) {PHOOL_VIRTUAL_WARN("set_zmax(const double)");}
  virtual void set_nscint(const int i) {PHOOL_VIRTUAL_WARN("set_nscint(const int)"); return;}
  virtual void set_tiltangle(const double i) {PHOOL_VIRTUAL_WARN("set_tiltangle(const double)"); return;}
  virtual void set_phi_slat_zero (const double phi) {PHOOL_VIRTUAL_WARN("set_phi_slat_zero (const double)"); return;}

  virtual void find_segment_center(const int segment_z_bin, const int segment_phi_bin, double location[]){PHOOL_VIRTUAL_WARN("find_sensor_center"); return;}
  virtual void find_strip_center(const int segment_z_bin, const int segment_phi_bin, const int strip_column, const int strip_index, double location[]){PHOOL_VIRTUAL_WARN("find_strip_center"); return;}
 virtual void find_strip_index_values(const int segment_z_bin, const double ypos, const double zpos, int &strip_y_index, int &strip_z_index){PHOOL_VIRTUAL_WARN("find_strip_index_values"); return;}
  virtual void find_strip_center_local_coords(const int segment_z_bin, const int strip_y_index, const int strip_z_index, double location[]){PHOOL_VIRTUAL_WARN("find_strip_center_localcoords"); return;}

  virtual double get_strip_y_spacing() const {PHOOL_VIRTUAL_WARN("get_strip_y_spacing"); return NAN;}
  virtual double get_strip_z_spacing() const {PHOOL_VIRTUAL_WARN("get_strip_z_spacing"); return NAN;}
  virtual double get_strip_tilt() const {PHOOL_VIRTUAL_WARN("get_strip_tilt"); return NAN;}

  virtual int get_N_strip_columns() const {PHOOL_VIRTUAL_WARN("get_N_strip_columns"); return -9999;}
  virtual int get_N_strips_per_column() const {PHOOL_VIRTUAL_WARN("get_N_strips_per_column"); return -9999;}
  virtual int get_N_sensors_in_layer() const {PHOOL_VIRTUAL_WARN("get_N_sensors_in_layer"); return -9999;}

  virtual double get_pixel_z() const {PHOOL_VIRTUAL_WARN("get_pixel_z"); return NAN;}
  virtual double get_pixel_x() const {PHOOL_VIRTUAL_WARN("get_pixel_x"); return NAN;}
  virtual double get_pixel_thickness() const {PHOOL_VIRTUAL_WARN("get_pixel_thickness"); return NAN;}

  //! load parameters from PHParameters, which interface to Database/XML/ROOT files
  virtual void ImportParameters(const PHParameters & param) {return ;}

 protected:
  PHG4CylinderGeom() {}

  ClassDef(PHG4CylinderGeom,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4CylinderGeom", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4CylinderGeom_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4CylinderGeom_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4CylinderGeom_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4CylinderGeom_Dict() {
  TriggerDictionaryInitialization_PHG4CylinderGeom_Dict_Impl();
}
