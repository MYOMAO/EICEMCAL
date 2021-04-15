// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4CylinderGeomv4_Dict
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
#include "PHG4CylinderGeomv4.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4CylinderGeomv4(void *p = 0);
   static void *newArray_PHG4CylinderGeomv4(Long_t size, void *p);
   static void delete_PHG4CylinderGeomv4(void *p);
   static void deleteArray_PHG4CylinderGeomv4(void *p);
   static void destruct_PHG4CylinderGeomv4(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderGeomv4*)
   {
      ::PHG4CylinderGeomv4 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderGeomv4 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderGeomv4", ::PHG4CylinderGeomv4::Class_Version(), "", 15,
                  typeid(::PHG4CylinderGeomv4), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderGeomv4::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderGeomv4) );
      instance.SetNew(&new_PHG4CylinderGeomv4);
      instance.SetNewArray(&newArray_PHG4CylinderGeomv4);
      instance.SetDelete(&delete_PHG4CylinderGeomv4);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderGeomv4);
      instance.SetDestructor(&destruct_PHG4CylinderGeomv4);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderGeomv4*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderGeomv4*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderGeomv4*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderGeomv4::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderGeomv4::Class_Name()
{
   return "PHG4CylinderGeomv4";
}

//______________________________________________________________________________
const char *PHG4CylinderGeomv4::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeomv4*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderGeomv4::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeomv4*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderGeomv4::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeomv4*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderGeomv4::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeomv4*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4CylinderGeomv4::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderGeomv4.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderGeomv4::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderGeomv4::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4CylinderGeomv4(void *p) {
      return  p ? new(p) ::PHG4CylinderGeomv4 : new ::PHG4CylinderGeomv4;
   }
   static void *newArray_PHG4CylinderGeomv4(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4CylinderGeomv4[nElements] : new ::PHG4CylinderGeomv4[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4CylinderGeomv4(void *p) {
      delete ((::PHG4CylinderGeomv4*)p);
   }
   static void deleteArray_PHG4CylinderGeomv4(void *p) {
      delete [] ((::PHG4CylinderGeomv4*)p);
   }
   static void destruct_PHG4CylinderGeomv4(void *p) {
      typedef ::PHG4CylinderGeomv4 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderGeomv4

namespace {
  void TriggerDictionaryInitialization_PHG4CylinderGeomv4_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4CylinderGeomv4_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4CylinderGeomv4;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4CylinderGeomv4_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERGEOMV4_H
#define G4DETECTORS_PHG4CYLINDERGEOMV4_H

#include "PHG4CylinderGeom.h"

#include <iostream>            // for cout, ostream

class PHG4CylinderGeomv4: public PHG4CylinderGeom
{
 public:
  PHG4CylinderGeomv4();
  PHG4CylinderGeomv4(const int lnsensors,
		    const int lnz,
		    const int nspc,
		    int nsc,
		    const int nstag,
		    const double lr,
		    const double rs,
		    const double szs,
		    const double sps,
		    const double sxo,
		    double syo,
		    const double szsp,
		    const double sys,
		    const double tck,
		    const double st) :
    N_sensors_in_layer(lnsensors),
    layer(-1),
    layer_radius(lr),
    radius_stagger(rs),
    layer_NZ(lnz),
    segment_z_step(szs),
    segment_phi_step(sps),
    sensor_x_offset(sxo),
    sensor_y_offset(syo),
    N_strip_columns(nsc),
    N_strips_per_column(nspc),
    N_staggers(nstag),
    strip_z_spacing(szsp),
    strip_y_spacing(sys),
    thickness(tck),
    strip_tilt(st)
  {}

  virtual ~PHG4CylinderGeomv4() {}

  void identify(std::ostream& os = std::cout) const;
  void set_layer(const int i) {layer = i;}
  int get_layer() const {return layer;}
  double get_radius() const {return layer_radius;}

  void find_segment_center(const int segment_z_bin, const int segment_phi_bin, double location[]);
  void find_strip_center(const int segment_z_bin, const int segment_phi_bin, const int strip_column, const int strip_index, double location[]);

  double get_thickness() const {return thickness;}
  double get_strip_y_spacing() const {return strip_y_spacing;}
  double get_strip_z_spacing() const {return strip_z_spacing;}
  double get_strip_tilt() const {return strip_tilt;}  
  double get_sensor_x_offset() const {return sensor_x_offset;}  
  double get_sensor_y_offset() const {return sensor_y_offset;}  
  int get_N_strip_columns() const {return N_strip_columns;}
  int get_N_strips_per_column() const {return N_strips_per_column;}
  int get_N_sensors_in_layer() const {return N_sensors_in_layer;}

protected:

  int N_sensors_in_layer;
  int layer;

  // finding the center of a sensor given ladder_segment_z and ladder_
  double layer_radius;
  double radius_stagger;
  int layer_NZ;
  double segment_z_step;
  double segment_phi_step;
  double sensor_x_offset;
  double sensor_y_offset;

  // navigation within a sensor
  //double strip_x_offset;
  // double strip_y_offset;
  int N_strip_columns;
  int N_strips_per_column;
  int N_staggers;
  double strip_z_spacing;
  double strip_y_spacing;
  double thickness;
  double strip_tilt;
  
  ClassDef(PHG4CylinderGeomv4,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4CylinderGeomv4", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4CylinderGeomv4_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4CylinderGeomv4_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4CylinderGeomv4_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4CylinderGeomv4_Dict() {
  TriggerDictionaryInitialization_PHG4CylinderGeomv4_Dict_Impl();
}
