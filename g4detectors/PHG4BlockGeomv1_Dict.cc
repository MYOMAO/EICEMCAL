// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4BlockGeomv1_Dict
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
#include "PHG4BlockGeomv1.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4BlockGeomv1(void *p = 0);
   static void *newArray_PHG4BlockGeomv1(Long_t size, void *p);
   static void delete_PHG4BlockGeomv1(void *p);
   static void deleteArray_PHG4BlockGeomv1(void *p);
   static void destruct_PHG4BlockGeomv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4BlockGeomv1*)
   {
      ::PHG4BlockGeomv1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4BlockGeomv1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4BlockGeomv1", ::PHG4BlockGeomv1::Class_Version(), "", 15,
                  typeid(::PHG4BlockGeomv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4BlockGeomv1::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4BlockGeomv1) );
      instance.SetNew(&new_PHG4BlockGeomv1);
      instance.SetNewArray(&newArray_PHG4BlockGeomv1);
      instance.SetDelete(&delete_PHG4BlockGeomv1);
      instance.SetDeleteArray(&deleteArray_PHG4BlockGeomv1);
      instance.SetDestructor(&destruct_PHG4BlockGeomv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4BlockGeomv1*)
   {
      return GenerateInitInstanceLocal((::PHG4BlockGeomv1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4BlockGeomv1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4BlockGeomv1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4BlockGeomv1::Class_Name()
{
   return "PHG4BlockGeomv1";
}

//______________________________________________________________________________
const char *PHG4BlockGeomv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockGeomv1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4BlockGeomv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockGeomv1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4BlockGeomv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockGeomv1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4BlockGeomv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockGeomv1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4BlockGeomv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4BlockGeomv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4BlockGeomv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4BlockGeomv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4BlockGeomv1(void *p) {
      return  p ? new(p) ::PHG4BlockGeomv1 : new ::PHG4BlockGeomv1;
   }
   static void *newArray_PHG4BlockGeomv1(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4BlockGeomv1[nElements] : new ::PHG4BlockGeomv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4BlockGeomv1(void *p) {
      delete ((::PHG4BlockGeomv1*)p);
   }
   static void deleteArray_PHG4BlockGeomv1(void *p) {
      delete [] ((::PHG4BlockGeomv1*)p);
   }
   static void destruct_PHG4BlockGeomv1(void *p) {
      typedef ::PHG4BlockGeomv1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4BlockGeomv1

namespace {
  void TriggerDictionaryInitialization_PHG4BlockGeomv1_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4BlockGeomv1_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4BlockGeomv1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4BlockGeomv1_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BLOCKGEOMV1_H
#define G4DETECTORS_PHG4BLOCKGEOMV1_H

#include "PHG4BlockGeom.h"

#include <iostream>         // for cout, ostream

class PHG4BlockGeomv1: public PHG4BlockGeom
{
 public:
  PHG4BlockGeomv1();
  PHG4BlockGeomv1( const int layer,
                   const double sizex, const double sizey, const double sizez, 
                   const double centerx, const double centery, const double centerz,
                   const double zrot );

  virtual ~PHG4BlockGeomv1() {}

  void identify(std::ostream& os = std::cout) const;

  int get_layer() const {return _layer;}
  double get_width() const {return _size[0];}
  double get_thickness() const {return _size[1];}
  double get_length() const {return _size[2];}
  double get_center_x() const {return _center[0];}
  double get_center_y() const {return _center[1];}
  double get_center_z() const {return _center[2];}
  double get_z_rot() const {return _rotation_z;}

  double get_size_x() const {return _size[0];}
  double get_size_y() const {return _size[1];}
  double get_size_z() const {return _size[2];}

  double get_rot_matrix(const int i, const int j) const {return _rot_matrix[i][j];}

  void set_layer(const int i) {_layer = i;}

  // size in local coordinates
  void set_size(const double sizex, const double sizey, const double sizez);

  void set_center(const double centerx, const double centery, const double centerz);
  void set_z_rot(const double zrot) {_build_rot_matrix(); _rotation_z = zrot;}

  void convert_local_to_global(double, double, double,
                               double &, double &, double &) const;
  void convert_global_x_to_local(double, double, double, 
                                 double &, double &, double &) const;

 protected:
  int _layer;
  double _size[3];
  double _center[3];
  double _rotation_z;

  void _build_rot_matrix();
  double _rot_matrix[3][3];  // global -> local coordinates rotation matrix

  ClassDef(PHG4BlockGeomv1,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4BlockGeomv1", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4BlockGeomv1_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4BlockGeomv1_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4BlockGeomv1_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4BlockGeomv1_Dict() {
  TriggerDictionaryInitialization_PHG4BlockGeomv1_Dict_Impl();
}
