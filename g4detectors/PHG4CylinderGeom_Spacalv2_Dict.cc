// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4CylinderGeom_Spacalv2_Dict
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
#include "PHG4CylinderGeom_Spacalv2.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4CylinderGeom_Spacalv2(void *p = 0);
   static void *newArray_PHG4CylinderGeom_Spacalv2(Long_t size, void *p);
   static void delete_PHG4CylinderGeom_Spacalv2(void *p);
   static void deleteArray_PHG4CylinderGeom_Spacalv2(void *p);
   static void destruct_PHG4CylinderGeom_Spacalv2(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderGeom_Spacalv2*)
   {
      ::PHG4CylinderGeom_Spacalv2 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderGeom_Spacalv2 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderGeom_Spacalv2", ::PHG4CylinderGeom_Spacalv2::Class_Version(), "", 27,
                  typeid(::PHG4CylinderGeom_Spacalv2), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderGeom_Spacalv2::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderGeom_Spacalv2) );
      instance.SetNew(&new_PHG4CylinderGeom_Spacalv2);
      instance.SetNewArray(&newArray_PHG4CylinderGeom_Spacalv2);
      instance.SetDelete(&delete_PHG4CylinderGeom_Spacalv2);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderGeom_Spacalv2);
      instance.SetDestructor(&destruct_PHG4CylinderGeom_Spacalv2);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderGeom_Spacalv2*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderGeom_Spacalv2*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv2*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderGeom_Spacalv2::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderGeom_Spacalv2::Class_Name()
{
   return "PHG4CylinderGeom_Spacalv2";
}

//______________________________________________________________________________
const char *PHG4CylinderGeom_Spacalv2::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv2*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderGeom_Spacalv2::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv2*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderGeom_Spacalv2::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv2*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderGeom_Spacalv2::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv2*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4CylinderGeom_Spacalv2::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderGeom_Spacalv2.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderGeom_Spacalv2::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderGeom_Spacalv2::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4CylinderGeom_Spacalv2(void *p) {
      return  p ? new(p) ::PHG4CylinderGeom_Spacalv2 : new ::PHG4CylinderGeom_Spacalv2;
   }
   static void *newArray_PHG4CylinderGeom_Spacalv2(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4CylinderGeom_Spacalv2[nElements] : new ::PHG4CylinderGeom_Spacalv2[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4CylinderGeom_Spacalv2(void *p) {
      delete ((::PHG4CylinderGeom_Spacalv2*)p);
   }
   static void deleteArray_PHG4CylinderGeom_Spacalv2(void *p) {
      delete [] ((::PHG4CylinderGeom_Spacalv2*)p);
   }
   static void destruct_PHG4CylinderGeom_Spacalv2(void *p) {
      typedef ::PHG4CylinderGeom_Spacalv2 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderGeom_Spacalv2

namespace {
  void TriggerDictionaryInitialization_PHG4CylinderGeom_Spacalv2_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4CylinderGeom_Spacalv2_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4CylinderGeom_Spacalv2;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4CylinderGeom_Spacalv2_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $$Id: PHG4CylinderGeom_Spacalv2.h,v 1.3 2014/08/28 22:18:35 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.3 $$
 * \date $$Date: 2014/08/28 22:18:35 $$
 */
#ifndef G4DETECTORS_PHG4CYLINDERGEOMSPACALV2_H
#define G4DETECTORS_PHG4CYLINDERGEOMSPACALV2_H

#include "PHG4CylinderGeom_Spacalv1.h"

#include <cmath>
#include <iostream>                     // for cout, ostream

class PHParameters;

class PHG4CylinderGeom_Spacalv2 : public PHG4CylinderGeom_Spacalv1
{
public:
  PHG4CylinderGeom_Spacalv2();

  virtual
  ~PHG4CylinderGeom_Spacalv2()
  {
  }

  virtual void
  identify(std::ostream& os = std::cout) const;
  virtual void
  Print(Option_t* option = "") const;
  virtual void
  SetDefault();

  //! load parameters from PHParameters, which interface to Database/XML/ROOT files
  virtual void ImportParameters(const PHParameters & param);

  virtual
  int
  get_azimuthal_n_sec() const;

  virtual
  void
  set_azimuthal_n_sec(int azimuthalNSec);

  //! azimuthal tilt in rad
  double
  get_azimuthal_tilt() const
  {
    return azimuthal_tilt;
  }

  //! azimuthal tilt in rad
  void
  set_azimuthal_tilt(double azimuthalTilt)
  {
    azimuthal_tilt = azimuthalTilt;
  }

  virtual
  bool is_azimuthal_seg_visible() const;

  virtual
  void set_azimuthal_seg_visible(bool b = true);

  double get_polar_taper_ratio() const
    {
      return polar_taper_ratio;
    }

  void
  set_polar_taper_ratio(double polarTaperRatio)
  {
    polar_taper_ratio = polarTaperRatio;
  }

  double get_half_polar_taper_angle() const;

  double
  get_sec_azimuthal_width() const;

  double
  get_sec_depth() const
  {
    const double available_depth = get_thickness()
        - (sqrt(
            (get_radius()) * (get_radius())
                + (get_sec_azimuthal_width() / 2)
                    * (get_sec_azimuthal_width() / 2)) - get_radius()) - get_assembly_spacing();
    if (available_depth < get_sec_azimuthal_width())
      return NAN;
    else
      return sqrt(
          available_depth * available_depth
              - get_sec_azimuthal_width() * get_sec_azimuthal_width());
  }

  double
  get_block_width() const
  {
    return get_sec_azimuthal_width() - 0.0001;
  }

  double
  get_block_depth() const
  {
    return sqrt(get_sec_depth()*get_sec_depth() - get_polar_taper_ratio()* get_polar_taper_ratio()* get_block_width()*get_block_width()) - get_assembly_spacing();
  }


  double get_assembly_spacing() const
    {
      return assembly_spacing;
    }

  void
  set_assembly_spacing(double assemblySpacing)
  {
    assembly_spacing = assemblySpacing;
  }

  //! regulated fiber distance in the tapering direction
  double
  get_reg_fiber_grid_distance_taper() const;

  //! regulated fiber distance in the non-tapering direction
  double
  get_reg_fiber_grid_distance_nontaper() const;


protected:

  int azimuthal_n_sec;

  //! azimuthal tilt in rad
  double azimuthal_tilt;
  bool azimuthal_seg_visible;
  double polar_taper_ratio;
  double assembly_spacing;

  ClassDef(PHG4CylinderGeom_Spacalv2,2)


};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4CylinderGeom_Spacalv2", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4CylinderGeom_Spacalv2_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4CylinderGeom_Spacalv2_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4CylinderGeom_Spacalv2_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4CylinderGeom_Spacalv2_Dict() {
  TriggerDictionaryInitialization_PHG4CylinderGeom_Spacalv2_Dict_Impl();
}
