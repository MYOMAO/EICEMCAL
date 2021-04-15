// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4CylinderGeom_Spacalv1_Dict
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
#include "PHG4CylinderGeom_Spacalv1.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4CylinderGeom_Spacalv1(void *p = 0);
   static void *newArray_PHG4CylinderGeom_Spacalv1(Long_t size, void *p);
   static void delete_PHG4CylinderGeom_Spacalv1(void *p);
   static void deleteArray_PHG4CylinderGeom_Spacalv1(void *p);
   static void destruct_PHG4CylinderGeom_Spacalv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderGeom_Spacalv1*)
   {
      ::PHG4CylinderGeom_Spacalv1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderGeom_Spacalv1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderGeom_Spacalv1", ::PHG4CylinderGeom_Spacalv1::Class_Version(), "", 28,
                  typeid(::PHG4CylinderGeom_Spacalv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderGeom_Spacalv1::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderGeom_Spacalv1) );
      instance.SetNew(&new_PHG4CylinderGeom_Spacalv1);
      instance.SetNewArray(&newArray_PHG4CylinderGeom_Spacalv1);
      instance.SetDelete(&delete_PHG4CylinderGeom_Spacalv1);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderGeom_Spacalv1);
      instance.SetDestructor(&destruct_PHG4CylinderGeom_Spacalv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderGeom_Spacalv1*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderGeom_Spacalv1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderGeom_Spacalv1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderGeom_Spacalv1::Class_Name()
{
   return "PHG4CylinderGeom_Spacalv1";
}

//______________________________________________________________________________
const char *PHG4CylinderGeom_Spacalv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderGeom_Spacalv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderGeom_Spacalv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderGeom_Spacalv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4CylinderGeom_Spacalv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderGeom_Spacalv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderGeom_Spacalv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderGeom_Spacalv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4CylinderGeom_Spacalv1(void *p) {
      return  p ? new(p) ::PHG4CylinderGeom_Spacalv1 : new ::PHG4CylinderGeom_Spacalv1;
   }
   static void *newArray_PHG4CylinderGeom_Spacalv1(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4CylinderGeom_Spacalv1[nElements] : new ::PHG4CylinderGeom_Spacalv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4CylinderGeom_Spacalv1(void *p) {
      delete ((::PHG4CylinderGeom_Spacalv1*)p);
   }
   static void deleteArray_PHG4CylinderGeom_Spacalv1(void *p) {
      delete [] ((::PHG4CylinderGeom_Spacalv1*)p);
   }
   static void destruct_PHG4CylinderGeom_Spacalv1(void *p) {
      typedef ::PHG4CylinderGeom_Spacalv1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderGeom_Spacalv1

namespace ROOT {
   static TClass *maplEintcOdoublegR_Dictionary();
   static void maplEintcOdoublegR_TClassManip(TClass*);
   static void *new_maplEintcOdoublegR(void *p = 0);
   static void *newArray_maplEintcOdoublegR(Long_t size, void *p);
   static void delete_maplEintcOdoublegR(void *p);
   static void deleteArray_maplEintcOdoublegR(void *p);
   static void destruct_maplEintcOdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,double>*)
   {
      map<int,double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,double>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,double>", -2, "map", 100,
                  typeid(map<int,double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,double>) );
      instance.SetNew(&new_maplEintcOdoublegR);
      instance.SetNewArray(&newArray_maplEintcOdoublegR);
      instance.SetDelete(&delete_maplEintcOdoublegR);
      instance.SetDeleteArray(&deleteArray_maplEintcOdoublegR);
      instance.SetDestructor(&destruct_maplEintcOdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,double> >()));

      ::ROOT::AddClassAlternate("map<int,double>","std::map<int, double, std::less<int>, std::allocator<std::pair<int const, double> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,double>*)0x0)->GetClass();
      maplEintcOdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,double> : new map<int,double>;
   }
   static void *newArray_maplEintcOdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,double>[nElements] : new map<int,double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOdoublegR(void *p) {
      delete ((map<int,double>*)p);
   }
   static void deleteArray_maplEintcOdoublegR(void *p) {
      delete [] ((map<int,double>*)p);
   }
   static void destruct_maplEintcOdoublegR(void *p) {
      typedef map<int,double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,double>

namespace {
  void TriggerDictionaryInitialization_PHG4CylinderGeom_Spacalv1_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4CylinderGeom_Spacalv1_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4CylinderGeom_Spacalv1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4CylinderGeom_Spacalv1_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $$Id: PHG4CylinderGeom_Spacalv1.h,v 1.2 2014/08/12 03:49:12 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.2 $$
 * \date $$Date: 2014/08/12 03:49:12 $$
 */
#ifndef G4DETECTORS_PHG4CYLINDERGEOMSPACALV1_H
#define G4DETECTORS_PHG4CYLINDERGEOMSPACALV1_H

#include "PHG4CylinderGeomv2.h"

#include <iostream>              // for cout, ostream
#include <map>
#include <string>

class PHParameters;

class PHG4CylinderGeom_Spacalv1 : public PHG4CylinderGeomv2
{
 public:
  /** @name Ctor DTor and IDs
   */
  ///@{
  PHG4CylinderGeom_Spacalv1();

  virtual ~PHG4CylinderGeom_Spacalv1()
  {
    sector_map.clear();
  }

  virtual void
  identify(std::ostream &os = std::cout) const;

  virtual void
  Print(Option_t *option = "") const;

  virtual void
  SetDefault();

  //! load parameters from PHParameters, which interface to Database/XML/ROOT files
  virtual void ImportParameters(const PHParameters &param);

  ///@}

  /** @name Set Cylinder Geometry
   */
  ///@{
  double
  get_max_radius() const
  {
    return get_radius() + get_thickness();
  }

  double
  get_half_radius() const
  {
    return get_radius() + get_thickness() / 2.;
  }

  double
  get_length() const
  {
    return get_zmax() - get_zmin();
  }

  double
  get_xpos() const
  {
    return xpos;
  }

  void
  set_xpos(double xpos)
  {
    this->xpos = xpos;
  }

  double
  get_ypos() const
  {
    return ypos;
  }

  void
  set_ypos(double ypos)
  {
    this->ypos = ypos;
  }

  double
  get_zpos() const
  {
    return zpos;
  }

  void
  set_zpos(double zpos)
  {
    this->zpos = zpos;
  }

  ///@}

  /** @name Azimuthal segments
   */
  ///@{
  virtual int
  get_azimuthal_n_sec() const;

  virtual double
  get_azimuthal_distance() const;

  virtual double
  get_z_distance() const;

  //! sector map sector_ID -> azimuthal rotation.
  typedef std::map<int, double> sector_map_t;

  //! sector map sector_ID -> azimuthal rotation.
  const sector_map_t &
  get_sector_map() const
  {
    return sector_map;
  }

  //! sector map sector_ID -> azimuthal rotation.
  sector_map_t &
  get_sector_map()
  {
    return sector_map;
  }

  //! load a default map that populate all the sectors
  void
  init_default_sector_map();

  ///@}

  /** @name Fiber geometry
   */
  ///@{
  double
  get_fiber_outer_r() const
  {
    return get_fiber_clading_thickness() + get_fiber_core_diameter() / 2;
  }

  double
  get_fiber_clading_thickness() const
  {
    return fiber_clading_thickness;
  }

  void
  set_fiber_clading_thickness(double fiberCladingThickness)
  {
    fiber_clading_thickness = fiberCladingThickness;
  }

  double
  get_fiber_core_diameter() const
  {
    return fiber_core_diameter;
  }

  void
  set_fiber_core_diameter(double fiberCoreDiameter)
  {
    fiber_core_diameter = fiberCoreDiameter;
  }

  double
  get_fiber_distance() const
  {
    return fiber_distance;
  }

  void
  set_fiber_distance(double fiberDistance)
  {
    fiber_distance = fiberDistance;
  }

  ///@}

  /** @name Materials
   */
  ///@{
  std::string
  get_absorber_mat() const
  {
    return absorber_mat;
  }

  void
  set_absorber_mat(const std::string &absorberMat)
  {
    absorber_mat = absorberMat;
  }

  std::string
  get_fiber_clading_mat() const
  {
    return fiber_clading_mat;
  }

  void
  set_fiber_clading_mat(const std::string &fiberCladingMat)
  {
    fiber_clading_mat = fiberCladingMat;
  }

  std::string
  get_fiber_core_mat() const
  {
    return fiber_core_mat;
  }

  void
  set_fiber_core_mat(const std::string &fiberCoreMat)
  {
    fiber_core_mat = fiberCoreMat;
  }

  ///@}

  /** @name General options
   */
  ///@{

  double
  get_fiber_core_step_size() const
  {
    return get_fiber_core_diameter() / 10.;
  }

  enum config_t
  {

    //! fiber always placed radially
    kNonProjective = 0,
    //! alias of above, more explicit
    k1DProjectiveSpacal = kNonProjective,

    //! Fully projective spacal with 2D tapered modules
    kFullProjective_2DTaper = 2,

    //! Fully projective spacal with 2D tapered modules. To speed up construction, same-length fiber is used cross one tower
    kFullProjective_2DTaper_SameLengthFiberPerTower = 3,

    //! Fully projective spacal with 2D tapered modules and allow azimuthal tilts
    kFullProjective_2DTaper_Tilted = 4,

    //! Fully projective spacal with 2D tapered modules and allow azimuthal tilts. To speed up construction, same-length fiber is used cross one tower
    kFullProjective_2DTaper_Tilted_SameLengthFiberPerTower = 5,

    //! alias of above, the default 2D-projective SPACAL
    k2DProjectiveSpacal = kFullProjective_2DTaper_Tilted_SameLengthFiberPerTower,

    //! max allowed value, for boundary cross check
    kInvalidSpacalConfig
  };

  config_t
  get_config() const
  {
    return config;
  }

  void
  set_config(config_t config)
  {
    this->config = config;
  }

  bool
  is_virualize_fiber() const
  {
    return virualize_fiber;
  }

  virtual bool is_azimuthal_seg_visible() const
  {
    return false;
  }

  void
  set_virualize_fiber(bool virualizeFiber = true)
  {
    virualize_fiber = virualizeFiber;
  }

  int get_construction_verbose() const
  {
    return construction_verbose;
  }

  void
  set_construction_verbose(int constructionVerbose)
  {
    construction_verbose = constructionVerbose;
  }

  ///@}

 protected:
  std::string absorber_mat;
  std::string fiber_core_mat;
  std::string fiber_clading_mat;
  double xpos;
  double ypos;
  double zpos;
  double fiber_core_diameter;
  double fiber_clading_thickness;
  double fiber_distance;
  config_t config;
  bool virualize_fiber;
  int construction_verbose;

  //! sector map sector_ID -> azimuthal rotation.
  sector_map_t sector_map;

  ClassDef(PHG4CylinderGeom_Spacalv1, 2)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4CylinderGeom_Spacalv1", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4CylinderGeom_Spacalv1_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4CylinderGeom_Spacalv1_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4CylinderGeom_Spacalv1_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4CylinderGeom_Spacalv1_Dict() {
  TriggerDictionaryInitialization_PHG4CylinderGeom_Spacalv1_Dict_Impl();
}
