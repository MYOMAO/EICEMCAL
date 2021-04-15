// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4CylinderCellGeom_Spacalv1_Dict
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
#include "PHG4CylinderCellGeom_Spacalv1.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4CylinderCellGeom_Spacalv1(void *p = 0);
   static void *newArray_PHG4CylinderCellGeom_Spacalv1(Long_t size, void *p);
   static void delete_PHG4CylinderCellGeom_Spacalv1(void *p);
   static void deleteArray_PHG4CylinderCellGeom_Spacalv1(void *p);
   static void destruct_PHG4CylinderCellGeom_Spacalv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderCellGeom_Spacalv1*)
   {
      ::PHG4CylinderCellGeom_Spacalv1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderCellGeom_Spacalv1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderCellGeom_Spacalv1", ::PHG4CylinderCellGeom_Spacalv1::Class_Version(), "", 30,
                  typeid(::PHG4CylinderCellGeom_Spacalv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderCellGeom_Spacalv1::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderCellGeom_Spacalv1) );
      instance.SetNew(&new_PHG4CylinderCellGeom_Spacalv1);
      instance.SetNewArray(&newArray_PHG4CylinderCellGeom_Spacalv1);
      instance.SetDelete(&delete_PHG4CylinderCellGeom_Spacalv1);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderCellGeom_Spacalv1);
      instance.SetDestructor(&destruct_PHG4CylinderCellGeom_Spacalv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderCellGeom_Spacalv1*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderCellGeom_Spacalv1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderCellGeom_Spacalv1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderCellGeom_Spacalv1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderCellGeom_Spacalv1::Class_Name()
{
   return "PHG4CylinderCellGeom_Spacalv1";
}

//______________________________________________________________________________
const char *PHG4CylinderCellGeom_Spacalv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellGeom_Spacalv1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderCellGeom_Spacalv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellGeom_Spacalv1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderCellGeom_Spacalv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellGeom_Spacalv1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderCellGeom_Spacalv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellGeom_Spacalv1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4CylinderCellGeom_Spacalv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderCellGeom_Spacalv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderCellGeom_Spacalv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderCellGeom_Spacalv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4CylinderCellGeom_Spacalv1(void *p) {
      return  p ? new(p) ::PHG4CylinderCellGeom_Spacalv1 : new ::PHG4CylinderCellGeom_Spacalv1;
   }
   static void *newArray_PHG4CylinderCellGeom_Spacalv1(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4CylinderCellGeom_Spacalv1[nElements] : new ::PHG4CylinderCellGeom_Spacalv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4CylinderCellGeom_Spacalv1(void *p) {
      delete ((::PHG4CylinderCellGeom_Spacalv1*)p);
   }
   static void deleteArray_PHG4CylinderCellGeom_Spacalv1(void *p) {
      delete [] ((::PHG4CylinderCellGeom_Spacalv1*)p);
   }
   static void destruct_PHG4CylinderCellGeom_Spacalv1(void *p) {
      typedef ::PHG4CylinderCellGeom_Spacalv1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderCellGeom_Spacalv1

namespace ROOT {
   static TClass *maplEintcOpairlEdoublecOdoublegRsPgR_Dictionary();
   static void maplEintcOpairlEdoublecOdoublegRsPgR_TClassManip(TClass*);
   static void *new_maplEintcOpairlEdoublecOdoublegRsPgR(void *p = 0);
   static void *newArray_maplEintcOpairlEdoublecOdoublegRsPgR(Long_t size, void *p);
   static void delete_maplEintcOpairlEdoublecOdoublegRsPgR(void *p);
   static void deleteArray_maplEintcOpairlEdoublecOdoublegRsPgR(void *p);
   static void destruct_maplEintcOpairlEdoublecOdoublegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,pair<double,double> >*)
   {
      map<int,pair<double,double> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,pair<double,double> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,pair<double,double> >", -2, "map", 100,
                  typeid(map<int,pair<double,double> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOpairlEdoublecOdoublegRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,pair<double,double> >) );
      instance.SetNew(&new_maplEintcOpairlEdoublecOdoublegRsPgR);
      instance.SetNewArray(&newArray_maplEintcOpairlEdoublecOdoublegRsPgR);
      instance.SetDelete(&delete_maplEintcOpairlEdoublecOdoublegRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOpairlEdoublecOdoublegRsPgR);
      instance.SetDestructor(&destruct_maplEintcOpairlEdoublecOdoublegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,pair<double,double> > >()));

      ::ROOT::AddClassAlternate("map<int,pair<double,double> >","std::map<int, std::pair<double, double>, std::less<int>, std::allocator<std::pair<int const, std::pair<double, double> > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,pair<double,double> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOpairlEdoublecOdoublegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,pair<double,double> >*)0x0)->GetClass();
      maplEintcOpairlEdoublecOdoublegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOpairlEdoublecOdoublegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOpairlEdoublecOdoublegRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,pair<double,double> > : new map<int,pair<double,double> >;
   }
   static void *newArray_maplEintcOpairlEdoublecOdoublegRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,pair<double,double> >[nElements] : new map<int,pair<double,double> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOpairlEdoublecOdoublegRsPgR(void *p) {
      delete ((map<int,pair<double,double> >*)p);
   }
   static void deleteArray_maplEintcOpairlEdoublecOdoublegRsPgR(void *p) {
      delete [] ((map<int,pair<double,double> >*)p);
   }
   static void destruct_maplEintcOpairlEdoublecOdoublegRsPgR(void *p) {
      typedef map<int,pair<double,double> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,pair<double,double> >

namespace ROOT {
   static TClass *maplEintcOintgR_Dictionary();
   static void maplEintcOintgR_TClassManip(TClass*);
   static void *new_maplEintcOintgR(void *p = 0);
   static void *newArray_maplEintcOintgR(Long_t size, void *p);
   static void delete_maplEintcOintgR(void *p);
   static void deleteArray_maplEintcOintgR(void *p);
   static void destruct_maplEintcOintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,int>*)
   {
      map<int,int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,int>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,int>", -2, "map", 100,
                  typeid(map<int,int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOintgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,int>) );
      instance.SetNew(&new_maplEintcOintgR);
      instance.SetNewArray(&newArray_maplEintcOintgR);
      instance.SetDelete(&delete_maplEintcOintgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOintgR);
      instance.SetDestructor(&destruct_maplEintcOintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,int> >()));

      ::ROOT::AddClassAlternate("map<int,int>","std::map<int, int, std::less<int>, std::allocator<std::pair<int const, int> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,int>*)0x0)->GetClass();
      maplEintcOintgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,int> : new map<int,int>;
   }
   static void *newArray_maplEintcOintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,int>[nElements] : new map<int,int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOintgR(void *p) {
      delete ((map<int,int>*)p);
   }
   static void deleteArray_maplEintcOintgR(void *p) {
      delete [] ((map<int,int>*)p);
   }
   static void destruct_maplEintcOintgR(void *p) {
      typedef map<int,int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,int>

namespace {
  void TriggerDictionaryInitialization_PHG4CylinderCellGeom_Spacalv1_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4CylinderCellGeom_Spacalv1_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4CylinderCellGeom_Spacalv1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4CylinderCellGeom_Spacalv1_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $Id: $                                                                                             

/*!
 * \file PHG4CylinderCellGeomSpacalv1.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4DETECTORS_PHG4CYLINDERCELLGEOMSPACALV1_H
#define G4DETECTORS_PHG4CYLINDERCELLGEOMSPACALV1_H

#include "PHG4CylinderCellGeom.h"

#include <iostream>                // for cout, ostream
#include <map>
#include <utility>                 // for pair

/*!
 * \brief PHG4CylinderCellGeom_Spacalv1
 */
class PHG4CylinderCellGeom_Spacalv1 : public PHG4CylinderCellGeom
{
public:
  PHG4CylinderCellGeom_Spacalv1();
  virtual
  ~PHG4CylinderCellGeom_Spacalv1();
  void
  identify(std::ostream& os = std::cout) const;

  virtual std::pair<double, double>
  get_zbounds(const int ibin) const;
  virtual std::pair<double, double>
  get_etabounds(const int ibin) const;

  virtual double
  get_etacenter(const int ibin) const;
  virtual double
  get_zcenter(const int ibin) const;

  virtual int
  get_etabin(const double eta) const;
  virtual int
  get_zbin(const double z) const;

  void
  set_zbounds(const int ibin, const std::pair<double, double> & bounds);
  void
  set_etabounds(const int ibin, const std::pair<double, double> & bounds);

  typedef std::pair<double, double> bound_t;
  typedef std::map<int, bound_t> bound_map_t;

  const bound_map_t &
  get_eta_bound_map() const
  {
    map_consistency_check();
    return eta_bound_map;
  }

  void
  set_eta_bound_map(const bound_map_t & etaBoundMap)
  {
    eta_bound_map = etaBoundMap;
  }

  const bound_map_t &
  get_z_bound_map() const
  {
    map_consistency_check();
    return z_bound_map;
  }

  void
  set_z_bound_map(const bound_map_t & boundMap)
  {
    z_bound_map = boundMap;
  }

  //! map tower_z_ID -> eta_bin number
  typedef std::map<int, int> tower_z_ID_eta_bin_map_t;

  //! map tower_z_ID -> eta_bin number
  const tower_z_ID_eta_bin_map_t &
  get_tower_z_ID_eta_bin_map() const
  {
    return tower_z_ID_eta_bin_map;
  }

  virtual int
  get_etabin_block(const int tower_z_ID) const;

  //! map tower_z_ID -> eta_bin number for blocks
  void
  set_tower_z_ID_eta_bin_map(const tower_z_ID_eta_bin_map_t & m)
  {
    tower_z_ID_eta_bin_map = m;
  }

protected:

  void
  map_consistency_check() const;

  bound_map_t z_bound_map;
  bound_map_t eta_bound_map;

  //! map tower_z_ID -> eta_bin number for blocks
  tower_z_ID_eta_bin_map_t tower_z_ID_eta_bin_map;

ClassDef(PHG4CylinderCellGeom_Spacalv1,2)

};

#endif /* PHG4CYLINDERCELLGEOMSPACALV1_H_ */

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4CylinderCellGeom_Spacalv1", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4CylinderCellGeom_Spacalv1_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4CylinderCellGeom_Spacalv1_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4CylinderCellGeom_Spacalv1_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4CylinderCellGeom_Spacalv1_Dict() {
  TriggerDictionaryInitialization_PHG4CylinderCellGeom_Spacalv1_Dict_Impl();
}
