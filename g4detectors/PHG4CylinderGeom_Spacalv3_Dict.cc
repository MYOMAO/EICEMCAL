// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4CylinderGeom_Spacalv3_Dict
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
#include "PHG4CylinderGeom_Spacalv3.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4CylinderGeom_Spacalv3(void *p = 0);
   static void *newArray_PHG4CylinderGeom_Spacalv3(Long_t size, void *p);
   static void delete_PHG4CylinderGeom_Spacalv3(void *p);
   static void deleteArray_PHG4CylinderGeom_Spacalv3(void *p);
   static void destruct_PHG4CylinderGeom_Spacalv3(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderGeom_Spacalv3*)
   {
      ::PHG4CylinderGeom_Spacalv3 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderGeom_Spacalv3 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderGeom_Spacalv3", ::PHG4CylinderGeom_Spacalv3::Class_Version(), "", 29,
                  typeid(::PHG4CylinderGeom_Spacalv3), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderGeom_Spacalv3::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderGeom_Spacalv3) );
      instance.SetNew(&new_PHG4CylinderGeom_Spacalv3);
      instance.SetNewArray(&newArray_PHG4CylinderGeom_Spacalv3);
      instance.SetDelete(&delete_PHG4CylinderGeom_Spacalv3);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderGeom_Spacalv3);
      instance.SetDestructor(&destruct_PHG4CylinderGeom_Spacalv3);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderGeom_Spacalv3*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderGeom_Spacalv3*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv3*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PHG4CylinderGeom_Spacalv3cLcLgeom_tower(void *p = 0);
   static void *newArray_PHG4CylinderGeom_Spacalv3cLcLgeom_tower(Long_t size, void *p);
   static void delete_PHG4CylinderGeom_Spacalv3cLcLgeom_tower(void *p);
   static void deleteArray_PHG4CylinderGeom_Spacalv3cLcLgeom_tower(void *p);
   static void destruct_PHG4CylinderGeom_Spacalv3cLcLgeom_tower(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderGeom_Spacalv3::geom_tower*)
   {
      ::PHG4CylinderGeom_Spacalv3::geom_tower *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderGeom_Spacalv3::geom_tower >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderGeom_Spacalv3::geom_tower", ::PHG4CylinderGeom_Spacalv3::geom_tower::Class_Version(), "", 114,
                  typeid(::PHG4CylinderGeom_Spacalv3::geom_tower), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderGeom_Spacalv3::geom_tower::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderGeom_Spacalv3::geom_tower) );
      instance.SetNew(&new_PHG4CylinderGeom_Spacalv3cLcLgeom_tower);
      instance.SetNewArray(&newArray_PHG4CylinderGeom_Spacalv3cLcLgeom_tower);
      instance.SetDelete(&delete_PHG4CylinderGeom_Spacalv3cLcLgeom_tower);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderGeom_Spacalv3cLcLgeom_tower);
      instance.SetDestructor(&destruct_PHG4CylinderGeom_Spacalv3cLcLgeom_tower);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderGeom_Spacalv3::geom_tower*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderGeom_Spacalv3::geom_tower*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv3::geom_tower*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_PHG4CylinderGeom_Spacalv3cLcLscint_id_coder(void *p);
   static void deleteArray_PHG4CylinderGeom_Spacalv3cLcLscint_id_coder(void *p);
   static void destruct_PHG4CylinderGeom_Spacalv3cLcLscint_id_coder(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderGeom_Spacalv3::scint_id_coder*)
   {
      ::PHG4CylinderGeom_Spacalv3::scint_id_coder *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderGeom_Spacalv3::scint_id_coder >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderGeom_Spacalv3::scint_id_coder", ::PHG4CylinderGeom_Spacalv3::scint_id_coder::Class_Version(), "", 227,
                  typeid(::PHG4CylinderGeom_Spacalv3::scint_id_coder), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderGeom_Spacalv3::scint_id_coder::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderGeom_Spacalv3::scint_id_coder) );
      instance.SetDelete(&delete_PHG4CylinderGeom_Spacalv3cLcLscint_id_coder);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderGeom_Spacalv3cLcLscint_id_coder);
      instance.SetDestructor(&destruct_PHG4CylinderGeom_Spacalv3cLcLscint_id_coder);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderGeom_Spacalv3::scint_id_coder*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderGeom_Spacalv3::scint_id_coder*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv3::scint_id_coder*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderGeom_Spacalv3::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderGeom_Spacalv3::Class_Name()
{
   return "PHG4CylinderGeom_Spacalv3";
}

//______________________________________________________________________________
const char *PHG4CylinderGeom_Spacalv3::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv3*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderGeom_Spacalv3::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv3*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderGeom_Spacalv3::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv3*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderGeom_Spacalv3::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv3*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderGeom_Spacalv3::geom_tower::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderGeom_Spacalv3::geom_tower::Class_Name()
{
   return "PHG4CylinderGeom_Spacalv3::geom_tower";
}

//______________________________________________________________________________
const char *PHG4CylinderGeom_Spacalv3::geom_tower::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv3::geom_tower*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderGeom_Spacalv3::geom_tower::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv3::geom_tower*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderGeom_Spacalv3::geom_tower::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv3::geom_tower*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderGeom_Spacalv3::geom_tower::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv3::geom_tower*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderGeom_Spacalv3::scint_id_coder::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderGeom_Spacalv3::scint_id_coder::Class_Name()
{
   return "PHG4CylinderGeom_Spacalv3::scint_id_coder";
}

//______________________________________________________________________________
const char *PHG4CylinderGeom_Spacalv3::scint_id_coder::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv3::scint_id_coder*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderGeom_Spacalv3::scint_id_coder::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv3::scint_id_coder*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderGeom_Spacalv3::scint_id_coder::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv3::scint_id_coder*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderGeom_Spacalv3::scint_id_coder::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeom_Spacalv3::scint_id_coder*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4CylinderGeom_Spacalv3::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderGeom_Spacalv3.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderGeom_Spacalv3::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderGeom_Spacalv3::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4CylinderGeom_Spacalv3(void *p) {
      return  p ? new(p) ::PHG4CylinderGeom_Spacalv3 : new ::PHG4CylinderGeom_Spacalv3;
   }
   static void *newArray_PHG4CylinderGeom_Spacalv3(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4CylinderGeom_Spacalv3[nElements] : new ::PHG4CylinderGeom_Spacalv3[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4CylinderGeom_Spacalv3(void *p) {
      delete ((::PHG4CylinderGeom_Spacalv3*)p);
   }
   static void deleteArray_PHG4CylinderGeom_Spacalv3(void *p) {
      delete [] ((::PHG4CylinderGeom_Spacalv3*)p);
   }
   static void destruct_PHG4CylinderGeom_Spacalv3(void *p) {
      typedef ::PHG4CylinderGeom_Spacalv3 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderGeom_Spacalv3

//______________________________________________________________________________
void PHG4CylinderGeom_Spacalv3::geom_tower::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderGeom_Spacalv3::geom_tower.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderGeom_Spacalv3::geom_tower::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderGeom_Spacalv3::geom_tower::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4CylinderGeom_Spacalv3cLcLgeom_tower(void *p) {
      return  p ? new(p) ::PHG4CylinderGeom_Spacalv3::geom_tower : new ::PHG4CylinderGeom_Spacalv3::geom_tower;
   }
   static void *newArray_PHG4CylinderGeom_Spacalv3cLcLgeom_tower(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4CylinderGeom_Spacalv3::geom_tower[nElements] : new ::PHG4CylinderGeom_Spacalv3::geom_tower[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4CylinderGeom_Spacalv3cLcLgeom_tower(void *p) {
      delete ((::PHG4CylinderGeom_Spacalv3::geom_tower*)p);
   }
   static void deleteArray_PHG4CylinderGeom_Spacalv3cLcLgeom_tower(void *p) {
      delete [] ((::PHG4CylinderGeom_Spacalv3::geom_tower*)p);
   }
   static void destruct_PHG4CylinderGeom_Spacalv3cLcLgeom_tower(void *p) {
      typedef ::PHG4CylinderGeom_Spacalv3::geom_tower current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderGeom_Spacalv3::geom_tower

//______________________________________________________________________________
void PHG4CylinderGeom_Spacalv3::scint_id_coder::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderGeom_Spacalv3::scint_id_coder.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderGeom_Spacalv3::scint_id_coder::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderGeom_Spacalv3::scint_id_coder::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PHG4CylinderGeom_Spacalv3cLcLscint_id_coder(void *p) {
      delete ((::PHG4CylinderGeom_Spacalv3::scint_id_coder*)p);
   }
   static void deleteArray_PHG4CylinderGeom_Spacalv3cLcLscint_id_coder(void *p) {
      delete [] ((::PHG4CylinderGeom_Spacalv3::scint_id_coder*)p);
   }
   static void destruct_PHG4CylinderGeom_Spacalv3cLcLscint_id_coder(void *p) {
      typedef ::PHG4CylinderGeom_Spacalv3::scint_id_coder current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderGeom_Spacalv3::scint_id_coder

namespace ROOT {
   static TClass *maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR_Dictionary();
   static void maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR_TClassManip(TClass*);
   static void *new_maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR(void *p = 0);
   static void *newArray_maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR(Long_t size, void *p);
   static void delete_maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR(void *p);
   static void deleteArray_maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR(void *p);
   static void destruct_maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,PHG4CylinderGeom_Spacalv3::geom_tower>*)
   {
      map<int,PHG4CylinderGeom_Spacalv3::geom_tower> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,PHG4CylinderGeom_Spacalv3::geom_tower>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,PHG4CylinderGeom_Spacalv3::geom_tower>", -2, "map", 100,
                  typeid(map<int,PHG4CylinderGeom_Spacalv3::geom_tower>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,PHG4CylinderGeom_Spacalv3::geom_tower>) );
      instance.SetNew(&new_maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR);
      instance.SetNewArray(&newArray_maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR);
      instance.SetDelete(&delete_maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR);
      instance.SetDeleteArray(&deleteArray_maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR);
      instance.SetDestructor(&destruct_maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,PHG4CylinderGeom_Spacalv3::geom_tower> >()));

      ::ROOT::AddClassAlternate("map<int,PHG4CylinderGeom_Spacalv3::geom_tower>","std::map<int, PHG4CylinderGeom_Spacalv3::geom_tower, std::less<int>, std::allocator<std::pair<int const, PHG4CylinderGeom_Spacalv3::geom_tower> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,PHG4CylinderGeom_Spacalv3::geom_tower>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,PHG4CylinderGeom_Spacalv3::geom_tower>*)0x0)->GetClass();
      maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,PHG4CylinderGeom_Spacalv3::geom_tower> : new map<int,PHG4CylinderGeom_Spacalv3::geom_tower>;
   }
   static void *newArray_maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,PHG4CylinderGeom_Spacalv3::geom_tower>[nElements] : new map<int,PHG4CylinderGeom_Spacalv3::geom_tower>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR(void *p) {
      delete ((map<int,PHG4CylinderGeom_Spacalv3::geom_tower>*)p);
   }
   static void deleteArray_maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR(void *p) {
      delete [] ((map<int,PHG4CylinderGeom_Spacalv3::geom_tower>*)p);
   }
   static void destruct_maplEintcOPHG4CylinderGeom_Spacalv3cLcLgeom_towergR(void *p) {
      typedef map<int,PHG4CylinderGeom_Spacalv3::geom_tower> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,PHG4CylinderGeom_Spacalv3::geom_tower>

namespace {
  void TriggerDictionaryInitialization_PHG4CylinderGeom_Spacalv3_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4CylinderGeom_Spacalv3_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4CylinderGeom_Spacalv3;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4CylinderGeom_Spacalv3_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $$Id: PHG4CylinderGeom_Spacalv3.h,v 1.3 2014/08/28 22:18:35 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.3 $$
 * \date $$Date: 2014/08/28 22:18:35 $$
 */
#ifndef G4DETECTORS_PHG4CYLINDERGEOMSPACALV3_H
#define G4DETECTORS_PHG4CYLINDERGEOMSPACALV3_H

#include "PHG4CylinderGeom_Spacalv2.h"

#include <iostream>                     // for operator<<, basic_ostream::op...
#include <map>
#include <string>
#include <utility>  // std::pair, std::make_pair

class PHParameters;

class PHG4CylinderGeom_Spacalv3 : public PHG4CylinderGeom_Spacalv2
{
 public:
  PHG4CylinderGeom_Spacalv3();

  virtual ~PHG4CylinderGeom_Spacalv3();

  virtual void
  identify(std::ostream& os = std::cout) const;
  virtual void
  Print(Option_t* option = "") const;
  virtual void
  SetDefault();

  //! load parameters from PHParameters, which interface to Database/XML/ROOT files
  virtual void
  ImportParameters(const PHParameters& param);

  double
  get_sidewall_outer_torr() const
  {
    return sidewall_outer_torr;
  }

  void
  set_sidewall_outer_torr(double sidewallOuterTorr)
  {
    sidewall_outer_torr = sidewallOuterTorr;
  }

  double
  get_sidewall_thickness() const
  {
    return sidewall_thickness;
  }

  void
  set_sidewall_thickness(double sidewallThickness)
  {
    sidewall_thickness = sidewallThickness;
  }

  std::string
  get_sidewall_mat() const
  {
    return sidewall_mat;
  }

  void
  set_sidewall_mat(const std::string &absorberMat)
  {
    sidewall_mat = absorberMat;
  }

  int get_max_phi_bin_in_sec() const
  {
    return max_phi_bin_in_sec;
  }

  void
  set_max_phi_bin_in_sec(int maxPhiBinInSec)
  {
    max_phi_bin_in_sec = maxPhiBinInSec;
  }

  const std::string& get_divider_mat() const
  {
    return divider_mat;
  }

  void set_divider_mat(const std::string& dividerMat)
  {
    divider_mat = dividerMat;
  }

  double get_divider_width() const
  {
    return divider_width;
  }

  void set_divider_width(double dividerWidth)
  {
    divider_width = dividerWidth;
  }

  class geom_tower
  {
   public:
    geom_tower();
    virtual ~geom_tower()
    {
    }

    int id;
    double pDz;

    double pDy1;
    double pDx1;
    double pDx2;
    double pDy2;
    double pDx3;
    double pDx4;

    double pTheta;
    double pPhi;
    double pAlp1;
    double pAlp2;

    double pRotationAngleX;
    double centralX;
    double centralY;
    double centralZ;

    double ModuleSkinThickness;

    //! number of fiber along final azimuthal direction
    int NFiberX;
    //! number of fiber along final polar direction
    int NFiberY;

    //! number of fiber along final azimuthal direction
    int NSubtowerX;
    //! number of fiber along final polar direction
    int NSubtowerY;

    //! fiber layout (2D index of 0...NFiberX/NFiberY) -> fiber_id
    int
    compose_fiber_id(int index_x, int index_y) const;
    //! fiber_id -> sub tower ID x/azimuthal direction: 0 ... NSubtowerX -1
    int
    get_sub_tower_ID_x(int fiber_id) const;
    //! fiber_id -> sub tower ID y/polar direction: 0 ... NSubtowerY -1
    int
    get_sub_tower_ID_y(int fiber_id) const;
    //! fiber_id -> fraction position in sub tower ID in the x/azimuthal direction, [0-1]
    double
    get_position_fraction_x_in_sub_tower(int fiber_id) const;
    //! fiber_id -> fraction position in sub tower ID in the y/polar direction, [0-1]
    double
    get_position_fraction_y_in_sub_tower(int fiber_id) const;

    //! height of light guide
    double LightguideHeight;
    //! edge length ratio, narrow / wide
    double LightguideTaperRatio;
    //! edge length ratio, narrow / wide
    std::string LightguideMaterial;

    virtual void
    identify(std::ostream& os = std::cout) const;

    //! read via PHParameters
    void
    ImportParameters(const PHParameters& param,
                     const std::string& param_prefix);

    ClassDef(PHG4CylinderGeom_Spacalv3::geom_tower, 3)
  };
  typedef std::map<int, geom_tower> tower_map_t;

  void
  load_demo_sector_tower_map1();
  void
  load_demo_sector_tower_map2();
  void
  load_demo_sector_tower_map_2015_Chris_Cullen_2D_spacal();
  void
  load_demo_sector_tower_map4();

  const tower_map_t&
  get_sector_tower_map() const
  {
    return sector_tower_map;
  }

  //! get approximate radial position of tower
  double
  get_tower_radial_position(const geom_tower& tower) const;

  //! check that all towers has consistent sub-tower divider
  void
  subtower_consistency_check() const;
  //! sub-tower divider along the polar direction
  int get_n_subtower_eta() const;
  //! sub-tower divider along the azimuthal direction
  int get_n_subtower_phi() const;
  //! max tolerance needed for the light guide
  double
  get_max_lightguide_height() const;

  //
  //  void
  //  set_geom_super_tower_map(geom_super_tower_map_t geomSuperTowerMap)
  //  {
  //    geom_super_tower_map = geomSuperTowerMap;
  //  }

  //! compact ID of each fiber in 32bit PHG4Hit::set_scint_id(). Buffer the result for repeated use.
  class scint_id_coder
  {
   public:
    explicit scint_id_coder(int scint_id);
    scint_id_coder(int sector_id, int tower_id, int fiber_id);
    virtual ~scint_id_coder()
    {
    }

    virtual void
    identify(std::ostream& os = std::cout) const
    {
      os << "scint_id_coder with "
         << "scint_ID(" << scint_ID << ") = "
         << "sector_ID(" << sector_ID << "), "
         << "tower_ID(" << tower_ID
         << "), "
         << "fiber_ID(" << fiber_ID << ")" << std::endl;
    }

    int scint_ID;
    int sector_ID;
    int tower_ID;
    int fiber_ID;

    static const int kfiber_bit = 13;  // max 8192 fiber per tower
    static const int ktower_bit = 11;  // max 2048 towers per sector
    static const int ksector_bit = 8;  // max 256 sectors

    ClassDef(PHG4CylinderGeom_Spacalv3::scint_id_coder, 1)
  };

  //! convert tower_ID + sector ID to eta and z bins as in other cylindrical calorimeters
  //! @return: a std::pair of zbin and phibin number
  virtual std::pair<int, int>
  get_tower_z_phi_ID(const int tower_ID, const int sector_ID) const;

 protected:
  double sidewall_thickness;
  double sidewall_outer_torr;
  std::string sidewall_mat;
  int max_phi_bin_in_sec;

  tower_map_t sector_tower_map;

  //! wdith along the approximate radial direction
  double divider_width;
  //! material for divider
  std::string divider_mat;

  ClassDef(PHG4CylinderGeom_Spacalv3, 4)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4CylinderGeom_Spacalv3", payloadCode, "@",
"PHG4CylinderGeom_Spacalv3::geom_tower", payloadCode, "@",
"PHG4CylinderGeom_Spacalv3::scint_id_coder", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4CylinderGeom_Spacalv3_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4CylinderGeom_Spacalv3_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4CylinderGeom_Spacalv3_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4CylinderGeom_Spacalv3_Dict() {
  TriggerDictionaryInitialization_PHG4CylinderGeom_Spacalv3_Dict_Impl();
}
