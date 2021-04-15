// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4CylinderCellv1_Dict
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
#include "PHG4CylinderCellv1.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4CylinderCellv1(void *p = 0);
   static void *newArray_PHG4CylinderCellv1(Long_t size, void *p);
   static void delete_PHG4CylinderCellv1(void *p);
   static void deleteArray_PHG4CylinderCellv1(void *p);
   static void destruct_PHG4CylinderCellv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderCellv1*)
   {
      ::PHG4CylinderCellv1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderCellv1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderCellv1", ::PHG4CylinderCellv1::Class_Version(), "", 22,
                  typeid(::PHG4CylinderCellv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderCellv1::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderCellv1) );
      instance.SetNew(&new_PHG4CylinderCellv1);
      instance.SetNewArray(&newArray_PHG4CylinderCellv1);
      instance.SetDelete(&delete_PHG4CylinderCellv1);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderCellv1);
      instance.SetDestructor(&destruct_PHG4CylinderCellv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderCellv1*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderCellv1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderCellv1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderCellv1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderCellv1::Class_Name()
{
   return "PHG4CylinderCellv1";
}

//______________________________________________________________________________
const char *PHG4CylinderCellv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellv1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderCellv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellv1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderCellv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellv1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderCellv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellv1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4CylinderCellv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderCellv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderCellv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderCellv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4CylinderCellv1(void *p) {
      return  p ? new(p) ::PHG4CylinderCellv1 : new ::PHG4CylinderCellv1;
   }
   static void *newArray_PHG4CylinderCellv1(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4CylinderCellv1[nElements] : new ::PHG4CylinderCellv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4CylinderCellv1(void *p) {
      delete ((::PHG4CylinderCellv1*)p);
   }
   static void deleteArray_PHG4CylinderCellv1(void *p) {
      delete [] ((::PHG4CylinderCellv1*)p);
   }
   static void destruct_PHG4CylinderCellv1(void *p) {
      typedef ::PHG4CylinderCellv1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderCellv1

namespace ROOT {
   static TClass *maplEintcOfloatgR_Dictionary();
   static void maplEintcOfloatgR_TClassManip(TClass*);
   static void *new_maplEintcOfloatgR(void *p = 0);
   static void *newArray_maplEintcOfloatgR(Long_t size, void *p);
   static void delete_maplEintcOfloatgR(void *p);
   static void deleteArray_maplEintcOfloatgR(void *p);
   static void destruct_maplEintcOfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,float>*)
   {
      map<int,float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,float>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,float>", -2, "map", 100,
                  typeid(map<int,float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,float>) );
      instance.SetNew(&new_maplEintcOfloatgR);
      instance.SetNewArray(&newArray_maplEintcOfloatgR);
      instance.SetDelete(&delete_maplEintcOfloatgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOfloatgR);
      instance.SetDestructor(&destruct_maplEintcOfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,float> >()));

      ::ROOT::AddClassAlternate("map<int,float>","std::map<int, float, std::less<int>, std::allocator<std::pair<int const, float> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,float>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,float>*)0x0)->GetClass();
      maplEintcOfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,float> : new map<int,float>;
   }
   static void *newArray_maplEintcOfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,float>[nElements] : new map<int,float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOfloatgR(void *p) {
      delete ((map<int,float>*)p);
   }
   static void deleteArray_maplEintcOfloatgR(void *p) {
      delete [] ((map<int,float>*)p);
   }
   static void destruct_maplEintcOfloatgR(void *p) {
      typedef map<int,float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,float>

namespace ROOT {
   static TClass *maplEULong64_tcOfloatgR_Dictionary();
   static void maplEULong64_tcOfloatgR_TClassManip(TClass*);
   static void *new_maplEULong64_tcOfloatgR(void *p = 0);
   static void *newArray_maplEULong64_tcOfloatgR(Long_t size, void *p);
   static void delete_maplEULong64_tcOfloatgR(void *p);
   static void deleteArray_maplEULong64_tcOfloatgR(void *p);
   static void destruct_maplEULong64_tcOfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<ULong64_t,float>*)
   {
      map<ULong64_t,float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<ULong64_t,float>));
      static ::ROOT::TGenericClassInfo 
         instance("map<ULong64_t,float>", -2, "map", 100,
                  typeid(map<ULong64_t,float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEULong64_tcOfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(map<ULong64_t,float>) );
      instance.SetNew(&new_maplEULong64_tcOfloatgR);
      instance.SetNewArray(&newArray_maplEULong64_tcOfloatgR);
      instance.SetDelete(&delete_maplEULong64_tcOfloatgR);
      instance.SetDeleteArray(&deleteArray_maplEULong64_tcOfloatgR);
      instance.SetDestructor(&destruct_maplEULong64_tcOfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<ULong64_t,float> >()));

      ::ROOT::AddClassAlternate("map<ULong64_t,float>","std::map<unsigned long long, float, std::less<unsigned long long>, std::allocator<std::pair<unsigned long long const, float> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<ULong64_t,float>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEULong64_tcOfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<ULong64_t,float>*)0x0)->GetClass();
      maplEULong64_tcOfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEULong64_tcOfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEULong64_tcOfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<ULong64_t,float> : new map<ULong64_t,float>;
   }
   static void *newArray_maplEULong64_tcOfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<ULong64_t,float>[nElements] : new map<ULong64_t,float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEULong64_tcOfloatgR(void *p) {
      delete ((map<ULong64_t,float>*)p);
   }
   static void deleteArray_maplEULong64_tcOfloatgR(void *p) {
      delete [] ((map<ULong64_t,float>*)p);
   }
   static void destruct_maplEULong64_tcOfloatgR(void *p) {
      typedef map<ULong64_t,float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<ULong64_t,float>

namespace {
  void TriggerDictionaryInitialization_PHG4CylinderCellv1_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4CylinderCellv1_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4CylinderCellv1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4CylinderCellv1_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERCELLV1_H
#define G4DETECTORS_PHG4CYLINDERCELLV1_H

#include "PHG4CylinderCell.h"

#include "PHG4Cell.h"              // for PHG4Cell::EdepMap, PHG4Cell::Showe...
#include "PHG4CylinderCellDefs.h"  // for keytype

#include <g4main/PHG4HitDefs.h>

#include <iostream>                // for cout, ostream
#include <map>
#include <utility>                 // for make_pair

class PHG4CylinderCellv1 : public PHG4CylinderCell
{
 public:

  PHG4CylinderCellv1();
  virtual ~PHG4CylinderCellv1(){}

  void identify(std::ostream& os = std::cout) const;

  EdepConstRange get_g4hits()
  {return std::make_pair(edeps.begin(), edeps.end());}
  void add_edep(const PHG4HitDefs::keytype g4hitid, const float edep);
  void add_edep(const PHG4HitDefs::keytype g4hitid, const float edep, const float light_yield);
  
  ShowerEdepConstRange get_g4showers()
  {return std::make_pair(showeredeps.begin(), showeredeps.end());}
  void add_shower_edep(const int g4showerid, const float edep);
  
  void set_cell_id(const PHG4CylinderCellDefs::keytype id) {cellid = id;}
  void set_layer(const unsigned int i) {layer = i;}
  double get_edep() const;
  unsigned int get_layer() const {return layer;}
  PHG4CylinderCellDefs::keytype get_cell_id() const {return cellid;}
  int get_binz() const {return binz;}
  int get_binphi() const {return binphi;}
  int get_bineta() const {return get_binz();}
  float get_light_yield() const  {return light_yield;}


  void set_zbin(const int i) {binz = i;}
  void set_etabin(const int i) {set_zbin(i);}
  void set_phibin(const int i) {binphi = i;}
  void set_light_yield(const float lightYield)  {    light_yield = lightYield;  }

 protected:

  unsigned int layer;
  PHG4CylinderCellDefs::keytype cellid;
  int binz;
  int binphi;
  EdepMap edeps;
  ShowerEdepMap showeredeps;
  float light_yield;
   
  ClassDef(PHG4CylinderCellv1,2)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4CylinderCellv1", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4CylinderCellv1_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4CylinderCellv1_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4CylinderCellv1_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4CylinderCellv1_Dict() {
  TriggerDictionaryInitialization_PHG4CylinderCellv1_Dict_Impl();
}
