// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4BlockCellGeomContainer_Dict
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
#include "PHG4BlockCellGeomContainer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4BlockCellGeomContainer(void *p = 0);
   static void *newArray_PHG4BlockCellGeomContainer(Long_t size, void *p);
   static void delete_PHG4BlockCellGeomContainer(void *p);
   static void deleteArray_PHG4BlockCellGeomContainer(void *p);
   static void destruct_PHG4BlockCellGeomContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4BlockCellGeomContainer*)
   {
      ::PHG4BlockCellGeomContainer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4BlockCellGeomContainer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4BlockCellGeomContainer", ::PHG4BlockCellGeomContainer::Class_Version(), "", 19,
                  typeid(::PHG4BlockCellGeomContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4BlockCellGeomContainer::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4BlockCellGeomContainer) );
      instance.SetNew(&new_PHG4BlockCellGeomContainer);
      instance.SetNewArray(&newArray_PHG4BlockCellGeomContainer);
      instance.SetDelete(&delete_PHG4BlockCellGeomContainer);
      instance.SetDeleteArray(&deleteArray_PHG4BlockCellGeomContainer);
      instance.SetDestructor(&destruct_PHG4BlockCellGeomContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4BlockCellGeomContainer*)
   {
      return GenerateInitInstanceLocal((::PHG4BlockCellGeomContainer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4BlockCellGeomContainer*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4BlockCellGeomContainer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4BlockCellGeomContainer::Class_Name()
{
   return "PHG4BlockCellGeomContainer";
}

//______________________________________________________________________________
const char *PHG4BlockCellGeomContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockCellGeomContainer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4BlockCellGeomContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockCellGeomContainer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4BlockCellGeomContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockCellGeomContainer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4BlockCellGeomContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockCellGeomContainer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4BlockCellGeomContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4BlockCellGeomContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4BlockCellGeomContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4BlockCellGeomContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4BlockCellGeomContainer(void *p) {
      return  p ? new(p) ::PHG4BlockCellGeomContainer : new ::PHG4BlockCellGeomContainer;
   }
   static void *newArray_PHG4BlockCellGeomContainer(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4BlockCellGeomContainer[nElements] : new ::PHG4BlockCellGeomContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4BlockCellGeomContainer(void *p) {
      delete ((::PHG4BlockCellGeomContainer*)p);
   }
   static void deleteArray_PHG4BlockCellGeomContainer(void *p) {
      delete [] ((::PHG4BlockCellGeomContainer*)p);
   }
   static void destruct_PHG4BlockCellGeomContainer(void *p) {
      typedef ::PHG4BlockCellGeomContainer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4BlockCellGeomContainer

namespace ROOT {
   static TClass *maplEintcOPHG4BlockCellGeommUgR_Dictionary();
   static void maplEintcOPHG4BlockCellGeommUgR_TClassManip(TClass*);
   static void *new_maplEintcOPHG4BlockCellGeommUgR(void *p = 0);
   static void *newArray_maplEintcOPHG4BlockCellGeommUgR(Long_t size, void *p);
   static void delete_maplEintcOPHG4BlockCellGeommUgR(void *p);
   static void deleteArray_maplEintcOPHG4BlockCellGeommUgR(void *p);
   static void destruct_maplEintcOPHG4BlockCellGeommUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,PHG4BlockCellGeom*>*)
   {
      map<int,PHG4BlockCellGeom*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,PHG4BlockCellGeom*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,PHG4BlockCellGeom*>", -2, "map", 100,
                  typeid(map<int,PHG4BlockCellGeom*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOPHG4BlockCellGeommUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,PHG4BlockCellGeom*>) );
      instance.SetNew(&new_maplEintcOPHG4BlockCellGeommUgR);
      instance.SetNewArray(&newArray_maplEintcOPHG4BlockCellGeommUgR);
      instance.SetDelete(&delete_maplEintcOPHG4BlockCellGeommUgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOPHG4BlockCellGeommUgR);
      instance.SetDestructor(&destruct_maplEintcOPHG4BlockCellGeommUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,PHG4BlockCellGeom*> >()));

      ::ROOT::AddClassAlternate("map<int,PHG4BlockCellGeom*>","std::map<int, PHG4BlockCellGeom*, std::less<int>, std::allocator<std::pair<int const, PHG4BlockCellGeom*> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,PHG4BlockCellGeom*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOPHG4BlockCellGeommUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,PHG4BlockCellGeom*>*)0x0)->GetClass();
      maplEintcOPHG4BlockCellGeommUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOPHG4BlockCellGeommUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOPHG4BlockCellGeommUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,PHG4BlockCellGeom*> : new map<int,PHG4BlockCellGeom*>;
   }
   static void *newArray_maplEintcOPHG4BlockCellGeommUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,PHG4BlockCellGeom*>[nElements] : new map<int,PHG4BlockCellGeom*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOPHG4BlockCellGeommUgR(void *p) {
      delete ((map<int,PHG4BlockCellGeom*>*)p);
   }
   static void deleteArray_maplEintcOPHG4BlockCellGeommUgR(void *p) {
      delete [] ((map<int,PHG4BlockCellGeom*>*)p);
   }
   static void destruct_maplEintcOPHG4BlockCellGeommUgR(void *p) {
      typedef map<int,PHG4BlockCellGeom*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,PHG4BlockCellGeom*>

namespace {
  void TriggerDictionaryInitialization_PHG4BlockCellGeomContainer_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4BlockCellGeomContainer_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4BlockCellGeomContainer;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4BlockCellGeomContainer_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BLOCKCELLGEOMCONTAINER_H
#define G4DETECTORS_PHG4BLOCKCELLGEOMCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream
#include <map>
#include <utility>           // for make_pair, pair

class PHG4BlockCellGeom;

class PHG4BlockCellGeomContainer: public PHObject
{
 public:
  typedef std::map<int, PHG4BlockCellGeom*> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  PHG4BlockCellGeomContainer(){}
  virtual ~PHG4BlockCellGeomContainer();

  void identify(std::ostream& os = std::cout) const;

  int AddLayerCellGeom(const int i, PHG4BlockCellGeom *mygeom);
  int AddLayerCellGeom(PHG4BlockCellGeom *mygeom);
  PHG4BlockCellGeom *GetLayerCellGeom(const int i);
  int get_NLayers() const {return layergeoms.size();}
  std::pair<std::map<int,PHG4BlockCellGeom *>::const_iterator, std::map<int,PHG4BlockCellGeom *>::const_iterator> get_begin_end() const {return std::make_pair(layergeoms.begin(), layergeoms.end());}

 protected:
  std::map<int,PHG4BlockCellGeom *> layergeoms ;
  ClassDef(PHG4BlockCellGeomContainer,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4BlockCellGeomContainer", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4BlockCellGeomContainer_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4BlockCellGeomContainer_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4BlockCellGeomContainer_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4BlockCellGeomContainer_Dict() {
  TriggerDictionaryInitialization_PHG4BlockCellGeomContainer_Dict_Impl();
}
