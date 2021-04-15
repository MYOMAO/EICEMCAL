// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4BlockGeomContainer_Dict
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
#include "PHG4BlockGeomContainer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4BlockGeomContainer(void *p = 0);
   static void *newArray_PHG4BlockGeomContainer(Long_t size, void *p);
   static void delete_PHG4BlockGeomContainer(void *p);
   static void deleteArray_PHG4BlockGeomContainer(void *p);
   static void destruct_PHG4BlockGeomContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4BlockGeomContainer*)
   {
      ::PHG4BlockGeomContainer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4BlockGeomContainer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4BlockGeomContainer", ::PHG4BlockGeomContainer::Class_Version(), "", 19,
                  typeid(::PHG4BlockGeomContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4BlockGeomContainer::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4BlockGeomContainer) );
      instance.SetNew(&new_PHG4BlockGeomContainer);
      instance.SetNewArray(&newArray_PHG4BlockGeomContainer);
      instance.SetDelete(&delete_PHG4BlockGeomContainer);
      instance.SetDeleteArray(&deleteArray_PHG4BlockGeomContainer);
      instance.SetDestructor(&destruct_PHG4BlockGeomContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4BlockGeomContainer*)
   {
      return GenerateInitInstanceLocal((::PHG4BlockGeomContainer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4BlockGeomContainer*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4BlockGeomContainer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4BlockGeomContainer::Class_Name()
{
   return "PHG4BlockGeomContainer";
}

//______________________________________________________________________________
const char *PHG4BlockGeomContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockGeomContainer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4BlockGeomContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockGeomContainer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4BlockGeomContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockGeomContainer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4BlockGeomContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockGeomContainer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4BlockGeomContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4BlockGeomContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4BlockGeomContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4BlockGeomContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4BlockGeomContainer(void *p) {
      return  p ? new(p) ::PHG4BlockGeomContainer : new ::PHG4BlockGeomContainer;
   }
   static void *newArray_PHG4BlockGeomContainer(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4BlockGeomContainer[nElements] : new ::PHG4BlockGeomContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4BlockGeomContainer(void *p) {
      delete ((::PHG4BlockGeomContainer*)p);
   }
   static void deleteArray_PHG4BlockGeomContainer(void *p) {
      delete [] ((::PHG4BlockGeomContainer*)p);
   }
   static void destruct_PHG4BlockGeomContainer(void *p) {
      typedef ::PHG4BlockGeomContainer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4BlockGeomContainer

namespace ROOT {
   static TClass *maplEintcOPHG4BlockGeommUgR_Dictionary();
   static void maplEintcOPHG4BlockGeommUgR_TClassManip(TClass*);
   static void *new_maplEintcOPHG4BlockGeommUgR(void *p = 0);
   static void *newArray_maplEintcOPHG4BlockGeommUgR(Long_t size, void *p);
   static void delete_maplEintcOPHG4BlockGeommUgR(void *p);
   static void deleteArray_maplEintcOPHG4BlockGeommUgR(void *p);
   static void destruct_maplEintcOPHG4BlockGeommUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,PHG4BlockGeom*>*)
   {
      map<int,PHG4BlockGeom*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,PHG4BlockGeom*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,PHG4BlockGeom*>", -2, "map", 100,
                  typeid(map<int,PHG4BlockGeom*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOPHG4BlockGeommUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,PHG4BlockGeom*>) );
      instance.SetNew(&new_maplEintcOPHG4BlockGeommUgR);
      instance.SetNewArray(&newArray_maplEintcOPHG4BlockGeommUgR);
      instance.SetDelete(&delete_maplEintcOPHG4BlockGeommUgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOPHG4BlockGeommUgR);
      instance.SetDestructor(&destruct_maplEintcOPHG4BlockGeommUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,PHG4BlockGeom*> >()));

      ::ROOT::AddClassAlternate("map<int,PHG4BlockGeom*>","std::map<int, PHG4BlockGeom*, std::less<int>, std::allocator<std::pair<int const, PHG4BlockGeom*> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,PHG4BlockGeom*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOPHG4BlockGeommUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,PHG4BlockGeom*>*)0x0)->GetClass();
      maplEintcOPHG4BlockGeommUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOPHG4BlockGeommUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOPHG4BlockGeommUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,PHG4BlockGeom*> : new map<int,PHG4BlockGeom*>;
   }
   static void *newArray_maplEintcOPHG4BlockGeommUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,PHG4BlockGeom*>[nElements] : new map<int,PHG4BlockGeom*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOPHG4BlockGeommUgR(void *p) {
      delete ((map<int,PHG4BlockGeom*>*)p);
   }
   static void deleteArray_maplEintcOPHG4BlockGeommUgR(void *p) {
      delete [] ((map<int,PHG4BlockGeom*>*)p);
   }
   static void destruct_maplEintcOPHG4BlockGeommUgR(void *p) {
      typedef map<int,PHG4BlockGeom*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,PHG4BlockGeom*>

namespace {
  void TriggerDictionaryInitialization_PHG4BlockGeomContainer_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4BlockGeomContainer_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4BlockGeomContainer;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4BlockGeomContainer_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BLOCKGEOMCONTAINER_H
#define G4DETECTORS_PHG4BLOCKGEOMCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream
#include <map>
#include <utility>           // for make_pair, pair

class PHG4BlockGeom;

class PHG4BlockGeomContainer: public PHObject
{
 public:
  typedef std::map<int,PHG4BlockGeom *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  PHG4BlockGeomContainer();
  virtual ~PHG4BlockGeomContainer();

  void identify(std::ostream& os = std::cout) const;

  int AddLayerGeom(const int i, PHG4BlockGeom *mygeom);
  int AddLayerGeom(PHG4BlockGeom *mygeom);
  PHG4BlockGeom *GetLayerGeom(const int i);
  int get_NLayers() const {return _layergeoms.size();}
  std::pair<std::map<int,PHG4BlockGeom *>::const_iterator, std::map<int,PHG4BlockGeom *>::const_iterator> get_begin_end() const {return std::make_pair(_layergeoms.begin(), _layergeoms.end());}

 protected:
  std::map<int,PHG4BlockGeom *> _layergeoms;
  float _magfield;

  ClassDef(PHG4BlockGeomContainer,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4BlockGeomContainer", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4BlockGeomContainer_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4BlockGeomContainer_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4BlockGeomContainer_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4BlockGeomContainer_Dict() {
  TriggerDictionaryInitialization_PHG4BlockGeomContainer_Dict_Impl();
}
