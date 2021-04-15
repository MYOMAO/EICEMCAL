// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4CylinderGeomContainer_Dict
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
#include "PHG4CylinderGeomContainer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4CylinderGeomContainer(void *p = 0);
   static void *newArray_PHG4CylinderGeomContainer(Long_t size, void *p);
   static void delete_PHG4CylinderGeomContainer(void *p);
   static void deleteArray_PHG4CylinderGeomContainer(void *p);
   static void destruct_PHG4CylinderGeomContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderGeomContainer*)
   {
      ::PHG4CylinderGeomContainer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderGeomContainer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderGeomContainer", ::PHG4CylinderGeomContainer::Class_Version(), "", 19,
                  typeid(::PHG4CylinderGeomContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderGeomContainer::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderGeomContainer) );
      instance.SetNew(&new_PHG4CylinderGeomContainer);
      instance.SetNewArray(&newArray_PHG4CylinderGeomContainer);
      instance.SetDelete(&delete_PHG4CylinderGeomContainer);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderGeomContainer);
      instance.SetDestructor(&destruct_PHG4CylinderGeomContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderGeomContainer*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderGeomContainer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderGeomContainer*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderGeomContainer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderGeomContainer::Class_Name()
{
   return "PHG4CylinderGeomContainer";
}

//______________________________________________________________________________
const char *PHG4CylinderGeomContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeomContainer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderGeomContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeomContainer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderGeomContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeomContainer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderGeomContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderGeomContainer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4CylinderGeomContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderGeomContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderGeomContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderGeomContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4CylinderGeomContainer(void *p) {
      return  p ? new(p) ::PHG4CylinderGeomContainer : new ::PHG4CylinderGeomContainer;
   }
   static void *newArray_PHG4CylinderGeomContainer(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4CylinderGeomContainer[nElements] : new ::PHG4CylinderGeomContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4CylinderGeomContainer(void *p) {
      delete ((::PHG4CylinderGeomContainer*)p);
   }
   static void deleteArray_PHG4CylinderGeomContainer(void *p) {
      delete [] ((::PHG4CylinderGeomContainer*)p);
   }
   static void destruct_PHG4CylinderGeomContainer(void *p) {
      typedef ::PHG4CylinderGeomContainer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderGeomContainer

namespace ROOT {
   static TClass *maplEintcOPHG4CylinderGeommUgR_Dictionary();
   static void maplEintcOPHG4CylinderGeommUgR_TClassManip(TClass*);
   static void *new_maplEintcOPHG4CylinderGeommUgR(void *p = 0);
   static void *newArray_maplEintcOPHG4CylinderGeommUgR(Long_t size, void *p);
   static void delete_maplEintcOPHG4CylinderGeommUgR(void *p);
   static void deleteArray_maplEintcOPHG4CylinderGeommUgR(void *p);
   static void destruct_maplEintcOPHG4CylinderGeommUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,PHG4CylinderGeom*>*)
   {
      map<int,PHG4CylinderGeom*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,PHG4CylinderGeom*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,PHG4CylinderGeom*>", -2, "map", 100,
                  typeid(map<int,PHG4CylinderGeom*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOPHG4CylinderGeommUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,PHG4CylinderGeom*>) );
      instance.SetNew(&new_maplEintcOPHG4CylinderGeommUgR);
      instance.SetNewArray(&newArray_maplEintcOPHG4CylinderGeommUgR);
      instance.SetDelete(&delete_maplEintcOPHG4CylinderGeommUgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOPHG4CylinderGeommUgR);
      instance.SetDestructor(&destruct_maplEintcOPHG4CylinderGeommUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,PHG4CylinderGeom*> >()));

      ::ROOT::AddClassAlternate("map<int,PHG4CylinderGeom*>","std::map<int, PHG4CylinderGeom*, std::less<int>, std::allocator<std::pair<int const, PHG4CylinderGeom*> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,PHG4CylinderGeom*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOPHG4CylinderGeommUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,PHG4CylinderGeom*>*)0x0)->GetClass();
      maplEintcOPHG4CylinderGeommUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOPHG4CylinderGeommUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOPHG4CylinderGeommUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,PHG4CylinderGeom*> : new map<int,PHG4CylinderGeom*>;
   }
   static void *newArray_maplEintcOPHG4CylinderGeommUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,PHG4CylinderGeom*>[nElements] : new map<int,PHG4CylinderGeom*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOPHG4CylinderGeommUgR(void *p) {
      delete ((map<int,PHG4CylinderGeom*>*)p);
   }
   static void deleteArray_maplEintcOPHG4CylinderGeommUgR(void *p) {
      delete [] ((map<int,PHG4CylinderGeom*>*)p);
   }
   static void destruct_maplEintcOPHG4CylinderGeommUgR(void *p) {
      typedef map<int,PHG4CylinderGeom*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,PHG4CylinderGeom*>

namespace {
  void TriggerDictionaryInitialization_PHG4CylinderGeomContainer_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4CylinderGeomContainer_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4CylinderGeomContainer;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4CylinderGeomContainer_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERGEOMCONTAINER_H
#define G4DETECTORS_PHG4CYLINDERGEOMCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream
#include <map>
#include <utility>           // for make_pair, pair

class PHG4CylinderGeom;

class PHG4CylinderGeomContainer: public PHObject
{
 public:
  typedef std::map<int,PHG4CylinderGeom *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  PHG4CylinderGeomContainer();
  virtual ~PHG4CylinderGeomContainer();

  void identify(std::ostream& os = std::cout) const;

  int AddLayerGeom(const int i, PHG4CylinderGeom *mygeom);
  int AddLayerGeom(PHG4CylinderGeom *mygeom);
  PHG4CylinderGeom *GetLayerGeom(const int i);
  PHG4CylinderGeom *GetFirstLayerGeom();
  int get_NLayers() const {return layergeoms.size();}
  std::pair<std::map<int,PHG4CylinderGeom *>::const_iterator, std::map<int,PHG4CylinderGeom *>::const_iterator> get_begin_end() const {return std::make_pair(layergeoms.begin(), layergeoms.end());}

 protected:
  std::map<int,PHG4CylinderGeom *> layergeoms ;
  float magfield;
  ClassDef(PHG4CylinderGeomContainer,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4CylinderGeomContainer", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4CylinderGeomContainer_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4CylinderGeomContainer_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4CylinderGeomContainer_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4CylinderGeomContainer_Dict() {
  TriggerDictionaryInitialization_PHG4CylinderGeomContainer_Dict_Impl();
}
