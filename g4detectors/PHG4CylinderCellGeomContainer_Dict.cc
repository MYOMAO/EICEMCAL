// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4CylinderCellGeomContainer_Dict
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
#include "PHG4CylinderCellGeomContainer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4CylinderCellGeomContainer(void *p = 0);
   static void *newArray_PHG4CylinderCellGeomContainer(Long_t size, void *p);
   static void delete_PHG4CylinderCellGeomContainer(void *p);
   static void deleteArray_PHG4CylinderCellGeomContainer(void *p);
   static void destruct_PHG4CylinderCellGeomContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderCellGeomContainer*)
   {
      ::PHG4CylinderCellGeomContainer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderCellGeomContainer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderCellGeomContainer", ::PHG4CylinderCellGeomContainer::Class_Version(), "", 19,
                  typeid(::PHG4CylinderCellGeomContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderCellGeomContainer::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderCellGeomContainer) );
      instance.SetNew(&new_PHG4CylinderCellGeomContainer);
      instance.SetNewArray(&newArray_PHG4CylinderCellGeomContainer);
      instance.SetDelete(&delete_PHG4CylinderCellGeomContainer);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderCellGeomContainer);
      instance.SetDestructor(&destruct_PHG4CylinderCellGeomContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderCellGeomContainer*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderCellGeomContainer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderCellGeomContainer*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderCellGeomContainer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderCellGeomContainer::Class_Name()
{
   return "PHG4CylinderCellGeomContainer";
}

//______________________________________________________________________________
const char *PHG4CylinderCellGeomContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellGeomContainer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderCellGeomContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellGeomContainer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderCellGeomContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellGeomContainer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderCellGeomContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellGeomContainer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4CylinderCellGeomContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderCellGeomContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderCellGeomContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderCellGeomContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4CylinderCellGeomContainer(void *p) {
      return  p ? new(p) ::PHG4CylinderCellGeomContainer : new ::PHG4CylinderCellGeomContainer;
   }
   static void *newArray_PHG4CylinderCellGeomContainer(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4CylinderCellGeomContainer[nElements] : new ::PHG4CylinderCellGeomContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4CylinderCellGeomContainer(void *p) {
      delete ((::PHG4CylinderCellGeomContainer*)p);
   }
   static void deleteArray_PHG4CylinderCellGeomContainer(void *p) {
      delete [] ((::PHG4CylinderCellGeomContainer*)p);
   }
   static void destruct_PHG4CylinderCellGeomContainer(void *p) {
      typedef ::PHG4CylinderCellGeomContainer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderCellGeomContainer

namespace ROOT {
   static TClass *maplEintcOPHG4CylinderCellGeommUgR_Dictionary();
   static void maplEintcOPHG4CylinderCellGeommUgR_TClassManip(TClass*);
   static void *new_maplEintcOPHG4CylinderCellGeommUgR(void *p = 0);
   static void *newArray_maplEintcOPHG4CylinderCellGeommUgR(Long_t size, void *p);
   static void delete_maplEintcOPHG4CylinderCellGeommUgR(void *p);
   static void deleteArray_maplEintcOPHG4CylinderCellGeommUgR(void *p);
   static void destruct_maplEintcOPHG4CylinderCellGeommUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,PHG4CylinderCellGeom*>*)
   {
      map<int,PHG4CylinderCellGeom*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,PHG4CylinderCellGeom*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,PHG4CylinderCellGeom*>", -2, "map", 100,
                  typeid(map<int,PHG4CylinderCellGeom*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOPHG4CylinderCellGeommUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,PHG4CylinderCellGeom*>) );
      instance.SetNew(&new_maplEintcOPHG4CylinderCellGeommUgR);
      instance.SetNewArray(&newArray_maplEintcOPHG4CylinderCellGeommUgR);
      instance.SetDelete(&delete_maplEintcOPHG4CylinderCellGeommUgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOPHG4CylinderCellGeommUgR);
      instance.SetDestructor(&destruct_maplEintcOPHG4CylinderCellGeommUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,PHG4CylinderCellGeom*> >()));

      ::ROOT::AddClassAlternate("map<int,PHG4CylinderCellGeom*>","std::map<int, PHG4CylinderCellGeom*, std::less<int>, std::allocator<std::pair<int const, PHG4CylinderCellGeom*> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,PHG4CylinderCellGeom*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOPHG4CylinderCellGeommUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,PHG4CylinderCellGeom*>*)0x0)->GetClass();
      maplEintcOPHG4CylinderCellGeommUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOPHG4CylinderCellGeommUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOPHG4CylinderCellGeommUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,PHG4CylinderCellGeom*> : new map<int,PHG4CylinderCellGeom*>;
   }
   static void *newArray_maplEintcOPHG4CylinderCellGeommUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,PHG4CylinderCellGeom*>[nElements] : new map<int,PHG4CylinderCellGeom*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOPHG4CylinderCellGeommUgR(void *p) {
      delete ((map<int,PHG4CylinderCellGeom*>*)p);
   }
   static void deleteArray_maplEintcOPHG4CylinderCellGeommUgR(void *p) {
      delete [] ((map<int,PHG4CylinderCellGeom*>*)p);
   }
   static void destruct_maplEintcOPHG4CylinderCellGeommUgR(void *p) {
      typedef map<int,PHG4CylinderCellGeom*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,PHG4CylinderCellGeom*>

namespace {
  void TriggerDictionaryInitialization_PHG4CylinderCellGeomContainer_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4CylinderCellGeomContainer_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4CylinderCellGeomContainer;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4CylinderCellGeomContainer_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERCELLGEOMCONTAINER_H
#define G4DETECTORS_PHG4CYLINDERCELLGEOMCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream
#include <map>
#include <utility>           // for make_pair, pair

class PHG4CylinderCellGeom;

class PHG4CylinderCellGeomContainer: public PHObject
{
 public:
  typedef std::map<int, PHG4CylinderCellGeom*> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  PHG4CylinderCellGeomContainer(){}
  virtual ~PHG4CylinderCellGeomContainer();

  void identify(std::ostream& os = std::cout) const;

  int AddLayerCellGeom(const int i, PHG4CylinderCellGeom *mygeom);
  int AddLayerCellGeom(PHG4CylinderCellGeom *mygeom);
  PHG4CylinderCellGeom *GetLayerCellGeom(const int i);
  PHG4CylinderCellGeom *GetFirstLayerCellGeom();
  int get_NLayers() const {return layergeoms.size();}
  std::pair<std::map<int,PHG4CylinderCellGeom *>::const_iterator, std::map<int,PHG4CylinderCellGeom *>::const_iterator> get_begin_end() const {return std::make_pair(layergeoms.begin(), layergeoms.end());}

 protected:
  std::map<int,PHG4CylinderCellGeom *> layergeoms ;
  ClassDef(PHG4CylinderCellGeomContainer,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4CylinderCellGeomContainer", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4CylinderCellGeomContainer_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4CylinderCellGeomContainer_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4CylinderCellGeomContainer_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4CylinderCellGeomContainer_Dict() {
  TriggerDictionaryInitialization_PHG4CylinderCellGeomContainer_Dict_Impl();
}
