// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4CylinderCellContainer_Dict
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
#include "PHG4CylinderCellContainer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4CylinderCellContainer(void *p = 0);
   static void *newArray_PHG4CylinderCellContainer(Long_t size, void *p);
   static void delete_PHG4CylinderCellContainer(void *p);
   static void deleteArray_PHG4CylinderCellContainer(void *p);
   static void destruct_PHG4CylinderCellContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderCellContainer*)
   {
      ::PHG4CylinderCellContainer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderCellContainer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderCellContainer", ::PHG4CylinderCellContainer::Class_Version(), "", 22,
                  typeid(::PHG4CylinderCellContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderCellContainer::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderCellContainer) );
      instance.SetNew(&new_PHG4CylinderCellContainer);
      instance.SetNewArray(&newArray_PHG4CylinderCellContainer);
      instance.SetDelete(&delete_PHG4CylinderCellContainer);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderCellContainer);
      instance.SetDestructor(&destruct_PHG4CylinderCellContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderCellContainer*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderCellContainer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderCellContainer*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderCellContainer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderCellContainer::Class_Name()
{
   return "PHG4CylinderCellContainer";
}

//______________________________________________________________________________
const char *PHG4CylinderCellContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellContainer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderCellContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellContainer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderCellContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellContainer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderCellContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellContainer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4CylinderCellContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderCellContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderCellContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderCellContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4CylinderCellContainer(void *p) {
      return  p ? new(p) ::PHG4CylinderCellContainer : new ::PHG4CylinderCellContainer;
   }
   static void *newArray_PHG4CylinderCellContainer(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4CylinderCellContainer[nElements] : new ::PHG4CylinderCellContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4CylinderCellContainer(void *p) {
      delete ((::PHG4CylinderCellContainer*)p);
   }
   static void deleteArray_PHG4CylinderCellContainer(void *p) {
      delete [] ((::PHG4CylinderCellContainer*)p);
   }
   static void destruct_PHG4CylinderCellContainer(void *p) {
      typedef ::PHG4CylinderCellContainer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderCellContainer

namespace ROOT {
   static TClass *setlEintgR_Dictionary();
   static void setlEintgR_TClassManip(TClass*);
   static void *new_setlEintgR(void *p = 0);
   static void *newArray_setlEintgR(Long_t size, void *p);
   static void delete_setlEintgR(void *p);
   static void deleteArray_setlEintgR(void *p);
   static void destruct_setlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const set<int>*)
   {
      set<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(set<int>));
      static ::ROOT::TGenericClassInfo 
         instance("set<int>", -2, "set", 94,
                  typeid(set<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &setlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(set<int>) );
      instance.SetNew(&new_setlEintgR);
      instance.SetNewArray(&newArray_setlEintgR);
      instance.SetDelete(&delete_setlEintgR);
      instance.SetDeleteArray(&deleteArray_setlEintgR);
      instance.SetDestructor(&destruct_setlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Insert< set<int> >()));

      ::ROOT::AddClassAlternate("set<int>","std::set<int, std::less<int>, std::allocator<int> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const set<int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *setlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const set<int>*)0x0)->GetClass();
      setlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void setlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_setlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<int> : new set<int>;
   }
   static void *newArray_setlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<int>[nElements] : new set<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_setlEintgR(void *p) {
      delete ((set<int>*)p);
   }
   static void deleteArray_setlEintgR(void *p) {
      delete [] ((set<int>*)p);
   }
   static void destruct_setlEintgR(void *p) {
      typedef set<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class set<int>

namespace ROOT {
   static TClass *maplEunsignedsPintcOPHG4CylinderCellmUgR_Dictionary();
   static void maplEunsignedsPintcOPHG4CylinderCellmUgR_TClassManip(TClass*);
   static void *new_maplEunsignedsPintcOPHG4CylinderCellmUgR(void *p = 0);
   static void *newArray_maplEunsignedsPintcOPHG4CylinderCellmUgR(Long_t size, void *p);
   static void delete_maplEunsignedsPintcOPHG4CylinderCellmUgR(void *p);
   static void deleteArray_maplEunsignedsPintcOPHG4CylinderCellmUgR(void *p);
   static void destruct_maplEunsignedsPintcOPHG4CylinderCellmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<unsigned int,PHG4CylinderCell*>*)
   {
      map<unsigned int,PHG4CylinderCell*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<unsigned int,PHG4CylinderCell*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<unsigned int,PHG4CylinderCell*>", -2, "map", 100,
                  typeid(map<unsigned int,PHG4CylinderCell*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEunsignedsPintcOPHG4CylinderCellmUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<unsigned int,PHG4CylinderCell*>) );
      instance.SetNew(&new_maplEunsignedsPintcOPHG4CylinderCellmUgR);
      instance.SetNewArray(&newArray_maplEunsignedsPintcOPHG4CylinderCellmUgR);
      instance.SetDelete(&delete_maplEunsignedsPintcOPHG4CylinderCellmUgR);
      instance.SetDeleteArray(&deleteArray_maplEunsignedsPintcOPHG4CylinderCellmUgR);
      instance.SetDestructor(&destruct_maplEunsignedsPintcOPHG4CylinderCellmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<unsigned int,PHG4CylinderCell*> >()));

      ::ROOT::AddClassAlternate("map<unsigned int,PHG4CylinderCell*>","std::map<unsigned int, PHG4CylinderCell*, std::less<unsigned int>, std::allocator<std::pair<unsigned int const, PHG4CylinderCell*> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<unsigned int,PHG4CylinderCell*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEunsignedsPintcOPHG4CylinderCellmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<unsigned int,PHG4CylinderCell*>*)0x0)->GetClass();
      maplEunsignedsPintcOPHG4CylinderCellmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEunsignedsPintcOPHG4CylinderCellmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEunsignedsPintcOPHG4CylinderCellmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<unsigned int,PHG4CylinderCell*> : new map<unsigned int,PHG4CylinderCell*>;
   }
   static void *newArray_maplEunsignedsPintcOPHG4CylinderCellmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<unsigned int,PHG4CylinderCell*>[nElements] : new map<unsigned int,PHG4CylinderCell*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEunsignedsPintcOPHG4CylinderCellmUgR(void *p) {
      delete ((map<unsigned int,PHG4CylinderCell*>*)p);
   }
   static void deleteArray_maplEunsignedsPintcOPHG4CylinderCellmUgR(void *p) {
      delete [] ((map<unsigned int,PHG4CylinderCell*>*)p);
   }
   static void destruct_maplEunsignedsPintcOPHG4CylinderCellmUgR(void *p) {
      typedef map<unsigned int,PHG4CylinderCell*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<unsigned int,PHG4CylinderCell*>

namespace {
  void TriggerDictionaryInitialization_PHG4CylinderCellContainer_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4CylinderCellContainer_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4CylinderCellContainer;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4CylinderCellContainer_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERCELLCONTAINER_H
#define G4DETECTORS_PHG4CYLINDERCELLCONTAINER_H

#include "PHG4CylinderCellDefs.h"  // for keytype

#include <phool/PHObject.h>

#include <iostream>                // for cout, ostream
#include <map>
#include <set>
#include <utility>                 // for pair, make_pair

class PHG4CylinderCell;

class PHG4CylinderCellContainer: public PHObject
{

  public:
  typedef std::map<PHG4CylinderCellDefs::keytype,PHG4CylinderCell *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  typedef std::set<int>::const_iterator LayerIter;
  typedef std::pair<LayerIter, LayerIter> LayerRange;

  PHG4CylinderCellContainer(){}

  virtual ~PHG4CylinderCellContainer() {}

  void Reset();

  void identify(std::ostream& os = std::cout) const;

  ConstIterator AddCylinderCell(const unsigned int detid, PHG4CylinderCell *newcylinderCell);
  ConstIterator AddCylinderCellSpecifyKey(const PHG4CylinderCellDefs::keytype key, PHG4CylinderCell *newcylinderCell);
  
  //! preferred removal method, key is currently the cell id
  void RemoveCylinderCell(PHG4CylinderCellDefs::keytype key) {
    cellmap.erase(key);
  }

  //! inefficent, use key where possible instead
  void RemoveCylinderCell(PHG4CylinderCell *cell)
  {
    Iterator its = cellmap.begin();
    Iterator last = cellmap.end();
    for (; its != last;)
      {
	if (its->second == cell)
	  {
	    cellmap.erase(its++);
	  }
	else
	  {
	    ++its;
	  }
      }
  }


  Iterator findOrAddCylinderCell(PHG4CylinderCellDefs::keytype key);

  PHG4CylinderCellDefs::keytype genkey(const unsigned int detid);

  //! return all cylinderCells matching a given detid
  ConstRange getCylinderCells(const unsigned int detid) const;

  //! return all hist
  ConstRange getCylinderCells( void ) const;

  PHG4CylinderCell* findCylinderCell(PHG4CylinderCellDefs::keytype key);

  unsigned int size( void ) const
  { return cellmap.size(); }
  unsigned int num_layers(void) const
  { return layers.size(); }
  LayerRange getLayers() const
  { return make_pair(layers.begin(), layers.end());}

  double getTotalEdep() const;

 protected:
  Map cellmap;
  std::set<int> layers; // layers is not reset since layers must not change event by event

  ClassDef(PHG4CylinderCellContainer,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4CylinderCellContainer", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4CylinderCellContainer_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4CylinderCellContainer_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4CylinderCellContainer_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4CylinderCellContainer_Dict() {
  TriggerDictionaryInitialization_PHG4CylinderCellContainer_Dict_Impl();
}
