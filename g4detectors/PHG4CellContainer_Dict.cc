// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4CellContainer_Dict
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
#include "PHG4CellContainer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4CellContainer(void *p = 0);
   static void *newArray_PHG4CellContainer(Long_t size, void *p);
   static void delete_PHG4CellContainer(void *p);
   static void deleteArray_PHG4CellContainer(void *p);
   static void destruct_PHG4CellContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CellContainer*)
   {
      ::PHG4CellContainer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CellContainer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CellContainer", ::PHG4CellContainer::Class_Version(), "", 21,
                  typeid(::PHG4CellContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CellContainer::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CellContainer) );
      instance.SetNew(&new_PHG4CellContainer);
      instance.SetNewArray(&newArray_PHG4CellContainer);
      instance.SetDelete(&delete_PHG4CellContainer);
      instance.SetDeleteArray(&deleteArray_PHG4CellContainer);
      instance.SetDestructor(&destruct_PHG4CellContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CellContainer*)
   {
      return GenerateInitInstanceLocal((::PHG4CellContainer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CellContainer*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4CellContainer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CellContainer::Class_Name()
{
   return "PHG4CellContainer";
}

//______________________________________________________________________________
const char *PHG4CellContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CellContainer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CellContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CellContainer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CellContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CellContainer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CellContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CellContainer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4CellContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CellContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CellContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CellContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4CellContainer(void *p) {
      return  p ? new(p) ::PHG4CellContainer : new ::PHG4CellContainer;
   }
   static void *newArray_PHG4CellContainer(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4CellContainer[nElements] : new ::PHG4CellContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4CellContainer(void *p) {
      delete ((::PHG4CellContainer*)p);
   }
   static void deleteArray_PHG4CellContainer(void *p) {
      delete [] ((::PHG4CellContainer*)p);
   }
   static void destruct_PHG4CellContainer(void *p) {
      typedef ::PHG4CellContainer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CellContainer

namespace ROOT {
   static TClass *maplEunsignedsPlongcOPHG4CellmUgR_Dictionary();
   static void maplEunsignedsPlongcOPHG4CellmUgR_TClassManip(TClass*);
   static void *new_maplEunsignedsPlongcOPHG4CellmUgR(void *p = 0);
   static void *newArray_maplEunsignedsPlongcOPHG4CellmUgR(Long_t size, void *p);
   static void delete_maplEunsignedsPlongcOPHG4CellmUgR(void *p);
   static void deleteArray_maplEunsignedsPlongcOPHG4CellmUgR(void *p);
   static void destruct_maplEunsignedsPlongcOPHG4CellmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<unsigned long,PHG4Cell*>*)
   {
      map<unsigned long,PHG4Cell*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<unsigned long,PHG4Cell*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<unsigned long,PHG4Cell*>", -2, "map", 100,
                  typeid(map<unsigned long,PHG4Cell*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEunsignedsPlongcOPHG4CellmUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<unsigned long,PHG4Cell*>) );
      instance.SetNew(&new_maplEunsignedsPlongcOPHG4CellmUgR);
      instance.SetNewArray(&newArray_maplEunsignedsPlongcOPHG4CellmUgR);
      instance.SetDelete(&delete_maplEunsignedsPlongcOPHG4CellmUgR);
      instance.SetDeleteArray(&deleteArray_maplEunsignedsPlongcOPHG4CellmUgR);
      instance.SetDestructor(&destruct_maplEunsignedsPlongcOPHG4CellmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<unsigned long,PHG4Cell*> >()));

      ::ROOT::AddClassAlternate("map<unsigned long,PHG4Cell*>","std::map<unsigned long, PHG4Cell*, std::less<unsigned long>, std::allocator<std::pair<unsigned long const, PHG4Cell*> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<unsigned long,PHG4Cell*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEunsignedsPlongcOPHG4CellmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<unsigned long,PHG4Cell*>*)0x0)->GetClass();
      maplEunsignedsPlongcOPHG4CellmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEunsignedsPlongcOPHG4CellmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEunsignedsPlongcOPHG4CellmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<unsigned long,PHG4Cell*> : new map<unsigned long,PHG4Cell*>;
   }
   static void *newArray_maplEunsignedsPlongcOPHG4CellmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<unsigned long,PHG4Cell*>[nElements] : new map<unsigned long,PHG4Cell*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEunsignedsPlongcOPHG4CellmUgR(void *p) {
      delete ((map<unsigned long,PHG4Cell*>*)p);
   }
   static void deleteArray_maplEunsignedsPlongcOPHG4CellmUgR(void *p) {
      delete [] ((map<unsigned long,PHG4Cell*>*)p);
   }
   static void destruct_maplEunsignedsPlongcOPHG4CellmUgR(void *p) {
      typedef map<unsigned long,PHG4Cell*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<unsigned long,PHG4Cell*>

namespace {
  void TriggerDictionaryInitialization_PHG4CellContainer_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4CellContainer_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4CellContainer;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4CellContainer_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CELLCONTAINER_H
#define G4DETECTORS_PHG4CELLCONTAINER_H

#include "PHG4CellDefs.h"    // for keytype

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream
#include <map>
#include <utility>           // for pair

class PHG4Cell;

class PHG4CellContainer: public PHObject
{

  public:
  typedef std::map<PHG4CellDefs::keytype,PHG4Cell *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  PHG4CellContainer();

  virtual ~PHG4CellContainer() {}

  void Reset();

  void identify(std::ostream& os = std::cout) const;

  ConstIterator AddCell(PHG4Cell *newCell);
  ConstIterator AddCellSpecifyKey(const PHG4CellDefs::keytype key, PHG4Cell *newCell);
  
  //! preferred removal method, key is currently the cell id
  void RemoveCell(PHG4CellDefs::keytype key) {
    cellmap.erase(key);
  }

  //! inefficent, use key where possible instead
  void RemoveCell(PHG4Cell *cell)
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


  Iterator findOrAddCell(PHG4CellDefs::keytype key);

  //! return all Cells matching a given detid
  ConstRange getCells(const unsigned short int detid) const;

  //! return all hist
  ConstRange getCells( void ) const;

  PHG4Cell* findCell(PHG4CellDefs::keytype key);

  unsigned int size( void ) const
  { return cellmap.size(); }

  double getTotalEdep() const;

 protected:
  Map cellmap;
  ClassDef(PHG4CellContainer,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4CellContainer", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4CellContainer_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4CellContainer_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4CellContainer_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4CellContainer_Dict() {
  TriggerDictionaryInitialization_PHG4CellContainer_Dict_Impl();
}
