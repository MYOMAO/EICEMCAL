// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4ScintillatorSlatContainer_Dict
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
#include "PHG4ScintillatorSlatContainer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4ScintillatorSlatContainer(void *p = 0);
   static void *newArray_PHG4ScintillatorSlatContainer(Long_t size, void *p);
   static void delete_PHG4ScintillatorSlatContainer(void *p);
   static void deleteArray_PHG4ScintillatorSlatContainer(void *p);
   static void destruct_PHG4ScintillatorSlatContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4ScintillatorSlatContainer*)
   {
      ::PHG4ScintillatorSlatContainer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4ScintillatorSlatContainer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4ScintillatorSlatContainer", ::PHG4ScintillatorSlatContainer::Class_Version(), "", 22,
                  typeid(::PHG4ScintillatorSlatContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4ScintillatorSlatContainer::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4ScintillatorSlatContainer) );
      instance.SetNew(&new_PHG4ScintillatorSlatContainer);
      instance.SetNewArray(&newArray_PHG4ScintillatorSlatContainer);
      instance.SetDelete(&delete_PHG4ScintillatorSlatContainer);
      instance.SetDeleteArray(&deleteArray_PHG4ScintillatorSlatContainer);
      instance.SetDestructor(&destruct_PHG4ScintillatorSlatContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4ScintillatorSlatContainer*)
   {
      return GenerateInitInstanceLocal((::PHG4ScintillatorSlatContainer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4ScintillatorSlatContainer*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4ScintillatorSlatContainer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4ScintillatorSlatContainer::Class_Name()
{
   return "PHG4ScintillatorSlatContainer";
}

//______________________________________________________________________________
const char *PHG4ScintillatorSlatContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4ScintillatorSlatContainer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4ScintillatorSlatContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4ScintillatorSlatContainer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4ScintillatorSlatContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4ScintillatorSlatContainer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4ScintillatorSlatContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4ScintillatorSlatContainer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4ScintillatorSlatContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4ScintillatorSlatContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4ScintillatorSlatContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4ScintillatorSlatContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4ScintillatorSlatContainer(void *p) {
      return  p ? new(p) ::PHG4ScintillatorSlatContainer : new ::PHG4ScintillatorSlatContainer;
   }
   static void *newArray_PHG4ScintillatorSlatContainer(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4ScintillatorSlatContainer[nElements] : new ::PHG4ScintillatorSlatContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4ScintillatorSlatContainer(void *p) {
      delete ((::PHG4ScintillatorSlatContainer*)p);
   }
   static void deleteArray_PHG4ScintillatorSlatContainer(void *p) {
      delete [] ((::PHG4ScintillatorSlatContainer*)p);
   }
   static void destruct_PHG4ScintillatorSlatContainer(void *p) {
      typedef ::PHG4ScintillatorSlatContainer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4ScintillatorSlatContainer

namespace ROOT {
   static TClass *maplEunsignedsPintcOPHG4ScintillatorSlatmUgR_Dictionary();
   static void maplEunsignedsPintcOPHG4ScintillatorSlatmUgR_TClassManip(TClass*);
   static void *new_maplEunsignedsPintcOPHG4ScintillatorSlatmUgR(void *p = 0);
   static void *newArray_maplEunsignedsPintcOPHG4ScintillatorSlatmUgR(Long_t size, void *p);
   static void delete_maplEunsignedsPintcOPHG4ScintillatorSlatmUgR(void *p);
   static void deleteArray_maplEunsignedsPintcOPHG4ScintillatorSlatmUgR(void *p);
   static void destruct_maplEunsignedsPintcOPHG4ScintillatorSlatmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<unsigned int,PHG4ScintillatorSlat*>*)
   {
      map<unsigned int,PHG4ScintillatorSlat*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<unsigned int,PHG4ScintillatorSlat*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<unsigned int,PHG4ScintillatorSlat*>", -2, "map", 100,
                  typeid(map<unsigned int,PHG4ScintillatorSlat*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEunsignedsPintcOPHG4ScintillatorSlatmUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<unsigned int,PHG4ScintillatorSlat*>) );
      instance.SetNew(&new_maplEunsignedsPintcOPHG4ScintillatorSlatmUgR);
      instance.SetNewArray(&newArray_maplEunsignedsPintcOPHG4ScintillatorSlatmUgR);
      instance.SetDelete(&delete_maplEunsignedsPintcOPHG4ScintillatorSlatmUgR);
      instance.SetDeleteArray(&deleteArray_maplEunsignedsPintcOPHG4ScintillatorSlatmUgR);
      instance.SetDestructor(&destruct_maplEunsignedsPintcOPHG4ScintillatorSlatmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<unsigned int,PHG4ScintillatorSlat*> >()));

      ::ROOT::AddClassAlternate("map<unsigned int,PHG4ScintillatorSlat*>","std::map<unsigned int, PHG4ScintillatorSlat*, std::less<unsigned int>, std::allocator<std::pair<unsigned int const, PHG4ScintillatorSlat*> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<unsigned int,PHG4ScintillatorSlat*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEunsignedsPintcOPHG4ScintillatorSlatmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<unsigned int,PHG4ScintillatorSlat*>*)0x0)->GetClass();
      maplEunsignedsPintcOPHG4ScintillatorSlatmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEunsignedsPintcOPHG4ScintillatorSlatmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEunsignedsPintcOPHG4ScintillatorSlatmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<unsigned int,PHG4ScintillatorSlat*> : new map<unsigned int,PHG4ScintillatorSlat*>;
   }
   static void *newArray_maplEunsignedsPintcOPHG4ScintillatorSlatmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<unsigned int,PHG4ScintillatorSlat*>[nElements] : new map<unsigned int,PHG4ScintillatorSlat*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEunsignedsPintcOPHG4ScintillatorSlatmUgR(void *p) {
      delete ((map<unsigned int,PHG4ScintillatorSlat*>*)p);
   }
   static void deleteArray_maplEunsignedsPintcOPHG4ScintillatorSlatmUgR(void *p) {
      delete [] ((map<unsigned int,PHG4ScintillatorSlat*>*)p);
   }
   static void destruct_maplEunsignedsPintcOPHG4ScintillatorSlatmUgR(void *p) {
      typedef map<unsigned int,PHG4ScintillatorSlat*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<unsigned int,PHG4ScintillatorSlat*>

namespace {
  void TriggerDictionaryInitialization_PHG4ScintillatorSlatContainer_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4ScintillatorSlatContainer_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4ScintillatorSlatContainer;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4ScintillatorSlatContainer_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4SCINTILLATORSLATCONTAINER_H
#define G4DETECTORS_PHG4SCINTILLATORSLATCONTAINER_H

#include "PHG4ScintillatorSlatDefs.h"  // for keytype

#include <phool/PHObject.h>

#include <map>
#include <set>
#include <iostream>                    // for cout, ostream
#include <utility>                     // for pair

class PHG4ScintillatorSlat;

class PHG4ScintillatorSlatContainer: public PHObject
{

  public:
  typedef std::map<PHG4ScintillatorSlatDefs::keytype,PHG4ScintillatorSlat *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  typedef std::set<int>::const_iterator LayerIter;
  typedef std::pair<LayerIter, LayerIter> LayerRange;

  PHG4ScintillatorSlatContainer();

  virtual ~PHG4ScintillatorSlatContainer() {}

  void Reset();

  void identify(std::ostream& os = std::cout) const;

  ConstIterator AddScintillatorSlat(const PHG4ScintillatorSlatDefs::keytype key, PHG4ScintillatorSlat *newscintillatorSlat);
  
  //! preferred removal method, key is currently the slat id
  void RemoveScintillatorSlat(PHG4ScintillatorSlatDefs::keytype key) {
    slatmap.erase(key);
  }

  //! inefficent, use key where possible instead
  void RemoveScintillatorSlat(PHG4ScintillatorSlat *slat)
  {
    Iterator its = slatmap.begin();
    Iterator last = slatmap.end();
    for (; its != last;)
      {
	if (its->second == slat)
	  {
	    slatmap.erase(its++);
	  }
	else
	  {
	    ++its;
	  }
      }
  }


  //! return all scintillatorSlats matching a given detid
  ConstRange getScintillatorSlats(const short icolumn) const;

  //! return all hist
  ConstRange getScintillatorSlats( void ) const;

  PHG4ScintillatorSlat* findScintillatorSlat(PHG4ScintillatorSlatDefs::keytype key);

  unsigned int size( void ) const
  { return slatmap.size(); }

  double getTotalEdep() const;

 protected:
  Map slatmap;

  ClassDef(PHG4ScintillatorSlatContainer,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4ScintillatorSlatContainer", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4ScintillatorSlatContainer_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4ScintillatorSlatContainer_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4ScintillatorSlatContainer_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4ScintillatorSlatContainer_Dict() {
  TriggerDictionaryInitialization_PHG4ScintillatorSlatContainer_Dict_Impl();
}
