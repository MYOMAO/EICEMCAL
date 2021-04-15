// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4HitContainer_Dict
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
#include "PHG4HitContainer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4HitContainer(void *p = 0);
   static void *newArray_PHG4HitContainer(Long_t size, void *p);
   static void delete_PHG4HitContainer(void *p);
   static void deleteArray_PHG4HitContainer(void *p);
   static void destruct_PHG4HitContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4HitContainer*)
   {
      ::PHG4HitContainer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4HitContainer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4HitContainer", ::PHG4HitContainer::Class_Version(), "", 20,
                  typeid(::PHG4HitContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4HitContainer::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4HitContainer) );
      instance.SetNew(&new_PHG4HitContainer);
      instance.SetNewArray(&newArray_PHG4HitContainer);
      instance.SetDelete(&delete_PHG4HitContainer);
      instance.SetDeleteArray(&deleteArray_PHG4HitContainer);
      instance.SetDestructor(&destruct_PHG4HitContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4HitContainer*)
   {
      return GenerateInitInstanceLocal((::PHG4HitContainer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4HitContainer*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4HitContainer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4HitContainer::Class_Name()
{
   return "PHG4HitContainer";
}

//______________________________________________________________________________
const char *PHG4HitContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4HitContainer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4HitContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4HitContainer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4HitContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4HitContainer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4HitContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4HitContainer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4HitContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4HitContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4HitContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4HitContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4HitContainer(void *p) {
      return  p ? new(p) ::PHG4HitContainer : new ::PHG4HitContainer;
   }
   static void *newArray_PHG4HitContainer(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4HitContainer[nElements] : new ::PHG4HitContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4HitContainer(void *p) {
      delete ((::PHG4HitContainer*)p);
   }
   static void deleteArray_PHG4HitContainer(void *p) {
      delete [] ((::PHG4HitContainer*)p);
   }
   static void destruct_PHG4HitContainer(void *p) {
      typedef ::PHG4HitContainer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4HitContainer

namespace ROOT {
   static TClass *setlEunsignedsPintgR_Dictionary();
   static void setlEunsignedsPintgR_TClassManip(TClass*);
   static void *new_setlEunsignedsPintgR(void *p = 0);
   static void *newArray_setlEunsignedsPintgR(Long_t size, void *p);
   static void delete_setlEunsignedsPintgR(void *p);
   static void deleteArray_setlEunsignedsPintgR(void *p);
   static void destruct_setlEunsignedsPintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const set<unsigned int>*)
   {
      set<unsigned int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(set<unsigned int>));
      static ::ROOT::TGenericClassInfo 
         instance("set<unsigned int>", -2, "set", 94,
                  typeid(set<unsigned int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &setlEunsignedsPintgR_Dictionary, isa_proxy, 0,
                  sizeof(set<unsigned int>) );
      instance.SetNew(&new_setlEunsignedsPintgR);
      instance.SetNewArray(&newArray_setlEunsignedsPintgR);
      instance.SetDelete(&delete_setlEunsignedsPintgR);
      instance.SetDeleteArray(&deleteArray_setlEunsignedsPintgR);
      instance.SetDestructor(&destruct_setlEunsignedsPintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Insert< set<unsigned int> >()));

      ::ROOT::AddClassAlternate("set<unsigned int>","std::set<unsigned int, std::less<unsigned int>, std::allocator<unsigned int> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const set<unsigned int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *setlEunsignedsPintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const set<unsigned int>*)0x0)->GetClass();
      setlEunsignedsPintgR_TClassManip(theClass);
   return theClass;
   }

   static void setlEunsignedsPintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_setlEunsignedsPintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<unsigned int> : new set<unsigned int>;
   }
   static void *newArray_setlEunsignedsPintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<unsigned int>[nElements] : new set<unsigned int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_setlEunsignedsPintgR(void *p) {
      delete ((set<unsigned int>*)p);
   }
   static void deleteArray_setlEunsignedsPintgR(void *p) {
      delete [] ((set<unsigned int>*)p);
   }
   static void destruct_setlEunsignedsPintgR(void *p) {
      typedef set<unsigned int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class set<unsigned int>

namespace ROOT {
   static TClass *maplEULong64_tcOPHG4HitmUgR_Dictionary();
   static void maplEULong64_tcOPHG4HitmUgR_TClassManip(TClass*);
   static void *new_maplEULong64_tcOPHG4HitmUgR(void *p = 0);
   static void *newArray_maplEULong64_tcOPHG4HitmUgR(Long_t size, void *p);
   static void delete_maplEULong64_tcOPHG4HitmUgR(void *p);
   static void deleteArray_maplEULong64_tcOPHG4HitmUgR(void *p);
   static void destruct_maplEULong64_tcOPHG4HitmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<ULong64_t,PHG4Hit*>*)
   {
      map<ULong64_t,PHG4Hit*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<ULong64_t,PHG4Hit*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<ULong64_t,PHG4Hit*>", -2, "map", 100,
                  typeid(map<ULong64_t,PHG4Hit*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEULong64_tcOPHG4HitmUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<ULong64_t,PHG4Hit*>) );
      instance.SetNew(&new_maplEULong64_tcOPHG4HitmUgR);
      instance.SetNewArray(&newArray_maplEULong64_tcOPHG4HitmUgR);
      instance.SetDelete(&delete_maplEULong64_tcOPHG4HitmUgR);
      instance.SetDeleteArray(&deleteArray_maplEULong64_tcOPHG4HitmUgR);
      instance.SetDestructor(&destruct_maplEULong64_tcOPHG4HitmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<ULong64_t,PHG4Hit*> >()));

      ::ROOT::AddClassAlternate("map<ULong64_t,PHG4Hit*>","std::map<unsigned long long, PHG4Hit*, std::less<unsigned long long>, std::allocator<std::pair<unsigned long long const, PHG4Hit*> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<ULong64_t,PHG4Hit*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEULong64_tcOPHG4HitmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<ULong64_t,PHG4Hit*>*)0x0)->GetClass();
      maplEULong64_tcOPHG4HitmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEULong64_tcOPHG4HitmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEULong64_tcOPHG4HitmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<ULong64_t,PHG4Hit*> : new map<ULong64_t,PHG4Hit*>;
   }
   static void *newArray_maplEULong64_tcOPHG4HitmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<ULong64_t,PHG4Hit*>[nElements] : new map<ULong64_t,PHG4Hit*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEULong64_tcOPHG4HitmUgR(void *p) {
      delete ((map<ULong64_t,PHG4Hit*>*)p);
   }
   static void deleteArray_maplEULong64_tcOPHG4HitmUgR(void *p) {
      delete [] ((map<ULong64_t,PHG4Hit*>*)p);
   }
   static void destruct_maplEULong64_tcOPHG4HitmUgR(void *p) {
      typedef map<ULong64_t,PHG4Hit*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<ULong64_t,PHG4Hit*>

namespace {
  void TriggerDictionaryInitialization_PHG4HitContainer_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4HitContainer_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4HitContainer;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4HitContainer_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4HITCONTAINER_H
#define G4MAIN_PHG4HITCONTAINER_H

#include "PHG4HitDefs.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>

class PHG4Hit;

class PHG4HitContainer: public PHObject
{

  public:
  typedef std::map<PHG4HitDefs::keytype, PHG4Hit *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  typedef std::set<unsigned int>::const_iterator LayerIter;

  PHG4HitContainer(); //< used only by ROOT for DST readback
  PHG4HitContainer(const std::string &nodename);

  virtual ~PHG4HitContainer() {}

  void Reset();

  void identify(std::ostream& os = std::cout) const;

  //! container ID should follow definition of PHG4HitDefs::get_volume_id(DST nodename)
  void SetID(int i) {id = i;}
  int GetID() const {return id;}
  
  ConstIterator AddHit(PHG4Hit *newhit);

  ConstIterator AddHit(const unsigned int detid, PHG4Hit *newhit);
  
  Iterator findOrAddHit(PHG4HitDefs::keytype key);

  PHG4Hit* findHit(PHG4HitDefs::keytype key );

  PHG4HitDefs::keytype genkey(const unsigned int detid);

  //! return all hits matching a given detid
  ConstRange getHits(const unsigned int detid) const;

  //! return all hist
  ConstRange getHits( void ) const;

  unsigned int size( void ) const
  { return hitmap.size(); }
  unsigned int num_layers(void) const
  { return layers.size(); }
  std::pair<LayerIter, LayerIter> getLayers() const
     { return make_pair(layers.begin(), layers.end());} 
  void AddLayer(const unsigned int ilayer) {layers.insert(ilayer);}
  void RemoveZeroEDep();
  PHG4HitDefs::keytype getmaxkey(const unsigned int detid);

 protected:

  int id; //< unique identifier from hash of node name. Defined following PHG4HitDefs::get_volume_id
  Map hitmap;
  std::set<unsigned int> layers; // layers is not reset since layers must not change event by event

  ClassDef(PHG4HitContainer,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4HitContainer", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4HitContainer_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4HitContainer_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4HitContainer_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4HitContainer_Dict() {
  TriggerDictionaryInitialization_PHG4HitContainer_Dict_Impl();
}
