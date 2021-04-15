// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4InEvent_Dict
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
#include "PHG4InEvent.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4InEvent(void *p = 0);
   static void *newArray_PHG4InEvent(Long_t size, void *p);
   static void delete_PHG4InEvent(void *p);
   static void deleteArray_PHG4InEvent(void *p);
   static void destruct_PHG4InEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4InEvent*)
   {
      ::PHG4InEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4InEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4InEvent", ::PHG4InEvent::Class_Version(), "", 17,
                  typeid(::PHG4InEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4InEvent::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4InEvent) );
      instance.SetNew(&new_PHG4InEvent);
      instance.SetNewArray(&newArray_PHG4InEvent);
      instance.SetDelete(&delete_PHG4InEvent);
      instance.SetDeleteArray(&deleteArray_PHG4InEvent);
      instance.SetDestructor(&destruct_PHG4InEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4InEvent*)
   {
      return GenerateInitInstanceLocal((::PHG4InEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4InEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4InEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4InEvent::Class_Name()
{
   return "PHG4InEvent";
}

//______________________________________________________________________________
const char *PHG4InEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4InEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4InEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4InEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4InEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4InEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4InEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4InEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4InEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4InEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4InEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4InEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4InEvent(void *p) {
      return  p ? new(p) ::PHG4InEvent : new ::PHG4InEvent;
   }
   static void *newArray_PHG4InEvent(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4InEvent[nElements] : new ::PHG4InEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4InEvent(void *p) {
      delete ((::PHG4InEvent*)p);
   }
   static void deleteArray_PHG4InEvent(void *p) {
      delete [] ((::PHG4InEvent*)p);
   }
   static void destruct_PHG4InEvent(void *p) {
      typedef ::PHG4InEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4InEvent

namespace ROOT {
   static TClass *multimaplEintcOPHG4ParticlemUgR_Dictionary();
   static void multimaplEintcOPHG4ParticlemUgR_TClassManip(TClass*);
   static void *new_multimaplEintcOPHG4ParticlemUgR(void *p = 0);
   static void *newArray_multimaplEintcOPHG4ParticlemUgR(Long_t size, void *p);
   static void delete_multimaplEintcOPHG4ParticlemUgR(void *p);
   static void deleteArray_multimaplEintcOPHG4ParticlemUgR(void *p);
   static void destruct_multimaplEintcOPHG4ParticlemUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const multimap<int,PHG4Particle*>*)
   {
      multimap<int,PHG4Particle*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(multimap<int,PHG4Particle*>));
      static ::ROOT::TGenericClassInfo 
         instance("multimap<int,PHG4Particle*>", -2, "map", 99,
                  typeid(multimap<int,PHG4Particle*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &multimaplEintcOPHG4ParticlemUgR_Dictionary, isa_proxy, 0,
                  sizeof(multimap<int,PHG4Particle*>) );
      instance.SetNew(&new_multimaplEintcOPHG4ParticlemUgR);
      instance.SetNewArray(&newArray_multimaplEintcOPHG4ParticlemUgR);
      instance.SetDelete(&delete_multimaplEintcOPHG4ParticlemUgR);
      instance.SetDeleteArray(&deleteArray_multimaplEintcOPHG4ParticlemUgR);
      instance.SetDestructor(&destruct_multimaplEintcOPHG4ParticlemUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< multimap<int,PHG4Particle*> >()));

      ::ROOT::AddClassAlternate("multimap<int,PHG4Particle*>","std::multimap<int, PHG4Particle*, std::less<int>, std::allocator<std::pair<int const, PHG4Particle*> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const multimap<int,PHG4Particle*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *multimaplEintcOPHG4ParticlemUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const multimap<int,PHG4Particle*>*)0x0)->GetClass();
      multimaplEintcOPHG4ParticlemUgR_TClassManip(theClass);
   return theClass;
   }

   static void multimaplEintcOPHG4ParticlemUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_multimaplEintcOPHG4ParticlemUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) multimap<int,PHG4Particle*> : new multimap<int,PHG4Particle*>;
   }
   static void *newArray_multimaplEintcOPHG4ParticlemUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) multimap<int,PHG4Particle*>[nElements] : new multimap<int,PHG4Particle*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_multimaplEintcOPHG4ParticlemUgR(void *p) {
      delete ((multimap<int,PHG4Particle*>*)p);
   }
   static void deleteArray_multimaplEintcOPHG4ParticlemUgR(void *p) {
      delete [] ((multimap<int,PHG4Particle*>*)p);
   }
   static void destruct_multimaplEintcOPHG4ParticlemUgR(void *p) {
      typedef multimap<int,PHG4Particle*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class multimap<int,PHG4Particle*>

namespace ROOT {
   static TClass *maplEintcOPHG4VtxPointmUgR_Dictionary();
   static void maplEintcOPHG4VtxPointmUgR_TClassManip(TClass*);
   static void *new_maplEintcOPHG4VtxPointmUgR(void *p = 0);
   static void *newArray_maplEintcOPHG4VtxPointmUgR(Long_t size, void *p);
   static void delete_maplEintcOPHG4VtxPointmUgR(void *p);
   static void deleteArray_maplEintcOPHG4VtxPointmUgR(void *p);
   static void destruct_maplEintcOPHG4VtxPointmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,PHG4VtxPoint*>*)
   {
      map<int,PHG4VtxPoint*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,PHG4VtxPoint*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,PHG4VtxPoint*>", -2, "map", 100,
                  typeid(map<int,PHG4VtxPoint*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOPHG4VtxPointmUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,PHG4VtxPoint*>) );
      instance.SetNew(&new_maplEintcOPHG4VtxPointmUgR);
      instance.SetNewArray(&newArray_maplEintcOPHG4VtxPointmUgR);
      instance.SetDelete(&delete_maplEintcOPHG4VtxPointmUgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOPHG4VtxPointmUgR);
      instance.SetDestructor(&destruct_maplEintcOPHG4VtxPointmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,PHG4VtxPoint*> >()));

      ::ROOT::AddClassAlternate("map<int,PHG4VtxPoint*>","std::map<int, PHG4VtxPoint*, std::less<int>, std::allocator<std::pair<int const, PHG4VtxPoint*> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,PHG4VtxPoint*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOPHG4VtxPointmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,PHG4VtxPoint*>*)0x0)->GetClass();
      maplEintcOPHG4VtxPointmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOPHG4VtxPointmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOPHG4VtxPointmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,PHG4VtxPoint*> : new map<int,PHG4VtxPoint*>;
   }
   static void *newArray_maplEintcOPHG4VtxPointmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,PHG4VtxPoint*>[nElements] : new map<int,PHG4VtxPoint*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOPHG4VtxPointmUgR(void *p) {
      delete ((map<int,PHG4VtxPoint*>*)p);
   }
   static void deleteArray_maplEintcOPHG4VtxPointmUgR(void *p) {
      delete [] ((map<int,PHG4VtxPoint*>*)p);
   }
   static void destruct_maplEintcOPHG4VtxPointmUgR(void *p) {
      typedef map<int,PHG4VtxPoint*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,PHG4VtxPoint*>

namespace ROOT {
   static TClass *maplEPHG4ParticlemUcOintgR_Dictionary();
   static void maplEPHG4ParticlemUcOintgR_TClassManip(TClass*);
   static void *new_maplEPHG4ParticlemUcOintgR(void *p = 0);
   static void *newArray_maplEPHG4ParticlemUcOintgR(Long_t size, void *p);
   static void delete_maplEPHG4ParticlemUcOintgR(void *p);
   static void deleteArray_maplEPHG4ParticlemUcOintgR(void *p);
   static void destruct_maplEPHG4ParticlemUcOintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<PHG4Particle*,int>*)
   {
      map<PHG4Particle*,int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<PHG4Particle*,int>));
      static ::ROOT::TGenericClassInfo 
         instance("map<PHG4Particle*,int>", -2, "map", 100,
                  typeid(map<PHG4Particle*,int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEPHG4ParticlemUcOintgR_Dictionary, isa_proxy, 0,
                  sizeof(map<PHG4Particle*,int>) );
      instance.SetNew(&new_maplEPHG4ParticlemUcOintgR);
      instance.SetNewArray(&newArray_maplEPHG4ParticlemUcOintgR);
      instance.SetDelete(&delete_maplEPHG4ParticlemUcOintgR);
      instance.SetDeleteArray(&deleteArray_maplEPHG4ParticlemUcOintgR);
      instance.SetDestructor(&destruct_maplEPHG4ParticlemUcOintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<PHG4Particle*,int> >()));

      ::ROOT::AddClassAlternate("map<PHG4Particle*,int>","std::map<PHG4Particle*, int, std::less<PHG4Particle*>, std::allocator<std::pair<PHG4Particle* const, int> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<PHG4Particle*,int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEPHG4ParticlemUcOintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<PHG4Particle*,int>*)0x0)->GetClass();
      maplEPHG4ParticlemUcOintgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEPHG4ParticlemUcOintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEPHG4ParticlemUcOintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<PHG4Particle*,int> : new map<PHG4Particle*,int>;
   }
   static void *newArray_maplEPHG4ParticlemUcOintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<PHG4Particle*,int>[nElements] : new map<PHG4Particle*,int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEPHG4ParticlemUcOintgR(void *p) {
      delete ((map<PHG4Particle*,int>*)p);
   }
   static void deleteArray_maplEPHG4ParticlemUcOintgR(void *p) {
      delete [] ((map<PHG4Particle*,int>*)p);
   }
   static void destruct_maplEPHG4ParticlemUcOintgR(void *p) {
      typedef map<PHG4Particle*,int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<PHG4Particle*,int>

namespace {
  void TriggerDictionaryInitialization_PHG4InEvent_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4InEvent_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4InEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4InEvent_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4INEVENT_H
#define G4MAIN_PHG4INEVENT_H

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <utility>

class PHG4Particle;
class PHG4VtxPoint;

class PHG4InEvent: public PHObject
{
 public:
  PHG4InEvent() {}
  virtual ~PHG4InEvent();

  virtual void identify(std::ostream& os = std::cout) const;
  void Reset();
  // this one is for HepMC records where we want to keep the HepMC vertex numbering
  int AddVtxHepMC(const int id, const double x, const double y, const double z, const double t);
  int AddVtx(const double x, const double y, const double z, const double t);
  int AddVtx(const int id,const PHG4VtxPoint &);
  int AddParticle(const int vtxid, PHG4Particle *particle);
  void AddEmbeddedParticle(PHG4Particle *particle, int flag) {embedded_particlelist.insert(std::make_pair(particle,flag));}

  //  PHG4VtxPoint *GetVtx() {return vtxlist.begin()->second;}
  std::pair< std::map<int, PHG4VtxPoint *>::const_iterator, std::map<int, PHG4VtxPoint *>::const_iterator > GetVertices() const;
  std::pair< std::multimap<int,PHG4Particle *>::const_iterator, std::multimap<int,PHG4Particle *>::const_iterator > GetParticles(const int vtxid) const;
  std::pair< std::multimap<int,PHG4Particle *>::const_iterator, std::multimap<int,PHG4Particle *>::const_iterator > GetParticles() const;
  std::pair< std::multimap<int,PHG4Particle *>::iterator, std::multimap<int,PHG4Particle *>::iterator > GetParticles_Modify();
  std::pair< std::map<PHG4Particle *,int>::const_iterator, std::map<PHG4Particle *,int>::const_iterator> GetEmbeddedParticles() const
    {return std::make_pair(embedded_particlelist.begin(), embedded_particlelist.end());}
  int isEmbeded(PHG4Particle *) const;
  int GetNEmbedded() const {return embedded_particlelist.size();}
  int GetNVtx() const {return vtxlist.size();}
  void DeleteParticle(std::multimap<int,PHG4Particle *>::iterator &iter);

 protected:

  int AddVtxCommon(PHG4VtxPoint *newvtx);
  std::map<int,PHG4VtxPoint *> vtxlist;
  std::multimap<int,PHG4Particle *> particlelist;
  std::map<PHG4Particle *,int> embedded_particlelist;

  ClassDef(PHG4InEvent,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4InEvent", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4InEvent_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4InEvent_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4InEvent_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4InEvent_Dict() {
  TriggerDictionaryInitialization_PHG4InEvent_Dict_Impl();
}
