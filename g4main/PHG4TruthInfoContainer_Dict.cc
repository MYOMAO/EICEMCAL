// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4TruthInfoContainer_Dict
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
#include "PHG4TruthInfoContainer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4TruthInfoContainer(void *p = 0);
   static void *newArray_PHG4TruthInfoContainer(Long_t size, void *p);
   static void delete_PHG4TruthInfoContainer(void *p);
   static void deleteArray_PHG4TruthInfoContainer(void *p);
   static void destruct_PHG4TruthInfoContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4TruthInfoContainer*)
   {
      ::PHG4TruthInfoContainer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4TruthInfoContainer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4TruthInfoContainer", ::PHG4TruthInfoContainer::Class_Version(), "", 19,
                  typeid(::PHG4TruthInfoContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4TruthInfoContainer::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4TruthInfoContainer) );
      instance.SetNew(&new_PHG4TruthInfoContainer);
      instance.SetNewArray(&newArray_PHG4TruthInfoContainer);
      instance.SetDelete(&delete_PHG4TruthInfoContainer);
      instance.SetDeleteArray(&deleteArray_PHG4TruthInfoContainer);
      instance.SetDestructor(&destruct_PHG4TruthInfoContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4TruthInfoContainer*)
   {
      return GenerateInitInstanceLocal((::PHG4TruthInfoContainer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4TruthInfoContainer*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4TruthInfoContainer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4TruthInfoContainer::Class_Name()
{
   return "PHG4TruthInfoContainer";
}

//______________________________________________________________________________
const char *PHG4TruthInfoContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4TruthInfoContainer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4TruthInfoContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4TruthInfoContainer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4TruthInfoContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4TruthInfoContainer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4TruthInfoContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4TruthInfoContainer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4TruthInfoContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4TruthInfoContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4TruthInfoContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4TruthInfoContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4TruthInfoContainer(void *p) {
      return  p ? new(p) ::PHG4TruthInfoContainer : new ::PHG4TruthInfoContainer;
   }
   static void *newArray_PHG4TruthInfoContainer(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4TruthInfoContainer[nElements] : new ::PHG4TruthInfoContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4TruthInfoContainer(void *p) {
      delete ((::PHG4TruthInfoContainer*)p);
   }
   static void deleteArray_PHG4TruthInfoContainer(void *p) {
      delete [] ((::PHG4TruthInfoContainer*)p);
   }
   static void destruct_PHG4TruthInfoContainer(void *p) {
      typedef ::PHG4TruthInfoContainer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4TruthInfoContainer

namespace ROOT {
   static TClass *maplEintcOintgR_Dictionary();
   static void maplEintcOintgR_TClassManip(TClass*);
   static void *new_maplEintcOintgR(void *p = 0);
   static void *newArray_maplEintcOintgR(Long_t size, void *p);
   static void delete_maplEintcOintgR(void *p);
   static void deleteArray_maplEintcOintgR(void *p);
   static void destruct_maplEintcOintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,int>*)
   {
      map<int,int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,int>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,int>", -2, "map", 100,
                  typeid(map<int,int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOintgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,int>) );
      instance.SetNew(&new_maplEintcOintgR);
      instance.SetNewArray(&newArray_maplEintcOintgR);
      instance.SetDelete(&delete_maplEintcOintgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOintgR);
      instance.SetDestructor(&destruct_maplEintcOintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,int> >()));

      ::ROOT::AddClassAlternate("map<int,int>","std::map<int, int, std::less<int>, std::allocator<std::pair<int const, int> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,int>*)0x0)->GetClass();
      maplEintcOintgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,int> : new map<int,int>;
   }
   static void *newArray_maplEintcOintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,int>[nElements] : new map<int,int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOintgR(void *p) {
      delete ((map<int,int>*)p);
   }
   static void deleteArray_maplEintcOintgR(void *p) {
      delete [] ((map<int,int>*)p);
   }
   static void destruct_maplEintcOintgR(void *p) {
      typedef map<int,int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,int>

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
   static TClass *maplEintcOPHG4ShowermUgR_Dictionary();
   static void maplEintcOPHG4ShowermUgR_TClassManip(TClass*);
   static void *new_maplEintcOPHG4ShowermUgR(void *p = 0);
   static void *newArray_maplEintcOPHG4ShowermUgR(Long_t size, void *p);
   static void delete_maplEintcOPHG4ShowermUgR(void *p);
   static void deleteArray_maplEintcOPHG4ShowermUgR(void *p);
   static void destruct_maplEintcOPHG4ShowermUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,PHG4Shower*>*)
   {
      map<int,PHG4Shower*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,PHG4Shower*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,PHG4Shower*>", -2, "map", 100,
                  typeid(map<int,PHG4Shower*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOPHG4ShowermUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,PHG4Shower*>) );
      instance.SetNew(&new_maplEintcOPHG4ShowermUgR);
      instance.SetNewArray(&newArray_maplEintcOPHG4ShowermUgR);
      instance.SetDelete(&delete_maplEintcOPHG4ShowermUgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOPHG4ShowermUgR);
      instance.SetDestructor(&destruct_maplEintcOPHG4ShowermUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,PHG4Shower*> >()));

      ::ROOT::AddClassAlternate("map<int,PHG4Shower*>","std::map<int, PHG4Shower*, std::less<int>, std::allocator<std::pair<int const, PHG4Shower*> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,PHG4Shower*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOPHG4ShowermUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,PHG4Shower*>*)0x0)->GetClass();
      maplEintcOPHG4ShowermUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOPHG4ShowermUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOPHG4ShowermUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,PHG4Shower*> : new map<int,PHG4Shower*>;
   }
   static void *newArray_maplEintcOPHG4ShowermUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,PHG4Shower*>[nElements] : new map<int,PHG4Shower*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOPHG4ShowermUgR(void *p) {
      delete ((map<int,PHG4Shower*>*)p);
   }
   static void deleteArray_maplEintcOPHG4ShowermUgR(void *p) {
      delete [] ((map<int,PHG4Shower*>*)p);
   }
   static void destruct_maplEintcOPHG4ShowermUgR(void *p) {
      typedef map<int,PHG4Shower*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,PHG4Shower*>

namespace ROOT {
   static TClass *maplEintcOPHG4ParticlemUgR_Dictionary();
   static void maplEintcOPHG4ParticlemUgR_TClassManip(TClass*);
   static void *new_maplEintcOPHG4ParticlemUgR(void *p = 0);
   static void *newArray_maplEintcOPHG4ParticlemUgR(Long_t size, void *p);
   static void delete_maplEintcOPHG4ParticlemUgR(void *p);
   static void deleteArray_maplEintcOPHG4ParticlemUgR(void *p);
   static void destruct_maplEintcOPHG4ParticlemUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,PHG4Particle*>*)
   {
      map<int,PHG4Particle*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,PHG4Particle*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,PHG4Particle*>", -2, "map", 100,
                  typeid(map<int,PHG4Particle*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOPHG4ParticlemUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,PHG4Particle*>) );
      instance.SetNew(&new_maplEintcOPHG4ParticlemUgR);
      instance.SetNewArray(&newArray_maplEintcOPHG4ParticlemUgR);
      instance.SetDelete(&delete_maplEintcOPHG4ParticlemUgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOPHG4ParticlemUgR);
      instance.SetDestructor(&destruct_maplEintcOPHG4ParticlemUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,PHG4Particle*> >()));

      ::ROOT::AddClassAlternate("map<int,PHG4Particle*>","std::map<int, PHG4Particle*, std::less<int>, std::allocator<std::pair<int const, PHG4Particle*> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,PHG4Particle*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOPHG4ParticlemUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,PHG4Particle*>*)0x0)->GetClass();
      maplEintcOPHG4ParticlemUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOPHG4ParticlemUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOPHG4ParticlemUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,PHG4Particle*> : new map<int,PHG4Particle*>;
   }
   static void *newArray_maplEintcOPHG4ParticlemUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,PHG4Particle*>[nElements] : new map<int,PHG4Particle*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOPHG4ParticlemUgR(void *p) {
      delete ((map<int,PHG4Particle*>*)p);
   }
   static void deleteArray_maplEintcOPHG4ParticlemUgR(void *p) {
      delete [] ((map<int,PHG4Particle*>*)p);
   }
   static void destruct_maplEintcOPHG4ParticlemUgR(void *p) {
      typedef map<int,PHG4Particle*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,PHG4Particle*>

namespace {
  void TriggerDictionaryInitialization_PHG4TruthInfoContainer_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4TruthInfoContainer_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4TruthInfoContainer;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4TruthInfoContainer_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4TRUTHINFOCONTAINER_H
#define G4MAIN_PHG4TRUTHINFOCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>
#include <iterator>  // for distance
#include <map>
#include <utility>

class PHG4Shower;
class PHG4Particle;
class PHG4VtxPoint;

class PHG4TruthInfoContainer : public PHObject
{
 public:
  typedef std::map<int, PHG4Particle*> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  typedef std::map<int, PHG4VtxPoint*> VtxMap;
  typedef VtxMap::iterator VtxIterator;
  typedef VtxMap::const_iterator ConstVtxIterator;
  typedef std::pair<VtxIterator, VtxIterator> VtxRange;
  typedef std::pair<ConstVtxIterator, ConstVtxIterator> ConstVtxRange;

  typedef std::map<int, PHG4Shower*> ShowerMap;
  typedef ShowerMap::iterator ShowerIterator;
  typedef ShowerMap::const_iterator ConstShowerIterator;
  typedef std::pair<ShowerIterator, ShowerIterator> ShowerRange;
  typedef std::pair<ConstShowerIterator, ConstShowerIterator> ConstShowerRange;

  PHG4TruthInfoContainer();
  virtual ~PHG4TruthInfoContainer();

  void Reset();
  void identify(std::ostream& os = std::cout) const;

  // --- particle storage ------------------------------------------------------

  //! Add a particle that the user has created
  ConstIterator AddParticle(const int particleid, PHG4Particle* newparticle);
  void delete_particle(Iterator piter);

  PHG4Particle* GetParticle(const int trackid);
  PHG4Particle* GetPrimaryParticle(const int trackid);

  bool is_primary(const PHG4Particle* p) const;

  //! Get a range of iterators covering the entire container
  Range GetParticleRange() { return Range(particlemap.begin(), particlemap.end()); }
  ConstRange GetParticleRange() const { return ConstRange(particlemap.begin(), particlemap.end()); }

  Range GetPrimaryParticleRange() { return Range(particlemap.upper_bound(0), particlemap.end()); }
  ConstRange GetPrimaryParticleRange() const { return ConstRange(particlemap.upper_bound(0), particlemap.end()); }

  Range GetSecondaryParticleRange() { return Range(particlemap.begin(), particlemap.upper_bound(0)); }
  ConstRange GetSecondaryParticleRange() const { return ConstRange(particlemap.begin(), particlemap.upper_bound(0)); }

  //! track -> particle map size
  unsigned int size(void) const { return particlemap.size(); }
  int GetNumPrimaryVertexParticles()
  {
    return std::distance(particlemap.upper_bound(0), particlemap.end());
  }

  //! Get the Particle Map storage
  const Map& GetMap() const { return particlemap; }

  int maxtrkindex() const;
  int mintrkindex() const;

  //! Retrieve the embedding ID for the HepMC subevent or track to be analyzed.
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  std::pair<std::map<int, int>::const_iterator,
            std::map<int, int>::const_iterator>
  GetEmbeddedTrkIds() const
  {
    return std::make_pair(particle_embed_flags.begin(), particle_embed_flags.end());
  }

  //! Set the embedding ID for the HepMC subevent or track to be analyzed.
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  void AddEmbededTrkId(const int id, const int flag)
  {
    particle_embed_flags.insert(std::make_pair(id, flag));
  }

  //! Retrieve the embedding ID for the HepMC subevent or track to be analyzed.
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int isEmbeded(const int trackid) const;

  // --- vertex storage --------------------------------------------------------

  //! Add a vertex and return an iterator to the user
  ConstVtxIterator AddVertex(const int vtxid, PHG4VtxPoint* vertex);
  void delete_vtx(VtxIterator viter);

  PHG4VtxPoint* GetVtx(const int vtxid);
  PHG4VtxPoint* GetPrimaryVtx(const int vtxid);

  bool is_primary_vtx(const PHG4VtxPoint* v) const;

  //! Get a range of iterators covering the entire vertex container
  VtxRange GetVtxRange() { return VtxRange(vtxmap.begin(), vtxmap.end()); }
  ConstVtxRange GetVtxRange() const { return ConstVtxRange(vtxmap.begin(), vtxmap.end()); }

  VtxRange GetPrimaryVtxRange() { return VtxRange(vtxmap.upper_bound(0), vtxmap.end()); }
  ConstVtxRange GetPrimaryVtxRange() const { return ConstVtxRange(vtxmap.upper_bound(0), vtxmap.end()); }

  VtxRange GetSecondaryVtxRange() { return VtxRange(vtxmap.begin(), vtxmap.upper_bound(0)); }
  ConstVtxRange GetSecondaryVtxRange() const { return ConstVtxRange(vtxmap.begin(), vtxmap.upper_bound(0)); }

  //! Get the number of vertices stored
  unsigned int GetNumVertices() const { return vtxmap.size(); }

  //! Get the Vertex Map storage
  const VtxMap& GetVtxMap() const { return vtxmap; }

  int maxvtxindex() const;
  int minvtxindex() const;

  //! Return ID of the truth primary vertex with highest embedding ID.
  //! For vertex with identical embedding ID, return first one simulated in Geant4.
  int GetPrimaryVertexIndex() const;

  //! Retrieve the embedding ID for the HepMC subevent or track to be analyzed.
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  std::pair<std::map<int, int>::const_iterator,
            std::map<int, int>::const_iterator>
  GetEmbeddedVtxIds() const
  {
    return std::make_pair(vertex_embed_flags.begin(), vertex_embed_flags.end());
  }

  //! Set the embedding ID for the HepMC subevent or track to be analyzed.
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  void AddEmbededVtxId(const int id, const int flag)
  {
    vertex_embed_flags.insert(std::make_pair(id, flag));
  }

  //! Retrieve the embedding ID for the HepMC subevent or track to be analyzed.
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int isEmbededVtx(const int vtxid) const;

  // --- shower storage ------------------------------------------------------

  //! Add a shower that the user has created
  ConstShowerIterator AddShower(const int showerid, PHG4Shower* newshower);
  void delete_shower(ShowerIterator piter);

  PHG4Shower* GetShower(const int showerid);
  PHG4Shower* GetPrimaryShower(const int showerid);

  //! Get a range of iterators covering the entire container
  ShowerRange GetShowerRange() { return ShowerRange(showermap.begin(), showermap.end()); }
  ConstShowerRange GetShowerRange() const { return ConstShowerRange(showermap.begin(), showermap.end()); }

  ShowerRange GetPrimaryShowerRange() { return ShowerRange(showermap.upper_bound(0), showermap.end()); }
  ConstShowerRange GetPrimaryShowerRange() const { return ConstShowerRange(showermap.upper_bound(0), showermap.end()); }

  ShowerRange GetSecondaryShowerRange() { return ShowerRange(showermap.begin(), showermap.upper_bound(0)); }
  ConstShowerRange GetSecondaryShowerRange() const { return ConstShowerRange(showermap.begin(), showermap.upper_bound(0)); }

  //! shower size
  unsigned int shower_size(void) const { return showermap.size(); }

  //! Get the Shower Map storage
  const ShowerMap& GetShowerMap() const { return showermap; }

  int maxshowerindex() const;
  int minshowerindex() const;

 private:
  /// particle storage map format description:
  /// primary particles are appended in the positive direction
  /// secondary particles are appended in the negative direction
  /// +N   primary particle id => particle*
  /// +N-1
  /// ...
  /// +1   primary particle id => particle*
  /// 0    no entry
  /// -1   secondary particle id => particle*
  /// ...
  /// -M+1
  /// -M   secondary particle id => particle*
  Map particlemap;

  /// vertex storage map format description:
  /// primary vertexes are appended in the positive direction
  /// secondary vertexes are appended in the negative direction
  /// +N   primary vertex id => vertex*
  /// +N-1
  /// ...
  /// +1   primary vertex id => vertex*
  /// 0    no entry
  /// -1   secondary vertex id => vertex*
  /// ...
  /// -M+1
  /// -M   secondary vertex id => vertex*
  VtxMap vtxmap;

  /// shower map
  /// showers encapsulate the secondaries and hits from a primary particle
  ShowerMap showermap;

  // embed flag storage, will typically be set for only a few entries or none at all
  std::map<int, int> particle_embed_flags;  //< trackid => embed flag
  std::map<int, int> vertex_embed_flags;    //< vtxid => embed flag

  ClassDef(PHG4TruthInfoContainer, 1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4TruthInfoContainer", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4TruthInfoContainer_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4TruthInfoContainer_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4TruthInfoContainer_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4TruthInfoContainer_Dict() {
  TriggerDictionaryInitialization_PHG4TruthInfoContainer_Dict_Impl();
}
