// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4Shower_Dict
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
#include "PHG4Shower.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_PHG4Shower(void *p);
   static void deleteArray_PHG4Shower(void *p);
   static void destruct_PHG4Shower(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4Shower*)
   {
      ::PHG4Shower *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4Shower >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4Shower", ::PHG4Shower::Class_Version(), "", 17,
                  typeid(::PHG4Shower), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4Shower::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4Shower) );
      instance.SetDelete(&delete_PHG4Shower);
      instance.SetDeleteArray(&deleteArray_PHG4Shower);
      instance.SetDestructor(&destruct_PHG4Shower);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4Shower*)
   {
      return GenerateInitInstanceLocal((::PHG4Shower*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4Shower*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4Shower::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4Shower::Class_Name()
{
   return "PHG4Shower";
}

//______________________________________________________________________________
const char *PHG4Shower::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Shower*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4Shower::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Shower*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4Shower::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Shower*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4Shower::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Shower*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4Shower::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4Shower.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4Shower::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4Shower::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PHG4Shower(void *p) {
      delete ((::PHG4Shower*)p);
   }
   static void deleteArray_PHG4Shower(void *p) {
      delete [] ((::PHG4Shower*)p);
   }
   static void destruct_PHG4Shower(void *p) {
      typedef ::PHG4Shower current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4Shower

namespace {
  void TriggerDictionaryInitialization_PHG4Shower_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4Shower_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4Shower;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4Shower_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4SHOWER_H
#define G4MAIN_PHG4SHOWER_H

#include "PHG4HitDefs.h"

#include <phool/PHObject.h>

#include <cmath>  // for NAN def
#include <iostream>
#include <map>
#include <set>

class PHG4Shower : public PHObject
{
 public:
  typedef std::set<int> ParticleIdSet;
  typedef ParticleIdSet::iterator ParticleIdIter;
  typedef ParticleIdSet::const_iterator ParticleIdConstIter;

  typedef std::set<int> VertexIdSet;
  typedef VertexIdSet::iterator VertexIdIter;
  typedef VertexIdSet::const_iterator VertexIdConstIter;

  typedef std::map<int, std::set<PHG4HitDefs::keytype> > HitIdMap;
  typedef HitIdMap::iterator HitIdIter;
  typedef HitIdMap::const_iterator HitIdConstIter;

  virtual ~PHG4Shower() {}

  // PHObject virtual overloads

  virtual void identify(std::ostream& os = std::cout) const { os << "PHG4Shower base class" << std::endl; }
  virtual PHG4Shower* CloneMe() const { return nullptr; }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }

  // shower info

  virtual int get_id() const { return 0; }
  virtual void set_id(int id) {}

  virtual int get_parent_particle_id() const { return 0; }
  virtual void set_parent_particle_id(int parent_particle_id) {}

  virtual int get_parent_shower_id() const { return 0; }
  virtual void set_parent_shower_id(int parent_shower_id) {}

  virtual float get_x() const { return NAN; }
  virtual void set_x(float x) {}

  virtual float get_y() const { return NAN; }
  virtual void set_y(float y) {}

  virtual float get_z() const { return NAN; }
  virtual void set_z(float x) {}

  virtual float get_position(unsigned int coor) const { return NAN; }
  virtual void set_position(unsigned int coor, float xi) {}

  virtual float get_covar(unsigned int i, unsigned int j) const { return NAN; }
  virtual void set_covar(unsigned int i, unsigned int j, float entry) {}

  virtual unsigned int get_nhits(int volume) const { return 0; }
  virtual void set_nhits(int volume, unsigned int nhits) {}

  virtual double get_edep() const { return NAN; }
  virtual float get_edep(int volume) const { return NAN; }
  virtual void set_edep(int volume, float edep) {}

  virtual double get_eion() const { return NAN; }
  virtual float get_eion(int volume) const { return NAN; }
  virtual void set_eion(int volume, float eion) {}

  virtual float get_light_yield(int volume) const { return NAN; }
  virtual void set_light_yield(int volume, float light_yield) {}

  virtual float get_eh_ratio(int volume) const { return NAN; }
  virtual void set_eh_ratio(int volume, float eh_ratio) {}

  virtual bool empty_g4particle_id() const { return true; }
  virtual size_t size_g4particle_id() const { return 0; }
  virtual void add_g4particle_id(int id) {}
  virtual ParticleIdIter begin_g4particle_id() { return ParticleIdSet().end(); }
  virtual ParticleIdConstIter begin_g4particle_id() const { return ParticleIdSet().end(); }
  virtual ParticleIdIter end_g4particle_id() { return ParticleIdSet().end(); }
  virtual ParticleIdConstIter end_g4particle_id() const { return ParticleIdSet().end(); }
  virtual size_t remove_g4particle_id(int id) { return 0; }
  virtual void clear_g4particle_id() {}

  virtual bool empty_g4vertex_id() const { return true; }
  virtual size_t size_g4vertex_id() const { return 0; }
  virtual void add_g4vertex_id(int id) {}
  virtual VertexIdIter begin_g4vertex_id() { return VertexIdSet().end(); }
  virtual VertexIdConstIter begin_g4vertex_id() const { return VertexIdSet().end(); }
  virtual VertexIdIter end_g4vertex_id() { return VertexIdSet().end(); }
  virtual VertexIdConstIter end_g4vertex_id() const { return VertexIdSet().end(); }
  virtual size_t remove_g4vertex_id(int id) { return 0; }
  virtual void clear_g4vertex_id() {}

  virtual bool empty_g4hit_id() const { return true; }
  virtual size_t size_g4hit_id() const { return 0; }
  virtual void add_g4hit_id(int volume, PHG4HitDefs::keytype id) {}
  virtual HitIdIter begin_g4hit_id() { return HitIdMap().end(); }
  virtual HitIdConstIter begin_g4hit_id() const { return HitIdMap().end(); }
  virtual HitIdIter find_g4hit_id(int volume) { return HitIdMap().end(); }
  virtual HitIdConstIter find_g4hit_id(int volume) const { return HitIdMap().end(); }
  virtual HitIdIter end_g4hit_id() { return HitIdMap().end(); }
  virtual HitIdConstIter end_g4hit_id() const { return HitIdMap().end(); }
  virtual size_t remove_g4hit_id(int volume, PHG4HitDefs::keytype id) { return 0; }
  virtual size_t remove_g4hit_volume(int volume) { return 0; }
  virtual void clear_g4hit_id() {}

 protected:
  PHG4Shower() {}

 private:
  ClassDef(PHG4Shower, 1);
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4Shower", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4Shower_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4Shower_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4Shower_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4Shower_Dict() {
  TriggerDictionaryInitialization_PHG4Shower_Dict_Impl();
}
