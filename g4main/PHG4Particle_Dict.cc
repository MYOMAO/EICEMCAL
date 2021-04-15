// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4Particle_Dict
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
#include "PHG4Particle.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4Particle(void *p = 0);
   static void *newArray_PHG4Particle(Long_t size, void *p);
   static void delete_PHG4Particle(void *p);
   static void deleteArray_PHG4Particle(void *p);
   static void destruct_PHG4Particle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4Particle*)
   {
      ::PHG4Particle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4Particle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4Particle", ::PHG4Particle::Class_Version(), "", 14,
                  typeid(::PHG4Particle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4Particle::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4Particle) );
      instance.SetNew(&new_PHG4Particle);
      instance.SetNewArray(&newArray_PHG4Particle);
      instance.SetDelete(&delete_PHG4Particle);
      instance.SetDeleteArray(&deleteArray_PHG4Particle);
      instance.SetDestructor(&destruct_PHG4Particle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4Particle*)
   {
      return GenerateInitInstanceLocal((::PHG4Particle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4Particle*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4Particle::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4Particle::Class_Name()
{
   return "PHG4Particle";
}

//______________________________________________________________________________
const char *PHG4Particle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Particle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4Particle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Particle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4Particle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Particle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4Particle::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Particle*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4Particle::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4Particle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4Particle::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4Particle::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4Particle(void *p) {
      return  p ? new(p) ::PHG4Particle : new ::PHG4Particle;
   }
   static void *newArray_PHG4Particle(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4Particle[nElements] : new ::PHG4Particle[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4Particle(void *p) {
      delete ((::PHG4Particle*)p);
   }
   static void deleteArray_PHG4Particle(void *p) {
      delete [] ((::PHG4Particle*)p);
   }
   static void destruct_PHG4Particle(void *p) {
      typedef ::PHG4Particle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4Particle

namespace {
  void TriggerDictionaryInitialization_PHG4Particle_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4Particle_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4Particle;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4Particle_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PARTICLE_H
#define G4MAIN_PHG4PARTICLE_H

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>
#include <string>

class PHG4Particle : public PHObject
{
 public:
  PHG4Particle() {}
  virtual ~PHG4Particle() {}
  virtual bool isIon() const { return false; }
  virtual int get_pid() const { return 0; }
  virtual std::string get_name() const { return "NONE"; }
  virtual double get_px() const { return NAN; }
  virtual double get_py() const { return NAN; }
  virtual double get_pz() const { return NAN; }
  virtual double get_e() const { return NAN; }

  virtual int get_track_id() const { return -9999; }
  virtual int get_vtx_id() const { return -9999; }
  virtual int get_parent_id() const { return -9999; }
  virtual int get_primary_id() const { return 0xFFFFFFFF; }

  virtual int get_barcode() const { return 0xFFFFFFFF; }

  virtual int get_A() const { return 0xFFFFFFFF; }
  virtual int get_Z() const { return 0xFFFFFFFF; }
  virtual double get_IonCharge() const { return NAN; }
  virtual double get_ExcitEnergy() const { return NAN; }

  virtual void set_track_id(const int i) { return; }
  virtual void set_vtx_id(const int i) { return; }
  virtual void set_parent_id(const int i) { return; }
  virtual void set_primary_id(const int i) { return; }
  virtual void set_name(const std::string &name) { return; }
  virtual void set_pid(const int i) { return; }
  virtual void set_px(const double x) { return; }
  virtual void set_py(const double x) { return; }
  virtual void set_pz(const double x) { return; }
  virtual void set_e(const double e) { return; }

  virtual void set_barcode(const int barcode) { return; }

  virtual void set_A(const int a) { return; }
  virtual void set_Z(const int z) { return; }
  virtual void set_NumCharge(const int c) { return; }
  virtual void set_IonCharge(const double e) { return; }
  virtual void set_ExcitEnergy(const double e) { return; }

  void identify(std::ostream &os = std::cout) const;

  bool operator==(const PHG4Particle &p) const;

 protected:
  ClassDef(PHG4Particle, 1)
};

#endif  // G4MAIN_PHG4PARTICLE_H

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4Particle", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4Particle_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4Particle_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4Particle_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4Particle_Dict() {
  TriggerDictionaryInitialization_PHG4Particle_Dict_Impl();
}
