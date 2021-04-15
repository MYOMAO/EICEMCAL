// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4Hitv1_Dict
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
#include "PHG4Hitv1.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4Hitv1(void *p = 0);
   static void *newArray_PHG4Hitv1(Long_t size, void *p);
   static void delete_PHG4Hitv1(void *p);
   static void deleteArray_PHG4Hitv1(void *p);
   static void destruct_PHG4Hitv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4Hitv1*)
   {
      ::PHG4Hitv1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4Hitv1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4Hitv1", ::PHG4Hitv1::Class_Version(), "", 17,
                  typeid(::PHG4Hitv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4Hitv1::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4Hitv1) );
      instance.SetNew(&new_PHG4Hitv1);
      instance.SetNewArray(&newArray_PHG4Hitv1);
      instance.SetDelete(&delete_PHG4Hitv1);
      instance.SetDeleteArray(&deleteArray_PHG4Hitv1);
      instance.SetDestructor(&destruct_PHG4Hitv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4Hitv1*)
   {
      return GenerateInitInstanceLocal((::PHG4Hitv1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4Hitv1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4Hitv1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4Hitv1::Class_Name()
{
   return "PHG4Hitv1";
}

//______________________________________________________________________________
const char *PHG4Hitv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Hitv1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4Hitv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Hitv1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4Hitv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Hitv1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4Hitv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Hitv1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4Hitv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4Hitv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4Hitv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4Hitv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4Hitv1(void *p) {
      return  p ? new(p) ::PHG4Hitv1 : new ::PHG4Hitv1;
   }
   static void *newArray_PHG4Hitv1(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4Hitv1[nElements] : new ::PHG4Hitv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4Hitv1(void *p) {
      delete ((::PHG4Hitv1*)p);
   }
   static void deleteArray_PHG4Hitv1(void *p) {
      delete [] ((::PHG4Hitv1*)p);
   }
   static void destruct_PHG4Hitv1(void *p) {
      typedef ::PHG4Hitv1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4Hitv1

namespace ROOT {
   static TClass *maplEunsignedsPcharcOunsignedsPintgR_Dictionary();
   static void maplEunsignedsPcharcOunsignedsPintgR_TClassManip(TClass*);
   static void *new_maplEunsignedsPcharcOunsignedsPintgR(void *p = 0);
   static void *newArray_maplEunsignedsPcharcOunsignedsPintgR(Long_t size, void *p);
   static void delete_maplEunsignedsPcharcOunsignedsPintgR(void *p);
   static void deleteArray_maplEunsignedsPcharcOunsignedsPintgR(void *p);
   static void destruct_maplEunsignedsPcharcOunsignedsPintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<unsigned char,unsigned int>*)
   {
      map<unsigned char,unsigned int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<unsigned char,unsigned int>));
      static ::ROOT::TGenericClassInfo 
         instance("map<unsigned char,unsigned int>", -2, "map", 100,
                  typeid(map<unsigned char,unsigned int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEunsignedsPcharcOunsignedsPintgR_Dictionary, isa_proxy, 0,
                  sizeof(map<unsigned char,unsigned int>) );
      instance.SetNew(&new_maplEunsignedsPcharcOunsignedsPintgR);
      instance.SetNewArray(&newArray_maplEunsignedsPcharcOunsignedsPintgR);
      instance.SetDelete(&delete_maplEunsignedsPcharcOunsignedsPintgR);
      instance.SetDeleteArray(&deleteArray_maplEunsignedsPcharcOunsignedsPintgR);
      instance.SetDestructor(&destruct_maplEunsignedsPcharcOunsignedsPintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<unsigned char,unsigned int> >()));

      ::ROOT::AddClassAlternate("map<unsigned char,unsigned int>","std::map<unsigned char, unsigned int, std::less<unsigned char>, std::allocator<std::pair<unsigned char const, unsigned int> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<unsigned char,unsigned int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEunsignedsPcharcOunsignedsPintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<unsigned char,unsigned int>*)0x0)->GetClass();
      maplEunsignedsPcharcOunsignedsPintgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEunsignedsPcharcOunsignedsPintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEunsignedsPcharcOunsignedsPintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<unsigned char,unsigned int> : new map<unsigned char,unsigned int>;
   }
   static void *newArray_maplEunsignedsPcharcOunsignedsPintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<unsigned char,unsigned int>[nElements] : new map<unsigned char,unsigned int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEunsignedsPcharcOunsignedsPintgR(void *p) {
      delete ((map<unsigned char,unsigned int>*)p);
   }
   static void deleteArray_maplEunsignedsPcharcOunsignedsPintgR(void *p) {
      delete [] ((map<unsigned char,unsigned int>*)p);
   }
   static void destruct_maplEunsignedsPcharcOunsignedsPintgR(void *p) {
      typedef map<unsigned char,unsigned int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<unsigned char,unsigned int>

namespace {
  void TriggerDictionaryInitialization_PHG4Hitv1_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4Hitv1_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4Hitv1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4Hitv1_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4HITV1_H
#define G4MAIN_PHG4HITV1_H

#include "PHG4Hit.h"
#include "PHG4HitDefs.h"

#include <climits>  // for INT_MIN, ULONG_LONG_MAX
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>

class PHG4Hitv1 : public PHG4Hit
{
 public:
  PHG4Hitv1() = default;
  explicit PHG4Hitv1(const PHG4Hit* g4hit);
  virtual ~PHG4Hitv1() = default;
  void identify(std::ostream& os = std::cout) const;
  void Reset();

  // The indices here represent the entry and exit points of the particle
  float get_x(const int i) const { return x[i]; }
  float get_y(const int i) const { return y[i]; }
  float get_z(const int i) const { return z[i]; }
  float get_t(const int i) const { return t[i]; }
  float get_edep() const { return edep; }
  PHG4HitDefs::keytype get_hit_id() const { return hitid; }
  int get_detid() const;
  int get_shower_id() const { return showerid; }
  int get_trkid() const { return trackid; }

  void set_x(const int i, const float f) { x[i] = f; }
  void set_y(const int i, const float f) { y[i] = f; }
  void set_z(const int i, const float f) { z[i] = f; }
  void set_t(const int i, const float f) { t[i] = f; }
  void set_edep(const float f) { edep = f; }
  void set_hit_id(const PHG4HitDefs::keytype i) { hitid = i; }
  void set_shower_id(const int i) { showerid = i; }
  void set_trkid(const int i) { trackid = i; }

  virtual void print() const;

  bool has_property(const PROPERTY prop_id) const;
  float get_property_float(const PROPERTY prop_id) const;
  int get_property_int(const PROPERTY prop_id) const;
  unsigned int get_property_uint(const PROPERTY prop_id) const;
  void set_property(const PROPERTY prop_id, const float value);
  void set_property(const PROPERTY prop_id, const int value);
  void set_property(const PROPERTY prop_id, const unsigned int value);

  virtual float get_px(const int i) const;
  virtual float get_py(const int i) const;
  virtual float get_pz(const int i) const;
  virtual float get_local_x(const int i) const;
  virtual float get_local_y(const int i) const;
  virtual float get_local_z(const int i) const;
  virtual float get_eion() const { return get_property_float(prop_eion); }
  virtual float get_light_yield() const { return get_property_float(prop_light_yield); }
  virtual float get_path_length() const { return get_property_float(prop_path_length); }
  virtual unsigned int get_layer() const { return get_property_uint(prop_layer); }
  virtual int get_scint_id() const { return get_property_int(prop_scint_id); }
  virtual int get_row() const { return get_property_int(prop_row); }
  virtual int get_strip_z_index() const { return get_property_int(prop_strip_z_index); }
  virtual int get_strip_y_index() const { return get_property_int(prop_strip_y_index); }
  virtual int get_ladder_z_index() const { return get_property_int(prop_ladder_z_index); }
  virtual int get_ladder_phi_index() const { return get_property_int(prop_ladder_phi_index); }
  virtual int get_index_i() const { return get_property_int(prop_index_i); }
  virtual int get_index_j() const { return get_property_int(prop_index_j); }
  virtual int get_index_k() const { return get_property_int(prop_index_k); }
  virtual int get_index_l() const { return get_property_int(prop_index_l); }
  virtual int get_hit_type() const { return get_property_int(prop_hit_type); }

  virtual void set_px(const int i, const float f);
  virtual void set_py(const int i, const float f);
  virtual void set_pz(const int i, const float f);
  virtual void set_local_x(const int i, const float f);
  virtual void set_local_y(const int i, const float f);
  virtual void set_local_z(const int i, const float f);
  virtual void set_eion(const float f) { set_property(prop_eion, f); }
  virtual void set_light_yield(const float f) { set_property(prop_light_yield, f); }
  virtual void set_path_length(const float f) { set_property(prop_path_length, f); }
  virtual void set_layer(const unsigned int i) { set_property(prop_layer, i); }
  virtual void set_scint_id(const int i) { set_property(prop_scint_id, i); }
  virtual void set_row(const int i) { set_property(prop_row, i); }
  virtual void set_strip_z_index(const int i) { set_property(prop_strip_z_index, i); }
  virtual void set_strip_y_index(const int i) { set_property(prop_strip_y_index, i); }
  virtual void set_ladder_z_index(const int i) { set_property(prop_ladder_z_index, i); }
  virtual void set_ladder_phi_index(const int i) { set_property(prop_ladder_phi_index, i); }
  virtual void set_index_i(const int i) { set_property(prop_index_i, i); }
  virtual void set_index_j(const int i) { set_property(prop_index_j, i); }
  virtual void set_index_k(const int i) { set_property(prop_index_k, i); }
  virtual void set_index_l(const int i) { set_property(prop_index_l, i); }
  virtual void set_hit_type(const int i) { set_property(prop_hit_type, i); }

 protected:
  unsigned int get_property_nocheck(const PROPERTY prop_id) const;
  void set_property_nocheck(const PROPERTY prop_id, const unsigned int ui) { prop_map[prop_id] = ui; }
  // Store both the entry and exit points of the particle
  // Remember, particles do not always enter on the inner edge!
  float x[2] = {NAN, NAN};
  float y[2] = {NAN, NAN};
  float z[2] = {NAN, NAN};
  float t[2] = {NAN, NAN};
  PHG4HitDefs::keytype hitid = ULONG_LONG_MAX;
  int trackid = INT_MIN;
  int showerid = INT_MIN;
  float edep = NAN;

  //! storage types for additional property
  typedef uint8_t prop_id_t;
  typedef uint32_t prop_storage_t;
  typedef std::map<prop_id_t, prop_storage_t> prop_map_t;

  //! convert between 32bit inputs and storage type prop_storage_t
  union u_property {
    float fdata;
    int32_t idata;
    uint32_t uidata;

    u_property(int32_t in)
      : idata(in)
    {
    }
    u_property(uint32_t in)
      : uidata(in)
    {
    }
    u_property(float in)
      : fdata(in)
    {
    }
    u_property()
      : uidata(0)
    {
    }

    prop_storage_t get_storage() const { return uidata; }
  };

  //! container for additional property
  prop_map_t prop_map;

  ClassDef(PHG4Hitv1, 2)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4Hitv1", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4Hitv1_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4Hitv1_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4Hitv1_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4Hitv1_Dict() {
  TriggerDictionaryInitialization_PHG4Hitv1_Dict_Impl();
}
