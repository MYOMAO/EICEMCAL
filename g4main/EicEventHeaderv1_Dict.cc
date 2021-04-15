// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME EicEventHeaderv1_Dict
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
#include "EicEventHeaderv1.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_EicEventHeaderv1(void *p = 0);
   static void *newArray_EicEventHeaderv1(Long_t size, void *p);
   static void delete_EicEventHeaderv1(void *p);
   static void deleteArray_EicEventHeaderv1(void *p);
   static void destruct_EicEventHeaderv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EicEventHeaderv1*)
   {
      ::EicEventHeaderv1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::EicEventHeaderv1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("EicEventHeaderv1", ::EicEventHeaderv1::Class_Version(), "", 14,
                  typeid(::EicEventHeaderv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::EicEventHeaderv1::Dictionary, isa_proxy, 4,
                  sizeof(::EicEventHeaderv1) );
      instance.SetNew(&new_EicEventHeaderv1);
      instance.SetNewArray(&newArray_EicEventHeaderv1);
      instance.SetDelete(&delete_EicEventHeaderv1);
      instance.SetDeleteArray(&deleteArray_EicEventHeaderv1);
      instance.SetDestructor(&destruct_EicEventHeaderv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EicEventHeaderv1*)
   {
      return GenerateInitInstanceLocal((::EicEventHeaderv1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::EicEventHeaderv1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr EicEventHeaderv1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *EicEventHeaderv1::Class_Name()
{
   return "EicEventHeaderv1";
}

//______________________________________________________________________________
const char *EicEventHeaderv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EicEventHeaderv1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int EicEventHeaderv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EicEventHeaderv1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *EicEventHeaderv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EicEventHeaderv1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *EicEventHeaderv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EicEventHeaderv1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void EicEventHeaderv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class EicEventHeaderv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(EicEventHeaderv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(EicEventHeaderv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_EicEventHeaderv1(void *p) {
      return  p ? new(p) ::EicEventHeaderv1 : new ::EicEventHeaderv1;
   }
   static void *newArray_EicEventHeaderv1(Long_t nElements, void *p) {
      return p ? new(p) ::EicEventHeaderv1[nElements] : new ::EicEventHeaderv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_EicEventHeaderv1(void *p) {
      delete ((::EicEventHeaderv1*)p);
   }
   static void deleteArray_EicEventHeaderv1(void *p) {
      delete [] ((::EicEventHeaderv1*)p);
   }
   static void destruct_EicEventHeaderv1(void *p) {
      typedef ::EicEventHeaderv1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::EicEventHeaderv1

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
  void TriggerDictionaryInitialization_EicEventHeaderv1_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "EicEventHeaderv1_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class EicEventHeaderv1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "EicEventHeaderv1_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICEVENTHEADERV1_H
#define EICEVENTHEADERV1_H

#include "EicEventHeader.h"

#include <cstdint>
#include <iostream>
#include <map>

class EicEventHeaderv1 : public EicEventHeader
{
 public:
  EicEventHeaderv1() {}
  explicit EicEventHeaderv1(const EicEventHeader *eicevt);
  virtual ~EicEventHeaderv1() {}

  //  void identify(std::ostream& os  = std::cout) const;
  void Reset();

  // property relates methods
  bool has_property(const PROPERTY prop_id) const;
  float get_property_float(const PROPERTY prop_id) const;
  int get_property_int(const PROPERTY prop_id) const;
  unsigned int get_property_uint(const PROPERTY prop_id) const;
  void set_property(const PROPERTY prop_id, const float value);
  void set_property(const PROPERTY prop_id, const int value);
  void set_property(const PROPERTY prop_id, const unsigned int value);

  // Generator specific values
  void set_eventgenerator_type(const int i) { set_property(prop_eventgen, i); }
  int get_eventgenerator_type() const { return get_property_int(prop_eventgen); }

  // Milou
  void set_milou_weight(const float val) { set_property(prop_milou_weight, val); }
  float get_milou_weight() const { return get_property_float(prop_milou_weight); }
  void set_milou_trueX(const float val) { set_property(prop_milou_truex, val); }
  float get_milou_trueX() const { return get_property_float(prop_milou_truex); }
  void set_milou_trueQ2(const float val) { set_property(prop_milou_trueq2, val); }
  float get_milou_trueQ2() const { return get_property_float(prop_milou_trueq2); }

 protected:
  unsigned int get_property_nocheck(const PROPERTY prop_id) const;
  void set_property_nocheck(const PROPERTY prop_id, const unsigned int ui) { prop_map[prop_id] = ui; }

  //! storage types for properties
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

  //! container for  properties
  prop_map_t prop_map;

  ClassDef(EicEventHeaderv1, 1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"EicEventHeaderv1", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("EicEventHeaderv1_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_EicEventHeaderv1_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_EicEventHeaderv1_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_EicEventHeaderv1_Dict() {
  TriggerDictionaryInitialization_EicEventHeaderv1_Dict_Impl();
}
