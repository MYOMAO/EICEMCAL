// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME EicEventHeader_Dict
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
#include "EicEventHeader.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_EicEventHeader(void *p = 0);
   static void *newArray_EicEventHeader(Long_t size, void *p);
   static void delete_EicEventHeader(void *p);
   static void deleteArray_EicEventHeader(void *p);
   static void destruct_EicEventHeader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EicEventHeader*)
   {
      ::EicEventHeader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::EicEventHeader >(0);
      static ::ROOT::TGenericClassInfo 
         instance("EicEventHeader", ::EicEventHeader::Class_Version(), "", 17,
                  typeid(::EicEventHeader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::EicEventHeader::Dictionary, isa_proxy, 4,
                  sizeof(::EicEventHeader) );
      instance.SetNew(&new_EicEventHeader);
      instance.SetNewArray(&newArray_EicEventHeader);
      instance.SetDelete(&delete_EicEventHeader);
      instance.SetDeleteArray(&deleteArray_EicEventHeader);
      instance.SetDestructor(&destruct_EicEventHeader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EicEventHeader*)
   {
      return GenerateInitInstanceLocal((::EicEventHeader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::EicEventHeader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr EicEventHeader::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *EicEventHeader::Class_Name()
{
   return "EicEventHeader";
}

//______________________________________________________________________________
const char *EicEventHeader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EicEventHeader*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int EicEventHeader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EicEventHeader*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *EicEventHeader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EicEventHeader*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *EicEventHeader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EicEventHeader*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void EicEventHeader::Streamer(TBuffer &R__b)
{
   // Stream an object of class EicEventHeader.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(EicEventHeader::Class(),this);
   } else {
      R__b.WriteClassBuffer(EicEventHeader::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_EicEventHeader(void *p) {
      return  p ? new(p) ::EicEventHeader : new ::EicEventHeader;
   }
   static void *newArray_EicEventHeader(Long_t nElements, void *p) {
      return p ? new(p) ::EicEventHeader[nElements] : new ::EicEventHeader[nElements];
   }
   // Wrapper around operator delete
   static void delete_EicEventHeader(void *p) {
      delete ((::EicEventHeader*)p);
   }
   static void deleteArray_EicEventHeader(void *p) {
      delete [] ((::EicEventHeader*)p);
   }
   static void destruct_EicEventHeader(void *p) {
      typedef ::EicEventHeader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::EicEventHeader

namespace ROOT {
   static TClass *maplEstringcOdoublegR_Dictionary();
   static void maplEstringcOdoublegR_TClassManip(TClass*);
   static void *new_maplEstringcOdoublegR(void *p = 0);
   static void *newArray_maplEstringcOdoublegR(Long_t size, void *p);
   static void delete_maplEstringcOdoublegR(void *p);
   static void deleteArray_maplEstringcOdoublegR(void *p);
   static void destruct_maplEstringcOdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,double>*)
   {
      map<string,double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,double>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,double>", -2, "map", 100,
                  typeid(map<string,double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,double>) );
      instance.SetNew(&new_maplEstringcOdoublegR);
      instance.SetNewArray(&newArray_maplEstringcOdoublegR);
      instance.SetDelete(&delete_maplEstringcOdoublegR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOdoublegR);
      instance.SetDestructor(&destruct_maplEstringcOdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,double> >()));

      ::ROOT::AddClassAlternate("map<string,double>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, double, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, double> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<string,double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<string,double>*)0x0)->GetClass();
      maplEstringcOdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,double> : new map<string,double>;
   }
   static void *newArray_maplEstringcOdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,double>[nElements] : new map<string,double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOdoublegR(void *p) {
      delete ((map<string,double>*)p);
   }
   static void deleteArray_maplEstringcOdoublegR(void *p) {
      delete [] ((map<string,double>*)p);
   }
   static void destruct_maplEstringcOdoublegR(void *p) {
      typedef map<string,double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<string,double>

namespace {
  void TriggerDictionaryInitialization_EicEventHeader_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "EicEventHeader_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class EicEventHeader;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "EicEventHeader_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICEVENTHEADER_H
#define EICEVENTHEADER_H

#include <phool/PHObject.h>

#include <climits>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <utility>

class EicEventHeader : public PHObject
{
 public:
  EicEventHeader() {}
  virtual ~EicEventHeader();

  virtual void identify(std::ostream& os = std::cout) const;
  virtual void CopyFrom(const PHObject* phobj);

  virtual void Reset();

  virtual void set_eventgenerator_type(const int i) { return; }
  virtual int get_eventgenerator_type() const { return -99999; }

  // Milou

  virtual void set_milou_weight(const float val) { return; }
  virtual float get_milou_weight() const { return NAN; }
  virtual void set_milou_trueX(const float val) { return; }
  virtual float get_milou_trueX() const { return NAN; }
  virtual void set_milou_trueQ2(const float val) { return; }
  virtual float get_milou_trueQ2() const { return NAN; }

  //! Procedure to add a new PROPERTY tag:
  //! 1.add new tag below with unique value,
  //! 2.add a short name to EicEventHeader::get_property_info
  enum PROPERTY
  {  //
    prop_eventgen = 1,
    prop_milou_weight = 2,
    prop_milou_truex = 3,
    prop_milou_trueq2 = 4,

    //! max limit in order to fit into 8 bit unsigned number
    prop_MAX_NUMBER = UCHAR_MAX
  };

  enum PROPERTY_TYPE
  {  //
    type_int = 1,
    type_uint = 2,
    type_float = 3,
    type_unknown = -1
  };

  enum EvtGen
  {
    Milou = 1
  };

  virtual bool has_property(const PROPERTY prop_id) const { return false; }
  virtual float get_property_float(const PROPERTY prop_id) const { return NAN; }
  virtual int get_property_int(const PROPERTY prop_id) const { return INT_MIN; }
  virtual unsigned int get_property_uint(const PROPERTY prop_id) const { return UINT_MAX; }
  virtual void set_property(const PROPERTY prop_id, const float value) { return; }
  virtual void set_property(const PROPERTY prop_id, const int value) { return; }
  virtual void set_property(const PROPERTY prop_id, const unsigned int value) { return; }
  static std::pair<const std::string, PROPERTY_TYPE> get_property_info(PROPERTY prop_id);
  static bool check_property(const PROPERTY prop_id, const PROPERTY_TYPE prop_type);
  static std::string get_property_type(const PROPERTY_TYPE prop_type);

 protected:
  virtual unsigned int get_property_nocheck(const PROPERTY prop_id) const { return UINT_MAX; }
  virtual void set_property_nocheck(const PROPERTY prop_id, const unsigned int) { return; }

  std::map<std::string, double> evInfo;

  ClassDef(EicEventHeader, 1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"EicEventHeader", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("EicEventHeader_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_EicEventHeader_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_EicEventHeader_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_EicEventHeader_Dict() {
  TriggerDictionaryInitialization_EicEventHeader_Dict_Impl();
}
