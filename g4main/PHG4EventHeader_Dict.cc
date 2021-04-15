// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4EventHeader_Dict
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
#include "PHG4EventHeader.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4EventHeader(void *p = 0);
   static void *newArray_PHG4EventHeader(Long_t size, void *p);
   static void delete_PHG4EventHeader(void *p);
   static void deleteArray_PHG4EventHeader(void *p);
   static void destruct_PHG4EventHeader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4EventHeader*)
   {
      ::PHG4EventHeader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4EventHeader >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4EventHeader", ::PHG4EventHeader::Class_Version(), "", 15,
                  typeid(::PHG4EventHeader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4EventHeader::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4EventHeader) );
      instance.SetNew(&new_PHG4EventHeader);
      instance.SetNewArray(&newArray_PHG4EventHeader);
      instance.SetDelete(&delete_PHG4EventHeader);
      instance.SetDeleteArray(&deleteArray_PHG4EventHeader);
      instance.SetDestructor(&destruct_PHG4EventHeader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4EventHeader*)
   {
      return GenerateInitInstanceLocal((::PHG4EventHeader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4EventHeader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4EventHeader::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4EventHeader::Class_Name()
{
   return "PHG4EventHeader";
}

//______________________________________________________________________________
const char *PHG4EventHeader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4EventHeader*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4EventHeader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4EventHeader*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4EventHeader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4EventHeader*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4EventHeader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4EventHeader*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4EventHeader::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4EventHeader.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4EventHeader::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4EventHeader::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4EventHeader(void *p) {
      return  p ? new(p) ::PHG4EventHeader : new ::PHG4EventHeader;
   }
   static void *newArray_PHG4EventHeader(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4EventHeader[nElements] : new ::PHG4EventHeader[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4EventHeader(void *p) {
      delete ((::PHG4EventHeader*)p);
   }
   static void deleteArray_PHG4EventHeader(void *p) {
      delete [] ((::PHG4EventHeader*)p);
   }
   static void destruct_PHG4EventHeader(void *p) {
      typedef ::PHG4EventHeader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4EventHeader

namespace {
  void TriggerDictionaryInitialization_PHG4EventHeader_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4EventHeader_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4EventHeader;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4EventHeader_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4EVENTHEADER_H
#define G4MAIN_PHG4EVENTHEADER_H

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <cmath>
#include <iostream>

///
class PHG4EventHeader: public PHObject
{
 public:

  /// dtor
  virtual ~PHG4EventHeader() {}

  /// Clear Event
  virtual void Reset()
    {
      std::cout << PHWHERE << "ERROR Reset() not implemented by daughter class" << std::endl;
      return;
    }

  /** identify Function from PHObject
      @param os Output Stream 
   */
  virtual void identify(std::ostream& os = std::cout) const
    {
      os << "identify yourself: virtual PHG4EventHeader Object" << std::endl;
      return;
    }

  /// isValid returns non zero if object contains valid data
  virtual int isValid() const
    {
      std::cout << PHWHERE << "isValid not implemented by daughter class" << std::endl;
      return 0;
    }

  /// get Event Number
  virtual int get_EvtSequence() const {return -9999;}
  /// set Event Number
  virtual void set_EvtSequence(const int /*ival*/) {return;}

  virtual float get_ImpactParameter() const {return NAN;}
  virtual void set_ImpactParameter(const float) {return;}

  virtual float get_EventPlaneAngle() const {return NAN;}
  virtual void set_EventPlaneAngle(const float) {return;}


 private: // prevent doc++ from showing ClassDef
  ClassDef(PHG4EventHeader,1)

};

#endif




#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4EventHeader", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4EventHeader_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4EventHeader_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4EventHeader_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4EventHeader_Dict() {
  TriggerDictionaryInitialization_PHG4EventHeader_Dict_Impl();
}
