// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4CylinderCellv3_Dict
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
#include "PHG4CylinderCellv3.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4CylinderCellv3(void *p = 0);
   static void *newArray_PHG4CylinderCellv3(Long_t size, void *p);
   static void delete_PHG4CylinderCellv3(void *p);
   static void deleteArray_PHG4CylinderCellv3(void *p);
   static void destruct_PHG4CylinderCellv3(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderCellv3*)
   {
      ::PHG4CylinderCellv3 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderCellv3 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderCellv3", ::PHG4CylinderCellv3::Class_Version(), "", 15,
                  typeid(::PHG4CylinderCellv3), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderCellv3::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderCellv3) );
      instance.SetNew(&new_PHG4CylinderCellv3);
      instance.SetNewArray(&newArray_PHG4CylinderCellv3);
      instance.SetDelete(&delete_PHG4CylinderCellv3);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderCellv3);
      instance.SetDestructor(&destruct_PHG4CylinderCellv3);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderCellv3*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderCellv3*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderCellv3*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderCellv3::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderCellv3::Class_Name()
{
   return "PHG4CylinderCellv3";
}

//______________________________________________________________________________
const char *PHG4CylinderCellv3::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellv3*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderCellv3::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellv3*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderCellv3::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellv3*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderCellv3::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellv3*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4CylinderCellv3::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderCellv3.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderCellv3::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderCellv3::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4CylinderCellv3(void *p) {
      return  p ? new(p) ::PHG4CylinderCellv3 : new ::PHG4CylinderCellv3;
   }
   static void *newArray_PHG4CylinderCellv3(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4CylinderCellv3[nElements] : new ::PHG4CylinderCellv3[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4CylinderCellv3(void *p) {
      delete ((::PHG4CylinderCellv3*)p);
   }
   static void deleteArray_PHG4CylinderCellv3(void *p) {
      delete [] ((::PHG4CylinderCellv3*)p);
   }
   static void destruct_PHG4CylinderCellv3(void *p) {
      typedef ::PHG4CylinderCellv3 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderCellv3

namespace {
  void TriggerDictionaryInitialization_PHG4CylinderCellv3_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4CylinderCellv3_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4CylinderCellv3;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4CylinderCellv3_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERCELLV3_H
#define G4DETECTORS_PHG4CYLINDERCELLV3_H

#include "PHG4CylinderCellv1.h"

#include <iostream>

class PHG4CylinderCellv3 : public PHG4CylinderCellv1
{
 public:

  PHG4CylinderCellv3();
  virtual ~PHG4CylinderCellv3(){}

  void identify(std::ostream& os = std::cout) const;

  void set_j_index(const int i) {j_index = i;}
  int get_j_index() const {return j_index;}

  void set_k_index(const int i) {k_index = i;}
  int get_k_index() const {return k_index;}

  void set_l_index(const int i) {l_index = i;}
  int get_l_index() const {return l_index;}
  
 protected:

  int j_index;
  int k_index; 
  int l_index; 
  
  ClassDef(PHG4CylinderCellv3,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4CylinderCellv3", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4CylinderCellv3_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4CylinderCellv3_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4CylinderCellv3_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4CylinderCellv3_Dict() {
  TriggerDictionaryInitialization_PHG4CylinderCellv3_Dict_Impl();
}
