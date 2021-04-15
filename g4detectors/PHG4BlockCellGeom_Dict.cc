// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4BlockCellGeom_Dict
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
#include "PHG4BlockCellGeom.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4BlockCellGeom(void *p = 0);
   static void *newArray_PHG4BlockCellGeom(Long_t size, void *p);
   static void delete_PHG4BlockCellGeom(void *p);
   static void deleteArray_PHG4BlockCellGeom(void *p);
   static void destruct_PHG4BlockCellGeom(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4BlockCellGeom*)
   {
      ::PHG4BlockCellGeom *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4BlockCellGeom >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4BlockCellGeom", ::PHG4BlockCellGeom::Class_Version(), "", 17,
                  typeid(::PHG4BlockCellGeom), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4BlockCellGeom::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4BlockCellGeom) );
      instance.SetNew(&new_PHG4BlockCellGeom);
      instance.SetNewArray(&newArray_PHG4BlockCellGeom);
      instance.SetDelete(&delete_PHG4BlockCellGeom);
      instance.SetDeleteArray(&deleteArray_PHG4BlockCellGeom);
      instance.SetDestructor(&destruct_PHG4BlockCellGeom);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4BlockCellGeom*)
   {
      return GenerateInitInstanceLocal((::PHG4BlockCellGeom*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4BlockCellGeom*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4BlockCellGeom::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4BlockCellGeom::Class_Name()
{
   return "PHG4BlockCellGeom";
}

//______________________________________________________________________________
const char *PHG4BlockCellGeom::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockCellGeom*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4BlockCellGeom::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockCellGeom*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4BlockCellGeom::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockCellGeom*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4BlockCellGeom::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4BlockCellGeom*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4BlockCellGeom::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4BlockCellGeom.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4BlockCellGeom::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4BlockCellGeom::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4BlockCellGeom(void *p) {
      return  p ? new(p) ::PHG4BlockCellGeom : new ::PHG4BlockCellGeom;
   }
   static void *newArray_PHG4BlockCellGeom(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4BlockCellGeom[nElements] : new ::PHG4BlockCellGeom[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4BlockCellGeom(void *p) {
      delete ((::PHG4BlockCellGeom*)p);
   }
   static void deleteArray_PHG4BlockCellGeom(void *p) {
      delete [] ((::PHG4BlockCellGeom*)p);
   }
   static void destruct_PHG4BlockCellGeom(void *p) {
      typedef ::PHG4BlockCellGeom current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4BlockCellGeom

namespace {
  void TriggerDictionaryInitialization_PHG4BlockCellGeom_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4BlockCellGeom_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4BlockCellGeom;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4BlockCellGeom_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BLOCKCELLGEOM_H
#define G4DETECTORS_PHG4BLOCKCELLGEOM_H

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream
#include <string>
#include <utility>           // for pair

class PHG4BlockCellGeom: public PHObject
{
 public:
  PHG4BlockCellGeom();

  virtual ~PHG4BlockCellGeom() {}

  void identify(std::ostream& os = std::cout) const;
  int get_layer() const {return _layer;}
  int get_binning() const {return _binning;}

  int get_zbins() const;
  double get_zstep() const;
  double get_zmin() const;

  int get_xbins() const;
  double get_xstep() const;
  double get_xmin() const;

  int get_etabins() const;
  double get_etastep() const;
  double get_etamin() const;

  std::pair<double, double> get_zbounds(const int ibin) const;
  std::pair<double, double> get_etabounds(const int ibin) const;
  std::pair<double, double> get_xbounds(const int ibin) const;

  double get_etacenter(const int ibin) const;
  double get_zcenter(const int ibin) const;
  double get_xcenter(const int ibin) const;

  int get_etabin(const double eta) const;
  int get_zbin(const double z) const;
  int get_xbin(const double x) const;

  void set_layer(const int i) {_layer = i;}
  void set_binning(const int i) {_binning = i;}

  void set_zbins(const int i);
  void set_zmin(const double z);
  void set_zstep(const double z);
  
  void set_xbins(const int i);
  void set_xstep(const double x);
  void set_xmin(const double x);
  
  void set_etabins(const int i);
  void set_etamin(const double z);
  void set_etastep(const double z);
  
 protected:
  void check_binning_method(const int i) const;
  void check_binning_method_eta(const std::string & src = "") const;
  void check_binning_method_x(const std::string & src = "") const;
  std::string methodname(const int i) const;

  int _layer;
  int _binning;

  int _nzbins;
  double _zmin;
  double _zstep;

  int _nxbins;
  double _xmin;
  double _xstep;

  ClassDef(PHG4BlockCellGeom,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4BlockCellGeom", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4BlockCellGeom_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4BlockCellGeom_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4BlockCellGeom_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4BlockCellGeom_Dict() {
  TriggerDictionaryInitialization_PHG4BlockCellGeom_Dict_Impl();
}
