// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4CylinderCellGeom_Dict
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
#include "PHG4CylinderCellGeom.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4CylinderCellGeom(void *p = 0);
   static void *newArray_PHG4CylinderCellGeom(Long_t size, void *p);
   static void delete_PHG4CylinderCellGeom(void *p);
   static void deleteArray_PHG4CylinderCellGeom(void *p);
   static void destruct_PHG4CylinderCellGeom(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderCellGeom*)
   {
      ::PHG4CylinderCellGeom *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderCellGeom >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderCellGeom", ::PHG4CylinderCellGeom::Class_Version(), "", 17,
                  typeid(::PHG4CylinderCellGeom), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderCellGeom::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderCellGeom) );
      instance.SetNew(&new_PHG4CylinderCellGeom);
      instance.SetNewArray(&newArray_PHG4CylinderCellGeom);
      instance.SetDelete(&delete_PHG4CylinderCellGeom);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderCellGeom);
      instance.SetDestructor(&destruct_PHG4CylinderCellGeom);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderCellGeom*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderCellGeom*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderCellGeom*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderCellGeom::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderCellGeom::Class_Name()
{
   return "PHG4CylinderCellGeom";
}

//______________________________________________________________________________
const char *PHG4CylinderCellGeom::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellGeom*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderCellGeom::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellGeom*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderCellGeom::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellGeom*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderCellGeom::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCellGeom*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4CylinderCellGeom::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderCellGeom.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderCellGeom::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderCellGeom::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4CylinderCellGeom(void *p) {
      return  p ? new(p) ::PHG4CylinderCellGeom : new ::PHG4CylinderCellGeom;
   }
   static void *newArray_PHG4CylinderCellGeom(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4CylinderCellGeom[nElements] : new ::PHG4CylinderCellGeom[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4CylinderCellGeom(void *p) {
      delete ((::PHG4CylinderCellGeom*)p);
   }
   static void deleteArray_PHG4CylinderCellGeom(void *p) {
      delete [] ((::PHG4CylinderCellGeom*)p);
   }
   static void destruct_PHG4CylinderCellGeom(void *p) {
      typedef ::PHG4CylinderCellGeom current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderCellGeom

namespace {
  void TriggerDictionaryInitialization_PHG4CylinderCellGeom_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4CylinderCellGeom_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4CylinderCellGeom;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4CylinderCellGeom_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERCELLGEOM_H
#define G4DETECTORS_PHG4CYLINDERCELLGEOM_H

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream
#include <string>
#include <utility>           // for pair

class PHG4CylinderCellGeom: public PHObject
{
 public:
  PHG4CylinderCellGeom();

  virtual ~PHG4CylinderCellGeom() {}

  void identify(std::ostream& os = std::cout) const;
  int get_layer() const {return layer;}
  double get_radius() const {return radius;}
  double get_thickness() const {return thickness;}
  int get_binning() const {return binning;}
  int get_zbins() const;
  int get_phibins() const;
  double get_zmin() const;
  double get_phistep() const;
  double get_phimin() const;
  double get_zstep() const;
  int get_etabins() const;
  double get_etastep() const;
  double get_etamin() const;

  virtual std::pair<double, double> get_zbounds(const int ibin) const;
  virtual std::pair<double, double> get_phibounds(const int ibin) const;
  virtual std::pair<double, double> get_etabounds(const int ibin) const;
  virtual double get_etacenter(const int ibin) const;
  virtual double get_zcenter(const int ibin) const;
  virtual double get_phicenter(const int ibin) const;

  virtual int get_etabin(const double eta) const;
  virtual int get_zbin(const double z) const;
  virtual int get_phibin(const double phi) const;

   void set_layer(const int i) {layer = i;}
   void set_binning(const int i) {binning = i;}
   void set_radius(const double r) {radius = r;}
   void set_thickness(const double t) {thickness = t;}
   void set_zbins(const int i);
   void set_zmin(const double z);
   void set_zstep(const double z);
   void set_phibins(const int i);
   void set_phistep(const double phi);
   void set_phimin(const double phi);
   void set_etabins(const int i);
   void set_etamin(const double z);
   void set_etastep(const double z);
  
 protected:
  void check_binning_method(const int i) const;
  void check_binning_method_eta(const std::string & src = "") const;
  void check_binning_method_phi(const std::string & src = "") const;
  std::string methodname(const int i) const;
  int layer;
  int binning;
  double radius;
  int nzbins;
  double zmin;
  double zstep;
  int nphibins;
  double phimin;
  double phistep;
  double thickness;

  ClassDef(PHG4CylinderCellGeom,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4CylinderCellGeom", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4CylinderCellGeom_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4CylinderCellGeom_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4CylinderCellGeom_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4CylinderCellGeom_Dict() {
  TriggerDictionaryInitialization_PHG4CylinderCellGeom_Dict_Impl();
}
