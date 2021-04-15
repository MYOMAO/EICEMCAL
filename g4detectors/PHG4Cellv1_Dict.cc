// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4Cellv1_Dict
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
#include "PHG4Cellv1.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4Cellv1(void *p = 0);
   static void *newArray_PHG4Cellv1(Long_t size, void *p);
   static void delete_PHG4Cellv1(void *p);
   static void deleteArray_PHG4Cellv1(void *p);
   static void destruct_PHG4Cellv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4Cellv1*)
   {
      ::PHG4Cellv1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4Cellv1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4Cellv1", ::PHG4Cellv1::Class_Version(), "", 23,
                  typeid(::PHG4Cellv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4Cellv1::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4Cellv1) );
      instance.SetNew(&new_PHG4Cellv1);
      instance.SetNewArray(&newArray_PHG4Cellv1);
      instance.SetDelete(&delete_PHG4Cellv1);
      instance.SetDeleteArray(&deleteArray_PHG4Cellv1);
      instance.SetDestructor(&destruct_PHG4Cellv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4Cellv1*)
   {
      return GenerateInitInstanceLocal((::PHG4Cellv1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4Cellv1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4Cellv1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4Cellv1::Class_Name()
{
   return "PHG4Cellv1";
}

//______________________________________________________________________________
const char *PHG4Cellv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Cellv1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4Cellv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Cellv1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4Cellv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Cellv1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4Cellv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Cellv1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4Cellv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4Cellv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4Cellv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4Cellv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4Cellv1(void *p) {
      return  p ? new(p) ::PHG4Cellv1 : new ::PHG4Cellv1;
   }
   static void *newArray_PHG4Cellv1(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4Cellv1[nElements] : new ::PHG4Cellv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4Cellv1(void *p) {
      delete ((::PHG4Cellv1*)p);
   }
   static void deleteArray_PHG4Cellv1(void *p) {
      delete [] ((::PHG4Cellv1*)p);
   }
   static void destruct_PHG4Cellv1(void *p) {
      typedef ::PHG4Cellv1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4Cellv1

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

namespace ROOT {
   static TClass *maplEintcOfloatgR_Dictionary();
   static void maplEintcOfloatgR_TClassManip(TClass*);
   static void *new_maplEintcOfloatgR(void *p = 0);
   static void *newArray_maplEintcOfloatgR(Long_t size, void *p);
   static void delete_maplEintcOfloatgR(void *p);
   static void deleteArray_maplEintcOfloatgR(void *p);
   static void destruct_maplEintcOfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,float>*)
   {
      map<int,float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,float>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,float>", -2, "map", 100,
                  typeid(map<int,float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,float>) );
      instance.SetNew(&new_maplEintcOfloatgR);
      instance.SetNewArray(&newArray_maplEintcOfloatgR);
      instance.SetDelete(&delete_maplEintcOfloatgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOfloatgR);
      instance.SetDestructor(&destruct_maplEintcOfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,float> >()));

      ::ROOT::AddClassAlternate("map<int,float>","std::map<int, float, std::less<int>, std::allocator<std::pair<int const, float> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,float>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,float>*)0x0)->GetClass();
      maplEintcOfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,float> : new map<int,float>;
   }
   static void *newArray_maplEintcOfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,float>[nElements] : new map<int,float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOfloatgR(void *p) {
      delete ((map<int,float>*)p);
   }
   static void deleteArray_maplEintcOfloatgR(void *p) {
      delete [] ((map<int,float>*)p);
   }
   static void destruct_maplEintcOfloatgR(void *p) {
      typedef map<int,float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,float>

namespace ROOT {
   static TClass *maplEULong64_tcOfloatgR_Dictionary();
   static void maplEULong64_tcOfloatgR_TClassManip(TClass*);
   static void *new_maplEULong64_tcOfloatgR(void *p = 0);
   static void *newArray_maplEULong64_tcOfloatgR(Long_t size, void *p);
   static void delete_maplEULong64_tcOfloatgR(void *p);
   static void deleteArray_maplEULong64_tcOfloatgR(void *p);
   static void destruct_maplEULong64_tcOfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<ULong64_t,float>*)
   {
      map<ULong64_t,float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<ULong64_t,float>));
      static ::ROOT::TGenericClassInfo 
         instance("map<ULong64_t,float>", -2, "map", 100,
                  typeid(map<ULong64_t,float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEULong64_tcOfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(map<ULong64_t,float>) );
      instance.SetNew(&new_maplEULong64_tcOfloatgR);
      instance.SetNewArray(&newArray_maplEULong64_tcOfloatgR);
      instance.SetDelete(&delete_maplEULong64_tcOfloatgR);
      instance.SetDeleteArray(&deleteArray_maplEULong64_tcOfloatgR);
      instance.SetDestructor(&destruct_maplEULong64_tcOfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<ULong64_t,float> >()));

      ::ROOT::AddClassAlternate("map<ULong64_t,float>","std::map<unsigned long long, float, std::less<unsigned long long>, std::allocator<std::pair<unsigned long long const, float> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<ULong64_t,float>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEULong64_tcOfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<ULong64_t,float>*)0x0)->GetClass();
      maplEULong64_tcOfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEULong64_tcOfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEULong64_tcOfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<ULong64_t,float> : new map<ULong64_t,float>;
   }
   static void *newArray_maplEULong64_tcOfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<ULong64_t,float>[nElements] : new map<ULong64_t,float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEULong64_tcOfloatgR(void *p) {
      delete ((map<ULong64_t,float>*)p);
   }
   static void deleteArray_maplEULong64_tcOfloatgR(void *p) {
      delete [] ((map<ULong64_t,float>*)p);
   }
   static void destruct_maplEULong64_tcOfloatgR(void *p) {
      typedef map<ULong64_t,float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<ULong64_t,float>

namespace {
  void TriggerDictionaryInitialization_PHG4Cellv1_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4Cellv1_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4Cellv1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4Cellv1_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CELLV1_H
#define G4DETECTORS_PHG4CELLV1_H

#include "PHG4Cell.h"
#include "PHG4CellDefs.h"

#include <g4main/PHG4HitDefs.h>  // for keytype

#include <cstdint>
#include <iostream>
#include <map>
#include <type_traits>           // for __decay_and_strip<>::__type
#include <utility>               // for make_pair


class PHG4Cellv1: public PHG4Cell
{
 public:
  PHG4Cellv1();
  PHG4Cellv1(const PHG4CellDefs::keytype g4cellid);
  virtual ~PHG4Cellv1();

  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset();

  void set_cellid(const PHG4CellDefs::keytype i) {cellid = i;}

  PHG4CellDefs::keytype get_cellid() const {return cellid;}
  bool has_binning(const PHG4CellDefs::CellBinning binning) const;
  short int get_detid() const;

  void add_edep(const PHG4HitDefs::keytype g4hitid, const float edep);
  void add_shower_edep(const int g4showerid, const float edep);

  EdepConstRange get_g4hits() {
    return std::make_pair(hitedeps.begin(), hitedeps.end());
  }
  
  ShowerEdepConstRange get_g4showers() {
    return std::make_pair(showeredeps.begin(),showeredeps.end());
  } 


  void add_edep(const float f) {add_property(prop_edep,f);}
  double get_edep() const {return get_property_float(prop_edep);}

  void add_eion(const float f) {add_property(prop_eion,f);}
  double get_eion() const {return get_property_float(prop_eion);}

  void add_light_yield(const float f) {add_property(prop_light_yield,f);}
  float get_light_yield() const {return get_property_float(prop_light_yield);}

  void set_chip_index(const int i) {set_property(prop_chip_index,i);}
  int get_chip_index() const {return get_property_int(prop_chip_index);}

  void set_half_stave_index(const int i) {set_property(prop_half_stave_index,i);}
  int get_half_stave_index() const {return get_property_int(prop_half_stave_index);}

  void set_ladder_phi_index(const int i) {set_property(prop_ladder_phi_index,i);}
  int get_ladder_phi_index() const {return get_property_int(prop_ladder_phi_index);}

  void set_ladder_z_index(const int i) {set_property(prop_ladder_z_index,i);}
  int get_ladder_z_index() const {return get_property_int(prop_ladder_z_index);}

  void set_module_index(const int i) {set_property(prop_module_index,i);}
  int get_module_index() const {return get_property_int(prop_module_index);}

  void set_phibin(const int i) {set_property(prop_phibin,i);}
  int get_phibin() const {return get_property_int(prop_phibin);}

  void set_pixel_index(const int i) {set_property(prop_pixel_index,i);}
  int get_pixel_index() const {return get_property_int(prop_pixel_index);}

  void set_stave_index(const int i) {set_property(prop_stave_index,i);}
  int get_stave_index() const {return get_property_int(prop_stave_index);}

  void set_zbin(const int i) {set_property(prop_zbin,i);}
  int get_zbin() const {return get_property_int(prop_zbin);}

  void print() const;

  bool  has_property(const PROPERTY prop_id) const;
  float get_property_float(const PROPERTY prop_id) const;
  int   get_property_int(const PROPERTY prop_id) const;
  unsigned int   get_property_uint(const PROPERTY prop_id) const;
  void  add_property(const PROPERTY prop_id, const float value);
  void  add_property(const PROPERTY prop_id, const int value);
  void  add_property(const PROPERTY prop_id, const unsigned int value);
  void  set_property(const PROPERTY prop_id, const float value);
  void  set_property(const PROPERTY prop_id, const int value);
  void  set_property(const PROPERTY prop_id, const unsigned int value);


 protected:
  unsigned int get_property_nocheck(const PROPERTY prop_id) const;
  void set_property_nocheck(const PROPERTY prop_id,const unsigned int ui) {prop_map[prop_id]=ui;}

  PHG4CellDefs::keytype cellid;
  EdepMap hitedeps;
  ShowerEdepMap showeredeps;

  //! storage types for additional property
  typedef uint8_t prop_id_t;
  typedef uint32_t prop_storage_t;
  typedef std::map<prop_id_t, prop_storage_t> prop_map_t;

  //! convert between 32bit inputs and storage type prop_storage_t
  union u_property{
    float fdata;
    int32_t idata;
    uint32_t uidata;

    u_property(int32_t in): idata(in) {}
    u_property(uint32_t in): uidata(in) {}
    u_property(float in): fdata(in) {}
    u_property(): uidata(0) {}

    prop_storage_t get_storage() const {return uidata;}
  };

  //! container for additional property
  prop_map_t prop_map;


  ClassDef(PHG4Cellv1,3)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4Cellv1", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4Cellv1_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4Cellv1_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4Cellv1_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4Cellv1_Dict() {
  TriggerDictionaryInitialization_PHG4Cellv1_Dict_Impl();
}
