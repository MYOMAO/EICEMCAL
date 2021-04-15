// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4Cell_Dict
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
#include "PHG4Cell.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_PHG4Cell(void *p);
   static void deleteArray_PHG4Cell(void *p);
   static void destruct_PHG4Cell(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4Cell*)
   {
      ::PHG4Cell *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4Cell >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4Cell", ::PHG4Cell::Class_Version(), "", 24,
                  typeid(::PHG4Cell), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4Cell::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4Cell) );
      instance.SetDelete(&delete_PHG4Cell);
      instance.SetDeleteArray(&deleteArray_PHG4Cell);
      instance.SetDestructor(&destruct_PHG4Cell);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4Cell*)
   {
      return GenerateInitInstanceLocal((::PHG4Cell*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4Cell*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4Cell::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4Cell::Class_Name()
{
   return "PHG4Cell";
}

//______________________________________________________________________________
const char *PHG4Cell::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Cell*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4Cell::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Cell*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4Cell::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Cell*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4Cell::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Cell*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4Cell::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4Cell.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4Cell::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4Cell::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PHG4Cell(void *p) {
      delete ((::PHG4Cell*)p);
   }
   static void deleteArray_PHG4Cell(void *p) {
      delete [] ((::PHG4Cell*)p);
   }
   static void destruct_PHG4Cell(void *p) {
      typedef ::PHG4Cell current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4Cell

namespace {
  void TriggerDictionaryInitialization_PHG4Cell_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4Cell_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4Cell;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4Cell_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CELL_H
#define G4DETECTORS_PHG4CELL_H

#include "PHG4CellDefs.h"

#include <g4main/PHG4HitDefs.h>  // for keytype

#include <phool/PHObject.h>

#include <cmath>
#include <climits>
#include <iostream>              // for ostream, cout, operator<<, endl, bas...
#include <map>
#include <string>                // for string
#include <utility>               // for pair, make_pair

class PHG4Cell: public PHObject
{
 public:
  typedef std::map<PHG4HitDefs::keytype, float> EdepMap;
  typedef EdepMap::iterator EdepIterator;
  typedef EdepMap::const_iterator EdepConstIterator;
  typedef std::pair<EdepIterator, EdepIterator> EdepRange;
  typedef std::pair<EdepConstIterator, EdepConstIterator> EdepConstRange;

  typedef std::map<int, float> ShowerEdepMap;
  typedef ShowerEdepMap::iterator ShowerEdepIterator;
  typedef ShowerEdepMap::const_iterator ShowerEdepConstIterator;
  typedef std::pair<ShowerEdepIterator, ShowerEdepIterator> ShowerEdepRange;
  typedef std::pair<ShowerEdepConstIterator, ShowerEdepConstIterator> ShowerEdepConstRange;

  virtual ~PHG4Cell() {}

  virtual void identify(std::ostream& os = std::cout) const;
  virtual void CopyFrom(const PHObject *phobj);
  friend std::ostream &operator<<(std::ostream & stream, const PHG4Cell * cell);
  virtual void Reset();

  // all methods connected to the cell id (encoding/decoding
  virtual void set_cellid(const PHG4CellDefs::keytype i) {return;}

  virtual PHG4CellDefs::keytype get_cellid() const {return ~0x0;}
  virtual bool has_binning(const PHG4CellDefs::CellBinning) const {return false;}

  // this adds hits to the g4 hit list map
  virtual void add_edep(const PHG4HitDefs::keytype g4hitid, const float edep) {return;}
  virtual void add_edep(const PHG4HitDefs::keytype g4hitid, const float edep, const float light_yield) {return;}
  virtual void add_edep(const PHG4HitDefs::keytype g4hitid, const int tbin, const float edep) {return;}
  // this adds showers to the shower map
  virtual void add_shower_edep(const int g4showerid, const float edep) {return;}

  virtual EdepConstRange get_g4hits();

  virtual ShowerEdepConstRange get_g4showers();

  virtual short int get_detid() const {return -1;}
  // for backward compatibility, layers and detector ids are identical
  short int get_layer() const {return get_detid();}

  virtual void add_edep(const float f) {return;}
  virtual double get_edep() const {return NAN;}

  virtual void add_eion(const float f) {return;}
  virtual double get_eion() const {return NAN;}

  virtual void add_light_yield(const float lightYield){return;}
  virtual float get_light_yield() const {return NAN;}

  // get/set methodes - PLEASE add those ALPHABETICALLY

  virtual void set_chip_index(const int i) {return;}
  virtual int get_chip_index() const {return ~0x0;}

  virtual void set_half_stave_index(const int i) {return;}
  virtual int get_half_stave_index() const {return ~0x0;}

  virtual void set_ladder_phi_index(const int i) {return;}
  virtual int get_ladder_phi_index() const {return ~0x0;}

  virtual void set_ladder_z_index(const int i) {return;}
  virtual int get_ladder_z_index() const {return ~0x0;}

  virtual void set_module_index(const int i) {return;}
  virtual int get_module_index() const {return ~0x0;}

  virtual void set_phibin(const int i) {return;}
  virtual int get_phibin() const {return ~0x0;}

  virtual void set_pixel_index(const int i) {return;}
  virtual int get_pixel_index() const {return ~0x0;}

  virtual void set_stave_index(const int i) {return;}
  virtual int get_stave_index() const {return ~0x0;}

//  virtual tpctod* get_train_of_digits() {return 0;}

  virtual void set_zbin(const int i) {return;}
  virtual int get_zbin() const {return ~0x0;}


  virtual void print() const {std::cout<<"virtual PHG4Cell"<<std::endl;}

  //! Procedure to add a new PROPERTY tag:
  //! 1.add new tag below with unique value,
  //! 2.add a short name to PHG4Cell::get_property_info
  enum PROPERTY 
  {//
    // first various coordinates 1-20
    //! Maps coordinates
    prop_stave_index = 1,
    prop_half_stave_index = 2,
    prop_module_index = 3,
    prop_chip_index = 4,
    prop_pixel_index = 5,
    prop_phibin = 6,
    prop_zbin = 7,
    prop_ladder_z_index = 8,
    prop_ladder_phi_index = 9,
    //-- summed energy:  - 20-30  --
    //! deposited energy
    prop_edep = 21,
    //! ionizing energy loss
    prop_eion = 22,

    //! for scintillation detectors, the amount of light produced
    prop_light_yield = 23,

    //! max limit in order to fit into 8 bit unsigned number
    prop_MAX_NUMBER = UCHAR_MAX
  };

  enum PROPERTY_TYPE 
  {//
    type_int = 1,
    type_uint = 2,
    type_float = 3,
    type_unknown = -1
  };

  virtual bool  has_property(const PROPERTY prop_id) const {return false;}
  virtual float get_property_float(const PROPERTY prop_id) const {return NAN;}
  virtual int   get_property_int(const PROPERTY prop_id) const {return INT_MIN;}
  virtual unsigned int   get_property_uint(const PROPERTY prop_id) const {return UINT_MAX;}
  virtual void  set_property(const PROPERTY prop_id, const float value) {return;}
  virtual void  set_property(const PROPERTY prop_id, const int value) {return;}
  virtual void  set_property(const PROPERTY prop_id, const unsigned int value) {return;}
  static std::pair<const std::string,PROPERTY_TYPE> get_property_info(PROPERTY prop_id);
  static bool check_property(const PROPERTY prop_id, const PROPERTY_TYPE prop_type);
  static std::string get_property_type(const PROPERTY_TYPE prop_type);

 protected:
  PHG4Cell() {}
  virtual unsigned int get_property_nocheck(const PROPERTY prop_id) const {return UINT_MAX;}
  virtual void set_property_nocheck(const PROPERTY prop_id,const unsigned int) {return;}
  ClassDef(PHG4Cell,2)
};


#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4Cell", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4Cell_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4Cell_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4Cell_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4Cell_Dict() {
  TriggerDictionaryInitialization_PHG4Cell_Dict_Impl();
}
