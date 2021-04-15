// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4Hit_Dict
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
#include "PHG4Hit.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4Hit(void *p = 0);
   static void *newArray_PHG4Hit(Long_t size, void *p);
   static void delete_PHG4Hit(void *p);
   static void deleteArray_PHG4Hit(void *p);
   static void destruct_PHG4Hit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4Hit*)
   {
      ::PHG4Hit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4Hit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4Hit", ::PHG4Hit::Class_Version(), "", 23,
                  typeid(::PHG4Hit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4Hit::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4Hit) );
      instance.SetNew(&new_PHG4Hit);
      instance.SetNewArray(&newArray_PHG4Hit);
      instance.SetDelete(&delete_PHG4Hit);
      instance.SetDeleteArray(&deleteArray_PHG4Hit);
      instance.SetDestructor(&destruct_PHG4Hit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4Hit*)
   {
      return GenerateInitInstanceLocal((::PHG4Hit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4Hit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4Hit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4Hit::Class_Name()
{
   return "PHG4Hit";
}

//______________________________________________________________________________
const char *PHG4Hit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Hit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4Hit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Hit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4Hit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Hit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4Hit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Hit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4Hit::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4Hit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4Hit::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4Hit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4Hit(void *p) {
      return  p ? new(p) ::PHG4Hit : new ::PHG4Hit;
   }
   static void *newArray_PHG4Hit(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4Hit[nElements] : new ::PHG4Hit[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4Hit(void *p) {
      delete ((::PHG4Hit*)p);
   }
   static void deleteArray_PHG4Hit(void *p) {
      delete [] ((::PHG4Hit*)p);
   }
   static void destruct_PHG4Hit(void *p) {
      typedef ::PHG4Hit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4Hit

namespace {
  void TriggerDictionaryInitialization_PHG4Hit_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4Hit_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4Hit;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4Hit_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4HIT_H
#define G4MAIN_PHG4HIT_H

#include "PHG4HitDefs.h"

#include <phool/PHObject.h>

#include <cmath>
#include <climits>
#include <iostream> 
#include <string>
#include <utility>

#include "TFile.h"
#include "TTree.h"



class PHG4Hit: public PHObject
{
 public:
  PHG4Hit() {}
  virtual ~PHG4Hit() {}

  virtual void identify(std::ostream& os = std::cout) const;
  virtual void CopyFrom(const PHObject *phobj);
  friend std::ostream &operator<<(std::ostream & stream, const PHG4Hit * hit);
  virtual void Reset();

  // The indices here represent the entry and exit points of the particle
  virtual float get_x(const int i) const {return NAN;}
  virtual float get_y(const int i) const {return NAN;}
  virtual float get_z(const int i) const {return NAN;}
  virtual float get_px(const int i) const {return NAN;}
  virtual float get_py(const int i) const {return NAN;}
  virtual float get_pz(const int i) const {return NAN;}
  virtual float get_local_x(const int i) const {return NAN;}
  virtual float get_local_y(const int i) const {return NAN;}
  virtual float get_local_z(const int i) const {return NAN;}
  virtual float get_t(const int i) const {return NAN;}
  virtual float get_edep() const {return NAN;}
  virtual float get_eion() const {return NAN;}
  virtual float get_light_yield() const {return NAN;}
  virtual float get_path_length() const {return NAN;}
  virtual unsigned int get_layer() const {return UINT_MAX;}
  virtual PHG4HitDefs::keytype get_hit_id() const {return ULONG_LONG_MAX;}
  virtual int get_detid() const {return INT_MIN;}
  virtual int get_shower_id() const {return INT_MIN;}
  virtual int get_scint_id() const {return INT_MIN;}
  virtual int get_row() const {return INT_MIN;}
  virtual int get_trkid() const {return INT_MIN;}
  virtual int get_strip_z_index() const {return INT_MIN;}
  virtual int get_strip_y_index() const {return INT_MIN;}
  virtual int get_ladder_z_index() const {return INT_MIN;}
  virtual int get_ladder_phi_index() const {return INT_MIN;}
  virtual int get_index_i() const {return INT_MIN;}
  virtual int get_index_j() const {return INT_MIN;}
  virtual int get_index_k() const {return INT_MIN;}
  virtual int get_index_l() const {return INT_MIN;}
  virtual int get_hit_type() const {return INT_MIN;}

  virtual void set_x(const int i, const float f) {return;}
  virtual void set_y(const int i, const float f) {return;}
  virtual void set_z(const int i, const float f) {return;}
  virtual void set_px(const int i, const float f) {return;}
  virtual void set_py(const int i, const float f) {return;}
  virtual void set_pz(const int i, const float f) {return;}
  virtual void set_local_x(const int i, const float f) {return;}
  virtual void set_local_y(const int i, const float f) {return;}
  virtual void set_local_z(const int i, const float f) {return;}
  virtual void set_t(const int i, const float f) {return;}
  virtual void set_edep(const float f) {return;}
  virtual void set_eion(const float f) {return;}
  virtual void set_light_yield(const float lightYield){return;}
  virtual void set_path_length(const float pathLength){return;}
  virtual void set_layer(const unsigned int i) {return;}
  virtual void set_hit_id(const PHG4HitDefs::keytype i) {return;}
  virtual void set_shower_id(const int i) {return;}
  virtual void set_scint_id(const int i) {return;}
  virtual void set_row(const int i) {return;}
  virtual void set_trkid(const int i) {return;}
  virtual void set_strip_z_index(const int i) {return;}
  virtual void set_strip_y_index(const int i) {return;}
  virtual void set_ladder_z_index(const int i) {return;}
  virtual void set_ladder_phi_index(const int i) {return;}
  virtual void set_index_i(const int i) {return;}
  virtual void set_index_j(const int i) {return;}
  virtual void set_index_k(const int i) {return;}
  virtual void set_index_l(const int i) {return;}
  virtual void set_hit_type(const int i) {return;}

  virtual float get_avg_x() const;
  virtual float get_avg_y() const;
  virtual float get_avg_z() const;
  virtual float get_avg_t() const;

  virtual void print() const {std::cout<<"PHG4Hit base class - print() not implemented"<<std::endl;}


  //! Procedure to add a new PROPERTY tag:
  //! 1.add new tag below with unique value,
  //! 2.add a short name to PHG4Hit::get_property_info
  enum PROPERTY 
  {//

    //-- hit properties: 1 - 10  --
    //! ionizing energy loss
    prop_eion = 1,

    //! for scintillation detectors, the amount of light produced
    prop_light_yield = 2,

    //-- track properties: 10 - 20  --

    //! momentum
    prop_px_0 = 10,
    prop_px_1 = 11,
    prop_py_0 = 12,
    prop_py_1 = 13,
    prop_pz_0 = 14,
    prop_pz_1 = 15,

    //! pathlength
    prop_path_length = 16,

    //! local coordinate
    prop_local_x_0 = 20,
    prop_local_x_1 = 21,
    prop_local_y_0 = 22,
    prop_local_y_1 = 23,
    prop_local_z_0 = 24,
    prop_local_z_1 = 25,

    //-- detector specific IDs: 100+ --

    //! layer ID
    prop_layer = 101,
    //! scintillator ID
    prop_scint_id = 102,
    //! row (mother volume or steel plate id)
    prop_row = 103,

    //! SVX stuff
    prop_strip_z_index = 110,
    prop_strip_y_index = 111,
    prop_ladder_z_index = 112,
    prop_ladder_phi_index = 113,

    // MAPS stuff 
    prop_stave_index = 114,
    prop_half_stave_index = 115,
    prop_module_index = 116,
    prop_chip_index = 117,

    prop_local_pos_x_0 = 118,
    prop_local_pos_y_0 = 119,
    prop_local_pos_z_0 = 120,

    //! generic indexes
    prop_index_i = 121,
    prop_index_j = 122,
    prop_index_k = 123,
    prop_index_l = 124,

    //! hit type
    prop_hit_type = 125,

    prop_local_pos_x_1 = 128,
    prop_local_pos_y_1 = 126,
    prop_local_pos_z_1 = 127,

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
  virtual unsigned int get_property_nocheck(const PROPERTY prop_id) const {return UINT_MAX;}
  virtual void set_property_nocheck(const PROPERTY prop_id,const unsigned int) {return;}
  ClassDef(PHG4Hit,1)
};


inline float PHG4Hit::get_avg_x() const { return 0.5*(get_x(0)+get_x(1)); }
inline float PHG4Hit::get_avg_y() const { return 0.5*(get_y(0)+get_y(1)); }
inline float PHG4Hit::get_avg_z() const { return 0.5*(get_z(0)+get_z(1)); }
inline float PHG4Hit::get_avg_t() const { return 0.5*(get_t(0)+get_t(1)); }


inline  TFile * fout; 
inline TTree * ShowerInfo;
inline int evtid;
inline  float X1;
inline float X2;
inline float Y1;
inline  float Y2;
inline  float Z1;
inline  float Z2;
inline  float T1;
inline  float T2;
inline  float BHBLedep;
inline  float BHFLedep;
inline  float	FEMCALedep;
inline  int hitid;
inline  int trkid;
inline  int showerid;
inline  int layer;	
inline float ETotalDep;
inline float FEMCALLYDep;
inline float FEMCALEIonDep;
inline float AbsorberDep;
inline float ERECO;
//inline int NClusters;

//Photon Info

inline float GPX1;
inline float GPY1;
inline float GPZ1;
inline float GE1;


inline float GPX2;
inline float GPY2;
inline float GPZ2;
inline float GE2;
inline int NumOfE;

inline float EOut;

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4Hit", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4Hit_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4Hit_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4Hit_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4Hit_Dict() {
  TriggerDictionaryInitialization_PHG4Hit_Dict_Impl();
}
