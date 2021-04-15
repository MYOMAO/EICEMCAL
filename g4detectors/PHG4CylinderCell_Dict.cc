// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4CylinderCell_Dict
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
#include "PHG4CylinderCell.h"

// Header files passed via #pragma extra_include

namespace PHG4CylinderCellDefs {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *PHG4CylinderCellDefs_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("PHG4CylinderCellDefs", 0 /*version*/, "PHG4CylinderCellDefs.h", 6,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &PHG4CylinderCellDefs_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *PHG4CylinderCellDefs_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace ROOT {
   static void delete_PHG4CylinderCell(void *p);
   static void deleteArray_PHG4CylinderCell(void *p);
   static void destruct_PHG4CylinderCell(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4CylinderCell*)
   {
      ::PHG4CylinderCell *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4CylinderCell >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4CylinderCell", ::PHG4CylinderCell::Class_Version(), "", 22,
                  typeid(::PHG4CylinderCell), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4CylinderCell::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4CylinderCell) );
      instance.SetDelete(&delete_PHG4CylinderCell);
      instance.SetDeleteArray(&deleteArray_PHG4CylinderCell);
      instance.SetDestructor(&destruct_PHG4CylinderCell);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4CylinderCell*)
   {
      return GenerateInitInstanceLocal((::PHG4CylinderCell*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4CylinderCell*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4CylinderCell::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4CylinderCell::Class_Name()
{
   return "PHG4CylinderCell";
}

//______________________________________________________________________________
const char *PHG4CylinderCell::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCell*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4CylinderCell::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCell*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4CylinderCell::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCell*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4CylinderCell::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4CylinderCell*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4CylinderCell::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4CylinderCell.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4CylinderCell::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4CylinderCell::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PHG4CylinderCell(void *p) {
      delete ((::PHG4CylinderCell*)p);
   }
   static void deleteArray_PHG4CylinderCell(void *p) {
      delete [] ((::PHG4CylinderCell*)p);
   }
   static void destruct_PHG4CylinderCell(void *p) {
      typedef ::PHG4CylinderCell current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4CylinderCell

namespace {
  void TriggerDictionaryInitialization_PHG4CylinderCell_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4CylinderCell_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4CylinderCell;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4CylinderCell_Dict dictionary payload"

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  #define CGAL_DISABLE_ROUNDING_MATH_CHECK 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERCELL_H
#define G4DETECTORS_PHG4CYLINDERCELL_H

#include "PHG4Cell.h"

#include "PHG4CylinderCellDefs.h"

#include <g4main/PHG4HitDefs.h>

#include <phool/PHObject.h>

#include <cmath>
#include <map>

class PHG4CylinderCell : public PHG4Cell
{
 public:
/*
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
*/
  
  virtual ~PHG4CylinderCell(){}

  virtual void identify(std::ostream& os = std::cout) const {
    os << "PHG4CylinderCell base class" << std::endl;
  }
  
  
  
  virtual void set_cell_id(const PHG4CylinderCellDefs::keytype id) {return;}
  virtual void set_layer(const unsigned int i) {return;}

  virtual unsigned int get_layer() const {return 0xFFFFFFFF;}
  virtual PHG4CylinderCellDefs::keytype get_cell_id() const {return 0xFFFFFFFF;}
  virtual int get_binz() const {return -1;}
  virtual int get_binphi() const {return -1;}
  virtual int get_bineta() const {return -1;}

  virtual void set_etabin(const int i) {return;}
  virtual void set_light_yield(float lightYield)  {   return;  }

  virtual void set_fiber_ID(int fiberId) {   return;  }
  virtual int get_fiber_ID() const {return -1;}

  virtual void set_sensor_index(const std::string &si) {return;}
  virtual std::string get_sensor_index() const {return "";}

  virtual void set_ladder_phi_index(const int i) {return;}
  virtual int get_ladder_phi_index() const {return -9999;}

  virtual void set_ladder_z_index(const int i) {return;}
  virtual int get_ladder_z_index() const {return -9999;}

  virtual int get_j_index() const {return -9999;}
  virtual void set_j_index(const int i) {return;}
  virtual int get_k_index() const {return -9999;}
  virtual void set_k_index(const int i) {return;}
  virtual int get_l_index() const {return -9999;}
  virtual void set_l_index(const int i) {return;}

   protected:

  PHG4CylinderCell() {}
  ClassDef(PHG4CylinderCell,2)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4CylinderCell", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4CylinderCell_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4CylinderCell_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4CylinderCell_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4CylinderCell_Dict() {
  TriggerDictionaryInitialization_PHG4CylinderCell_Dict_Impl();
}
