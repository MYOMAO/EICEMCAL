// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PHG4Showerv1_Dict
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
#include "PHG4Showerv1.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PHG4Showerv1(void *p = 0);
   static void *newArray_PHG4Showerv1(Long_t size, void *p);
   static void delete_PHG4Showerv1(void *p);
   static void deleteArray_PHG4Showerv1(void *p);
   static void destruct_PHG4Showerv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PHG4Showerv1*)
   {
      ::PHG4Showerv1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PHG4Showerv1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PHG4Showerv1", ::PHG4Showerv1::Class_Version(), "", 17,
                  typeid(::PHG4Showerv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PHG4Showerv1::Dictionary, isa_proxy, 4,
                  sizeof(::PHG4Showerv1) );
      instance.SetNew(&new_PHG4Showerv1);
      instance.SetNewArray(&newArray_PHG4Showerv1);
      instance.SetDelete(&delete_PHG4Showerv1);
      instance.SetDeleteArray(&deleteArray_PHG4Showerv1);
      instance.SetDestructor(&destruct_PHG4Showerv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PHG4Showerv1*)
   {
      return GenerateInitInstanceLocal((::PHG4Showerv1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PHG4Showerv1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PHG4Showerv1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PHG4Showerv1::Class_Name()
{
   return "PHG4Showerv1";
}

//______________________________________________________________________________
const char *PHG4Showerv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Showerv1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PHG4Showerv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PHG4Showerv1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PHG4Showerv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Showerv1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PHG4Showerv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PHG4Showerv1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PHG4Showerv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class PHG4Showerv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PHG4Showerv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(PHG4Showerv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PHG4Showerv1(void *p) {
      return  p ? new(p) ::PHG4Showerv1 : new ::PHG4Showerv1;
   }
   static void *newArray_PHG4Showerv1(Long_t nElements, void *p) {
      return p ? new(p) ::PHG4Showerv1[nElements] : new ::PHG4Showerv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_PHG4Showerv1(void *p) {
      delete ((::PHG4Showerv1*)p);
   }
   static void deleteArray_PHG4Showerv1(void *p) {
      delete [] ((::PHG4Showerv1*)p);
   }
   static void destruct_PHG4Showerv1(void *p) {
      typedef ::PHG4Showerv1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PHG4Showerv1

namespace ROOT {
   static TClass *setlEintgR_Dictionary();
   static void setlEintgR_TClassManip(TClass*);
   static void *new_setlEintgR(void *p = 0);
   static void *newArray_setlEintgR(Long_t size, void *p);
   static void delete_setlEintgR(void *p);
   static void deleteArray_setlEintgR(void *p);
   static void destruct_setlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const set<int>*)
   {
      set<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(set<int>));
      static ::ROOT::TGenericClassInfo 
         instance("set<int>", -2, "set", 94,
                  typeid(set<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &setlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(set<int>) );
      instance.SetNew(&new_setlEintgR);
      instance.SetNewArray(&newArray_setlEintgR);
      instance.SetDelete(&delete_setlEintgR);
      instance.SetDeleteArray(&deleteArray_setlEintgR);
      instance.SetDestructor(&destruct_setlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Insert< set<int> >()));

      ::ROOT::AddClassAlternate("set<int>","std::set<int, std::less<int>, std::allocator<int> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const set<int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *setlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const set<int>*)0x0)->GetClass();
      setlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void setlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_setlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<int> : new set<int>;
   }
   static void *newArray_setlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<int>[nElements] : new set<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_setlEintgR(void *p) {
      delete ((set<int>*)p);
   }
   static void deleteArray_setlEintgR(void *p) {
      delete [] ((set<int>*)p);
   }
   static void destruct_setlEintgR(void *p) {
      typedef set<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class set<int>

namespace ROOT {
   static TClass *setlEULong64_tgR_Dictionary();
   static void setlEULong64_tgR_TClassManip(TClass*);
   static void *new_setlEULong64_tgR(void *p = 0);
   static void *newArray_setlEULong64_tgR(Long_t size, void *p);
   static void delete_setlEULong64_tgR(void *p);
   static void deleteArray_setlEULong64_tgR(void *p);
   static void destruct_setlEULong64_tgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const set<ULong64_t>*)
   {
      set<ULong64_t> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(set<ULong64_t>));
      static ::ROOT::TGenericClassInfo 
         instance("set<ULong64_t>", -2, "set", 94,
                  typeid(set<ULong64_t>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &setlEULong64_tgR_Dictionary, isa_proxy, 0,
                  sizeof(set<ULong64_t>) );
      instance.SetNew(&new_setlEULong64_tgR);
      instance.SetNewArray(&newArray_setlEULong64_tgR);
      instance.SetDelete(&delete_setlEULong64_tgR);
      instance.SetDeleteArray(&deleteArray_setlEULong64_tgR);
      instance.SetDestructor(&destruct_setlEULong64_tgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Insert< set<ULong64_t> >()));

      ::ROOT::AddClassAlternate("set<ULong64_t>","std::set<unsigned long long, std::less<unsigned long long>, std::allocator<unsigned long long> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const set<ULong64_t>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *setlEULong64_tgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const set<ULong64_t>*)0x0)->GetClass();
      setlEULong64_tgR_TClassManip(theClass);
   return theClass;
   }

   static void setlEULong64_tgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_setlEULong64_tgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<ULong64_t> : new set<ULong64_t>;
   }
   static void *newArray_setlEULong64_tgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<ULong64_t>[nElements] : new set<ULong64_t>[nElements];
   }
   // Wrapper around operator delete
   static void delete_setlEULong64_tgR(void *p) {
      delete ((set<ULong64_t>*)p);
   }
   static void deleteArray_setlEULong64_tgR(void *p) {
      delete [] ((set<ULong64_t>*)p);
   }
   static void destruct_setlEULong64_tgR(void *p) {
      typedef set<ULong64_t> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class set<ULong64_t>

namespace ROOT {
   static TClass *maplEintcOunsignedsPintgR_Dictionary();
   static void maplEintcOunsignedsPintgR_TClassManip(TClass*);
   static void *new_maplEintcOunsignedsPintgR(void *p = 0);
   static void *newArray_maplEintcOunsignedsPintgR(Long_t size, void *p);
   static void delete_maplEintcOunsignedsPintgR(void *p);
   static void deleteArray_maplEintcOunsignedsPintgR(void *p);
   static void destruct_maplEintcOunsignedsPintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,unsigned int>*)
   {
      map<int,unsigned int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,unsigned int>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,unsigned int>", -2, "map", 100,
                  typeid(map<int,unsigned int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOunsignedsPintgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,unsigned int>) );
      instance.SetNew(&new_maplEintcOunsignedsPintgR);
      instance.SetNewArray(&newArray_maplEintcOunsignedsPintgR);
      instance.SetDelete(&delete_maplEintcOunsignedsPintgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOunsignedsPintgR);
      instance.SetDestructor(&destruct_maplEintcOunsignedsPintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,unsigned int> >()));

      ::ROOT::AddClassAlternate("map<int,unsigned int>","std::map<int, unsigned int, std::less<int>, std::allocator<std::pair<int const, unsigned int> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,unsigned int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOunsignedsPintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,unsigned int>*)0x0)->GetClass();
      maplEintcOunsignedsPintgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOunsignedsPintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOunsignedsPintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,unsigned int> : new map<int,unsigned int>;
   }
   static void *newArray_maplEintcOunsignedsPintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,unsigned int>[nElements] : new map<int,unsigned int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOunsignedsPintgR(void *p) {
      delete ((map<int,unsigned int>*)p);
   }
   static void deleteArray_maplEintcOunsignedsPintgR(void *p) {
      delete [] ((map<int,unsigned int>*)p);
   }
   static void destruct_maplEintcOunsignedsPintgR(void *p) {
      typedef map<int,unsigned int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,unsigned int>

namespace ROOT {
   static TClass *maplEintcOsetlEULong64_tgRsPgR_Dictionary();
   static void maplEintcOsetlEULong64_tgRsPgR_TClassManip(TClass*);
   static void *new_maplEintcOsetlEULong64_tgRsPgR(void *p = 0);
   static void *newArray_maplEintcOsetlEULong64_tgRsPgR(Long_t size, void *p);
   static void delete_maplEintcOsetlEULong64_tgRsPgR(void *p);
   static void deleteArray_maplEintcOsetlEULong64_tgRsPgR(void *p);
   static void destruct_maplEintcOsetlEULong64_tgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,set<ULong64_t> >*)
   {
      map<int,set<ULong64_t> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,set<ULong64_t> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,set<ULong64_t> >", -2, "map", 100,
                  typeid(map<int,set<ULong64_t> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOsetlEULong64_tgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,set<ULong64_t> >) );
      instance.SetNew(&new_maplEintcOsetlEULong64_tgRsPgR);
      instance.SetNewArray(&newArray_maplEintcOsetlEULong64_tgRsPgR);
      instance.SetDelete(&delete_maplEintcOsetlEULong64_tgRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOsetlEULong64_tgRsPgR);
      instance.SetDestructor(&destruct_maplEintcOsetlEULong64_tgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,set<ULong64_t> > >()));

      ::ROOT::AddClassAlternate("map<int,set<ULong64_t> >","std::map<int, std::set<unsigned long long, std::less<unsigned long long>, std::allocator<unsigned long long> >, std::less<int>, std::allocator<std::pair<int const, std::set<unsigned long long, std::less<unsigned long long>, std::allocator<unsigned long long> > > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,set<ULong64_t> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOsetlEULong64_tgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,set<ULong64_t> >*)0x0)->GetClass();
      maplEintcOsetlEULong64_tgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOsetlEULong64_tgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOsetlEULong64_tgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,set<ULong64_t> > : new map<int,set<ULong64_t> >;
   }
   static void *newArray_maplEintcOsetlEULong64_tgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,set<ULong64_t> >[nElements] : new map<int,set<ULong64_t> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOsetlEULong64_tgRsPgR(void *p) {
      delete ((map<int,set<ULong64_t> >*)p);
   }
   static void deleteArray_maplEintcOsetlEULong64_tgRsPgR(void *p) {
      delete [] ((map<int,set<ULong64_t> >*)p);
   }
   static void destruct_maplEintcOsetlEULong64_tgRsPgR(void *p) {
      typedef map<int,set<ULong64_t> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,set<ULong64_t> >

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

namespace {
  void TriggerDictionaryInitialization_PHG4Showerv1_Dict_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PHG4Showerv1_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class PHG4Showerv1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PHG4Showerv1_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4SHOWERV1_H
#define G4MAIN_PHG4SHOWERV1_H

#include "PHG4Shower.h"

#include "PHG4HitDefs.h"

#include <cstddef>       // for size_t
#include <iostream>
#include <map>
#include <set>

class PHG4Showerv1 : public PHG4Shower
{
 public:
  PHG4Showerv1();
  virtual ~PHG4Showerv1() {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const;
  PHG4Shower* CloneMe() const { return (new PHG4Showerv1(*this)); }
  void Reset() { *this = PHG4Showerv1(); }
  int isValid() const;

  // shower info

  int get_id() const { return _id; }
  void set_id(int id) { _id = id; }

  int get_parent_particle_id() const { return _parent_particle_id; }
  void set_parent_particle_id(int parent_particle_id) { _parent_particle_id = parent_particle_id; }

  int get_parent_shower_id() const { return _parent_shower_id; }
  void set_parent_shower_id(int parent_shower_id) { _parent_shower_id = parent_shower_id; }

  float get_x() const { return _pos[0]; }
  void set_x(float x) { _pos[0] = x; }

  float get_y() const { return _pos[1]; }
  void set_y(float y) { _pos[1] = y; }

  float get_z() const { return _pos[2]; }
  void set_z(float z) { _pos[2] = z; }

  float get_position(unsigned int coor) const { return _pos[coor]; }
  void set_position(unsigned int coor, float xi) { _pos[coor] = xi; }

  float get_covar(unsigned int i, unsigned int j) const;
  void set_covar(unsigned int i, unsigned int j, float entry);

  unsigned int get_nhits(int volume) const;
  void set_nhits(int volume, unsigned int nhits) { _nhits[volume] = nhits; }

  double get_edep() const;
  float get_edep(int volume) const;
  void set_edep(int volume, float edep) { _edep[volume] = edep; }

  double get_eion() const;
  float get_eion(int volume) const;
  void set_eion(int volume, float eion) { _eion[volume] = eion; }

  float get_light_yield(int volume) const;
  void set_light_yield(int volume, float light_yield) { _light_yield[volume] = light_yield; }

  float get_eh_ratio(int volume) const;
  void set_eh_ratio(int volume, float eh_ratio) { _eh_ratio[volume] = eh_ratio; }

  // container methods for ids
  bool empty_g4particle_id() const { return _g4particle_ids.empty(); }
  size_t size_g4particle_id() const { return _g4particle_ids.size(); }
  void add_g4particle_id(int id) { _g4particle_ids.insert(id); }
  PHG4Shower::ParticleIdIter begin_g4particle_id() { return _g4particle_ids.begin(); }
  PHG4Shower::ParticleIdConstIter begin_g4particle_id() const { return _g4particle_ids.begin(); }
  PHG4Shower::ParticleIdIter end_g4particle_id() { return _g4particle_ids.end(); }
  PHG4Shower::ParticleIdConstIter end_g4particle_id() const { return _g4particle_ids.end(); }
  size_t remove_g4particle_id(int id) { return _g4particle_ids.erase(id); }
  void clear_g4particle_id() { return _g4particle_ids.clear(); }

  bool empty_g4vertex_id() const { return _g4vertex_ids.empty(); }
  size_t size_g4vertex_id() const { return _g4vertex_ids.size(); }
  void add_g4vertex_id(int id) { _g4vertex_ids.insert(id); }
  PHG4Shower::VertexIdIter begin_g4vertex_id() { return _g4vertex_ids.begin(); }
  PHG4Shower::VertexIdConstIter begin_g4vertex_id() const { return _g4vertex_ids.begin(); }
  PHG4Shower::VertexIdIter end_g4vertex_id() { return _g4vertex_ids.end(); }
  PHG4Shower::VertexIdConstIter end_g4vertex_id() const { return _g4vertex_ids.end(); }
  size_t remove_g4vertex_id(int id) { return _g4vertex_ids.erase(id); }
  void clear_g4vertex_id() { return _g4vertex_ids.clear(); }

  bool empty_g4hit_id() const { return _g4hit_ids.empty(); }
  size_t size_g4hit_id() const { return _g4hit_ids.size(); }
  void add_g4hit_id(int volume, PHG4HitDefs::keytype id) { _g4hit_ids[volume].insert(id); }
  PHG4Shower::HitIdIter begin_g4hit_id() { return _g4hit_ids.begin(); }
  PHG4Shower::HitIdConstIter begin_g4hit_id() const { return _g4hit_ids.begin(); }
  PHG4Shower::HitIdIter find_g4hit_id(int volume) { return _g4hit_ids.find(volume); }
  PHG4Shower::HitIdConstIter find_g4hit_id(int volume) const { return _g4hit_ids.find(volume); }
  PHG4Shower::HitIdIter end_g4hit_id() { return _g4hit_ids.end(); }
  PHG4Shower::HitIdConstIter end_g4hit_id() const { return _g4hit_ids.end(); }
  size_t remove_g4hit_id(int volume, PHG4HitDefs::keytype id) { return _g4hit_ids[volume].erase(id); }
  size_t remove_g4hit_volume(int volume) { return _g4hit_ids.erase(volume); }
  void clear_g4hit_id() { return _g4hit_ids.clear(); }

 private:
  unsigned int covar_index(unsigned int i, unsigned int j) const;

  int _id;                             //< unique identifier within container
  int _parent_particle_id;             //< association of shower to parent particle id
  int _parent_shower_id;               //< association of shower to parent shower id
  float _pos[3];                       //< mean position of the shower hits
  float _covar[6];                     //< covariance of shower hits
  std::map<int, unsigned int> _nhits;  //< number of hits in different volumes
  std::map<int, float> _edep;          //< energy deposit in different volumes
  std::map<int, float> _eion;          //< ionization energy in different volumes
  std::map<int, float> _light_yield;   //< light yield in different volumes
  std::map<int, float> _eh_ratio;      //< electron/hadron ratio of energy in different volumes

  // these containers are cleared during dst reduction, but are available in full dsts
  std::set<int> _g4particle_ids;                            //< contained secondary particle ids
  std::set<int> _g4vertex_ids;                              //< contained secondary vertex ids
  std::map<int, std::set<unsigned long long> > _g4hit_ids;  //< contained hit ids

  ClassDef(PHG4Showerv1, 1);
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PHG4Showerv1", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PHG4Showerv1_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PHG4Showerv1_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PHG4Showerv1_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PHG4Showerv1_Dict() {
  TriggerDictionaryInitialization_PHG4Showerv1_Dict_Impl();
}
