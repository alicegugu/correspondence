//*****************************************************************************
// Non-Intrusive Simple SmartPtr (RefCount based)
//*****************************************************************************

#ifndef E_SMARTPTR_H
#define E_SMARTPTR_H

#include <string>
#include <assert.h>

namespace cwc
{
   template< class T > class SmartPtr;

   template< class T >
   class SmartPtr 
   {
   public:
      SmartPtr( T* ptr = 0 ) : _ref( new int( 1 ) ), _ptr( ptr ) {}
      SmartPtr( const SmartPtr& sp )
      {
         Copy(sp);
      }
      template< class X > friend class SmartPtr;
      template< class X >
      SmartPtr( const SmartPtr< X >& sp ) 
      {
         ++*sp._ref;
         _ptr = sp._ptr;
         _ref = sp._ref;
      }
      ~SmartPtr() 
      {
         Decrease();
      }

   public:
      SmartPtr& operator=( const SmartPtr& sp ) 
      {
         if ( this != &sp ) 
         {
	         Decrease();
	         Copy(sp);
         }
         return *this;
      }
      
      T* operator->() 
      {
         return _ptr;
      }
      
      const T* operator->() const 
      {
         return _ptr;
      }

      operator bool() const 
      { 
         return _ptr != 0; 
      }

   private:
      void Copy( const SmartPtr& sp ) 
      {
         ++*sp._ref;
         _ref = sp._ref;
         _ptr = sp._ptr;
      }
      
      void Decrease() 
      {
         --*_ref;
         if ( 0 == *_ref ) 
         {
	         delete _ref;
	         delete _ptr; 
	         _ptr = 0;
         }
      }

   private:
      int* _ref;
      T* _ptr;
   };
}

#endif