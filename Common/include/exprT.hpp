#ifndef SU2_EXPRT
#define SU2_EXPRT

#include <iostream>
#include <type_traits>
#include "su2_assert.hpp"

namespace Common {

  /**
    * Definition of a wrapper base expression class ExprT: derived objects from ExprT use the CRTP technique
    * The expression template accepts two template parameters:
    * 1- the derived expression
    * 2- the return type of the epxression
    * \author G. Orlando
  */

  template <class Derived, typename T>
  struct ExprT {

    using Type = T;

    /*!
      *\brief Cast operator to derived class (const version)
    */
    inline operator const Derived&() const {
      return static_cast<const Derived&>(*this);
    }

    /*!
     *\brief  Cast operator to derived class (non const version)
    */
    inline operator Derived&() {
      return static_cast<Derived>(*this);
    }

   /*!
    * \brief Access operator (const version)
   */
   inline const Type& operator[](std::size_t i) const {
     return this->operator const Derived&().operator[](i);
   }

   /*!
    * \brief Access operator (non const version)
   */
   inline Type& operator[](std::size_t i) {
     return this->operator Derived&().operator[](i);
   }

   /*!
    * \brief Access operator (const version)
   */
   static inline Type& at(std::size_t i) {
     return Derived::at(i);
   }

   /*!
    * \brief Access operator (non const version)
   */
   /*
   inline const Type& at(std::size_t i) const {
     return this->operator const Derived&().at(i);
   }
   */

   /*!
    * \brief Size of the derived object
   */
   std::size_t size() const {
     return this->operator const Derived&().size();
   }

 }; /*--- End of class ExprT ----*/

  /**
   * Definition of an expression template class for basic binary operations.
   * A macro is used to save code duplication.
   * The expression template accepts two parameters:
   * 1. left operand type
   * 2. right operand type.
   *
   * \author G. Orlando
   */

   #define TYPE(a) typename a::Type
   #define EETYPE(a) ExprT<a,TYPE(a)>

   #define ET_BINARY(OpName,operation) \
   template <class Left, class Right>			\
   struct OpName : public ExprT<OpName<Left,Right>, TYPE(Left)> {	\
     static_assert(std::is_convertible<TYPE(Left),TYPE(Right)>::value,"The types are not convertible"); \
                          \
     OpName(EETYPE(Left) l, EETYPE(Right) r) :	    \
     ExprT<OpName<Left,Right>, TYPE(Left)>(this), e1(l), e2(r) {}	\
    									\
     inline TYPE(Left)& operator[](std::size_t i) { \
       return operation;\
     } \
                        \
     inline TYPE(Left)& operator[](std::size_t i) const { \
       return operation;\
     } \
                      \
     static inline TYPE(Left)& at(std::size_t i) { \
       return operation;\
     } \
                      \
     /*TYPE(Left)& at(std::size_t i) const {*/ \
    /*   return operation;*/\
    /* }*/ \
									\
    std::size_t size() const { \
      SU2_Assert(e1.size() == e2.size(),"The size of vectors is not the same"); \
      return e1.size();\
    }		\
                          \
   private:								\
      const EETYPE(Left)& e1;							\
      const EETYPE(Right)& e2;							\
  };

  ET_BINARY(Add, Left::at(i)+Right::at(i))
  ET_BINARY(Sub, Left::at(i)-Right::at(i))
  ET_BINARY(Mul, Left::at(i)*Right::at(i))
  ET_BINARY(Div, Left::at(i)/Right::at(i))
  //ET_BINARY(Max, std::max(e1.at(i),e2.at(i)))
  //ET_BINARY(Min, std::min(e1.at(i),e2.at(i)))
  #undef ET_BINARY

  template <class Left, class Right>
  inline Add<Left,Right> operator+(const Left& l, const Right& r){return  Add<Left,Right>(l,r);}

  template <class Left, class Right>
  inline Sub<Left,Right> operator-(const Left& l, const Right& r){return  Sub<Left,Right>(l,r);}

  template <class Left, class Right>
  inline Mul<Left,Right> operator*(const Left& l, const Right& r){return  Mul<Left,Right>(l,r);}

  template <class Left>
  inline Mul<Left,double> operator*=(const Left& l, const double& r){return  Mul<Left,double>(l,r);}

  template <class Left, class Right>
  inline Div<Left,Right> operator/(const Left& l, const Right& r){return  Div<Left,Right>(l,r);}

  template <class Left>
  inline Div<Left,double> operator/=(const Left& l, const double& r){return  Div<Left,double>(l,r);}


}   /*--- End of namespace Common ----*/

#endif
