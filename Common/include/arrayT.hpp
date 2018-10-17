#ifndef SU2__ARRAYT_HPP
#define SU2_ARRAYT_HPP

#include <cmath>
#include <algorithm>
#include <cassert>
#include <numeric>
//#include "MathTools/MacrosET.hh"

namespace Common {


#define ETVEC(a,b,c,d)   ArrayT<a<b<c,d>, ETuple<c,d,0> > >
#define ETMAT(a,b,c,d,e) ArrayT<a<b<c,d,e>, ETuple<c,d,e> > >
#define EDATAPTR(a)      this->getData()->ptr()[a]


/** Definition of a class ArrayT that implements an expression template technique
  * with statically defined size
  * @author G.Orlando
  */
template <typename Derived>
class ArrayT : public Derived {
public:

  typedef typename Derived::Tuple::Type Type;

  /*
   * \brief Default constructor
   */
  ArrayT(typename Derived::Ptr* der): Derived(der) {}

  /*
   * \brief Default destructor
  */
  virtual ~ArrayT() {}

  /*!
   * \brief Overloading of operator [] for accessing
   * \param[in] iElem - Index of the element.
   */
  inline Type& operator[] (std::size_t iElem) {
    return EDATAPTR(iElem);
  }

  /*!
   * \brief Overloading of operator [] for accessing (const version)
   * \param[in] iElem - Index of the element.
   */
   inline const Type& operator[] (size_t iElem) const {
     return EDATAPTR(iElem);
   }

  /*!
   * \brief Function for accessing element
   * \param[in] iElem - Index of the element.
   */
  inline Type& at(std::size_t iElem) {
    static_assert(iElem < size(),"The index exceeds the size");
    return EDATAPTR(iElem);
  }

  /*!
   * \brief Function for accessing element (const version)
   * \param[in] iElem - Index of the element.
   */
  inline const Type& at(std::size_t iElem) const {
    static_assert(iElem < size(),"The index exceeds the size");
    return EDATAPTR(iElem);
  }

  /*!
   * \brief Return the size of the array
   */
  inline std::size_t size() const {
    return this->size();
  }

  /*!
   * \brief Return the maximum within the array
   */
  inline Type max() const  {
    static_assert(size() > 0,"Empty array");
    return std::max_element(this->getData(),this->getData() + size());
  }

  /*!
   * \brief Return the minimum within the array
   */
  inline Type min() const  {
    static_assert(size() > 0,"Empty array");
    return std::min_element(this->getData(),this->getData() + size());
  }

  /*!
   * \brief Return the sum of all the elements within the array
   */
  inline Type sum() const  {
    static_assert(size() > 0,"Empty array");
    return std::accumulate(this->getData(),this->getData() + size(),0.0,std::plus<Type>());
  }

  /*!
   * \brief Return the sum of all the elements within a certain range
   */
  inline Type partialSum(std::size_t init,std::size_t end) const {
    static_assert(size()>=end,"Not enough elements");
    return std::accumulate(this->getData() + init,this->getData() + end,0.0,std::plus<Type>());
  }

  /*!
   * \brief Return the norm2
   */
  inline Type norm2() const {
    return std::sqrt(sqrNorm());
  }

  /*!
   * \brief Return the norm1
   */
  inline Type norm1() const {
    static_assert(size() > 0,"Empty array");
    Type norm = std::abs(EDATAPTR(0));
    for (std::size_t i = 1; i < size(); ++i)
      norm += std::abs(EDATAPTR(i));
    return norm;
  }

  /*!
   * \brief Return the infinity norm
   */
  inline Type normInf() const {
    return std::max(std::abs(max()),std::abs(min()));
  }

  /*!
   * \brief Return the squred L2 norm
   */
   inline Type sqrNorm() const  {
    static_assert(size() > 0,"Empty array");
    return std::inner_product(this->getData(),this->getData()+size(),this->getData(),0.0);
  }

protected:

  /*!
   * \brief Return the infinity norm
   */
  inline void initialize(const Type& value) {
    std::fill(this->getData(),this->getData()+size(),value);
};

}; /*--- End of class ArrayT ----*/

} /*--- End of namespace Common ----*/

#endif
