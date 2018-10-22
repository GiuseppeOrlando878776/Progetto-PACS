#ifndef SU2_VECTORT_HPP
#define SU2_VECTORT_HPP

#include "exprT.hpp"
#include <vector>
#include <exception>

namespace Common {

 /*!
  * Definition of a class Su2Vec that implements an expression template technique
  * \author G.Orlando
  */

  template <typename T,unsigned ID_ENUM,char* ID_FILE>
  class SU2Vec : public ExprT<SU2Vec<T,ID_ENUM,ID_FILE>,T> {
  public:

    using Type = T;

    /*!
      * \brief Default constructor
    */
    SU2Vec() : N() {
      data = NULL;
    }

    /*!
      * \brief Class constructor
      * \param[in] _N - size of the vector
      * \param[in] value - value to initialize (default = 0.0)
    */
    SU2Vec(unsigned _N,Type value = 0.0);

    /*!
      * \brief Class destructor
    */
    ~SU2Vec();

    /*!
      * \brief Copy constructor
      * \param[in] su2vec - SU2Vec to copy from
    */
    SU2Vec(const SU2Vec& su2vec);

    /*!
      * \brief Move constructor
      * \param[in] su2vec - SU2Vec to copy from
    */
    SU2Vec(SU2Vec&& su2vec);

    /*!
      * \brief Copy assignment operator
      * \param[in] su2vec - SU2Vec to copy from
    */
    SU2Vec& operator=(const SU2Vec& su2vec);

    /*!
      * \brief Move assignment operator
      * \param[in] su2vec - SU2Vec to copy from
    */
    SU2Vec& operator=(SU2Vec&& su2vec);

    /*!
      * \brief Copy Constructor from standard library vector
      * \param[in] v - std::vector to copy from
    */
    explicit SU2Vec(const std::vector<Type>& v);

    /*!
      * \brief Move Constructor from standard library vector
      * \param[in] v - std::vector to copy from
    */
    explicit SU2Vec(std::vector<Type>&& v);

    /*!
      * \brief Construct from ExprT;
      * \param[in] e - ExprT to copy from
     */
    template <class Derived>
    SU2Vec& operator=(const ExprT<Derived,Type>& e);

    /*!
     * \brief Returns i-th element (non const version)
     * \param[in] i - index of the element
    */
    Type& operator[](std::size_t i);

    /*!
     * \brief Returns i-th element (const version)
     * \param[in] i - index of the element
    */
    const Type& operator[](std::size_t i) const;

    /*!
     * \brief Returns i-th element (non const version)
     * \param[in] i - index of the element
    */
    static inline Type& at(std::size_t i);

    /*!
     * \brief Returns i-th element (const version)
     * \param[in] i - index of the element
    */
    //const Type& at(std::size_t i) const;

    /*!
     * \brief Returns the size of the vector
    */
    std::size_t size() const {
      return N;
    }

    /*!
     * \brief Cast to standard library vector (const version)
    */
    operator const std::vector<Type>& () const;

    /*!
     * \brief Cast to standard library vector (non const version)
    */
    operator std::vector<Type>& ();

  private:

    std::size_t N;  /*!< \brief Size of the vector. */
    static Type* data;  /*!< \brief Stored elements. */

  }; /*--- End of class SU2Vec ---*/

} /*--- End of namespace Common ---*/

//////////////////////////////////////////////////////////////////////////////

#endif
