#ifndef SU2_VECTORT_HPP
#define SU2_VECTORT_HPP

#include "exprT.hpp"
#include <vector>

namespace Common {

 /*!
  * Definition of a class Su2Vec that implements an expression template technique
  * \author G.Orlando
  */

  template <typename T,unsigned ID_ENUM,typename S,S ID_FILE>
  class SU2Vec : public ExprT<SU2Vec<T,ID_ENUM,S,ID_FILE>,T> {
  public:

    using Type = T;

    /*!
      * \brief Default constructor
    */
    SU2Vec() : N() {
      m_data = NULL;
    }

    /*!
      * \brief Class constructor
      * \param[in] _N - size of the vector
      * \param[in] init - value to initialize (default provided)
    */
    SU2Vec(std::size_t _N,Type init = Type());

    /*!
      * \brief Class constructor from already allocated memory
      * \param[in] _N - size of the vector
      * \param[in] values - value to initialize (default = 0.0)
    */
    SU2Vec(std::size_t _N,Type* values);

    /*!
      * \brief Class destructor
    */
    ~SU2Vec() {
      free();
    }

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
      * \brief Assignment opertor from ExprT;
      * \param[in] e - ExprT to copy from
     */
    template <class Derived>
    inline void operator=(const ExprT<Derived,Type>& e) {
      for(auto i = 0; i < N; ++i)
        m_data[i] = Derived::at(i);
    }

    /*!
     * \brief Returns i-th element (non const version)
     * \param[in] i - index of the element
    */
    inline Type& operator[](std::size_t i) {
      return m_data[i];
    }

    /*!
     * \brief Returns i-th element (const version)
     * \param[in] i - index of the element
    */
    const Type& operator[](std::size_t i) const {
      return m_data[i];
    }

    /*!
     * \brief Returns i-th element (non const version)
     * \param[in] i - index of the element
    */
    static inline Type& at(std::size_t i) {
      SU2_Assert(i < size(),"Index is beyond the size of the vector");
      return m_data[i];
    }

    /*!
     * \brief Returns i-th element (const version)
     * \param[in] i - index of the element
    */
    //const Type& at(std::size_t i) const;

    /*!
     * \brief Returns the size of the vector
    */
    inline unsigned id() const {
      return ID_ENUM;
    }

    /*!
     * \brief Returns the size of the vector
    */
    inline std::size_t size() const {
      return N;
    }

    /*!
     * \brief Resize of the vector
     * \param[in] _N - new size of the vector
     * \param[in] init - value to initialize (default provided)
    */
    void resize(std::size_t _N, T init = T());

    /*!
     * \brief Returns the underlined pointer
    */
    Type* data() const {
      return m_data;
    }

    /*!
     * \brief Cast to standard library vector (const version)
    */
    operator const std::vector<Type>& () const {
      return static_cast<const std::vector<Type>&>(*this);
    }

    /*!
     * \brief Cast to standard library vector (non const version)
    */
    operator std::vector<Type>& () {
      return static_cast<std::vector<Type>&>(*this);
    }

  private:

    std::size_t N;  /*!< \brief Size of the vector. */
    static Type* m_data;  /*!< \brief Stored elements. */

    /*!
     * \brief Helper function to allocate memory
    */
    inline void allocate(void) {
      if (size() > 0)
        m_data = new Type[N];
    }

    /*!
     * \brief Helper function to initialize vector
    */
    inline void initialize(Type value) {
      for(auto i = 0; i < size(); ++i)
        m_data[i] = value;
    }

    /*!
     * \brief Helper function to free memory
    */
    inline void free(void) {
      if (m_data!= NULL) {
        delete [] m_data;
        m_data = NULL;
        N = 0;
      }
    }

  }; /*--- End of class declaration SU2Vec ---*/

  //
  //
  /*--- SU2Vec value constructor---*/
  template<typename T,unsigned ID_ENUM,typename S,S ID_FILE>
  SU2Vec<T,ID_ENUM,S,ID_FILE>::SU2Vec(std::size_t _N,T init):N(_N) {
    allocate();
    initiliaze(init);
  }

  //
  //
  /*--- SU2Vec preallocated memory constructor---*/
  template<typename T,unsigned ID_ENUM,typename S,S ID_FILE>
  SU2Vec<T,ID_ENUM,S,ID_FILE>::SU2Vec(std::size_t _N,T* values):N(_N) {
    m_data = values;
  }

  //
  //
  /*--- SU2Vec copy constructor---*/
  template<typename T,unsigned ID_ENUM,typename S,S ID_FILE>
  SU2Vec<T,ID_ENUM,S,ID_FILE>::SU2Vec(const SU2Vec& su2vec):N(su2vec.N) {
    allocate();
    for(auto i = 0; i<size(); ++i)
      m_data[i] = su2vec.m_data[i];
  }

  //
  //
  /*--- SU2Vec move constructor---*/
  template<typename T,unsigned ID_ENUM,typename S,S ID_FILE>
  SU2Vec<T,ID_ENUM,S,ID_FILE>::SU2Vec(SU2Vec&& su2vec):N(su2vec.N) {
    m_data = su2vec.m_data;
    su2vec.m_data = NULL;
    su2vec.N = 0;
  }

  //
  //
  /*--- SU2Vec copy assignment operator ---*/
  template<typename T,unsigned ID_ENUM,typename S,S ID_FILE>
  SU2Vec<T,ID_ENUM,S,ID_FILE>& SU2Vec<T,ID_ENUM,S,ID_FILE>::operator=(const SU2Vec& su2vec) {
      if(this!= su2vec) {
          N = su2vec.N;
          delete[] m_data;
          m_data = new T(N);
          for(auto i = 0; i<size(); ++i)
            m_data[i] = su2vec.m_data[i];
      }
      return *this;
  }

  //
  //
  /*--- SU2Vec move assignment operator ---*/
  template<typename T,unsigned ID_ENUM,typename S,S ID_FILE>
  SU2Vec<T,ID_ENUM,S,ID_FILE>& SU2Vec<T,ID_ENUM,S,ID_FILE>::operator=(SU2Vec&& su2vec) {
      if(this!= su2vec) {
          N = su2vec.N;
          delete[] m_data;
          m_data = su2vec.m_data;
          su2vec.m_data = NULL;
          su2vec.N = 0;
      }
      return *this;
  }

  //
  //
  /*--- SU2Vec copy constructor from standard library vector ---*/
  template<typename T,unsigned ID_ENUM,typename S,S ID_FILE>
  SU2Vec<T,ID_ENUM,S,ID_FILE>::SU2Vec(const std::vector<T>& v):N(v.size()) {
    allocate();
    for(auto i = 0; i<size(); ++i)
      m_data[i] = v[i];
  }

  //
  //
  /*--- SU2Vec move constructor from standard library vector ---*/
  template<typename T,unsigned ID_ENUM,typename S,S ID_FILE>
  SU2Vec<T,ID_ENUM,S,ID_FILE>::SU2Vec(std::vector<T>&& v):N(v.size()) {
    auto ptr = v.data();
    m_data = ptr;
    ptr = NULL;
  }

  //
  //
  /*--- Resize container ---*/
  template<typename T,unsigned ID_ENUM,typename S,S ID_FILE>
  void SU2Vec<T,ID_ENUM,S,ID_FILE>::resize(std::size_t _N,T init) {
    free();
    N = _N;
    allocate();
    initiliaze(init);
  }

  //
  //
  /*--- SU2Vec begin iterator (non const version)---*/
  template<typename T,unsigned ID_ENUM,typename S,S ID_FILE>
  inline auto begin(SU2Vec<T,ID_ENUM,S,ID_FILE>& a)->decltype(std::declval<std::vector<T>>().begin()) {
    return static_cast<std::vector<T>&>(a).begin();
  }

  //
  //
  /*--- SU2Vec end iterator (non cosnt version)---*/
  template<typename T,unsigned ID_ENUM,typename S,S ID_FILE>
  inline auto end(SU2Vec<T,ID_ENUM,S,ID_FILE>& a)->decltype(std::declval<std::vector<T>>().end()) {
    return static_cast<std::vector<T>&>(a).end();
  }

  //
  //
  /*--- SU2Vec begin iterator (const version)---*/
  template<typename T,unsigned ID_ENUM,typename S,S ID_FILE>
  inline auto cbegin(SU2Vec<T,ID_ENUM,S,ID_FILE>& a)->decltype(std::declval<std::vector<T>>().cbegin()) {
    return static_cast<std::vector<T>&>(a).cbegin();
  }

  //
  //
  /*--- SU2Vec end iterator (const version)---*/
  template<typename T,unsigned ID_ENUM,typename S,S ID_FILE>
  inline auto cend(SU2Vec<T,ID_ENUM,S,ID_FILE>& a)->decltype(std::declval<std::vector<T>>().cend()) {
    return static_cast<std::vector<T>&>(a).cend();
  }

  #define RealVec SU2Vec<double,__LINE__,std::string*, __FILE__>

} /*--- End of namespace Common ---*/

#endif
