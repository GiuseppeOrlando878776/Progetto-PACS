#ifndef SU2_VECTORT_HPP
#define SU2_VECTORT_HPP

#include "exprT.hpp"
#include <iomanip>
#include <algorithm>

namespace Common {

  /*!
   * Forward declaration of SU2Vec class
   */
  template<typename T> class SU2Vec;

  /*!
   * Declaration of output stream operator
   */
  template<typename T> std::ostream& operator<<(std::ostream& out,const SU2Vec<T>& v);

  /*!
   * Definition of a class SU2Vec that implements an expression template technique
   * \author G.Orlando
  */

  template <typename T>
  class SU2Vec : public ExprT<SU2Vec<T>,T> {
  public:

    using Type = T;
    using value_type = T;
    typedef typename ExprT<SU2Vec<T>,T>::size_type size_type;

    /*!
      * \brief Default constructor
    */
    SU2Vec(): N(), m_data(NULL) {}

    /*!
      * \brief Class constructor
      * \param[in] _N - size of the vector
      * \param[in] init - value to initialize (default provided)
    */
    explicit SU2Vec(size_type _N,Type init = Type());

    /*!
      * \brief Class destructor
    */
    ~SU2Vec() {
      free_mem();
    }

    /*!
      * \brief Copy constructor
      * \param[in] su2vec - SU2Vec to copy from
    */
    SU2Vec(const SU2Vec& su2vec);

    /*!
      * \brief Move constructor
      * \param[in] su2vec - SU2Vec to move from
    */
    SU2Vec(SU2Vec&& su2vec) noexcept;

    /*!
      * \brief Copy assignment operator
      * \param[in] su2vec - SU2Vec to copy from
    */
    SU2Vec& operator=(const SU2Vec& su2vec);

    /*!
      * \brief Move assignment operator
      * \param[in] su2vec - SU2Vec to move from
    */
    SU2Vec& operator=(SU2Vec&& su2vec) noexcept;

    /*!
      * \brief Copy Constructor with range
      * \param[in] init - iterator to first datum
      * \param[in] end - iterator to one past last datum
    */
    template<typename Iterator>
    SU2Vec(Iterator b,Iterator e) {
      N = std::distance(b,e);
      m_data = new Type[N];
      for(size_type i = 0; i<N; ++i)
        m_data[i] = *b++;
    }

    /*!
      * \brief Copy Constructor from standard library vector
      * \param[in] v - std::vector to copy from
    */
    explicit SU2Vec(const std::vector<Type>& v);

    /*!
      * \brief Copy constructor from ExprT;
      * \param[in] e - ExprT to copy from
     */
    template <class Derived>
    SU2Vec(const ExprT<Derived,Type>& e):N(e.size()) {
      //const Derived& et(e);
      allocate();
      for(size_type i = 0; i < N; ++i)
        m_data[i] = e[i];
    }

    /*!
      * \brief Assignment operator from ExprT;
      * \param[in] e - ExprT to assign
     */

    template <class Derived>
    SU2Vec& operator=(ExprT<Derived,Type>& e) {
      //const Derived& et(e);
      N = e.size();
      if(m_data!=NULL)
        delete[] m_data;
      allocate();
      for(size_type i = 0; i < N; ++i)
        m_data[i] = e[i];
      return *this;
    }

    /*!
      * \brief Overloaded assignment operator with constants (macro used);
     */

    #define SU2VEC_ASSIGN_OP_CONST(__op__) \
    SU2Vec& operator __op__ (Type value) { \
      SU2_Assert(strcmp(#__op__ ,"/=") && std::abs(value)>0,"You can't divide by zero"); \
      for(size_type i = 0; i < size(); ++i) \
        m_data[i] __op__ value; \
      return *this; \
    }
    SU2VEC_ASSIGN_OP_CONST(+=)
    SU2VEC_ASSIGN_OP_CONST(-=)
    SU2VEC_ASSIGN_OP_CONST(*=)
    SU2VEC_ASSIGN_OP_CONST(/=)
    #undef SU2VEC_ASSIGN_OP_CONST

    /*!
     * \brief Overloaded assignment operator with ExprT (macro used);
    */
    #define SU2VEC_ASSIGN_OP_EXPR(__op__) \
    template<class Derived> \
    SU2Vec& operator __op__ (const ExprT<Derived,Type>& e) { \
      SU2_Assert(this->size() == e.size(),"The size of the vectors is not the same and so you can't operate with them"); \
      for(size_type i = 0; i < size(); ++i) \
        m_data[i] __op__ e[i]; \
      return *this; \
    }
    //SU2VEC_ASSIGN_OP_EXPR(=)
    SU2VEC_ASSIGN_OP_EXPR(+=)
    SU2VEC_ASSIGN_OP_EXPR(-=)
    SU2VEC_ASSIGN_OP_EXPR(*=)
    //SU2VEC_ASSIGN_OP_CONST(/=)
    #undef SU2VEC_ASSIGN_OP_EXPR

    /*!
     * \brief Overloaded /= operator with ExprT (checks if rhs contains zeros);
    */
    template<class Derived>
    SU2Vec& operator /= (const ExprT<Derived,Type>& e) {
      SU2_Assert(this->size() == e.size(),"The size of the vectors is not the same and so you can't operate with them");
      SU2Vec tmp(e);
      SU2_Assert(tmp.check_notzero(),"The vector on rhs contains zeros and so you can't divide by it");
      for(size_type i = 0; i < size(); ++i) \
        m_data[i] /= tmp[i]; \
      return *this;
    }


    /*!
     * \brief Returns i-th element (non const version)
     * \param[in] i - index of the element
    */
    inline Type& operator[](size_type i) {
      return m_data[i];
    }

    /*!
     * \brief Returns i-th element (const version)
     * \param[in] i - index of the element
    */
    //inline Type operator[](size_type i) const
    inline const Type& operator[](size_type i) const {
      return m_data[i];
    }

    /*!
     * \brief Returns i-th element (non const version)
     * \param[in] i - index of the element
    */
    inline Type& at(size_type i) {
      SU2_Assert(i < size(),"Index is beyond the size of the vector");
      return m_data[i];
    }

    /*!
     * \brief Returns i-th element (const version)
     * \param[in] i - index of the element
    */
    //inline Type at(size_type i) const
    inline const Type& at(size_type i) const {
      SU2_Assert(i < size(),"Index is beyond the size of the vector");
      return m_data[i];
    }

    /*!
     * \brief Returns the size of the vector
    */
    inline size_type size(void) const {
      return N;
    }

    /*!
     * \brief Clears the vector
    */
    inline void clear(void) {
      free_mem();
    }

    /*!
     * \brief Resize of the vector
     * \param[in] _N - new size of the vector
     * \param[in] init - value to initialize (default provided)
    */
    void resize(size_type _N, Type init = Type());

    /*!
     * \brief Add an element at the end of the vector
     * \param[in] value - value to add
    */
    void push_back(const Type& value);

    /*!
     * \brief Check if there are zeros in the vector
    */
    inline bool check_notzero(void) {
      return std::none_of(this->cbegin(), this->cend(), [&](const Type& x){return x == Type();});
    }

    /*!
     * \brief Returns the underlined pointer
    */
    inline Type* data(void) const {
      return m_data;
    }

    /*!
     * \brief Cast to standard library vector
    */
    inline operator std::vector<Type> () const {
      return std::vector<Type>(m_data,m_data + N);
    }

    /*!
     * \brief Begin iterator (non const version)
    */
    inline auto begin()->decltype(std::declval<std::vector<Type>>().begin()) {
      return static_cast<decltype(std::declval<std::vector<Type>>().begin())>(m_data);
    }

    /*!
     * \brief End iterator (non const version)
    */
    inline auto end()->decltype(std::declval<std::vector<Type>>().end()) {
      return static_cast<decltype(std::declval<std::vector<Type>>().end())>(m_data + N);
    }

    /*!
     * \brief Begin iterator (const version)
    */
    inline decltype(std::declval<std::vector<Type>>().cbegin()) cbegin() const {
      return static_cast<decltype(std::declval<std::vector<Type>>().cbegin())>(m_data);
    }

    /*!
     * \brief End iterator (const version)
    */
    inline decltype(std::declval<std::vector<Type>>().cend()) cend() const {
      return static_cast<decltype(std::declval<std::vector<Type>>().cend())>(m_data + N);
    }

    /*!
     * \brief Overloading of output stream operator
    */
    friend std::ostream& operator<< TF (std::ostream& out,const SU2Vec<Type>& v);

  private:

    size_type N;  /*!< \brief Size of the vector. */
    Type* m_data;  /*!< \brief Stored elements. */

    /*!
     * \brief Helper function to allocate memory
    */
    void allocate(void);

    /*!
     * \brief Helper function to initialize vector
     * \param[in] value - Value initialization
    */
    void initialize(Type value);

    /*!
     * \brief Helper function to free memory
    */
    void free_mem(void);

  }; /*--- End of class declaration SU2Vec ---*/

  //
  //
  /*--- Allocate memory---*/
  template<typename T>
  void SU2Vec<T>::allocate(void) {
    if (size() > 0)
      m_data = new T[N];
  }

  //
  //
  /*--- Initialize vector ---*/
  template<typename T>
  void SU2Vec<T>::initialize(T value) {
    for(size_type i = 0; i < size(); ++i)
      m_data[i] = value;
  }

  //
  //
  /*--- free memory ---*/
  template<typename T>
  void SU2Vec<T>::free_mem(void) {
    if (m_data!= NULL) {
      delete[] m_data;
      m_data = NULL;
      N = 0;
    }
  }

  //
  //
  /*--- SU2Vec value constructor---*/
  template<typename T>
  SU2Vec<T>::SU2Vec(size_type _N,T init):N(_N) {
    allocate();
    initialize(init);
  }

  //
  //
  /*--- SU2Vec copy constructor---*/
  template<typename T>
  SU2Vec<T>::SU2Vec(const SU2Vec& su2vec):N(su2vec.N) {
    allocate();
    for(size_type i = 0; i < N; ++i)
      m_data[i] = su2vec.m_data[i];
  }

  //
  //
  /*--- SU2Vec move constructor---*/
  template<typename T>
  SU2Vec<T>::SU2Vec(SU2Vec&& su2vec) noexcept:N(su2vec.N) {
    m_data = su2vec.m_data;
    su2vec.m_data = NULL;
    su2vec.N = 0;
  }

  //
  //
  /*--- SU2Vec copy assignment operator ---*/
  template<typename T>
  SU2Vec<T>& SU2Vec<T>::operator=(const SU2Vec& su2vec) {
      if(this!= &su2vec) {
          N = su2vec.N;
          if(m_data!= NULL)
            delete[] m_data;
          m_data = new T[N];
          for(size_type i = 0; i < N; ++i)
            m_data[i] = su2vec.m_data[i];
      }
      return *this;
  }

  //
  //
  /*--- SU2Vec move assignment operator ---*/

  template<typename T>
  SU2Vec<T>& SU2Vec<T>::operator=(SU2Vec&& su2vec) noexcept {
      if(this!= &su2vec) {
        N = std::move(su2vec.N);
        if(m_data!= NULL)
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
  template<typename T>
  SU2Vec<T>::SU2Vec(const std::vector<T>& v):N(v.size()) {
    allocate();
    for(size_type i = 0; i<size(); ++i)
      m_data[i] = v[i];
  }

  //
  //
  /*--- Resize container ---*/
  template<typename T>
  void SU2Vec<T>::resize(size_type _N,T init) {
    if(_N == size())
      return;
    if(_N < size()) {
      N = _N;
      auto tmp = new T[size()];
      for(size_type i = 0; i < size(); ++i)
        tmp[i] = m_data[i];
      delete[] m_data;
      m_data = tmp;
    }
    else {
      auto tmp = new T[_N];
      size_type i;
      std::fill(tmp + size(),tmp + _N,init);
      for(i = 0; i < size(); ++i)
        tmp[i] = m_data[i];
      N = _N;
      delete[] m_data;
      m_data = tmp;
    }
  }

  //
  //
  /*--- Add an element at the end of the container ---*/
  template<typename T>
  void SU2Vec<T>::push_back(const T& value) {
    T* tmp = new T[++N];
    size_type i;
    for(i = 0; i < size() - 1; ++i)
      tmp[i] = m_data[i];
    tmp[i] = value;
    if(m_data!=NULL)
      delete[] m_data;
    m_data = tmp;
  }

  //
  //
  /*--- Stream output operator ---*/
  template <typename T>
  std::ostream& operator<<(std::ostream& out,const SU2Vec<T>& v) {
    std::cout<<std::fixed<<std::setprecision(4);
    for(std::size_t i = 0; i < v.size(); ++i)
      out<<std::right<<v.m_data[i]<<std::endl;
    //out<<'\n';
    return out;
  }

  //
  //
  /*--- Operator / between two vectors ---*/
  template<class Derived_Left,class Derived_Right,typename T>
  SU2Vec<T> operator/(const ExprT<Derived_Left,T>& l, const ExprT<Derived_Right,T>& r) {
    SU2_Assert(l.size() == r.size(),"The size of the vectors is not the same and so you can't operate with them");
    SU2Vec<T> result(r);
    SU2_Assert(result.check_notzero(),"The vector on rhs contains zero and so you can't divide by it");
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = l[i]/result[i];
    return result;
  }

  //
  //
  /*--- Operator / between a constant and a vector ---*/
  template<class Derived,typename T>
  SU2Vec<T> operator/(double l, const ExprT<Derived,T>& r) {
    static_assert(std::is_convertible<T,double>::value,"The type of the vector is not convertible to double"); \
    SU2Vec<T> result(r);
    SU2_Assert(result.check_notzero(),"The vector on rhs contains zero and so you can't compute its reciprocal");
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = l/result[i];
    return result;
  }

  using RealVec = SU2Vec<double>;

} /*--- End of namespace Common ---*/

#endif
