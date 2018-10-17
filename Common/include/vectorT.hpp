#ifndef SU2_VECTORT_HPP
#define SU2_VECTORT_HPP

#include "arrayT.hpp"
//#include "MathTools/CFVecSlice.hh"
//#include "MathTools/MathFunctions.hh"


namespace Common {

  template  <typename T, int N> class RealVec;

  template <typename T, int N> std::ostream& operator<< (std::ostream& out, const RealVec<T,N>& v);
  template <typename T, int N> std::istream& operator>> (std::istream& in, RealVec<T,N>& v);
  template <typename T, int N> bool operator== (const RealVec<T,N>& v1, const RealVec<T,N>& v2);
  template <typename T, int N> bool operator== (const RealVec<T,N>& v1, T value);
  template <typename T, int N> bool operator!= (const RealVec<T,N>& v1, const RealVec<T,N>& v2);
  template <typename T, int N> bool operator!= (const RealVec<T,N>& v1, T value);
  template <typename T, int N> bool operator<  (const RealVec<T,N>& v1, const RealVec<T,N>& v2);
  template <typename T, int N> bool operator<  (const RealVec<T,N>& v1, T value);
  template <typename T, int N> bool operator>  (const RealVec<T,N>& v1, const RealVec<T,N>& v2);
  template <typename T, int N> bool operator>  (const RealVec<T,N>& v1, T value);
  template <typename T, int N> bool operator<= (const RealVec<T,N>& v1, const RealVec<T,N>& v2);
  template <typename T, int N> bool operator<= (const RealVec<T,N>& v1, T value);
  template <typename T, int N> bool operator>= (const RealVec<T,N>& v1, const RealVec<T,N>& v2);
  template <typename T, int N> bool operator>= (const RealVec<T,N>& v1, T value);
  template <typename T, int N> void copy       (const RealVec<T,N>& v1, RealVec<T,N>& v2);


/*
 * Definition of a class RealVec that implements an expression template technique
 * with statically defined size
 * @author G.Orlando
 */
template <typename T, int N = 0>
class RealVec : public ETVEC(ExprT,RealVec,T,N) {
public:

  /// Default constructor
  RealVec() : ETVEC(ExprT,RealVec,T,N)(this) {this->initialize(T());}

  /// Copy constructor from expression
  template <typename EXPR>
  RealVec(EETYPEV(EXPR) expr) : ETVEC(ExprT,RealVec,T,N)(this) {EVECLOOP(i, 0, N, m_data[i] = expr.at(i));}

  /// Copy constructor from RealVec
   RealVec(const RealVec<T,N>& orig) : ETVEC(ExprT,RealVec,T,N)(this) {copy(orig,*this);}

  /// Default destructor
   ~RealVec() {}

  /// Overloading of assignment operator(s)
#define VEC_EASSIGN_OP(__op__) \
  template <typename EXPR>						\
   const RealVec<T,N>& operator __op__ (EETYPEV(EXPR) expr) {EVECLOOP(i, 0, N, m_data[i] __op__ expr.at(i)); return *this;}
VEC_EASSIGN_OP(=)
VEC_EASSIGN_OP(+=)
VEC_EASSIGN_OP(-=)
VEC_EASSIGN_OP(*=)
VEC_EASSIGN_OP(/=)
#undef VEC_EASSIGN_OP

  /// Overloading of assignment operator(s) with constants
#define VEC_EASSIGN_OP_CONST(__op__) \
   const RealVec<T,N>& operator __op__ (T expr) {EVECLOOP(i, 0, N, m_data[i] __op__ expr); return *this;}
VEC_EASSIGN_OP_CONST(=)
VEC_EASSIGN_OP_CONST(+=)
VEC_EASSIGN_OP_CONST(-=)
VEC_EASSIGN_OP_CONST(*=)
VEC_EASSIGN_OP_CONST(/=)
#undef VEC_EASSIGN_OP_CONST

  /// copy content of another array
  template <typename ARRAY>
   void copyFrom(const ARRAY& in) {EVECLOOP(i, 0, N, m_data[i] = in[i]);}

  /// copy content of another array
  template <typename ARRAY>
   void copyTo(ARRAY& in) {EVECLOOP(i, 0, N, in[i] = m_data[i]);}

  /// @return a vector slice with fixed size
  template <int NV>
   RealVecSlice<T,NV> slice(std::size_t start) {return RealVecSlice<T,NV>(&m_data[start]);}

  /// @return a vector slice
   RealVecSlice<T,0> slice(std::size_t start, std::size_t ns) {return RealVecSlice<T,0>(&m_data[start], ns);}

  /// @return the raw data
   T* ptr() {return &m_data[0];}

  /// return the array size
   std::size_t size() const {return N;}

  /// Overloading of the stream operator "<<" for the output.
  /// "\n"ine introduced at the end of every line of the matrix.
  friend std::ostream& operator<<
    LTGT (std::ostream& out, const RealVec<T,N>& v);

  /// Overloading of the stream operator ">>" for the input
  friend std::istream& operator>>
    LTGT (std::istream& in, RealVec<T,N>& v);

  /// Overloading of the "==" operator.
  /// @return true if all elements are equal elementwise
   friend bool operator== LTGT (const RealVec<T,N>& v1, const RealVec<T,N>& v2);

  /// Overloading of the "==" operator.
  /// @param v given array
  /// @param value value for the comparison
  /// @return true if all elements are equal to value
   friend bool operator== LTGT (const RealVec<T,N>& v1, T value);

  /// Overloading of the "!=" operator.
  /// @return true if all elements are not equal elementwise
   friend bool operator!= LTGT (const RealVec<T,N>& v1, const RealVec<T,N>& v2);

  /// Overloading of the "!=" operator.
  /// @param v given array
  /// @param value value for the comparison
  /// @return true if at least one element is not equal to value
   friend bool operator!= LTGT (const RealVec<T,N>& v1, T value);

#ifndef CF_HAVE_CUDA
  /// Overloading of the "<" operator.
  /// @return true if the norm of the first RealVec is < than the norm of the second RealVec.
   friend bool operator< LTGT (const RealVec<T,N>& v1, const RealVec<T,N>& v2);

  /// Overloading of the "<" operator.
  /// @return true if all entries of the first RealVec are < than the given value
   friend bool operator< LTGT (const RealVec<T,N>& v1, T value);
#endif

  /// Overloading of the ">" operator.
  /// @return true if the norm of the first RealVec is > than the norm of the second RealVec.
   friend bool operator> LTGT (const RealVec<T,N>& v1, const RealVec<T,N>& v2);

  /// Overloading of the ">" operator.
  /// @return true if all entries of the first RealVec are > than the given value
   friend bool operator> LTGT (const RealVec<T,N>& v1, T value);

  /// Overloading of the "<=" operator.
  /// @return true if the norm of the first RealVec is <= than the norm of the second RealVec.
   friend bool operator<= LTGT (const RealVec<T,N>& v1, const RealVec<T,N>& v2);

  /// Overloading of the "<=" operator.
  /// @return true if all entries of the first RealVec are <= than the given value
   friend bool operator<= LTGT (const RealVec<T,N>& v1, T value);

  /// Overloading of the ">=" operator.
  /// @return true if the norm of the first RealVec is >= than the norm of the second RealVec.
   friend bool operator>= LTGT (const RealVec<T,N>& v1, const RealVec<T,N>& v2);

  /// Overloading of the ">=" operator.
  /// @return true if all entries of the first RealVec are >= than the given value
   friend bool operator>= LTGT (const RealVec<T,N>& v1, T value);

  /// Copy one RealVec into another one
  /// @pre v1.size() == v2.size()
  /// @param v1 source vector
  /// @param v2 destination vector
   friend void copy LTGT (const RealVec<T,N>& orig, RealVec<T,N>& dest);

  /// Normalizes this vector so that the norm is equal to one.
  /// \f$ v = v/\|v\|\f$
   void normalize() {*this /= this->norm2();}

  /// Projection of one vector onto another.
  /// @param v1 1st CFRealVector
  /// @param v2 2nd CFRealVector
  /// @post this RealVec contains the projected vector of v2 onto v1.
  template <typename V1, typename V2>
   void proj(const V1& v1, const V2& v2)
  {
    assert(v1.size() = v2.size());
    T scale = (MathFunctions::innerProd(v1, v2) / v1.sqrNorm());
    *this = v1 * scale;
  }

  /// Projection of this CFRealVector onto a given one.
  /// @pre this and the given object must be of the same size.
  /// @param v1 the other CFRealVector
  /// @post this CFRealVector contains the projected vector of itself onto v1.
  template <typename V1>
   void proj(const V1& v1) { proj(*this,v1);}

private:

  /// array data
  T m_data[N];
};

//////////////////////////////////////////////////////////////////////////////

/// Definition of a class RealVec that implements an expression template technique
/// with dynamical size
/// @author Andrea Lani
template <typename T>
class RealVec<T,0> : public ETVEC(ExprT,RealVec,T,0) {
public:

  /// Constructor
   RealVec(T value, std::size_t ns) : ETVEC(ExprT,RealVec,T,0)(this), m_owner(true), m_size(ns), m_data(NULL)
  {
    allocate();
    this->initialize(value);
  }

  /// Constructor
   RealVec(std::size_t ns = 0) : ETVEC(ExprT,RealVec,T,0)(this), m_owner(true), m_size(ns), m_data(NULL)
  {
   if (ns > 0) {
      allocate();
      this->initialize(T());
    }
  }

  // Constructor from preallocated memory
   RealVec(std::size_t ns, T* data) :
    ETVEC(ExprT,RealVec,T,0)(this), m_owner(false), m_size(ns), m_data(data) {}

  /// Copy constructor from expression
  template <typename EXPR>
   RealVec(EETYPEV(EXPR) expr) : ETVEC(ExprT,RealVec,T,0)(this)
  {
    m_owner = true; m_size = expr.size(); allocate();
    EVECLOOP(i, 0, EGETSIZE1(ENMAX(0,EXPR::SIZE1)), m_data[i] = expr.at(i));
  }

  /// Copy constructor from RealVec
   RealVec(const RealVec<T,0>& orig) :
    ETVEC(ExprT,RealVec,T,0)(this), m_owner(orig.m_owner), m_size(orig.m_size)
  {
    if (m_owner) {allocate(); copy(orig,*this);}
    else {m_data = orig.m_data;}
  }

  /// Default destructor
   ~RealVec() {free();}

  /// Tell if the array allocate its own memory
  bool isMemoryOwner() const {return m_owner;}

  /// Overloading of assignment operator(s)
#define VEC0_EASSIGN_OP(__op__) \
  template <typename EXPR>						\
   const RealVec<T,0>& operator __op__ (EETYPEV(EXPR) expr)	\
    {EVECLOOP(i, 0, EGETSIZE1(ENMAX(0,EXPR::SIZE1)), m_data[i] __op__ expr.at(i)); return *this;}
VEC0_EASSIGN_OP(=)
VEC0_EASSIGN_OP(+=)
VEC0_EASSIGN_OP(-=)
VEC0_EASSIGN_OP(*=)
VEC0_EASSIGN_OP(/=)
#undef VEC0_EASSIGN_OP

   /// Overloading of assignment operator(s) with constants
#define VEC0_EASSIGN_OP_CONST(__op__) \
   const RealVec<T,0>& operator __op__ (T expr) {EVECLOOP(i, 0, size(), m_data[i] __op__ expr); return *this;}
VEC0_EASSIGN_OP_CONST(=)
VEC0_EASSIGN_OP_CONST(+=)
VEC0_EASSIGN_OP_CONST(-=)
VEC0_EASSIGN_OP_CONST(*=)
VEC0_EASSIGN_OP_CONST(/=)
#undef VEC0_EASSIGN_OP_CONST

  /// Overloading of the assignment operator "="
  /// @pre the assignee is supposed to have same size and ownership
  /// as the given RealVec (as in std::valarray).
   const RealVec<T,0>& operator= (const RealVec<T,0>& other)
  {
    // cf_assert(&other != this);
    copy(other,*this);
    return *this;
  }

  // Constructor from preallocated memory
   void wrap(std::size_t ns, T* data) {m_owner = false; m_size = ns; m_data = data;}

  /// return the array size
   std::size_t size() const {return m_size;}

  /// Overloading of the stream operator "<<" for the output.
  /// "\n"ine introduced at the end of every line of the matrix.
  friend std::ostream& operator<<
    LTGT (std::ostream& out, const RealVec<T,0>& v);

  /// Overloading of the stream operator ">>" for the input
  friend std::istream& operator>>
    LTGT (std::istream& in, RealVec<T,0>& v);

  /// Overloading of the "==" operator.
  /// @return true if all elements are equal elementwise
   friend bool operator== LTGT (const RealVec<T,0>& v1, const RealVec<T,0>& v2);

  /// Overloading of the "==" operator.
  /// @param v given array
  /// @param value value for the comparison
  /// @return true if all elements are equal to value
   friend bool operator== LTGT (const RealVec<T,0>& v1, T value);

  /// Overloading of the "!=" operator.
  /// @return true if all elements are not equal elementwise
   friend bool operator!= LTGT (const RealVec<T,0>& v1, const RealVec<T,0>& v2);

  /// Overloading of the "!=" operator.
  /// @param v given array
  /// @param value value for the comparison
  /// @return true if at least one element is not equal to value
   friend bool operator!= LTGT (const RealVec<T,0>& v1, T value);

#ifndef CF_HAVE_CUDA
  /// Overloading of the "<" operator.
  /// @return true if the norm of the first RealVec is < than the norm of the second RealVec.
   friend bool operator< LTGT (const RealVec<T,0>& v1, const RealVec<T,0>& v2);

  /// Overloading of the "<" operator.
  /// @return true if all entries of the first RealVec are < than the given value
   friend bool operator< LTGT (const RealVec<T,0>& v1, T value);
#endif

  /// Overloading of the ">" operator.
  /// @return true if the norm of the first RealVec is > than the norm of the second RealVec.
   friend bool operator> LTGT (const RealVec<T,0>& v1, const RealVec<T,0>& v2);

  /// Overloading of the ">" operator.
  /// @return true if all entries of the first RealVec are > than the given value
   friend bool operator> LTGT (const RealVec<T,0>& v1, T value);

  /// Overloading of the "<=" operator.
  /// @return true if the norm of the first RealVec is <= than the norm of the second RealVec.
   friend bool operator<= LTGT (const RealVec<T,0>& v1, const RealVec<T,0>& v2);

  /// Overloading of the "<=" operator.
  /// @return true if all entries of the first RealVec are <= than the given value
   friend bool operator<= LTGT (const RealVec<T,0>& v1, T value);

  /// Overloading of the ">=" operator.
  /// @return true if the norm of the first RealVec is >= than the norm of the second RealVec.
   friend bool operator>= LTGT (const RealVec<T,0>& v1, const RealVec<T,0>& v2);

  /// Overloading of the ">=" operator.
  /// @return true if all entries of the first RealVec are >= than the given value
   friend bool operator>= LTGT (const RealVec<T,0>& v1, T value);

  /// Copy one RealVec into another one
  /// @pre v1.size() == v2.size()
  /// @param v1 source vector
  /// @param v2 destination vector
   friend void copy LTGT (const RealVec<T,0>& orig, RealVec<T,0>& dest);

  /// @return a vector slice with fixed size
  template <int NV>
   RealVecSlice<T,NV> slice(std::size_t start) {return RealVecSlice<T,NV>(&m_data[start]);}

  /// @return a vector slice
   RealVecSlice<T,0> slice(std::size_t start, std::size_t ns) {return RealVecSlice<T,0>(&m_data[start], ns);}

  /// @return the raw data
   T* ptr() {return m_data;}

  /// resize the array
   void resize(std::size_t ns, T init = T())
  {free(); m_size = ns; allocate(); this->initialize(init);}

  /// Normalizes this vector so that the norm is equal to one.
  /// \f$ v = v/\|v\|\f$
   void normalize() {*this /= this->norm2();}

  /// Projection of one vector onto another.
  /// @param v1 1st CFRealVector
  /// @param v2 2nd CFRealVector
  /// @post this RealVec contains the projected vector of v2 onto v1.
  template <typename V1, typename V2>
   void proj(const V1& v1, const V2& v2)
  {
    assert(v1.size() == v2.size());
    T scale = (MathFunctions::innerProd(v1, v2) / v1.sqrNorm());
    *this = v1 * scale;
  }

  /// Projection of this CFRealVector onto a given one.
  /// @pre this and the given object must be of the same size.
  /// @param v1 the other CFRealVector
  /// @post this CFRealVector contains the projected vector of itself onto v1.
  template <typename V1>
   void proj(const V1& v1) { proj(*this,v1); }

private: // helper functions

  /// allocate the memory
   void allocate() {if (size() > 0) {assert(m_owner); m_data = new T[m_size];}}

  /// free the memory
   void free() {if (m_owner && m_size > 0) {delete [] m_data; m_data = NULL; m_size = 0;}}

private:

  /// boolean flag to know if RealVec has ownership on the data
  bool m_owner;

  /// array size
  std::size_t m_size;

  /// array data
  T* m_data;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
std::ostream& operator<< (std::ostream& out, const RealVec<T,N>& v)
{
  const std::size_t size = v.size();
  for (std::size_t i = 0; i < size; ++i)
    out << v.m_data[i] << " " ;
  return out;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
std::istream& operator>> (std::istream& in, RealVec<T,N>& v)
{
  const std::size_t size = v.size();
  for (std::size_t i = 0; i < size; ++i)
    in >> v.m_data[i];
  return in;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
bool operator== (const RealVec<T,N>& v1, const RealVec<T,N>& v2)
{
  assert(v1.size() == v2.size());
  const std::size_t size = v1.size();
  for (std::size_t i = 0; i < size; ++i) {
    if (v1.m_data[i] != v2.m_data[i]) {
      return false;
    }
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
bool operator== (const RealVec<T,N>& v, T value)
{
  const std::size_t size = v.size();
  for (std::size_t i = 0; i < size; ++i) {
    if (v[i] != value)  return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
bool operator!= (const RealVec<T,N>& v1, const RealVec<T,N>& v2)
{
  return !(v1 == v2);
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
bool operator!= (const RealVec<T,N>& v, T value)
{
  return !(v == value);
}

//////////////////////////////////////////////////////////////////////////////

#ifndef CF_HAVE_CUDA
template <typename T, int N>
bool operator< (const RealVec<T,N>& v1, const RealVec<T,N>& v2)
{
  return (v1.sqrNorm() < v2.sqrNorm());
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
bool operator< (const RealVec<T,N>& v1, const T value)
{
 const std::size_t size = v1.size();
 for (std::size_t i = 0; i < size; ++i) {
   if (v1[i] >= value) return false;
 }
 return true;
}
#endif

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
bool operator> (const RealVec<T,N>& v1, const RealVec<T,N>& v2)
{
  return v1.sqrNorm() > v2.sqrNorm();
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
bool operator> (const RealVec<T,N>& v1, T value)
{
  const std::size_t size = v1.size();
  for (std::size_t i = 0; i < size; ++i) {
    if (v1[i] <= value) return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
bool operator<= (const RealVec<T,N>& v1, const RealVec<T,N>& v2)
{
  return !(v1.sqrNorm() > v2.sqrNorm());
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
bool operator<= (const RealVec<T,N>& v1, T value)
{
  const std::size_t size = v1.size();
  for (std::size_t i = 0; i < size; ++i) {
    if (v1[i] > value) return false;
  }
 return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
bool operator>= (const RealVec<T,N>& v1, const RealVec<T,N>& v2)
{
  return !(v1.sqrNorm() < v2.sqrNorm());
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
bool operator>= (const RealVec<T,N>& v1, T value)
{
  const std::size_t size = v1.size();
  for (std::size_t i = 0; i < size; ++i) {
    if (v1[i] < value) return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
void copy (const RealVec<T,N>& orig, RealVec<T,N>& dest)
{
  static_assert(orig.size() == dest.size(),"Length of two vectors for copy doesn't match");
  const std::size_t size = orig.size();
  for (std::size_t i = 0; i < size; ++i)
    dest.m_data[i] = orig.m_data[i];
}

//////////////////////////////////////////////////////////////////////////////

}; /*--- End of class RealVec ---*/

} /*--- End of namespace Common ---*/

//////////////////////////////////////////////////////////////////////////////

#endif
